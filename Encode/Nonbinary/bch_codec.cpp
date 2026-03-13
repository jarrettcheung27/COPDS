#include "bch_codec.h"

#include <stdexcept>
#include <algorithm>

namespace
{
    std::vector<uint8_t> poly_mod2_div_high(const std::vector<uint8_t> &dividend, const std::vector<uint8_t> &divisor)
    {
        if (divisor.empty() || divisor[0] == 0)
        {
            throw std::runtime_error("invalid divisor polynomial");
        }
        std::vector<uint8_t> work = dividend;
        const int n = static_cast<int>(dividend.size());
        const int m = static_cast<int>(divisor.size());
        for (int i = 0; i <= n - m; ++i)
        {
            if (work[i])
            {
                for (int j = 0; j < m; ++j)
                {
                    work[i + j] ^= divisor[j];
                }
            }
        }
        return std::vector<uint8_t>(work.end() - (m - 1), work.end());
    }
}

BinaryBCHCodec::BinaryBCHCodec(int n, int k)
    : n_short_(n), k_short_(k), n_full_(0), k_full_(0), m_(0), t_(0), shorten_(0), parity_bits_(0), prim_poly_(0)
{
    if (n <= 0 || k <= 0 || k >= n)
    {
        throw std::runtime_error("invalid BCH (n, k)");
    }

    m_ = deduce_m(n, k);
    n_full_ = (1 << m_) - 1;
    parity_bits_ = n - k;
    k_full_ = n_full_ - parity_bits_;
    shorten_ = n_full_ - n;
    if (k_full_ - shorten_ != k_short_)
    {
        throw std::runtime_error("unsupported shortened BCH parameters");
    }

    t_ = deduce_t(m_, n, k);
    prim_poly_ = primitive_poly(m_);

    init_field();
    init_generator();
}

int BinaryBCHCodec::deduce_m(int n, int k)
{
    const int parity = n - k;
    for (int m = 2; m <= 15; ++m)
    {
        const int nfull = (1 << m) - 1;
        const int kfull = nfull - parity;
        const int shorten = nfull - n;
        if (shorten >= 0 && kfull - shorten == k)
        {
            return m;
        }
    }
    throw std::runtime_error("cannot deduce BCH order m from (n, k)");
}

int BinaryBCHCodec::deduce_t(int m, int n, int k)
{
    const int parity = n - k;
    if (m == 5 && parity == 15)
    {
        return 3;
    }
    if (m == 9 && parity == 99)
    {
        return 11;
    }
    const int t = parity / m;
    if (t <= 0)
    {
        throw std::runtime_error("cannot deduce BCH correction capability t");
    }
    return t;
}

int BinaryBCHCodec::primitive_poly(int m)
{
    switch (m)
    {
    case 5:
        return 0x25; // x^5 + x^2 + 1
    case 9:
        return 0x211; // x^9 + x^4 + 1
    default:
        throw std::runtime_error("unsupported BCH field order");
    }
}

void BinaryBCHCodec::init_field()
{
    const int nn = n_full_;
    alpha_to_.assign(nn * 2 + 1, 0);
    index_of_.assign(nn + 1, -1);

    alpha_to_[0] = 1;
    for (int i = 1; i < nn; ++i)
    {
        alpha_to_[i] = alpha_to_[i - 1] << 1;
        if (alpha_to_[i] & (1 << m_))
        {
            alpha_to_[i] ^= prim_poly_;
        }
        alpha_to_[i] &= nn;
    }

    for (int i = 0; i < nn; ++i)
    {
        index_of_[alpha_to_[i]] = i;
    }
    index_of_[0] = -1;

    for (int i = nn; i <= 2 * nn; ++i)
    {
        alpha_to_[i] = alpha_to_[i - nn];
    }
}

int BinaryBCHCodec::gf_mul(int a, int b) const
{
    if (a == 0 || b == 0)
    {
        return 0;
    }
    const int nn = n_full_;
    return alpha_to_[(index_of_[a] + index_of_[b]) % nn];
}

int BinaryBCHCodec::gf_div(int a, int b) const
{
    if (b == 0)
    {
        throw std::runtime_error("GF division by zero");
    }
    if (a == 0)
    {
        return 0;
    }
    const int nn = n_full_;
    int idx = index_of_[a] - index_of_[b];
    idx %= nn;
    if (idx < 0)
    {
        idx += nn;
    }
    return alpha_to_[idx];
}

int BinaryBCHCodec::gf_pow_alpha(int p) const
{
    const int nn = n_full_;
    int idx = p % nn;
    if (idx < 0)
    {
        idx += nn;
    }
    return alpha_to_[idx];
}

void BinaryBCHCodec::init_generator()
{
    const int nn = n_full_;
    std::vector<int> include(nn, 0);
    for (int e = 1; e <= 2 * t_; ++e)
    {
        int cur = e % nn;
        do
        {
            include[cur] = 1;
            cur = (cur * 2) % nn;
        } while (cur != (e % nn));
    }

    std::vector<int> poly(1, 1);
    for (int e = 1; e < nn; ++e)
    {
        if (!include[e])
        {
            continue;
        }
        const int root = gf_pow_alpha(e);
        std::vector<int> next(poly.size() + 1, 0);
        for (size_t i = 0; i < poly.size(); ++i)
        {
            next[i] = gf_add(next[i], gf_mul(poly[i], root));
            next[i + 1] = gf_add(next[i + 1], poly[i]);
        }
        poly.swap(next);
    }

    if (static_cast<int>(poly.size()) != parity_bits_ + 1)
    {
        throw std::runtime_error("generator polynomial degree mismatch");
    }

    std::vector<uint8_t> gen_low(poly.size(), 0);
    for (size_t i = 0; i < poly.size(); ++i)
    {
        if (poly[i] != 0 && poly[i] != 1)
        {
            throw std::runtime_error("generator polynomial has non-binary coefficients");
        }
        gen_low[i] = static_cast<uint8_t>(poly[i]);
    }

    generator_high_.assign(gen_low.rbegin(), gen_low.rend());
}

std::vector<uint8_t> BinaryBCHCodec::encode_full_systematic(const std::vector<uint8_t> &msg_full) const
{
    if (static_cast<int>(msg_full.size()) != k_full_)
    {
        throw std::runtime_error("full-message length mismatch");
    }

    std::vector<uint8_t> dividend(n_full_, 0);
    for (int i = 0; i < k_full_; ++i)
    {
        dividend[i] = static_cast<uint8_t>(msg_full[i] & 1);
    }

    auto remainder = poly_mod2_div_high(dividend, generator_high_);

    std::vector<uint8_t> codeword(n_full_, 0);
    for (int i = 0; i < k_full_; ++i)
    {
        codeword[i] = static_cast<uint8_t>(msg_full[i] & 1);
    }
    for (int i = 0; i < parity_bits_; ++i)
    {
        codeword[k_full_ + i] = remainder[i];
    }
    return codeword;
}

std::vector<int> BinaryBCHCodec::compute_syndromes(const std::vector<uint8_t> &rx_full) const
{
    std::vector<int> syndromes(2 * t_ + 1, 0);
    for (int s = 1; s <= 2 * t_; ++s)
    {
        int Sj = 0;
        for (int i = 0; i < n_full_; ++i)
        {
            if (rx_full[i])
            {
                const int power = s * (n_full_ - 1 - i);
                Sj = gf_add(Sj, gf_pow_alpha(power));
            }
        }
        syndromes[s] = Sj;
    }
    return syndromes;
}

std::vector<int> BinaryBCHCodec::berlekamp_massey(const std::vector<int> &syndromes) const
{
    const int max_deg = 2 * t_;
    std::vector<int> C(max_deg + 1, 0);
    std::vector<int> B(max_deg + 1, 0);
    C[0] = 1;
    B[0] = 1;

    int L = 0;
    int m = 1;
    int b = 1;

    for (int n = 0; n < max_deg; ++n)
    {
        int d = syndromes[n + 1];
        for (int i = 1; i <= L; ++i)
        {
            d = gf_add(d, gf_mul(C[i], syndromes[n + 1 - i]));
        }

        if (d == 0)
        {
            ++m;
            continue;
        }

        std::vector<int> T = C;
        int coef = gf_div(d, b);
        for (int i = 0; i + m <= max_deg; ++i)
        {
            C[i + m] = gf_add(C[i + m], gf_mul(coef, B[i]));
        }

        if (2 * L <= n)
        {
            L = n + 1 - L;
            B = T;
            b = d;
            m = 1;
        }
        else
        {
            ++m;
        }
    }

    C.resize(L + 1);
    return C;
}

std::vector<int> BinaryBCHCodec::chien_search(const std::vector<int> &sigma) const
{
    std::vector<int> positions;
    const int degree = static_cast<int>(sigma.size()) - 1;
    if (degree <= 0)
    {
        return positions;
    }

    for (int pos = 0; pos < n_full_; ++pos)
    {
        int xinv = gf_pow_alpha(-pos);
        int xpow = 1;
        int val = sigma[0];
        for (int i = 1; i <= degree; ++i)
        {
            xpow = gf_mul(xpow, xinv);
            val = gf_add(val, gf_mul(sigma[i], xpow));
        }
        if (val == 0)
        {
            positions.push_back(pos);
        }
    }
    return positions;
}

std::pair<std::vector<uint8_t>, int> BinaryBCHCodec::decode_full_message(const std::vector<uint8_t> &rx_full) const
{
    if (static_cast<int>(rx_full.size()) != n_full_)
    {
        throw std::runtime_error("full-codeword length mismatch");
    }

    auto syndromes = compute_syndromes(rx_full);
    bool has_error = false;
    for (int s = 1; s <= 2 * t_; ++s)
    {
        if (syndromes[s] != 0)
        {
            has_error = true;
            break;
        }
    }

    std::vector<uint8_t> corrected = rx_full;
    int flag = 0;

    if (has_error)
    {
        auto sigma = berlekamp_massey(syndromes);
        const int sigma_deg = static_cast<int>(sigma.size()) - 1;
        auto positions = chien_search(sigma);

        if (sigma_deg <= 0 || sigma_deg > t_ || static_cast<int>(positions.size()) != sigma_deg)
        {
            flag = -1;
        }
        else
        {
            for (int p : positions)
            {
                const int idx = n_full_ - 1 - p;
                if (idx >= 0 && idx < n_full_)
                {
                    corrected[idx] ^= 1;
                }
            }
            flag = static_cast<int>(positions.size());

            auto verify = compute_syndromes(corrected);
            for (int s = 1; s <= 2 * t_; ++s)
            {
                if (verify[s] != 0)
                {
                    flag = -1;
                    break;
                }
            }
        }
    }

    std::vector<uint8_t> msg(k_full_, 0);
    for (int i = 0; i < k_full_; ++i)
    {
        msg[i] = corrected[i];
    }
    return {msg, flag};
}

std::vector<uint8_t> BinaryBCHCodec::encode(const std::vector<uint8_t> &msg) const
{
    if (static_cast<int>(msg.size()) != k_short_)
    {
        throw std::runtime_error("message length mismatch for shortened BCH encode");
    }

    std::vector<uint8_t> msg_full(k_full_, 0);
    for (int i = 0; i < k_short_; ++i)
    {
        msg_full[shorten_ + i] = static_cast<uint8_t>(msg[i] & 1);
    }

    auto code_full = encode_full_systematic(msg_full);
    std::vector<uint8_t> code_short(n_short_, 0);
    for (int i = 0; i < n_short_; ++i)
    {
        code_short[i] = code_full[shorten_ + i];
    }
    return code_short;
}

std::pair<std::vector<uint8_t>, int> BinaryBCHCodec::decode(const std::vector<uint8_t> &rx) const
{
    if (static_cast<int>(rx.size()) != n_short_)
    {
        throw std::runtime_error("codeword length mismatch for shortened BCH decode");
    }

    std::vector<uint8_t> rx_full(n_full_, 0);
    for (int i = 0; i < n_short_; ++i)
    {
        rx_full[shorten_ + i] = static_cast<uint8_t>(rx[i] & 1);
    }

    auto dec_full = decode_full_message(rx_full);
    std::vector<uint8_t> msg_short(k_short_, 0);
    for (int i = 0; i < k_short_; ++i)
    {
        msg_short[i] = dec_full.first[shorten_ + i];
    }
    return {msg_short, dec_full.second};
}
