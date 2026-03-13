#pragma once

#include <vector>
#include <cstdint>

class BinaryBCHCodec
{
public:
    BinaryBCHCodec(int n, int k);

    int n() const { return n_short_; }
    int k() const { return k_short_; }

    std::vector<uint8_t> encode(const std::vector<uint8_t> &msg) const;
    std::pair<std::vector<uint8_t>, int> decode(const std::vector<uint8_t> &rx) const;

private:
    int n_short_;
    int k_short_;

    int n_full_;
    int k_full_;
    int m_;
    int t_;
    int shorten_;
    int parity_bits_;

    int prim_poly_;

    std::vector<int> alpha_to_;
    std::vector<int> index_of_;
    std::vector<uint8_t> generator_high_;

    static int deduce_m(int n, int k);
    static int deduce_t(int m, int n, int k);
    static int primitive_poly(int m);

    void init_field();
    void init_generator();

    int gf_add(int a, int b) const { return a ^ b; }
    int gf_mul(int a, int b) const;
    int gf_div(int a, int b) const;
    int gf_pow_alpha(int p) const;

    std::vector<uint8_t> encode_full_systematic(const std::vector<uint8_t> &msg_full) const;
    std::pair<std::vector<uint8_t>, int> decode_full_message(const std::vector<uint8_t> &rx_full) const;

    std::vector<int> compute_syndromes(const std::vector<uint8_t> &rx_full) const;
    std::vector<int> berlekamp_massey(const std::vector<int> &syndromes) const;
    std::vector<int> chien_search(const std::vector<int> &sigma) const;
};
