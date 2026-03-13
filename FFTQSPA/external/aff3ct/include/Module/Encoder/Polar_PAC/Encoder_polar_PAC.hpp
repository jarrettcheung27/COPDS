/*!
 * \file
 * \brief Class module::Encoder_polar_PAC.
 */
#ifndef ENCODER_POLAR_PAC_HPP_
#define ENCODER_POLAR_PAC_HPP_

#include <vector>

#include "Module/Encoder/Encoder.hpp"
#include "Tools/Interface/Interface_get_set_frozen_bits.hpp"

namespace aff3ct
{
namespace module
{
template<typename B = int>
class Encoder_polar_PAC
  : public Encoder<B>
  , public tools::Interface_get_set_frozen_bits
{
  protected:
    const int m;                   // log_2 of code length
    std::vector<bool> frozen_bits; // true means frozen, false means set to 0/1
    std::vector<B> X_N_tmp;
    std::vector<uint8_t> conv_reg;
    std::vector<B> curState;

  public:
    Encoder_polar_PAC(const int& K, const int& N, const std::vector<bool>& frozen_bits, const std::string& conv);
    virtual ~Encoder_polar_PAC() = default;

    virtual Encoder_polar_PAC<B>* clone() const;

    void light_encode(B* bits);

    bool is_codeword(const B* X_N);

    virtual const std::vector<bool>& get_frozen_bits() const;
    virtual void set_frozen_bits(const std::vector<bool>& frozen_bits);

  protected:
    virtual void _encode(const B* U_K, B* X_N, const size_t frame_id);
    void convert(const B* U_K, B* U_N);

    B conv1bitEnc(B cbit);
    void convEnc(B* X_N);
};
}
}

#endif // ENCODER_POLAR_PAC_HPP_
