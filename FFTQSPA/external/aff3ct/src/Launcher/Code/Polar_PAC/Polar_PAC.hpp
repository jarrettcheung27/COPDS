#ifndef LAUNCHER_POLAR_PAC_HPP_
#define LAUNCHER_POLAR_PAC_HPP_

#include <iostream>

#include "Factory/Tools/Codec/Polar_PAC/Codec_polar_PAC.hpp"

namespace aff3ct
{
namespace launcher
{
template<class L, typename B = int, typename R = float, typename Q = R>
class Polar_PAC : public L
{
  protected:
    factory::Codec_polar_PAC* params_cdc;

  public:
    Polar_PAC(const int argc, const char** argv, std::ostream& stream = std::cout);
    virtual ~Polar_PAC() = default;

  protected:
    virtual void get_description_args();
    virtual void store_args();
};
}
}

#endif /* LAUNCHER_POLAR_PAC_HPP_ */
