
/*!
 * \file
 * \brief Class tools::Frozenbits_generator_RM.
 */
#ifndef FROZENBITS_GENERATOR_RM_HPP_
#define FROZENBITS_GENERATOR_RM_HPP_

#include "Tools/Code/Polar/Frozenbits_generator/Frozenbits_generator.hpp"

namespace aff3ct
{
namespace tools
{
class Frozenbits_generator_RM : public Frozenbits_generator
{
  public:
    Frozenbits_generator_RM(const int K,
                            const int N,
                            const std::string& dump_channels_path = "",
                            const bool dump_channels_single_thread = true);

    ~Frozenbits_generator_RM();

    virtual Frozenbits_generator_RM* clone() const;

  private:
    void evaluate();
};
}
}

#endif /* FROZENBITS_GENERATOR_RM_HPP_ */
