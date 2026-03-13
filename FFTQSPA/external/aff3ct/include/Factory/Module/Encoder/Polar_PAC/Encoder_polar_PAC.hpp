/*!
 * \file
 * \brief Class factory::Encoder_polar_PAC.
 */
#ifndef FACTORY_ENCODER_POLAR_PAC_HPP
#define FACTORY_ENCODER_POLAR_PAC_HPP

#include <cli.hpp>
#include <map>
#include <string>
#include <vector>

#include "Factory/Module/Encoder/Encoder.hpp"
#include "Module/Encoder/Polar_PAC/Encoder_polar_PAC.hpp"
#include "Tools/Factory/Header.hpp"

namespace aff3ct
{
namespace factory
{
extern const std::string Encoder_polar_PAC_name;
extern const std::string Encoder_polar_PAC_prefix;
class Encoder_polar_PAC : public Encoder
{
  public:
    // ----------------------------------------------------------------------------------------------------- PARAMETERS
    std::string conv = "NO";

    // -------------------------------------------------------------------------------------------------------- METHODS
    explicit Encoder_polar_PAC(const std::string& p = Encoder_polar_PAC_prefix);
    virtual ~Encoder_polar_PAC() = default;
    Encoder_polar_PAC* clone() const;

    // parameters construction
    void get_description(cli::Argument_map_info& args) const;
    void store(const cli::Argument_map_value& vals);
    void get_headers(std::map<std::string, tools::header_list>& headers, const bool full = true) const;

    // builder
    template<typename B = int>
    module::Encoder_polar_PAC<B>* build(const std::vector<bool>& frozen_bits) const;
};
}
}

#endif /* FACTORY_ENCODER_POLAR_PAC_HPP */
