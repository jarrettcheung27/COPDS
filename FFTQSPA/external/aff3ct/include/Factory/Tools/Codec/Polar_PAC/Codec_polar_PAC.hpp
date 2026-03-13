/*!
 * \file
 * \brief Class factory::Codec_polar_PAC.
 */
#ifndef FACTORY_CODEC_POLAR_PAC_HPP
#define FACTORY_CODEC_POLAR_PAC_HPP

#include <cli.hpp>
#include <map>
#include <string>
#include <vector>

#include "Factory/Tools/Code/Polar/Frozenbits_generator.hpp"
#include "Factory/Tools/Codec/Codec_SISO.hpp"
#include "Module/CRC/CRC.hpp"
#include "Tools/Codec/Polar_PAC/Codec_polar_PAC.hpp"
#include "Tools/Factory/Header.hpp"

namespace aff3ct
{
namespace factory
{
extern const std::string Codec_polar_PAC_name;
extern const std::string Codec_polar_PAC_prefix;
class Codec_polar_PAC : public Codec_SISO
{
  public:
    // ----------------------------------------------------------------------------------------------------- PARAMETERS
    // depending parameters
    tools::auto_cloned_unique_ptr<Frozenbits_generator> fbg;

    // -------------------------------------------------------------------------------------------------------- METHODS
    explicit Codec_polar_PAC(const std::string& p = Codec_polar_PAC_prefix);
    virtual ~Codec_polar_PAC() = default;
    Codec_polar_PAC* clone() const;

    virtual std::vector<std::string> get_names() const;
    virtual std::vector<std::string> get_short_names() const;
    virtual std::vector<std::string> get_prefixes() const;

    // parameters construction
    void get_description(cli::Argument_map_info& args) const;
    void store(const cli::Argument_map_value& vals);
    void get_headers(std::map<std::string, tools::header_list>& headers, const bool full = true) const;

    // builder
    template<typename B = int, typename Q = float>
    tools::Codec_polar_PAC<B, Q>* build(const module::CRC<B>* crc = nullptr) const;
};
}
}

#endif /* FACTORY_CODEC_POLAR_PAC_HPP */
