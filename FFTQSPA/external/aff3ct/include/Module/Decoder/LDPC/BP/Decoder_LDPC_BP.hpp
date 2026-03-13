/*!
 * \file
 * \brief Class module::Decoder_LDPC_BP.
 */
#ifndef DECODER_LDPC_BP_HPP_
#define DECODER_LDPC_BP_HPP_

#include "Tools/Algo/Matrix/Sparse_matrix/Sparse_matrix.hpp"
#include "Tools/Code/LDPC/Syndrome/LDPC_syndrome.hpp"

namespace aff3ct
{
namespace module
{
class Decoder_LDPC_BP
{
  protected:
    const int n_ite;
    const tools::Sparse_matrix H;
    const bool enable_syndrome;
    const int syndrome_depth;

    int cur_syndrome_depth;

  public:
    Decoder_LDPC_BP(const int K,
                    const int N,
                    const int n_ite,
                    const tools::Sparse_matrix& H,
                    const bool enable_syndrome = true,
                    const int syndrome_depth = 1);

    virtual ~Decoder_LDPC_BP() = default;

    template<typename R, class Syndrome_checker = tools::LDPC_syndrome>
    inline bool check_syndrome_soft(const R* Y_N);

    template<typename B, class Syndrome_checker = tools::LDPC_syndrome>
    inline bool check_syndrome_hard(const B* V_N);
};
}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Decoder/LDPC/BP/Decoder_LDPC_BP.hxx"
#endif

#endif /* DECODER_LDPC_BP_HPP_ */
