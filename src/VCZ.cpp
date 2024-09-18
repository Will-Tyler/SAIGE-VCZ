#include "VCZ.hpp"

#include <string>

#include "tensorstore/array.h"
#include "tensorstore/chunk_layout.h"
#include "tensorstore/index_space/dim_expression.h"
#include "tensorstore/open.h"
#include "tensorstore/spec.h"
#include "tensorstore/tensorstore.h"

using ::tensorstore::Index;
using ::tensorstore::TensorStore;

namespace VCZ {

class VczClass::Impl {
   public:
    tensorstore::TensorStore<int8_t, tensorstore::dynamic_rank,
                             tensorstore::ReadWriteMode::read>
        m_store;
    std::string m_chrom;
    int m_pos;
    int m_startPos;
    int m_endPos;
};

VczClass::VczClass(std::string t_vczFileName,
                   std::vector<std::string> t_SampleInModel) {
    std::string path = t_vczFileName + "/call_genotype";
    tensorstore::Spec input_spec =
        tensorstore::Spec::FromJson(
            {
                {"driver", "zarr"},
                {"kvstore", {{"driver", "file"}, {"path", path}}},
            })
            .value();
    pimpl->m_store = tensorstore::Open<int8_t, tensorstore::dynamic_rank,
                                       tensorstore::ReadWriteMode::read>(
                         input_spec, tensorstore::OpenMode::open,
                         tensorstore::ReadWriteMode::read)
                         .value();
    pimpl->m_chrom = "";
    pimpl->m_pos = 1;
    pimpl->m_startPos = 1;
    pimpl->m_endPos = std::numeric_limits<int>::max();
}

void VczClass::set_iterator(std::string& chrom, const int start_pos,
                            const int end_pos) {
    pimpl->m_chrom = chrom;
    pimpl->m_startPos = start_pos;
    pimpl->m_endPos = end_pos;
    pimpl->m_pos = start_pos;
}

void VczClass::move_forward_iterator(const int a) { pimpl->m_pos += a; }

bool VczClass::check_iterator_end() {
    // TODO
    return true;
}

}  // namespace VCZ
