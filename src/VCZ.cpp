#include <string>

#include "tensorstore/array.h"
#include "tensorstore/chunk_layout.h"
#include "tensorstore/index_space/dim_expression.h"
#include "tensorstore/open.h"
#include "tensorstore/spec.h"
#include "tensorstore/tensorstore.h"
#include "VCZ.hpp"

namespace VCZ {

VczClass::VczClass(std::string t_vczFileName, std::vector<std::string> t_SampleInModel) {
    std::string path = t_vczFileName + "/call_genotype";
    tensorstore::Spec input_spec =
        tensorstore::Spec::FromJson(
            {
                {"driver", "zarr"},
                {"kvstore", {{"driver", "file"}, {"path", path}}},
            })
            .value();
}

}  // namespace VCZ
