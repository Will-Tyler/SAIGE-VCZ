#include "VCZ.hpp"

#include <cassert>
#include <memory>
#include <string>
#include <vector>

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
   private:
    tensorstore::TensorStore<int64_t, tensorstore::dynamic_rank,
                             tensorstore::ReadWriteMode::read>
        m_contig_length_store;
    tensorstore::SharedArray<int64_t, tensorstore::dynamic_rank,
                             tensorstore::offset_origin>
        m_contig_length_array;

    tensorstore::TensorStore<int8_t, tensorstore::dynamic_rank,
                             tensorstore::ReadWriteMode::read>
        m_variant_contig_store;
    tensorstore::SharedArray<int8_t, tensorstore::dynamic_rank,
                             tensorstore::offset_origin>
        m_variant_contig_array;

    tensorstore::TensorStore<int8_t, tensorstore::dynamic_rank,
                             tensorstore::ReadWriteMode::read>
        m_variant_position_store;
    tensorstore::SharedArray<int8_t, tensorstore::dynamic_rank,
                             tensorstore::offset_origin>
        m_variant_position_array;

    tensorstore::TensorStore<int8_t, tensorstore::dynamic_rank,
                             tensorstore::ReadWriteMode::read>
        m_call_genotype_store;
    tensorstore::SharedArray<int8_t, tensorstore::dynamic_rank,
                             tensorstore::offset_origin>
        m_call_genotype_array;

    tensorstore::span<const Index> m_genotype_shape;
    std::string m_chrom;
    Index m_marker_index;
    int m_startPos;
    int m_endPos;

   public:
    Impl(std::string& t_vczFileName) {
        tensorstore::Spec input_spec =
            tensorstore::Spec::FromJson(
                {
                    {"driver", "zarr"},
                    {"kvstore",
                     {{"driver", "file"},
                      {"path", t_vczFileName + "/call_genotype"}}},
                })
                .value();
        m_call_genotype_store =
            tensorstore::Open<int8_t, tensorstore::dynamic_rank,
                              tensorstore::ReadWriteMode::read>(
                input_spec, tensorstore::OpenMode::open,
                tensorstore::ReadWriteMode::read)
                .value();
        m_genotype_shape = m_call_genotype_store.domain().shape();
        m_call_genotype_array =
            tensorstore::Read(m_call_genotype_store).value();

        input_spec = tensorstore::Spec::FromJson(
                         {
                             {"driver", "zarr"},
                             {"kvstore",
                              {{"driver", "file"},
                               {"path", t_vczFileName + "/contig_length"}}},
                         })
                         .value();
        m_contig_length_store =
            tensorstore::Open<int64_t, tensorstore::dynamic_rank,
                              tensorstore::ReadWriteMode::read>(
                input_spec, tensorstore::OpenMode::open,
                tensorstore::ReadWriteMode::read)
                .value();
        m_contig_length_array =
            tensorstore::Read(m_contig_length_store).value();

        input_spec = tensorstore::Spec::FromJson(
                         {
                             {"driver", "zarr"},
                             {"kvstore",
                              {{"driver", "file"},
                               {"path", t_vczFileName + "/variant_contig"}}},
                         })
                         .value();
        m_variant_contig_store =
            tensorstore::Open<int8_t, tensorstore::dynamic_rank,
                              tensorstore::ReadWriteMode::read>(
                input_spec, tensorstore::OpenMode::open,
                tensorstore::ReadWriteMode::read)
                .value();
        m_variant_contig_array =
            tensorstore::Read(m_variant_contig_store).value();

        input_spec = tensorstore::Spec::FromJson(
                         {
                             {"driver", "zarr"},
                             {"kvstore",
                              {{"driver", "file"},
                               {"path", t_vczFileName + "/variant_position"}}},
                         })
                         .value();
        m_variant_position_store =
            tensorstore::Open<int8_t, tensorstore::dynamic_rank,
                              tensorstore::ReadWriteMode::read>(
                input_spec, tensorstore::OpenMode::open,
                tensorstore::ReadWriteMode::read)
                .value();
        m_variant_position_array =
            tensorstore::Read(m_variant_position_store).value();

        m_chrom = "";
        m_marker_index = 0;
        m_startPos = 1;
        m_endPos = 1000;
    }

    void set_iterator(std::string& chrom, const int start_pos,
                      const int end_pos) {
        assert(start_pos > 0);
        assert(end_pos >= start_pos);
        m_chrom = chrom;
        m_startPos = start_pos;
        m_endPos = std::min(end_pos, (int)m_genotype_shape[0]);
        m_marker_index = start_pos - 1;
    }

    void move_forward_iterator(const int a) {
        assert(a >= 0);
        m_marker_index += a;
    }

    bool check_iterator_end() { return m_marker_index >= m_endPos; }

    void getOneMarker(std::string& t_ref, std::string& t_alt,
                      std::string& t_marker, uint32_t& t_pos,
                      std::string& t_chr, double& t_altFreq,
                      double& t_altCounts, double& t_missingRate,
                      double& t_imputeInfo, bool t_isOutputIndexForMissing,
                      std::vector<uint32_t>& t_indexForMissingforOneMarker,
                      bool t_isOnlyOutputNonZero,
                      std::vector<uint32_t>& t_indexForNonZero,
                      bool& t_isBoolRead, std::vector<double>& dosages,
                      bool t_isImputation) {
        if (check_iterator_end()) {
            t_isBoolRead = false;
            return;
        }
        
        t_ref = "A";
        t_alt = "B";
        t_marker = "ID";
        t_pos = m_variant_position_array({m_marker_index});
        t_chr = "1";
        t_altFreq = 0.0;
        t_altCounts = 0.0;
        t_imputeInfo = 1.0;

        size_t missing_count = 0;
        dosages.clear();
        dosages.resize(m_genotype_shape[1]);

        for (Index sample_index = 0; sample_index < m_genotype_shape[1];
             sample_index++) {
            const int8_t a =
                m_call_genotype_array({m_marker_index, sample_index, 0});
            const int8_t b =
                m_call_genotype_array({m_marker_index, sample_index, 1});

            // TODO: use m_posSampleInModel

            if (a >= 0 && b >= 0) {
                const double dosage = a + b;
                dosages[sample_index] = dosage;
                t_altCounts += dosage;
            } else {
                missing_count += 1;
            }
        }

        t_isBoolRead = true;
    }
};

VczClass::VczClass(std::string t_vczFileName,
                   std::vector<std::string> t_SampleInModel) {
    pimpl = std::make_unique<Impl>(t_vczFileName);
}

void VczClass::set_iterator(std::string& chrom, const int start_pos,
                            const int end_pos) {
    pimpl->set_iterator(chrom, start_pos, end_pos);
}

void VczClass::move_forward_iterator(const int a) {
    assert(a >= 0);
    pimpl->move_forward_iterator(a);
}

bool VczClass::check_iterator_end() { return pimpl->check_iterator_end(); }

void VczClass::getOneMarker(
    std::string& t_ref, std::string& t_alt, std::string& t_marker,
    uint32_t& t_pos, std::string& t_chr, double& t_altFreq, double& t_altCounts,
    double& t_missingRate, double& t_imputeInfo, bool t_isOutputIndexForMissing,
    std::vector<uint32_t>& t_indexForMissingforOneMarker,
    bool t_isOnlyOutputNonZero, std::vector<uint32_t>& t_indexForNonZero,
    bool& t_isBoolRead, std::vector<double>& dosages, bool t_isImputation) {
    pimpl->getOneMarker(
        t_ref, t_alt, t_marker, t_pos, t_chr, t_altFreq, t_altCounts,
        t_missingRate, t_imputeInfo, t_isOutputIndexForMissing,
        t_indexForMissingforOneMarker, t_isOnlyOutputNonZero, t_indexForNonZero,
        t_isBoolRead, dosages, t_isImputation);
}

}  // namespace VCZ
