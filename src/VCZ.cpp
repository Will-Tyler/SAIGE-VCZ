#include "VCZ.hpp"

#include <cassert>
#include <memory>
#include <string>
#include <utility>
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

template <typename ElementType = void>
class ChunkCacher {
   public:
    ChunkCacher() = default;

    ChunkCacher(tensorstore::Spec input_spec) {
        store_ = tensorstore::Open<ElementType, tensorstore::dynamic_rank,
                                   tensorstore::ReadWriteMode::read>(input_spec)
                     .value();
        cache_first_chunk();
    }

    ChunkCacher(
        tensorstore::TensorStore<ElementType, tensorstore::dynamic_rank,
                                 tensorstore::ReadWriteMode::read>&& store) {
        store_ = store;
        cache_first_chunk();
    }

    template <size_t N>
    ElementType& operator()(const Index (&indices)[N]) {
        if (tensorstore::Contains(array_.domain(), indices))
            return array_(indices);

        auto rank = store_.rank();
        const auto shape = store_.domain().shape();
        auto chunk_shape = store_.chunk_layout().value().read_chunk().shape();
        std::vector<Index> chunk_starts(rank);
        std::vector<Index> chunk_ends(rank);

        for (Index dim_index = 0; dim_index < rank; dim_index++) {
            const Index chunk_index = indices[dim_index] / chunk_shape[dim_index];
            chunk_starts[dim_index] = chunk_index * chunk_shape[dim_index];
            chunk_ends[dim_index] = chunk_starts[dim_index] + chunk_shape[dim_index];
            chunk_ends[dim_index] = std::min(chunk_ends[dim_index], shape[dim_index]);
        }

        array_ = tensorstore::Read(store_ | tensorstore::AllDims().HalfOpenInterval(chunk_starts, chunk_ends)).value();
        return array_(indices);
    }

   private:
    tensorstore::TensorStore<ElementType, tensorstore::dynamic_rank,
                             tensorstore::ReadWriteMode::read>
        store_;
    tensorstore::SharedArray<ElementType, tensorstore::dynamic_rank,
                             tensorstore::offset_origin>
        array_;

    void cache_first_chunk() {
        const auto rank = store_.rank();
        const auto shape = store_.domain().shape();
        const auto chunk_shape = store_.chunk_layout().value().read_chunk().shape();
        const std::vector<Index> chunk_starts(rank);
        std::vector<Index> chunk_ends(rank);

        for (Index dim_index = 0; dim_index < rank; dim_index++)
            chunk_ends[dim_index] = std::min(shape[dim_index], chunk_shape[dim_index]);

        array_ = tensorstore::Read(store_ | tensorstore::AllDims().HalfOpenInterval(chunk_starts, chunk_ends)).value();
    }
};

class VczClass::Impl {
   public:
    Impl(std::string& t_vczFileName,
         std::vector<std::string>& t_sampleInModel) {
        m_model_sample_count = t_sampleInModel.size();
        tensorstore::Spec input_spec;
        input_spec = tensorstore::Spec::FromJson(
                         {
                             {"driver", "zarr"},
                             {"kvstore",
                              {{"driver", "file"},
                               {"path", t_vczFileName + "/call_genotype"}}},
                         })
                         .value();
        auto call_genotype_store =
            tensorstore::Open<int8_t, tensorstore::dynamic_rank,
                              tensorstore::ReadWriteMode::read>(
                input_spec, tensorstore::OpenMode::open,
                tensorstore::ReadWriteMode::read)
                .value();
        m_genotype_shape = call_genotype_store.domain().shape();
        m_call_genotype_array =
            ChunkCacher<int8_t>(std::move(call_genotype_store));

        input_spec = tensorstore::Spec::FromJson(
                         {
                             {"driver", "zarr"},
                             {"kvstore",
                              {{"driver", "file"},
                               {"path", t_vczFileName + "/contig_length"}}},
                         })
                         .value();
        m_contig_length_array = ChunkCacher<int64_t>(input_spec);

        input_spec = tensorstore::Spec::FromJson(
                         {
                             {"driver", "zarr"},
                             {"kvstore",
                              {{"driver", "file"},
                               {"path", t_vczFileName + "/variant_contig"}}},
                         })
                         .value();
        m_variant_contig_array = ChunkCacher<int8_t>(input_spec);

        input_spec = tensorstore::Spec::FromJson(
                         {
                             {"driver", "zarr"},
                             {"kvstore",
                              {{"driver", "file"},
                               {"path", t_vczFileName + "/variant_position"}}},
                         })
                         .value();
        m_variant_position_array = ChunkCacher<int8_t>(input_spec);

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
        const Index sample_count =
            std::min(m_genotype_shape[1], (Index)m_model_sample_count);
        dosages.clear();
        dosages.resize(sample_count);

        for (Index sample_index = 0; sample_index < sample_count;
             sample_index++) {
            const int8_t a =
                m_call_genotype_array({m_marker_index, sample_index, 0});
            const int8_t b =
                m_call_genotype_array({m_marker_index, sample_index, 1});

            if (a >= 0 && b >= 0) {
                const double dosage = a + b;
                dosages[sample_index] = dosage;
                t_altCounts += dosage;

                if (dosage > 0) t_indexForNonZero.push_back(sample_index);
            } else {
                dosages[sample_index] = -1;
                t_indexForMissingforOneMarker.push_back(sample_index);
                missing_count += 1;
            }
        }

        t_missingRate = (double)missing_count / sample_count;

        if (missing_count == sample_count) {
            t_altFreq = 0.0;
        } else {
            t_altFreq = t_altCounts / 2 / (sample_count - missing_count);
        }

        t_isBoolRead = true;
    }

   private:
    ChunkCacher<int64_t> m_contig_length_array;
    ChunkCacher<int8_t> m_variant_contig_array;
    ChunkCacher<int8_t> m_variant_position_array;
    ChunkCacher<int8_t> m_call_genotype_array;

    tensorstore::span<const Index> m_genotype_shape;
    std::size_t m_model_sample_count;
    std::string m_chrom;
    Index m_marker_index;
    int m_startPos;
    int m_endPos;
};

VczClass::VczClass(std::string t_vczFileName,
                   std::vector<std::string> t_SampleInModel) {
    pimpl = std::make_unique<Impl>(t_vczFileName, t_SampleInModel);
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
