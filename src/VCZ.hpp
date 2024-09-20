#ifndef VCZ_HPP
#define VCZ_HPP

#include <string>
#include <vector>

namespace VCZ {

class VczClass {
   private:
    class Impl;
    std::unique_ptr<Impl> pimpl;

   public:
    VczClass(std::string t_vczFileName,
             std::vector<std::string> t_SampleInModel);
    void set_iterator(std::string& chrom, const int start_pos,
                      const int end_pos);
    void move_forward_iterator(int a);
    bool check_iterator_end();
    void getOneMarker(std::string& t_ref, std::string& t_alt,
                      std::string& t_marker, uint32_t& t_pos,
                      std::string& t_chr, double& t_altFreq,
                      double& t_altCounts, double& t_missingRate,
                      double& t_imputeInfo, bool t_isOutputIndexForMissing,
                      std::vector<uint32_t>& t_indexForMissingforOneMarker,
                      bool t_isOnlyOutputNonZero,
                      std::vector<uint32_t>& t_indexForNonZero,
                      bool& t_isBoolRead, std::vector<double>& dosages,
                      bool t_isImputation);
};

}  // namespace VCZ

#endif
