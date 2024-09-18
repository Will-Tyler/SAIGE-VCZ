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
    void set_iterator(std::string& chrom, const int start_pos, const int end_pos);
    void move_forward_iterator(int a);
    bool check_iterator_end();
};

}  // namespace VCZ

#endif
