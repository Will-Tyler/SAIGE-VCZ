#ifndef VCZ_HPP
#define VCZ_HPP

#include <string>
#include <vector>

namespace VCZ {

class VczClass {
   public:
    VczClass(std::string t_vczFileName, std::vector<std::string> t_SampleInModel);
    bool setVczObj(const std::string &t_vczFileName);
};

}  // namespace VCZ

#endif
