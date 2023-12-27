#include <math3d/opflags.h>
#include <iostream>

namespace math3d {
//

bool CHECK(OpResult r) { return r.status == SUCCESS; }

OpResult INFO(const OpResult res) {
  std::cerr << res.status << " :: " << res.file_name
            << " :: " << res.line_info
            << " :: " << res.fn_name
            << " :: " << res.call_name << std::endl;
  return res;
}
OpResult INFO_VERBOSE(OpResult res) {
  if (res.status == SUCCESS) {
    std::cerr << "SUCCESS "
              << " :: " << res.file_name
              << " :: " << res.fn_name
              << " :: " << res.call_name
              << " :: " << res.duration_info
              << " microseconds" << std::endl;
    return res;
  }
  return INFO(res);
}

} // namespace math3d
