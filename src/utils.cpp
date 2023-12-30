#include <iostream>
#include <math3d/utils.h>
#include <sstream>

namespace math3d {

bool CHECK(OpResult r) { return r.status == SUCCESS; }

std::stringstream to_xml(const OpResult &res) {
  std::stringstream s;
  s << "<OpResult status='" << res.status << "' function='"
    << res.fn_name << "' file='" << res.file_name
    << "' callsite='" << res.call_name << "' line='"
    << res.line_info << "' duration='" << res.duration_info
    << "'"
    << " durationUnit='microsecond'/>";
  return s;
}
OpResult INFO(const OpResult &res) {
  std::stringstream s = to_xml(res);
  std::cerr << s.str() << std::endl;
  return res;
}
OpResult INFO_VERBOSE(const OpResult &res) {
  return INFO(res);
}

} // namespace math3d

//
std::ostream &operator<<(std::ostream &out,
                         math3d::opstatus_t flag) {
  switch (flag) {
  case math3d::SUCCESS: {
    out << "SUCCESS";
    break;
  }
  case math3d::SIZE_ERROR: {
    out << "SIZE_ERROR";
    break;
  }
  case math3d::INDEX_ERROR: {
    out << "INDEX_ERROR";
    break;
  }
  case math3d::ARG_ERROR: {
    out << "ARG_ERROR";
    break;
  }
  case math3d::NOT_IMPLEMENTED: {
    out << "NOT_IMPLEMENTED";
    break;
  }
  case math3d::LU_ERROR: {
    out << "LU_DECOMPOSITION_ERROR";
    break;
  }
  }
  return out;
}
