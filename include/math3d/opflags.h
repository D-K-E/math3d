#ifndef OPFLAGS_H
#define OPFLAGS_H

#include <cstdint>

namespace math3d {

enum opstatus_t : uint8_t {
  SUCCESS = 1,
  INDEX_ERROR = 2,
  SIZE_ERROR = 3,
  ARG_ERROR = 4,
  LU_ERROR = 5,
  NOT_IMPLEMENTED = 6
};

struct OpResult {
  //
  size_t line_info = 0;
  const char *file_name{nullptr};
  const char *fn_name{nullptr};
  const char *call_name{nullptr};
  const char *duration_info{nullptr};

  opstatus_t status = NOT_IMPLEMENTED;
  bool success = false;

  OpResult() = default;
  ~OpResult() = default;

  OpResult(size_t line, const char *fname,
           const char *funcname, const char *cname,
           opstatus_t op)
      : line_info(line), file_name(fname),
        fn_name(funcname), call_name(cname), status(op),
        success(op == SUCCESS) {}
};

} // namespace math3d

#endif
