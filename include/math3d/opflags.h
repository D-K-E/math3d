#ifndef OPFLAGS_H
#define OPFLAGS_H

namespace math3d {

enum opstatus_t : std::uint_least8_t {
  SUCCESS = 1,
  INDEX_ERROR = 2,
  SIZE_ERROR = 3,
  ARG_ERROR = 4,
  LU_ERROR = 5,
  NOT_IMPLEMENTED = 6
};

struct OpResult {
  //
  unsigned int line_info;
  std::string file_name;
  std::string fn_name;
  std::string call_name;
  std::string duration_info;

  opstatus_t status;
  bool success = false;

  OpResult() = delete;
  ~OpResult() = default;

  OpResult(unsigned int line, const std::string &fname,
           const std::string &funcname,
           const std::string &cname, opstatus_t op)
      : line_info(line), file_name(fname),
        fn_name(funcname), call_name(cname), status(op),
        success(op == SUCCESS) {}
  OpResult(unsigned int line, const char *fname,
           const char *funcname, const char *cname,
           opstatus_t op)
      : line_info(line), file_name(fname),
        fn_name(funcname), call_name(cname), status(op),
        success(op == SUCCESS) {}
};

} // namespace math3d
#endif
