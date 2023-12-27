#ifndef OPFLAGS_H
#define OPFLAGS_H

#include <cstdint>
#include <string>

namespace math3d {

enum opstatus_t : std::uint8_t {
  SUCCESS = 1,
  INDEX_ERROR = 2,
  SIZE_ERROR = 3,
  ARG_ERROR = 4,
  LU_ERROR = 5,
  NOT_IMPLEMENTED = 6
};

struct OpResult {
  //
  std::size_t line_info;
  std::string file_name;
  std::string fn_name;
  std::string call_name;
  std::string duration_info;

  opstatus_t status;
  bool success = false;

  OpResult() = delete;
  ~OpResult() = default;

  OpResult(std::size_t line, const std::string &fname,
           const std::string &funcname,
           const std::string &cname, opstatus_t op)
      : line_info(line), file_name(fname),
        fn_name(funcname), call_name(cname), status(op),
        success(op == SUCCESS) {}
  OpResult(std::size_t line, const char *fname,
           const char *funcname, const char *cname,
           opstatus_t op)
      : line_info(line), file_name(fname),
        fn_name(funcname), call_name(cname), status(op),
        success(op == SUCCESS) {}
};

bool CHECK(OpResult r);

OpResult INFO(const OpResult res);

OpResult INFO_VERBOSE(OpResult res);

} // namespace math3d
/**Checks if operation was successful*/
#define CHECK_MATH3D(call, res)                            \
  do {                                                     \
    res = call;                                            \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
  } while (0)

/**Prints information about the operation*/
#define INFO_MATH3D(call, res)                             \
  do {                                                     \
    res = call;                                            \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
    if (res.status != SUCCESS) {                           \
      res = math3d::INFO(res);                             \
    }                                                      \
  } while (0)

/**Prints everything about operation result*/
#define INFO_VERBOSE_MATH3D(call, res)                     \
  do {                                                     \
    auto start = std::chrono::steady_clock::now();         \
    res = call;                                            \
    auto stop = std::chrono::steady_clock::now();          \
    auto duration = std::chrono::duration_cast<            \
        std::chrono::microseconds>(stop - start);          \
    res.duration_info =                                    \
        to_string(static_cast<int>(duration.count()));     \
    res.call_name = #call;                                 \
    res.line_info = __LINE__;                              \
    res.file_name = __FILE__;                              \
    math3d::INFO_VERBOSE(res);                             \
  } while (0)

#endif
