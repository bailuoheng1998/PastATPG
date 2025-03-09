// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_TIMEUTIL_H
#define ICTEST_TIMEUTIL_H
#include <chrono>
#include <iomanip>
#include <iostream>
#include "common/Define.h"
#include "util/Log.h"

#ifdef INSTRUMENT_TIMER
#define FUNC_TIME_REPORT_ICT Timer timer(__FUNCTION__);
#define FUNC_TIME_REPORT_WITH_FILE_ICT \
  Timer timer(                         \
      fmt::format("FILE :{} LINE: {} {} ", __FILE__, __LINE__, __FUNCTION__));
#else
#define FUNC_TIME_REPORT_ICT
#define FUNC_TIME_REPORT_WITH_FILE_ICT
#endif
namespace ictest {
class Timer {
 public:
  Timer() = delete;

  explicit Timer(const std::string_view FuncName) {
    function_name_ = FuncName;
    std::string res = fmt::format("-------{:->30} {} {:-<15}-------",
                                  get_current_time_and_date(), function_name_,
                                  "start running");
    std::cout << res << std::endl;
    start_ = std::chrono::steady_clock::now();
  }

  ~Timer() {
    std::chrono::duration<double> duration =
        std::chrono::steady_clock::now() - start_;
    std::string res = fmt::format("-------{:->30} cost {}s-----------",
                                  function_name_, duration.count());
    std::cout << res << std::endl;
  }

  static std::string get_current_time_and_date() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    return ss.str();
  }

 private:
  std::chrono::time_point<std::chrono::steady_clock> start_;
  std::string function_name_;
};
}  // namespace ictest
#endif  // ICTEST_TIMEUTIL_H
