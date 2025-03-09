// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_LOG_H
#define ICTEST_LOG_H
#include <iostream>
#include <string>

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/rotating_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

namespace ictest {

#define LOG_DEFAULT_LOGGER_NAME "ictest_logger"
// #define LOG_CONFIG_RELEASE

#if defined(LOG_PLATFORM_WINDOWS)
#define LOG_BREAK __debugbreak();
#elif defined(LOG_PLATFORM_MAC)
#define LOG_BREAK __builtin_debugtrap();
#else
#define LOG_BREAK __builtin_trap();
#endif

#ifndef LOG_CONFIG_RELEASE
#define LOG_TRACE(...)                                        \
  if (spdlog::get(LOG_DEFAULT_LOGGER_NAME) != nullptr) {      \
    spdlog::get(LOG_DEFAULT_LOGGER_NAME)->trace(__VA_ARGS__); \
  }
#define LOG_DEBUG(...)                                        \
  if (spdlog::get(LOG_DEFAULT_LOGGER_NAME) != nullptr) {      \
    spdlog::get(LOG_DEFAULT_LOGGER_NAME)->debug(__VA_ARGS__); \
  }
#define LOG_INFO(...)                                        \
  if (spdlog::get(LOG_DEFAULT_LOGGER_NAME) != nullptr) {     \
    spdlog::get(LOG_DEFAULT_LOGGER_NAME)->info(__VA_ARGS__); \
  }
#define LOG_WARN(...)                                        \
  if (spdlog::get(LOG_DEFAULT_LOGGER_NAME) != nullptr) {     \
    spdlog::get(LOG_DEFAULT_LOGGER_NAME)->warn(__VA_ARGS__); \
  }
#define LOG_ERROR(...)                                        \
  if (spdlog::get(LOG_DEFAULT_LOGGER_NAME) != nullptr) {      \
    spdlog::get(LOG_DEFAULT_LOGGER_NAME)->error(__VA_ARGS__); \
  }
#define LOG_FATAL(...)                                           \
  if (spdlog::get(LOG_DEFAULT_LOGGER_NAME) != nullptr) {         \
    spdlog::get(LOG_DEFAULT_LOGGER_NAME)->critical(__VA_ARGS__); \
  }
#define LOG_ASSERT(x, msg)                                                \
  if ((x)) {                                                              \
  } else {                                                                \
    LOG_FATAL("ASSERT - {}\n\t{}\n\tin file: {}\n\ton line: {}", #x, msg, \
              __FILE__, __LINE__);                                        \
    LOG_BREAK                                                             \
  }
#else
// Disable logging for release builds
#define LOG_TRACE(...) (void)0
#define LOG_DEBUG(...) (void)0
#define LOG_INFO(...) (void)0
#define LOG_WARN(...) (void)0
#define LOG_ERROR(...) (void)0
#define LOG_FATAL(...) (void)0
#define LOG_ASSERT(x, msg) (void)0
#endif

class Log {
 public:
  //  static bool InitSpdlog(const std::string& logger_name,
  //                         const std::string& logger_file_path);
  Log() = default;
  ~Log() = default;
  void Init(const std::string& log_path = "");
  void ShutDown();
};

}  // namespace ictest

#endif  // ICTEST_LOG_H
