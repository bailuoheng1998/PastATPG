// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3
#include "util/Log.h"

// bool Log::InitSpdlog(const std::string& logger_name,
//                      const std::string& logger_file_path) {
//   try {
//     auto file_logger = spdlog::rotating_logger_mt(logger_name,
//     logger_file_path,
//                                                   1024 * 1024 * 5,  // 5MB
//                                                   10);              // 10
//                                                   file
//     // spdlog::register_logger(file_logger);
//     spdlog::set_default_logger(file_logger);
//   } catch (const spdlog::spdlog_ex& ex) {
//     std::cerr << "Log initialization failed: " << ex.what() << std::endl;
//   }
//   return false;
// }

namespace ictest {

void Log::Init(const std::string& log_path) {
  // to stdout
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] %v%$");
  console_sink->set_level(spdlog::level::trace);
  std::vector<spdlog::sink_ptr> sinks{console_sink};
  // to file
  if (!log_path.empty() && log_path.size() > 0) {
    auto daily_file_sink =
        std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_path, true);
    daily_file_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] %v%$");
    daily_file_sink->set_level(spdlog::level::trace);
    sinks.emplace_back(daily_file_sink);
  }

  auto logger = std::make_shared<spdlog::logger>(LOG_DEFAULT_LOGGER_NAME,
                                                 sinks.begin(), sinks.end());
  logger->set_level(spdlog::level::trace);
  logger->flush_on(spdlog::level::trace);

  spdlog::register_logger(logger);
}

void Log::ShutDown() { spdlog::shutdown(); }
}  // namespace ictest