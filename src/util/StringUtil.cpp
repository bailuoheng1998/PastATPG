// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3
#include "util/StringUtil.h"

namespace ictest {
std::vector<std::string> StringUtil::SplitBySpace(std::string str) {
  std::istringstream str_stream(str);
  std::vector<std::string> tokens;
  std::string t;
  while (str_stream >> t) {
    tokens.push_back(t);
  }
  return tokens;
}

//void fast_split(const std::string& source_str, std::vector<std::string>& target,
//                const std::string& split_flag = " ") {
//  target.clear();
//  auto start =
//      source_str.find_first_not_of(split_flag, 0);  // 分割到的字符串的第一个字符
//  auto position = source_str.find_first_of(split_flag, start);  // 分隔符的位置
//  while (position != std::string::npos || start != std::string::npos) {
//    // [start, position) 为分割下来的字符串
//    target.emplace_back(std::move(source_str.substr(start, position - start)));
//    start = source_str.find_first_not_of(split_flag, position);
//    position = source_str.find_first_of(split_flag, start);
//  }
//}

void StringUtil::split(const std::string_view& source_str,
                       std::vector<std::string_view>& target, char split_flag) {
  target.clear();
  auto max_idx = source_str.size();
  auto start_pos = source_str.find_first_not_of(split_flag, 0);
  auto end_pos = source_str.find_first_of(split_flag, start_pos);
  while (start_pos != std::string_view::npos) {
    end_pos = std::string_view::npos == end_pos ? (max_idx) : end_pos;
    target.emplace_back(&source_str[start_pos], end_pos - start_pos);
    start_pos = source_str.find_first_not_of(split_flag, end_pos);
    end_pos = source_str.find_first_of(split_flag, start_pos);
  }
}

void StringUtil::TrimWinCRLine(std::string& line) {
  line.erase(std::find(line.begin(), line.end(), '\r'), line.end());
}

void StringUtil::ClearAllSpace(std::string& line) {
  size_t idx = 0;
  if (!line.empty()) {
    while ((idx = line.find(' ', idx)) != std::string::npos) {
      line.erase(idx, 1);
    }
  }
}

}  // namespace ictest
