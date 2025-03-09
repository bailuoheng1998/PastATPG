// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_STRINGUTIL_H
#define ICTEST_STRINGUTIL_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// #include "boost/algorithm/string.hpp"
namespace ictest {
class StringUtil {
 public:
  //  static void Split(std::vector<std::string>& res, const std::string& src,
  //  const std::string& sep) {
  //    boost::split(res, src, boost::is_any_of(sep));
  //  }
  //
  //  static std::string Join(std::vector<std::string>& src, const std::string&
  //  sep) {
  //    return boost::join(src, sep);
  //  }
  //
  //  static std::string ToUpperCopy(const std::string& str) {
  //    return boost::to_upper_copy(str);
  //  }
  //
  //  static std::string ToLowerCopy(const std::string& str) {
  //    return boost::to_lower_copy(str);
  //  }
  //
  //  static bool StartsWith(const std::string& a, const std::string& b) {
  //    return boost::starts_with(a, b);
  //  }
  //
  //  static bool IgnoreStartsWith(const std::string& a, const std::string& b) {
  //    return boost::istarts_with(a, b);
  //  }
  //
  //  static bool EndsWith(const std::string& a, const std::string& b) {
  //    return boost::ends_with(a, b);
  //  }
  //
  //  static bool IgnoreEndsWith(const std::string& a, const std::string& b) {
  //    return boost::iends_with(a, b);
  //  }
  //
  //  static bool Contains(const std::string& a, const std::string& b) {
  //    return boost::contains(a, b);
  //  }
  //
  //  static bool IgnoreContains(const std::string& a, const std::string& b) {
  //    return boost::icontains(a, b);
  //  }
  //
  //  static bool Equals(const std::string& a, const std::string& b) {
  //    return boost::equals(a, b);
  //  }
  //
  //  static bool IgnoreEquals(const std::string& a, const std::string& b) {
  //    return boost::iequals(a, b);
  //  }
  //
  //  static bool IsAlNum(char c) {
  //    return boost::is_alnum()(c);
  //  }
  //
  //  static bool IsAlpha(char c) {
  //    return boost::is_alpha()(c);
  //  }
  //
  //  static bool IsDigit(char c) {
  //    return boost::is_digit()(c);
  //  }
  //
  //  static bool IsUpper(char c) {
  //    return boost::is_upper()(c);
  //  }
  //
  //  static bool IsLower(char c) {
  //    return boost::is_lower()(c);
  //  }
  //
  //  static bool IsAnyOf(const std::string& str,  char c) {
  //    return boost::is_any_of(str)(c);
  //  }
  //
  //  static void Trim(std::string& str) {
  //    boost::trim(str);
  //  }
  //
  //  static std::string TrimCopy(const std::string& str) {
  //    return boost::trim_copy(str);
  //  }
  //
  //  static void TrimLeft(std::string& str) {
  //    boost::trim_left(str);
  //  }
  //
  //  static std::string TrimLeftCopy(const std::string& str) {
  //    return boost::trim_left_copy(str);
  //  }
  //
  //  static void TrimRight(std::string& str) {
  //    boost::trim_right(str);
  //  }
  //
  //  static std::string TrimRightCopy(const std::string& str) {
  //    return boost::trim_right_copy(str);
  //  }

  static std::vector<std::string> SplitBySpace(std::string str);
  /**
   * 处理win \n\r和linux \name换行符不同的问题
   * @param line
   */
  static void TrimWinCRLine(std::string& line);

  /**
   * @brief 清空字符串中所有的空格
   * @param line
   */
  static void ClearAllSpace(std::string& line);

  static void split(const std::string_view& source_str,
                  std::vector<std::string_view>& target,
                  char split_flag) ;
};
}  // namespace ictest

#endif  // ICTEST_STRINGUTIL_H
