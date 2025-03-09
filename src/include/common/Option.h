// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_SIMOPTION_H
#define ICTEST_SIMOPTION_H
#include <cstdint>

namespace ictest {

struct Option {
  /// 设置故障控制选项
  /// 故障来源 ICTest or TMax(Tessent): 1. internal(default) 2. external
  std::string flist_source_ = "internal";
  /// 指定故障格式: 1. tmax(default) 2. tessent
  std::string fault_format_ = "tmax";
  /// 如果flist_source是external，指定故障输入文件，默认为空
  std::string external_flist_path_ = "";
  /// 指定是否采用故障压缩技术
  bool is_eqv_flist_ = false;
  /// 指定故障模型: 1. stuck-at(default) 2. transition
  std::string fault_model_ = "stuck-at";

  /// 设置向量控制选项
  /// 向量来源 atpg or stil file: 1.internal(default) 2. external
  std::string pattern_source_ = "internal";
  /// 如果向量来自于external，指定向量文件路径，默认为空
  std::string external_pattern_path_ = "";

  /// 设置atpg控制选项
  /// 选择ATPG引擎，要和对应的电路类型，向量类型保持一致
  /// 类型包括: 1.comb(default) 2.basic_scan 3.fast_seq 4.full_seq
  std::string atpg_engine_type_ = "comb";
  /// 设置atpg abort次数，默认为10
  int32_t abort_limit_ = 10;
  /// 设置ATPG每次仿真向量数，默认为64
  int32_t pattern_per_sim_ = 64;


  /// 设置仿真控制选项
  /// 选择仿真引擎，要和对应的电路类型，向量类型保持一致
  /// 类型包括: 1.comb(default) 2.basic_scan 3.fast_seq 4.full_seq
  std::string sim_engine_type_ = "comb";
  /// 是否比较PO仿真直，默认为false
  bool compare_po_value_ = false;
  /// 是否比较PPO仿真值，默认为false
  /// 1. serial pattern 2. parallel pattern(default)
  std::string pattern_bit_mode_ = "parallel";

  bool dump_good_sim_res = false;

  // measure po/ppo for gsim
  std::string out_file_path;
};
}  // namespace ictest
#endif  // ICTEST_SIMOPTION_H
