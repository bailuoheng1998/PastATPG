// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_PATTERN_H
#define ICTEST_PATTERN_H
#include <string>
#include <vector>

#include "common/Define.h"
#include "common/Option.h"
#include "netlist/PrimNetlist.h"
#include "simulation/SimVal.h"
#include "util/Log.h"

namespace ictest {

enum class PatternCycleType {
  UNKNOWN_PCT,
  LOAD_UNLOAD_PCT,
  CAPTURE_PCT,
};

enum class PatternType {
  UNKNOWN_PT,
  COMB_PT,
  BASIC_SCAN_PT,
  FAST_SEQ_PT,
  FULL_SEQ_PT,
};

class PatternCycle {
 public:
  PatternCycle()
      : pctype_(PatternCycleType::UNKNOWN_PCT), last_parallel_bit_idx_(0) {}

  inline size_t GetSize() const { return input_.size(); }

  inline int32_t& LastParallelBitIdx() { return last_parallel_bit_idx_; }

  inline void SetPatternCycleType(PatternCycleType pctype) {
    this->pctype_ = pctype;
  }
  inline PatternCycleType GetPatternCycleType() { return pctype_; }

  // size: scan chain nums or frame nums
  inline void AllocInputMem(int32_t size) { input_.resize(size); }
  // size: scan chain nums or frame nums
  inline void AllocOutputMem(int32_t size) { output_.resize(size); }

  void AddInputVal(int32_t idx, const val64_t& logic_val) {
    input_[idx].emplace_back(logic_val);
  }
  void AddOutputVal(int32_t idx, const val64_t& logic_val) {
    output_[idx].emplace_back(logic_val);
  }

  void ResetInputValWithX() {
    for (int i = 0; i < input_.size(); i++) {
      for (int j = 0; j < input_[i].size(); j++) {
        input_[i][j].SetX();
      }
    }
  }

  void ResetOutputValWithX() {
    for (int i = 0; i < output_.size(); i++) {
      for (int j = 0; j < output_[i].size(); j++) {
        output_[i][j].SetX();
      }
    }
  }

  std::vector<std::vector<val64_t>>& GetInputVal() { return input_; }
  std::vector<std::vector<val64_t>>& GetOutputVal() { return output_; }

 private:
  int32_t last_parallel_bit_idx_{0};
  PatternCycleType pctype_{PatternCycleType::UNKNOWN_PCT};
  // multi scan chain load val64_t in scan chain order or pi val64_t
  std::vector<std::vector<val64_t>> input_;
  // multi scan chain unload val64_t in scan chain order or po val64_t
  std::vector<std::vector<val64_t>> output_;
};

class Pattern {
 public:
  void Clear() {
    load_unload_nums_ = 0;
    capture_nums_ = 0;
    pattern.clear();
  }

  void SetPatternCycleNums();

  inline int32_t GetLoadUnloadCycleNums() { return load_unload_nums_; }
  inline int32_t GetCaptureCycleNums() { return capture_nums_; }

  PatternCycle& GetPatternCycle(int32_t idx) {

    LOG_ASSERT(idx >= 0 && idx < pattern.size(),
               "ERROR: wrong pattern cycle idx");
    return pattern[idx];
  }

  void AddPattern(const PatternCycle& pcycle) { pattern.emplace_back(pcycle); }

  std::vector<PatternCycle>& GetPattern() { return pattern; }

  uint32_t GetPatternId() { return pattern_id_; }

  void SetPatternId(uint32_t pattern_id) { pattern_id_ = pattern_id; }

  /**
   * 保存当前向量的格式，但是将当前向量中的每个cycle的input和output的值都恢复为X。
   */
  void KeepPatternFormatAndInitialWithX() {
    for (auto& patternCycle : pattern) {
      patternCycle.ResetInputValWithX();
      patternCycle.ResetOutputValWithX();
      patternCycle.LastParallelBitIdx() = 0;
    }
  }

  /**
   * 复制template_pattern的格式到当前向量，并使用X初始化当前向量。
   * @param template_pattern 模板向量
   * @attention 该函数只复制template_pattern向量的格式，并不复制其里面的值。
   */
  void CopyPatternFormatAndInitialWithX(const Pattern& template_pattern) {
    pattern.clear();
    for (PatternCycle patternCycle : template_pattern.pattern) {
      patternCycle.ResetInputValWithX();
      patternCycle.ResetOutputValWithX();
      pattern.push_back(patternCycle);
    }
    SetPatternCycleNums();
  }

  bool IsSamePatternFormat(const Pattern& template_pattern) {
    return template_pattern.load_unload_nums_ == load_unload_nums_
    && template_pattern.capture_nums_ == capture_nums_;
  }

  std::string FormatToString(){
    std::string format = "Pattern " + std::to_string(pattern_id_) + " Format is: ";
    for (auto& pac : pattern) {
      if (pac.GetPatternCycleType() == PatternCycleType::LOAD_UNLOAD_PCT) {
        format += "loun_unlod\t";
      } else if (pac.GetPatternCycleType() == PatternCycleType::CAPTURE_PCT) {
        format += "capture\t";
      } else {
        format += "Unknown\t";
      }
    }
    return format + "\n";
  }

 private:
  int32_t load_unload_nums_ = 0;
  int32_t capture_nums_ = 0;
  uint32_t pattern_id_ = 0;
  std::vector<PatternCycle> pattern;
};

class PatternParser {
 public:
  inline void SetupPrimNetlist(PrimNetlist* pnlist) {
    LOG_ASSERT(pnlist, "ERROR: prim netlist is nullptr");
    pnlist_ = pnlist;
  }

  std::vector<Pattern>& GetSerialPatternSet() { return serial_pattern_sets_; }
  std::vector<Pattern>& GetParallelPatternSet() {
    return parallel_pattern_sets_;
  }

  // set pattern type from different circuit type
  inline void SetPatternType(PatternType ptype) { ptype_ = ptype; }
  inline PatternType GetPatternType() { return ptype_; }

  // get internal com patterns from ATPG
  void LoadInternalComPatternForSAF(
      std::vector<std::vector<int>>& serial_pi_pats,
      std::vector<Pattern>& parallel_pi_pats, bool set_x_2_rand01 = false);



  void LoadInternalFullScanPatternForSAF(
      std::vector<std::vector<std::vector<int>>>& serial_pi_pats,
      std::vector<std::vector<std::vector<int>>>& serial_sc_pats,
      std::vector<Pattern>& parallel_patterns, bool set_x_2_rand01 = false);


  // get external patterns from .ictest.pat file
  bool ReadPatternFile(const std::string& pattern_file);
  void GetLoadUnloadStringLine(const std::vector<std::string>& lines,
                               PatternCycle& pc);
  void GetCaptureStringLine(const std::vector<std::string>& lines,
                            PatternCycle& pc);



  //  transform serial pattern to parallel pattern
  void TransformSerialToParallelForSAF(std::vector<Pattern>& serial_patterns,
                                       std::vector<Pattern>& parallel_patterns,
                                       bool is_external);


  // fill pattern with X val
  void InitPatternWithX(Pattern& sp);

  void InitPatternWithX(Pattern& sp, int capture_size, int load_unload_size);
  void InitPatternCycleWithX(PatternCycle& patternCycle,
                             int frameSize_or_scanChainSize,
                             int piSize_or_siSize, int poSize_or_soSize);

  // fill 64 serial patterns into a packet
  void FillPatternPacketForSAF(int pattern_id, Pattern& sp,
                               std::vector<Pattern>& serial_patterns,
                               bool is_external);



  val64_t GetOneBit(const val64_t& source, int32_t bit_index);

  void TransformParallelToSerial(std::vector<Pattern>& parallel_patterns,
                                 std::vector<Pattern>& serial_patterns);

  // report输出到终端
  void ReportPattern();
  // dump向量到文件
  void DumpPatternFile(std::vector<Pattern>& serial_pattern_set,
                       const std::string& pattern_file);

 private:
  PrimNetlist* pnlist_{nullptr};
  PatternType ptype_{PatternType::UNKNOWN_PT};
  std::vector<Pattern> serial_pattern_sets_;
  // only support com and full scan circuit
  std::vector<Pattern> parallel_pattern_sets_;
};
}  // namespace ictest
#endif  // ICTEST_PATTERN_H
