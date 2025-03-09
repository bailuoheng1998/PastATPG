// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "pattern/Pattern.h"

#include <fstream>
#include <sstream>

#include "atpg/ATPGDefine.h"
#include "util/StringUtil.h"
#include "util/TimeUtil.h"
namespace ictest {

void Pattern::SetPatternCycleNums() {
  load_unload_nums_ = 0;
  capture_nums_ = 0;
  for (int i = 0; i < pattern.size(); i++) {
    auto pcycle = pattern[i];
    if (pcycle.GetPatternCycleType() == PatternCycleType::LOAD_UNLOAD_PCT) {
      load_unload_nums_++;
    } else if (pcycle.GetPatternCycleType() == PatternCycleType::CAPTURE_PCT) {
      capture_nums_++;
    } else {
      LOG_ASSERT(false, "ERROR:unknown pattern cycle type");
    }
  }
}

void PatternParser::LoadInternalComPatternForSAF(
    std::vector<std::vector<int>>& serial_pi_pats,
    std::vector<Pattern>& parallel_pi_pats, bool set_x_2_rand01) {
  val64_t value;
  std::vector<Pattern> serial_patterns;
  for (int pat_id = 0; pat_id < serial_pi_pats.size(); pat_id++) {
    Pattern sp;
    PatternCycle pc;
    pc.SetPatternCycleType(PatternCycleType::CAPTURE_PCT);
    pc.AllocInputMem(1);
    pc.AllocOutputMem(1);
    for (int pi_id = 0; pi_id < pnlist_->NumPIs(); pi_id++) {
      if (serial_pi_pats[pat_id][pi_id] == LOGIC_0 ||
          serial_pi_pats[pat_id][pi_id] == LOGIC_S0) {
        value.Set0();
      } else if (serial_pi_pats[pat_id][pi_id] == LOGIC_1 ||
                 serial_pi_pats[pat_id][pi_id] == LOGIC_S1) {
        value.Set1();
      } else if (serial_pi_pats[pat_id][pi_id] == LOGIC_x) {
        if (set_x_2_rand01) {
          if (rand() % 2 == 0) {
            value.Set0();
            serial_pi_pats[pat_id][pi_id] = LOGIC_0;
          } else {
            value.Set1();
            serial_pi_pats[pat_id][pi_id] = LOGIC_1;
          }
        } else {
          value.SetX();
        }
      } else {
        LOG_ASSERT(false, "ERROR: not support logic value");
      }
      pc.AddInputVal(0, value);
    }
    sp.AddPattern(pc);
    sp.SetPatternCycleNums();
    serial_patterns.emplace_back(sp);
  }
  TransformSerialToParallelForSAF(serial_patterns, parallel_pi_pats, false);
}


void PatternParser::LoadInternalFullScanPatternForSAF(
    std::vector<std::vector<std::vector<int>>>& serial_pi_pats,
    std::vector<std::vector<std::vector<int>>>& serial_sc_pats,
    std::vector<Pattern>& parallel_patterns, bool set_x_2_rand01) {
  val64_t value;
  std::vector<Pattern> serial_patterns;
  int pat_nums = serial_pi_pats.size();

  LOG_ASSERT(serial_sc_pats.size() == pat_nums,
             "ERROR: pi pattern nums not equal scan pattern nums");
  for (int pat_id = 0; pat_id < pat_nums; pat_id++) {
    Pattern sp;
    // scan chain load_unload pattern cycle
    PatternCycle sc_pc;
    sc_pc.SetPatternCycleType(PatternCycleType::LOAD_UNLOAD_PCT);
    sc_pc.AllocInputMem(pnlist_->ScanChains().size());
    sc_pc.AllocOutputMem(pnlist_->ScanChains().size());

    // scan chain nums
    for (int i = 0; i < pnlist_->ScanChains().size(); i++) {
      // scan cell size in single chain
      for (int j = 0; j < pnlist_->ScanChains()[i]->GetChainLength(); j++) {
        if (serial_sc_pats[pat_id][i][j] == LOGIC_0 ||
            serial_sc_pats[pat_id][i][j] == LOGIC_S0) {
          value.Set0();
        } else if (serial_sc_pats[pat_id][i][j] == LOGIC_1 ||
                   serial_sc_pats[pat_id][i][j] == LOGIC_S1) {
          value.Set1();
        } else if (serial_sc_pats[pat_id][i][j] == LOGIC_x) {
          if (set_x_2_rand01) {
            if (rand() % 2 == 0) {
              value.Set0();
              serial_sc_pats[pat_id][i][j] = LOGIC_0;
            } else {
              value.Set1();
              serial_sc_pats[pat_id][i][j] = LOGIC_1;
            }
          } else {
            value.SetX();
          }
        } else {
          LOG_ASSERT(false, "ERROR: not support scan chain values");
        }
        sc_pc.AddInputVal(i, value);
      }
    }
    sp.AddPattern(sc_pc);
    // pi capture pattern cycle
    PatternCycle pi_pc;
    pi_pc.SetPatternCycleType(PatternCycleType::CAPTURE_PCT);
    pi_pc.AllocInputMem(MAX_FRAME_NUM);
    pi_pc.AllocOutputMem(MAX_FRAME_NUM);
    for (int i = 0; i < pnlist_->NumPIs(); i++) {
      for (int j = 0; j < MAX_FRAME_NUM; j++) {
        if (serial_pi_pats[pat_id][i][j] == LOGIC_0 ||
            serial_pi_pats[pat_id][i][j] == LOGIC_S0) {
          value.Set0();
        } else if (serial_pi_pats[pat_id][i][j] == LOGIC_1 ||
                   serial_pi_pats[pat_id][i][j] == LOGIC_S1) {
          value.Set1();
        } else if (serial_pi_pats[pat_id][i][j] == LOGIC_x) {
          if (j == FRAME_0) {
            if (set_x_2_rand01) {
              if (rand() % 2 == 0) {
                value.Set0();
                serial_pi_pats[pat_id][i][j] = LOGIC_0;
              } else {
                value.Set1();
                serial_pi_pats[pat_id][i][j] = LOGIC_1;
              }
            } else {
              value.SetX();
            }
          } else {
            serial_pi_pats[pat_id][i][j] = serial_pi_pats[pat_id][i][j - 1];
          }
        } else {
          LOG_ASSERT(false, "ERROR: not support pi values");
        }
        pi_pc.AddInputVal(j, value);
      }
    }
    sp.AddPattern(pi_pc);
    sp.SetPatternCycleNums();
    serial_patterns.emplace_back(sp);
    serial_pattern_sets_.emplace_back(sp);
  }
  TransformSerialToParallelForSAF(serial_patterns, parallel_patterns, false);
}


bool PatternParser::ReadPatternFile(const std::string& pattern_file) {
  std::ifstream fin(pattern_file);
  if (!fin.is_open()) {
    LOG_ERROR("ERROR:failed open pattern file : " + pattern_file);
    return false;
  }
  FUNC_TIME_REPORT_WITH_FILE_ICT
  std::string line;
  Pattern sp;
  std::vector<std::string> load_unload;
  std::vector<std::string> capture;
  int32_t idx = 0;
  while (getline(fin, line)) {
    StringUtil::TrimWinCRLine(line);
    if (line.empty()) {
      continue;
    } else if ("[BEGIN]" == line) {
      continue;
    } else if ("[PATTERN_START]" == line) {
      continue;
    } else if ("[LOAD_UNLOAD_START]" == line) {
      while (getline(fin, line)) {
        StringUtil::TrimWinCRLine(line);
        if ("[LOAD_UNLOAD_END]" == line) {
          PatternCycle pc;
          GetLoadUnloadStringLine(load_unload, pc);
          sp.AddPattern(pc);
          load_unload.clear();
          break;
        } else {
          if (!line.empty()) {
            load_unload.emplace_back(line);
          }
        }
      }
    } else if ("[CAPTURE_START]" == line) {
      while (getline(fin, line)) {
        StringUtil::TrimWinCRLine(line);
        if ("[CAPTURE_END]" == line) {
          PatternCycle pc;
          this->GetCaptureStringLine(capture, pc);
          sp.AddPattern(pc);
          capture.clear();
          break;
        } else {
          if (!line.empty()) {
            capture.emplace_back(line);
          }
        }
      }
    } else if ("[PATTERN_END]" == line) {
      sp.SetPatternCycleNums();
      sp.SetPatternId(idx++);
      serial_pattern_sets_.emplace_back(sp);
      sp.Clear();
    } else if ("[END]" == line) {
      break;
    }
  }
  return true;
}

void PatternParser::ReportPattern() {
  LOG_INFO("{:^55}", "## Report Pattern ##");
  LOG_INFO("{:^20}:     {:<}", "Pattern Nums", serial_pattern_sets_.size());
}

void PatternParser::DumpPatternFile(std::vector<Pattern>& serial_pattern_set,
                                    const std::string& pattern_file) {
  std::ofstream file(pattern_file);
  file << "[BEGIN]";
  file << "\n";
  for (int i = 0; i < serial_pattern_sets_.size(); i++) {
    file << "[PATTERN_START]";
    file << "\n";
    auto& single_pattern = serial_pattern_sets_[i];
    auto& pattern_cycles = single_pattern.GetPattern();
    for (int j = 0; j < pattern_cycles.size(); j++) {
      auto& pc = pattern_cycles[j];
      if (pc.GetPatternCycleType() == PatternCycleType::LOAD_UNLOAD_PCT) {
        file << "[LOAD_UNLOAD_START]";
        file << "\n";
        auto& inputs = pc.GetInputVal();
        auto& outputs = pc.GetOutputVal();
        std::stringstream ss;
        for (int k = 0; k < inputs.size(); k++) {
          ss << "test_si ";
          auto& si_vals = inputs[k];
          for (auto& val : si_vals) {
            ss << val.GetSerialValue();
          }
          ss << " test_so ";
          auto& so_vals = outputs[k];
          for (auto& val : so_vals) {
            char c = val.GetSerialValue()[0];
            switch (c) {
              case '0':
                ss << "L";
                break;
              case '1':
                ss << "H";
                break;
              case 'X':
                ss << "X";
                break;
              case 'Z':
                ss << "Z";
                break;
              default:
                break;
            }
          }
          file << ss.str();
          file << "\n";
          ss.str("");
        }
        file << "[LOAD_UNLOAD_END]\n";
        file << "\n";
      } else if (pc.GetPatternCycleType() == PatternCycleType::CAPTURE_PCT) {
        file << "[CAPTURE_START]";
        file << "\n";
        auto& inputs = pc.GetInputVal();
        auto& outputs = pc.GetOutputVal();
        std::stringstream ss;
        LOG_ASSERT(inputs.size() == outputs.size(),
                   "ERROR: capture pattern input size not equal output size");
        for (int k = 0; k < inputs.size(); k++) {
          ss << "_pi ";
          auto& pi_vals = inputs[k];
          for (auto& val : pi_vals) {
            ss << val.GetSerialValue();
          }

          auto& po_vals = outputs[k];
          if (!po_vals.empty()) {
            ss << " _po ";
          }
          for (auto& val : po_vals) {
            char c = val.GetSerialValue()[0];
            switch (c) {
              case '0':
                ss << "L";
                break;
              case '1':
                ss << "H";
                break;
              case 'X':
                ss << "X";
                break;
              case 'Z':
                ss << "Z";
                break;
              default:
                break;
            }
          }
          file << ss.str();
          file << "\n";
          ss.str("");
        }
        file << "[CAPTURE_END]\n";
      }
    }
    file << "[PATTERN_END]\n";
    file << "\n";
  }
  file << "[END]";
}


void PatternParser::GetLoadUnloadStringLine(
    const std::vector<std::string>& lines, PatternCycle& pc) {
  pc.AllocInputMem(lines.size());
  pc.AllocOutputMem(lines.size());
  pc.SetPatternCycleType(PatternCycleType::LOAD_UNLOAD_PCT);
//    std::vector<std::string_view> load_unload_view;

  for (int i = 0; i < lines.size(); i++) {
    std::istringstream iss(lines[i]);
    std::string tmp;
    std::vector<std::string> load_unload;
    while (iss >> tmp) {
      load_unload.emplace_back(tmp);
    }
//    StringUtil::split(lines[i],load_unload_view,' ');
    // test_si si_val test_so so_val
    //    LOG_ASSERT(load_unload.size() == 4, "ERROR: wrong load unload line
    //    size");
    // si_val
    for (auto elem : load_unload[1]) {
      val64_t val;
      if (elem == '0') {
        val.Set0();
      } else if (elem == '1') {
        val.Set1();
      } else if (elem == 'X') {
        val.SetX();
      } else {
        LOG_ASSERT(false, "ERROR:not support stil load vals:" + elem);
      }
      pc.AddInputVal(i, val);
    }
    // so_val
    for (auto elem : load_unload[3]) {
      val64_t val;
      if (elem == 'L' || elem == '0') {
        val.Set0();
      } else if (elem == 'H' || elem == '1') {
        val.Set1();
      } else if (elem == 'X') {
        val.SetX();
      } else if (elem == 'T') {
        val.SetX();
      } else {
        LOG_ASSERT(false, "ERROR:not support stil unload vals:" + elem);
      }
      pc.AddOutputVal(i, val);
    }
    iss.str("");
  }
}

void PatternParser::GetCaptureStringLine(const std::vector<std::string>& lines,
                                         PatternCycle& pc) {
  pc.AllocInputMem(lines.size());
  pc.AllocOutputMem(lines.size());
  pc.SetPatternCycleType(PatternCycleType::CAPTURE_PCT);
  for (int i = 0; i < lines.size(); i++) {
    std::istringstream iss(lines[i]);
    std::string tmp;
    std::vector<std::string> capture;
    while (iss >> tmp) {
      capture.emplace_back(tmp);
    }
    // pi pi_val po po_val
//    LOG_ASSERT(capture.size() >= 2, "ERROR: wrong capture line size");
    // pi_val
    for (int idx = 0; idx < pnlist_->NumPIs(); idx++) {
      val64_t val;
      char elem = capture[1][idx];
      if (elem == '0') {
        val.Set0();
      } else if (elem == '1') {
        val.Set1();
      } else if (elem == 'X') {
        val.SetX();
      } else if (elem == 'Z') {
        val.SetZ();
      } else {
        LOG_ASSERT(false,fmt::format("ERROR:not support stil pi vals:{}" , elem));
      }
      pc.AddInputVal(i, val);
    }
    // po_val
    if (capture.size() == 4) {
      for (auto elem : capture[3]) {
        val64_t val;
        if (elem == 'L' || elem == '0') {
          val.Set0();
        } else if (elem == 'H' || elem == '1') {
          val.Set1();
        } else if (elem == 'X') {
          val.SetX();
        } else if (elem == 'T') {
          val.SetX();
        } else {
          LOG_ASSERT(false, "ERROR:not support stil po vals:" + elem);
        }
        pc.AddOutputVal(i, val);
      }
    }
    iss.str("");
  }
}

void PatternParser::TransformSerialToParallelForSAF(
    std::vector<Pattern>& serial_patterns,
    std::vector<Pattern>& parallel_patterns, bool is_external) {
  parallel_patterns.clear();
  Pattern sp;
  InitPatternWithX(sp);
  int pattern_nums = serial_patterns.size();
  for (int i = 0; i < pattern_nums; i++) {
    if (i % PARALLEL_BIT_SIZE == PARALLEL_BIT_SIZE - 1) {
      FillPatternPacketForSAF(i, sp, serial_patterns, is_external);
      sp.SetPatternCycleNums();
      parallel_patterns.emplace_back(sp);
      sp.Clear();
      InitPatternWithX(sp);
    } else {
      FillPatternPacketForSAF(i, sp, serial_patterns, is_external);
    }
  }
  // the last group pattern which size < 64
  if (pattern_nums % PARALLEL_BIT_SIZE != 0) {
    sp.SetPatternCycleNums();
    parallel_patterns.emplace_back(sp);
    sp.Clear();
  }
  int begin_pattern_id = 0;
  for (auto i = 0; i < parallel_patterns.size(); ++i) {
    parallel_patterns[i].SetPatternId(begin_pattern_id);
    begin_pattern_id +=
        parallel_patterns[i].GetPatternCycle(0).LastParallelBitIdx();
  }
}


void PatternParser::InitPatternCycleWithX(PatternCycle& patternCycle,
                                          int frameSize_or_scanChainSize,
                                          int piSize_or_siSize,
                                          int poSize_or_soSize) {
  patternCycle.AllocInputMem(frameSize_or_scanChainSize);
  patternCycle.AllocOutputMem(frameSize_or_scanChainSize);
  for (int i = 0; i < frameSize_or_scanChainSize; i++) {
    if (patternCycle.GetPatternCycleType() ==
        PatternCycleType::LOAD_UNLOAD_PCT) {
      piSize_or_siSize = pnlist_->ScanChains()[i]->GetChainLength();
    }
    patternCycle.GetInputVal()[i].resize(piSize_or_siSize);
    patternCycle.GetOutputVal()[i].resize(poSize_or_soSize);
  }
}

void PatternParser::InitPatternWithX(Pattern& sp, int capture_size,
                                     int load_unload_size) {
  // todo
  if (ptype_ == PatternType::COMB_PT) {
    // com circuit only has pi po values
    PatternCycle pc;
    pc.SetPatternCycleType(PatternCycleType::CAPTURE_PCT);
    InitPatternCycleWithX(pc, 1, pnlist_->NumPIs(), pnlist_->NumPOs());
    sp.AddPattern(pc);
  } else if (ptype_ == PatternType::BASIC_SCAN_PT) {
    // full scan circuit has load_unload and pi_po values
    // 1. load_unload
    PatternCycle pc_load_unload;
    pc_load_unload.SetPatternCycleType(PatternCycleType::LOAD_UNLOAD_PCT);
    auto scanChains = pnlist_->ScanChains();
    // siSize and soSize will be fixed in function.
    InitPatternCycleWithX(pc_load_unload, scanChains.size(), 0, 0);
    sp.AddPattern(pc_load_unload);
    // 2. pi_po
    PatternCycle pc_capture;
    pc_capture.SetPatternCycleType(PatternCycleType::CAPTURE_PCT);
    InitPatternCycleWithX(pc_capture, MAX_FRAME_NUM, pnlist_->NumPIs(),
                          pnlist_->NumPOs());
    sp.AddPattern(pc_capture);
  }
}

void PatternParser::InitPatternWithX(Pattern& sp) {
  if (ptype_ == PatternType::COMB_PT) {
    // com circuit only has pi po values
    PatternCycle pc;
    pc.SetPatternCycleType(PatternCycleType::CAPTURE_PCT);
    pc.AllocInputMem(1);
    pc.AllocOutputMem(1);
    pc.GetInputVal()[0].resize(pnlist_->NumPIs());
    pc.GetOutputVal()[0].resize(pnlist_->NumPOs());
    sp.AddPattern(pc);
  } else if (ptype_ == PatternType::BASIC_SCAN_PT) {
    // full scan circuit has load_unload and pi_po values
    // 1. load_unload
    PatternCycle pc_load_unload;
    pc_load_unload.SetPatternCycleType(PatternCycleType::LOAD_UNLOAD_PCT);
    // the nums of scan chain
    int scan_chain_size = pnlist_->ScanChains().size();
    pc_load_unload.AllocInputMem(scan_chain_size);
    pc_load_unload.AllocOutputMem(scan_chain_size);
    for (int i = 0; i < scan_chain_size; i++) {
      int single_sc_len = pnlist_->ScanChains()[i]->GetChainLength();
      pc_load_unload.GetInputVal()[i].resize(single_sc_len);
      pc_load_unload.GetOutputVal()[i].resize(single_sc_len);
    }
    sp.AddPattern(pc_load_unload);
    // 2. pi_po
    PatternCycle pc_capture;
    pc_capture.SetPatternCycleType(PatternCycleType::CAPTURE_PCT);
    pc_capture.AllocInputMem(MAX_FRAME_NUM);
    pc_capture.AllocOutputMem(MAX_FRAME_NUM);
    for (int i = 0; i < pc_capture.GetSize(); i++) {
      pc_capture.GetInputVal()[i].resize(pnlist_->NumPIs());
      pc_capture.GetOutputVal()[i].resize(pnlist_->NumPOs());
    }
    sp.AddPattern(pc_capture);
  }
}

void PatternParser::FillPatternPacketForSAF(
    int pattern_id, Pattern& sp, std::vector<Pattern>& serial_patterns,
    bool is_external) {
  if (ptype_ == PatternType::COMB_PT) {
    // com circuit no 3-frames
    PatternCycle& pc_parallel = sp.GetPatternCycle(0);
    PatternCycle& pc_serial = serial_patterns[pattern_id].GetPatternCycle(0);
    // fill pi
    for (int pi_id = 0; pi_id < pnlist_->NumPIs(); pi_id++) {
      if (pc_serial.GetInputVal()[0][pi_id].v0_ == V64_ONE) {
        pc_parallel.GetInputVal()[0][pi_id].v0_ =
            (pc_parallel.GetInputVal()[0][pi_id].v0_ << 1) + 1;
      } else {
        pc_parallel.GetInputVal()[0][pi_id].v0_ =
            (pc_parallel.GetInputVal()[0][pi_id].v0_ << 1) + 0;
      }
      if (pc_serial.GetInputVal()[0][pi_id].v1_ == V64_ONE) {
        pc_parallel.GetInputVal()[0][pi_id].v1_ =
            (pc_parallel.GetInputVal()[0][pi_id].v1_ << 1) + 1;
      } else {
        pc_parallel.GetInputVal()[0][pi_id].v1_ =
            (pc_parallel.GetInputVal()[0][pi_id].v1_ << 1) + 0;
      }
    }
    if (is_external) {
      // fill po, po only the first frame is needed,other frame is X
      for (int po_id = 0; po_id < pnlist_->NumPOs(); po_id++) {
        if (pc_serial.GetOutputVal()[0][po_id].v0_ == V64_ONE) {
          pc_parallel.GetOutputVal()[0][po_id].v0_ =
              (pc_parallel.GetOutputVal()[0][po_id].v0_ << 1) + 1;
        } else {
          pc_parallel.GetOutputVal()[0][po_id].v0_ =
              (pc_parallel.GetOutputVal()[0][po_id].v0_ << 1) + 0;
        }
        if (pc_serial.GetOutputVal()[0][po_id].v1_ == V64_ONE) {
          pc_parallel.GetOutputVal()[0][po_id].v1_ =
              (pc_parallel.GetOutputVal()[0][po_id].v1_ << 1) + 1;
        } else {
          pc_parallel.GetOutputVal()[0][po_id].v1_ =
              (pc_parallel.GetOutputVal()[0][po_id].v1_ << 1) + 0;
        }
      }
    }
    // update the last bit position
    pc_parallel.LastParallelBitIdx()++;
  } else if (ptype_ == PatternType::BASIC_SCAN_PT) {
    // todo
    PatternCycle& pc_parallel_load_unload = sp.GetPatternCycle(0);
    PatternCycle& pc_serial_load_unload =
        serial_patterns[pattern_id].GetPatternCycle(0);
    int scan_chain_size = pnlist_->ScanChains().size();
    ;
    // fill pi
    for (int i = 0; i < scan_chain_size; i++) {
      int single_sc_len = pnlist_->ScanChains()[i]->GetChainLength();
      for (int pi_id = 0; pi_id < single_sc_len; pi_id++) {
        if (pc_serial_load_unload.GetInputVal()[i][pi_id].v0_ == V64_ONE) {
          pc_parallel_load_unload.GetInputVal()[i][pi_id].v0_ =
              (pc_parallel_load_unload.GetInputVal()[i][pi_id].v0_ << 1) + 1;
        } else {
          pc_parallel_load_unload.GetInputVal()[i][pi_id].v0_ =
              (pc_parallel_load_unload.GetInputVal()[i][pi_id].v0_ << 1) + 0;
        }

        if (pc_serial_load_unload.GetInputVal()[i][pi_id].v1_ == V64_ONE) {
          pc_parallel_load_unload.GetInputVal()[i][pi_id].v1_ =
              (pc_parallel_load_unload.GetInputVal()[i][pi_id].v1_ << 1) + 1;
        } else {
          pc_parallel_load_unload.GetInputVal()[i][pi_id].v1_ =
              (pc_parallel_load_unload.GetInputVal()[i][pi_id].v1_ << 1) + 0;
        }
      }
      if (is_external) {
        // fill po, po only the first frame is needed,other frame is X
        for (int po_id = 0; po_id < single_sc_len; po_id++) {
          if (pc_serial_load_unload.GetOutputVal()[i][po_id].v0_ == V64_ONE) {
            pc_parallel_load_unload.GetOutputVal()[i][po_id].v0_ =
                (pc_parallel_load_unload.GetOutputVal()[i][po_id].v0_ << 1) + 1;
          } else {
            pc_parallel_load_unload.GetOutputVal()[i][po_id].v0_ =
                (pc_parallel_load_unload.GetOutputVal()[i][po_id].v0_ << 1) + 0;
          }

          if (pc_serial_load_unload.GetOutputVal()[i][po_id].v1_ == V64_ONE) {
            pc_parallel_load_unload.GetOutputVal()[i][po_id].v1_ =
                (pc_parallel_load_unload.GetOutputVal()[i][po_id].v1_ << 1) + 1;
          } else {
            pc_parallel_load_unload.GetOutputVal()[i][po_id].v1_ =
                (pc_parallel_load_unload.GetOutputVal()[i][po_id].v1_ << 1) + 0;
          }
        }
      }
    }
    // update the last bit position
    pc_parallel_load_unload.LastParallelBitIdx()++;
    PatternCycle& pc_parallel_capture = sp.GetPatternCycle(1);
    PatternCycle& pc_serial_capture =
        serial_patterns[pattern_id].GetPatternCycle(1);
    int frame_size = pc_serial_capture.GetSize();
    if (frame_size == 1) {
      for (int pi_id = 0; pi_id < pnlist_->NumPIs(); pi_id++) {
        // frame-0
        if (pc_serial_capture.GetInputVal()[0][pi_id].v0_ == V64_ONE) {
          pc_parallel_capture.GetInputVal()[0][pi_id].v0_ =
              (pc_parallel_capture.GetInputVal()[0][pi_id].v0_ << 1) + 1;
        } else {
          pc_parallel_capture.GetInputVal()[0][pi_id].v0_ =
              (pc_parallel_capture.GetInputVal()[0][pi_id].v0_ << 1) + 0;
        }
        if (pc_serial_capture.GetInputVal()[0][pi_id].v1_ == V64_ONE) {
          pc_parallel_capture.GetInputVal()[0][pi_id].v1_ =
              (pc_parallel_capture.GetInputVal()[0][pi_id].v1_ << 1) + 1;
        } else {
          pc_parallel_capture.GetInputVal()[0][pi_id].v1_ =
              (pc_parallel_capture.GetInputVal()[0][pi_id].v1_ << 1) + 0;
        }
        // frame-1
        if (pc_serial_capture.GetInputVal()[0][pi_id].v0_ == V64_ONE) {
          pc_parallel_capture.GetInputVal()[1][pi_id].v0_ =
              (pc_parallel_capture.GetInputVal()[1][pi_id].v0_ << 1) + 1;
        } else {
          pc_parallel_capture.GetInputVal()[1][pi_id].v0_ =
              (pc_parallel_capture.GetInputVal()[1][pi_id].v0_ << 1) + 0;
        }
        if (pc_serial_capture.GetInputVal()[0][pi_id].v1_ == V64_ONE) {
          pc_parallel_capture.GetInputVal()[1][pi_id].v1_ =
              (pc_parallel_capture.GetInputVal()[1][pi_id].v1_ << 1) + 1;
        } else {
          pc_parallel_capture.GetInputVal()[1][pi_id].v1_ =
              (pc_parallel_capture.GetInputVal()[1][pi_id].v1_ << 1) + 0;
        }
        // frame-2
        if (pc_serial_capture.GetInputVal()[0][pi_id].v0_ == V64_ONE) {
          pc_parallel_capture.GetInputVal()[2][pi_id].v0_ =
              (pc_parallel_capture.GetInputVal()[2][pi_id].v0_ << 1) + 1;
        } else {
          pc_parallel_capture.GetInputVal()[2][pi_id].v0_ =
              (pc_parallel_capture.GetInputVal()[2][pi_id].v0_ << 1) + 0;
        }
        if (pc_serial_capture.GetInputVal()[0][pi_id].v1_ == V64_ONE) {
          pc_parallel_capture.GetInputVal()[2][pi_id].v1_ =
              (pc_parallel_capture.GetInputVal()[2][pi_id].v1_ << 1) + 1;
        } else {
          pc_parallel_capture.GetInputVal()[2][pi_id].v1_ =
              (pc_parallel_capture.GetInputVal()[2][pi_id].v1_ << 1) + 0;
        }
        //                pc_parallel_capture.GetInputVal()[1][pi_id].v0_ =
        //                        pc_parallel_capture.GetInputVal()[1][pi_id].v0_;
        //                pc_parallel_capture.GetInputVal()[2][pi_id].v0_ =
        //                        pc_parallel_capture.GetInputVal()[2][pi_id].v0_;
        //
        //                pc_parallel_capture.GetInputVal()[1][pi_id].v1_ =
        //                        pc_parallel_capture.GetInputVal()[1][pi_id].v1_;
        //                pc_parallel_capture.GetInputVal()[2][pi_id].v1_ =
        //                        pc_parallel_capture.GetInputVal()[2][pi_id].v1_;
      }
    } else if (frame_size == MAX_FRAME_NUM) {
      for (int frame_id = 0; frame_id < frame_size; frame_id++) {
        for (int pi_id = 0; pi_id < pnlist_->NumPIs(); pi_id++) {
          if (pc_serial_capture.GetInputVal()[frame_id][pi_id].v0_ == V64_ONE) {
            pc_parallel_capture.GetInputVal()[frame_id][pi_id].v0_ =
                (pc_parallel_capture.GetInputVal()[frame_id][pi_id].v0_ << 1) +
                1;
          } else {
            pc_parallel_capture.GetInputVal()[frame_id][pi_id].v0_ =
                (pc_parallel_capture.GetInputVal()[frame_id][pi_id].v0_ << 1) +
                0;
          }
          if (pc_serial_capture.GetInputVal()[frame_id][pi_id].v1_ == V64_ONE) {
            pc_parallel_capture.GetInputVal()[frame_id][pi_id].v1_ =
                (pc_parallel_capture.GetInputVal()[frame_id][pi_id].v1_ << 1) +
                1;
          } else {
            pc_parallel_capture.GetInputVal()[frame_id][pi_id].v1_ =
                (pc_parallel_capture.GetInputVal()[frame_id][pi_id].v1_ << 1) +
                0;
          }
        }
      }
    }
    if (is_external) {
      // fill po, po only the first frame is needed,other frame is X
      for (int po_id = 0; po_id < pnlist_->NumPOs(); po_id++) {
        if (pc_serial_capture.GetOutputVal()[0][po_id].v0_ == V64_ONE) {
          pc_parallel_capture.GetOutputVal()[0][po_id].v0_ =
              (pc_parallel_capture.GetOutputVal()[0][po_id].v0_ << 1) + 1;
        } else {
          pc_parallel_capture.GetOutputVal()[0][po_id].v0_ =
              (pc_parallel_capture.GetOutputVal()[0][po_id].v0_ << 1) + 0;
        }
        if (pc_serial_capture.GetOutputVal()[0][po_id].v1_ == V64_ONE) {
          pc_parallel_capture.GetOutputVal()[0][po_id].v1_ =
              (pc_parallel_capture.GetOutputVal()[0][po_id].v1_ << 1) + 1;
        } else {
          pc_parallel_capture.GetOutputVal()[0][po_id].v1_ =
              (pc_parallel_capture.GetOutputVal()[0][po_id].v1_ << 1) + 0;
        }
      }
    }
    pc_parallel_capture.LastParallelBitIdx()++;
  }
}


void PatternParser::TransformParallelToSerial(
    std::vector<Pattern>& parallel_patterns,
    std::vector<Pattern>& serial_patterns) {
  serial_patterns.clear();
  serial_patterns.reserve(parallel_patterns.size() * PARALLEL_BIT_SIZE);
  for (int parall_id = 0; parall_id < parallel_patterns.size();
       ++parall_id) {                                     // all parall_pattern
    Pattern& curr_parall = parallel_patterns[parall_id];  // a parall pattern
    std::vector<PatternCycle>& pcs = curr_parall.GetPattern();
    int pattern_len = pcs[0].LastParallelBitIdx();

    for (int p_idx = 0; p_idx < pattern_len; ++p_idx) {
      Pattern sp_temp;
      for (int i = 0; i < pcs.size(); ++i) {
        PatternCycle& pc = pcs[i];
        std::vector<std::vector<val64_t>>& input = pc.GetInputVal();
        std::vector<std::vector<val64_t>>& output = pc.GetOutputVal();
        PatternCycle pc_temp;
        // init date member
        pc_temp.SetPatternCycleType(pc.GetPatternCycleType());
        pc_temp.LastParallelBitIdx() = 0;
        pc_temp.AllocInputMem(pc.GetSize());
        pc_temp.AllocOutputMem(pc.GetSize());

        if (pc.GetPatternCycleType() == PatternCycleType::LOAD_UNLOAD_PCT) {
          for (int sc_id = 0; sc_id < input.size(); ++sc_id) {
            for (int sc_dff_id = 0; sc_dff_id < input[sc_id].size();
                 ++sc_dff_id) {
              val64_t parall_loadval = input[sc_id][sc_dff_id];
              val64_t serial_loadval =
                  GetOneBit(parall_loadval, pattern_len - 1 - p_idx);
              pc_temp.AddInputVal(sc_id, serial_loadval);
            }
            for (int sc_dff_id = 0; sc_dff_id < output[sc_id].size();
                 ++sc_dff_id) {
              val64_t parall_unloadval = output[sc_id][sc_dff_id];
              val64_t serial_unloadval =
                  GetOneBit(parall_unloadval, pattern_len - 1 - p_idx);
              pc_temp.AddOutputVal(sc_id, serial_unloadval);
            }
          }
        } else if (pc.GetPatternCycleType() == PatternCycleType::CAPTURE_PCT) {
          for (int frame_id = 0; frame_id < input.size(); ++frame_id) {
            for (int pi_id = 0; pi_id < input[frame_id].size(); ++pi_id) {
              val64_t parall_pival = input[frame_id][pi_id];
              val64_t serial_pival =
                  GetOneBit(parall_pival, pattern_len - 1 - p_idx);
              pc_temp.AddInputVal(frame_id, serial_pival);
            }
            for (int po_id = 0; po_id < output[frame_id].size(); ++po_id) {
              val64_t parall_poval = output[frame_id][po_id];
              val64_t serial_poval =
                  GetOneBit(parall_poval, pattern_len - 1 - p_idx);
              pc_temp.AddOutputVal(frame_id, serial_poval);
            }
          }
        } else {
          LOG_ASSERT(false, "ERROR:not support pattern cycle type");
        }
        sp_temp.AddPattern(pc_temp);
      }
      sp_temp.SetPatternCycleNums();
      serial_patterns.emplace_back(sp_temp);
    }
  }
  for (auto i = 0; i < serial_patterns.size(); ++i) {
    serial_patterns[i].SetPatternId(i);
  }
}

val64_t PatternParser::GetOneBit(const val64_t& source, int32_t bit_index) {
  const uint64_t right1 = 1;
  uint64_t b0 = source.v0_ & (right1 << bit_index);
  uint64_t b1 = source.v1_ & (right1 << bit_index);

  if (b0 > 0) {
    b0 = 0xFFFFFFFFFFFFFFFF >> (64 - PARALLEL_BIT_SIZE);
  }
  if (b1 > 0) {
    b1 = 0xFFFFFFFFFFFFFFFF >> (64 - PARALLEL_BIT_SIZE);
  }

  return val64_t(b0, b1);
}

}  // namespace ictest