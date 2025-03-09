// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "simulation/GoodSimulator.h"
namespace ictest {

/// 分配保存逻辑值的数组空间 建立dff slave==>master的虚拟映射
void GoodSimulator::AllocAllMem() {
  auto& dffs = prim_netlist_->GetDFFGates();
  slave2master_.resize(prim_netlist_->NumGates());
  for (int i = 0; i < dffs.size(); ++i) {
    slave2master_[dffs[i]->GetGId()] = i;
  }
  good_machine_vals_.resize(MAX_FRAME_NUM,
                            std::vector<val64_t>(prim_netlist_->NumGates()));
  good_master_vals_.resize(MAX_FRAME_NUM,
                           std::vector<val64_t>(prim_netlist_->NumDFFs()));
}

/// 重置保存逻辑仿真逻辑值的数组的所有值为X
///@attention 同时重置了虚拟master的值
void GoodSimulator::ResetGoodMachine() {
  std::fill(good_master_vals_[FRAME_0].begin(),
            good_master_vals_[FRAME_0].end(), val64_t(0, 0));
  std::fill(good_machine_vals_[FRAME_0].begin(),
            good_machine_vals_[FRAME_0].end(), val64_t(0, 0));
}
/**
 * @brief 重置所有锁存器的逻辑值
 * 将该锁存器及其扇出加入事件队列等待重新计算
 */
void GoodSimulator::ResetDLatch(const val64_t& init_val, int frame_id) {
  // latch
  for (int i = 0; i < prim_netlist_->NumDLATs(); ++i) {
    Gate& latch_ptr = *prim_netlist_->GetDLATGates()[i];

    if (latch_ptr.GetDPI() >= 0) {
      bool is_update = UpdateGoodVal(init_val, latch_ptr.GetGId(), frame_id);
      if (is_update) {
        AddGateToEventQueue(latch_ptr);
        event_queue_.push(latch_ptr.GetGId(), latch_ptr.GetDPI());
      }
    }
  }
}
/**
 * @brief 重置所有非扫描单元dff的值
 */
void GoodSimulator::ResetNonScanDFF(const val64_t& init_val, int frame_id) {
  // latch and non scan dff may could ctrl
  const auto& nsc_dffs = prim_netlist_->GetNonScanDFFGates();
  for (int i = 0; i < nsc_dffs.size(); ++i) {
    auto& dff_gptr = *nsc_dffs[i];
    int dff_id = dff_gptr.GetGId();
    int master_dff_id = slave2master_[dff_id];

    if (good_master_vals_[frame_id][master_dff_id] != init_val) {
      good_master_vals_[frame_id][master_dff_id] = init_val;
      event_queue_.push(dff_id, dff_gptr.GetDPI());
    }

    bool is_update = UpdateGoodVal(init_val, dff_id, frame_id);
    if (is_update) {
      AddGateToEventQueue(dff_gptr);
      event_queue_.push(dff_id, dff_gptr.GetDPI());
    }
  }
}
/**
 * @brief 重置所有dff的值
 */
void GoodSimulator::ResetAllDFF(const val64_t& init_val, int frame_id) {
  const auto& dffs = prim_netlist_->GetDFFGates();
  for (int i = 0; i < dffs.size(); ++i) {
    auto& dff_gptr = *dffs[i];
    int dff_id = dff_gptr.GetGId();
    int master_dff_id = slave2master_[dff_id];

    if (good_master_vals_[frame_id][master_dff_id] != init_val) {
      good_master_vals_[frame_id][master_dff_id] = init_val;
      event_queue_.push(dff_id, dff_gptr.GetDPI());
    }

    bool is_update = UpdateGoodVal(init_val, dff_id, frame_id);
    if (is_update) {
      AddGateToEventQueue(dff_gptr);
      event_queue_.push(dff_id, dff_gptr.GetDPI());
    }
  }
}
/**
 * @brief 逻辑仿真主函数 支持组合电路/全扫描电路仿真
 */
void GoodSimulator::GoodSimulation(Pattern& sp) {
  // include load_unload & capture procedure
  // init simulation
  ResetDLatch(val64_t(0, 0), FRAME_0);
  ResetNonScanDFF(val64_t(0, 0), FRAME_0);
  auto& last_pc = sp.GetPattern().back();
  int load_pc_id = -1;
  for (int pc_id = 0; pc_id < sp.GetPattern().size(); ++pc_id) {
    auto& pattern_cycle = sp.GetPattern()[pc_id];

    InitTieVal(pattern_cycle.LastParallelBitIdx(), 0);
    // init scan dff load_unload value
    PatternCycleType pct = pattern_cycle.GetPatternCycleType();
    switch (pct) {
      case PatternCycleType::LOAD_UNLOAD_PCT:
        load_pc_id = pc_id;
        Load(pattern_cycle.GetInputVal(), 0);
        break;
      case PatternCycleType::CAPTURE_PCT:
        Capture(pattern_cycle);
        break;
      default:
        LOG_ASSERT(false,
                   "Not support pct : " + std::to_string(static_cast<int>(pct)))
    }
  }
}


///load 扫描链的值 只支持并行load
void GoodSimulator::Load(std::vector<std::vector<val64_t>>& load_vals,
                         int shift_nums) {
  const int max_shift_cycle_nums = simdata_->max_chain_lenth_;

  // 给扫描链dff直接赋值
  ParallelAssignDffs(load_vals, shift_nums);
}

/**
 * @brief load过程中用于并行赋值的函数部分
 * @param shift_nums 需要仿真的shift过程次数，取值区间为[0~最大扫描链长度]
 * @param load_vals 给扫描链赋值的二维数组向量 load_vals[i]赋值给第i条扫描链
 * @attention 当前只支持复制到frame 0，后续有需要时再做拓展
 */
void GoodSimulator::ParallelAssignDffs(
    const std::vector<std::vector<val64_t>>& load_vals, int shift_nums) {
  for (int chain_id = 0; chain_id < simdata_->GetChainNums(); ++chain_id) {
    int max_assign_idx =
        prim_netlist_->ScanChains()[chain_id]->GetChainLength() - shift_nums;
    const auto& sc_dffs =
        prim_netlist_->ScanChains()[chain_id]->ScanChainDFFs();

    for (int inpux_idx = 0; inpux_idx < max_assign_idx; ++inpux_idx) {
      // 因为后续需要做shift过程 所以这里赋值不能直接赋值到最终对应的触发器
      // 需要有shift_nums的偏移
      int dff_idx = inpux_idx + shift_nums;
      Gate& slave_dff_gptr = *sc_dffs[dff_idx];
      int master_dff_id = slave2master_[slave_dff_gptr.GetGId()];
      const val64_t& dff_val =prim_netlist_->GetDffInvFlag()[slave_dff_gptr.GetGId()]?~ load_vals[chain_id][inpux_idx]: load_vals[chain_id][inpux_idx];
      // update dff master value and slave value
      if (good_master_vals_[FRAME_0][master_dff_id] != dff_val) {
        good_master_vals_[FRAME_0][master_dff_id] = dff_val;
        event_queue_.push(slave_dff_gptr.GetGId(), slave_dff_gptr.GetDPI());
      }
      bool is_update = UpdateGoodVal(dff_val, slave_dff_gptr.GetGId(), FRAME_0);
      if (is_update) {
        AddGateToEventQueue(slave_dff_gptr);
        event_queue_.push(slave_dff_gptr.GetGId(), slave_dff_gptr.GetDPI());
      }
    }
  }
}
/**
 * @brief 只支持纯并行扫描链赋值的逻辑仿真函数 可以完全用load 函数代替
 * @attention 后续可能会删除该函数 应减少使用
 */
void GoodSimulator::LoadScanDFFs(PatternCycle& pc, int frame_id) {
  // init scan chain
  assert(pc.GetSize() == prim_netlist_->ScanChains().size());

  for (int chain_id = 0; chain_id < prim_netlist_->ScanChains().size();
       chain_id++) {
    auto& single_chain = prim_netlist_->ScanChains()[chain_id];
    std::vector<Gate*>& single_chain_dffs = single_chain->ScanChainDFFs();
    const auto& load_vals = pc.GetInputVal()[chain_id];

    for (int i = 0; i < single_chain_dffs.size(); i++) {
      auto& dff_gptr = *single_chain_dffs[i];
      int dff_id = dff_gptr.GetGId();
      const auto& good_val = prim_netlist_->GetDffInvFlag()[dff_id]?~load_vals[i]:load_vals[i];
      int master_dff_id = slave2master_[dff_id];

      if (good_master_vals_[frame_id][master_dff_id] != good_val) {
        good_master_vals_[frame_id][master_dff_id] = good_val;
        event_queue_.push(dff_id, dff_gptr.GetDPI());
      }
      bool is_update = UpdateGoodVal(good_val, dff_id, frame_id);
      if (is_update) {
        AddGateToEventQueue(dff_gptr);
        event_queue_.push(dff_id, dff_gptr.GetDPI());
      }
    }
  }
}

/**
 * @brief Capture过程主函数 包括force pi/pulse clk/measure po
 * @param pattern_cycle capture向量
 * @note measure po用于跟外部向量保存的结果作对比调试逻辑仿真的正确性
 * 当前仅支持在frame 0 measure，在全局参数option中设置compare_po_value_为true时对比
 */
void GoodSimulator::Capture(PatternCycle& pattern_cycle) {
  for (int frame_id = 0; frame_id < pattern_cycle.GetSize(); frame_id++) {
    if (frame_id == FRAME_0) {
      ForcePIs(pattern_cycle.GetInputVal()[frame_id], 0);
      EventDriveSim(frame_id);
    } else {
      // frame-1 and frame-2 simulation
      good_machine_vals_[frame_id] = good_machine_vals_[frame_id - 1];
      good_master_vals_[frame_id] = good_master_vals_[frame_id - 1];
      PulseCaptureClock(pattern_cycle, frame_id);
      EventDriveSim(frame_id);
    }
  }  // end 3-frame
}

/**
 * @brief 给PI端口赋值的三个重载函数之一 给pi和tie端口赋值
 * 并将值发生变化的pi和tie端口的扇出加入事件对列
 * @param pattern_cycle capture向量
 * @param frame_id
 * 取capture向量的frame_id对应的输入值赋值到电路的frame_id对应的值数组
 * @attention 将对tie端口的赋值和pi的赋值集成到了一个函数 可能不太好
 * 后面准备弃用
 */
void GoodSimulator::ForcePIs(PatternCycle& pattern_cycle, int frame_id) {
  auto& pi_vals = pattern_cycle.GetInputVal()[frame_id];

  auto& tie_gates = prim_netlist_->GetTieGates();

  auto& pi_gates = prim_netlist_->GetPIGates();

  // force pi
  for (int idx = 0; idx < pi_gates.size(); idx++) {
    auto& pi_gptr = *pi_gates[idx];
    const auto& pi_val = pi_vals[idx];

    bool is_update = UpdateGoodVal(pi_val, pi_gptr.GetGId(), frame_id);
    if (is_update) {
      AddGateToEventQueue(pi_gptr);
    }
  }

  // init tie gate values
  constexpr uint64_t all_one = ~0;
  int useful_bits = pattern_cycle.LastParallelBitIdx() == 0
                        ? 64
                        : pattern_cycle.LastParallelBitIdx();

  const int nouseful_bits = 64 - useful_bits;

  const uint64_t pc_all_one = all_one >> nouseful_bits;

  val64_t tie_val;
  for (int i = 0; i < tie_gates.size(); i++) {
    auto& tie_gate = *tie_gates[i];
    auto logic_type = tie_gate.GetGType();

    switch (logic_type) {
      case GType::G_TIE0:
        tie_val = tie_val.MaskTo0(pc_all_one);
        break;
      case GType::G_TIE1:
        tie_val = tie_val.MaskTo1(pc_all_one);
        break;
      case GType::G_TIEX:
        tie_val = tie_val.MaskToX(pc_all_one);
        break;
      default:
        LOG_ASSERT(false, "not support tie gate type")
    }

    bool is_update = UpdateGoodVal(tie_val, tie_gate.GetGId(), frame_id);
    if (is_update) {
      AddGateToEventQueue(tie_gate);
    }
  }
}

/**
 * @brief 给PI端口赋值的三个重载函数之一 给pi端口赋值
 * 并将值发生变化的pi端口的扇出加入事件对列
 * @param pi_vals 长度和pi端口数量对应相等 表示给pi端口赋值的向量基本输入单位
 * @param frame_id 赋值到电路的frame_id对应的值数组
 * @attention 将对tie端口的赋值拆分出来 只对pi端口进行处理
 */
void GoodSimulator::ForcePIs(const std::vector<val64_t>& pi_vals,
                             int frame_id) {
  const auto& pi_gates = prim_netlist_->GetPIGates();
  // force pi
  for (int idx = 0; idx < pi_gates.size(); idx++) {
    auto& pi_gptr = *pi_gates[idx];
    const auto& pi_val = pi_vals[idx];

    bool is_update = UpdateGoodVal(pi_val, pi_gptr.GetGId(), frame_id);
    if (is_update) {
      AddGateToEventQueue(pi_gptr);
    }
  }
}
/**
 * @brief 给PI端口赋值的三个重载函数之一 给特定pi端口赋值
 * @param pi2vals pi的逻辑门id==>pi端口初始值的映射数组
 * @param frame_id 赋值到电路的frame_id对应的值数组
 */
void GoodSimulator::ForcePIs(std::vector<std::pair<int, val64_t>>& pi2vals,
                             int frame_id) {
  const auto& netlist = prim_netlist_->GetPrimNetlist();
  for (const auto& hold_sig : pi2vals) {
    Gate& active_gptr = *netlist[hold_sig.first];
    bool is_update = UpdateGoodVal(hold_sig.second, hold_sig.first, frame_id);
    if (is_update) {
      AddGateToEventQueue(active_gptr);
    }
  }
}

/**
 * @brief 初始化电路tie端口
 * @param useful_bits 取值范围[0~64] 一般不会取0 除非初始化所有位为X
 * 初始化一个val64_t中的多少个比特位 从最低位开始计算 未填充的高位用X填充
 * @param frame_id 初始化到电路的frame_id对应的值数组
 */
void GoodSimulator::InitTieVal(int useful_bits, int frame_id) {
  const auto& tie_gates = prim_netlist_->GetTieGates();
  // init tie gate values
  constexpr uint64_t all_one = ~0;
  const int nouseful_bits = 64 - useful_bits;
  const uint64_t pc_all_one = all_one >> nouseful_bits;

  val64_t tie_val;
  for (int i = 0; i < tie_gates.size(); i++) {
    auto& tie_gate = *tie_gates[i];
    auto logic_type = tie_gate.GetGType();

    switch (logic_type) {
      case GType::G_TIE0:
        tie_val = tie_val.MaskTo0(pc_all_one);
        break;
      case GType::G_TIE1:
        tie_val = tie_val.MaskTo1(pc_all_one);
        break;
      case GType::G_TIEX:
        tie_val = tie_val.MaskToX(pc_all_one);
        break;
      default:
        LOG_ASSERT(false, "not support tie gate type")
    }

    bool is_update = UpdateGoodVal(tie_val, tie_gate.GetGId(), frame_id);
    if (is_update) {
      AddGateToEventQueue(tie_gate);
    }
  }
}

void GoodSimulator::EventDriveSim(int frame_id) {
  auto& netlist = prim_netlist_->GetPrimNetlist();
  event_queue_.InitCurrLevel(1);

  while (!event_queue_.empty()) {
    int gate_id = event_queue_.pop();
    events_count_++;
    auto& active_gptr = *netlist[gate_id];
    auto logic_type = active_gptr.GetGType();

    // calculate new out put
    SetGateInputVal(active_gptr, frame_id);
    auto output = evaluate_.Evaluate(logic_type);

    // update value
    if (logic_type == GType::G_DFF) {
      good_master_vals_[frame_id][slave2master_[gate_id]] =
          evaluate_.GetDffMasterVal();
    }

    bool is_update = UpdateGoodVal(output, gate_id, frame_id);

    // update event queue
    if (is_update) {
      AddGateToEventQueue(active_gptr);
    }
  }
}

/// 找到时钟所在的pi value 输入到电路中
///@attention 输入向量应该保证在一个capture过程中
///电路只有时钟（包括set/reset）端口的值可以发生改变
void GoodSimulator::PulseCaptureClock(PatternCycle& pc, int frame_id) {
  for (auto& clk_pair : simdata_->clk_to_pi_id_) {
    int clk_id = clk_pair.first;
    int clk_pi_idx = clk_pair.second;
    const auto& clk_val = pc.GetInputVal()[frame_id][clk_pi_idx];
    bool is_pulse = UpdateGoodVal(clk_val, clk_id, frame_id);
    if (is_pulse) {
      AddGateToEventQueue(*prim_netlist_->GetPrimNetlist()[clk_id]);
    }
  }
}

///////////////////////////////////////////////////////////
// set and get value
/**
 * @brief 计算之前向evaluate类赋值
 * @param gptr 当前待计算的逻辑门
 * @param frame_id 待计算的逻辑门的取值来自于值数组的哪个frame
 * @attention dff和dlatch只支持四输入 如果输入数量不为4则报错
 */
inline void GoodSimulator::SetGateInputVal(Gate& gptr, int frame_id) {
  size_t fanin_size = gptr.FaninSize();

  LOG_ASSERT(fanin_size != 0, "Wrong fanin size")

  auto logic_type = gptr.GetGType();
  auto* fanin_data = gptr.FaninGates().data();

  evaluate_.SetFaninSize(fanin_size);
  for (int i = 0; i < fanin_size; ++i) {
    evaluate_.SetPinValue(
        i, good_machine_vals_[frame_id][fanin_data[i]->GetGId()]);
  }

  switch (logic_type) {
    case GType::G_DFF:
      LOG_ASSERT(gptr.FaninSize() == 4, "DFF fanin size ERROR")
      evaluate_.SetDffMasterVal(
          good_master_vals_[frame_id][slave2master_[gptr.GetGId()]]);
      evaluate_.SetDffSlaveVal(good_machine_vals_[frame_id][gptr.GetGId()]);
      break;

    case GType::G_DLAT:
      LOG_ASSERT(gptr.FaninSize() == 4, "G_DLAT fanin size ERROR")
      evaluate_.SetLatchHoldVal(good_machine_vals_[frame_id][gptr.GetGId()]);
      break;

    default:
      break;
  }
}

bool GoodSimulator::UpdateGoodVal(const val64_t& new_value, gid_t gid,
                                  int frame_id) {
  if (new_value != good_machine_vals_[frame_id][gid]) {
    good_machine_vals_[frame_id][gid] = new_value;
    return true;
  }
  return false;
}

void GoodSimulator::AddGateToEventQueue(Gate& gptr) {
  for (const auto& gate : gptr.FanoutGates()) {
    event_queue_.push(gate->GetGId(), gate->GetDPI());
  }
}

bool GoodSimulator::EvalAndUpdate(Gate& gptr, int frame_id) {
  SetGateInputVal(gptr, frame_id);
  val64_t output = evaluate_.Evaluate(gptr.GetGType());
  bool is_update = UpdateGoodVal(output, gptr.GetGId(), frame_id);
  return is_update;
}

void GoodSimulator::ScanChainAnalyse() {
//  FUNC_TIME_REPORT_ICT
  ResetGoodMachine();
  std::vector<GateClassify>& marks = simdata_->gate_classify_mark_;
  marks.resize(prim_netlist_->NumGates());
  // 标记全局信号的传播路径
  std::fill(marks.begin(), marks.end(), GateClassify::UNKNOW);
  TraceAndMarkHoldSigPath(simdata_->tm_sigs_, GateClassify::TM_SIG);
  TraceAndMarkHoldSigPath(simdata_->set_reset_sigs_,
                          GateClassify::SET_RESET_SIG);
  TraceAndMarkHoldSigPath(simdata_->se_sigs_, GateClassify::SE_SIG);
  TraceAndMarkHoldSigPath(simdata_->tie_sigs_, GateClassify::TIE_SIG);
  TraceAndMarkMasterClkPath(simdata_->shift_clk_off_state_);
  TraceScanChainsSo2Si();
  // 标记shift过程中活跃的non scan dff
  MarkActiveNSC();
  // 统计上升沿下降沿数量
  StatisticChains();
  ResetGoodMachine();
}

void GoodSimulator::MarkActiveNSC() {
  auto& marks = simdata_->gate_classify_mark_;
  const auto logic_0 = val64_t(0, ~0);
  for (auto nsc_dff : prim_netlist_->GetNonScanDFFGates()) {
    Gate& nsc_gptr = *nsc_dff;
    if (nsc_gptr.FaninSize() == 0) {
      LOG_INFO("has float nsc dff {} ", nsc_gptr.GetInstName())
      continue;
    }
    if (nsc_gptr.FaninSize() != 4) {
      LOG_INFO("not float nsc dff fanin size {} is not support : {} ",
               nsc_gptr.FaninSize(), nsc_gptr.GetInstName())
      continue;
    }
    auto const data = nsc_gptr.FaninGates().data();
    const auto& set_val = good_machine_vals_[0][data[0]->GetGId()];
    const auto& reset_val = good_machine_vals_[0][data[1]->GetGId()];
    const auto& clk_val = good_machine_vals_[0][data[2]->GetGId()];
    if (set_val == logic_0 && reset_val == logic_0) {
      if (clk_val == val64_t(2, 1)) {  // 10
        marks[nsc_dff->GetGId()] = GateClassify::ACT_RISING_NSC;
      } else if (clk_val == val64_t(1, 2)) {  // 01
        marks[nsc_dff->GetGId()] = GateClassify::ACT_FALLING_NSC;
      }
    } else {
      LOG_INFO("HAS NOT ACT NSC DFF : {}", nsc_gptr.GetInstName())
    }
  }
}

void GoodSimulator::StatisticChains() {
  const auto& marks = simdata_->gate_classify_mark_;
  int count_sc_falling = 0;
  int count_sc_rising = 0;
  int active_nsc_rising_gate = 0;
  int active_nsc_falling_gate = 0;
  auto& dffs = prim_netlist_->GetDFFGates();
  for (int i = 0; i < dffs.size(); ++i) {
    auto gate_mark = marks[dffs[i]->GetGId()];
    if (gate_mark == GateClassify::FALLING_SC_DFF) {
      count_sc_falling++;
#if SIM_DEBUG_LEVEL>0
      LOG_INFO("FALLING SC DFF {}", dffs[i]->GetInstName())
#endif
    } else if (gate_mark == GateClassify::RISING_SC_DFF) {
      count_sc_rising++;
    } else if (gate_mark == GateClassify::ACT_FALLING_NSC) {
      active_nsc_falling_gate++;
#if SIM_DEBUG_LEVEL>0
      LOG_INFO("ACT_FALLING_NSC {}", dffs[i]->GetInstName())
#endif
    } else if (gate_mark == GateClassify::ACT_RISING_NSC) {
      active_nsc_rising_gate++;
#if SIM_DEBUG_LEVEL>0
      LOG_INFO("ACT_RISING_NSC {}", dffs[i]->GetInstName())
#endif
    }
  }
  simdata_->falling_edge_sc_dff_nums = count_sc_falling;
  simdata_->rising_edge_sc_dff_nums = count_sc_rising;
  simdata_->act_falling_edge_nsc_dff_nums = active_nsc_falling_gate;
  simdata_->act_rising_edge_nsc_dff_nums = active_nsc_rising_gate;
}

/**
 * @brief 通过对clk信号端口赋值 标记shift master clk 信号路径
 */
void GoodSimulator::TraceAndMarkMasterClkPath(
    const std::vector<std::pair<int, val64_t>>& sigs) {
  const auto& netlist = prim_netlist_->GetPrimNetlist();
  std::vector<GateClassify>& marks = simdata_->gate_classify_mark_;
  const val64_t clk_val_01 = val64_t(1, 2);
  const val64_t clk_val_10 = val64_t(2, 1);
  val64_t clk_val;
  // clk port input
  for (const auto& sig : sigs) {
    Gate& active_gptr = *netlist[sig.first];
    if (sig.second.ExpandByBitIdx(0) == val64_t(0, ~0)) {
      clk_val = clk_val_10;
    } else if (sig.second.ExpandByBitIdx(0) == val64_t(~0, 0)) {
      clk_val = clk_val_01;
    } else {
      LOG_ASSERT(false, "wrong pulse init val")
    }
    bool is_update = UpdateGoodVal(clk_val, sig.first, FRAME_0);
    if (is_update) {
      if (clk_val != val64_t()) {
        assert(marks[sig.first] == GateClassify::UNKNOW);
        marks[sig.first] = GateClassify::CLK_SIG;
      }
      for (const auto& fanout : active_gptr.FanoutGates()) {
        event_queue_.push(fanout->GetGId(), fanout->GetDPI());
      }
    }
  }
  while (!event_queue_.empty()) {
    int gid = event_queue_.pop();
    Gate& active_gptr = *netlist[gid];
    events_count_++;

    SetGateInputVal(active_gptr, FRAME_0);
    val64_t output = evaluate_.Evaluate(active_gptr.GetGType());
    bool is_update = UpdateGoodVal(output, gid, FRAME_0);

    if (is_update) {
      if (output == clk_val_01 || output == clk_val_10) {
        assert(marks[gid] == GateClassify::UNKNOW);
        marks[gid] = GateClassify::CLK_SIG;
      }
      for (const auto& fanout : active_gptr.FanoutGates()) {
        event_queue_.push(fanout->GetGId(), fanout->GetDPI());
      }
    }
  }
}

void GoodSimulator::TraceAndMarkHoldSigPath(
    const std::vector<std::pair<int, val64_t>>& sigs, GateClassify mark_val) {
  const auto& netlist = prim_netlist_->GetPrimNetlist();
  std::vector<GateClassify>& marks = simdata_->gate_classify_mark_;
  for (const auto& sig : sigs) {
    Gate& active_gptr = *netlist[sig.first];
    bool is_update = UpdateGoodVal(sig.second, sig.first, FRAME_0);
    if (is_update) {
      if (sig.second != val64_t()) {
        //        assert(marks[sig.first] == GateClassify::UNKNOW);
        marks[sig.first] = mark_val;
      }
      for (const auto& fanout : active_gptr.FanoutGates()) {
        event_queue_.push(fanout->GetGId(), fanout->GetDPI());
      }
    }
  }
  EventDrivenPathMark(mark_val);
}

/**
 * @brief 从so端口到si端口反响追踪扫描链 并根据遇到扫描单元dff的顺序给扫描链分段
 */
void GoodSimulator::TraceScanChainsSo2Si() {
  auto& netlist = prim_netlist_->GetPrimNetlist();
  simdata_->chain_gate_2_sigment_.resize(prim_netlist_->NumGates(),
                                         std::make_pair(-1, -1));
  auto& marks = simdata_->gate_classify_mark_;
  for (int chain_id = 0; chain_id < simdata_->chain_nums_; ++chain_id) {
    auto sc_out = prim_netlist_->ScanChains()[chain_id]->GetScanOutGate();
    int curr_seg_idx = 0;
    auto curr_gate = sc_out;
    GType logic_type;
    while (curr_gate->FaninSize() > 0) {
      auto temp_curr_gate = curr_gate;
      logic_type = curr_gate->GetGType();
      switch (logic_type) {
        case GType::G_DFF: {
          for (int sig_id = 0; sig_id < 2; ++sig_id) {
            int sig_in = curr_gate->FaninGates()[sig_id]->GetGId();
            if (good_machine_vals_[FRAME_0][sig_in] != val64_t(0, ~0)) {
              LOG_ASSERT(false, "Chain DFF SET/RESET wrong ");
            }
          }
          val64_t clk_val =
              good_machine_vals_[FRAME_0][curr_gate->FaninGates()[2]->GetGId()];
          if (clk_val == val64_t(1, 2)) {  // 01
            marks[curr_gate->GetGId()] = GateClassify::FALLING_SC_DFF;
          } else if (clk_val == val64_t(2, 1)) {  // 10
            marks[curr_gate->GetGId()] = GateClassify::RISING_SC_DFF;
          } else {
            LOG_ASSERT(false, "Chain DFF CLK wrong ");
          }
        }
          curr_gate = curr_gate->FaninGates()[3];
          break;
        case GType::G_DLAT:
          LOG_INFO("DLATCH is in Scan Chain")
          for (int sig_id = 0; sig_id < 2; ++sig_id) {
            int sig_in = curr_gate->FaninGates()[sig_id]->GetGId();
            if (good_machine_vals_[FRAME_0][sig_in] != val64_t(0, ~0)) {
              LOG_ASSERT(false, "Chain dlatch SET/RESET wrong ");
            }
          }
          marks[curr_gate->GetGId()] = GateClassify::SCAN_DLATCH;
          curr_gate = curr_gate->FaninGates()[3];
          break;
        case GType::G_MUX: {
          marks[curr_gate->GetGId()] = GateClassify::SE_MUX;
          int sel = curr_gate->FaninGates()[0]->GetGId();
          if (good_machine_vals_[FRAME_0][sel] == val64_t(0, ~0)) {
            curr_gate = curr_gate->FaninGates()[1];
          } else if (good_machine_vals_[FRAME_0][sel] == val64_t(~0, 0)) {
            curr_gate = curr_gate->FaninGates()[2];
          } else {
            LOG_ASSERT(false, "Chain MUX SEL wrong ");
          }
        } break;
        default:
          if (curr_gate->FaninSize() == 1) {
            marks[curr_gate->GetGId()] = GateClassify::ONLY_SI_PATH;
            curr_gate = curr_gate->FaninGates()[0];
          } else {
            LOG_ASSERT(false, "Chain Gate Type wrong")
          }
      }
      simdata_->chain_gate_2_sigment_[temp_curr_gate->GetGId()] =
          std::make_pair(chain_id, curr_seg_idx);
      if (logic_type == GType::G_DFF) {
        curr_seg_idx++;
      }
    }
    // SI
    if (curr_gate != nullptr && curr_gate->FaninSize() == 0) {
      marks[curr_gate->GetGId()] = GateClassify::ONLY_SI_PATH;
    } else {
      LOG_ASSERT(false, "Chain Gate Type wrong")
    }
  }
}

void GoodSimulator::EventDrivenPathMark(GateClassify marks_val) {
  auto& marks = simdata_->gate_classify_mark_;
  while (!event_queue_.empty()) {
    int gid = event_queue_.pop();
    Gate& active_gptr = *prim_netlist_->GetPrimNetlist()[gid];
    events_count_++;
    SetGateInputVal(active_gptr, FRAME_0);
    val64_t output = evaluate_.Evaluate(active_gptr.GetGType());
    bool is_update = UpdateGoodVal(output, gid, FRAME_0);
    if (is_update) {
      if (output != val64_t()) {
        //        assert(marks[gid] == GateClassify::UNKNOW);
        marks[gid] = marks_val;
      }
      for (const auto& fanout : active_gptr.FanoutGates()) {
        event_queue_.push(fanout->GetGId(), fanout->GetDPI());
      }
    }
  }
}

}  // namespace ictest
