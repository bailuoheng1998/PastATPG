// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "simulation/SAFSimulator.h"
namespace ictest {

// only frame_0 used in combined circuit
int SAFSimulator::BasicComPPFSim(Pattern& sp, std::vector<SAF*>& sa_flist) {
  LOG_ASSERT(prim_netlist_->GetNetlistType() == NetlistType::COMB_CIRCUIT,
             "circuit is not comb ,BasicComPPFSim only support comb")
  int useful_bits = sp.GetPatternCycle(0).LastParallelBitIdx() == 0
                        ? 64
                        : sp.GetPatternCycle(0).LastParallelBitIdx();
  gsim_.GoodSimulation(sp);

  int dt_count = 0;

  for (int fault_idx = 0; fault_idx < sa_flist.size(); ++fault_idx) {
    auto* fptr = sa_flist[fault_idx];
    auto fault_status = fptr->GetSAFStatus();

    if (fault_status == SAFStatus::DS || fault_status == SAFStatus::UU ||
        fault_status == SAFStatus::TI || fault_status == SAFStatus::BL ||
        fault_status == SAFStatus::DI) {
      continue;
    }

    // inject fault
    Gate& fgate = *fptr->GetSAFGate();
    auto fval = GenSAFaultValue(fptr, useful_bits);
    uint64_t inv_mask =
        ((fval ^ gsim_.good_machine_vals_[FRAME_0][fgate.GetGId()])
             .GetBit1Pos());
    fval = gsim_.good_machine_vals_[FRAME_0][fgate.GetGId()].MaskInv(inv_mask);
    bool is_diff = InitFaultyEvent(fgate, fval, FRAME_0, -1);

    if (is_diff) {
      FaultEventDriveSim(FRAME_0, -1);
      auto status = MeasurePOs(FRAME_0);
      if (status == SAFStatus::DS) {
        fptr->SetSAFStatus(SAFStatus::DS);
        dt_count++;
      } else if (status == SAFStatus::PT) {
        fptr->SetSAFStatus(SAFStatus::PT);
      }
    }
    ResetFaultMachineToGood();
  }

  return dt_count;
}

int SAFSimulator::BasicScan(Pattern& sp, std::vector<SAF*>& sa_flist) {
  gsim_.GoodSimulation(sp);
  int dt_count = 0;
  for (int fault_idx = 0; fault_idx < sa_flist.size(); ++fault_idx) {
    auto* fptr = sa_flist[fault_idx];
    auto fault_status = fptr->GetSAFStatus();

    if (fault_status == SAFStatus::DS || fault_status == SAFStatus::UU ||
        fault_status == SAFStatus::TI || fault_status == SAFStatus::BL ||
        fault_status == SAFStatus::DI) {
      continue;
    }
    auto f_status = PPSFDT(*fptr, sp);

    if (f_status == SAFStatus::DS) {
      ++dt_count;
      fptr->SetSAFStatus(SAFStatus::DS);
    } else if (f_status == SAFStatus::PT) {
      fptr->SetSAFStatus(SAFStatus::PT);
    }
    ResetFaultMachineToGood();
  }
  return dt_count;
}

SAFStatus SAFSimulator::PPSFDT(SAF& fptr, Pattern& sp) {
  auto& fgate = *fptr.GetSAFGate();
  bool is_pt = false;
  auto& pc = sp.GetPatternCycle(0);
  auto dt_status = fptr.GetSAFStatus();
  bool is_active = false;
  int useful_bits = pc.LastParallelBitIdx() == 0 ? 64 : pc.LastParallelBitIdx();
  int fgid = fgate.GetGId();
  val64_t fval = GenSAFaultValue(&fptr, useful_bits);
  has_fault_[fgid] = 1;
  for (int frame_id = 0; frame_id < MAX_FRAME_NUM; frame_id++) {
    if (frame_id == FRAME_0) {
      is_active = UpdateFaultVal(fval, fgid, frame_id);
      if (is_active) {
        AddActiveGateToEvents(fgate);
        fevent_queue_.InitCurrLevel(0);
        FaultEventDriveSim(frame_id, -1);
        dt_status = MeasurePOs(frame_id);
        if (dt_status == SAFStatus::DS) {
          break;
        } else if (dt_status == SAFStatus::PT) {
          is_pt = true;
        }
      }
    } else {
      // assum clk sig is fault free
      UpdateMemEventsBasicScan(frame_id);
      if (frame_id != 2) {
        InitFaultyEvent(fgate, fval, frame_id, -1);
        if (is_active) {
          AddActiveGateToEvents(fgate);
        }
      }

      fevent_queue_.InitCurrLevel(0);
        FaultEventDriveSim(frame_id, -1);
      if (frame_id == 2) {  // check at frame_2
          dt_status = MeasureSOs(frame_id);
        if (dt_status == SAFStatus::DS) {
          break;
        } else if (dt_status == SAFStatus::PT) {
          is_pt = true;
        }
      }
    }  // end fram_1 or fram_2
  }    // end 3-frame
  has_fault_[fgid] = 0;
  if (dt_status == SAFStatus::DS) {
    return SAFStatus::DS;
  }
  if (is_pt) {
    return SAFStatus::PT;
  }
  return SAFStatus::UO;
}

void SAFSimulator::UpdateMemEventsBasicScan(int frame_id) {
  for (int dff_id : fault_dff_ids_) {
        int master_latch_id = gsim_.slave2master_[dff_id];
        fault_path_mark_[frame_id][dff_id] = 1;
        fault_vals_[frame_id][dff_id] = fault_vals_[frame_id - 1][dff_id];
        fault_master_vals_[frame_id][master_latch_id] =
            fault_master_vals_[frame_id - 1][master_latch_id];
        auto& dff_gptr = prim_netlist_->GetPrimNetlist()[dff_id];
        fevent_queue_.push(dff_id, dff_gptr->GetDPI());
  }
    for (const auto& latch : prim_netlist_->GetDLATGates()) {
      int latch_id = latch->GetGId();
      if (!fault_path_mark_[frame_id - 1][latch_id]) {
        continue;
      }
      fault_path_mark_[frame_id][latch_id] = 1;
      fault_vals_[frame_id][latch_id] = fault_vals_[frame_id - 1][latch_id];
      //      AddActiveGateToEvents(*latch);
      fevent_queue_.push(latch->GetGId(), latch->GetDPI());
      AddActiveGateToEvents(*latch);
    }
}

SAFStatus SAFSimulator::MeasureSOs(int frame_id) {
  bool is_dt_possiable = false;
  const auto& scan_chain_mask = simdata_->gate_classify_mark_;
  for (auto dff_id : fault_dff_ids_) {
      if (scan_chain_mask[dff_id] != GateClassify::RISING_SC_DFF &&
          scan_chain_mask[dff_id] != GateClassify::FALLING_SC_DFF) {
          continue;
      }
    const Gate& dff_gptr = *prim_netlist_->GetPrimNetlist()[dff_id];
    // only for full scan no non scan cell
    if (!fault_path_mark_[frame_id][dff_id]) continue;
    val64_t faulty_val =
        fault_vals_[frame_id][dff_id];  // GetGateValue(frame_id, slave_id);
    val64_t good_val = gsim_.good_machine_vals_[frame_id][dff_id];
    val64_t res = (faulty_val ^ good_val).MaskTo0(good_val.GetBitXPos());
    if (res.GetBit1Pos()) {
      return SAFStatus::DS;
    }
    if (!is_dt_possiable) {
      is_dt_possiable = res.GetBitXPos();
    }
  }
  if (is_dt_possiable) {
    return SAFStatus::PT;
  }
  return SAFStatus::UO;
}


val64_t SAFSimulator::GetFaultVal(int frame_id, int gid) const {
  return fault_path_mark_[frame_id][gid]
             ? fault_vals_[frame_id][gid]
             : gsim_.good_machine_vals_[frame_id][gid];
}

val64_t SAFSimulator::GetFaultMasterVal(int frame_id, int dff_gid) const {
  const int master_id = gsim_.slave2master_[dff_gid];
  return fault_path_mark_[frame_id][dff_gid]
             ? fault_master_vals_[frame_id][master_id]
             : gsim_.good_master_vals_[frame_id][master_id];
}

void SAFSimulator::FaultEventDriveSim(int frame_id, int stop_flag) {
  while (!fevent_queue_.empty()) {
    int gate_id = fevent_queue_.pop();
    if (has_fault_[gate_id]) {
      // new_val= GetFaultVal(frame_id,gate_id);//only for temp
      continue;
    }
    fault_event_count_++;
    Gate& active_gptr = *prim_netlist_->GetPrimNetlist()[gate_id];
    auto logic_type = active_gptr.GetGType();

    SetInputLineVal(active_gptr, frame_id);
    auto new_val = evaluate_.Evaluate(logic_type);

    bool is_diff = SetFaultVal(active_gptr, new_val, frame_id);

    if (stop_flag == gate_id) {
      break;
    }
    if (is_diff) {
      AddActiveGateToEvents(active_gptr);
    }
  }
}

val64_t SAFSimulator::GenSAFaultValue(SAF* fptr, int useful_bits) const {
  val64_t fval;
  constexpr uint64_t all_one = ~(uint64_t(0));
  const int nouseful_bits = 64 - useful_bits;
  const uint64_t pc_all_one = all_one >> nouseful_bits;

  auto fault_type = fptr->GetSAFType();
  if (fault_type == SAFType::SA0) {
    fval = fval.MaskTo0(pc_all_one);
  } else if (fault_type == SAFType::SA1) {
    fval = fval.MaskTo1(pc_all_one);
  } else {
    std::string msg = "ERROR:not support fault type";
    LOG_ASSERT(false, msg)
  }
  return fval;
}

bool SAFSimulator::InitFaultyEvent(Gate& gptr, const val64_t& faulty_val,
                                   int frame_id, int ffr_head_id) {
  int gid = gptr.GetGId();
  bool is_diff = UpdateFaultVal(faulty_val, gid, frame_id);
  if (is_diff) {
    if (gptr.GetGType() == GType::G_PO) {
      fault_po_ids_.emplace_back(gid);
    }
    if (gid != ffr_head_id) {
      AddActiveGateToEvents(gptr);
    }
  }
  return is_diff;
}

void SAFSimulator::AddActiveGateToEvents(
    Gate& gptr) {  // may be better in event queue
  const size_t fan_out_size = gptr.FanoutGates().size();
  const auto& fanout_datas = gptr.FanoutGates().data();

  for (int i = 0; i < fan_out_size; ++i) {
    int gate_id = fanout_datas[i]->GetGId();
    int dpi = fanout_datas[i]->GetDPI();
    fevent_queue_.push(gate_id, dpi);
  }
}

void SAFSimulator::SetInputLineVal(Gate& gptr, int frame_id) {
  size_t fanin_size = gptr.FaninSize();

  auto logic_type = gptr.GetGType();
  auto* fanin_data = gptr.FaninGates().data();

  evaluate_.SetFaninSize(fanin_size);

  switch (logic_type) {
    case GType::G_DFF: {
      for (int i = 0; i < fanin_size; ++i) {
        evaluate_.SetPinValue(i,
                              GetFaultVal(frame_id, fanin_data[i]->GetGId()));
        if (i < 3 && fault_path_mark_[frame_id][fanin_data[i]->GetGId()]) {
          fault_to_clk = true;
        }
      }
      int gate_id = gptr.GetGId();
      int master_latch_id = gsim_.slave2master_[gate_id];

      if (fault_path_mark_[frame_id][gate_id]) {
        evaluate_.SetDffMasterVal(
            fault_master_vals_[frame_id][master_latch_id]);
        evaluate_.SetDffSlaveVal(fault_vals_[frame_id][gate_id]);
      } else {
        //当前逻辑门未被标记时 其初始状态应该是上一个frame的值
        //所以这里取上一个frmae的值作为其初始值
        if (frame_id > 0) frame_id--;
        evaluate_.SetDffMasterVal(
            gsim_.good_master_vals_[frame_id][master_latch_id]);
        evaluate_.SetDffSlaveVal(gsim_.good_machine_vals_[frame_id][gate_id]);
      }
    } break;
    case GType::G_DLAT: {
      for (int i = 0; i < fanin_size; ++i) {
        evaluate_.SetPinValue(i,
                              GetFaultVal(frame_id, fanin_data[i]->GetGId()));
      }
      evaluate_.SetLatchHoldVal(GetFaultVal(frame_id, gptr.GetGId()));
    } break;
    default:
      for (int i = 0; i < fanin_size; ++i) {
        evaluate_.SetPinValue(i,
                              GetFaultVal(frame_id, fanin_data[i]->GetGId()));
      }
      break;
  }
}

bool SAFSimulator::SetFaultVal(Gate& gptr, const val64_t& new_val,
                               int frame_id) {
  auto logic_type = gptr.GetGType();
  int gate_id = gptr.GetGId();
  bool is_diff = false;

  switch (logic_type) {
    case GType::G_PO:
      if (new_val != GetFaultVal(frame_id, gate_id)) {
        fault_vals_[frame_id][gate_id] = new_val;
        is_diff = true;
        if (!fault_path_mark_[frame_id][gate_id]) {
          fault_path_mark_[frame_id][gate_id] = 1;
          fault_gids_.emplace_back(gate_id);
          fault_po_ids_.emplace_back(gate_id);
        }
      }
      break;

    case GType::G_DFF: {
      int master_latch_id = gsim_.slave2master_[gate_id];

      if (!fault_path_mark_[frame_id][gate_id]) {  // path not marked
        is_diff = new_val != gsim_.good_machine_vals_[frame_id][gate_id];
        fault_master_vals_[frame_id][master_latch_id] =
            evaluate_.GetDffMasterVal();
        fault_vals_[frame_id][gate_id] = new_val;

        fault_path_mark_[frame_id][gate_id] = 1;
        fault_gids_.emplace_back(gate_id);
        fault_dff_ids_.emplace_back(gate_id);
      } else {  // path marked
        is_diff = fault_vals_[frame_id][gate_id] != new_val;
        if (is_diff) {
          fault_vals_[frame_id][gate_id] = new_val;
        }
        if (fault_master_vals_[frame_id][master_latch_id] !=
            evaluate_.GetDffMasterVal()) {
          fault_master_vals_[frame_id][master_latch_id] =
              evaluate_.GetDffMasterVal();
        }
      }
    } break;
    case GType::G_DLAT:
      is_diff = UpdateFaultVal(new_val, gate_id, frame_id);  // temp
      break;
    default:
      is_diff = UpdateFaultVal(new_val, gate_id, frame_id);
  }
  return is_diff;
}

bool SAFSimulator::UpdateFaultVal(const val64_t& new_value, int gate_id,
                                  int frame_id) {
  const auto old_val = GetFaultVal(frame_id, gate_id);
  if (new_value != old_val) {
    fault_vals_[frame_id][gate_id] = new_value;
    if (!fault_path_mark_[frame_id][gate_id]) {
      fault_path_mark_[frame_id][gate_id] = 1;
      fault_gids_.emplace_back(gate_id);
    }
    return true;
  }
  return false;
}

SAFStatus SAFSimulator::MeasurePOs(int frame_id) {
  uint64_t is_dt_possiable = 0;
  for (auto po_id : fault_po_ids_) {
    if (!fault_path_mark_[frame_id][po_id]) continue;
    auto good_val = gsim_.good_machine_vals_[frame_id][po_id];
    auto faulty_val = fault_vals_[frame_id][po_id];
    auto res = (faulty_val ^ good_val).MaskTo0(good_val.GetBitXPos());
    if (res.GetBit1Pos()) {
      return SAFStatus::DS;
    }
    if (!is_dt_possiable) {
      is_dt_possiable = res.GetBitXPos();
    }
  }
  if (is_dt_possiable) {
    return SAFStatus::PT;
  }
  return SAFStatus::UO;
}

uint64_t SAFSimulator::GetEventsCount() { return fault_event_count_; }

// init func
void SAFSimulator::AllocAllMem() {
  fault_master_vals_.resize(MAX_FRAME_NUM,
                            std::vector<val64_t>(prim_netlist_->NumDFFs()));
  fault_vals_.resize(MAX_FRAME_NUM,
                     std::vector<val64_t>(prim_netlist_->NumGates()));
  fault_path_mark_.resize(MAX_FRAME_NUM,
                          std::vector<uint8_t>(prim_netlist_->NumGates()));
  has_fault_.resize(prim_netlist_->NumGates());
  _mask.resize(prim_netlist_->NumGates());
}

void SAFSimulator::ResetFaultMachineToGood() {
  for (int frame_id = 0; frame_id < MAX_FRAME_NUM; frame_id++) {
    for (int i = 0; i < fault_gids_.size(); ++i) {
      int gate_id = fault_gids_[i];
      fault_path_mark_[frame_id][gate_id] = 0;
    }
  }
  fault_po_ids_.clear();
  fault_gids_.clear();
  fault_dff_ids_.clear();
}

}  // namespace ictest
