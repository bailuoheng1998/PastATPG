// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "netlist/PrimNetlist.h"
#include <fstream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include "common/Define.h"
#include "util/Log.h"

namespace ictest {

void PrimNetlist::ClearPrimNetlist() {
  for (auto& gate : prim_netlist_) {
    if (gate != nullptr) {
      delete gate;
      gate = nullptr;
    }
  }
  prim_netlist_.clear();
  gate_status_.clear();
  pi_gates_.clear();
  po_gates_.clear();
  ppi_gates_.clear();
  ppo_gates_.clear();
  tie_gates_.clear();
  stem_gates_.clear();
  dlat_gates_.clear();
  dff_gates_.clear();
  non_scan_dffs_.clear();
  clk_gates_.clear();
  scan_enable_gates_.clear();
  test_mode_gates_.clear();
  constraint_gates_.clear();
  level_prim_netlist_.clear();
  scan_chains_.clear();

  net2names_.clear();
  pin2names_.clear();
  pin2branch_.clear();
  netname2pipogptr_.clear();
  clk_off_state_.clear();
  se2values_.clear();
  tm2values_.clear();
  ct2values_.clear();
}

void PrimNetlist::InitPrimNetlist() {
  top_module_name_ = "";
  num_gates_ = 0;
  num_pis_ = 0;
  num_pos_ = 0;
  max_level_ = 0;
}

void PrimNetlist::ResizePrimNetlist() {
  prim_netlist_.resize(--num_gates_);
  prim_netlist_.shrink_to_fit();
}

bool PrimNetlist::ReadPrimVY(const std::string& vy_file) {
  std::ifstream fin(vy_file);
  if (!fin.is_open()) {
    LOG_ERROR("ERROR:cannot open primitive vy file:" + vy_file);
    return false;
  }
  std::stringstream is(vy_file);
  std::string tmp;
  std::vector<std::string> tmp_vec;
  while (getline(is, tmp, '/'))
    ;
  //  int pos_ = tmp.find('_');
  // todo write topmodule name in vy
  size_t pos_ = tmp.find_last_of('_');
  LOG_ASSERT(pos_ < tmp.size(), "find top module name err ")
  top_module_name_ = tmp.substr(0, pos_);
  //  LOG_INFO("gaojun top_module_name_ : " + top_module_name_);

  std::vector<std::string> vstr;
  std::string line;
  while (fin.good()) {
    getline(fin, line);
    vstr.push_back(line);
  }
  // current gate's output net id
  int gate_idx = 0;
  Gate* gate = nullptr;

  // initial pi, inout, po
  num_gates_ = stoi(vstr[0]);
  prim_netlist_.resize(3 * num_gates_);
  num_pis_ = stoi(vstr[1]);
  num_pos_ = stoi(vstr[5]);
  std::istringstream pi_str(vstr[2]), po_str(vstr[6]);
  while (pi_str >> gate_idx) {
    //    gate = new Gate(gate_idx, GType::G_PI, "");
    gate = new Gate(GType::G_PI, "");
    gate->SetGId(gate_idx);
    prim_netlist_[gate_idx] = gate;
    pi_gates_.push_back(gate);
  }
  while (po_str >> gate_idx) {
    // po gate's logic id and name is defined in the following code
    if (prim_netlist_[gate_idx] != nullptr) {
      po_gates_.push_back(prim_netlist_[gate_idx]);
    } else {
      gate = new Gate(GType::G_UNKNOWN, "");
      gate->SetGId(gate_idx);
      prim_netlist_[gate_idx] = gate;
      po_gates_.push_back(gate);
    }
  }

  std::string str;
  int count_flag = 0;
  // initial other logic gates
  for (int i = 7; i < vstr.size(); i++) {
    std::vector<std::string> vtmp;
    std::istringstream ifline(vstr[i]);
    while (ifline >> str) {
      vtmp.push_back(str);
    }
//    std::cout<<vtmp[0]<<std::endl;
    if (vtmp.size() > 0 && vtmp.size() == 2) {
      int id = stoi(vtmp[0]);
      if (-1 == id) {
//        return false;
        count_flag++;
      }
      if (1 == count_flag) {
        int net_id = id;
        net2names_[net_id] = vtmp[1];
      } else if (2 == count_flag) {
        int pin_id = id;
        pin2names_[pin_id] = vtmp[1];
      }
//      std::cout<<name<<std::endl;
    }
    else if (vtmp.size() > 2) {

      int logic_id = stoi(vtmp[0]);
      std::string name = vtmp[1];

      gate_idx = stoi(vtmp[3]);
      int succ_pin_id = stoi(vtmp[4]);
      int num_inout = stoi(vtmp[5]);
      int num_in = stoi(vtmp[6]);
//      std::cout<<name<< 122 <<std::endl;
      int index = 6;
      // po gate
      if (prim_netlist_[gate_idx] != nullptr) {
        Gate* po = prim_netlist_[gate_idx];
        po->SetGType(static_cast<GType>(logic_id));
        po->SetInstName(name);
        po->SuccPinId().push_back(succ_pin_id);

        for (size_t j = 0; j < num_in; j++) {
          int pred_net_id = stoi(vtmp[++index]);
          int pred_pin_id = stoi(vtmp[++index]);
          po->AddToPredNetId(pred_net_id);
          po->PredPinId().push_back(pred_pin_id);
        }

      } else {
        // other gate
        // gate = new Gate(logic_id, name);
        gate = new Gate(static_cast<GType>(logic_id), name);
        gate->SetGId(gate_idx);
        gate->SuccPinId().push_back(succ_pin_id);

        for (size_t j = 0; j < num_in; j++) {
          int pred_net_id = stoi(vtmp[++index]);
          int pred_pin_id = stoi(vtmp[++index]);
          gate->AddToPredNetId(pred_net_id);
          gate->PredPinId().push_back(pred_pin_id);
        }

        if (gate->GetGType() == GType::G_DFF) {
          dff_gates_.push_back(gate);
        }
        if (gate->GetGType() == GType::G_DLAT) {
          dlat_gates_.push_back(gate);
        }
        prim_netlist_[gate_idx] = gate;
      }
    }

  }
  // set pi name
  for (int i = 0; i < pi_gates_.size(); i++) {
    Gate* gptr = pi_gates_[i];
    LOG_ASSERT(gptr, "gate is nullptr");
    int pi_idx = gptr->GetGId();
    gptr->SetInstName(net2names_[pi_idx]);
    netname2pipogptr_[gptr->GetInstName()] = gptr;
  }
  fin.close();
  return true;
}

void PrimNetlist::ConnectPrimNetlist() {
  std::vector<bool> visited(num_gates_, false);
  std::map<int, bool> ste_in;
  bool supply0_flag = false, supply1_flag = false, supplyX_flag = false;

  for (int i = 0; i < num_pos_; i++) {
    Gate* gate = po_gates_[i];
    LOG_ASSERT(gate, "gate is nullptr");
    std::queue<Gate*> gate_que;
    if (!visited[gate->GetGId()]) {
      gate_que.push(gate);
      visited[gate->GetGId()] = true;
    }

    while (!gate_que.empty()) {
      Gate* tmp = gate_que.front();
      gate_que.pop();
      std::vector<int>& pred = tmp->PredNetId();
      if (tmp->GetGType() == GType::G_DLAT) {
        //        LOG_INFO(tmp->GetInstName());
        //        LOG_INFO(pred.size());
      }
      for (int j = 0; j < pred.size(); j++) {
        // processing supply0 & supply1
        if (net2names_[pred[j]] == "supply0" && !supply0_flag) {
          Gate* supply0_gate = new Gate(GType::G_TIE0, "supply0");
          supply0_gate->SetGId(pred[j]);
          //          supply0_gate->SetValue(LOGIC_0);  // atpg value
          prim_netlist_[pred[j]] = supply0_gate;
          tie_gates_.push_back(supply0_gate);
          // SetGlobalGate(supply0_gate->GetId());
          // SetSpecialGate(supply0_gate->GetGId(), SUPPLY_GATE);  // no use
          supply0_flag = true;
        }
        if (net2names_[pred[j]] == "supply1" && !supply1_flag) {
          Gate* supply1_gate = new Gate(GType::G_TIE1, "supply1");
          supply1_gate->SetGId(pred[j]);
          // supply1_gate->SetValue(LOGIC_1);  // atpg value
          prim_netlist_[pred[j]] = supply1_gate;
          tie_gates_.push_back(supply1_gate);
          // SetGlobalGate(supply1_gate->GetId());
          // SetSpecialGate(supply1_gate->GetGId(), SUPPLY_GATE);  // no use
          supply1_flag = true;
        }
        if (net2names_[pred[j]] == "supplyX" && !supplyX_flag) {
          Gate* supplyX_gate = new Gate(GType::G_TIEX, "supplyX");
          supplyX_gate->SetGId(pred[j]);
          //          supplyX_gate->SetValue(LOGIC_x);  // atpg value
          prim_netlist_[pred[j]] = supplyX_gate;
          tie_gates_.push_back(supplyX_gate);
          // SetGlobalGate(supply1_gate->GetId());
          // SetSpecialGate(supplyX_gate->GetGId(), SUPPLY_GATE);  // no use
          supplyX_flag = true;
        }
        // connection fanin & fanouts
        Gate* fanin = prim_netlist_[pred[j]];
        tmp->AddToFaninGates(fanin);
        fanin->AddToFanoutGates(tmp);
        if (!visited[fanin->GetGId()]) {
          gate_que.push(fanin);
          visited[fanin->GetGId()] = true;
        }
      }
    }
  }
}

void PrimNetlist::AddPOGates() {
  Gate* gptr = nullptr;
  for (int i = 0; i < po_gates_.size(); i++) {
    gptr = po_gates_[i];
    LOG_ASSERT(gptr, "gate is nullptr");
    num_gates_ += 1;
    int id = num_gates_ - 2;
    int po_idx = gptr->GetGId();
    Gate* po = new Gate(GType::G_PO, net2names_[po_idx]);
    po->SetGId(id);
    prim_netlist_[id] = po;
    gptr->AddToFanoutGates(po);
    po->AddToFaninGates(gptr);
    po->AddToPredNetId(gptr->GetGId());
    //    po->AddToPredNetId().push_back(gptr->GetGId());
    po_gates_[i] = po;
    netname2pipogptr_[po->GetInstName()] = po->FaninGates()[0];
  }
}

void PrimNetlist::AddBranchGates() {
  Gate* stem = nullptr;
  Gate* branch = nullptr;
  Gate* new_gate = nullptr;
  stem_gates_.clear();
  for (int i = 0; i < num_gates_; i++) {
    Gate* stem = prim_netlist_[i];
    if (stem != nullptr) {
      if (stem->FanoutSize() > 1) {
        stem_gates_.push_back(stem);
      }
    }
  }
  for (int m = 0; m < stem_gates_.size(); m++) {
    stem = stem_gates_[m];
    LOG_ASSERT(stem, "gate is nullptr");
    for (int i = 0; i < stem->FanoutSize(); i++) {
      branch = stem->FanoutGates()[i];
      num_gates_ += 1;
      int id = num_gates_ - 2;
      new_gate = new Gate(GType::G_BRH, "STEM_" + std::to_string(id));
      new_gate->SetGId(id);
      prim_netlist_[id] = new_gate;
      new_gate->AddToFaninGates(stem);
      new_gate->AddToFanoutGates(branch);
      new_gate->AddToPredNetId(stem->GetGId());
      stem->FanoutGates()[i] = new_gate;

      for (int j = 0; j < branch->FaninSize(); j++) {
        Gate* cur_fanin = branch->FaninGates()[j];
        if (cur_fanin->GetGId() == stem->GetGId()) {
          branch->PredNetId()[j] = new_gate->GetGId();
          branch->FaninGates()[j] = new_gate;
          branch->PredNetId()[j] = new_gate->GetGId();

          if (branch->PredPinId().size() > 0)
            //            new_gate->AddSuccPin(branch->g_pred_pin_id_[j]);
            new_gate->SuccPinId().push_back(branch->PredPinId()[j]);
          break;
        }
      }
    }
  }
}

void PrimNetlist::Levelize() {
  SetDPI();
  SetLevelPrimNetlist();
  SetDPO();
}
///@todo need detected feed back combine loop
void PrimNetlist::SetDPI() {
  // topo sort
  std::queue<Gate*> gq;
  std::vector<int> fanin_count(num_gates_, 0);
  std::vector<uint8_t> visited(num_gates_, 0);

  for (auto& pi_gate : pi_gates_) {
    gq.push(pi_gate);
    visited[pi_gate->GetGId()] = 1;
  }
  for (auto& tie_gate : tie_gates_) {
    gq.push(tie_gate);
    visited[tie_gate->GetGId()] = 1;
  }
  for (auto& dff_gate : dff_gates_) {
    gq.push(dff_gate);
    visited[dff_gate->GetGId()] = 1;
  }
  gq.push(nullptr);  // next level sign
  int curr_level = 0;
  while (!gq.empty()) {
    Gate* gptr = gq.front();
    gq.pop();
    if (gptr == nullptr) {
      if (!gq.empty()) {
        curr_level++;
        gq.push(nullptr);
      }
      continue;
    }
    gptr->SetDPI(curr_level);
    for (auto out : gptr->FanoutGates()) {
      if (visited[out->GetGId()]) {
        continue;
      }
      fanin_count[out->GetGId()]++;
      if (fanin_count[out->GetGId()] == out->FaninSize()) {
        gq.push(out);
        visited[out->GetGId()] = 1;
      }
    }
  }
  max_level_ = curr_level + 1;
}

void PrimNetlist::SetDPO() {
  // from po to pi
  Gate* gptr = nullptr;
  int depth = 0;

  for (int i = 0; i < po_gates_.size(); i++) {
    gptr = po_gates_[i];
    LOG_ASSERT(gptr, "gptr is nullptr");
    gptr->SetDPO(0);
  }
  // TODO: support non-scan-cell dff
  // if dff is non scan cell, maybe the dff pin Q can not propagate to po
  // the DFF is unused the dff's input gate also unused
  for (auto dff : dff_gates_) {
    for (auto fanin : dff->FaninGates()) {
      fanin->SetDPO(0);
    }
  }

  for (int i = max_level_; i >= 0; i--) {
    for (int j = 0; j < level_prim_netlist_[i].size(); j++) {
      gptr = level_prim_netlist_[i][j];
      LOG_ASSERT(gptr, "gptr is nullptr");
      if (gptr->FanoutSize() == 0) {
        continue;
      }
      int count = 0;
      int fanout_count = 0;
      depth = INT_MAX;
      for (auto fanout : gptr->FanoutGates()) {
        if (fanout->GetGType() != GType::G_DFF) {
          fanout_count++;
          if (fanout->GetDPI() >= 0) {
            count++;
            depth = std::min(depth, fanout->GetDPO());
          }
        }
      }
      if (fanout_count == 0) {
        depth = 0;
      } else {
        if (count > 0) {
          depth += 1;
        } else {
          depth = -1;
        }
      }
      gptr->SetDPO(depth);
    }
  }
}

void PrimNetlist::SetLevelPrimNetlist() {
  // initialize _levels
  level_prim_netlist_.resize(max_level_ + 1);
  for (int i = 0; i < num_gates_; i++) {
    Gate* gptr = prim_netlist_[i];
    if (gptr != nullptr && gptr->GetDPI() >= 0) {
      level_prim_netlist_[gptr->GetDPI()].push_back(gptr);
    }
  }
}

void PrimNetlist::AddPin2BranchGates() {
  pin2branch_.resize(pin2names_.size());
  for (int i = 0; i < stem_gates_.size(); ++i) {
    Gate* tmp_gate = stem_gates_[i];
    for (int j = 0; j < tmp_gate->FanoutSize(); ++j) {
      Gate* branch_gate = tmp_gate->FanoutGates()[j];
      if (branch_gate->SuccPinId().size() == 1) {
        int pin_id=branch_gate->SuccPinId()[0];
        pin2branch_[pin_id] = branch_gate;
      }
    }
  }
}

// todo
bool PrimNetlist::CheckPrimNetlistValidation() {
  bool exist_warning = false;

  for (auto& gate : pi_gates_) {
    if (gate == nullptr) {
      LOG_ASSERT(gate, "ERROR: exist nullptr in pi gate");
    }
  }
  for (auto& gate : po_gates_) {
    if (gate == nullptr) {
      LOG_ASSERT(gate, "ERROR: exist nullptr in po gate");
    }
  }
  for(uint32_t gate_id=0;gate_id<prim_netlist_.size();++gate_id){
    auto& gate=prim_netlist_[gate_id];
        if (gate == nullptr) {
      LOG_ERROR("ERROR: exist nullptr in primitive netlist gate idx {} ",gate_id);
    }
  }
  return exist_warning;
}

void PrimNetlist::ReportPI() {
  // 输入pi数量
  LOG_INFO("{:^55}", "## Report PI ##")
  LOG_INFO("{:^20}:     {:<}", "PI Nums", pi_gates_.size());
}

void PrimNetlist::ReportPO() {
  // 输出po数量
  LOG_INFO("{:^55}", "## Report PO ##")
  LOG_INFO("{:^20}:     {:<}", "PO Nums", po_gates_.size());
}

void PrimNetlist::ReportCLK() {
  LOG_INFO("{:^55}", "## Report CLK ##")
  // 打印时钟信号的id，名称和offstate值
  for (auto clk : clk_gates_) {
    auto clk_id = clk->GetGId();
    auto clk_name = clk->GetInstName();
    auto clk_offstate = clk_off_state_[clk_id];
    // 将offstate值从枚举类转换为数值
    int offstate = 0;
    if (clk_offstate == LogicValue::LOGIC_0) {
      offstate = 0;
    } else if (clk_offstate == LogicValue::LOGIC_1) {
      offstate = 1;
    }
    LOG_INFO("{:^20}:     {:<}", "CLK ID", clk_id);
    LOG_INFO("{:^20}:     {:<}", "CLK Name", clk_name);
    LOG_INFO("{:^20}:     {:<}", "CLK Offstate", offstate);
  }
}

void PrimNetlist::ReportScanChain() {
  LOG_INFO("{:^55}", "## Report ScanChain ##");
  LOG_INFO("{:^20}:     {:<}", "Netlist Type", NetlistType2Str());
  LOG_INFO("{:^20}:     {:<}", "Chain Nums", scan_chains_.size());
  for (auto sc : scan_chains_) {
    auto chain_name = sc->GetChainName();
    auto chain_len = sc->GetChainLength();
    auto si_name = sc->GetScanInGate()->GetInstName();
    auto so_name = sc->GetScanOutGate()->GetInstName();
    LOG_INFO("{:^20}:     {:<}", "Chain Name", chain_name);
    LOG_INFO("{:^20}:     {:<}", "Chain Length", chain_len);
    LOG_INFO("{:^20}:     {:<}", "SI Name", si_name);
    LOG_INFO("{:^20}:     {:<}", "SO Name", so_name);
  }
}

void PrimNetlist::ReportNetlist() {
  LOG_INFO("{:^55}", "## Report Netlist ##");
  LOG_INFO("{:^20}:     {:<}", "Top Module", top_module_name_);
  //  LOG_INFO("{:^20}:     {:<}", "Netlist Type", NetlistType2Str());
  LOG_INFO("{:^20}:     {:<}", "Gate Nums", num_gates_);
  LOG_INFO("{:^20}:     {:<}", "PI Nums", num_pis_);
  LOG_INFO("{:^20}:     {:<}", "PO Nums", num_pos_);
  LOG_INFO("{:^20}:     {:<}", "Tie Nums", tie_gates_.size());
  LOG_INFO("{:^20}:     {:<}", "Latch Nums", dlat_gates_.size());
  LOG_INFO("{:^20}:     {:<}", "All DFF Nums", dff_gates_.size());
  //  LOG_INFO("{:^20}:     {:<}", "Scan Chain nums", scan_chains_.size());
  //  LOG_INFO("{:^20}:     {:<}", "Non-Scan DFF nums", non_scan_dffs_.size());
  LOG_INFO("{:^20}:     {:<}", "Max DPI Level", max_level_);
}

// output Gate status to file_path for comparison
void PrimNetlist::DumpPrimNetlist(const std::string& file_path) {
  std::ofstream ofs(file_path, std::ios::out);
  ofs << "Netlist "
         "Data============================================================"
      << std::endl;
  ofs << "\tTop module:   " << top_module_name_ << std::endl;
  ofs << "\tNetlist type: " << NetlistType2Str() << std::endl;
  ofs << "\tGate nums:    " << num_gates_ << std::endl;
  ofs << "\tMax level:    " << max_level_ << std::endl;
  int max_width = 0;
  GetMaxWidthNums(max_width);
  ofs << "\tMax width:    " << max_width << std::endl;
  ofs << "\tPI nums:      " << num_pis_ << std::endl;
  ofs << "\tPO nums:      " << num_pos_ << std::endl;
  ofs << "\tPPI nums:     " << ppi_gates_.size() << std::endl;
  ofs << "\tPPO nums:     " << ppo_gates_.size() << std::endl;
  ofs << "\tTie nums:     " << tie_gates_.size() << std::endl;
  ofs << "\tLatch nums:   " << dlat_gates_.size() << std::endl;
  ofs << "\tDFF nums:     " << dff_gates_.size() << std::endl;
  ofs << "\tNon Scan DFF nums:" << non_scan_dffs_.size() << std::endl;
  ofs << "\tScan Chain nums:  " << scan_chains_.size() << std::endl;

  ofs << "\tRising edge DFF:  " << has_rising_edge_dffs_ << std::endl;
  ofs << "\tFalling edge DFF: " << has_falling_edge_dffs_ << std::endl;
  ofs << "====================================================================="
         "==="
      << std::endl;

  ofs << "" << std::endl;
  ofs << "Detailed "
         "Data==========================================================="
      << std::endl;

  for (int i = 0; i < GetPrimNetlist().size() && GetPrimNetlist()[i] != nullptr;
       ++i) {
    Gate* gate = GetPrimNetlist()[i];
    LOG_ASSERT(gate, "gate is nullptr");
    ofs << "GateId:\t" << gate->GetGId() << " " << std::endl;
    ofs << "\tLogicType:\t" << (int)gate->GetGType() << " " << std::endl;
    ofs << "\tName:\t" << gate->GetInstName() << " " << std::endl;
    ofs << "\tDpi:\t" << gate->GetDPI() << " " << std::endl;
    ofs << "\tDpo:\t" << gate->GetDPO() << " " << std::endl;
    ofs << "\tis cell border\t" << gate->IsCellBorder() << " " << std::endl;
    ofs << "\tis fb loop\t" << gate->IsFBLoop() << " " << std::endl;
    ofs << "\tdff master id\t" << gate->GetDffMasterId() << " " << std::endl;

    int output_size = gate->FanoutSize();
    ofs << "\toutput_size:\t" << output_size << " " << std::endl;
    ofs << "\t\toutput:\t\t";
    for (int i = 0; i < output_size; i++) {
      ofs << gate->FanoutGates()[i]->GetGId() << " ";
    }
    ofs << std::endl;
    int input_size = gate->FaninSize();
    ofs << "\tinput_size:\t" << input_size << " " << std::endl;
    ofs << "\t\tinput:\t\t";
    for (int i = 0; i < input_size; i++) {
      ofs << gate->FaninGates()[i]->GetGId() << " ";
    }
    ofs << std::endl;
    ofs << "-------------------------------------------------------------------"
           "-----"
        << std::endl;
  }

  ofs << std::endl;

  ofs.close();
}

std::string PrimNetlist::NetlistType2Str() {
  switch (netlist_type_) {
    case NetlistType::COMB_CIRCUIT:
      return "Comb Circuit";
    case NetlistType::FULL_SCAN_CIRCUIT:
      return "Full Scan Circuit";
    case NetlistType::PART_SCAN_CIRCUIT:
      return "Partial Scan Circuit";
    case NetlistType::SEQU_CIRCUIT:
      return "Seq Circuit";
    default:
      return "Unknown Circuit";
  }
}

void PrimNetlist::GetMaxWidthNums(int& max_width) {
  max_width = 0;
  for (int i = 0; i <= max_level_; i++) {
    int width_num = level_prim_netlist_[i].size();
    if (width_num > max_width) {
      max_width = width_num;
    }
  }
}

void PrimNetlist::GetMaxFaninNums(int& max_fanin_nums) {
  num_gates_--;

  max_fanin_nums = 0;
  for (int i = 0; i < num_gates_; i++) {
    Gate& gptr = *prim_netlist_[i];
    int fanin_nums = gptr.FaninSize();
    if (fanin_nums > max_fanin_nums) {
      max_fanin_nums = fanin_nums;
    }
  }
}

void PrimNetlist::GetBufferNums(int& total_buffer, int& buffer_cas) {
  total_buffer = 0;
  buffer_cas = 0;
  for (int i = 0; i < num_gates_; i++) {
    Gate& gptr = *prim_netlist_[i];
    if (gptr.GetGType() == GType::G_BUF) {
      total_buffer++;
      if (gptr.FaninSize() > 0 && gptr.FanoutSize() > 0) {
        //        int fanin_id = gptr.GetFaninId(0);
        int fanin_id = gptr.FaninGates()[0]->GetGId();
        //        int fanout_id = gptr.GetFanoutId(0);
        int fanout_id = gptr.FanoutGates()[0]->GetGId();

        Gate& fanin = *prim_netlist_[fanin_id];
        Gate& fanout = *prim_netlist_[fanout_id];
        if (fanin.GetGType() != GType::G_BUF &&
            fanout.GetGType() == GType::G_BUF) {
          buffer_cas++;
        }
      }
    }
  }
}

void PrimNetlist::GetInverterNums(int& total_inv, int& inv_cas) {
  total_inv = 0;
  inv_cas = 0;
  for (int i = 0; i < num_gates_; i++) {
    Gate& gptr = *prim_netlist_[i];
    if (gptr.GetGType() == GType::G_NOT) {
      total_inv++;
      if (gptr.FanoutSize() > 0) {
        int fanout_id = gptr.FanoutGates()[0]->GetGId();
        Gate& fanout = *prim_netlist_[fanout_id];
        if (fanout.GetGType() != GType::G_NOT) {
          inv_cas++;
        }
      }
    }
  }
}

void PrimNetlist::GetStemNums(int& stem_nums) {
  stem_nums = 0;
  for (int i = 0; i < num_gates_; i++) {
    Gate& gptr = *prim_netlist_[i];
    if (gptr.GetGType() == GType::G_BRH) {
      stem_nums++;
    }
  }
}


void PrimNetlist::AnalyseNetlistType() {
  if (GetDFFGates().size() == 0) {
    SetNetlistType(NetlistType::COMB_CIRCUIT);
  } else if (GetNonScanDFFGates().size() == 0 && ScanChains().size() > 0) {
    SetNetlistType(NetlistType::FULL_SCAN_CIRCUIT);
  } else if (GetNonScanDFFGates().size() > 0 && ScanChains().size() > 0) {
    SetNetlistType(NetlistType::PART_SCAN_CIRCUIT);
  } else if (GetNonScanDFFGates().size() > 0 && ScanChains().size() == 0) {
    SetNetlistType(NetlistType::SEQU_CIRCUIT);
  } else {
    SetNetlistType(NetlistType::UNKNOWN_CIRCUIT);
  }
}

}  // namespace ictest
