// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3
#include "fault/SAFList.h"

#include "util/Log.h"
#include "util/TimeUtil.h"

namespace ictest {
void SAFList::SetupNetlist(PrimNetlist* primitive_netlist,
                           CellNetlist* cell_netlist) {
  LOG_ASSERT(primitive_netlist, "ERROR: primitive netlist is nullptr");
  LOG_ASSERT(cell_netlist, "ERROR: cell netlist is nullptr");
  pnlist_ = primitive_netlist;
  cnlist_ = cell_netlist;
}

void SAFList::SetGate2SAFs() { gate2SAFs_.resize(pnlist_->NumGates()); }

void SAFList::SetSAFs() {
  int fault_id = 0;
  for (auto& group : *eqv_safs_) {
    for (int i = 0; i < group.size(); i++) {
      auto fptr = group[i];
      if (i != 0) {
        fptr->SetSAFEqv("--");

      } else {
        collapsed_safs_.emplace_back(fptr);
      }
      fptr->SetSAFId(fault_id);
      fault_id++;
      uncollapsed_safs_.emplace_back(fptr);
    }
  }
}

void SAFList::SetInternalSAFs() {
  GenerateEqvSAFs(internal_safs_);

  eqv_safs_ = &internal_eqv_safs_;
  SetSAFs();

}

void SAFList::GenerateInternalSAFs() {
  Gate* gptr = nullptr;
  SAF* fptr = nullptr;
  auto pin_name_map = cnlist_->GetPin2Names();
  // 遍历所有的cell netlist，在cell的输入和输出引脚上生成故障
  for (auto cell_ptr : cnlist_->GetCellNetlist()) {
    // 标识该cell是否是floating的（故障无法传播到po或pi的故障无法传播到该cell）
    bool floating = true;
    std::vector<int> pins = cell_ptr->GetOutputPins();
    std::vector<Gate*> pin_gates = cell_ptr->GetFanoutGates();
    if (pin_gates.size() != pins.size()) {
      LOG_ERROR("ERROR: pin size not equal: " + cell_ptr->GetInstanceName());
    }

    // 如果assert报错，则说明cell netlist有错误，应当检查vy图是否正确
    LOG_ASSERT(pin_gates.size() == pins.size(), "ERROR: pin size not equal");

    // 遍历输出引脚，生成故障
    for (int i = 0; i < pins.size(); i++) {
      // DPI >= 0代表pi能够连接到该cell
      gptr = pin_gates[i];
      int pin_id = pins[i];
      // 如果一个cell至少存在一个引脚的dpo>=0且dpi>0，则说明该cell不是floating的
      // 不在assign buf上生成故障
      // TODO: maybe dpo<0 can generate fault
      // tmax不在floating的net上生成故障
      // 但tessent生成
      if (gptr->GetDPI() >= 0 &&
          gptr->GetDPO() >= 0
          && gptr->GetGType() != GType::G_ABUF) {
        floating = false;
        LOG_ASSERT(pin_name_map.count(pin_id), "ERROR: can not find pi");
        // 使用pin path name作为故障名称
        std::string fault_name = pin_name_map[pin_id];
        // 在引脚上生成sa0和sa1故障
        GenerateInternalCellSAFs(fault_name, gptr, cell_ptr, false, i);
      }
    }
    // 如果cell的所有输出引脚都不是floating的
    // 在输入引脚上生成故障
    if (!floating) {
      pins = cell_ptr->GetInputPins();
      pin_gates = cell_ptr->GetFaninGates();
      if (pin_gates.size() != pins.size()) {
        LOG_ERROR("ERROR: pin size not equal: " + cell_ptr->GetInstanceName());
      }

      LOG_ASSERT(pin_gates.size() == pins.size(), "ERROR: pin size not equal");
      for (int i = 0; i < pins.size(); i++) {
        gptr = pin_gates[i];
        int pin_id = pins[i];

        LOG_ASSERT(pin_name_map.count(pin_id), "ERROR: can not find pi");
        // 使用pin path name作为故障名称
        std::string fault_name = pin_name_map[pin_id];
        // 在引脚上生成sa0和sa1故障
        GenerateInternalCellSAFs(fault_name, gptr, cell_ptr, true, i);
      }
    }
  }

  // 生成pi和po的故障
  // 遍历pi，生成故障
  for (auto gptr : pnlist_->GetPIGates()) {
    // 在引脚上生成sa0和sa1故障
    GenerateInternalCellSAFs(gptr->GetInstName(), gptr, nullptr, false, -1);
  }
  // 遍历po，生成故障
  for (auto gptr : pnlist_->GetPOGates()) {
    GenerateInternalCellSAFs(gptr->GetInstName(), gptr->FaninGates()[0],
                             nullptr, false, -1);
  }
}
//

void SAFList::GenerateInternalCellSAFs(const std::string& fault_name,
                                       Gate* gptr, Cell* cell_ptr,
                                       bool input_pin, int cell_pin_index) {
  SAF* fptr;
  // 生成位于引脚上的故障
  fptr = GenerateInternalCellPinSAF(fault_name, SAFType::SA0, gptr);
  if (cell_ptr != nullptr) {
    if (input_pin) {
      cell_ptr->AddSAFsAtInputPin(fptr, cell_pin_index);

    } else {
      cell_ptr->AddSAFsAtOutputPin(fptr, cell_pin_index);
    }
  }
  // 加入故障列表
  internal_safs_.emplace_back(fptr);
  name2saf_[fault_name + "sa0"] = fptr;

  fptr = GenerateInternalCellPinSAF(fault_name, SAFType::SA1, gptr);
  if (cell_ptr != nullptr) {
    if (input_pin) {
      cell_ptr->AddSAFsAtInputPin(fptr, cell_pin_index);

    } else {
      cell_ptr->AddSAFsAtOutputPin(fptr, cell_pin_index);
    }
  }
  internal_safs_.emplace_back(fptr);
  name2saf_[fault_name + "sa1"] = fptr;
}

SAF* SAFList::GenerateInternalCellPinSAF(std::string fault_name, SAFType sa,
                                         Gate* gptr) {
  auto fptr = new SAF(gptr, sa, fault_name);
  gate2SAFs_[gptr->GetGId()].emplace_back(fptr);
  return fptr;
}



void SAFList::GenerateEqvSAFs(std::vector<SAF*>& flist) {
  Gate* gptr;
  // 层级化遍历所有的gate，调用CheckEqvRelation函数处理该门上的故障的等价关系
  for (int v = pnlist_->GetLevelPrimNetlist().size() - 1; v >= 0; v--) {
    for (int idx = 0; idx < pnlist_->GetLevelPrimNetlist()[v].size(); idx++) {
      gptr = pnlist_->GetLevelPrimNetlist()[v][idx];
      LOG_ASSERT(gptr, "ERROR: gate is nullptr");
      for (int sa = 0; sa < 2; sa++) {
        CheckEqvSAFRelation(gptr, sa, false);
      }
    }
  }

  // 遍历flist，根据计算好的等价关系设置internal_eqv_safs_和saf_eqv_group_idx_
  internal_eqv_safs_.resize(eqv_group_.size(), std::vector<SAF*>());
  int fault_idx = 0;
  int sa = 0;

  for (auto fptr : flist) {
    auto gptr = fptr->GetSAFGate();
    auto ftype = fptr->GetSAFType();
    if (ftype == SAFType::SA1) {
      sa = 1;
    } else {
      sa = 0;
    }
    // 利用gptr和sa故障得出唯一的key
    int eqv_group_id = eqv_group_idx_[2 * gptr->GetGId() + sa];

    internal_eqv_safs_[eqv_group_id].emplace_back(fptr);
    fptr->SetSAFEqvGroupId(eqv_group_id);
  }
}
void SAFList::CheckEqvSAFRelation(Gate* gptr, int sa,
                                  bool optimize_gate_flag = false) {
  // 根据gate id和sa得出唯一的key，利用该key作为eqv_group_idx_的id
  int hash_key = 2 * gptr->GetGId() + sa;

  // eqv_group_idx_是hash_key到group_id的映射，代表(gptr, sa)对属于哪个等价故障组
  // eqv_group_存储了该组内所有等价的故障
  // 如果eqv_group_idx_不存在hash_key，说明(gptr,sa)对未被访问过
  if (!eqv_group_idx_.count(hash_key)) {
    // 则创建新的等价故障组
    eqv_group_.push_back(std::vector<int>());
    // 维护等价故障组id和(gptr, sa) hash_key之间的关系
    int new_group = eqv_group_.size() - 1;
    eqv_group_idx_[hash_key] = new_group;
    eqv_group_[new_group].push_back(hash_key);
  }

  int id = eqv_group_idx_[hash_key];
  int new_key;

  // 根据gptr的类型来判断输入门和输出门之间的故障等价关系
  switch (gptr->GetGType()) {
    case GType::G_BUF:
      // buf的输入sa0/sa1和输出sa0/sa1是等价故障
      new_key = 2 * gptr->FaninGates()[0]->GetGId() + sa;
      eqv_group_idx_[new_key] = id;
      eqv_group_[id].push_back(new_key);
      break;
    case GType::G_ABUF:
      // 和buf一样
      new_key = 2 * gptr->FaninGates()[0]->GetGId() + sa;
      eqv_group_idx_[new_key] = id;
      eqv_group_[id].push_back(new_key);
      break;
    case GType::G_BRH:
      break;
    case GType::G_NOT:
      // not的输入sa1/sa0和输出sa0/sa1是等价故障
      new_key = 2 * gptr->FaninGates()[0]->GetGId() + (1 - sa);
      eqv_group_idx_[new_key] = id;
      eqv_group_[id].push_back(new_key);
      break;
    case GType::G_AND:
      // AND的输出sa0和输入sa0是等价故障
      if (sa == 0) {
        for (int fi = 0; fi < gptr->FaninSize(); fi++) {
          new_key = 2 * gptr->FaninGates()[fi]->GetGId();
          eqv_group_idx_[new_key] = id;
          eqv_group_[id].push_back(new_key);
        }
      }
      break;
    case GType::G_NAND:
      // NAND的输出sa1和输入sa0是等价故障
      if (sa == 1) {
        for (int fi = 0; fi < gptr->FaninSize(); fi++) {
          new_key = 2 * gptr->FaninGates()[fi]->GetGId();
          eqv_group_idx_[new_key] = id;
          eqv_group_[id].push_back(new_key);
        }
      }
      break;
    case GType::G_OR:
      // OR输出sa1和输入sa1是等价故障
      if (sa == 1) {
        for (int fi = 0; fi < gptr->FaninSize(); fi++) {
          new_key = 2 * gptr->FaninGates()[fi]->GetGId() + 1;
          eqv_group_idx_[new_key] = id;
          eqv_group_[id].push_back(new_key);
        }
      }
      break;
    case GType::G_NOR:
      // NOR输出sa0和输入sa1是等价故障
      if (sa == 0) {
        for (int fi = 0; fi < gptr->FaninSize(); fi++) {
          new_key = 2 * gptr->FaninGates()[fi]->GetGId() + 1;
          eqv_group_idx_[new_key] = id;
          eqv_group_[id].push_back(new_key);
        }
      }
      break;
    case GType::G_PI:
    case GType::G_PO:
    case GType::G_XOR:
    case GType::G_XNOR:
    case GType::G_MUX:
    case GType::G_DFF:
    case GType::G_DLAT:
    case GType::G_TIE0:
    case GType::G_TIE1:
    case GType::G_TIEX: {
      break;
    }
    default: {
      LOG_ASSERT(false, "ERROR: not supported logic type: " +
                            std::to_string(static_cast<int>(gptr->GetGType())));
    }
  }
}


void SAFList::FindScanEnablesDISAF(std::unordered_set<int>& gid_visited_map) {
  // 用于蕴涵的queue
  std::queue<Gate*> q;
  // 从scan_enable向前蕴涵，直到传播到了类型不是NOT,INV,BRH,ABUF的门。
  // 如果gate的值为LOGIC 1，则sa0是DI，反之亦然
  for (auto& scan_enable_gptr : pnlist_->GetSEGates()) {
    pnlist_->GetGateStatus()[scan_enable_gptr->GetGId()].logic_value_ =
        pnlist_->GetSE2Values()[scan_enable_gptr->GetGId()];
    q.push(scan_enable_gptr);
    gid_visited_map.insert(scan_enable_gptr->GetGId());
  }

  // 蕴涵
  while (!q.empty()) {
    auto gptr = q.front();
    q.pop();

    LogicValue gptr_v = pnlist_->GetGateStatus()[gptr->GetGId()].logic_value_;
    assert(gptr_v == LogicValue::LOGIC_0 || gptr_v == LogicValue::LOGIC_1);
    for (auto fanout : gptr->FanoutGates()) {
      // imply one input gate
      if ((fanout->GetGType() == GType::G_ABUF ||
           fanout->GetGType() == GType::G_BUF ||
           fanout->GetGType() == GType::G_BRH ||
           fanout->GetGType() == GType::G_NOT) &&
          gid_visited_map.find(fanout->GetGId()) == gid_visited_map.end()) {
        LogicValue fanout_v = gptr_v;
        if (fanout->GetGType() == GType::G_NOT) {
          fanout_v = NOT[LogicV2IntV(gptr_v)];
        }
        pnlist_->GetGateStatus()[fanout->GetGId()].logic_value_ = fanout_v;
        q.push(fanout);
        gid_visited_map.insert(gptr->GetGId());
      }
    }

    // 如果sa故障在scan_cell的SE引脚上，我们应检查scan_cell的D引脚和SI引脚是否只通过
    // inv/buf门连接到了同一个scan_cell上
    bool have_DI_SAF = true;
    if (gptr->FanoutGates()[0]->GetGType() == GType::G_MUX &&
        pnlist_->GetGateStatus()[gptr->GetGId()].logic_value_ !=
            LogicValue::LOGIC_X) {
      auto mux_gptr = gptr->FanoutGates()[0];
      assert(gptr->FanoutSize() == 1);

      // 从scan_cell的D引脚找到只通过inv/buf门连接的DFF
      // 如果没找到，pre_scanchain_dff_ptr应为null
      Gate *pre_scanchain_dff_ptr[2] = {nullptr, nullptr},
           *fanin_d_or_si_ptr[2] = {nullptr, nullptr};
      // 获取SI和D gate，我们不关心顺序，只关心对应的pre_scanchain_dff_ptr是否相等
      for (int i = 0; i < 2; i++) {
        fanin_d_or_si_ptr[i] = mux_gptr->FaninGates()[i + 1];
      }
      // 找到只通过inv/buf/brh门连接的DFF
      for (int i = 0; i < 2; i++) {
        while (fanin_d_or_si_ptr[i] && fanin_d_or_si_ptr[i]->FaninSize() == 1 ||
               fanin_d_or_si_ptr[i]->GetGType() == GType::G_DFF) {
          if (fanin_d_or_si_ptr[i]->GetGType() == GType::G_DFF) {
            pre_scanchain_dff_ptr[i] = fanin_d_or_si_ptr[i];
            break;
          }
          fanin_d_or_si_ptr[i] = fanin_d_or_si_ptr[i]->FaninGates()[0];
        }
      }
      if (pre_scanchain_dff_ptr[0] && pre_scanchain_dff_ptr[1] &&
          pre_scanchain_dff_ptr[0] == pre_scanchain_dff_ptr[1]) {
        have_DI_SAF = false;
      }
    }
    // 设置故障为DI
    if (have_DI_SAF) {
      for (auto fptr : gate2SAFs_[gptr->GetGId()]) {
        if ((gptr_v == LogicValue::LOGIC_0 &&
             fptr->GetSAFType() == SAFType::SA1) ||
            (gptr_v == LogicValue::LOGIC_1 &&
             fptr->GetSAFType() == SAFType::SA0) &&
                fptr->GetSAFStatus() == SAFStatus::UC) {
          fptr->SetSAFStatus(SAFStatus::DI);
        }
        AddSAFToChainTestSAFs(fptr);
      }
    } else {
      // 回溯，设置所有的DI故障为NC
      auto not_have_DI_fanin_gptr = gptr->FaninGates()[0];
      while (not_have_DI_fanin_gptr &&
             not_have_DI_fanin_gptr->FanoutSize() == 1) {
        for (auto fptr : gate2SAFs_[not_have_DI_fanin_gptr->GetGId()]) {
          if (fptr->GetSAFStatus() == SAFStatus::DI) {
            fptr->SetSAFStatus(SAFStatus::UC);
            AddSAFToChainTestSAFs(fptr);
          }
        }
        not_have_DI_fanin_gptr = not_have_DI_fanin_gptr->FaninGates()[0];
      }
    }
  }
}
void SAFList::FindScanChainDISAF(std::unordered_set<int>& gid_visited_map,
                                 std::queue<Gate*>& test_se_imply_q) {
  // 设置所有的scan_cell的SI Q引脚的sa0/sa1为DI
  // 设置test_si/test_so的sa0 sa1为DI
  // 将clock/set/reset加入队列
  std::queue<Gate*> q;
  for (auto& scan_chain : pnlist_->ScanChains()) {
    for (auto& gptr : scan_chain->ScanChainGates()) {
      if (gptr->GetGType() == GType::G_DFF) {
        gid_visited_map.insert(gptr->GetGId());
        // 给logic_value_赋值为1则代表sa1故障应为DI
        // logic_x代表clock的sa0和sa1都应为DI
        pnlist_->GetGateStatus()[gptr->FaninGates()[0]->GetGId()].logic_value_ =
            LogicValue::LOGIC_1;
        // reset
        pnlist_->GetGateStatus()[gptr->FaninGates()[1]->GetGId()].logic_value_ =
            LogicValue::LOGIC_1;
        // clock
        pnlist_->GetGateStatus()[gptr->FaninGates()[2]->GetGId()].logic_value_ =
            LogicValue::LOGIC_X;
        for (int i = 0; i <= 2; i++) {
          if (gid_visited_map.find(gptr->FaninGates()[i]->GetGId()) ==
              gid_visited_map.end()) {
            gid_visited_map.insert(gptr->FaninGates()[i]->GetGId());
            q.push(gptr->FaninGates()[i]);
          }
        }
      }
      for (auto fptr : gate2SAFs_[gptr->GetGId()]) {
        if (fptr->GetSAFStatus() == SAFStatus::UC) {
          fptr->SetSAFStatus(SAFStatus::DI);
        }
      }
    }
  }

  // 从clock/set/reset向后蕴涵到pi
  // set/reset在遇到非G_BUF_ASSGIN/G_BUF/G_NOT/G_BRH门时停止蕴涵
  // 通过区分logic_value_来判断门是连接到set/reset上还是clock上
  // 连接到clock门上的门的logic_value_应为X
  while (!q.empty()) {
    auto gptr = q.front();
    q.pop();
    auto gptr_v = pnlist_->GetGateStatus()[gptr->GetGId()].logic_value_;
    assert(gptr_v == LogicValue::LOGIC_0 || gptr_v == LogicValue::LOGIC_1 ||
           gptr_v == LogicValue::LOGIC_X);

    // 检查gptr的safs是否是DI
    for (auto fptr : gate2SAFs_[gptr->GetGId()]) {
      if ((gptr_v == LogicValue::LOGIC_0 &&
           fptr->GetSAFType() == SAFType::SA0) ||
          (gptr_v == LogicValue::LOGIC_1 &&
           fptr->GetSAFType() == SAFType::SA1) ||
          gptr_v == LogicValue::LOGIC_X &&
              fptr->GetSAFStatus() == SAFStatus::UC) {
        fptr->SetSAFStatus(SAFStatus::DI);
      }
      AddSAFToChainTestSAFs(fptr);
    }

    // 把gptr的fanin加入队列
    if (gptr->GetGType() == GType::G_ABUF || gptr->GetGType() == GType::G_BUF ||
        gptr->GetGType() == GType::G_BRH || gptr->GetGType() == GType::G_NOT) {
      auto fanin_gptr = gptr->FaninGates()[0];
      if (gid_visited_map.find(fanin_gptr->GetGId()) != gid_visited_map.end()) {
        continue;
      }
      gid_visited_map.insert(fanin_gptr->GetGId());

      auto fanin_v = gptr_v;
      if (gptr->GetGType() == GType::G_NOT) {
        fanin_v = NOT[LogicV2IntV(gptr_v)];
      }
      pnlist_->GetGateStatus()[fanin_gptr->GetGId()].logic_value_ = fanin_v;
      q.push(fanin_gptr);
    }

    // set/reset在遇到G_MUX时停止蕴涵
    // 但clock继续
    if (gptr->GetGType() == GType::G_MUX &&
        pnlist_->GetGateStatus()[gptr->GetGId()].logic_value_ ==
            LogicValue::LOGIC_X) {
      auto sel_gptr = gptr->FaninGates()[0];
      // 把select门加入test_se_imply_q
      test_se_imply_q.push(sel_gptr);

      Gate* fanin_gptr = nullptr;
      SAFType sel_saf_type;
      if (pnlist_->GetGateStatus()[sel_gptr->GetGId()].ct_value_ ==
          ConstraintValue::CT_0) {
        fanin_gptr = gptr->FaninGates()[1];
        sel_saf_type = SAFType::SA1;
      } else if (pnlist_->GetGateStatus()[sel_gptr->GetGId()].ct_value_ ==
                 ConstraintValue::CT_1) {
        fanin_gptr = gptr->FaninGates()[2];
        sel_saf_type = SAFType::SA0;
      } else {
        assert(false);
      }
      for (auto fptr : gate2SAFs_[sel_gptr->GetGId()]) {
        if (fptr->GetSAFType() == sel_saf_type &&
            fptr->GetSAFStatus() == SAFStatus::UC) {
          fptr->SetSAFStatus(SAFStatus::DI);
        }
        AddSAFToChainTestSAFs(fptr);
      }
      if (fanin_gptr) {
        pnlist_->GetGateStatus()[fanin_gptr->GetGId()].logic_value_ =
            LogicValue::LOGIC_X;
        q.push(fanin_gptr);
        gid_visited_map.insert(fanin_gptr->GetGId());
      }
    }
  }
}
void SAFList::FindTestModeDISAF(std::unordered_set<int>& gid_visited_map,
                                std::queue<Gate*>& test_se_imply_q) {
  // if test_mode equal Logic_1 select test mode to run atpg
  // then test_mode sa0 fault can be classified to NC or DI
  // 如果test_mode等于LOGIC_1，选择test_mode跑atpg
  // 那么test_mode的sa0应为NC或DI
  // 如果test_mode只控制clock信号，则test_mode的sa0应为DI
  std::unordered_set<int> test_mode_visited_set;
  Gate* test_mode_gptr = nullptr;
  // 从mux sel引脚向后蕴涵到test_mode
  while (!test_se_imply_q.empty()) {
    auto gptr = test_se_imply_q.front();
    test_se_imply_q.pop();

    SAFType saf_type;
    if (pnlist_->GetGateStatus()[gptr->GetGId()].ct_value_ ==
        ConstraintValue::CT_0) {
      saf_type = SAFType::SA0;
    } else {
      saf_type = SAFType::SA1;
    }
    for (auto fptr : gate2SAFs_[gptr->GetGId()]) {
      if (fptr->GetSAFType() == saf_type &&
          fptr->GetSAFStatus() == SAFStatus::UC) {
        fptr->SetSAFStatus(SAFStatus::DI);
      }
    }

    if (gptr->GetGType() == GType::G_ABUF || gptr->GetGType() == GType::G_BUF ||
        gptr->GetGType() == GType::G_NOT || gptr->GetGType() == GType::G_BRH) {
      auto fanin_gptr = gptr->FaninGates()[0];
      if (gptr->GetGType() == GType::G_BRH) {
        LOG_ASSERT(fanin_gptr->GetGType() == GType::G_PI, "");
        //        assert();
      }
      // 如果gptr的fanin是test_mode，把gptr加入test_mode_visited_set
      if (fanin_gptr->GetGType() == GType::G_PI) {
        assert(pnlist_->GetGateStatus()[fanin_gptr->GetGId()].ct_value_ !=
               ConstraintValue::CT_FREE);
        test_mode_visited_set.insert(gptr->GetGId());
        test_mode_gptr = fanin_gptr;
        continue;
      }
      if (gid_visited_map.find(fanin_gptr->GetGId()) != gid_visited_map.end()) {
        test_se_imply_q.push(fanin_gptr);
        gid_visited_map.insert(fanin_gptr->GetGId());
      }
    } else {
      assert(false);
    }
  }

  // 如果test_mode_visited_set.size等于test_mode fanout size
  // 表示test_mode只控制clock信号
  // 把test_mode对应的sa故障设为DI
  if (test_mode_gptr != nullptr &&
      test_mode_gptr->FanoutSize() == test_mode_visited_set.size()) {
    SAFType saf_type;
    if (pnlist_->GetGateStatus()[test_mode_gptr->GetGId()].ct_value_ ==
        ConstraintValue::CT_0) {
      saf_type = SAFType::SA1;
    } else {
      saf_type = SAFType::SA0;
    }
    for (auto fptr : gate2SAFs_[test_mode_gptr->GetGId()]) {
      if (fptr->GetSAFType() == saf_type) {
        if (fptr->GetSAFStatus() == SAFStatus::UC) {
          fptr->SetSAFStatus(SAFStatus::DI);
        }
      } else {
        AddSAFToChainTestSAFs(fptr);
      }
    }
  }
}
void SAFList::FindDISAFs() {
  std::unordered_set<int> gid_visited_map;
  std::queue<Gate*> test_se_imply_q;
  FindScanEnablesDISAF(gid_visited_map);
  FindScanChainDISAF(gid_visited_map, test_se_imply_q);
  FindTestModeDISAF(gid_visited_map, test_se_imply_q);

  // 如果有一个branch是UB故障，则把scan_enable设为UNDETECTED
  for (auto& gptr : pnlist_->GetSEGates()) {
    for (auto fanout : gptr->FanoutGates()) {
      if (pnlist_->GetGateStatus()[fanout->GetGId()].is_blocked_) {
        for (auto fptr : gate2SAFs_[gptr->GetGId()]) {
          if (fptr->GetSAFStatus() == SAFStatus::DI) {
            fptr->SetSAFStatus(SAFStatus::UC);
            AddSAFToChainTestSAFs(fptr);
          }
        }
        break;
      }
    }
  }

  std::queue<Gate*> q;
  for (auto tm_gptr : pnlist_->GetTMGates()) {
    if (tm_gptr->GetGType() != GType::G_PI) {
      continue;
    }
    assert(q.size() == 0);
    q.push(tm_gptr);
    while (!q.empty()) {
      auto gptr = q.front();
      q.pop();

      for (auto fptr : gate2SAFs_[gptr->GetGId()]) {
        AddSAFToChainTestSAFs(fptr);
      }
      for (auto fanout : gptr->FanoutGates()) {
        if (fanout->GetGType() != GType::G_DFF) {
          q.push(fanout);
        }
      }
    }
  }
}


// statistic fault type count
void SAFList::SetCountSAFs() {
  count_safs_.DS_ = 0;
  count_safs_.DI_ = 0;
  count_safs_.DT_CHAIN_TEST_ = 0;
  count_safs_.PT_ = 0;
  count_safs_.PU_ = 0;
  count_safs_.AU_ = 0;
  count_safs_.ATPG_ABORT_ = 0;
  count_safs_.UC_ = 0;
  count_safs_.UO_ = 0;
  count_safs_.UU_ = 0;
  count_safs_.TI_ = 0;
  count_safs_.BL_ = 0;
  count_safs_.RE_ = 0;
  int total_fault_nums = uncollapsed_safs_.size();
  for (auto& fptr : uncollapsed_safs_) {
    auto status = fptr->GetSAFStatus();
    if (status == SAFStatus::DS) {
      ++count_safs_.DS_;
    } else if (status == SAFStatus::DI) {
      ++count_safs_.DI_;
    } else if (status == SAFStatus::DT_CHAIN_TEST) {
      ++count_safs_.DT_CHAIN_TEST_;
    } else if (status == SAFStatus::PT) {
      ++count_safs_.PT_;
    } else if (status == SAFStatus::PU) {
      ++count_safs_.PU_;
    } else if (status == SAFStatus::AU) {
      ++count_safs_.AU_;
    } else if (status == SAFStatus::ATPG_ABORT) {
      ++count_safs_.ATPG_ABORT_;
    } else if (status == SAFStatus::UC) {
      ++count_safs_.UC_;
    } else if (status == SAFStatus::UO) {
      ++count_safs_.UO_;
    } else if (status == SAFStatus::UU) {
      ++count_safs_.UU_;
    } else if (status == SAFStatus::TI) {
      ++count_safs_.TI_;
    } else if (status == SAFStatus::BL) {
      ++count_safs_.BL_;
    } else if (status == SAFStatus::RE) {
      ++count_safs_.RE_;
    } else {
      LOG_ASSERT(false, "ERROR: not support fault status ");
    }
  }
}

int32_t SAFList::NumDS() { return count_safs_.DS_; }
int32_t SAFList::NumDI() { return count_safs_.DI_; }
int32_t SAFList::NumDT_CHAIN_TEST() { return count_safs_.DT_CHAIN_TEST_; }
int32_t SAFList::NumPT() { return count_safs_.PT_; }
int32_t SAFList::NumPU() { return count_safs_.PU_; }
int32_t SAFList::NumAU() { return count_safs_.AU_; }
int32_t SAFList::NumATPG_ABORT() { return count_safs_.ATPG_ABORT_; }
int32_t SAFList::NumUC() { return count_safs_.UC_; }
int32_t SAFList::NumUO() { return count_safs_.UO_; }
int32_t SAFList::NumUU() { return count_safs_.UU_; }
int32_t SAFList::NumTI() { return count_safs_.TI_; }
int32_t SAFList::NumBL() { return count_safs_.BL_; }
int32_t SAFList::NumRE() { return count_safs_.RE_; }

void SAFList::ReportSAF() {
  LOG_INFO("{:^55}", "## Report Stuck-At Fault ##")
  LOG_INFO("{:^20}:     {:<}", "Stuck-At Fault Nums", uncollapsed_safs_.size());
}

void SAFList::GenerateExternalSAFs(const std::string& file_path) {
  std::ifstream fin(file_path);

  std::string fault_str;

  std::string fault_name;
  std::string fault_class = "";
  SAFType ftype;
  SAF* fptr;
  SAFStatus group_represent_type;

  // find  the faults that ICTest's  more than Tmax's
  int group_id = -1;
  bool have_last_fptr = true;
  while (getline(fin, fault_str)) {
    // faults[0] represent sa0 or sa1
    // faults[1] represent fault status
    // faults[2] represent fault name
    std::vector<std::string> faults;

    faults = StringUtil::SplitBySpace(fault_str);
    if (faults.size() < 3) {
      continue;
    }

    fault_name = TmaxFormat2ICTest(faults[2]);
    if (faults[0] == "sa1") {
      ftype = SAFType::SA1;
    } else {
      ftype = SAFType::SA0;
    }
    // fptr = InjectFault(fault_name, ftype);

    if (name2saf_.find(fault_name + faults[0]) != name2saf_.end()) {
      fptr = name2saf_[fault_name + faults[0]];
    } else {
      // the external fault is not in internal fault list
      // skip it
      fptr = nullptr;
      LOG_INFO("ERROR: the external fault is not in internal fault list : " +
               fault_name);
      //            assert(false);
      // because skip a fault, the eqv relation is change
      // like
      // a sa0 NC
      // e sa0 --
      // b sa0 NC
      // d sa1 --
      // if b is skip
      // then should new a group eqv faults
      if (faults[1] != "--") {
        group_id++;
        fault_class = "NC";
        external_eqv_safs_.push_back(std::vector<SAF*>());
      }
      have_last_fptr = false;
      continue;
    }
    if (faults[1] != "--") {
      if (faults[1] == "DI") {
        group_represent_type = SAFStatus::DI;
      } else {
        group_represent_type = SAFStatus::UC;
      }
      fault_class = "NC";
      group_id++;
      external_eqv_safs_.push_back(std::vector<SAF*>());
    } else {
      if (have_last_fptr || fault_class != "NC") {
        fault_class = "--";
      }
    }
    // store fault in external_flist
    fptr->SetSAFStatus(group_represent_type);
    fptr->SetSAFEqv(fault_class);
    external_safs_.emplace_back(fptr);
    external_eqv_safs_[group_id].push_back(fptr);
    fptr->SetSAFEqvGroupId(group_id);
    have_last_fptr = true;
  }
}
std::string SAFList::TmaxFormat2ICTest(const std::string& fault_name) {
  // replace / with ==

  int count = 0, len = fault_name.size();
  for (char c : fault_name) {
    if (c == '/') {
      count++;
    }
  }
  std::string new_fname;
  if (count > 0) {
    new_fname.resize(len + count + 2);
    new_fname[0] = '=';
    new_fname[1] = '=';
    for (int i = 0, j = 2; i < fault_name.size() && j < new_fname.size();
         i++, j++) {
      if (fault_name[i] != '/') {
        new_fname[j] = fault_name[i];
      } else {
        new_fname[j] = '=';
        j++;
        new_fname[j] = '=';
      }
    }
  } else {
    new_fname = fault_name;
  }
  return new_fname;
}

std::string SAFList::ICTestFormat2TMax(const std::string& fault_name) {
  std::string res;
  int begin;
  int len = 0;
  if (fault_name[0] != '=') {
    res = fault_name;
  } else {
    for (int i = 0; i < fault_name.size(); i++) {
      if (fault_name[i] != '=') {
        if (len == 0) {
          begin = i;
        }
        len++;
      } else {
        if (len != 0) {
          if (res.size() != 0) {
            res += "/";
          }
          res += fault_name.substr(begin, len);
          len = 0;
        }
      }
    }
    if (len != 0) {
      if (res.size() != 0) {
        res += "/";
      }
      res += fault_name.substr(begin, len);
      len = 0;
    }
  }
  return res;
}

void SAFList::SetExternalSAFs() {
  //  GenerateDomSAFTree();
  eqv_safs_ = &external_eqv_safs_;
  SetSAFs();
  //  sort(uncollapsed_SAFs_.begin(), uncollapsed_SAFs_.end(), SAF_LEVEL_CMP());
}



}  // namespace ictest