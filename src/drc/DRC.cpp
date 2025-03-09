// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3
#include "drc/DRC.h"

#include <queue>
#include <unordered_set>
#include <vector>

#include "util/Log.h"

namespace ictest {

void DRC::RunDRC() {
  IdentifyPIs();
  IdentifyPOs();
  IdentifyCLKs();
  IdentifySE();
  IdentifyTM();
  IdentifyCTs();
  IdentifyScanChains();
  IdentifyNonScanDFFs();
  IdentifyPPIs();
  IdentifyPPOs();
  Markchaingates();
}

void DRC::SortByNames(std::vector<std::string>& names,
                      std::vector<Gate*>& gates) {
  LOG_ASSERT(names.size() == gates.size(),
             "ERROR: pin names size not match gates size");
  std::unordered_map<std::string, Gate*> index_map;
  // build name to gptr map
  for (auto& gptr : gates) {
    LOG_ASSERT(gptr, "ERROR: gate a is nullptr");
    index_map[gptr->GetInstName()] = gptr;
  }
  // swap gptr order
  for (int i = 0; i < names.size(); i++) {
    std::string& name = names[i];
    Gate* gptr = index_map[name];
    if (gptr != gates[i]) {
      gates[i] = gptr;
    }
  }
}

void DRC::IdentifyPIs() {
  SortByNames(stil_parser_->PINames(), pnlist_->GetPIGates());
}

void DRC::IdentifyPOs() {
  SortByNames(stil_parser_->PONames(), pnlist_->GetPOGates());
}

void DRC::IdentifyPPIs() {
    for (int i = 0; i < pnlist_->ScanChains().size(); ++i) {
        std::vector<Gate*> single_ppi;
        for (int j = 0; j < pnlist_->ScanChains()[i]->ScanChainDFFs().size(); ++j) {
            Gate* dff = pnlist_->ScanChains()[i]->ScanChainDFFs()[j];
            assert(dff->FanoutSize() == 1 || dff->FanoutSize() == 3);
            if (dff->FanoutSize() == 1) {
                dff->FanoutGates()[0]->if_ppi_ = true;
                single_ppi.emplace_back(dff->FanoutGates()[0]);
            } else if (dff->FanoutSize() == 3) {
                for (int k = 0; k < dff->FanoutGates().size(); ++k) {
                    if (dff->FanoutGates()[k]->FanoutGates()[0]->GetGType() != GType::G_NOT) {
                        dff->FanoutGates()[k]->if_ppi_ = true;
                        single_ppi.emplace_back(dff->FanoutGates()[0]);
                    } else {
                        dff->FanoutGates()[k]->if_ppi_ = true;
                    }
                }
            }
        }
        pnlist_->GetPPIGates().emplace_back(single_ppi);
        assert(single_ppi.size() == pnlist_->ScanChains()[i]->ScanChainDFFs().size());
    }

}

void DRC::IdentifyPPOs() {
    for (int i = 0; i < pnlist_->ScanChains().size(); ++i) {
        std::vector<Gate*> single_ppo;
        for (int j = 0; j < pnlist_->ScanChains()[i]->ScanChainDFFs().size(); ++j) {
            Gate* dff = pnlist_->ScanChains()[i]->ScanChainDFFs()[j];
            dff->FaninGates()[3]->if_ppo_ = true;
            single_ppo.emplace_back(dff->FaninGates()[3]);
        }
        pnlist_->GetPPOGates().emplace_back(single_ppo);
        assert(single_ppo.size() == pnlist_->ScanChains()[i]->ScanChainDFFs().size());
    }
}

// TODO(wsl): default value?
void DRC::IdentifyCTs() {
  auto& ct_values = stil_parser_->CT2Value();
  for (auto& ct : ct_values) {
    auto& ct_name = ct.first;
    auto& ct_val = ct.second;
    auto& pipo_gates = pnlist_->GetNetName2PIPO();
    Gate* ct_gate = pipo_gates[ct_name];
    LOG_ASSERT(ct_gate, "ERROR: ct gate is nullptr");
    pnlist_->GetCTGates().emplace_back(ct_gate);
    auto gid = ct_gate->GetGId();
    if (ct_val == "CT0") {
      pnlist_->GetCT2Values()[gid] = ConstraintValue::CT_0;
    } else if (ct_val == "CT1") {
      pnlist_->GetCT2Values()[gid] = ConstraintValue::CT_1;
    } else if (ct_val == "CTX") {
      pnlist_->GetCT2Values()[gid] = ConstraintValue::CT_X;
    } else if (ct_val == "CTN0") {
      pnlist_->GetCT2Values()[gid] = ConstraintValue::CT_N0;
    } else if (ct_val == "CTN1") {
      pnlist_->GetCT2Values()[gid] = ConstraintValue::CT_N1;
    } else if (ct_val == "CTFREE") {
      pnlist_->GetCT2Values()[gid] = ConstraintValue::CT_FREE;
    } else {
      LOG_ASSERT(false, "ERROR: unsupported .ictest.spf constraint value");
    }
  }
}

void DRC::IdentifyTM() {
  auto& tm_values = stil_parser_->TM2Value();
  for (auto& tm : tm_values) {
    auto& tm_name = tm.first;
    auto& tm_val = tm.second;
    auto& pipo_gates = pnlist_->GetNetName2PIPO();
    Gate* tm_gate = pipo_gates[tm_name];
    LOG_ASSERT(tm_gate, "ERROR: tm gate is nullptr");
    pnlist_->GetTMGates().emplace_back(tm_gate);
    auto gid = tm_gate->GetGId();
    if (tm_val == "1") {
      pnlist_->GetTM2Values()[gid] = LogicValue::LOGIC_1;
    } else if (tm_val == "0") {
      pnlist_->GetTM2Values()[gid] = LogicValue::LOGIC_0;
    } else {
      LOG_ASSERT(false, "ERROR: unsupported .ictest.spf test mode value");
    }
  }
}

void DRC::IdentifySE() {
  auto& se_values = stil_parser_->SE2Value();
  for (auto& se : se_values) {
    auto& se_name = se.first;
    auto& se_val = se.second;
    auto& pipo_gates = pnlist_->GetNetName2PIPO();
    Gate* se_gate = pipo_gates[se_name];
    LOG_ASSERT(se_gate, "ERROR: se gate is nullptr");
    pnlist_->GetSEGates().emplace_back(se_gate);
    auto gid = se_gate->GetGId();
    if (se_val == "1") {
      pnlist_->GetSE2Values()[gid] = LogicValue::LOGIC_1;
    } else if (se_val == "0") {
      pnlist_->GetSE2Values()[gid] = LogicValue::LOGIC_0;
    } else {
      LOG_ASSERT(false, "ERROR: unsupported .ictest.spf scan enable value");
    }
  }
}

void DRC::IdentifyCLKs() {
  auto& clk2offstate = stil_parser_->CLK2Offstate();
  for (auto& clk : clk2offstate) {
    auto& clk_name = clk.first;
    auto& offstate = clk.second;
    auto& pipo_gates = pnlist_->GetNetName2PIPO();
    Gate* clk_gate = pipo_gates[clk_name];
    LOG_ASSERT(clk_gate, "ERROR: clk gate is nullptr");
    pnlist_->GetClkGates().emplace_back(clk_gate);
    auto gid = clk_gate->GetGId();
    if (offstate == "1") {
      pnlist_->GetClkOffState()[gid] = LogicValue::LOGIC_1;
    } else if (offstate == "0") {
      pnlist_->GetClkOffState()[gid] = LogicValue::LOGIC_0;
    } else {
      LOG_ASSERT(false, "ERROR: unsupported .ictest.spf clock value");
    }
  }
}

void DRC::IdentifyScanChains() {
  std::vector<uint8_t> inv_dff_flag(pnlist_->GetPrimNetlist().size(),0);
  auto& scan_chain_defs = stil_parser_->ScanChainDefs();
  auto& pipo_gates = pnlist_->GetNetName2PIPO();
  auto& scan_chains = pnlist_->ScanChains();
  int inv_nums = 0;
  bool is_mux_a = false;
  for (const auto& scan_chain_def : scan_chain_defs) {
    Gate* si_gate = pipo_gates[scan_chain_def.si_name_];
    //    Gate* so_gate = pipo_gates[scan_chain_def.so_name_];
    Gate* so_gate = nullptr;
    for (auto gptr : pnlist_->GetPOGates()) {
      if (gptr->GetInstName() == scan_chain_def.so_name_) {
        so_gate = gptr;
        break;
      }
    }
    ScanChain* scan_chain =
        new ScanChain(scan_chain_def.chain_name_, si_gate, so_gate);

    std::queue<Gate*> gate_queue;
    gate_queue.push(so_gate);

    std::unordered_set<Gate*> visit;
    while (!gate_queue.empty()) {
      auto gate = gate_queue.front();
      gate_queue.pop();

      Gate* fanin;
      LOG_ASSERT(gate != nullptr,
                 "ERROR: trace scan chain, gate pointer cannot be nullptr");

//      if (gate->GetGType() == GType::G_DFF) {
//        if(inv_nums & 1){
//          LOG_ERROR("ERROR: not support odd nums INV gate between 2 dff si path")
//        }
//        inv_nums = 0;
//      }
      if (gate->FaninSize() == 0) {
        if (gate->GetGType() == GType::G_PI &&
            gate->GetInstName() == scan_chain_def.si_name_) {
          scan_chain->ScanChainGates().emplace_back(gate);
        } else {
          LOG_ASSERT(false, "ERROR: trace scan chain, wrong si pin");
        }
        // TODO(wsl): test_so -> multi-test_si ?
        continue;
      } else if (gate->FaninSize() == 1) {
        fanin = gate->FaninGates()[0];
        scan_chain->ScanChainGates().emplace_back(gate);
        if (gate->GetGType() == GType::G_NOT) {
          inv_nums++;
        }
      } else if (gate->GetGType() == GType::G_DFF) {
        LOG_ASSERT(gate->FaninSize() == 4,
                   "ERROR: trace scan chain, DFF gate fanin size != 4");
        fanin = gate->FaninGates()[3];
        scan_chain->ScanChainDFFs().emplace_back(gate);
        scan_chain->ScanChainGates().emplace_back(gate);
      } else if (gate->GetGType() == GType::G_MUX) {
        LOG_ASSERT(gate->FaninSize() == 3,
                   "ERROR: trace scan chain, MUX gate fanin size != 3");
        uint8_t input_index = is_mux_a ? 1 : 2;
        fanin = gate->FaninGates()[input_index];
        scan_chain->ScanChainGates().emplace_back(gate);
      } else if (gate->GetGType() == GType::G_DLAT) {
        LOG_ASSERT(gate->FaninSize() == 4,
                   "ERROR: trace scan chain, DLAT gate fanin size != 4");
        fanin = gate->FaninGates()[3];
        scan_chain->ScanChainGates().emplace_back(gate);
      } else {
        LOG_ASSERT(
            false,
            "ERROR: trace scan chain, gate has multiple input but not MUX/DFF");
      }

      if (visit.find(fanin) == visit.end()) {
        gate_queue.push(fanin);
        visit.insert(fanin);
      }
    }
    inv_nums=0;
    uint32_t sc_dff_nums=scan_chain->ScanChainDFFs().size();
    scan_chain->GetScanChainsInvFlag().resize(sc_dff_nums);

    //invter count
    for(uint32_t i=0;i<scan_chain->ScanChainGates().size();++i){
      uint32_t from_si_idx=scan_chain->ScanChainGates().size()-1-i;
      Gate* gptr=scan_chain->ScanChainGates()[from_si_idx];
      if(gptr->GetGType()==GType::G_NOT){
        inv_nums++;
      }else if(gptr->GetGType()==GType::G_DFF){
        LOG_ASSERT(inv_dff_flag[gptr->GetGId()]==0, "ERROR: gate pointer aready inv in other scan chain");
        if(inv_nums & 1){
          inv_dff_flag[gptr->GetGId()]=1;
        }
      }else if(gptr->GetGType()==GType::G_PO){
        LOG_ASSERT(from_si_idx==0,"SO is not last gate in scan chain");
        if(inv_nums & 1){
          inv_dff_flag[gptr->GetGId()]=1;
        }
      }
    }
    pnlist_->GetDffInvFlag()=inv_dff_flag;
    for(uint32_t i=0;i<sc_dff_nums;++i){
      Gate* gptr=scan_chain->ScanChainDFFs()[i];
      if(inv_dff_flag[gptr->GetGId()]){
        scan_chain->GetScanChainsInvFlag()[i]=true;
      }else{
        scan_chain->GetScanChainsInvFlag()[i]=false;
      }
    }
    scan_chains.emplace_back(scan_chain);
  }
}

void DRC::IdentifyNonScanDFFs() {
  auto& dff_gates = pnlist_->GetDFFGates();
  auto& non_dff_gates = pnlist_->GetNonScanDFFGates();
  auto& scan_chains = pnlist_->ScanChains();

  std::unordered_set<Gate*> visited;
  for (const auto scan_chain : scan_chains) {
    for (const auto dff : scan_chain->ScanChainDFFs()) {
      LOG_ASSERT(dff != nullptr, "ERROR: gate pointer cannot be nullptr");
      visited.insert(dff);
    }
  }

  for (const auto& dff_gate : dff_gates) {
    LOG_ASSERT(dff_gate != nullptr, "ERROR: gate pointer cannot be nullptr");
    if (visited.find(dff_gate) == visited.end()) {
      non_dff_gates.emplace_back(dff_gate);
    }
  }
}

void DRC::DumpScanChain(const std::string& scan_chain_file) {
  std::ofstream ofs(scan_chain_file);
  if (!ofs.is_open()) {
    LOG_ERROR("ERROR: failed open file scan chain dump file: " +
              scan_chain_file);
    return;
  }
  for (auto& scan_chain : pnlist_->ScanChains()) {
    ofs << scan_chain->GetChainName() << "\n";
    auto& chain_dff = scan_chain->ScanChainDFFs();
    auto iter = chain_dff.rbegin();
    for (; iter != chain_dff.rend(); iter++) {
      auto dff_name = (*iter)->GetInstName();
      ofs << dff_name << "\n";
    }
  }
  ofs.close();
}

void DRC::Markchaingates() {
    for (int i = 0; i < pnlist_->ScanChains().size(); ++i) {
        for (int j = 0; j < pnlist_->ScanChains()[i]->ScanChainGates().size(); ++j) {
            if (!pnlist_->ScanChains()[i]->ScanChainGates()[j]->if_ppo_ &&
            !pnlist_->ScanChains()[i]->ScanChainGates()[j]->if_ppi_) {
                pnlist_->ScanChains()[i]->ScanChainGates()[j]->if_chain_ = true;
            }
        }
    }
    for (int i = 0; i < pnlist_->GetSEGates().size(); ++i) {
        pnlist_->GetSEGates()[i]->if_chain_ = true;
    }
    for (int i = 0; i < pnlist_->GetDLATGates().size(); ++i) {
        pnlist_->GetDLATGates()[i]->if_chain_ = true;
    }
    for (int i = 0; i < pnlist_->GetTMGates().size(); ++i) {
        pnlist_->GetTMGates()[i]->if_chain_ = true;
    }
    for (int i = 0; i < pnlist_->GetClkGates().size(); ++i) {
        pnlist_->GetClkGates()[i]->if_chain_ = true;
    }
}

}  // namespace ictest