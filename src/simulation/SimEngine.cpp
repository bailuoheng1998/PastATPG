// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "simulation/SimEngine.h"
namespace ictest {

SimEngine::SimEngine(PrimNetlist *pnlist, ictest::Option &sim_option)
    : gsim_(pnlist), fsim_(pnlist, sim_option, gsim_), sim_option_(&sim_option),
      pnlist_(pnlist) {
  BuildSimData();
  gsim_.SetSimData(simData_);
  fsim_.SetSimData(simData_);
  gsim_.ScanChainAnalyse();
  total_pat_num_ = 0;
}
SimEngine::~SimEngine() { delete simData_; }
/**
 * @brief 执行逻辑仿真
 * @param patterns 需要仿真的输入向量集合
 */
void SimEngine::RunGoodSimulation(std::vector<Pattern> &patterns) {
  sim_record_.clear();
  sim_record_.reserve(20);
  std::ofstream ofs;
  FUNC_TIME_REPORT_WITH_FILE_ICT
  int total_pat_nums = 0;
  auto start = std::chrono::steady_clock::now();
  gsim_.ResetGoodMachine();
  for (int i = 0; i < patterns.size(); ++i) {
    gsim_.GoodSimulation(patterns[i]);
    total_pat_nums += patterns[i].GetPatternCycle(0).LastParallelBitIdx();
  }
  auto end = std::chrono::steady_clock::now();
  auto cost_time = std::chrono::duration<double>(end - start);
  sim_record_.emplace_back("TopModuleName", pnlist_->GetTopModuleName());
  sim_record_.emplace_back("Circuit Type",
                           NetlistType2Name[pnlist_->GetNetlistType()]);
  sim_record_.emplace_back("PatternNums", std::to_string(total_pat_nums));
  auto gsim_key = "GSim time";
  sim_record_.emplace_back(gsim_key, std::to_string(cost_time.count()) + "s");
  sim_record_.emplace_back("GSim report at",
                           Timer::get_current_time_and_date());
}
/**
 * @brief 故障仿真函数，支持组合电路，全扫描电路和单capture的部分扫描电路 。
 * @param patterns 需要仿真的输入向量集合
 * @param saf_list 待仿真的故障列表
 */
void SimEngine::RunSAFSimulation(std::vector<Pattern> &patterns,
                                 std::vector<SAF *> &saf_list) {
  sim_record_.clear();
  sim_record_.reserve(20);
  int dt_nums = 0;

  for (auto &fptr : saf_list) {
    if (fptr->GetSAFStatus() == SAFStatus::DI ||
        fptr->GetSAFStatus() == SAFStatus::DS) {
      dt_nums++;
    }
  }
  int dt_count = 0;
  int total_pat_nums = 0;
  auto netlist_type = pnlist_->GetNetlistType();
  uint32_t epoch_flist_dt_count;
  LOG_INFO("{}...", "fault sim start")
  LOG_INFO("{:^15}{:^15}{:^15}{:^15}", "total_pat_nums", "detected", "coverage",
           "cost time")
  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < patterns.size(); ++i) {
    auto pattern_start = std::chrono::steady_clock::now();
    switch (netlist_type) {
    case NetlistType::COMB_CIRCUIT:
      epoch_flist_dt_count = fsim_.BasicComPPFSim(patterns[i], saf_list);
      break;
    case NetlistType::PART_SCAN_CIRCUIT:
    case NetlistType::FULL_SCAN_CIRCUIT:
      epoch_flist_dt_count = fsim_.BasicScan(patterns[i], saf_list);
      break;
    case NetlistType::SEQU_CIRCUIT:
      LOG_ASSERT(false, "NOT SUPPORT SEQU_CIRCUIT")
      break;
    case NetlistType::UNKNOWN_CIRCUIT:
      LOG_ASSERT(false, "NOT SUPPORT UNKNOWN TYPE CIRCUIT")
      break;
    }
    auto pattern_end = std::chrono::steady_clock::now();
    auto pattern_cost =
        std::chrono::duration<double>(pattern_end - pattern_start).count();
    dt_count += epoch_flist_dt_count;
    int curr_pattern_nums = patterns[i].GetPatternCycle(0).LastParallelBitIdx();
    total_pat_nums += (curr_pattern_nums == 0 ? 1 : curr_pattern_nums);
    auto coverage = (double)(dt_count + dt_nums) / saf_list.size();
    LOG_INFO("{0:^15}{1:^15}{2:^15}{3:^15}", total_pat_nums,
             epoch_flist_dt_count, fmt::format("{:.4f}%", coverage * 100),
             fmt::format("{:.4f}s", pattern_cost))
  }

  auto end = std::chrono::steady_clock::now();
  auto cost_time = std::chrono::duration<double>(end - start);

  sim_record_.emplace_back("TopModuleName", pnlist_->GetTopModuleName());
  sim_record_.emplace_back("SimType",
                           NetlistType2Name[pnlist_->GetNetlistType()]);
  sim_record_.emplace_back("PatternNums", std::to_string(total_pat_nums));
  sim_record_.emplace_back("FSim time (new netlist)",
                           std::to_string(cost_time.count()) + "s");
  StatisticSAFList(saf_list);
  sim_record_.emplace_back("Total faults", std::to_string(saf_list.size()));
  uint64_t event_count = fsim_.GetEventsCount();
  sim_record_.emplace_back("Total events", std::to_string(event_count));
  sim_record_.emplace_back("FSim report at",
                           Timer::get_current_time_and_date());


}

void SimEngine::RunSAFSimulation(std::vector<Pattern> &patterns,
                                 std::vector<SAF *> &saf_list, bool is_report) {
  sim_record_.clear();
  sim_record_.reserve(20);
  int dt_nums = 0;

  for (auto &fptr : saf_list) {
    if (fptr->GetSAFStatus() == SAFStatus::DI ||
        fptr->GetSAFStatus() == SAFStatus::DS) {
      dt_nums++;
    }
  }

  int dt_count = 0;
  int total_pat_nums = 0;
  auto netlist_type = pnlist_->GetNetlistType();
  uint32_t epoch_flist_dt_count;

  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < patterns.size(); ++i) {
    auto pattern_start = std::chrono::steady_clock::now();
    switch (netlist_type) {
    case NetlistType::COMB_CIRCUIT:
      epoch_flist_dt_count = fsim_.BasicComPPFSim(patterns[i], saf_list);
      break;
    case NetlistType::PART_SCAN_CIRCUIT:
    case NetlistType::FULL_SCAN_CIRCUIT:
      epoch_flist_dt_count = fsim_.BasicScan(patterns[i], saf_list);
      break;
    case NetlistType::SEQU_CIRCUIT:
      LOG_ASSERT(false, "NOT SUPPORT SEQU_CIRCUIT")
      break;
    case NetlistType::UNKNOWN_CIRCUIT:
      LOG_ASSERT(false, "NOT SUPPORT UNKNOWN TYPE CIRCUIT")
      break;
    }
    if (is_report) {
      auto pattern_end = std::chrono::steady_clock::now();
      auto pattern_cost =
          std::chrono::duration<double>(pattern_end - pattern_start).count();
      dt_count += epoch_flist_dt_count;
      int curr_pattern_nums =
          patterns[i].GetPatternCycle(0).LastParallelBitIdx();
      total_pat_nums += (curr_pattern_nums == 0 ? 1 : curr_pattern_nums);
      total_pat_num_ += (curr_pattern_nums == 0 ? 1 : curr_pattern_nums);
      auto coverage = (double)(dt_count + dt_nums) / saf_list.size();
      LOG_INFO("{0:^15}{1:^15}{2:^15}", total_pat_num_,
               epoch_flist_dt_count, fmt::format("{:.4f}%", coverage * 100));
    }
  }
  if (is_report) {
    auto end = std::chrono::steady_clock::now();
    auto cost_time = std::chrono::duration<double>(end - start);

    sim_record_.emplace_back("TopModuleName", pnlist_->GetTopModuleName());
    sim_record_.emplace_back("SimType",
                             NetlistType2Name[pnlist_->GetNetlistType()]);
    sim_record_.emplace_back("PatternNums", std::to_string(total_pat_nums));
    sim_record_.emplace_back("FSim cost time",
                             std::to_string(cost_time.count()) + "s");
    StatisticSAFList(saf_list);
    sim_record_.emplace_back("Total faults", std::to_string(saf_list.size()));
    uint64_t event_count = fsim_.GetEventsCount();
    sim_record_.emplace_back("Total events", std::to_string(event_count));
    sim_record_.emplace_back("FSim report at",
                             Timer::get_current_time_and_date());
  }
}


/**
 * @brief 初始化仿真所需要的数据
 */
void SimEngine::BuildSimData() {
  simData_ = new SimData;
  SimData &simdata = *simData_;
  val64_t logic_0 = val64_t(0, ~0);
  //  print(logic_0)
  auto &scan_chains = pnlist_->ScanChains();
  auto &netlist = pnlist_->GetPrimNetlist();
  simdata.chain_nums_ = pnlist_->ScanChains().size();
  simdata.si_so_port_.resize(simdata.chain_nums_);
  // init clk 2 pi idx
  for (int i = 0; i < pnlist_->GetPIGates().size(); ++i) {
    for (auto &clk : pnlist_->GetClkGates()) {
      if (clk->GetGId() == pnlist_->GetPIGates()[i]->GetGId()) {
        simdata.clk_to_pi_id_.emplace_back(clk->GetGId(), i);
      }
    }
  }
  for (int sc_id = 0; sc_id < pnlist_->ScanChains().size(); ++sc_id) {
    int lenth = scan_chains[sc_id]->ScanChainDFFs().size();
    simdata.max_chain_lenth_ = std::max(lenth, simdata.max_chain_lenth_);
  }

  for (int sc_id = 0; sc_id < pnlist_->ScanChains().size(); ++sc_id) {
    int si_id = scan_chains[sc_id]->GetScanInGate()->GetGId();
    int so_id = scan_chains[sc_id]->GetScanOutGate()->GetGId();
    simdata.SetChainSi(sc_id, si_id);
    simdata.SetChainSo(sc_id, so_id);
  }
  val64_t init_val;
  for (auto &tm_pair : pnlist_->GetTM2Values()) {
    Gate &gptr = *netlist[tm_pair.first];
    if (tm_pair.second == LogicValue::LOGIC_0) {
      init_val = logic_0;
    } else if (tm_pair.second == LogicValue::LOGIC_1) {
      init_val = ~logic_0;
    }
    simdata.AddTmSigs(tm_pair.first, init_val);
  }
  for (auto &se_pair : pnlist_->GetSE2Values()) {
    Gate &gptr = *netlist[se_pair.first];
    if (se_pair.second == LogicValue::LOGIC_0) {
      init_val = logic_0;
    } else if (se_pair.second == LogicValue::LOGIC_1) {
      init_val = ~logic_0;
    }
    simdata.AddSESigs(se_pair.first, init_val);
  }
  for (auto &clk_pair : pnlist_->GetClkOffState()) {
    Gate &gptr = *netlist[clk_pair.first];
    if (clk_pair.second == LogicValue::LOGIC_0) {
      init_val = logic_0;
    } else if (clk_pair.second == LogicValue::LOGIC_1) {
      init_val = ~logic_0;
    } else {
      LOG_ASSERT(false, "clk has X off state")
    }
    // todo specify master clock in drc
    // temp for find master clk
    if (pnlist_->GetClkOffState().size() == 1) {
      simdata.AddShiftPulseClks(clk_pair.first, init_val);
    } else {
      const auto &clk_name = netlist[clk_pair.first]->GetInstName();
      if (clk_name.find("clk") != std::string::npos) {
        simdata.AddShiftPulseClks(clk_pair.first, init_val);
      } else {
        simdata.set_reset_sigs_.emplace_back(clk_pair.first, init_val);
      }
    }
  }

  for (auto &tie : pnlist_->GetTieGates()) {
    Gate &gptr = *tie;
    if (gptr.GetGType() == GType::G_TIE0) {
      init_val = logic_0;
    } else if (gptr.GetGType() == GType::G_TIE1) {
      init_val = ~logic_0;
    } else if (gptr.GetGType() == GType::G_TIEX) {
      init_val = val64_t(0, 0);
    } else {
      LOG_ASSERT(false, "Not support Tie Gate Type")
    }
    simdata.tie_sigs_.emplace_back(gptr.GetGId(), init_val);
  }
}

} // namespace ictest