// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_SIMREPOTER_H
#define ICTEST_SIMREPOTER_H
#include "fault/SAF.h"
namespace ictest{

struct SimReporter {
  std::vector<std::pair<std::string, std::string>> sim_record_;  // key value
  double cpu_time_{0};
  uint32_t di_nums_{0};
  uint32_t ds_nums_{0};
  uint32_t undetectable_nums_{0};
  uint32_t totle_fault_nums_{0};
  uint32_t total_pat_num_{0};

  [[nodiscard]] double GetFaultCoverage() const {
    if (totle_fault_nums_ == 0) {
      return 1;
    }
    return static_cast<double>(ds_nums_ + di_nums_) / totle_fault_nums_;
  }

  [[nodiscard]] double GetTestCoverage() const {
    uint32_t totle_testable = totle_fault_nums_ - undetectable_nums_;
    if (totle_testable == 0) {
      return 1;
    }
    return static_cast<double>(ds_nums_ + di_nums_) / totle_testable;
  }

  void ClearReport() {
    di_nums_ = 0;
    ds_nums_ = 0;
    undetectable_nums_ = 0;
    totle_fault_nums_ = 0;
    total_pat_num_ = 0;
    cpu_time_=0;
    sim_record_.clear();
  }

  void InitReporter(const std::vector<SAF*>& sa_flist) {
    ClearReport();
    for (auto& fptr : sa_flist) {
      if (fptr == nullptr) {
        continue;
      }
      totle_fault_nums_++;
      auto fault_status = fptr->GetSAFStatus();
      if (fault_status == SAFStatus::UU || fault_status == SAFStatus::TI ||
          fault_status == SAFStatus::BL) {
        undetectable_nums_++;
      } else if (fault_status == SAFStatus::DI) {
        di_nums_++;
      } else if (fault_status == SAFStatus::DS) {
        ds_nums_++;
      }
    }
  }


  void StatisticSAFList(const std::vector<SAF*>& sa_flist) {
    std::vector<int> statistic(20, 0);
    for (auto& fptr : sa_flist) {
      auto status = fptr->GetSAFStatus();
      int status_idx = static_cast<int>(status);
      assert(status_idx < statistic.size());
      statistic[status_idx]++;
    }
    std::unordered_map<std::string, int> res_map;
    for (auto& kv : SAFStatus2Tmax) {
      int sidx = static_cast<int>(kv.first);
      if (res_map.find(kv.second) != res_map.end()) {
        res_map[kv.second] += statistic[sidx];
      } else {
        res_map[kv.second] = statistic[sidx];
      }
    }
    int dt = statistic[static_cast<int>(SAFStatus::DS)] +
             statistic[static_cast<int>(SAFStatus::DI)];
    double coverage = (double)dt / (double)sa_flist.size();
    for (auto& kv : res_map) {
      sim_record_.emplace_back(kv.first, std::to_string(kv.second));
    }
    sim_record_.emplace_back("totle fault", std::to_string(sa_flist.size()));
    sim_record_.emplace_back("fault coverage", std::to_string(100 * coverage) + "%");
  }

};
}  // namespace ictest
#endif  // ICTEST_SIMREPOTER_H
