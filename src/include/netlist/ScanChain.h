// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_SCANCHAIN_H
#define ICTEST_SCANCHAIN_H
#include <cstdint>
#include <list>
#include <string>
#include <unordered_map>
#include <vector>

#include "Gate.h"

namespace ictest {
class ScanChain {
 public:
  ScanChain(const std::string& chain_name, Gate* si, Gate* so)
      : chain_name_(chain_name), scan_in_(si), scan_out_(so) {}

  // todo:
  void TraceChain() {}

  inline size_t GetChainLength() const { return chain_dffs_.size(); }

  inline void SetChainName(const std::string& name) {
    chain_name_ = chain_name_;
  }
  inline std::string GetChainName() const { return chain_name_; }

  inline int32_t GetScanInversion() const { return scan_inversion_; }

  inline Gate* GetScanInGate() const { return scan_in_; }

  inline Gate* GetScanOutGate() const { return scan_out_; }

  inline Gate* GetScanEnableGate() const { return scan_enable_; }

  inline Gate* GetScanClkGate() const { return scan_clk_; }

  std::vector<Gate*>& ScanChainDFFs() { return chain_dffs_; }

  std::vector<Gate*>& ScanChainGates() { return chain_gates_; }

  std::vector<bool>& GetScanChainsInvFlag() { return inv_flag_; }
 private:
  std::string chain_name_;
  int32_t scan_inversion_{0};
  Gate* scan_in_{nullptr};
  Gate* scan_out_{nullptr};
  Gate* scan_enable_{nullptr};
  Gate* scan_clk_{nullptr};
  std::vector<Gate*> chain_dffs_;  //
  std::vector<bool> inv_flag_;
  std::vector<Gate*> chain_gates_;
};
}  // namespace ictest
#endif  // ICTEST_SCANCHAIN_H
