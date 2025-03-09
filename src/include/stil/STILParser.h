// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_STILPARSER_H
#define ICTEST_STILPARSER_H

#include <string>
#include <unordered_map>
#include <vector>

#include "common/Define.h"

namespace ictest {

// scan chain define
struct ScanChainDef {
  std::string chain_name_;
  std::string si_name_;
  std::string so_name_;

  ScanChainDef(std::string chain_name, std::string si_name, std::string so_name)
      : chain_name_(std::move(chain_name)),
        si_name_(std::move(si_name)),
        so_name_(std::move(so_name)) {}
};

class STILParser {
 public:
  STILParser() = default;
  ~STILParser() = default;

  bool ReadICTestSPF(const std::string& spf_name);
  void DumpICTestSPF(const std::string& file_path);

  inline std::vector<std::string>& PINames() { return pi_names_; }
  inline std::vector<std::string>& PONames() { return po_names_; }
  inline std::vector<std::string>& SINames() { return si_names_; }
  inline std::vector<std::string>& SONames() { return so_names_; }
  inline std::vector<ScanChainDef>& ScanChainDefs() { return scan_chain_defs_; }

  inline std::unordered_map<std::string, std::string>& CLK2Offstate() {
    return clk2offstate_;
  }
  inline std::unordered_map<std::string, std::string>& SE2Value() {
    return se2value_;
  }
  inline std::unordered_map<std::string, std::string>& TM2Value() {
    return tm2value_;
  }
  inline std::unordered_map<std::string, std::string>& CT2Value() {
    return ct2value_;
  }

 private:
  void ReadInputNames(std::istringstream& iss, std::vector<std::string>& names);
  void ReadInputValue(std::istringstream& iss,
                      std::unordered_map<std::string, std::string>& names);
  void ReadScanChain(std::istringstream& iss,
                     std::vector<ScanChainDef>& scan_chains);

  void WriteNames(std::ofstream& file, std::string&& cmd,
                  std::vector<std::string>& names);
  void WriteValue(std::ofstream& file, std::string&& cmd,
                  std::unordered_map<std::string, std::string>& names);
  void WriteScanChain(std::ofstream& file, std::string&& cmd,
                      std::vector<ScanChainDef>& scan_chains);

 private:
  std::vector<std::string> pi_names_;
  std::vector<std::string> po_names_;
  std::vector<std::string> si_names_;
  std::vector<std::string> so_names_;
  std::unordered_map<std::string, std::string> clk2offstate_;
  std::unordered_map<std::string, std::string> se2value_;
  std::unordered_map<std::string, std::string> tm2value_;
  std::unordered_map<std::string, std::string> ct2value_;  // input constraint
  std::vector<ScanChainDef> scan_chain_defs_;
};

}  // namespace ictest
#endif  // ICTEST_STILPARSER_H
