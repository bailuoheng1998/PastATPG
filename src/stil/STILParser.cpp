// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "stil/STILParser.h"

#include <fstream>
#include <sstream>
#include <string>

#include "util/Log.h"
#include "util/StringUtil.h"

namespace ictest {
bool STILParser::ReadICTestSPF(const std::string& spf_name) {
  std::ifstream file(spf_name);
  if (!file.is_open()) {
    LOG_ERROR("ERROR:failed to open .ictest.spf file : " + spf_name);
    return false;
  }

  std::string line;
  while (std::getline(file, line)) {
    StringUtil::TrimWinCRLine(line);
    if (line.empty()) {
      continue;
    }

    std::istringstream iss(line);
    std::string cmd;
    if (!(iss >> cmd)) {
      break;  // file format error
    }
    if(cmd[0]=='#'){
      continue;
    }
    if (cmd == "add_pi") {
      ReadInputNames(iss, pi_names_);
    } else if (cmd == "add_po") {
      ReadInputNames(iss, po_names_);
    } else if (cmd == "add_si") {
      ReadInputNames(iss, si_names_);
    } else if (cmd == "add_so") {
      ReadInputNames(iss, so_names_);
    } else if (cmd == "add_se") {
      ReadInputValue(iss, se2value_);
    } else if (cmd == "add_tm") {
      ReadInputValue(iss, tm2value_);
    } else if (cmd == "add_clk") {
      ReadInputValue(iss, clk2offstate_);
    } else if (cmd == "add_input_constraint") {
      ReadInputValue(iss, ct2value_);
    } else if (cmd == "add_scan_chain") {
      ReadScanChain(iss, scan_chain_defs_);
    } else {
      LOG_ERROR("ERROR:unsupported .ictest.spf command");
    }
  }

  return true;
}

void STILParser::ReadInputNames(std::istringstream& iss,
                                std::vector<std::string>& names) {
  std::string name;
  while (iss >> name) {
    names.emplace_back(name);
  }
}

void STILParser::ReadInputValue(
    std::istringstream& iss,
    std::unordered_map<std::string, std::string>& names) {
  std::string name, value;
  iss >> name >> value;
  names[name] = value;
}

void STILParser::ReadScanChain(std::istringstream& iss,
                               std::vector<ScanChainDef>& scan_chains) {
  std::string chain_name, input, output;
  iss >> chain_name >> input >> output;
  scan_chains.emplace_back(chain_name, input, output);
}

void STILParser::DumpICTestSPF(const std::string& file_path) {
  std::ofstream file(file_path);

  // write pi
  WriteNames(file, "add_pi", pi_names_);
  // write po
  WriteNames(file, "add_po", po_names_);
  // write si
  WriteNames(file, "add_si", si_names_);
  // write so
  WriteNames(file, "add_so", so_names_);
  // write se
  WriteValue(file, "add_se", se2value_);
  // write tm
  WriteValue(file, "add_tm", tm2value_);
  // write clk
  WriteValue(file, "add_clk", clk2offstate_);
  // write input constraint
  WriteValue(file, "add_input_constraint", ct2value_);
  // write scan chain
  WriteScanChain(file, "add_scan_chain", scan_chain_defs_);
}

void STILParser::WriteNames(std::ofstream& file, std::string&& cmd,
                            std::vector<std::string>& names) {
  if (names.empty()) return;

  file << cmd;
  for (const auto& name : names) {
    file << " " << name;
  }
  file << "\n";
}

void STILParser::WriteValue(
    std::ofstream& file, std::string&& cmd,
    std::unordered_map<std::string, std::string>& names) {
  if (names.empty()) return;

  file << cmd << " ";
  for (const auto& [name, value] : names) {
    file << name << " " << value << "\n";
  }
}

void STILParser::WriteScanChain(std::ofstream& file, std::string&& cmd,
                                std::vector<ScanChainDef>& scan_chains) {
  if (scan_chains.empty()) return;

  file << cmd << " ";
  for (const auto& scan_chain : scan_chains) {
    file << scan_chain.chain_name_ << " " << scan_chain.si_name_ << " "
         << scan_chain.so_name_ << "\n";
  }
  //  file << "\n";
}
}  // namespace ictest