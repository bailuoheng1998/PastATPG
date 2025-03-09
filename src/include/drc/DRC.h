// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_DRC_H
#define ICTEST_DRC_H

#include <fstream>

#include "netlist/PrimNetlist.h"
#include "netlist/ScanChain.h"
#include "stil/STILParser.h"

namespace ictest {

class DRC {
 public:
  DRC() {}
  DRC(PrimNetlist* pnlist, STILParser* stil_parser)
      : pnlist_(pnlist), stil_parser_{stil_parser} {}

  // todo:
  void RunDRC();
  void IdentifyPIs();         // sort pi
  void IdentifyPOs();         // sort po
  void IdentifyPPIs();         // SET UP Ppi
  void IdentifyPPOs();         // SET UP Ppo
  void IdentifyCLKs();        // clk and offstate
  void IdentifySE();          // scan enable
  void IdentifyTM();          // test mode
  void IdentifyCTs();         // constraint
  void IdentifyScanChains();  // scan chain trace
  void IdentifyNonScanDFFs();
  void Markchaingates();

  void DumpScanChain(const std::string& scan_chain_file);

 private:
  void SortByNames(std::vector<std::string>& names, std::vector<Gate*>& gates);

  // Design Rule Check
  // void CheckClkRule_C0();
  // void CheckScanChainRule_S0();

 private:
  PrimNetlist* pnlist_{nullptr};
  STILParser* stil_parser_{nullptr};
};
}  // namespace ictest
#endif  // ICTEST_DRC_H
