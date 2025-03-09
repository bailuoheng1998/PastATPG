// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_SIMENGINE_H
#define ICTEST_SIMENGINE_H
#include "SAFSimulator.h"
#include "SimData.h"
#include "common/Option.h"
#include "SimReporter.h"
namespace ictest {
class SimEngine : public SimReporter{
 public:

  SimEngine(PrimNetlist* pnlist, ictest::Option& sim_option);
  ~SimEngine();

 public:
  // good simulation
  void RunGoodSimulation(std::vector<Pattern>& patterns);

  // fault simulation
  void RunSAFSimulation(std::vector<Pattern>& patterns,
                        std::vector<SAF*>& saf_list);

  void RunSAFSimulation(std::vector<Pattern>& patterns,
                        std::vector<SAF*>& saf_list, bool is_report);

  //for build new netlist
  void BuildSimData();



 public:
  ictest::SimData* simData_{nullptr};

 private:
  ictest::GoodSimulator gsim_;
  ictest::SAFSimulator fsim_;

  ictest::Option* sim_option_;
  PrimNetlist* pnlist_;

};

}  // namespace ictest
#endif  // ICTEST_SIMENGINE_H
