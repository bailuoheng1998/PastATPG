//
// Created by tianpengyu on 23-4-4.
//

#ifndef ICTEST_SATENGINE_H
#define ICTEST_SATENGINE_H

#include "atpg/CNFGenerator.h"
#include "fault/SAF.h"
#include "fault/SAFList.h"
#include "netlist/PrimNetlist.h"
#include "pattern/Pattern.h"
#include "simulation/SimEngine.h"
#include "util/TimeUtil.h"

class CNFGenerator;
class Solver;

namespace ictest {
class SATEngine {
 public:
  SATEngine(PrimNetlist* prim, SAFList* saf_list, SimEngine* sim,
            std::string dump_path) {
    prim_ = prim;
    saf_list_ = saf_list;
    sim_ = sim;
    dump_path_ = dump_path;
  }

  ~SATEngine() {}

  // L2 Logic System.
  void GenerateCompactTestSetForSAFListHard2DTfaults(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompaction(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompactionDynamic(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompactionDynamic_2(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompareMini(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompareKis(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompareKis_dynamic(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompareCal(int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompaction_czt(int pattern_per_sim = 64);
  int SATGeneratePat4Pf(SAF* primary_fault, std::vector<int> &new_pattern, double &cnf_gen_time, double &cnf_trans_time, double &solver_time);
  int SATGeneratePat4Sf(SAF* secondary_fault, std::vector<int> &new_pattern, double &cnf_gen_time, double &cnf_trans_time, double &solver_time);
  void GenerateCompactTestSetForSAFListHard2DTfaultsComparePotopara(int thread_num, int pattern_per_sim = 64);
  void GenerateCompactTestSetForSAFListHard2DTfaultsCompareCal_dynamic(int pattern_per_sim = 64);


  void CollectPoByDistance(Gate* gate, std::vector<Gate*>& ordered_pos,
                           std::vector<int>& fault_cone);

  void DumpSAFList();

 private:
  PrimNetlist* prim_;
  SAFList* saf_list_;

  SimEngine* sim_;
  std::string dump_path_;

};
}  // namespace ictest

#endif  // ICTEST_SATENGINE_H
