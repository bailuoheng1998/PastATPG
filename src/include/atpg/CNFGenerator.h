//
// Created by tianpengyu on 23-4-4.
//

#ifndef ICTEST_CNFGENERATOR_H
#define ICTEST_CNFGENERATOR_H

#include "atpg/ATPGDefine.h"
#include "fault/SAF.h"
#include "netlist/PrimNetlist.h"
#include "pattern/Pattern.h"
#include "simulation/SimEngine.h"
#include "util/TimeUtil.h"
#include <fstream>
#include <errno.h>
#include <zlib.h>

#include "zhangsat/System.h"
#include "zhangsat/ParseUtils.h"
#include "zhangsat/Options.h"
#include "zhangsat/Dimacs.h"
#include "zhangsat/SATSolver.h"

class ATPGEngine;

using namespace Minisat;

namespace ictest {

class Solver {
  public:
    Solver () {
        cnf_var_intv_.clear();
        var_vals_.clear();
        res.clear();
    }
    std::vector<int> cnf_var_intv_;
    std::vector<int> var_vals_;
    std::set<int> x_vars_;
    std::vector<int> gid2dchain_vars_;
    std::vector<int> gid_po2pi_;
    std::vector<int> gid_fsite2po_;
    std::vector<Gate*> gid_po2pi_gptr_;
    Gate* fault_site_;
    std::vector<std::string> res;
    int max_var_size_ = 0;
    int max_dpi_ = 0;
    double cnf_trans_time_ = 0;
    double solve_time_ = 0;

    void add(int var) {
       cnf_var_intv_.push_back(var);
    }
    int val(int var_id) {
       return var_vals_[var_id];
    }
    int solve();
    int solve_help_checkx();
    int solve(ATPGEngine* d_test_engine, SAF* fault);
    void RecordActFS2NEAR(vector<int> &var_old2new, vector<int> &var_new2old, Minisat::SATSolver& S, vector<double> &initial_activity);
    bool CheckXbits(vector<vector<int>> &cnf_new, Minisat::SATSolver &S);
};

class CNFGenerator {
 public:
  CNFGenerator(PrimNetlist* prim) {
    prim_ = prim;
    model_.assign(2 * prim->GetPrimNetlist().size(), LOGIC_x);
    //log_info_ = false;
      log_info_ = true;

    for (int i = 0; i < prim_->ScanChains().size(); ++i) {
      for (int j = 0; j < prim_->ScanChains()[i]->ScanChainDFFs().size(); ++j) {
        is_scan_cell.insert(prim_->ScanChains()[i]->ScanChainDFFs()[j]);
      }
    }
  }
  ~CNFGenerator() {}

  void LogOn() { log_info_ = true; }
  void LogOff() { log_info_ = false; }
  void SetSolver(Solver* solver) { solver_ = solver; }
  void ResetSolver() { solver_ = nullptr; }
  void WriteCNF(vector<int>& var_old2new, vector<int>& var_new2old);
  int ReadSATResultMini(vector<bool>& old_vars_signs, vector<int>& var_new2old);
  int ReadSATResultKis(vector<bool>& old_vars_signs, vector<int>& var_new2old);
  int ReadSATResultCal(vector<bool>& old_vars_signs, vector<int>& var_new2old);
  int ReadSATResultPotopara(vector<bool>& old_vars_signs, vector<int>& var_new2old);

  void AddClauseToSolver(std::vector<int> literals) {
    for (auto lit : literals) {
      solver_->add(lit);
    }
    solver_->add(0);
  }

  bool CollectDChainCircuitClauseSAFL2(SAF* saf_fault,
                                       std::vector<int>& fanin_mark);
  bool CollectDChainCircuitClauseSAFL2TPI(SAF* saf_fault,
                                          std::vector<int>& fanin_mark);


  void CollectGateL2(Gate* target, int out, std::vector<int> fanins);
  void CollectL2AND(Gate* target, int out, std::vector<int> fanins);
  void CollectL2OR(Gate* target, int out, std::vector<int> fanins);
  void CollectL2NAND(Gate* target, int out, std::vector<int> fanins);
  void CollectL2NOR(Gate* target, int out, std::vector<int> fanins);
  void CollectL2MUX2(Gate* target, int out, int sel, int in0, int in1);
  void CollectL2XOR(Gate* target, int out, std::vector<int> fanins);
  void CollectL2XNOR(Gate* target, int out, std::vector<int> fanins);
  void CollectL2XOR2(Gate* target, int out, int in0, int in1);
  void CollectL2XNOR2(Gate* target, int out, int in0, int in1);
  void CollectL2BUF(Gate* target, int out, int in);
  void CollectL2INV(Gate* target, int out, int in);
  int GetVarIdOfGoodGate(Gate* gate);
  int GetVarIdOfBadGate(Gate* gate);

  // L4 logic system with Z.
  void CollectGateL4(Gate* target, int out, std::vector<int> fanins,
                     std::vector<int>& x_cone);
  void CollectL4AND(Gate* target, int out, std::vector<int> fanins,
                    std::vector<int>& x_cone);
  void CollectL4OR(Gate* target, int out, std::vector<int> fanins,
                   std::vector<int>& x_cone);
  void CollectL4NAND(Gate* target, int out, std::vector<int> fanins,
                     std::vector<int>& x_cone);
  void CollectL4NAND2(Gate* target, int out, std::vector<int> fanins,
                      std::vector<int>& x_cone);
  void CollectL4NOR(Gate* target, int out, std::vector<int> fanins,
                    std::vector<int>& x_cone);
  void CollectL4MUX2(Gate* target, int out, int sel, int in0, int in1,
                     std::vector<int>& x_cone);
  void CollectL4XOR(Gate* target, int out, std::vector<int> fanins,
                    std::vector<int>& x_cone);
  void CollectL4XNOR(Gate* target, int out, std::vector<int> fanins,
                     std::vector<int>& x_cone);
  void CollectL4BUF(Gate* target, int out, int in, std::vector<int>& x_cone);
  void CollectL4INV(Gate* target, int out, int in, std::vector<int>& x_cone);

  // L4 logic system Z as X.
  void CollectGateL3(Gate* target, int out, std::vector<int> fanins,
                     std::vector<int>& x_cone);
  void CollectL3AND(Gate* target, int out, std::vector<int> fanins,
                    std::vector<int>& x_cone);
  void CollectL3OR(Gate* target, int out, std::vector<int> fanins,
                   std::vector<int>& x_cone);
  void CollectL3NAND(Gate* target, int out, std::vector<int> fanins,
                     std::vector<int>& x_cone);
  void CollectL3NOR(Gate* target, int out, std::vector<int> fanins,
                    std::vector<int>& x_cone);
  void CollectL3MUX2(Gate* target, int out, int sel, int in0, int in1,
                     std::vector<int>& x_cone);
  void CollectL3XOR(Gate* target, int out, std::vector<int> fanins,
                    std::vector<int>& x_cone);
  void CollectL3XNOR(Gate* target, int out, std::vector<int> fanins,
                     std::vector<int>& x_cone);
  void CollectL3BUF(Gate* target, int out, int in, std::vector<int>& x_cone);
  void CollectL3INV(Gate* target, int out, int in, std::vector<int>& x_cone);

  int GetVarIdOfL4SecondBit(int var_id);

  void PassValToModel(Gate* gate, int good_val, int bad_val) {
    model_[gate->GetGId()] = good_val;
    model_[gate->GetGId() + prim_->GetPrimNetlist().size()] = bad_val;
  }
  void DebugModel(SAF* target_fault, const std::vector<int>& fanin_mark) {
    // debug out circuit model.
    std::cout << "Good model" << std::endl;
    for (int l = 0; l < prim_->GetLevelPrimNetlist().size(); l++) {
      for (int n = 0; n < prim_->GetLevelPrimNetlist()[l].size(); n++) {
        auto gate = prim_->GetLevelPrimNetlist()[l][n];
        if (gate->GetDPI() < 0 || gate->GetDPO() < 0) {
          continue;
        }
        if (fanin_mark[gate->GetGId()]) {
          std::cout << gate->GetInstName() << ":(" << model_[gate->GetGId()]
                    << ") ";
        }
      }
      std::cout << std::endl;
    }

    std::cout << "Bad model" << std::endl;
    for (int l = 0; l < prim_->GetLevelPrimNetlist().size(); l++) {
      for (int n = 0; n < prim_->GetLevelPrimNetlist()[l].size(); n++) {
        auto gate = prim_->GetLevelPrimNetlist()[l][n];
        if (gate->GetDPI() < 0 || gate->GetDPO() < 0) {
          continue;
        }
        if (fanin_mark[gate->GetGId()]) {
          std::cout << gate->GetInstName() << ":("
                    << model_[prim_->GetPrimNetlist().size() + gate->GetGId()]
                    << ") ";
        }
      }
      std::cout << std::endl;
    }

    for (int i = 0; i < prim_->GetPrimNetlist().size(); i++) {
      auto gate = prim_->GetPrimNetlist()[i];
      if (gate->GetDPI() < 0 || gate->GetDPO() < 0) {
        continue;
      }
      if (fanin_mark[gate->GetGId()] &&
          model_[gate->GetGId()] != EvalGate(gate)) {
        std::cout << "Wrong Model at: " << gate->GetInstName() << std::endl;
        std::cout << "Gate type: " << static_cast<int>(gate->GetGType())
                  << std::endl;
        std::cout << "Fanin vals: " << std::endl;
        for (auto fanin : gate->FaninGates()) {
          std::cout << fanin->GetInstName() << "(" << model_[fanin->GetGId()]
                    << ") ";
        }
        std::cout << std::endl;

        std::cout << "Out val: " << model_[gate->GetGId()] << std::endl;
      }
    }
  }
  int EvalGate(Gate* gate);

 private:
  PrimNetlist* prim_;
  Solver* solver_;
  bool log_info_;

  int max_extra_var_id_;
  std::vector<int> model_;
  std::unordered_set<Gate*> is_scan_cell;
};



}  // namespace ictest

#endif  // ICTEST_CNFGENERATOR_H
