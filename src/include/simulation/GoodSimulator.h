// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_PPGOODSIMULATOR_H
#define ICTEST_PPGOODSIMULATOR_H

#include <util/TimeUtil.h>

#include "Evaluation.h"
#include "EventsQueue.h"
#include "SimData.h"
#include "SimVal.h"
#include "common/Option.h"
#include "netlist/PrimNetlist.h"
#include "pattern/Pattern.h"
#include <fstream>
#include <iomanip>
#include <queue>
#include <sstream>

namespace ictest {
class GoodSimulator {
public:
  GoodSimulator(const GoodSimulator &) = delete;

  GoodSimulator operator=(const GoodSimulator &) = delete;

  explicit GoodSimulator(PrimNetlist *prim_netlist)
      : prim_netlist_(prim_netlist),
        event_queue_(prim_netlist->MaxLevel() + 1, prim_netlist->NumGates()) {
    po_ids_.resize(prim_netlist->NumPOs());
    for (int i = 0; i < prim_netlist->NumPOs(); ++i) {
      po_ids_[i] = prim_netlist->GetPOGates()[i]->GetGId();
    }
    constexpr int max_fanin_nums = 8;

    evaluate_.AllocFaninMem(max_fanin_nums);

    AllocAllMem();
  }

  ~GoodSimulator() = default;

  ///////////////////////////////////////
  // main func
  void GoodSimulation(Pattern &sp);

  void ResetGoodMachine();

  // scan chain
  void ScanChainAnalyse();

  uint64_t GetEventsCount() const { return events_count_; }

  ////////////////////////////
  // option
  void SetSimData(SimData *simData) { simdata_ = simData; }

protected:
  void AllocAllMem();

  void ResetDLatch(const val64_t &init_val, int frame_id);

  void ResetNonScanDFF(const val64_t &init_val, int frame_id);

  void ResetAllDFF(const val64_t &init_val, int frame_id);

  //////////////////////////////////////////////////////////
  // load scan DFFs
  void LoadScanDFFs(PatternCycle &load_pc, int frame_id);

  void Load(std::vector<std::vector<val64_t>> &load_vals, int shift_nums);

  void ParallelAssignDffs(const std::vector<std::vector<val64_t>> &load_vals,
                          int shift_nums);

  //////////////////////////////////////////////////////////
  // capture
  void Capture(PatternCycle &pc);

  void ForcePIs(PatternCycle &pc, int frame_id);

  void ForcePIs(const std::vector<val64_t> &pi_vals, int frame_id);

  void InitTieVal(int useful_bits, int frame_id);

  void PulseCaptureClock(PatternCycle &pc, int frame_id);

  void EventDriveSim(int frame_id);

  void AddGateToEventQueue(Gate &gptr);

  void SetGateInputVal(Gate &gptr, int frame_id);

  /**
   * @brief 计算之后更新good machine 的逻辑值
   * @param gid 当前待更新的逻辑门索引
   * @param frame_id 更新到值数组的哪个frame
   */
  bool UpdateGoodVal(const val64_t &new_value, gid_t gid, int frame_id);

  bool EvalAndUpdate(Gate &gptr, int frame_id);

  /////////////////////////////////////////////////////////
  // scan chain check
  void ForcePIs(std::vector<std::pair<int, val64_t>> &pi2vals, int frame_id);

  void MarkActiveNSC();

  void StatisticChains();

  void
  TraceAndMarkMasterClkPath(const std::vector<std::pair<int, val64_t>> &sigs);

  void TraceAndMarkHoldSigPath(const std::vector<std::pair<int, val64_t>> &sigs,
                               GateClassify mark_val);

  void TraceScanChainsSo2Si();

  void EventDrivenPathMark(GateClassify marks_val);

private:
  EventsQueue event_queue_;
  Evaluation<val64_t> evaluate_;
  PrimNetlist *prim_netlist_{nullptr};

  std::vector<int> po_ids_;
  std::vector<int> slave2master_;
  std::vector<std::vector<val64_t>> good_machine_vals_;
  std::vector<std::vector<val64_t>> good_master_vals_;

  uint64_t events_count_ = 0;
  SimData *simdata_{nullptr};

  friend class SAFSimulator;
};
} // namespace ictest
#endif // ICTEST_PPGOODSIMULATOR_H
