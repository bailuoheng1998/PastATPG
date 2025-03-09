// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_PPFAULTSIMULATOR_H
#define ICTEST_PPFAULTSIMULATOR_H

#include "GoodSimulator.h"
#include "fault/SAFList.h"
namespace ictest {

class SAFSimulator {
 public:
  SAFSimulator() = delete;

  SAFSimulator(PrimNetlist* pnlist, Option& option, GoodSimulator& gsim)
      : prim_netlist_(pnlist),
        option_(&option),
        fevent_queue_(pnlist->MaxLevel() + 1, pnlist->NumGates()),
        gsim_(gsim) {
    constexpr int max_fanin_nums = 8;
    evaluate_.AllocFaninMem(max_fanin_nums);
    AllocAllMem();
  }
  /**
 * @brief 组合电路故障仿真主函数 未使用优化策略
   */
  int BasicComPPFSim(Pattern& sp, std::vector<SAF*>& sa_flist);

  /**
   * @brief 全扫描电路和单个capture部分扫描仿SA故障仿真主函数
   */
  int BasicScan( Pattern &sp,std::vector<SAF*> &sa_flist);

  /**
   * @brief 单故障单并行向量下的故障仿真
   */

  SAFStatus PPSFDT(SAF& fptr, Pattern &sp) ;

  ///将上一个frame计算的故障值传播到下一个frame
  void UpdateMemEventsBasicScan(int frame_id);
  /**
   * @brief 向网标注入故障
   */
  bool InitFaultyEvent(Gate& gptr, const val64_t& faulty_val, int frame_id,
                       int ffr_head_id);
  /**
   * @brief 故障仿真的事件驱动主函数
   */
  void FaultEventDriveSim(int frame_id, int stop_flag);

  /**
   * @brief 故障仿真在po端口对比逻辑值 返回检测之后的故障分类状态
   */
  SAFStatus MeasurePOs(int frame_id);
  /**
   * @brief 故障仿真在ppo端口对比逻辑值 返回检测之后的故障分类状态
   */
  SAFStatus MeasureSOs(int frame_id);
  ///获取当前故障下 某个frame下逻辑门id为gid的逻辑门逻辑值 当前仅支持frame范围为0～2
  val64_t GetFaultVal(int frame_id, int gid) const;
  //获取当前故障下 某个frame下逻辑门id为dff_gid的逻辑门虚拟master 锁存器逻辑值 当前仅支持frame范围为0～2
  val64_t GetFaultMasterVal(int frame_id, int dff_gid) const;
  uint64_t GetEventsCount();

  void SetSAFList(SAFList& saf_list){   // gaojun temp be del
    saf_list_=&saf_list;
  }

  void SetSimData(SimData* simdata){
    simdata_=simdata;
  }
  void SetCnlist(CellNetlist* cnlist){
    cnlist_=cnlist;
  }
 private:
  void AllocAllMem();
  void ResetFaultMachineToGood();
  /**
   * @brief 计算之前向evaluate类赋值
   */
  void SetInputLineVal(Gate& gptr, int frame_id);
  /**
   * @brief 计算之后更新fault machine 的逻辑值
   */
  bool SetFaultVal(Gate& gptr, const val64_t& new_val, int frame_id);
  // update  by assign val
  bool UpdateFaultVal(const val64_t& new_value, int gate_id, int frame_id);
  /// 将当前逻辑门的扇出加入fault事件对列
  void AddActiveGateToEvents(Gate& gptr);
  ///根据故障时生成一个64位故障值，useful_bits表示从最低位开始填充了多少位，未填充故障值的位用X填充
  val64_t GenSAFaultValue(SAF* fptr, int useful_bits) const;

 private:
  SAFList* saf_list_{nullptr};
  //for edt compact record whitch dff has X value
  ictest::GoodSimulator& gsim_;
  PrimNetlist* prim_netlist_{nullptr};
  CellNetlist*cnlist_{nullptr};
  EventsQueue fevent_queue_;
  Evaluation<val64_t> evaluate_;

  const Option* option_{nullptr};
  const SimData *simdata_{nullptr};
  std::vector<std::vector<val64_t>> fault_vals_;
  std::vector<std::vector<val64_t>> fault_master_vals_;
  std::vector<std::vector<uint8_t>> fault_path_mark_;
  ///标记某个逻辑门是否有故障，有故障的情况下不计算，避免覆盖掉故障值。
  std::vector<uint8_t> has_fault_;
  std::vector<int> fault_po_ids_;
  std::vector<int> fault_gids_;
  std::vector<int> fault_dff_ids_;
  uint64_t fault_event_count_{0};
  bool fault_to_clk{false};
  bool sim_falling_edge{true};
  bool sim_rising_edge{true};
  ///目前暂时用于debug
  std::vector<std::pair<SAF*,uint64_t>> _inject_faults;
  std::vector<std::vector<std::pair<SAF*,uint64_t>>> _mask;
  std::vector<std::string> debug_pins_;
  std::vector<uint32_t> debug_gate_ids_;
  bool measure_ppo_=false;
  bool measure_po_=false;
};

}  // namespace ictest

#endif  // ICTEST_PPFAULTSIMULATOR_H
