// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_PRIMNETLIST_H
#define ICTEST_PRIMNETLIST_H
#include <string>
#include <unordered_map>


#include "Gate.h"
#include "ScanChain.h"
#include "common/Define.h"
#include "fault/SAF.h"
namespace ictest {

// learning analysis
struct GateStatus {
  int32_t gid_{-1};
  LogicValue logic_value_{LogicValue::LOGIC_X};
  TieValue tie_value_{TieValue::TIE_FREE};
  ConstraintValue ct_value_{ConstraintValue::CT_FREE};
  bool is_blocked_{false};
  bool is_scheduled_{false};
  bool is_blocked_by_constraint_{false};
};

class PrimNetlist {
 public:
  PrimNetlist()
      : top_module_name_(""),
        num_gates_(0),
        num_pis_(0),
        num_pos_(0),
        max_level_(0),
        max_width_(0) {}

  ~PrimNetlist() { ClearPrimNetlist(); }
  void ClearPrimNetlist();
  void InitPrimNetlist();
  void ResizePrimNetlist();
  bool ReadPrimVY(const std::string& vy_file);
  void ConnectPrimNetlist();

  void AddPOGates();
  // stem's branch
  void AddBranchGates();
  void AddPin2BranchGates();

  // levelize call dpi and dpoInitPrimNetlist
  void Levelize();
  void SetDPI();
  void SetDPO();
  void SetLevelPrimNetlist();

  // check validation
  bool CheckPrimNetlistValidation();


  // set/get module name
  void SetTopModName(std::string name) { top_module_name_ = name; }
  inline std::string GetTopModuleName() const { return top_module_name_; }

  // set/get netlist type
  void AnalyseNetlistType();
  inline NetlistType GetNetlistType() const { return netlist_type_; }
  inline void SetNetlistType(NetlistType ntype) { netlist_type_ = ntype; }
  std::string NetlistType2Str();

  inline int32_t NumGates() const { return num_gates_; }
  inline int32_t NumPIs() const { return num_pis_; }
  inline int32_t NumPOs() const { return num_pos_; }
  inline int32_t NumDFFs() const { return dff_gates_.size(); }
  inline int32_t NumDLATs() const { return dlat_gates_.size(); }
  inline int32_t NumTies() const { return tie_gates_.size(); }
  inline int32_t MaxLevel() const { return max_level_; }
  inline int32_t MaxWidth() const { return max_width_; }

  std::vector<Gate*>& GetPrimNetlist() { return prim_netlist_; }

  void SetGate2Status() { gate_status_.resize(num_gates_); }
  std::vector<GateStatus>& GetGateStatus() { return gate_status_; }

  std::vector<Gate*>& GetPIGates() { return pi_gates_; }
  std::vector<Gate*>& GetPOGates() { return po_gates_; }
  std::vector<std::vector<Gate*>>& GetPPIGates() { return ppi_gates_; }
  std::vector<std::vector<Gate*>>& GetPPOGates() { return ppo_gates_; }
  std::vector<Gate*>& GetTieGates() { return tie_gates_; }
  std::vector<Gate*>& GetStemGates() { return stem_gates_; }
  std::vector<Gate*>& GetDLATGates() { return dlat_gates_; }
  std::vector<Gate*>& GetDFFGates() { return dff_gates_; }
  std::vector<Gate*>& GetNonScanDFFGates() { return non_scan_dffs_; }
  std::vector<Gate*>& GetClkGates() { return clk_gates_; }
  std::vector<Gate*>& GetSEGates() { return scan_enable_gates_; }
  std::vector<Gate*>& GetTMGates() { return test_mode_gates_; }
  std::vector<Gate*>& GetCTGates() { return constraint_gates_; }
  std::vector<std::vector<Gate*>>& GetLevelPrimNetlist() {
    return level_prim_netlist_;
  }

  std::vector<ScanChain*>& ScanChains() { return scan_chains_; }



  std::unordered_map<int32_t, std::string>& GetNet2Names() {
    return net2names_;
  }
  std::unordered_map<int32_t, std::string>& GetPin2Names() {
    return pin2names_;
  }
  std::vector<Gate*>& GetPin2BranchGates() { return pin2branch_; }
  std::unordered_map<std::string, Gate*>& GetNetName2PIPO() {
    return netname2pipogptr_;
  }

  std::unordered_map<int32_t, LogicValue>& GetClkOffState() {
    return clk_off_state_;
  }
  std::unordered_map<int32_t, LogicValue>& GetSE2Values() { return se2values_; }
  std::unordered_map<int32_t, LogicValue>& GetTM2Values() { return tm2values_; }
  std::unordered_map<int32_t, ConstraintValue>& GetCT2Values() {
    return ct2values_;
  }

  void GetMaxWidthNums(int32_t& max_width);
  void GetMaxFaninNums(int32_t& max_fanin_nums);
  void GetBufferNums(int32_t& total_buffer, int32_t& buffer_cas);
  void GetInverterNums(int32_t& total_inv, int32_t& inv_cas);
  void GetStemNums(int32_t& stem_nums);
  auto& GetDffInvFlag(){
    return dff_parallel_inv_flag_;
  }
  // report指令，方便shell调用
  void ReportPI();
  void ReportPO();
  void ReportCLK();
  void ReportScanChain();
  void ReportNetlist();
  // dump
  void DumpPrimNetlist(const std::string& file_path);
  void DumpVerilog(const std::string& file_path);

  // todo ParseNetlist()

 protected:
 private:
  std::string top_module_name_;

  NetlistType netlist_type_{NetlistType::UNKNOWN_CIRCUIT};
  int32_t num_gates_{0};
  int32_t num_pis_{0};
  int32_t num_pos_{0};
  int32_t max_level_{0};
  int32_t max_width_{0};

  std::vector<Gate*> prim_netlist_;
  std::vector<GateStatus> gate_status_;
  std::vector<Gate*> pi_gates_;
  std::vector<Gate*> po_gates_;
  std::vector<std::vector<Gate*>> ppi_gates_;
  std::vector<std::vector<Gate*>> ppo_gates_;
  std::vector<Gate*> tie_gates_;
  std::vector<Gate*> stem_gates_;
  std::vector<Gate*> dlat_gates_;
  /// include scan dffs and non-scan dffs
  std::vector<Gate*> dff_gates_;
  std::vector<Gate*> non_scan_dffs_;
  std::vector<Gate*> clk_gates_;
  std::vector<Gate*> scan_enable_gates_;
  std::vector<Gate*> test_mode_gates_;
  std::vector<Gate*> constraint_gates_;
  std::vector<std::vector<Gate*>> level_prim_netlist_;

  std::vector<ScanChain*> scan_chains_;

  std::vector<uint8_t> dff_parallel_inv_flag_;
  std::unordered_map<int32_t, std::string> net2names_;
  std::unordered_map<int32_t, std::string> pin2names_;
  std::vector<Gate*> pin2branch_;
  std::unordered_map<std::string, Gate*> netname2pipogptr_;
  // clk gate id, clk off state
  std::unordered_map<int32_t, LogicValue> clk_off_state_;
  // scan enable gate, enable signal value
  std::unordered_map<int32_t, LogicValue> se2values_;
  // test mode gate, signal value
  std::unordered_map<int32_t, LogicValue> tm2values_;
  // input constraint gate, constraint signal value
  std::unordered_map<int32_t, ConstraintValue> ct2values_;

  //
  bool has_rising_edge_dffs_{false};
  bool has_falling_edge_dffs_{false};
};
}  // namespace ictest
#endif  // ICTEST_PRIMNETLIST_H
