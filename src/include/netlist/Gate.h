// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_GATE_H
#define ICTEST_GATE_H
#include <cassert>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_set>

#include "common/Define.h"

namespace ictest {

class Gate {
 public:
  Gate() {}
  Gate(GType gtype, const std::string& inst_name)
      : g_id_(-1),
        g_dpi_(-1),
        g_dpo_(-1),
        g_type_(gtype),
        g_instance_name_(inst_name),
        g_is_cell_border_(false),
        g_is_fb_loop_(false),
        g_dff_master_id_(-1),
        if_ppi_(false),
        if_ppo_(false),
        if_chain_(false){}
  ~Gate() {}

  inline void SetGId(const int32_t& gid) { g_id_ = gid; }
  inline int32_t GetGId() const { return g_id_; }

  inline void SetGType(enum GType gtype) { g_type_ = gtype; }
  inline GType GetGType() const { return g_type_; }

  inline void SetDPI(const int32_t& dpi) { g_dpi_ = dpi; }
  inline int32_t GetDPI() const { return g_dpi_; }

  inline void SetDPO(const int32_t& dpo) { g_dpo_ = dpo; }
  inline int32_t GetDPO() const { return g_dpo_; }

  inline void SetInstName(const std::string& inst_name) {
    g_instance_name_ = inst_name;
  }
  inline std::string GetInstName() const { return g_instance_name_; }

  inline void SetCellBorder(bool is_cell_border) {
    g_is_cell_border_ = is_cell_border;
  }
  inline bool IsCellBorder() const { return g_is_cell_border_; }

  inline void SetFBLoop(bool is_fb_loop) { g_is_fb_loop_ = is_fb_loop; }
  inline bool IsFBLoop() const { return g_is_fb_loop_; }

  inline void SetDffMasterId(int32_t dff_master_id) {
    g_dff_master_id_ = dff_master_id;
  }
  inline int32_t GetDffMasterId() const { return g_dff_master_id_; }

  // set/get pred/succ pin
  std::vector<int32_t>& PredPinId() { return g_pred_pin_id_; }
  std::vector<int32_t>& SuccPinId() { return g_succ_pin_id_; }
  std::vector<int32_t>& PredNetId() { return g_pred_net_id_; }
  void AddToPredNetId(int32_t idx) { g_pred_net_id_.push_back(idx); }

  // fanins gate set/get
  std::vector<Gate*>& FaninGates() { return g_fanin_gates_; }
  void AddToFaninGates(Gate* gate) { g_fanin_gates_.push_back(gate); }

  // fanouts gate set/get
  std::vector<Gate*>& FanoutGates() { return g_fanout_gates_; }
  void AddToFanoutGates(Gate* gate) { g_fanout_gates_.push_back(gate); }

  inline int32_t FaninSize() const { return g_fanin_gates_.size(); }
  inline int32_t FanoutSize() const { return g_fanout_gates_.size(); }

  // gate ctr/non_ctr
  LogicValue GetCtrInVal() const;
  LogicValue GetCtrOutVal() const;
  LogicValue GetNonCtrInVal() const;
  LogicValue GetNonCtrOutVal() const;

 public:
    bool if_ppi_;
    bool if_ppo_;
    bool if_chain_;
 private:
  int32_t g_id_{-1};  // -1 invalid
  int32_t g_dpi_{-1};
  int32_t g_dpo_{-1};
  GType g_type_{GType::G_UNKNOWN};
  std::string g_instance_name_;

  std::vector<Gate*> g_fanin_gates_;
  std::vector<Gate*> g_fanout_gates_;

  // for connection
  std::vector<int32_t> g_pred_pin_id_;
  std::vector<int32_t> g_succ_pin_id_;
  std::vector<int32_t> g_pred_net_id_;

  //
  bool g_is_cell_border_{false};
  bool g_is_fb_loop_{false};
  int32_t g_dff_master_id_{-1};
};
}  // namespace ictest
#endif  // ICTEST_GATE_H
