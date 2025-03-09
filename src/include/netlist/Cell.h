// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_CELL_H
#define ICTEST_CELL_H
#include <cassert>
#include <string>
#include <vector>

#include "common/Define.h"
#include "fault/SAF.h"

namespace ictest {
class Gate;
class SAF;

class Cell {
 public:
  Cell() {}
  Cell(int32_t id, int32_t logic_type, const std::string& inst_name)
      : c_id_(id), c_type_(logic_type), c_instance_name_(inst_name) {}
  ~Cell() {}

  inline void SetCId(int32_t cid) { c_id_ = cid; }
  inline int32_t GetCId() const { return c_id_; }

  inline void SetCType(int32_t ctype) { c_type_ = ctype; }
  inline int32_t GetCType() const { return c_type_; }

  inline void SetInstanceName(const std::string& instance_name) {
    c_instance_name_ = instance_name;
  }
  inline std::string GetInstanceName() const { return c_instance_name_; }

  inline int32_t FaninSize() const { return c_fanin_cells_.size(); }
  inline int32_t FanoutSize() const { return c_fanout_cells_.size(); }

  // set/get
  std::vector<int32_t>& GetInputNets() { return c_input_nets_; }
  std::vector<int32_t>& GetInputPins() { return c_input_pins_; }
  std::vector<int32_t>& GetInoutNets() { return c_inout_nets_; }
  std::vector<int32_t>& GetInoutPins() { return c_inout_pins_; }
  std::vector<int32_t>& GetOutputNets() { return c_output_nets_; }
  std::vector<int32_t>& GetOutputPins() { return c_output_pins_; }
  std::vector<Cell*>& GetFaninCells() { return c_fanin_cells_; }
  std::vector<Cell*>& GetFanoutCells() { return c_fanout_cells_; }
  std::vector<Gate*>& GetFaninGates() { return c_fanin_gates_; }
  std::vector<Gate*>& GetFanoutGates() { return c_fanout_gates_; }
  std::vector<SAF*>& GetFaninSAFs() { return c_fanin_safs_; }
  std::vector<SAF*>& GetFanoutSAFs() { return c_fanout_safs_; }


  // fault

  void AddSAFsAtInputPin(SAF* fptr, int32_t pin_id);
  void AddSAFsAtOutputPin(SAF* fptr, int32_t pin_id);


 private:
  int32_t c_id_{-1};
  int32_t c_type_{-1};
  //  std::string c_module_name_;
  std::string c_instance_name_;

  std::vector<int32_t> c_input_nets_;
  std::vector<int32_t> c_input_pins_;

  std::vector<int32_t> c_inout_nets_;
  std::vector<int32_t> c_inout_pins_;

  std::vector<int32_t> c_output_nets_;
  std::vector<int32_t> c_output_pins_;

  std::vector<Cell*> c_fanin_cells_;  // connect cell
  std::vector<Cell*> c_fanout_cells_;

  std::vector<Gate*> c_fanin_gates_;
  std::vector<Gate*> c_fanout_gates_;

  std::vector<SAF*> c_fanin_safs_;
  std::vector<SAF*> c_fanout_safs_;

};
}  // namespace ictest
#endif  // ICTEST_CELL_H
