// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_SAF_H
#define ICTEST_SAF_H

#include <cstdint>
#include <string>

#include "common/Define.h"
#include "netlist/Gate.h"

namespace ictest {

class SAF {
 public:
  SAF(Gate* fault_gate, SAFType fault_type, const std::string& fault_name)
      : saf_id_(-1),
        saf_type_(fault_type),
        saf_gate_(fault_gate),
        saf_status_(SAFStatus::UC),
        saf_name_(fault_name),
        saf_eqv_("UC") {}

  inline void SetSAFId(int32_t fid) { saf_id_ = fid; }
  inline int32_t GetSAFId() const { return saf_id_; }

  inline void SetSAFGate(Gate* fault_gate) { saf_gate_ = fault_gate; }
  inline Gate* GetSAFGate() const { return saf_gate_; }

  inline void SetSAFStatus(SAFStatus fault_status) {
    saf_status_ = fault_status;
  }
  inline const SAFStatus& GetSAFStatus() const { return saf_status_; }

  inline void SetSAFName(const std::string& fault_name) {
    saf_name_ = fault_name;
  }
  inline std::string GetSAFName() const { return saf_name_; }

  inline void SetSAFEqv(const std::string& fault_eqv) { saf_eqv_ = fault_eqv; }
  inline std::string GetSAFEqv() const { return saf_eqv_; }

  inline void SetSAFEqvGroupId(int id) { saf_eqv_group_idx_ = id; }
  inline int GetSAFEqvGroupId() const { return saf_eqv_group_idx_; }

  inline void SetSAFType(SAFType fault_type) { saf_type_ = fault_type; }
  inline SAFType GetSAFType() const { return saf_type_; }

 private:
  int32_t saf_id_{-1};
  int32_t saf_eqv_group_idx_{-1};
  Gate* saf_gate_{nullptr};
  SAFType saf_type_{SAFType::UNKNOWN_TYPE};
  SAFStatus saf_status_{SAFStatus::UC};
  std::string saf_name_{""};   // fault name at cell pin
  std::string saf_eqv_{"UC"};  // TMax is "--"/ TESSENT is "EQ"
};
}  // namespace ictest
#endif  // ICTEST_SAF_H
