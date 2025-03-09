// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3
#include "netlist/Cell.h"

#include "util/Log.h"

namespace ictest {

void Cell::AddSAFsAtInputPin(SAF* fptr, int32_t pin_id) {
  LOG_ASSERT(pin_id < GetInputPins().size(),
             "ERROR: pin_id < GetInputPins().size()");
  if (c_fanin_safs_.empty()) {
    c_fanin_safs_.resize(2 * GetInputPins().size());
  }

  int32_t sa = 0;
  if (fptr->GetSAFType() == SAFType::SA1) {
    sa = 1;
  }
  c_fanin_safs_[pin_id * 2 + sa] = fptr;
}

void Cell::AddSAFsAtOutputPin(SAF* fptr, int32_t pin_id) {
  LOG_ASSERT(pin_id < GetOutputPins().size(),
             "ERROR: pin_id < GetOutputPins().size()");
  if (c_fanout_safs_.empty()) {
    c_fanout_safs_.resize(2 * GetOutputPins().size());
  }

  int32_t sa = 0;
  if (fptr->GetSAFType() == SAFType::SA1) {
    sa = 1;
  }
  c_fanout_safs_[pin_id * 2 + sa] = fptr;
}

}  // namespace ictest