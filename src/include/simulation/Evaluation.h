// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_EVALUATION_H
#define ICTEST_EVALUATION_H
#include "SimVal.h"
#include "common/Define.h"

namespace ictest {
template<class T_VAL_4>
class Evaluation {
  using base_bit_type = typename T_VAL_4::base_bit_type;

  // only for Evaluation
  enum PIN_DEFINE {
    // for latch and dff
    SET_PIN_IDX = 0,
    RESET_PIN_IDX = 1,
    CLK_PIN_IDX = 2,
    DATA_PIN_IDX = 3,
    // only for latch
    LATCH_HOLD_VAL_IDX = 4,
    // only for dff
    DFF_MASTER_VAL_IDX = 4,
    DFF_SLAVE_VAL_IDX = 5,
    DFF_MASTER_NEW_VAL_IDX = 6,
    // for mux
    MUX_SEL_IDX = 0,
    MUX_IN0_IDX = 1,
    MUX_IN1_IDX = 2,
  };

 public:
  Evaluation() {
    constexpr int min_vals_size = 8;
    AllocFaninMem(min_vals_size);
  };

  explicit Evaluation(int max_fanin_size) { AllocFaninMem(max_fanin_size); }

  void AllocFaninMem(int max_fanin_size) {
    constexpr int dff_input_nums = 8;
    int size = std::max(dff_input_nums, max_fanin_size);
    pin_vals_.resize(size);
  }

  ~Evaluation() noexcept {
//    delete[] pin_vals_;
  }

  Evaluation(const Evaluation& rhs) = delete;

  Evaluation(Evaluation&& rhs) = delete;

  Evaluation& operator=(const Evaluation& rhs) = delete;

  Evaluation& operator=(Evaluation&& rhs) = delete;

  // init input vals
  void SetFaninSize(size_t fanin_size) { fanin_size_ = fanin_size; }

  void SetPinValue(int pin_id, const T_VAL_4& val) { pin_vals_[pin_id] = val; }

  // for dlatch
  void SetLatchHoldVal(const T_VAL_4& hold_val) {
    pin_vals_[LATCH_HOLD_VAL_IDX] = hold_val;
  }

  // for dff
  void SetDffMasterVal(const T_VAL_4& master_val) {
    pin_vals_[DFF_MASTER_VAL_IDX] = master_val;
  }

  void SetDffSlaveVal(const T_VAL_4& slave_val) {
    pin_vals_[DFF_SLAVE_VAL_IDX] = slave_val;
  }

  // only for dff master
  T_VAL_4 GetDffMasterVal() const { return pin_vals_[DFF_MASTER_NEW_VAL_IDX]; }

  // interface func
  T_VAL_4 Evaluate(GType logic_type) {
    size_t fanin_size = fanin_size_;

    LOG_ASSERT(fanin_size > 0, "Fanin size wrong");
    // compute logic value according to gate's logic type
    switch (logic_type) {
      case GType::G_PO:
      case GType::G_BUF:
      case GType::G_ABUF:
      case GType::G_BRH: {
        return Eval_BUF();
      }
      case GType::G_NOT: {
        return Eval_NOT();
      }
      case GType::G_AND: {
        return Eval_AND();
      }
      case GType::G_NAND: {
        return Eval_NAND();
      }
      case GType::G_OR: {
        return Eval_OR();
      }
      case GType::G_NOR: {
        return Eval_NOR();
      }
      case GType::G_XOR: {
        return Eval_XOR();
      }
      case GType::G_XNOR: {
        return Eval_XNOR();
      }
      case GType::G_MUX: {
        return Eval_MUX();
      }
      case GType::G_DFF: {
        return Eval_DFF();
      }
      case GType::G_DLAT: {
        return Eval_DLATCH();
      }
      case GType::G_CMUX: {
        return Eval_CMUX();
      }
      default: {
        LOG_ASSERT(false, "ERROR:not support simulation gate type: " +
                              std::to_string(static_cast<int>(logic_type)))
      }
    }
  }

 private:
  // help func
  T_VAL_4 Select(int sel_idx, int in0_idx, int in1_idx) const {
    const T_VAL_4& sel = pin_vals_[sel_idx];
    const T_VAL_4& in0 = pin_vals_[in0_idx];
    const T_VAL_4& in1 = pin_vals_[in1_idx];

    T_VAL_4 same_input_bits = (~(in0 ^ in1)).MaskTo0((in0 ^ in1).GetBitXPos());
    T_VAL_4 new_value = ((sel & in1) | (~sel & in0));
    new_value = (same_input_bits & in0) | (~same_input_bits & new_value);
    return new_value;
  }
  static T_VAL_4 Select(const T_VAL_4& sel, const T_VAL_4& in0, const T_VAL_4& in1){

    T_VAL_4 same_input_bits = (~(in0 ^ in1)).MaskTo0((in0 ^ in1).GetBitXPos());
    T_VAL_4 new_value = ((sel & in1) | (~sel & in0));
    new_value = (same_input_bits & in0) | (~same_input_bits & new_value);
    return new_value;
  }
  // if triggle level==true high level triggle else low level triggle
  T_VAL_4 Latch(int data_pin_idx, int hold_val_idx, bool triggle_level) const {
    const T_VAL_4& set = pin_vals_[SET_PIN_IDX];
    // reset signal
    const T_VAL_4& reset = pin_vals_[RESET_PIN_IDX];


    if (!triggle_level) {
      std::swap(hold_val_idx, data_pin_idx);
    }
    T_VAL_4 new_value = Select(CLK_PIN_IDX, hold_val_idx, data_pin_idx);
    new_value= Select(set,new_value,T_VAL_4(~0,0));
    new_value= Select(reset,new_value,T_VAL_4(0,~0));
    // set and reset all 1,q=x
    base_bit_type conflict_bits=set.GetBit1Pos() & reset.GetBit1Pos();
    new_value = new_value.MaskToX(conflict_bits);
    return new_value;
  }

  // logic func
  T_VAL_4 Eval_BUF() const { return pin_vals_[0]; }

  T_VAL_4 Eval_NOT() const { return ~pin_vals_[0]; }

  T_VAL_4 Eval_OR() const {
    T_VAL_4 new_val = pin_vals_[0];
    for (int fi = 1; fi < fanin_size_; fi++) {
      new_val = new_val | pin_vals_[fi];
    }
    return new_val;
  }

  T_VAL_4 Eval_NOR() const { return ~Eval_OR(); }

  T_VAL_4 Eval_AND() const {
    T_VAL_4 new_val = pin_vals_[0];
    for (int fi = 1; fi < fanin_size_; fi++) {
      new_val = new_val & pin_vals_[fi];
    }
    return new_val;
  }

  T_VAL_4 Eval_NAND() const { return ~Eval_AND(); }

  T_VAL_4 Eval_XOR() const {
    T_VAL_4 new_val = pin_vals_[0];
    for (int fi = 1; fi < fanin_size_; fi++) {
      new_val = new_val ^ pin_vals_[fi];
    }
    return new_val;
  }

  T_VAL_4 Eval_XNOR() const { return ~Eval_XOR(); }

  T_VAL_4 Eval_MUX() const {
    return Select(MUX_SEL_IDX, MUX_IN0_IDX, MUX_IN1_IDX);
  }

  T_VAL_4 Eval_CMUX() const {
    // only support 3-input mux
    const T_VAL_4& sel = pin_vals_[0];
    const T_VAL_4& in0 = pin_vals_[1];
    const T_VAL_4& in1 = pin_vals_[2];

    T_VAL_4 new_value = (sel & in1) | (~sel & in0);
    new_value = new_value.MaskToX(sel.GetBitXPos());

    return new_value;
  }

  T_VAL_4 Eval_DLATCH() const {
    return Latch(DATA_PIN_IDX, LATCH_HOLD_VAL_IDX, true);
  }

  T_VAL_4 Eval_DFF() {
    T_VAL_4 new_slave_value =
        Latch(DFF_MASTER_VAL_IDX, DFF_SLAVE_VAL_IDX, true);

    pin_vals_[DFF_MASTER_NEW_VAL_IDX] =
        Latch(DATA_PIN_IDX, DFF_MASTER_VAL_IDX, false);
    return new_slave_value;
  }

 private:
  //  std::vector<func> func_array;
  size_t fanin_size_ = 0;
//  T_VAL_4* pin_vals_{nullptr};
  std::vector<T_VAL_4> pin_vals_;
};
}  // namespace ictest

#endif  // ICTEST_EVALUATION_H
