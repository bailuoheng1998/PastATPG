// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_VALUE_H
#define ICTEST_VALUE_H

#include <bitset>
#include <cstdint>
#include <iostream>
#include <string>

#include "common/Define.h"
#include "util/Log.h"

namespace ictest {
template<class T>
struct SimVal {
 public:
  using base_bit_type = T;
  SimVal() : v0_(0), v1_(0) {}
  
  constexpr SimVal(base_bit_type v0, base_bit_type v1) : v0_(v0), v1_(v1) {}

  SimVal(const SimVal& rval) = default;

  SimVal& operator=(const SimVal& rval) = default;

//  SimVal& operator=(const SimVal& rval) {
////    assert(!rval.GetBitZPos());
////    *this = rval;
//    return *this;
//  }

  inline void Set0() {
    v0_ = 0;
    v1_ = -1;
  }
  inline void Set1() {
    v0_ = -1;
    v1_ = 0;
  }
  inline void SetX() {
    v0_ = 0;
    v1_ = 0;
  }
  inline void SetZ() {
    v0_ = -1;
    v1_ = -1;
  }

  // val1 == val2
  inline bool operator==(const SimVal& other) const {
    return (this->v0_ == other.v0_) && (this->v1_ == other.v1_);
  }
  static SimVal MergeByBitMask(const SimVal& merge_by_bit0, const SimVal& merge_by_bit1,base_bit_type bit_mask){
    base_bit_type v0=(~bit_mask & merge_by_bit0.v0_) | (bit_mask & merge_by_bit1.v0_);//bit 1 :merge_by_bit1
    base_bit_type v1=(~bit_mask & merge_by_bit0.v1_) | (bit_mask & merge_by_bit1.v1_);
    return SimVal(v0,v1);
  }
  static base_bit_type GetSameBits(const SimVal& first,const SimVal& other){
    return ~((first.v0_ ^ other.v0_) | (first.v1_ ^ other.v1_));
  }
  // val1 != val2
  bool operator!=(const SimVal& other) const { return !(*this == other); }

  // val1 & val2
  SimVal operator&(const SimVal& other) const {
    return SimVal(this->v0_ & other.v0_, this->v1_ | other.v1_);
  }
  // val1 | val2
  SimVal operator|(const SimVal& other) const {
    return SimVal(this->v0_ | other.v0_, this->v1_ & other.v1_);
  }
  // val1 ^ val2
  SimVal operator^(const SimVal& other) const {
    return SimVal((this->v1_ & other.v0_) | (this->v0_ & other.v1_),
                  (this->v1_ & other.v1_) | (this->v0_ & other.v0_));
  }
  // ~val
  SimVal operator~() const { return SimVal(this->v1_, this->v0_); }

  base_bit_type GetBit1Pos() const {  // only 1 is 1 in return uint_64
    base_bit_type mask = this->v0_ ^ this->v1_;
    return (this->v0_ & mask);
  }

  base_bit_type GetBit0Pos() const {  // only 1 is 1 in return uint_64
    base_bit_type mask = this->v0_ ^ this->v1_;
    return (this->v1_ & mask);
  }

  base_bit_type GetBitXPos() const {  // only x is 1 in return uint_64
    base_bit_type mask = ~(this->v0_ ^ this->v1_);
    return (~this->v1_) & mask;
  }

  base_bit_type GetBitZPos() const {  // only Z is 1 in return uint_64
    base_bit_type mask = ~(this->v0_ ^ this->v1_);
    return this->v0_ & mask;
  }

  //  base_bit_type GetBitXAndZPos() const {  // only x or z is 1 in return uint_64
  //    return ~(this->v0_ ^ this->v1_);
  //  }

  SimVal MaskToX(base_bit_type mask) const {
    return SimVal(this->v0_ & (~mask), this->v1_ & (~mask));
  }

  SimVal MaskToZ(base_bit_type mask) const {
    return SimVal(this->v0_ | mask, this->v1_ | mask);
  }

  SimVal MaskTo1(base_bit_type mask) const {
    return SimVal(this->v0_ | mask, this->v1_ & (~mask));
  }

  SimVal MaskTo0(base_bit_type mask) const {
    return SimVal(this->v0_ & (~mask), this->v1_ | mask);
  }

  SimVal MaskInv(base_bit_type mask) const { return *this ^ SimVal(mask, ~mask); }

  SimVal ExpandByBitIdx(int bit_idx) const {
    //*this :            X 1 0 1
    // bit idx            3 2 1 0
    // bit idx=3 return   x x x x
//    assert(bit_idx < 64);
    constexpr base_bit_type right1 = 1;
    base_bit_type bit0 = this->v0_ & (right1 << bit_idx);
    base_bit_type bit1 = this->v1_ & (right1 << bit_idx);
    if (bit0 > 0) {
      //            bit0 = 0xFFFFFFFFFFFFFFFF >> (64 - SIZE_OF_PACKET);
      bit0 = ~0;
    }
    if (bit1 > 0) {
      //            bit1 = 0xFFFFFFFFFFFFFFFF >> (64 - SIZE_OF_PACKET);
      bit1 = ~0;
    }
    return SimVal(bit0, bit1);
  }

  SimVal ExpandByBitMask(base_bit_type mask) const {
    base_bit_type bit0 = this->v0_ & (mask);
    base_bit_type bit1 = this->v1_ & (mask);
    if (bit0 > 0) {
      bit0 = ~0;
    }
    if (bit1 > 0) {
      bit1 = ~0;
    }
    return SimVal(bit0, bit1);
  }

  std::string GetSerialValue() const {
    std::string res;
    if (v0_ == V64_ONE && v1_ == V64_ZERO) {
      res = "1";
    } else if (v0_ == V64_ZERO && v1_ == V64_ONE) {
      res = "0";
    } else if (v0_ == V64_ZERO && v1_ == V64_ZERO) {
      res = "X";
    } else if (v0_ == V64_ONE && v1_ == V64_ONE) {
      res = "Z";
    } else {
      LOG_ASSERT(false, "ERROR: not supported logic value");
    }
    return res;
  }

  char ToChar(int bit_idx = 0) const {
    constexpr base_bit_type one = 1;
    base_bit_type v0 = ((one << bit_idx) & this->v0_);
    base_bit_type v1 = ((one << bit_idx) & this->v1_);
    if (v0 ^ v1) {///0 and 1
      if (v0) {
        return '1';
      } else {
        return '0';
      }
    } else if (v0 & v1) {
      return 'Z';
    } else {
      return 'X';
    }
  }

  std::string ToString(int lenth = sizeof (base_bit_type)*CHAR_BIT) const {
    constexpr int packet_size = sizeof (base_bit_type)*CHAR_BIT;

    std::string res;
    std::string bit0 = std::bitset<packet_size>(v0_).to_string();
    std::string bit1 = std::bitset<packet_size>(v1_).to_string();

    for (int i = packet_size - lenth; i < packet_size; ++i) {
      char b0 = bit0[i];
      char b1 = bit1[i];
      if (b0 == '0' && b1 == '1')
        res.push_back('0');
      else if (b0 == '1' && b1 == '0')
        res.push_back('1');
      else if (b0 == '0' && b1 == '0')
        res.push_back('X');
      else if (b0 == '1' && b1 == '1')
        res.push_back('Z');
    }
    return res;
  }

 public:
  base_bit_type v0_{0};
  base_bit_type v1_{0};
};

template<class T>
std::ostream & operator << (std::ostream & cout,const SimVal<T>& val){
  cout<<val.ToString();
  return cout;
}
using val64_t = SimVal<uint64_t>;
using val32_t = SimVal<uint32_t>;
using val16_t = SimVal<uint16_t>;
using val8_t = SimVal<uint8_t>;
}  // namespace ictest
#endif  // ICTEST_VALUE_H
