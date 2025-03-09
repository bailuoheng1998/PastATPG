// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_DEFINE_H
#define ICTEST_DEFINE_H

#include <dirent.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>

//#define DEBUG_MODE
//#define DC_MODE

namespace ictest {

////////////////////////////////////////////////////
// gate logic type
enum class GType : std::uint8_t {
  G_UNKNOWN = 0,
  G_PI = 1,
  G_PO = 2,
//  G_PPI = 3,
//  G_PPO = 4,
  G_AND = 6,
  G_NAND = 7,
  G_OR = 8,
  G_NOR = 9,
  G_NOT = 10,
  G_BUF = 11,
  G_BRH = 12, // stem's branch   // gaojun ---- G_STEM
  G_XOR = 101,
  G_XNOR = 102,
  G_ABUF = 180, // buf assign
  G_TIE0 = 195, // gaojun ---- supply
  G_TIE1 = 196,
  G_TIEX = 181,
  G_TIEZ = 182,
  G_MUX = 197,
  G_DLAT = 199,
  G_DFF = 198,
  G_CMUX = 187,
  G_EQUIV = 188,
  G_SW = 189,
  G_TSD = 190,
  G_WIRE = 191,
};

static std::unordered_map<GType, std::string> GType2Name{
    {GType::G_PI, "INPUT"},  {GType::G_PO, "OUTPUT"}, {GType::G_AND, "AND"},
    {GType::G_NAND, "NAND"},
    {GType::G_OR, "OR"},     {GType::G_NOR, "NOR"},   {GType::G_NOT, "NOT"},
    {GType::G_BUF, "BUF"},   {GType::G_BRH, "BRH"},   {GType::G_XOR, "XOR"},
    {GType::G_XNOR, "XNOR"}, {GType::G_ABUF, "ABUF"}, {GType::G_TIE0, "TIE0"},
    {GType::G_TIE1, "TIE1"}, {GType::G_TIEX, "TIEX"}, {GType::G_TIEZ, "TIEZ"},
    {GType::G_MUX, "MUX"},   {GType::G_DLAT, "DLAT"}, {GType::G_DFF, "DFF"},
};

////////////////////////////////////////////////////
// primitives netlist type
enum class NetlistType {
  COMB_CIRCUIT,
  FULL_SCAN_CIRCUIT,
  PART_SCAN_CIRCUIT,
  SEQU_CIRCUIT,
  UNKNOWN_CIRCUIT,
};

static std::unordered_map<NetlistType, std::string> NetlistType2Name{
    {NetlistType::COMB_CIRCUIT, "comb"},
    {NetlistType::FULL_SCAN_CIRCUIT, "full_scan"},
    {NetlistType::SEQU_CIRCUIT, "seq"},
    {NetlistType::PART_SCAN_CIRCUIT, "partial_scan"},
    {NetlistType::UNKNOWN_CIRCUIT, "unknown type circuit"},
};
////////////////////////////////////////////////////
// logic value
enum class LogicValue {
  // atpg
  LOGIC_0 = 0,
  LOGIC_1,
  LOGIC_X,
  LOGIC_Z,
  LOGIC_S1,
  LOGIC_S0,
  LOGIC_G0,
  LOGIC_F0,
  LOGIC_G1,
  LOGIC_F1,
  LOGIC_CF,
  LOGIC_UJ,
  UNKNOWN_LOGIC,
};

enum class ConstraintValue {
  // input constraint value
  CT_FREE = 0, // constraint free
  CT_0,        // constraint 0
  CT_1,        // constraint 1
  CT_X,        // constraint X
  CT_N0,       // constraint non 0
  CT_N1,       // constraint non 1
};

enum class TieValue {
  TIE_FREE = 0,
  TIE_0, // gate is always 0 after learning
  TIE_1,
  TIE_X,
  TIE_Z,
};
////////////////////////////////////////////////////
// stuck-at fault type
enum class SAFType : std::uint8_t {
  SA0,
  SA1,
  UNKNOWN_TYPE = 2,
};

// stuck-at fault type to name
static std::unordered_map<SAFType, std::string> SAFType2Name{
    {SAFType::SA0, "sa0"},
    {SAFType::SA1, "sa1"},
};

// stuck-at fault status
enum class SAFStatus : std::uint8_t {
  // detected
  DS, // detected by simulation
  DI, // detected by implication
  DT_CHAIN_TEST,
  // possible
  PT, // possible testable
  PU, // possible untestable
  // atpg untestable
  AU,         //
  ATPG_ABORT, // ATPG abort
  // undetected
  UC, // uncontrolled
  UO, // unobserved
  // untestable
  UU, // unused
  TI, // tied  ,UNTESTABLE_TIED
  BL, // blocked
  RE, // redundant
  CONFLICT_FREE
};

// stuck-at fault status to tmax format
static std::unordered_map<SAFStatus, std::string> SAFStatus2Tmax{
    {SAFStatus::DS, "DS"},
    {SAFStatus::DI, "DI"},
    {SAFStatus::DT_CHAIN_TEST, "DT_CHAIN_TEST"},
    {SAFStatus::PT, "AP"},
    {SAFStatus::PU, "NP"},
    {SAFStatus::AU, "AN"},
    {SAFStatus::ATPG_ABORT, "NC"},
    {SAFStatus::UC, "NC"},
    {SAFStatus::UO, "NO"},
    {SAFStatus::UU, "UU"},
    {SAFStatus::TI, "UT"},
    {SAFStatus::BL, "UB"},
    {SAFStatus::RE, "UR"},
    {SAFStatus::CONFLICT_FREE, "CONFLICT_FREE"},
};

// stuck-at fault status to tessent format
static std::unordered_map<SAFStatus, std::string> SAFStatus2Tessent{
    {SAFStatus::DS, "DS"},
    {SAFStatus::DI, "DI"},
    {SAFStatus::DT_CHAIN_TEST, "DT_CHAIN_TEST"},
    {SAFStatus::PT, "PT"},
    {SAFStatus::PU, "PU"},
    {SAFStatus::AU, "AU"},
    {SAFStatus::ATPG_ABORT, "UC"},
    {SAFStatus::UC, "UC"},
    {SAFStatus::UO, "UO"},
    {SAFStatus::UU, "UU"},
    {SAFStatus::TI, "TI"},
    {SAFStatus::BL, "BL"},
    {SAFStatus::RE, "RE"},
    {SAFStatus::CONFLICT_FREE, "CONFLICT_FREE"},

};



////////////////////////////////////////////////////
// Learning Bool Table

static LogicValue AND[3][3] = {
    //  0        1        x
    {LogicValue::LOGIC_0, LogicValue::LOGIC_0, LogicValue::LOGIC_0}, // 0
    {LogicValue::LOGIC_0, LogicValue::LOGIC_1, LogicValue::LOGIC_X}, // 1
    {LogicValue::LOGIC_0, LogicValue::LOGIC_X, LogicValue::LOGIC_X}, // x
};
static LogicValue OR[3][3] = {
    //  0        1        x
    {LogicValue::LOGIC_0, LogicValue::LOGIC_1, LogicValue::LOGIC_X}, // 0
    {LogicValue::LOGIC_1, LogicValue::LOGIC_1, LogicValue::LOGIC_1}, // 1
    {LogicValue::LOGIC_X, LogicValue::LOGIC_1, LogicValue::LOGIC_X}, // x
};
static LogicValue XOR[3][3] = {
    //  0        1        x
    {LogicValue::LOGIC_0, LogicValue::LOGIC_1, LogicValue::LOGIC_X}, // 0
    {LogicValue::LOGIC_1, LogicValue::LOGIC_0, LogicValue::LOGIC_X}, // 1
    {LogicValue::LOGIC_X, LogicValue::LOGIC_X, LogicValue::LOGIC_X}, // x
};
static LogicValue NOT[3] = {
    //  0        1        x
    LogicValue::LOGIC_1, LogicValue::LOGIC_0, LogicValue::LOGIC_X};
static LogicValue BUF[3] = {
    //  0        1        x
    LogicValue::LOGIC_0, LogicValue::LOGIC_1, LogicValue::LOGIC_X};

// TODO
////////////////////////////////////////////////////
// SIM Define
#define PARALLEL_BIT_SIZE 64
#define V64_ZERO 0X0
#define V64_ONE 0XFFFFFFFFFFFFFFFF

#define MAX_FRAME_NUM 3
// frame id
#define FRAME_0 0
#define FRAME_1 1
#define FRAME_2 2

#define SIM_DEBUG_LEVEL                                                        \
  1 // 0 debug off ;1 report a few debug info ;2 report more debug info ...need
    // more define
// #define INSTRUMENT_TIMER //report time
} // namespace ictest

#endif // ICTEST_DEFINE_H
