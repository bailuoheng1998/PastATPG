// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "netlist/Gate.h"
namespace ictest {

LogicValue Gate::GetCtrInVal() const {
  switch (g_type_) {
    case GType::G_AND:
    case GType::G_NAND:
      return LogicValue::LOGIC_0;
    case GType::G_OR:
    case GType::G_NOR:
      return LogicValue::LOGIC_1;
    default:
      break;
  }

  return LogicValue::LOGIC_X;
}

LogicValue Gate::GetCtrOutVal() const {
  switch (g_type_) {
    case GType::G_AND:
    case GType::G_NOR:
      return LogicValue::LOGIC_0;
    case GType::G_NAND:
    case GType::G_OR:
      return LogicValue::LOGIC_1;
    default:
      break;
  }

  return LogicValue::LOGIC_X;
}

LogicValue Gate::GetNonCtrInVal() const {
  switch (g_type_) {
    case GType::G_AND:
    case GType::G_NAND:
      return LogicValue::LOGIC_1;
    case GType::G_OR:
    case GType::G_NOR:
      return LogicValue::LOGIC_0;
    default:
      break;
  }

  return LogicValue::LOGIC_X;
}

LogicValue Gate::GetNonCtrOutVal() const {
  switch (g_type_) {
    case GType::G_AND:
    case GType::G_NOR:
      return LogicValue::LOGIC_1;
    case GType::G_NAND:
    case GType::G_OR:
      return LogicValue::LOGIC_0;
    default:
      break;
  }

  return LogicValue::LOGIC_X;
}

}  // namespace ictest