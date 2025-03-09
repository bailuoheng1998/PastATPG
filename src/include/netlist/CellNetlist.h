// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_CELLNETLIST_H
#define ICTEST_CELLNETLIST_H
#include <string>
#include <unordered_map>
#include <vector>

#include "Cell.h"
#include "PrimNetlist.h"

namespace ictest {

class CellNetlist {
 public:
  CellNetlist() {}

  ~CellNetlist() { ClearCellNetlist(); }

  void ClearCellNetlist();
  void InitCellNetlist();
  bool ReadCellVy(const std::string& vy_path);
  void ConnectCellNetlist(PrimNetlist* prim_netlist);
  void SetCellBorder();
  bool CheckCellNetlistValidation();

  inline int32_t NumCells() const { return num_cells_; }
  std::vector<Cell*>& GetCellNetlist() { return cell_netlist_; }
  std::unordered_map<int32_t, std::string>& GetNet2Names() {
    return net2names_;
  }
  std::unordered_map<int32_t, std::string>& GetPin2Names() {
    return pin2names_;
  }
  std::unordered_map<int32_t, std::string>& GetCType2Names() {
    return ctype2names_;
  }
  std::unordered_map<int32_t, int32_t>& GetOutNet2Cell() {
    return outnet2cell_;
  }
  std::unordered_map<int32_t, int32_t>& GetInNet2Cell() { return innet2cell_; }
  std::unordered_map<std::string, Cell*>& GetCellName2Cell() {
    return cellname2cell_;
  }

  void DumpCellNetlist(const std::string& file_path);

 private:
  int32_t num_cells_{0};
  std::vector<Cell*> cell_netlist_;
  std::unordered_map<int32_t, std::string> net2names_;
  std::unordered_map<int32_t, std::string> pin2names_;
  std::unordered_map<int32_t, std::string> ctype2names_;
  std::unordered_map<int32_t, int32_t> outnet2cell_;
  std::unordered_map<int32_t, int32_t> innet2cell_;
  // cell instance name to cell ptr
  std::unordered_map<std::string, Cell*> cellname2cell_;
};
}  // namespace ictest
#endif  // ICTEST_CELLNETLIST_H
