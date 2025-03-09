// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "netlist/CellNetlist.h"

#include <fstream>
#include <sstream>

#include "util/Log.h"

namespace ictest {

void CellNetlist::ClearCellNetlist() {
  for (auto& cell : cell_netlist_) {
    if (cell != nullptr) {
      delete cell;
      cell = nullptr;
    }
  }
  cell_netlist_.clear();
  net2names_.clear();
  pin2names_.clear();
  ctype2names_.clear();
  outnet2cell_.clear();
  innet2cell_.clear();
  cellname2cell_.clear();
}

void CellNetlist::InitCellNetlist() {
  // cell data
  num_cells_ = 0;
  cell_netlist_.clear();
  net2names_.clear();
  pin2names_.clear();
  ctype2names_.clear();
  outnet2cell_.clear();
  innet2cell_.clear();
  cellname2cell_.clear();
}

bool CellNetlist::ReadCellVy(const std::string& vy_path) {
  std::ifstream fin(vy_path);
  if (!fin.is_open()) {
    LOG_ERROR("ERROR:cannot open cell vy file:" + vy_path);
    return false;
  }
  std::vector<std::string> vstr;
  std::string line;
  while (fin.good()) {
    getline(fin, line);
    vstr.push_back(line);
  }
  int count = 0;
  int count_cell = 0;
  std::string str;
  for (int i = 7; i < vstr.size(); i++) {
    std::vector<std::string> elems;
    std::istringstream ifline(vstr[i]);
    while (ifline >> str) {
      elems.push_back(str);
    }
    if (elems.size() > 2) {
      int m = 2;  // out id start from 2: 200 ==U9 1 9 0 0 3 10 1 11 2 2 3
      int cell_logic_id = stoi(elems[0]);
      std::string inst_name = elems[1];
      Cell* cell = new Cell(count_cell, cell_logic_id, inst_name);
      ++count_cell;
      int num_out = stoi(elems[m]);
      for (int i = 0; i < num_out; ++i) {
        int net_id = stoi(elems[++m]);
        int pin_id = stoi(elems[++m]);
        cell->GetOutputNets().push_back(net_id);
        outnet2cell_[net_id] = count_cell - 1;
        cell->GetOutputPins().push_back(pin_id);
      }
      int num_inout = stoi(elems[++m]);
      for (int i = 0; i < num_inout; ++i) {
        int net_id = stoi(elems[++m]);
        int pin_id = stoi(elems[++m]);
        cell->GetInoutNets().push_back(net_id);
        cell->GetInoutPins().push_back(pin_id);

        innet2cell_[net_id] = count_cell - 1;
      }
      int num_in = stoi(elems[++m]);
      for (int i = 0; i < num_in; ++i) {
        int net_id = stoi(elems[++m]);
        int pin_id = stoi(elems[++m]);
        cell->GetInputNets().push_back(net_id);
        cell->GetInputPins().push_back(pin_id);
      }
      cell_netlist_.push_back(cell);
    }


    else if (elems.size() > 0 && elems.size() == 2) {
      int id = stoi(elems[0]);
      std::string name = elems[1];
      if (count == 1) {
        net2names_[id] = name;
      } else if (count == 2) {
        pin2names_[id] = name;
      } else if (count == 3) {
        ctype2names_[id] = name;
      }
    } else {
      ++count;
    }
  }

//  for(auto& cell_ptr:cell_netlist_){
//    std::cout<<cell_ptr->GetInstanceName()<<"-------------------------"<<std::endl;
//    assert(cell_ptr->GetInputNets().size()==cell_ptr->GetInputPins().size());
//    for(int i=0;i<cell_ptr->GetInputNets().size();++i){
//      int net_id=cell_ptr->GetInputNets()[i];
//      int pin_id=cell_ptr->GetInputPins()[i];
//      std::cout<<net2names_[net_id]<<" + "<<pin2names_[pin_id]<<" ";
//    }
//    std::cout<<"\n";
//    assert(cell_ptr->GetOutputNets().size()==cell_ptr->GetOutputPins().size());
//    for(int i=0;i<cell_ptr->GetOutputNets().size();++i){
//      int net_id=cell_ptr->GetOutputNets()[i];
//      int pin_id=cell_ptr->GetOutputPins()[i];
//      std::cout<<net2names_[net_id]<<" + "<<pin2names_[pin_id]<<" ";
//    }
//    std::cout<<"\n";
//  }

  num_cells_ = count_cell;
  fin.close();
  return true;
}

void CellNetlist::ConnectCellNetlist(PrimNetlist* prim_netlist) {
  // add fanouts
  for (int i = 0; i < num_cells_; ++i) {
    Cell* tmp_cell = cell_netlist_[i];

    LOG_ASSERT(tmp_cell, "tmp_cell nullptr");
    cellname2cell_[tmp_cell->GetInstanceName()] = tmp_cell;

    std::vector<int> out_nets_id = tmp_cell->GetOutputNets();

    for (int j = 0; j < out_nets_id.size(); j++) {
      // find and add fanouts cell if this net is a fanin
      if (innet2cell_.find(out_nets_id[j]) != innet2cell_.end()) {
        int suc_cell_id = innet2cell_[out_nets_id[j]];
        LOG_ASSERT(cell_netlist_[suc_cell_id], "gate nullptr");

        tmp_cell->GetFanoutCells().push_back(cell_netlist_[suc_cell_id]);
      }

      Gate* out_gate = prim_netlist->GetPrimNetlist()[out_nets_id[j]];
      LOG_ASSERT(out_gate, "out_gate nullptr");
      // tmp_cell->AddFanoutGates(out_gate);
      tmp_cell->GetFanoutGates().push_back(out_gate);
    }
  }
  // add fanin
  for (int i = 0; i < num_cells_; ++i) {
    Cell* tmp_cell = cell_netlist_[i];
    LOG_ASSERT(tmp_cell, "tmp_cell nullptr");
    std::vector<int> in_nets_id = tmp_cell->GetInputNets();
    std::vector<int> in_pins_id = tmp_cell->GetInputPins();

    for (int j = 0; j < in_nets_id.size(); j++) {
      // find and add fanin cell if this net is a fanouts
      if (outnet2cell_.find(in_nets_id[j]) != outnet2cell_.end()) {
        int tmp_cell_id = outnet2cell_[in_nets_id[j]];
        LOG_ASSERT(cell_netlist_[tmp_cell_id], "cell nullptr");
        tmp_cell->GetFaninCells().push_back(cell_netlist_[tmp_cell_id]);
      }
      // find and add fanin gate if this net is cellboader
      Gate* in_gate = prim_netlist->GetPrimNetlist()[in_nets_id[j]];
      LOG_ASSERT(in_gate, "in_gate nullptr");

      if (in_gate->FanoutSize() <= 1) {
        tmp_cell->GetFaninGates().push_back(in_gate);
      } else {
        int cell_pin_id = in_pins_id[j];
        Gate* branch = prim_netlist->GetPin2BranchGates()[cell_pin_id];
        // why??? branch == null, which means pin's gate is a stem but do not
        // have branch?
        if (branch != nullptr) {  // actually FaninGate can be nullptr
          tmp_cell->GetFaninGates().push_back(branch);
        }
      }
    }
  }
}

void CellNetlist::SetCellBorder() {
  for (auto cell : cell_netlist_) {
    LOG_ASSERT(cell, "cell nullptr");

    std::vector<Gate*> fanins = cell->GetFaninGates();
    std::vector<Gate*> fanouts = cell->GetFanoutGates();

    for (auto gptr : fanins) {
      gptr->SetCellBorder(true);
      if (gptr->FanoutSize() > 1) {
        for (auto branch : gptr->FanoutGates()) {
          branch->SetCellBorder(true);
        }
      }
    }

    for (auto gptr : fanouts) {
      gptr->SetCellBorder(true);
      if (gptr->FanoutSize() > 1) {
        for (auto branch : gptr->FanoutGates()) {
          branch->SetCellBorder(true);
        }
      }
    }
  }
}

// todo
bool CellNetlist::CheckCellNetlistValidation() { return false; }

void CellNetlist::DumpCellNetlist(const std::string& file_path) {
  std::ofstream ofs(file_path, std::ios::out);
  ofs << "Cellist "
         "Data============================================================"
      << std::endl;

  for (int index = 0;
       index < cell_netlist_.size() && cell_netlist_[index] != nullptr;
       index++) {
    Cell* cell = cell_netlist_[index];
    LOG_ASSERT(cell, "cell is nullptr");
    ofs << "Id:\t" << cell->GetCId() << " " << std::endl;
    ofs << "\tLogicType:\t" << cell->GetCType() << " " << std::endl;
    ofs << "\tName:\t" << cell->GetInstanceName() << " " << std::endl;
    ofs << "\tFanin Cell Size:\t" << cell->GetFaninCells().size() << " "
        << std::endl;
    ofs << "\tFanout Cell Size:\t" << cell->GetFanoutCells().size() << " "
        << std::endl;
    ofs << "\tFanin Gate Size:\t" << cell->GetFaninGates().size() << " "
        << std::endl;
    ofs << "\tFanout Gate Size:\t" << cell->GetFanoutGates().size() << " "
        << std::endl;
    ofs << "\tFanin Faults Size:\t" << cell->GetFaninSAFs().size() << " "
        << std::endl;
    ofs << "\tFanout Faults Size:\t" << cell->GetFanoutSAFs().size() << " "
        << std::endl;

    int output_size = cell->GetOutputNets().size();
    ofs << "\toutput_size:\t" << output_size << " " << std::endl;
    ofs << "\t\toutput:\t";
    for (int i = 0; i < output_size; i++) {
      ofs << cell->GetOutputNets()[i] << " ";
      ofs << cell->GetOutputPins()[i] << " ";
    }
    ofs << std::endl;

    int inout_size = cell->GetInoutNets().size();
    ofs << "\tinout_size:\t" << inout_size << " " << std::endl;
    ofs << "\t\tinout:\t";
    for (int i = 0; i < inout_size; i++) {
      ofs << cell->GetInoutNets()[i] << " ";
      ofs << cell->GetInoutPins()[i] << " ";
    }
    ofs << std::endl;

    int input_size = cell->GetInputNets().size();
    ofs << "\tinput_size:\t" << input_size << " " << std::endl;
    ofs << "\t\tinput:\t";
    for (int i = 0; i < input_size; i++) {
      ofs << cell->GetInputNets()[i] << " ";
      ofs << cell->GetInputPins()[i] << " ";
    }

    ofs << std::endl;
    ofs << "-------------------------------------------------------------------"
           "-----"
        << std::endl;
  }
}

}  // namespace ictest