// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_SIMDATA_H
#define ICTEST_SIMDATA_H
#include "SimVal.h"
namespace ictest {
struct scan_info{
  int scan_lenth_;
  int scan_si;
  int scan_so;
  int scan_en;
  int scan_master_clock;
};

enum class GateClassify: std::uint8_t {
  UNKNOW,
  TM_SIG,
  SE_SIG,
  CLK_SIG,
  TIE_SIG,
  TM_MUX,
  SE_MUX,
  SET_RESET_SIG,
  ONLY_SI_PATH,//with out DFF MUX
  RISING_SC_DFF,//Define clk 1->0->1 falling 0->1->0 rising in this clk
  FALLING_SC_DFF,
  SCAN_DLATCH,
  ACT_RISING_NSC,
  ACT_FALLING_NSC,
  MIXED,
};

struct SimData {
 public:
  SimData()=default;
  ~SimData()=default;

  //scan chain
  void SetChainSo(int chain_id,int so_gid){
    si_so_port_[chain_id].second=so_gid;
  }

  void SetChainSi(int chain_id,int si_gid){
    si_so_port_[chain_id].first=si_gid;
  }

  int GetChainSo(int sc_id)const{
    return si_so_port_[sc_id].second;
  }

  int GetChainSi(int chain_id){
    return si_so_port_[chain_id].first;
  }
  void AddShiftPulseClks(int gid ,val64_t init_val){
    shift_clk_off_state_.emplace_back(gid,init_val);
  }

   int GetChainNums()const{
    return chain_nums_;
  }

  int GetChainSi(int sc_id)const{
    return si_so_port_[sc_id].first;
  }

  int GetMaxChainLenth()const{
    return max_chain_lenth_;
  }

  void AddTmSigs(int gid ,val64_t init_val){
    tm_sigs_.emplace_back(gid,init_val);
  }
  void AddSESigs(int gid , val64_t init_val){
    se_sigs_.emplace_back(gid,init_val);
  }

  std::vector<std::pair<int,int>>& GetAllChainGateSeg(){
    return chain_gate_2_sigment_;
  }
 public:
  // so<-buf<-dff<buf<-dff<-buf<-si  pair<chain_id,seg_idx>
  // 0    0    0   1    1    2   2
  std::vector<std::pair<int,int>> chain_gate_2_sigment_;
  std::vector<std::pair<int,int>> clk_to_pi_id_;
  std::vector<std::pair<int,int>> si_so_port_;
  std::vector<GateClassify> gate_classify_mark_;
  /////////////////////////////////////////////////////
  //signals pair==>  gate-id: val(off state)
  std::vector<std::pair<int,val64_t>> shift_clk_off_state_;
  std::vector<std::pair<int,val64_t>> tm_sigs_;
  std::vector<std::pair<int,val64_t>> se_sigs_;
  std::vector<std::pair<int,val64_t>> set_reset_sigs_;
  std::vector<std::pair<int,val64_t>> tie_sigs_;
  //chain_id ->val_idx in sc_pat->nsc_dff_id->inverter
  std::vector<std::tuple<int,int,int,bool>> shadow_dff;
  int max_chain_lenth_ =0;
  int chain_nums_ =0;
  int rising_edge_sc_dff_nums{0};
  int falling_edge_sc_dff_nums{0};
  int act_rising_edge_nsc_dff_nums{0};
  int act_falling_edge_nsc_dff_nums{0};
};
}  // namespace ictest
#endif  // ICTEST_SIMDATA_H
