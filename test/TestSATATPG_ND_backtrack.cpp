// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#include "atpg/SATEngine.h"
#include "drc/DRC.h"
#include "fault/SAFList.h"
#include "gtest/gtest.h"
#include "netlist/CellNetlist.h"
#include "netlist/PrimNetlist.h"
#include "pattern/Pattern.h"
#include "simulation/SimEngine.h"
#include "stil/STILParser.h"
#include "util/Log.h"

// Global namespace
namespace global {
std::string run_mode;
std::string prim_vy;
std::string cell_vy;
std::string spf;
std::string external_fl;
std::string dump_path;
}  // namespace global

// Global command line parameters.
class MyTestEnvironment : public testing::Environment {
 public:
  MyTestEnvironment(int argc, char** command_line_arg) {
    if (argc != 5 && argc != 6) {
      std::cout << argc << std::endl;
      std::cout << "Error! invalid numbers of params" << std::endl;
      assert(false);
    }
    if (command_line_arg[1] != nullptr) {
      int begin_idx = 1;
      global::prim_vy = command_line_arg[begin_idx];
      global::cell_vy = command_line_arg[++begin_idx];
      global::spf = command_line_arg[++begin_idx];
      if (argc == 6)global::external_fl = command_line_arg[++begin_idx];
      global::dump_path = command_line_arg[++begin_idx];
    }
  }
};


void RunSingleSAF(const std::string prim_vy_path,
                  const std::string cell_vy_path,
                  const std::string spf_path,
                  const std::string dump_path,
                  const std::string external_fl) {
  ictest::Log log;
  log.Init();

  SPDLOG_INFO("# start SAF ATPG test ");

  // 1. primitives netlist flow
  ictest::PrimNetlist* pnlist = new ictest::PrimNetlist();
  pnlist->InitPrimNetlist();
  pnlist->ReadPrimVY(prim_vy_path);
  pnlist->ConnectPrimNetlist();
  pnlist->AddPOGates();
  pnlist->AddBranchGates();
  pnlist->Levelize();
  pnlist->AddPin2BranchGates();
  pnlist->ResizePrimNetlist();
  pnlist->CheckPrimNetlistValidation();


  // 2. cell netlist flow
  ictest::CellNetlist* cnlist = new ictest::CellNetlist();
  cnlist->InitCellNetlist();
  cnlist->ReadCellVy(cell_vy_path);
  cnlist->ConnectCellNetlist(pnlist);
  cnlist->SetCellBorder();
  cnlist->CheckCellNetlistValidation();


  // 3. spf flow
  ictest::STILParser* stil_parser = new ictest::STILParser();
  stil_parser->ReadICTestSPF(spf_path);
  ictest::DRC* drc = new ictest::DRC(pnlist, stil_parser);
  //drc->RunDRC();
  pnlist->AnalyseNetlistType();

  // 4. init gate status
  pnlist->SetGate2Status();

  // 5. fault list flow
  ictest::SAFList* saf_flist = new ictest::SAFList();
  saf_flist->SetupNetlist(pnlist, cnlist);
  saf_flist->SetGate2SAFs();

  // 5.1 generate internal fault
  saf_flist->GenerateInternalSAFs();
    if (global::external_fl == "") {
        saf_flist->SetInternalSAFs();
    } else {
        // 5.2 generate external fault
        saf_flist->GenerateExternalSAFs(external_fl);
        saf_flist->SetExternalSAFs();
    }



  saf_flist->SetCountSAFs();
  saf_flist->FindDISAFs();

  // 6. set option
  ictest::Option* options = new ictest::Option();
  options->fault_model_ = "stuck-at";
  options->flist_source_ = "internal";
  options->sim_engine_type_ = "basic_scan";
  options->pattern_per_sim_ = 64;
  options->abort_limit_ = 10;

  // 7. set simulator
  ictest::SimEngine* sim = new ictest::SimEngine(pnlist, *options);

  // 8.stuck-at fault ATPG.
  // sat atpg
  auto* test_engine = new ictest::SATEngine(pnlist, saf_flist, sim, dump_path);

//  test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompactionDynamic_2();
    test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompaction_czt(1);
   //compaction_dynamic
    //test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompactionDynamic();
  // kissat
//  test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompareKis();
  //test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompareKis_dynamic();
  // minisat
//  test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompareMini();
  // cadical
  //test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompareCal();

  //test_engine->GenerateCompactTestSetForSAFListHard2DTfaultsCompareCal_dynamic();

  SPDLOG_INFO("# end test ");
  log.ShutDown();
}

void runFault() {
  RunSingleSAF(global::prim_vy, global::cell_vy, global::spf, global::dump_path, global::external_fl);
}
TEST(TestATPG, TestFault) { runFault(); }

int main(int argc, char** argv) {
  std::cout << "Start test all ATPG" << std::endl;
  testing::InitGoogleTest(&argc, argv);
  auto const& arg =
      testing::AddGlobalTestEnvironment(new MyTestEnvironment(argc, argv));
  return RUN_ALL_TESTS();
}