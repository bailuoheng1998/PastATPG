//
// Created by tianpengyu on 23-4-4.
//

#include "atpg/SATEngine.h"

namespace ictest {


void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompactionDynamic(int pattern_per_sim) {

        vector<int>  confirmed_bits;

        double pattern_xbits_rate_sum = 0;
        int pattern_count = 0;

        int test_found = 0, no_test = 0, abort = 0;
        int num_curr_pattern = 0;
        int total_pattern_num = 0;
        double sim_time = 0;
        double cnf_time = 0;
        double cnf_trans_time = 0;
        double solve_time = 0;
        int num_faults = saf_list_->GetUncollapsedSAFs().size();
        auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
        vector<std::pair<int, double>> single_time;
        std::unordered_set<SAF *> success_targets;
        std::vector<SAF *> rest_faults;
        std::vector<SAF *> success_pfs;
        for (auto fault : uncollapsed_saf_list) {
            if (fault->GetSAFStatus() == SAFStatus::UC) {
                rest_faults.emplace_back(fault);
            }
        }
        int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
        std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
        CNFGenerator saf_test_generator(prim_);

        std::vector<Pattern> patterns;
        std::vector<Pattern> total_patterns;
        std::vector<std::vector<int>> atpg_patterns;
        std::vector<std::vector<int>> validation_patterns;
        std::vector<int> new_pattern(prim_->NumPIs());
        PatternParser pp;
        pp.SetupPrimNetlist(prim_);
        pp.SetPatternType(PatternType::COMB_PT);
        assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);


        vector<vector<int>> compaction_patterns;

        vector<int> check_fault(num_faults,0);
        vector<int> new_detected_fault;

        vector<int> Necessary_fault_assignment(prim_->NumPIs(),0);


        auto start = std::chrono::steady_clock::now();


        //vector<double> time_for_each_fault;

        for (int pf_id = 0; pf_id < num_rest_faults; ++pf_id) {
            vector<int> determinate_bits;
            //int last_fault;
            auto* primary_fault = rest_faults[pf_id];
            /*if (pf_id > 0){
                last_fault = rest_faults[pf_id - 1]->GetSAFId();
            }*/

            if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
               && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
                continue;
            }

            /*std::cout << "At Fault: " << primary_fault->GetSAFName() << " and id is"<<primary_fault->GetSAFId()<<" "
                      << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

            std::cout<<"fault gate is : "<<primary_fault->GetSAFGate()->GetGId()<<" and type is : "<<int(primary_fault->GetSAFType())<<std::endl;*/

            auto solver = new Solver;
            solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
            solver->max_dpi_ = prim_->MaxLevel();
            solver->fault_site_ = primary_fault->GetSAFGate();
            saf_test_generator.SetSolver(solver);

            int test_res;
            int test_res_compaction;
            std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
            //auto cnf_start = std::chrono::steady_clock::now();

            if (saf_test_generator.CollectDChainCircuitClauseSAFL2TPI(primary_fault,
                                                                      fanin_mark)){
                //auto cnf_end = std::chrono::steady_clock::now();
                //cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
//            test_res = solver->solve(d_test_engine_, primary_fault);
                test_res = solver->solve_help_checkx();

                int xpi = 0;
                for (auto& xval: solver->x_vars_) {
                    if(xval < prim_->GetPrimNetlist().size() && prim_->GetPrimNetlist()[xval]->GetGType() == GType::G_PI) {
                        xpi++;
                    }
                }

            }
            else {
                test_res = 20;
            }

            if (test_res == 20) {
                primary_fault->SetSAFStatus(SAFStatus::RE);
                no_test++;
                //      std::cout << "No Test" << std::endl;
            }

            else if (test_res == 0) {
                primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
            }

            else if (test_res == 10){
                // Load test cube from ATPG circuit model.
                int determined_bits = 0;
                int x_bits = 0;


                for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                    auto pi_gate = prim_->GetPIGates()[pi_id];
                    if (!fanin_mark[pi_gate->GetGId()] ) {
                        new_pattern[pi_id] = LOGIC_x;
                        x_bits++;


                        //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                        //          ";
                        continue;
                    }


                    //add to determinate bits
                    else{
                        determinate_bits.emplace_back(pi_gate->GetGId());
                    }


                    assert(pi_gate->GetGType() == GType::G_PI);
                    auto pi_bit =
                            (solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
                    //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
                    determined_bits++;
                    //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
//                assert(pi_bit != LOGIC_x);
                    new_pattern[pi_id] = pi_bit;
                }





            //start dynamiic compaction

            /*for (int i = pf_id + 1; i < std::min(num_rest_faults,20); ++i) {

                CNFGenerator saf_test_generator_dynamic(prim_);

                auto *primary_fault = rest_faults[i];

                if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
                    && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
                    continue;
                }

                auto solver_2 = new Solver;
                solver_2->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
                solver_2->max_dpi_ = prim_->MaxLevel();
                solver_2->fault_site_ = primary_fault->GetSAFGate();
                saf_test_generator_dynamic.SetSolver(solver_2);

                std::vector<int> fanin_mark_2(prim_->GetPrimNetlist().size(), 0);

                int compaction_success = 0;
                if (saf_test_generator_dynamic.CollectDChainCircuitClauseSAFL2TPI(primary_fault,
                                                                          fanin_mark_2)) {


                    for (int j = 0; j < prim_->GetPIGates().size(); j++){
                        if (new_pattern[j] == 0){
                            solver_2->add(-1 * (prim_->GetPIGates()[j]->GetGId() + 1));
                            solver_2->add(0);
                        }else if (new_pattern[j] == 1){
                            solver_2->add(prim_->GetPIGates()[j]->GetGId() + 1);
                            solver_2->add(0);
                        }
                    }

                    *//*for (int j = 0; j < new_pattern.size(); ++j) {
                        if (new_pattern[j] == 1){
                            solver_2->add(j + 1);
                            solver_2->add(0);
                        }else if (new_pattern[j] == 0){
                            solver_2->add(-1*(j + 1));
                            solver_2->add(0);
                        }

                    }*//*

                    test_res_compaction = solver_2->solve_help_checkx();
                    //test_res = solver->solve();


                    int xpi = 0;
                    for (auto& xval: solver_2->x_vars_) {
                        if(xval < prim_->GetPrimNetlist().size() && prim_->GetPrimNetlist()[xval]->GetGType() == GType::G_PI) {
                            xpi++;
                        }
                    }

                    if (test_res_compaction == 10){
                        //compaction_success = 1;
                        //std::cout<<"dynamic compaction success!"<<std::endl;


                        *//*for (int i = 0; i < new_pattern.size(); ++i) {
                            if (new_pattern[i] != LOGIC_x){
                              fanin_mark_2[i] = 1;
                            }
                        }*//*

                    //int determined_bits_2 = 0;
                    //int x_bits_2 = 0;

                        for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                            auto pi_gate = prim_->GetPIGates()[pi_id];
                            *//*if (!fanin_mark_2[pi_gate->GetGId()] ) {
                                new_pattern[pi_id] = LOGIC_x;
                                //x_bits_2++;
                                //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                                //          ";
                                continue;
                            }*//*
                            assert(pi_gate->GetGType() == GType::G_PI);
                            auto pi_bit =
                                (solver_2->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver_2->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
                            //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
                            //determined_bits_2++;
                            //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
//                  assert(pi_bit != LOGIC_x)
                            if (new_pattern[pi_id] == LOGIC_x){
                                //new_pattern[pi_id] = pi_bit;
                            }

                        }
                       *//* std::cout<<"new pattern is:";
                        for (int j = 0; j < new_pattern.size(); ++j) {
                            std::cout<<" "<<new_pattern[j]<<" ";
                        }
                        std::cout<<std::endl;*//*

                    }


                }
                //fanin_mark_2.clear();
                delete solver_2;
            }*/
            }



            delete solver;

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);
            success_pfs.emplace_back(primary_fault);

            success_targets.insert(primary_fault);

            num_curr_pattern = static_cast<int>(atpg_patterns.size());
            if (num_curr_pattern) {
                //auto sim_start = std::chrono::steady_clock::now();
                total_pattern_num += num_curr_pattern;
                pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
                sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
                atpg_patterns.clear();
                total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
                patterns.clear();
                //auto sim_end = std::chrono::steady_clock::now();
                //sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
            }

        }


        /*std::cout<<"time for each fault is: ";
        for (int i = 0; i < time_for_each_fault.size(); ++i) {
            std::cout<<time_for_each_fault[i]<<std::endl;
        }*/
        //outputFile.close();



        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration<double>(end - start);

        int fault_detect = 0;
        int fault_redundant = 0;
        int fault_not_detect = 0;
        double fault_coverage = 0;
        double test_coverage = 0;
        DumpSAFList();
        for (auto fault : saf_list_->GetUncollapsedSAFs()) {
            if (fault->GetSAFStatus() == SAFStatus::DS ||
                fault->GetSAFStatus() == SAFStatus::DI ||
                fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
                fault_detect++;
            } else if (fault->GetSAFStatus() == SAFStatus::RE) {
                fault_redundant++;
            } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                       fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                       fault->GetSAFStatus() == SAFStatus::UO) {
                fault_not_detect++;
            }
        }
        fault_coverage =
                (double)fault_detect / (double)uncollapsed_saf_list.size();
        test_coverage =
                (double)(fault_detect + fault_redundant) /
                (double)uncollapsed_saf_list.size();

        std::cout << "Detected Faults: " << fault_detect << std::endl;
        std::cout << "No Test Faults: " << fault_redundant << std::endl;
        std::cout << "Abort Faults: " << fault_not_detect << std::endl;
        std::cout << "Pattern generated: " << total_pattern_num << std::endl;
        std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
                  << std::endl;
        std::cout << "Test Coverage: " << test_coverage * 100 << "%"
                  << std::endl;
        std::cout << "Test time: " << duration.count() << "s." << std::endl;
        std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
        std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
        std::cout << "Solver time: " << solve_time << "s." << std::endl;
        std::cout << "Sim time: " << sim_time << "s." << std::endl;


        std::cout<<"Average x_bits fill rate is used fanin : "<<pattern_xbits_rate_sum/pattern_count<<std::endl;


        // ==================== VALIDATION ===========================
        std::cout << std::endl;
        std::cout << "----------Validate test patterns-----------" << std::endl;
        for (auto fault : uncollapsed_saf_list) {
            fault->SetSAFStatus(SAFStatus::UC);
        }
        pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
        sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
        validation_patterns.clear();
        patterns.clear();

        /*for (int i = 0; i < uncollapsed_saf_list.size(); ++i) {
            if(int(uncollapsed_saf_list[i]->GetSAFStatus()) != 0){
                std::cout<<"undetected fault id is :"<<i<<std::endl;
                std::cout<<"type is :"<<int(uncollapsed_saf_list[i]->GetSAFStatus())<<std::endl;
            }
        }*/

        // ==================== VALIDATION END ==============
        // ==================== CHECK X BIT ===========================
        std::cout << std::endl;
        /*std::cout << "----------Check X bit in test patterns-----------" << std::endl;
        assert(success_pfs.size() == total_patterns.size());
        for (int i = 0; i < total_patterns.size(); ++i) {
            auto& pat_cycle_pis = total_patterns[i].GetPatternCycle(0).GetInputVal();
            assert(pat_cycle_pis.size() == 1);
            val64_t ori_val;
            int add_x_bit = 0;
            int ori_nox_bit = 0;
            for (int j = 0; j < pat_cycle_pis[0].size(); ++j) {
                auto& pi_v = pat_cycle_pis[0][j];
                val64_t x_v (0, 0);
                if (pi_v != x_v) {
                    ori_nox_bit++;
                    ori_val = pi_v;
                    pi_v.SetX();
                    std::vector<Pattern> one_pat;
                    one_pat.push_back(total_patterns[i]);
                    std::vector<SAF*> one_pf;
                    success_pfs[i]->SetSAFStatus(SAFStatus::UC);
                    one_pf.push_back(success_pfs[i]);
                    sim_->RunSAFSimulation(one_pat, one_pf, false);
                    if (one_pf[0]->GetSAFStatus() != SAFStatus::DS) {
                        pi_v = ori_val;
                        one_pat.clear();
                        one_pat.push_back(total_patterns[i]);
                        sim_->RunSAFSimulation(one_pat, one_pf, false);
                        assert(one_pf[0]->GetSAFStatus() == SAFStatus::DS);
                    } else {
                        add_x_bit++;
                    }
                }
            }
            std::cout << "pat " << i << ": add_x: " << add_x_bit <<"; total: " << ori_nox_bit << std::endl;
        }*/
        // ==================== CHECK X BIT END ==============
//    // ==================== CHECK OVERLAP ===========================
//    std::cout << std::endl;
//    std::cout << "----------Check overlap in test patterns-----------" << std::endl;
//    int overlap_bit = 0;
//    vector<int> overlap_bits;
//    for (int i = 0; i < total_patterns.size() - 1; ++i) {
//        for (int j = i + 1; j < total_patterns.size(); ++j) {
//            auto& pat_cycle_pis1 = total_patterns[i].GetPatternCycle(0).GetInputVal();
//            auto& pat_cycle_pis2 = total_patterns[j].GetPatternCycle(0).GetInputVal();
//            assert(pat_cycle_pis1.size() == pat_cycle_pis2.size());
//            for (int k = 0; k < pat_cycle_pis1[0].size(); ++k) {
//                auto& pi_v1 = pat_cycle_pis1[0][k];
//                auto& pi_v2 = pat_cycle_pis2[0][k];
//                if (pi_v1 != pi_v2) {
//                    overlap_bit++;
//                }
//            }
//            overlap_bits.push_back(overlap_bit);
//            overlap_bit = 0;
//        }
//    }
//    double avg = std::accumulate(overlap_bits.begin(), overlap_bits.end(), 0.0);
//    avg = avg / overlap_bits.size();
//    std::cout << "avg overlap: " << avg << std::endl;
//    std::sort(overlap_bits.begin(), overlap_bits.end());
//    double min_ = overlap_bits[0];
//    double max_ = overlap_bits.back();
//    std::cout << "max overlap: " << max_ << std::endl;
//    std::cout << "min overlap: " << min_ << std::endl;
//    std::cout << "mean overlap: " << overlap_bits[overlap_bits.size()/2] << std::endl;
//    // ==================== CHECK OVERLAP END ==============
        std::sort(single_time.begin(), single_time.end(),
                  [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                      return a.second > b.second;
                  });

        // 打印前十个最大的项
        for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
            std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
        }


       /* std::cout<<"compaction_pattern_size is: "<<compaction_patterns.size()<<std::endl;
        //start compaction step_1
        vector<vector<int>> new_compacted_patterns;
        new_compacted_patterns.clear();
        vector<int> mid_pattern;
        vector<int> fanin_mark(prim_->NumPIs(),0);
        vector<int> pattern_compcted_mark(compaction_patterns.size(),0);
        for (int i = 0; i < compaction_patterns.size(); ++i) {
            if (pattern_compcted_mark[i] == 0){
                pattern_compcted_mark[i] = 1;
                for (int j = 0; j < compaction_patterns[i].size(); ++j) {
                    mid_pattern.emplace_back(compaction_patterns[i][j]);
                }
                for (int j = i + 1; j < compaction_patterns.size(); ++j) {
                    int conflict_flag = 0;
                    for (int k = 0; k < compaction_patterns[j].size(); ++k) {
                        if (mid_pattern[k] != compaction_patterns[j][k] && mid_pattern[k] != 2 && compaction_patterns[j][k] != 2){
                            conflict_flag = 1;
                        }
                    }
                    if (conflict_flag == 0){
                        pattern_compcted_mark[j] = 1;
                        for (int k = 0; k < mid_pattern.size(); ++k) {
                            if (mid_pattern[k] == 2){
                                mid_pattern[k] = compaction_patterns[j][k];
                            }
                        }
                    }

                }
            }
            if (mid_pattern.size() > 0){
                new_compacted_patterns.emplace_back(mid_pattern);
            }

            *//*for (int j = 0; j < mid_pattern.size(); ++j) {
                std::cout<<" "<<mid_pattern[j]<<" ";
            }
            std::cout<<std::endl;*//*
            mid_pattern.clear();
        }


        std::cout<<"compacted pattern number is : "<<new_compacted_patterns.size()<<std::endl;*/



    }

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompaction(int pattern_per_sim) {

    std::string time_file_str = "/home/eda/CLionProjects/ICTest-ATPG-system/example/SAT_BENCH/b17/time_for_each_fault.txt";
    std::ofstream outputFile;
    outputFile.open(time_file_str);



    double pattern_xbits_rate_sum = 0;
        int pattern_count = 0;

        int test_found = 0, no_test = 0, abort = 0;
        int num_curr_pattern = 0;
        int total_pattern_num = 0;
        double sim_time = 0;
        double cnf_time = 0;
        double cnf_trans_time = 0;
        double solve_time = 0;
        int num_faults = saf_list_->GetUncollapsedSAFs().size();
        auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
        vector<std::pair<int, double>> single_time;
        std::unordered_set<SAF *> success_targets;
        std::vector<SAF *> rest_faults;
        std::vector<SAF *> success_pfs;
        for (auto fault : uncollapsed_saf_list) {
            if (fault->GetSAFStatus() == SAFStatus::UC) {
                rest_faults.emplace_back(fault);
            }
        }
        int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
        std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
        CNFGenerator saf_test_generator(prim_);

        std::vector<Pattern> patterns;
        std::vector<Pattern> total_patterns;
        std::vector<std::vector<int>> atpg_patterns;
        std::vector<std::vector<int>> validation_patterns;
        std::vector<int> new_pattern(prim_->NumPIs());
        PatternParser pp;
        pp.SetupPrimNetlist(prim_);
        pp.SetPatternType(PatternType::COMB_PT);
        assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);


        vector<vector<int>> compaction_patterns;

        vector<int> check_fault(num_faults,0);
        vector<int> new_detected_fault;

        vector<int> Necessary_fault_assignment(prim_->NumPIs(),0);


        auto start = std::chrono::steady_clock::now();


        vector<double> time_for_each_fault;

        for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {

            //time for each fault
            auto start_each_fault = std::chrono::steady_clock::now();

            auto *primary_fault = rest_faults[pf_id];



            //if (primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault)) {
            if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
                && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
//    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
                continue;
            }
            //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
            std::cout << "At Fault: " << primary_fault->GetSAFName() << " and id is"<<primary_fault->GetSAFId()<<" "
                      << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

            std::cout<<"fault gate is : "<<primary_fault->GetSAFGate()->GetGId()<<" and type is : "<<int(primary_fault->GetSAFType())<<std::endl;

            auto solver = new Solver;
            solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
            solver->max_dpi_ = prim_->MaxLevel();
            solver->fault_site_ = primary_fault->GetSAFGate();
            saf_test_generator.SetSolver(solver);





            int test_res;
            std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
            auto cnf_start = std::chrono::steady_clock::now();


            /*int used_fanin_size = 0;

            //get used fanin gates
            vector<int> visited2po(prim_->GetPrimNetlist().size(),0);
            vector<int> used_fanin_mark(prim_->GetPIGates().size(),0);
            vector<Gate*> que;
            que.emplace_back(primary_fault->GetSAFGate());
            visited2po[primary_fault->GetSAFGate()->GetGId()] = 1;
            while (!que.empty()){
                auto front = que.back();
                que.pop_back();
                for (int i = 0; i < front->FanoutSize(); ++i) {
                    if (visited2po[front->FanoutGates()[i]->GetGId()] == 0){
                        que.emplace_back(front->FanoutGates()[i]);
                        visited2po[front->FanoutGates()[i]->GetGId()] = 1;
                    }
                }
            }


            vector<Gate*> find_used_pi(prim_->GetPrimNetlist().size(),0);
            vector<int> visited2pi(prim_->GetPrimNetlist().size(),0);
            for (int i = 0; i < visited2po.size(); ++i) {
                if (visited2po[i] == 1){
                    find_used_pi.emplace_back(prim_->GetPrimNetlist()[i]);
                    visited2pi[i] = 1;
                }
            }*/



            /*while (!find_used_pi.empty()){
                auto front = find_used_pi.back();
                find_used_pi.pop_back();
                for (int i = 0; i < front->FaninSize(); ++i) {
                    if (visited2pi[front->FaninGates()[i]->GetGId()] == 0){
                        find_used_pi.emplace_back(front->FaninGates()[i]);
                        visited2pi[front->FaninGates()[i]->GetGId()] = 1;
                    }
                }
            }*/













            //std::cout<<"used pi number is : "<<used_pi_count<<std::endl;
            if (saf_test_generator.CollectDChainCircuitClauseSAFL2TPI(primary_fault,
                                                                      fanin_mark)) {


                auto cnf_end = std::chrono::steady_clock::now();
                cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
//            test_res = solver->solve(d_test_engine_, primary_fault);
                test_res = solver->solve_help_checkx();
                //test_res = solver->solve();




                int xpi = 0;
                for (auto& xval: solver->x_vars_) {
                    if(xval < prim_->GetPrimNetlist().size() && prim_->GetPrimNetlist()[xval]->GetGType() == GType::G_PI) {
                        xpi++;
                    }
                }



                std::cout << "xpi num: " << xpi << std::endl;
                cnf_trans_time += solver->cnf_trans_time_;
                solve_time += solver->solve_time_;
                single_time.push_back(std::make_pair(primary_fault->GetSAFId(), solver->solve_time_));
            } else {
                test_res = 20;
            }

            // No test.
            if (test_res == 20) {
                primary_fault->SetSAFStatus(SAFStatus::RE);
                no_test++;
                //      std::cout << "No Test" << std::endl;
            }
                // Over backtrack.
            else if (test_res == 0) {
                primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
            }
                // Test found.
            else if (test_res == 10) {
                // Load test cube from ATPG circuit model.
                int determined_bits = 0;
                int x_bits = 0;

                for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                    auto pi_gate = prim_->GetPIGates()[pi_id];
                    if (!fanin_mark[pi_gate->GetGId()] ) {
                        new_pattern[pi_id] = LOGIC_x;
                        x_bits++;
                        //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                        //          ";
                        continue;
                    }
                    assert(pi_gate->GetGType() == GType::G_PI);
                    auto pi_bit =
                            (solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
                    //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
                    determined_bits++;
                    //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
//                assert(pi_bit != LOGIC_x);
                    new_pattern[pi_id] = pi_bit;
                }




                std::cout<<"new pattern is : ";
                for (int i = 0; i < new_pattern.size(); ++i) {
                    std::cout<<" "<<new_pattern[i]<<" ";
                }
                std::cout<<std::endl;

                std::cout << "determined_bits: " << determined_bits << std::endl;
                std::cout << "x_bits: " << x_bits << std::endl;
                std::cout << "fill rate: " << double(determined_bits) / prim_->GetPIGates().size() << std::endl;

                int used_fanin = 0;
                int xbits_in_used_fanin = 0;
                for (int i = 0; i < prim_->NumPIs(); ++i) {
                    if (fanin_mark[i] == 1){
                        used_fanin++;
                        if (new_pattern[i] == 2){
                            xbits_in_used_fanin++;
                        }
                    }
                }
                double used_fanin_rate = double(xbits_in_used_fanin)/double (used_fanin);
                std::cout<<"x_bits in used fanin rate is :"<<used_fanin_rate<<std::endl;
                pattern_xbits_rate_sum = pattern_xbits_rate_sum + used_fanin_rate;
                pattern_count ++;



                atpg_patterns.emplace_back(new_pattern);
                validation_patterns.emplace_back(new_pattern);
                success_pfs.emplace_back(primary_fault);

                success_targets.insert(primary_fault);
            } else {
                std::cout << "Unexpected test gen res." << std::endl;
                exit(137);
            }



            saf_test_generator.ResetSolver();
            delete solver;
#ifdef DEBUG_MODE
            DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
            num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
            if (num_curr_pattern == 1) {
#else
            if (num_curr_pattern == 1) {
#endif
                auto sim_start = std::chrono::steady_clock::now();
                total_pattern_num += num_curr_pattern;



                pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);




                sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);






                for (int i = 0; i < check_fault.size(); ++i) {
                    //std::cout<<int(uncollapsed_saf_list[i]->GetSAFStatus())<<std::endl;
                    if (int(uncollapsed_saf_list[i]->GetSAFStatus()) == 0 && check_fault[i] == 0){
                        new_detected_fault.emplace_back(uncollapsed_saf_list[i]->GetSAFId());
                        check_fault.at(i) = 1;
                    }
                }

                std::cout<<"new detected fault is :";
                for (int i = 0; i < new_detected_fault.size(); ++i) {
                    std::cout<<new_detected_fault[i]<<" ";
                }
                std::cout<<std::endl;
                new_detected_fault.clear();



                compaction_patterns.emplace_back(atpg_patterns[0]);

                atpg_patterns.clear();
                total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
                patterns.clear();
                auto sim_end = std::chrono::steady_clock::now();
                sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
                int cnt = 0;
      for (int i = 0; i < DETECT_FAULT.size(); ++i) {
          if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
              std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
              cnt++;
          }
      }
      if (cnt != 0) {
          std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
      }
      std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
      DETECT_FAULT.clear();
      for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
          if ((*it)->GetSAFStatus() != SAFStatus::DS) {
              std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
              " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
              exit(1);
          }
      }
#endif
            }
            auto end_each_fault = std::chrono::steady_clock::now();
                auto fault_time = std::chrono::duration<double>(end_each_fault-start_each_fault).count();
                //time_for_each_fault.emplace_back(fault_time);
            if (outputFile.is_open()){
                outputFile<<fault_time<<std::endl;
            }else{
                std::cout<<"fail to open txt!"<<std::endl;
            }


        }
    /*std::cout<<"time for each fault is: ";
    for (int i = 0; i < time_for_each_fault.size(); ++i) {
        std::cout<<time_for_each_fault[i]<<std::endl;
    }*/
    outputFile.close();

        num_curr_pattern = static_cast<int>(atpg_patterns.size());
        if (num_curr_pattern) {
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
        }

        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration<double>(end - start);

        int fault_detect = 0;
        int fault_redundant = 0;
        int fault_not_detect = 0;
        double fault_coverage = 0;
        double test_coverage = 0;
        DumpSAFList();
        for (auto fault : saf_list_->GetUncollapsedSAFs()) {
            if (fault->GetSAFStatus() == SAFStatus::DS ||
                fault->GetSAFStatus() == SAFStatus::DI ||
                fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
                fault_detect++;
            } else if (fault->GetSAFStatus() == SAFStatus::RE) {
                fault_redundant++;
            } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                       fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                       fault->GetSAFStatus() == SAFStatus::UO) {
                fault_not_detect++;
            }
        }
        fault_coverage =
                (double)fault_detect / (double)uncollapsed_saf_list.size();
        test_coverage =
                (double)(fault_detect + fault_redundant) /
                (double)uncollapsed_saf_list.size();

        std::cout << "Detected Faults: " << fault_detect << std::endl;
        std::cout << "No Test Faults: " << fault_redundant << std::endl;
        std::cout << "Abort Faults: " << fault_not_detect << std::endl;
        std::cout << "Pattern generated: " << total_pattern_num << std::endl;
        std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
                  << std::endl;
        std::cout << "Test Coverage: " << test_coverage * 100 << "%"
                  << std::endl;
        std::cout << "Test time: " << duration.count() << "s." << std::endl;
        std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
        std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
        std::cout << "Solver time: " << solve_time << "s." << std::endl;
        std::cout << "Sim time: " << sim_time << "s." << std::endl;


        std::cout<<"Average x_bits fill rate is used fanin : "<<pattern_xbits_rate_sum/pattern_count<<std::endl;


        // ==================== VALIDATION ===========================
        std::cout << std::endl;
        std::cout << "----------Validate test patterns-----------" << std::endl;
        for (auto fault : uncollapsed_saf_list) {
            fault->SetSAFStatus(SAFStatus::UC);
        }
        pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
        sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
        validation_patterns.clear();
        patterns.clear();

        for (int i = 0; i < uncollapsed_saf_list.size(); ++i) {
            if(int(uncollapsed_saf_list[i]->GetSAFStatus()) != 0){
                std::cout<<"undetected fault id is :"<<i<<std::endl;
                std::cout<<"type is :"<<int(uncollapsed_saf_list[i]->GetSAFStatus())<<std::endl;
            }
        }

        // ==================== VALIDATION END ==============
        // ==================== CHECK X BIT ===========================
        std::cout << std::endl;
        std::cout << "----------Check X bit in test patterns-----------" << std::endl;
        assert(success_pfs.size() == total_patterns.size());
        for (int i = 0; i < total_patterns.size(); ++i) {
            auto& pat_cycle_pis = total_patterns[i].GetPatternCycle(0).GetInputVal();
            assert(pat_cycle_pis.size() == 1);
            val64_t ori_val;
            int add_x_bit = 0;
            int ori_nox_bit = 0;
            for (int j = 0; j < pat_cycle_pis[0].size(); ++j) {
                auto& pi_v = pat_cycle_pis[0][j];
                val64_t x_v (0, 0);
                if (pi_v != x_v) {
                    ori_nox_bit++;
                    ori_val = pi_v;
                    pi_v.SetX();
                    std::vector<Pattern> one_pat;
                    one_pat.push_back(total_patterns[i]);
                    std::vector<SAF*> one_pf;
                    success_pfs[i]->SetSAFStatus(SAFStatus::UC);
                    one_pf.push_back(success_pfs[i]);
                    sim_->RunSAFSimulation(one_pat, one_pf, false);
                    if (one_pf[0]->GetSAFStatus() != SAFStatus::DS) {
                        pi_v = ori_val;
                        one_pat.clear();
                        one_pat.push_back(total_patterns[i]);
                        sim_->RunSAFSimulation(one_pat, one_pf, false);
                        assert(one_pf[0]->GetSAFStatus() == SAFStatus::DS);
                    } else {
                        add_x_bit++;
                    }
                }
            }
            std::cout << "pat " << i << ": add_x: " << add_x_bit <<"; total: " << ori_nox_bit << std::endl;
        }
        // ==================== CHECK X BIT END ==============
//    // ==================== CHECK OVERLAP ===========================
//    std::cout << std::endl;
//    std::cout << "----------Check overlap in test patterns-----------" << std::endl;
//    int overlap_bit = 0;
//    vector<int> overlap_bits;
//    for (int i = 0; i < total_patterns.size() - 1; ++i) {
//        for (int j = i + 1; j < total_patterns.size(); ++j) {
//            auto& pat_cycle_pis1 = total_patterns[i].GetPatternCycle(0).GetInputVal();
//            auto& pat_cycle_pis2 = total_patterns[j].GetPatternCycle(0).GetInputVal();
//            assert(pat_cycle_pis1.size() == pat_cycle_pis2.size());
//            for (int k = 0; k < pat_cycle_pis1[0].size(); ++k) {
//                auto& pi_v1 = pat_cycle_pis1[0][k];
//                auto& pi_v2 = pat_cycle_pis2[0][k];
//                if (pi_v1 != pi_v2) {
//                    overlap_bit++;
//                }
//            }
//            overlap_bits.push_back(overlap_bit);
//            overlap_bit = 0;
//        }
//    }
//    double avg = std::accumulate(overlap_bits.begin(), overlap_bits.end(), 0.0);
//    avg = avg / overlap_bits.size();
//    std::cout << "avg overlap: " << avg << std::endl;
//    std::sort(overlap_bits.begin(), overlap_bits.end());
//    double min_ = overlap_bits[0];
//    double max_ = overlap_bits.back();
//    std::cout << "max overlap: " << max_ << std::endl;
//    std::cout << "min overlap: " << min_ << std::endl;
//    std::cout << "mean overlap: " << overlap_bits[overlap_bits.size()/2] << std::endl;
//    // ==================== CHECK OVERLAP END ==============
        std::sort(single_time.begin(), single_time.end(),
                  [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                      return a.second > b.second;
                  });

        // 打印前十个最大的项
        for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
            std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
        }


    std::cout<<"compaction_pattern_size is: "<<compaction_patterns.size()<<std::endl;
    //start compaction step_1
    vector<vector<int>> new_compacted_patterns;
    new_compacted_patterns.clear();
    vector<int> mid_pattern;
    vector<int> fanin_mark(prim_->NumPIs(),0);
    vector<int> pattern_compcted_mark(compaction_patterns.size(),0);
    for (int i = 0; i < compaction_patterns.size(); ++i) {
        if (pattern_compcted_mark[i] == 0){
            pattern_compcted_mark[i] = 1;
            for (int j = 0; j < compaction_patterns[i].size(); ++j) {
                mid_pattern.emplace_back(compaction_patterns[i][j]);
            }
            for (int j = i + 1; j < compaction_patterns.size(); ++j) {
                int conflict_flag = 0;
                for (int k = 0; k < compaction_patterns[j].size(); ++k) {
                    if (mid_pattern[k] != compaction_patterns[j][k] && mid_pattern[k] != 2 && compaction_patterns[j][k] != 2){
                        conflict_flag = 1;
                    }
                }
                if (conflict_flag == 0){
                    pattern_compcted_mark[j] = 1;
                    for (int k = 0; k < mid_pattern.size(); ++k) {
                        if (mid_pattern[k] == 2){
                            mid_pattern[k] = compaction_patterns[j][k];
                        }
                    }
                }

            }
        }
        if (mid_pattern.size() > 0){
            new_compacted_patterns.emplace_back(mid_pattern);
        }

        /*for (int j = 0; j < mid_pattern.size(); ++j) {
            std::cout<<" "<<mid_pattern[j]<<" ";
        }
        std::cout<<std::endl;*/
        mid_pattern.clear();
    }


    std::cout<<"compacted pattern number is : "<<new_compacted_patterns.size()<<std::endl;



    }

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompactionDynamic_2(int pattern_per_sim) {
    double pattern_xbits_rate_sum = 0;
    int pattern_count = 0;

    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    vector<std::pair<int, double>> single_time;
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    std::vector<SAF *> success_pfs;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;
    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);


    vector<vector<int>> compaction_patterns;

    vector<int> check_fault(num_faults,0);
    vector<int> new_detected_fault;

    vector<int> Necessary_fault_assignment(prim_->NumPIs(),0);


    auto start = std::chrono::steady_clock::now();


    vector<double> time_for_each_fault;

    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {

        //time for each fault
        auto start_each_fault = std::chrono::steady_clock::now();

        auto *primary_fault = rest_faults[pf_id];



        //if (primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault)) {
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
            //    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " and id is"<<primary_fault->GetSAFId()<<" "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        //std::cout<<"fault gate is : "<<primary_fault->GetSAFGate()->GetGId()<<" and type is : "<<int(primary_fault->GetSAFType())<<std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        solver->max_dpi_ = prim_->MaxLevel();
        solver->fault_site_ = primary_fault->GetSAFGate();
        saf_test_generator.SetSolver(solver);





        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();


        /*int used_fanin_size = 0;

        //get used fanin gates
        vector<int> visited2po(prim_->GetPrimNetlist().size(),0);
        vector<int> used_fanin_mark(prim_->GetPIGates().size(),0);
        vector<Gate*> que;
        que.emplace_back(primary_fault->GetSAFGate());
        visited2po[primary_fault->GetSAFGate()->GetGId()] = 1;
        while (!que.empty()){
            auto front = que.back();
            que.pop_back();
            for (int i = 0; i < front->FanoutSize(); ++i) {
                if (visited2po[front->FanoutGates()[i]->GetGId()] == 0){
                    que.emplace_back(front->FanoutGates()[i]);
                    visited2po[front->FanoutGates()[i]->GetGId()] = 1;
                }
            }
        }


        vector<Gate*> find_used_pi(prim_->GetPrimNetlist().size(),0);
        vector<int> visited2pi(prim_->GetPrimNetlist().size(),0);
        for (int i = 0; i < visited2po.size(); ++i) {
            if (visited2po[i] == 1){
                find_used_pi.emplace_back(prim_->GetPrimNetlist()[i]);
                visited2pi[i] = 1;
            }
        }*/



        /*while (!find_used_pi.empty()){
            auto front = find_used_pi.back();
            find_used_pi.pop_back();
            for (int i = 0; i < front->FaninSize(); ++i) {
                if (visited2pi[front->FaninGates()[i]->GetGId()] == 0){
                    find_used_pi.emplace_back(front->FaninGates()[i]);
                    visited2pi[front->FaninGates()[i]->GetGId()] = 1;
                }
            }
        }*/






        //std::cout<<"used pi number is : "<<used_pi_count<<std::endl;
        if (saf_test_generator.CollectDChainCircuitClauseSAFL2TPI(primary_fault,
                                                                  fanin_mark)) {


            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
            //            test_res = solver->solve(d_test_engine_, primary_fault);
            test_res = solver->solve_help_checkx();
            //test_res = solver->solve();




            int xpi = 0;
            for (auto& xval: solver->x_vars_) {
                if(xval < prim_->GetPrimNetlist().size() && prim_->GetPrimNetlist()[xval]->GetGType() == GType::G_PI) {
                    xpi++;
                }
            }



            //std::cout << "xpi num: " << xpi << std::endl;
            cnf_trans_time += solver->cnf_trans_time_;
            solve_time += solver->solve_time_;
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), solver->solve_time_));
        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
        // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
        // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            int determined_bits = 0;
            int x_bits = 0;

            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                auto pi_gate = prim_->GetPIGates()[pi_id];
                if (!fanin_mark[pi_gate->GetGId()] ) {
                    new_pattern[pi_id] = LOGIC_x;
                    x_bits++;
                    //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                    //          ";
                    continue;
                }
                assert(pi_gate->GetGType() == GType::G_PI);
                auto pi_bit =
                    (solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
                //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
                determined_bits++;
                //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
                //                assert(pi_bit != LOGIC_x);
                new_pattern[pi_id] = pi_bit;
            }



            /*std::cout << "determined_bits: " << determined_bits << std::endl;
            std::cout << "x_bits: " << x_bits << std::endl;
            std::cout << "fill rate: " << double(determined_bits) / prim_->GetPIGates().size() << std::endl;*/

            /*int used_fanin = 0;
            int xbits_in_used_fanin = 0;
            for (int i = 0; i < prim_->NumPIs(); ++i) {
                if (fanin_mark[i] == 1){
                    used_fanin++;
                    if (new_pattern[i] == 2){
                        xbits_in_used_fanin++;
                    }
                }
            }
            double used_fanin_rate = double(xbits_in_used_fanin)/double (used_fanin);
            std::cout<<"x_bits in used fanin rate is :"<<used_fanin_rate<<std::endl;
            pattern_xbits_rate_sum = pattern_xbits_rate_sum + used_fanin_rate;
            pattern_count ++;*/

            int compact_fail_count = 0;
            int continue_flag = 0;
            for (int i = pf_id + 1; i < num_rest_faults && continue_flag != 1; ++i) {

                continue_flag = 0;

                auto *primary_fault_compaction = rest_faults[i];
                if ((primary_fault_compaction->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault_compaction))
                    && primary_fault_compaction->GetSAFStatus() != SAFStatus::PT && primary_fault_compaction->GetSAFStatus() != SAFStatus::PU) {
                    //    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
                    continue;
                }
                auto solver_2 = new Solver;
                solver_2->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
                solver_2->max_dpi_ = prim_->MaxLevel();
                solver_2->fault_site_ = primary_fault_compaction->GetSAFGate();
                saf_test_generator.SetSolver(solver_2);

                int test_res_2;
                std::vector<int> fanin_mark_2(prim_->GetPrimNetlist().size(), 0);

                if (saf_test_generator.CollectDChainCircuitClauseSAFL2TPI(primary_fault_compaction,
                                                                          fanin_mark_2)){
                    for (int j = 0; j < new_pattern.size(); ++j) {
                        if (new_pattern[j] == 1){
                            solver_2->add(prim_->GetPIGates()[j]->GetGId() + 1);
                            solver_2->add(0);
                        }
                        else if(new_pattern[j] == 0){
                              solver_2->add(-1 * (prim_->GetPIGates()[j]->GetGId() + 1));
                              solver_2->add(0);
                        }

                    }
                    test_res_2 = solver_2->solve_help_checkx();

                    if (test_res_2 == 10){
                        for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                              auto pi_gate = prim_->GetPIGates()[pi_id];
                              assert(pi_gate->GetGType() == GType::G_PI);
                              auto pi_bit =
                                  (solver_2->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver_2->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
                              //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
                              determined_bits++;
                              //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
                              //                assert(pi_bit != LOGIC_x);
                              if (new_pattern[pi_id] == LOGIC_x){
                                new_pattern[pi_id] = pi_bit;
                              }

                        }

                    compact_fail_count = 0;
                    }else{
                        compact_fail_count ++;
                    }
                }

                if (compact_fail_count > 50){
                    continue_flag = 1;
                }

            delete solver_2;
            }

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);
            success_pfs.emplace_back(primary_fault);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }



        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
        if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;

            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);

            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);


            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
            for (int i = 0; i < DETECT_FAULT.size(); ++i) {
                if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
                    std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
                    cnt++;
                }
            }
            if (cnt != 0) {
                std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
            }
            std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
            DETECT_FAULT.clear();
            for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
                if ((*it)->GetSAFStatus() != SAFStatus::DS) {
                    std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
                        " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
                    exit(1);
                }
            }
#endif
        }
        auto end_each_fault = std::chrono::steady_clock::now();
        auto fault_time = std::chrono::duration<double>(end_each_fault-start_each_fault).count();
        //time_for_each_fault.emplace_back(fault_time);
        /*if (outputFile.is_open()){
            outputFile<<fault_time<<std::endl;
        }else{
            std::cout<<"fail to open txt!"<<std::endl;
        }*/


    }
    /*std::cout<<"time for each fault is: ";
    for (int i = 0; i < time_for_each_fault.size(); ++i) {
        std::cout<<time_for_each_fault[i]<<std::endl;
    }*/
    //outputFile.close();

    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
        (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
        (double)(fault_detect + fault_redundant) /
        (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;


    //std::cout<<"Average x_bits fill rate is used fanin : "<<pattern_xbits_rate_sum/pattern_count<<std::endl;


    std::ofstream outputfile;
    outputfile.open("/home/dell/Desktop/ICTest-ATPG-system/compaction output/output.txt");

    if (outputfile.is_open()){
        outputfile<<"Pattern generated: "<< total_pattern_num << std::endl;
        outputfile<< "Test Coverage: " << test_coverage * 100 << "%"<< std::endl;
    }
    outputfile.close();




    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }
    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();

    /*for (int i = 0; i < uncollapsed_saf_list.size(); ++i) {
        if(int(uncollapsed_saf_list[i]->GetSAFStatus()) != 0){
            std::cout<<"undetected fault id is :"<<i<<std::endl;
            std::cout<<"type is :"<<int(uncollapsed_saf_list[i]->GetSAFStatus())<<std::endl;
        }
    }*/

    // ==================== VALIDATION END ==============
    // ==================== CHECK X BIT ===========================
    std::cout << std::endl;
    std::cout << "----------Check X bit in test patterns-----------" << std::endl;
    assert(success_pfs.size() == total_patterns.size());
    for (int i = 0; i < total_patterns.size(); ++i) {
        auto& pat_cycle_pis = total_patterns[i].GetPatternCycle(0).GetInputVal();
        assert(pat_cycle_pis.size() == 1);
        val64_t ori_val;
        int add_x_bit = 0;
        int ori_nox_bit = 0;
        for (int j = 0; j < pat_cycle_pis[0].size(); ++j) {
            auto& pi_v = pat_cycle_pis[0][j];
            val64_t x_v (0, 0);
            if (pi_v != x_v) {
                ori_nox_bit++;
                ori_val = pi_v;
                pi_v.SetX();
                std::vector<Pattern> one_pat;
                one_pat.push_back(total_patterns[i]);
                std::vector<SAF*> one_pf;
                success_pfs[i]->SetSAFStatus(SAFStatus::UC);
                one_pf.push_back(success_pfs[i]);
                sim_->RunSAFSimulation(one_pat, one_pf, false);
                if (one_pf[0]->GetSAFStatus() != SAFStatus::DS) {
                    pi_v = ori_val;
                    one_pat.clear();
                    one_pat.push_back(total_patterns[i]);
                    sim_->RunSAFSimulation(one_pat, one_pf, false);
                    assert(one_pf[0]->GetSAFStatus() == SAFStatus::DS);
                } else {
                    add_x_bit++;
                }
            }
        }
        std::cout << "pat " << i << ": add_x: " << add_x_bit <<"; total: " << ori_nox_bit << std::endl;
    }
    // ==================== CHECK X BIT END ==============
    //    // ==================== CHECK OVERLAP ===========================
    //    std::cout << std::endl;
    //    std::cout << "----------Check overlap in test patterns-----------" << std::endl;
    //    int overlap_bit = 0;
    //    vector<int> overlap_bits;
    //    for (int i = 0; i < total_patterns.size() - 1; ++i) {
    //        for (int j = i + 1; j < total_patterns.size(); ++j) {
    //            auto& pat_cycle_pis1 = total_patterns[i].GetPatternCycle(0).GetInputVal();
    //            auto& pat_cycle_pis2 = total_patterns[j].GetPatternCycle(0).GetInputVal();
    //            assert(pat_cycle_pis1.size() == pat_cycle_pis2.size());
    //            for (int k = 0; k < pat_cycle_pis1[0].size(); ++k) {
    //                auto& pi_v1 = pat_cycle_pis1[0][k];
    //                auto& pi_v2 = pat_cycle_pis2[0][k];
    //                if (pi_v1 != pi_v2) {
    //                    overlap_bit++;
    //                }
    //            }
    //            overlap_bits.push_back(overlap_bit);
    //            overlap_bit = 0;
    //        }
    //    }
    //    double avg = std::accumulate(overlap_bits.begin(), overlap_bits.end(), 0.0);
    //    avg = avg / overlap_bits.size();
    //    std::cout << "avg overlap: " << avg << std::endl;
    //    std::sort(overlap_bits.begin(), overlap_bits.end());
    //    double min_ = overlap_bits[0];
    //    double max_ = overlap_bits.back();
    //    std::cout << "max overlap: " << max_ << std::endl;
    //    std::cout << "min overlap: " << min_ << std::endl;
    //    std::cout << "mean overlap: " << overlap_bits[overlap_bits.size()/2] << std::endl;
    //    // ==================== CHECK OVERLAP END ==============
    std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second > b.second;
              });

    // 打印前十个最大的项
    /*for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }*/


    /*std::cout<<"compaction_pattern_size is: "<<compaction_patterns.size()<<std::endl;
    //start compaction step_1
    vector<vector<int>> new_compacted_patterns;
    new_compacted_patterns.clear();
    vector<int> mid_pattern;
    vector<int> fanin_mark(prim_->NumPIs(),0);
    vector<int> pattern_compcted_mark(compaction_patterns.size(),0);
    for (int i = 0; i < compaction_patterns.size(); ++i) {
        if (pattern_compcted_mark[i] == 0){
            pattern_compcted_mark[i] = 1;
            for (int j = 0; j < compaction_patterns[i].size(); ++j) {
                mid_pattern.emplace_back(compaction_patterns[i][j]);
            }
            for (int j = i + 1; j < compaction_patterns.size(); ++j) {
                int conflict_flag = 0;
                for (int k = 0; k < compaction_patterns[j].size(); ++k) {
                    if (mid_pattern[k] != compaction_patterns[j][k] && mid_pattern[k] != 2 && compaction_patterns[j][k] != 2){
                        conflict_flag = 1;
                    }
                }
                if (conflict_flag == 0){
                    pattern_compcted_mark[j] = 1;
                    for (int k = 0; k < mid_pattern.size(); ++k) {
                        if (mid_pattern[k] == 2){
                            mid_pattern[k] = compaction_patterns[j][k];
                        }
                    }
                }

            }
        }
        if (mid_pattern.size() > 0){
            new_compacted_patterns.emplace_back(mid_pattern);
        }

        for (int j = 0; j < mid_pattern.size(); ++j) {
            std::cout<<" "<<mid_pattern[j]<<" ";
        }
        std::cout<<std::endl;
        mid_pattern.clear();
    }


    std::cout<<"compacted pattern number is : "<<new_compacted_patterns.size()<<std::endl;*/

}

int SATEngine::SATGeneratePat4Pf(SAF* primary_fault, std::vector<int> &new_pattern, double &cnf_gen_time, double &cnf_trans_time, double &solver_time) {
    int res_atpg_statu = -1;
    CNFGenerator saf_test_generator(prim_);
    auto solver = new Solver;
    solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
    solver->max_dpi_ = prim_->MaxLevel();
    solver->fault_site_ = primary_fault->GetSAFGate();
    saf_test_generator.SetSolver(solver);

    int test_res;
    std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
    auto cnf_start = std::chrono::steady_clock::now();

    std::cout << "name" << primary_fault->GetSAFName() << " type" << static_cast<bool>(primary_fault->GetSAFType())
              << std::endl;

    //std::cout<<"used pi number is : "<<used_pi_count<<std::endl;
    if (saf_test_generator.CollectDChainCircuitClauseSAFL2TPI(primary_fault,
                                                              fanin_mark)) {

        auto cnf_end = std::chrono::steady_clock::now();
        cnf_gen_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
        //            test_res = solver->solve(d_test_engine_, primary_fault);
//        test_res = solver->solve_help_checkx();
        test_res = solver->solve();

        int xpi = 0;
        for (auto& xval: solver->x_vars_) {
            if(xval < prim_->GetPrimNetlist().size() && prim_->GetPrimNetlist()[xval]->GetGType() == GType::G_PI) {
                xpi++;
            }
        }
        cnf_trans_time += solver->cnf_trans_time_;
        solver_time += solver->solve_time_;
    } else {
        test_res = 20;
    }

    // No test.
    if (test_res == 20) {
        primary_fault->SetSAFStatus(SAFStatus::RE);
        res_atpg_statu = NO_TEST;
        //      std::cout << "No Test" << std::endl;
    }
        // Over backtrack.
    else if (test_res == 0) {
        primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        res_atpg_statu = OVER_BACKTRACK;
    }
        // Test found.
    else if (test_res == 10) {
        // Load test cube from ATPG circuit model.
        int determined_bits = 0;
        int x_bits = 0;
        res_atpg_statu = TEST_FOUND;

        for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
            auto pi_gate = prim_->GetPIGates()[pi_id];
            if (!fanin_mark[pi_gate->GetGId()] ) {
                new_pattern[pi_id] = LOGIC_x;
                x_bits++;
                //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                //          ";
                continue;
            }
            assert(pi_gate->GetGType() == GType::G_PI);
            auto pi_bit =
                    (solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
            //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
            determined_bits++;
            //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
            //                assert(pi_bit != LOGIC_x);
            new_pattern[pi_id] = pi_bit;
        }
    }
    saf_test_generator.ResetSolver();
    delete solver;
    return res_atpg_statu;
}

int SATEngine::SATGeneratePat4Sf(SAF* secondary_fault, std::vector<int> &new_pattern, double &cnf_gen_time, double &cnf_trans_time, double &solver_time){
    int res_atpg_statu = -1;
    CNFGenerator saf_test_generator(prim_);
    auto solver = new Solver;
    solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
    solver->max_dpi_ = prim_->MaxLevel();
    solver->fault_site_ = secondary_fault->GetSAFGate();
    saf_test_generator.SetSolver(solver);

    int test_res;
    std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
    auto cnf_start = std::chrono::steady_clock::now();

    //std::cout<<"used pi number is : "<<used_pi_count<<std::endl;
    if (saf_test_generator.CollectDChainCircuitClauseSAFL2TPI(secondary_fault,
                                                              fanin_mark)) {
        for (int j = 0; j < new_pattern.size(); ++j) {
            if (new_pattern[j] == 1){
                solver->add(prim_->GetPIGates()[j]->GetGId() + 1);
                solver->add(0);
            }
            else if(new_pattern[j] == 0){
                solver->add(-1 * (prim_->GetPIGates()[j]->GetGId() + 1));
                solver->add(0);
            }

        }

        auto cnf_end = std::chrono::steady_clock::now();
        cnf_gen_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
        //            test_res = solver->solve(d_test_engine_, primary_fault);
//        test_res = solver->solve_help_checkx();
        test_res = solver->solve();

        int xpi = 0;
        for (auto& xval: solver->x_vars_) {
            if(xval < prim_->GetPrimNetlist().size() && prim_->GetPrimNetlist()[xval]->GetGType() == GType::G_PI) {
                xpi++;
            }
        }
        cnf_trans_time += solver->cnf_trans_time_;
        solver_time += solver->solve_time_;
    } else {
        test_res = 20;
    }

    // No test.
    if (test_res == 20) {
//        secondary_fault->SetSAFStatus(SAFStatus::RE);
        res_atpg_statu = NO_TEST;
        //      std::cout << "No Test" << std::endl;
    }
        // Over backtrack.
    else if (test_res == 0) {
//        secondary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        res_atpg_statu = OVER_BACKTRACK;
    }
        // Test found.
    else if (test_res == 10) {
        // Load test cube from ATPG circuit model.
        int determined_bits = 0;
        res_atpg_statu = TEST_FOUND;
        int x_bits = 0;

        for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
            auto pi_gate = prim_->GetPIGates()[pi_id];
            assert(pi_gate->GetGType() == GType::G_PI);
            auto pi_bit =
                    (solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
            //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
            determined_bits++;
            //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
            //                assert(pi_bit != LOGIC_x);
            if (new_pattern[pi_id] == LOGIC_x){
                new_pattern[pi_id] = pi_bit;
            }
        }
    }
    saf_test_generator.ResetSolver();
    delete solver;
    return res_atpg_statu;

}

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompaction_czt(int pattern_per_sim) {
    LOG_ASSERT(pattern_per_sim >= 1 && pattern_per_sim <= 64,
               "Illegal Simulated Patterns");

    int atpg_status = NO_TEST;
    int num_curr_pattern = 0;
    int num_total_pattern = 0;
    int max_fail_secondary = 50;
    std::unordered_set<SAF *> success_targets;
    auto &uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    std::vector<SAF *> rest_faults;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

    std::vector<std::vector<int>> atpg_pi_pats;
    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    double cnf_gen_time = 0, cnf_trans_time = 0, solver_time = 0;


    PatternParser pp;
    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    CNFGenerator saf_test_generator(prim_);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);

    LOG_INFO("{} ...", "ATPG Start")
    LOG_INFO("{:^15}{:^15}{:^15}", "total_pat_nums", "detected", "coverage")

    auto start = std::chrono::steady_clock::now();

    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {
//    int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
//    std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
//        if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        std::vector<int> new_pattern_sat(prim_->NumPIs());

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        solver->max_dpi_ = prim_->MaxLevel();
        solver->fault_site_ = primary_fault->GetSAFGate();
        saf_test_generator.SetSolver(solver);

        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();

        atpg_status = SATGeneratePat4Pf(primary_fault, new_pattern_sat, cnf_gen_time, cnf_trans_time, solver_time);
        if (atpg_status == OVER_BACKTRACK) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
            continue;
        }
        if (atpg_status == NO_TEST) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            continue;
        }

        assert(atpg_status == TEST_FOUND);

        success_targets.insert(primary_fault);
        int try_secondary_faults = 0;
        int fail_secondary_faults = 0;
        for (int sf_id = pf_id + 1; sf_id < num_rest_faults; sf_id++) {
            if (fail_secondary_faults > max_fail_secondary) {
                break;
            }
//        int srand_idx = rand() % (num_rest_faults - sf_id) + sf_id;
//        std::swap(rest_faults[sf_id], rest_faults[srand_idx]);
            auto *secondary_fault = rest_faults[sf_id];
            if ((secondary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(secondary_fault))
                && secondary_fault->GetSAFStatus() != SAFStatus::PT && secondary_fault->GetSAFStatus() != SAFStatus::PU) {
//            if (secondary_fault->GetSAFStatus() != SAFStatus::UC) {
                continue;
            }

            try_secondary_faults++;
            atpg_status = SATGeneratePat4Sf(secondary_fault, new_pattern_sat, cnf_gen_time, cnf_trans_time, solver_time);

            if (atpg_status == OVER_BACKTRACK || atpg_status == NO_TEST) {
                fail_secondary_faults++;
                continue;
            }
            assert(atpg_status == TEST_FOUND);
            success_targets.insert(secondary_fault);
        }


        atpg_pi_pats.emplace_back(new_pattern_sat);

        num_curr_pattern = static_cast<int>(atpg_pi_pats.size());
        if (num_curr_pattern == pattern_per_sim) {
            num_total_pattern += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_pi_pats, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list, true);
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            atpg_pi_pats.clear();
            patterns.clear();
        }
    }
    num_curr_pattern = static_cast<int>(atpg_pi_pats.size());
    if (num_curr_pattern) {
        num_total_pattern += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_pi_pats, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list, true);
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        atpg_pi_pats.clear();
        patterns.clear();
    }
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);
    DumpSAFList();

    LOG_INFO("#Pattern : {}", num_total_pattern)
    LOG_INFO("RunTime  : {}", fmt::format("{:.3f} s", duration.count()))
    LOG_INFO("{} ...", "ATPG End")
    LOG_INFO("")
    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;

    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
            (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
            (double)(fault_detect + fault_redundant) /
            (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << num_total_pattern << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_gen_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solver_time << "s." << std::endl;
    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::DS)
            fault->SetSAFStatus(SAFStatus::UC);
    }
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);

    for (int i = 0; i < uncollapsed_saf_list.size(); ++i) {
        if(uncollapsed_saf_list[i]->GetSAFStatus() != SAFStatus::DS && uncollapsed_saf_list[i]->GetSAFStatus() != SAFStatus::RE){
            std::cout << "undetected fault id is :" << i << " name" << uncollapsed_saf_list[i]->GetSAFName() << " type" << static_cast<bool>(uncollapsed_saf_list[i]->GetSAFType())
                      << std::endl;
            std::cout<<"status is :"<<int(uncollapsed_saf_list[i]->GetSAFStatus())<<std::endl;
        }
    }

    // ==================== VALIDATION END ==============
}

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaults(
        int pattern_per_sim) {

    double pattern_xbits_rate_sum = 0;
    int pattern_count = 0;

    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    vector<std::pair<int, double>> single_time;
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    std::vector<SAF *> success_pfs;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;
    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);



    vector<int> check_fault(num_faults,0);
    vector<int> new_detected_fault;



    //std::cout<<"fault 27 :gate is :"<<uncollapsed_saf_list[27]->GetSAFGate()->GetGId()<<",type is :"<<int(uncollapsed_saf_list[27]->GetSAFType())<<std::endl;
    //std::cout<<"fault 62 :gate is :"<<uncollapsed_saf_list[62]->GetSAFGate()->GetGId()<<",type is :"<<int(uncollapsed_saf_list[62]->GetSAFType())<<std::endl;


    /*for (int i = 0; i < prim_->NumGates(); ++i) {
        std::cout<<"gate id is: "<<prim_->GetPrimNetlist()[i]->GetGId()<<" type is : "<<int(prim_->GetPrimNetlist()[i]->GetGType())<<std::endl;
        for (int j = 0; j < prim_->GetPrimNetlist()[i]->FaninSize(); ++j) {
            std::cout<<"fan in is : "<<prim_->GetPrimNetlist()[i]->FaninGates()[j]->GetGId()<<std::endl;
        }
    }*/


    auto start = std::chrono::steady_clock::now();
    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {

        //std::cout<<"fault 62 ststus is :"<<int(uncollapsed_saf_list[62]->GetSAFStatus())<<std::endl;

        //int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
        //std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];



        //if (primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault)) {
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
//    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " and id is"<<primary_fault->GetSAFId()<<" "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        std::cout<<"fault gate is : "<<primary_fault->GetSAFGate()->GetGId()<<" and type is : "<<int(primary_fault->GetSAFType())<<std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        solver->max_dpi_ = prim_->MaxLevel();
        solver->fault_site_ = primary_fault->GetSAFGate();
        saf_test_generator.SetSolver(solver);





        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();


        /*int used_fanin_size = 0;

        //get used fanin gates
        vector<int> visited2po(prim_->GetPrimNetlist().size(),0);
        vector<int> used_fanin_mark(prim_->GetPIGates().size(),0);
        vector<Gate*> que;
        que.emplace_back(primary_fault->GetSAFGate());
        visited2po[primary_fault->GetSAFGate()->GetGId()] = 1;
        while (!que.empty()){
            auto front = que.back();
            que.pop_back();
            for (int i = 0; i < front->FanoutSize(); ++i) {
                if (visited2po[front->FanoutGates()[i]->GetGId()] == 0){
                    que.emplace_back(front->FanoutGates()[i]);
                    visited2po[front->FanoutGates()[i]->GetGId()] = 1;
                }
            }
        }


        vector<Gate*> find_used_pi(prim_->GetPrimNetlist().size(),0);
        vector<int> visited2pi(prim_->GetPrimNetlist().size(),0);
        for (int i = 0; i < visited2po.size(); ++i) {
            if (visited2po[i] == 1){
                find_used_pi.emplace_back(prim_->GetPrimNetlist()[i]);
                visited2pi[i] = 1;
            }
        }*/



        /*while (!find_used_pi.empty()){
            auto front = find_used_pi.back();
            find_used_pi.pop_back();
            for (int i = 0; i < front->FaninSize(); ++i) {
                if (visited2pi[front->FaninGates()[i]->GetGId()] == 0){
                    find_used_pi.emplace_back(front->FaninGates()[i]);
                    visited2pi[front->FaninGates()[i]->GetGId()] = 1;
                }
            }
        }*/













        //std::cout<<"used pi number is : "<<used_pi_count<<std::endl;
        if (saf_test_generator.CollectDChainCircuitClauseSAFL2TPI(primary_fault,
                                                               fanin_mark)) {


            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
//            test_res = solver->solve(d_test_engine_, primary_fault);
            test_res = solver->solve_help_checkx();
            //test_res = solver->solve();




            int xpi = 0;
            for (auto& xval: solver->x_vars_) {
                if(xval < prim_->GetPrimNetlist().size() && prim_->GetPrimNetlist()[xval]->GetGType() == GType::G_PI) {
                    xpi++;
                }
            }



            std::cout << "xpi num: " << xpi << std::endl;
            cnf_trans_time += solver->cnf_trans_time_;
            solve_time += solver->solve_time_;
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), solver->solve_time_));
        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
            // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
            // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            int determined_bits = 0;
            int x_bits = 0;

            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                auto pi_gate = prim_->GetPIGates()[pi_id];
                if (!fanin_mark[pi_gate->GetGId()] ) {
                    new_pattern[pi_id] = LOGIC_x;
                    x_bits++;
                    //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                    //          ";
                    continue;
                }
                assert(pi_gate->GetGType() == GType::G_PI);
                auto pi_bit =
                        (solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : ((solver->val(pi_gate->GetGId() + 1) == 0) ? LOGIC_0 : LOGIC_x);
                        //(solver->val(pi_gate->GetGId() + 1) > 0) ? LOGIC_1 : LOGIC_0;
                determined_bits++;
                //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
//                assert(pi_bit != LOGIC_x);
                new_pattern[pi_id] = pi_bit;
            }




            std::cout<<"new pattern is : ";
            for (int i = 0; i < new_pattern.size(); ++i) {
                std::cout<<" "<<new_pattern[i]<<" ";
            }
            std::cout<<std::endl;

            std::cout << "determined_bits: " << determined_bits << std::endl;
            std::cout << "x_bits: " << x_bits << std::endl;
            std::cout << "fill rate: " << double(determined_bits) / prim_->GetPIGates().size() << std::endl;

            int used_fanin = 0;
            int xbits_in_used_fanin = 0;
            for (int i = 0; i < prim_->NumPIs(); ++i) {
                if (fanin_mark[i] == 1){
                    used_fanin++;
                    if (new_pattern[i] == 2){
                        xbits_in_used_fanin++;
                    }
                }
            }
            double used_fanin_rate = double(xbits_in_used_fanin)/double (used_fanin);
            std::cout<<"x_bits in used fanin rate is :"<<used_fanin_rate<<std::endl;
            pattern_xbits_rate_sum = pattern_xbits_rate_sum + used_fanin_rate;
            pattern_count ++;



            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);
            success_pfs.emplace_back(primary_fault);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }



        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
        if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);




            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);


            for (int i = 0; i < check_fault.size(); ++i) {
                //std::cout<<int(uncollapsed_saf_list[i]->GetSAFStatus())<<std::endl;
                if (int(uncollapsed_saf_list[i]->GetSAFStatus()) == 0 && check_fault[i] == 0){
                    new_detected_fault.emplace_back(uncollapsed_saf_list[i]->GetSAFId());
                    check_fault.at(i) = 1;
                }
            }

            std::cout<<"new detected fault is :";
            for (int i = 0; i < new_detected_fault.size(); ++i) {
                std::cout<<new_detected_fault[i]<<" ";
            }
            std::cout<<std::endl;
            new_detected_fault.clear();




            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
      for (int i = 0; i < DETECT_FAULT.size(); ++i) {
          if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
              std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
              cnt++;
          }
      }
      if (cnt != 0) {
          std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
      }
      std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
      DETECT_FAULT.clear();
      for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
          if ((*it)->GetSAFStatus() != SAFStatus::DS) {
              std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
              " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
              exit(1);
          }
      }
#endif
        }
    }
    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
            (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
            (double)(fault_detect + fault_redundant) /
            (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;


    std::cout<<"Average x_bits fill rate is used fanin : "<<pattern_xbits_rate_sum/pattern_count<<std::endl;


    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }
    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();

    for (int i = 0; i < uncollapsed_saf_list.size(); ++i) {
        if(int(uncollapsed_saf_list[i]->GetSAFStatus()) != 0){
            std::cout<<"undetected fault id is :"<<i<<std::endl;
            std::cout<<"type is :"<<int(uncollapsed_saf_list[i]->GetSAFStatus())<<std::endl;
        }
    }

    // ==================== VALIDATION END ==============
    // ==================== CHECK X BIT ===========================
    std::cout << std::endl;
    std::cout << "----------Check X bit in test patterns-----------" << std::endl;
    assert(success_pfs.size() == total_patterns.size());
    for (int i = 0; i < total_patterns.size(); ++i) {
        auto& pat_cycle_pis = total_patterns[i].GetPatternCycle(0).GetInputVal();
        assert(pat_cycle_pis.size() == 1);
        val64_t ori_val;
        int add_x_bit = 0;
        int ori_nox_bit = 0;
        for (int j = 0; j < pat_cycle_pis[0].size(); ++j) {
            auto& pi_v = pat_cycle_pis[0][j];
            val64_t x_v (0, 0);
            if (pi_v != x_v) {
                ori_nox_bit++;
                ori_val = pi_v;
                pi_v.SetX();
                std::vector<Pattern> one_pat;
                one_pat.push_back(total_patterns[i]);
                std::vector<SAF*> one_pf;
                success_pfs[i]->SetSAFStatus(SAFStatus::UC);
                one_pf.push_back(success_pfs[i]);
                sim_->RunSAFSimulation(one_pat, one_pf, false);
                if (one_pf[0]->GetSAFStatus() != SAFStatus::DS) {
                    pi_v = ori_val;
                    one_pat.clear();
                    one_pat.push_back(total_patterns[i]);
                    sim_->RunSAFSimulation(one_pat, one_pf, false);
                    assert(one_pf[0]->GetSAFStatus() == SAFStatus::DS);
                } else {
                    add_x_bit++;
                }
            }
        }
        std::cout << "pat " << i << ": add_x: " << add_x_bit <<"; total: " << ori_nox_bit << std::endl;
    }
    // ==================== CHECK X BIT END ==============
//    // ==================== CHECK OVERLAP ===========================
//    std::cout << std::endl;
//    std::cout << "----------Check overlap in test patterns-----------" << std::endl;
//    int overlap_bit = 0;
//    vector<int> overlap_bits;
//    for (int i = 0; i < total_patterns.size() - 1; ++i) {
//        for (int j = i + 1; j < total_patterns.size(); ++j) {
//            auto& pat_cycle_pis1 = total_patterns[i].GetPatternCycle(0).GetInputVal();
//            auto& pat_cycle_pis2 = total_patterns[j].GetPatternCycle(0).GetInputVal();
//            assert(pat_cycle_pis1.size() == pat_cycle_pis2.size());
//            for (int k = 0; k < pat_cycle_pis1[0].size(); ++k) {
//                auto& pi_v1 = pat_cycle_pis1[0][k];
//                auto& pi_v2 = pat_cycle_pis2[0][k];
//                if (pi_v1 != pi_v2) {
//                    overlap_bit++;
//                }
//            }
//            overlap_bits.push_back(overlap_bit);
//            overlap_bit = 0;
//        }
//    }
//    double avg = std::accumulate(overlap_bits.begin(), overlap_bits.end(), 0.0);
//    avg = avg / overlap_bits.size();
//    std::cout << "avg overlap: " << avg << std::endl;
//    std::sort(overlap_bits.begin(), overlap_bits.end());
//    double min_ = overlap_bits[0];
//    double max_ = overlap_bits.back();
//    std::cout << "max overlap: " << max_ << std::endl;
//    std::cout << "min overlap: " << min_ << std::endl;
//    std::cout << "mean overlap: " << overlap_bits[overlap_bits.size()/2] << std::endl;
//    // ==================== CHECK OVERLAP END ==============
    std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second > b.second;
    });

    // 打印前十个最大的项
    for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }
}

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompareMini(
        int pattern_per_sim) {
    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    vector<std::pair<int, double>> single_time;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;

    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);

    auto start = std::chrono::steady_clock::now();
    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {
//        int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
//        std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];
        if (primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault)) {
//    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        std::vector<bool> old_vars_sign (solver->max_var_size_, false);
        saf_test_generator.SetSolver(solver);

        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();


        if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault,
                                                               fanin_mark)) {
            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
            vector<int> var_old2new, var_new2old;
            auto cnf_trans_start = std::chrono::steady_clock::now();
            saf_test_generator.WriteCNF(var_old2new, var_new2old);
            auto cnf_trans_end = std::chrono::steady_clock::now();
            cnf_trans_time += std::chrono::duration<double>(cnf_trans_end - cnf_trans_start).count();

            auto solve_start = std::chrono::steady_clock::now();
            system("../../compare_sat/minisat -verb=0 ../../compare_sat/cnf.txt ../../compare_sat/output.txt");
            test_res = saf_test_generator.ReadSATResultMini(old_vars_sign, var_new2old);
            auto solve_end = std::chrono::steady_clock::now();
            solve_time += std::chrono::duration<double>(solve_end - solve_start).count();
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), std::chrono::duration<double>(solve_end - solve_start).count()));




        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
            // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
            // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                auto pi_gate = prim_->GetPIGates()[pi_id];
                if (!fanin_mark[pi_gate->GetGId()] ) {
                    new_pattern[pi_id] = LOGIC_x;
                    //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                    //          ";
                    continue;
                }
                assert(pi_gate->GetGType() == GType::G_PI);
                auto pi_bit =
                        old_vars_sign[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;
                //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
                assert(pi_bit != LOGIC_x);
                new_pattern[pi_id] = pi_bit;
            }
            //      std::cout << std::endl;

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }

        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
        if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
  for (int i = 0; i < DETECT_FAULT.size(); ++i) {
      if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
          std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
          cnt++;
      }
  }
  if (cnt != 0) {
      std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
  }
  std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
  DETECT_FAULT.clear();
  for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
      if ((*it)->GetSAFStatus() != SAFStatus::DS) {
          std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
          " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
          exit(1);
      }
  }
#endif
        }
    }
    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
            (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
            (double)(fault_detect + fault_redundant) /
            (double)uncollapsed_saf_list.size();


    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;

    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }

    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();
    // ==================== VALIDATION END ==============
    std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second > b.second;
    });
    // 打印前十个最大的项
    for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }
}

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompareKis(
        int pattern_per_sim) {
    vector<vector<int>> compaction_patterns;
    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    vector<std::pair<int, double>> single_time;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;

    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);

    auto start = std::chrono::steady_clock::now();
    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {
//        int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
//        std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
//    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        std::vector<bool> old_vars_sign (solver->max_var_size_, false);
        saf_test_generator.SetSolver(solver);

        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();

        if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault,
                                                               fanin_mark)) {
            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
            vector<int> var_old2new, var_new2old;
            auto cnf_trans_start = std::chrono::steady_clock::now();
            saf_test_generator.WriteCNF(var_old2new, var_new2old);
            auto cnf_trans_end = std::chrono::steady_clock::now();
            cnf_trans_time += std::chrono::duration<double>(cnf_trans_end - cnf_trans_start).count();

            auto solve_start = std::chrono::steady_clock::now();
            system("../../compare_sat/kissat ../../compare_sat/cnf.txt > ../../compare_sat/output.txt");
            test_res = saf_test_generator.ReadSATResultKis(old_vars_sign, var_new2old);
            auto solve_end = std::chrono::steady_clock::now();
            solve_time += std::chrono::duration<double>(solve_end - solve_start).count();
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), std::chrono::duration<double>(solve_end - solve_start).count()));

        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
            // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
            // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                auto pi_gate = prim_->GetPIGates()[pi_id];
                if (!fanin_mark[pi_gate->GetGId()] ) {
                    new_pattern[pi_id] = LOGIC_x;
                    //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                    //          ";
                    continue;
                }
                assert(pi_gate->GetGType() == GType::G_PI);
                auto pi_bit =
                        old_vars_sign[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;
                //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
                assert(pi_bit != LOGIC_x);
                new_pattern[pi_id] = pi_bit;
            }
            //      std::cout << std::endl;

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }

        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
        if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);

            compaction_patterns.emplace_back(atpg_patterns[0]);

            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
for (int i = 0; i < DETECT_FAULT.size(); ++i) {
  if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
      std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
      cnt++;
  }
}
if (cnt != 0) {
  std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
}
std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
DETECT_FAULT.clear();
for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
  if ((*it)->GetSAFStatus() != SAFStatus::DS) {
      std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
      " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
      exit(1);
  }
}
#endif
        }
    }
    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
            (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
            (double)(fault_detect + fault_redundant) /
            (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;

    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }

    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();
    // ==================== VALIDATION END ==============
    std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second > b.second;
    });
    // 打印前十个最大的项
    for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }


    std::cout<<"compaction_pattern_size is: "<<compaction_patterns.size()<<std::endl;
    //start compaction step_1
    vector<vector<int>> new_compacted_patterns;
    new_compacted_patterns.clear();
    vector<int> mid_pattern;
    vector<int> fanin_mark(prim_->NumPIs(),0);
    vector<int> pattern_compcted_mark(compaction_patterns.size(),0);
    for (int i = 0; i < compaction_patterns.size(); ++i) {
        if (pattern_compcted_mark[i] == 0){
            pattern_compcted_mark[i] = 1;
            for (int j = 0; j < compaction_patterns[i].size(); ++j) {
  mid_pattern.emplace_back(compaction_patterns[i][j]);
            }
            for (int j = i + 1; j < compaction_patterns.size(); ++j) {
  int conflict_flag = 0;
  for (int k = 0; k < compaction_patterns[j].size(); ++k) {
      if (mid_pattern[k] != compaction_patterns[j][k] && mid_pattern[k] != 2 && compaction_patterns[j][k] != 2){
                        conflict_flag = 1;
      }
  }
  if (conflict_flag == 0){
      pattern_compcted_mark[j] = 1;
      for (int k = 0; k < mid_pattern.size(); ++k) {
                        if (mid_pattern[k] == 2){
                              mid_pattern[k] = compaction_patterns[j][k];
                        }
      }
  }

            }
        }
        if (mid_pattern.size() > 0){
            new_compacted_patterns.emplace_back(mid_pattern);
        }

        /*for (int j = 0; j < mid_pattern.size(); ++j) {
            std::cout<<" "<<mid_pattern[j]<<" ";
        }
        std::cout<<std::endl;*/
        mid_pattern.clear();
    }


    std::cout<<"compacted pattern number is : "<<new_compacted_patterns.size()<<std::endl;



}

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompareKis_dynamic(
    int pattern_per_sim) {
    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    vector<std::pair<int, double>> single_time;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;

    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);

    auto start = std::chrono::steady_clock::now();
    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {
        //        int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
        //        std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
            //    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        std::vector<bool> old_vars_sign (solver->max_var_size_, false);
        saf_test_generator.SetSolver(solver);

        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();

        if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault,
                                                               fanin_mark)) {
            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
            vector<int> var_old2new, var_new2old;
            auto cnf_trans_start = std::chrono::steady_clock::now();
            saf_test_generator.WriteCNF(var_old2new, var_new2old);
            auto cnf_trans_end = std::chrono::steady_clock::now();
            cnf_trans_time += std::chrono::duration<double>(cnf_trans_end - cnf_trans_start).count();

            auto solve_start = std::chrono::steady_clock::now();
            system("../../compare_sat/kissat ../../compare_sat/cnf.txt > ../../compare_sat/output.txt");
            test_res = saf_test_generator.ReadSATResultKis(old_vars_sign, var_new2old);
            auto solve_end = std::chrono::steady_clock::now();
            solve_time += std::chrono::duration<double>(solve_end - solve_start).count();
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), std::chrono::duration<double>(solve_end - solve_start).count()));

        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
        // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
        // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
  auto pi_gate = prim_->GetPIGates()[pi_id];
  if (!fanin_mark[pi_gate->GetGId()] ) {
      new_pattern[pi_id] = LOGIC_x;
      //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
      //          ";
      continue;
  }
  assert(pi_gate->GetGType() == GType::G_PI);
  auto pi_bit =
      old_vars_sign[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;
  //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
  assert(pi_bit != LOGIC_x);
  new_pattern[pi_id] = pi_bit;
            }

            //start dynaamic compaction
            //for (int i = pf_id + 1; i < num_rest_faults; ++i) {
            int compact_fail_count = 0;
            int continue_flag = 0;
  for (int i = pf_id + 1; i < num_rest_faults && continue_flag != 1; ++i) {

  continue_flag = 0;
  auto *primary_fault_2 = rest_faults[i];
  if ((primary_fault_2->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault_2))
      && primary_fault_2->GetSAFStatus() != SAFStatus::PT && primary_fault_2->GetSAFStatus() != SAFStatus::PU) {
      //    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
      continue;
  }

  int test_res_2;
  auto solver_2 = new Solver;
  solver_2->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
  std::vector<bool> old_vars_sign_2 (solver_2->max_var_size_, false);
  saf_test_generator.SetSolver(solver_2);

  std::vector<int> fanin_mark_2(prim_->GetPrimNetlist().size(), 0);


  if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault_2,
                                                         fanin_mark_2)) {

      for (int j = 0; j < new_pattern.size(); ++j) {
                        if (new_pattern[j] == 1){
                              solver_2->add(prim_->GetPIGates()[j]->GetGId() + 1);
                              solver_2->add(0);
                        }
                        else if(new_pattern[j] == 0){
                              solver_2->add(-1 * (prim_->GetPIGates()[j]->GetGId() + 1));
                              solver_2->add(0);
                        }

      }

      vector<int> var_old2new_2, var_new2old_2;
      saf_test_generator.WriteCNF(var_old2new_2, var_new2old_2);


      system("../../compare_sat/kissat ../../compare_sat/cnf.txt > ../../compare_sat/output.txt");
      test_res_2 = saf_test_generator.ReadSATResultKis(old_vars_sign_2, var_new2old_2);

      if (test_res_2 == 10){
                        for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                              auto pi_gate = prim_->GetPIGates()[pi_id];
                              assert(pi_gate->GetGType() == GType::G_PI);

                              auto pi_bit =
                                  old_vars_sign_2[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;

                              if (new_pattern[pi_id] == LOGIC_x && (pi_bit == LOGIC_0 || pi_bit == LOGIC_1)){
                                new_pattern[pi_id] = pi_bit;
                              }

                        }
                        compact_fail_count = 0;
      }else{
                        compact_fail_count++;
      }

      var_new2old_2.clear();
      var_old2new_2.clear();


  }


  if (compact_fail_count > 50){
      continue_flag = 1;
  }
            }







            //      std::cout << std::endl;

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }

        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
        if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
            for (int i = 0; i < DETECT_FAULT.size(); ++i) {
  if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
      std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
      cnt++;
  }
            }
            if (cnt != 0) {
  std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
            }
            std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
            DETECT_FAULT.clear();
            for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
  if ((*it)->GetSAFStatus() != SAFStatus::DS) {
      std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
          " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
      exit(1);
  }
            }
#endif
        }
    }
    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
        (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
        (double)(fault_detect + fault_redundant) /
        (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;

    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }

    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();
    // ==================== VALIDATION END ==============
    std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second > b.second;
              });
    // 打印前十个最大的项
    for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }
}

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompareCal(
        int pattern_per_sim) {
    vector<vector<int>> compaction_patterns;
    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    vector<std::pair<int, double>> single_time;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;

    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);

    auto start = std::chrono::steady_clock::now();
    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {
//        int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
//        std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
//    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        std::vector<bool> old_vars_sign (solver->max_var_size_, false);
        saf_test_generator.SetSolver(solver);

        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();

        if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault,
                                                               fanin_mark)) {
            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
            vector<int> var_old2new, var_new2old;
            auto cnf_trans_start = std::chrono::steady_clock::now();
            saf_test_generator.WriteCNF(var_old2new, var_new2old);
            auto cnf_trans_end = std::chrono::steady_clock::now();
            cnf_trans_time += std::chrono::duration<double>(cnf_trans_end - cnf_trans_start).count();

            auto solve_start = std::chrono::steady_clock::now();
            system("../../compare_sat/cadical ../../compare_sat/cnf.txt > ../../compare_sat/output.txt");
            test_res = saf_test_generator.ReadSATResultCal(old_vars_sign, var_new2old);
            auto solve_end = std::chrono::steady_clock::now();
            solve_time += std::chrono::duration<double>(solve_end - solve_start).count();
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), std::chrono::duration<double>(solve_end - solve_start).count()));

        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
            // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
            // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                auto pi_gate = prim_->GetPIGates()[pi_id];
                if (!fanin_mark[pi_gate->GetGId()] ) {
                    new_pattern[pi_id] = LOGIC_x;
                    //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                    //          ";
                    continue;
                }
                assert(pi_gate->GetGType() == GType::G_PI);
                auto pi_bit =
                        old_vars_sign[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;
                //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
                assert(pi_bit != LOGIC_x);
                new_pattern[pi_id] = pi_bit;
            }
            //      std::cout << std::endl;

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }

        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
            if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);

            compaction_patterns.emplace_back(atpg_patterns[0]);

            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
            for (int i = 0; i < DETECT_FAULT.size(); ++i) {
                if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
                    std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
                    cnt++;
                }
            }
            if (cnt != 0) {
                std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
            }
            std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
            DETECT_FAULT.clear();
            for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
                if ((*it)->GetSAFStatus() != SAFStatus::DS) {
                    std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
                              " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
                    exit(1);
                }
            }
#endif
        }
    }
    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
            (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
            (double)(fault_detect + fault_redundant) /
            (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;

    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }

    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();
    // ==================== VALIDATION END ==============
     std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second > b.second;
    });
    // 打印前十个最大的项
    for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }

    std::cout<<"compaction_pattern_size is: "<<compaction_patterns.size()<<std::endl;
    //start compaction step_1
    vector<vector<int>> new_compacted_patterns;
    new_compacted_patterns.clear();
    vector<int> mid_pattern;
    vector<int> fanin_mark(prim_->NumPIs(),0);
    vector<int> pattern_compcted_mark(compaction_patterns.size(),0);
    for (int i = 0; i < compaction_patterns.size(); ++i) {
        if (pattern_compcted_mark[i] == 0){
            pattern_compcted_mark[i] = 1;
            for (int j = 0; j < compaction_patterns[i].size(); ++j) {
                mid_pattern.emplace_back(compaction_patterns[i][j]);
            }
            for (int j = i + 1; j < compaction_patterns.size(); ++j) {
                int conflict_flag = 0;
                for (int k = 0; k < compaction_patterns[j].size(); ++k) {
                    if (mid_pattern[k] != compaction_patterns[j][k] && mid_pattern[k] != 2 && compaction_patterns[j][k] != 2){
                        conflict_flag = 1;
                    }
                }
                if (conflict_flag == 0){
                    pattern_compcted_mark[j] = 1;
                    for (int k = 0; k < mid_pattern.size(); ++k) {
                        if (mid_pattern[k] == 2){
                              mid_pattern[k] = compaction_patterns[j][k];
                        }
                    }
                }

            }
        }
        if (mid_pattern.size() > 0){
            new_compacted_patterns.emplace_back(mid_pattern);
        }

        /*for (int j = 0; j < mid_pattern.size(); ++j) {
            std::cout<<" "<<mid_pattern[j]<<" ";
        }
        std::cout<<std::endl;*/
        mid_pattern.clear();
    }


    std::cout<<"compacted pattern number is : "<<new_compacted_patterns.size()<<std::endl;
}


void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsComparePotopara(
        int thread_num, int pattern_per_sim) {
    vector<vector<int>> compaction_patterns;
    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    vector<std::pair<int, double>> single_time;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;

    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);

    auto start = std::chrono::steady_clock::now();
    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {
//        int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
//        std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
//    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        std::vector<bool> old_vars_sign (solver->max_var_size_, false);
        saf_test_generator.SetSolver(solver);

        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();

        if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault,
                                                               fanin_mark)) {
            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
            vector<int> var_old2new, var_new2old;
            auto cnf_trans_start = std::chrono::steady_clock::now();
            saf_test_generator.WriteCNF(var_old2new, var_new2old);
            auto cnf_trans_end = std::chrono::steady_clock::now();
            cnf_trans_time += std::chrono::duration<double>(cnf_trans_end - cnf_trans_start).count();

            auto solve_start = std::chrono::steady_clock::now();
            std::string command1 = "../../compare_sat/PRS ../../compare_sat/cnf.txt --clause_sharing=1 --DCE=0 --preprocessor=1 --nThreads=";
            std::string command2 = " > ../../compare_sat/output.txt";
            std::string command = command1 + std::to_string(thread_num) + command2;
            std::cout<< command << std::endl;
            system(command.c_str());
            test_res = saf_test_generator.ReadSATResultPotopara(old_vars_sign, var_new2old);
            auto solve_end = std::chrono::steady_clock::now();
            solve_time += std::chrono::duration<double>(solve_end - solve_start).count();
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), std::chrono::duration<double>(solve_end - solve_start).count()));

        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
            // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
            // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                auto pi_gate = prim_->GetPIGates()[pi_id];
                if (!fanin_mark[pi_gate->GetGId()] ) {
                    new_pattern[pi_id] = LOGIC_x;
                    //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                    //          ";
                    continue;
                }
                assert(pi_gate->GetGType() == GType::G_PI);
                auto pi_bit =
                        old_vars_sign[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;
                //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
                assert(pi_bit != LOGIC_x);
                new_pattern[pi_id] = pi_bit;
            }
            //      std::cout << std::endl;

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }

        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
        if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);

            compaction_patterns.emplace_back(atpg_patterns[0]);

            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
        for (int i = 0; i < DETECT_FAULT.size(); ++i) {
            if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
                std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
                cnt++;
            }
        }
        if (cnt != 0) {
            std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
        }
        std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
        DETECT_FAULT.clear();
        for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
            if ((*it)->GetSAFStatus() != SAFStatus::DS) {
                std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
                          " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
                exit(1);
            }
        }
#endif
        }
    }
    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
            (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
            (double)(fault_detect + fault_redundant) /
            (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;

    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }

    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();
    // ==================== VALIDATION END ==============
    std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second > b.second;
              });
    // 打印前十个最大的项
    for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }

    std::cout<<"compaction_pattern_size is: "<<compaction_patterns.size()<<std::endl;
    //start compaction step_1
    vector<vector<int>> new_compacted_patterns;
    new_compacted_patterns.clear();
    vector<int> mid_pattern;
    vector<int> fanin_mark(prim_->NumPIs(),0);
    vector<int> pattern_compcted_mark(compaction_patterns.size(),0);
    for (int i = 0; i < compaction_patterns.size(); ++i) {
        if (pattern_compcted_mark[i] == 0){
            pattern_compcted_mark[i] = 1;
            for (int j = 0; j < compaction_patterns[i].size(); ++j) {
                mid_pattern.emplace_back(compaction_patterns[i][j]);
            }
            for (int j = i + 1; j < compaction_patterns.size(); ++j) {
                int conflict_flag = 0;
                for (int k = 0; k < compaction_patterns[j].size(); ++k) {
                    if (mid_pattern[k] != compaction_patterns[j][k] && mid_pattern[k] != 2 && compaction_patterns[j][k] != 2){
                        conflict_flag = 1;
                    }
                }
                if (conflict_flag == 0){
                    pattern_compcted_mark[j] = 1;
                    for (int k = 0; k < mid_pattern.size(); ++k) {
                        if (mid_pattern[k] == 2){
                            mid_pattern[k] = compaction_patterns[j][k];
                        }
                    }
                }

            }
        }
        if (mid_pattern.size() > 0){
            new_compacted_patterns.emplace_back(mid_pattern);
        }

        /*for (int j = 0; j < mid_pattern.size(); ++j) {
            std::cout<<" "<<mid_pattern[j]<<" ";
        }
        std::cout<<std::endl;*/
        mid_pattern.clear();
    }


    std::cout<<"compacted pattern number is : "<<new_compacted_patterns.size()<<std::endl;
}

void SATEngine::GenerateCompactTestSetForSAFListHard2DTfaultsCompareCal_dynamic(
    int pattern_per_sim) {
    int test_found = 0, no_test = 0, abort = 0;
    int num_curr_pattern = 0;
    int total_pattern_num = 0;
    double sim_time = 0;
    double cnf_time = 0;
    double cnf_trans_time = 0;
    double solve_time = 0;
    int num_faults = saf_list_->GetUncollapsedSAFs().size();
    auto uncollapsed_saf_list = saf_list_->GetUncollapsedSAFs();
    std::unordered_set<SAF *> success_targets;
    std::vector<SAF *> rest_faults;
    vector<std::pair<int, double>> single_time;
    for (auto fault : uncollapsed_saf_list) {
        if (fault->GetSAFStatus() == SAFStatus::UC) {
            rest_faults.emplace_back(fault);
        }
    }
    int num_rest_faults = (int)rest_faults.size();

#ifdef DEBUG_MODE
    std::vector<std::pair<int, SAF*>> DETECT_FAULT;
#endif
    CNFGenerator saf_test_generator(prim_);

    std::vector<Pattern> patterns;
    std::vector<Pattern> total_patterns;
    std::vector<std::vector<int>> atpg_patterns;
    std::vector<std::vector<int>> validation_patterns;
    std::vector<int> new_pattern(prim_->NumPIs());
    PatternParser pp;

    pp.SetupPrimNetlist(prim_);
    pp.SetPatternType(PatternType::COMB_PT);
    assert(prim_->GetNetlistType() == NetlistType::COMB_CIRCUIT);

    auto start = std::chrono::steady_clock::now();
    for (int pf_id = 0; pf_id < num_rest_faults; pf_id++) {
        //        int rand_idx = rand() % (num_rest_faults - pf_id) + pf_id;
        //        std::swap(rest_faults[pf_id], rest_faults[rand_idx]);
        auto *primary_fault = rest_faults[pf_id];
        if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
            && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
            //    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
            continue;
        }
        //assert(primary_fault->GetSAFStatus() == SAFStatus::UC);
        std::cout << "At Fault: " << primary_fault->GetSAFName() << " "
                  << static_cast<int>(primary_fault->GetSAFType()) << std::endl;

        auto solver = new Solver;
        solver->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
        std::vector<bool> old_vars_sign (solver->max_var_size_, false);
        saf_test_generator.SetSolver(solver);

        int test_res;
        std::vector<int> fanin_mark(prim_->GetPrimNetlist().size(), 0);
        auto cnf_start = std::chrono::steady_clock::now();

        if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault,
                                                               fanin_mark)) {
            auto cnf_end = std::chrono::steady_clock::now();
            cnf_time += std::chrono::duration<double>(cnf_end - cnf_start).count();
            vector<int> var_old2new, var_new2old;
            auto cnf_trans_start = std::chrono::steady_clock::now();
            saf_test_generator.WriteCNF(var_old2new, var_new2old);
            auto cnf_trans_end = std::chrono::steady_clock::now();
            cnf_trans_time += std::chrono::duration<double>(cnf_trans_end - cnf_trans_start).count();

            auto solve_start = std::chrono::steady_clock::now();
            system("../../compare_sat/cadical ../../compare_sat/cnf.txt > ../../compare_sat/output.txt");
            test_res = saf_test_generator.ReadSATResultCal(old_vars_sign, var_new2old);
            auto solve_end = std::chrono::steady_clock::now();
            solve_time += std::chrono::duration<double>(solve_end - solve_start).count();
            single_time.push_back(std::make_pair(primary_fault->GetSAFId(), std::chrono::duration<double>(solve_end - solve_start).count()));

        } else {
            test_res = 20;
        }

        // No test.
        if (test_res == 20) {
            primary_fault->SetSAFStatus(SAFStatus::RE);
            no_test++;
            //      std::cout << "No Test" << std::endl;
        }
        // Over backtrack.
        else if (test_res == 0) {
            primary_fault->SetSAFStatus(SAFStatus::ATPG_ABORT);
        }
        // Test found.
        else if (test_res == 10) {
            // Load test cube from ATPG circuit model.
            for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                auto pi_gate = prim_->GetPIGates()[pi_id];
                if (!fanin_mark[pi_gate->GetGId()] ) {
                    new_pattern[pi_id] = LOGIC_x;
                    //          std::cout << pi_gate->GetInstName() << ":" << 'x' << "
                    //          ";
                    continue;
                }
                assert(pi_gate->GetGType() == GType::G_PI);
                auto pi_bit =
                    old_vars_sign[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;
                //        std::cout << pi_gate->GetInstName() << ":" << pi_bit << " ";
                assert(pi_bit != LOGIC_x);
                new_pattern[pi_id] = pi_bit;
            }
            //start dynaamic compaction
            int compact_fail_count = 0;
            int continue_flag = 0;
            for (int i = pf_id + 1; i < num_rest_faults && continue_flag != 1; ++i) {
            //for (int i = pf_id + 1; i < num_rest_faults; ++i) {

                auto *primary_fault_2 = rest_faults[i];
                if ((primary_fault->GetSAFStatus() != SAFStatus::UC || success_targets.count(primary_fault))
                    && primary_fault->GetSAFStatus() != SAFStatus::PT && primary_fault->GetSAFStatus() != SAFStatus::PU) {
                    //    if (primary_fault->GetSAFStatus() != SAFStatus::UC ) {
                    continue;
                }

                int test_res_2;
                auto solver_2 = new Solver;
                solver_2->max_var_size_ = prim_->GetPrimNetlist().size() * 3;
                std::vector<bool> old_vars_sign_2 (solver_2->max_var_size_, false);
                saf_test_generator.SetSolver(solver_2);

                std::vector<int> fanin_mark_2(prim_->GetPrimNetlist().size(), 0);


                if (saf_test_generator.CollectDChainCircuitClauseSAFL2(primary_fault_2,
                                                                       fanin_mark_2)) {

                    for (int j = 0; j < new_pattern.size(); ++j) {
                        if (new_pattern[j] == 1){
                              solver_2->add(prim_->GetPIGates()[j]->GetGId() + 1);
                              solver_2->add(0);
                        }
                        else if(new_pattern[j] == 0){
                              solver_2->add(-1 * (prim_->GetPIGates()[j]->GetGId() + 1));
                              solver_2->add(0);
                        }

                    }

                    vector<int> var_old2new_2, var_new2old_2;
                    saf_test_generator.WriteCNF(var_old2new_2, var_new2old_2);


                    system("../../compare_sat/cadical ../../compare_sat/cnf.txt > ../../compare_sat/output.txt");
                    test_res_2 = saf_test_generator.ReadSATResultKis(old_vars_sign_2, var_new2old_2);

                    if (test_res_2 == 10){
                        for (int pi_id = 0; pi_id < prim_->GetPIGates().size(); pi_id++) {
                              auto pi_gate = prim_->GetPIGates()[pi_id];
                              assert(pi_gate->GetGType() == GType::G_PI);

                              auto pi_bit =
                                  old_vars_sign_2[pi_gate->GetGId() + 1] ? LOGIC_1 : LOGIC_0;

                              if (new_pattern[pi_id] == LOGIC_x && (pi_bit == LOGIC_0 || pi_bit == LOGIC_1)){
                                new_pattern[pi_id] = pi_bit;
                              }

                        }
                        compact_fail_count = 0;
                    }else {
                        compact_fail_count++;
                    }




                }

                if (compact_fail_count > 50){
                    continue_flag = 1;
                }

            }

            //      std::cout << std::endl;

            atpg_patterns.emplace_back(new_pattern);
            validation_patterns.emplace_back(new_pattern);

            success_targets.insert(primary_fault);
        } else {
            std::cout << "Unexpected test gen res." << std::endl;
            exit(137);
        }

        saf_test_generator.ResetSolver();
        delete solver;
#ifdef DEBUG_MODE
        DETECT_FAULT.emplace_back(std::make_pair(pf_id, primary_fault));
#endif
        num_curr_pattern = atpg_patterns.size();
#ifdef DEBUG_MODE
        if (num_curr_pattern == 1) {
#else
        if (num_curr_pattern == 1) {
#endif
            auto sim_start = std::chrono::steady_clock::now();
            total_pattern_num += num_curr_pattern;
            pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
            sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
            atpg_patterns.clear();
            total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
            patterns.clear();
            auto sim_end = std::chrono::steady_clock::now();
            sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
#ifdef DEBUG_MODE
            int cnt = 0;
            for (int i = 0; i < DETECT_FAULT.size(); ++i) {
                if(DETECT_FAULT[i].second->GetSAFStatus() != SAFStatus::DS) {
                    std::cout << "" << DETECT_FAULT[i].second->GetSAFGate()->GetGId() << " " << DETECT_FAULT[i].second->GetSAFName() << " " << int(DETECT_FAULT[i].second->GetSAFType()) << std::endl;
                    cnt++;
                }
            }
            if (cnt != 0) {
                std::cout << "Gen pattern not detect ATPG detect faults number: " << cnt << std::endl;
            }
            std::cout << "Newly detected faults number: " << DETECT_FAULT.size() << std::endl;
            DETECT_FAULT.clear();
            for (auto it = success_targets.begin(); it != success_targets.end(); ++it) {
                if ((*it)->GetSAFStatus() != SAFStatus::DS) {
                    std::cout << "Gen pattern not detect ATPG detect fault: " << (*it)->GetSAFGate()->GetGId() <<
                        " " << (*it)->GetSAFName() << " " << int((*it)->GetSAFType()) << std::endl;
                    exit(1);
                }
            }
#endif
        }
    }
    num_curr_pattern = static_cast<int>(atpg_patterns.size());
    if (num_curr_pattern) {
        auto sim_start = std::chrono::steady_clock::now();
        total_pattern_num += num_curr_pattern;
        pp.LoadInternalComPatternForSAF(atpg_patterns, patterns, false);
        sim_->RunSAFSimulation(patterns, uncollapsed_saf_list);
        atpg_patterns.clear();
        total_patterns.insert(total_patterns.end(), patterns.begin(), patterns.end());
        patterns.clear();
        auto sim_end = std::chrono::steady_clock::now();
        sim_time += std::chrono::duration<double>(sim_end - sim_start).count();
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration<double>(end - start);

    int fault_detect = 0;
    int fault_redundant = 0;
    int fault_not_detect = 0;
    double fault_coverage = 0;
    double test_coverage = 0;
    DumpSAFList();
    for (auto fault : saf_list_->GetUncollapsedSAFs()) {
        if (fault->GetSAFStatus() == SAFStatus::DS ||
            fault->GetSAFStatus() == SAFStatus::DI ||
            fault->GetSAFStatus() == SAFStatus::DT_CHAIN_TEST) {
            fault_detect++;
        } else if (fault->GetSAFStatus() == SAFStatus::RE) {
            fault_redundant++;
        } else if (fault->GetSAFStatus() == SAFStatus::UC ||
                   fault->GetSAFStatus() == SAFStatus::ATPG_ABORT ||
                   fault->GetSAFStatus() == SAFStatus::UO) {
            fault_not_detect++;
        }
    }
    fault_coverage =
        (double)fault_detect / (double)uncollapsed_saf_list.size();
    test_coverage =
        (double)(fault_detect + fault_redundant) /
        (double)uncollapsed_saf_list.size();

    std::cout << "Detected Faults: " << fault_detect << std::endl;
    std::cout << "No Test Faults: " << fault_redundant << std::endl;
    std::cout << "Abort Faults: " << fault_not_detect << std::endl;
    std::cout << "Pattern generated: " << total_pattern_num << std::endl;
    std::cout << "Fault Coverage: " << fault_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test Coverage: " << test_coverage * 100 << "%"
              << std::endl;
    std::cout << "Test time: " << duration.count() << "s." << std::endl;
    std::cout << "CNF generation time: " << cnf_time << "s." << std::endl;
    std::cout << "CNF transformation time: " << cnf_trans_time << "s." << std::endl;
    std::cout << "Solver time: " << solve_time << "s." << std::endl;
    std::cout << "Sim time: " << sim_time << "s." << std::endl;

    // ==================== VALIDATION ===========================
    std::cout << std::endl;
    std::cout << "----------Validate test patterns-----------" << std::endl;
    for (auto fault : uncollapsed_saf_list) {
        fault->SetSAFStatus(SAFStatus::UC);
    }

    pp.LoadInternalComPatternForSAF(validation_patterns, patterns, false);
    sim_->RunSAFSimulation(total_patterns, uncollapsed_saf_list);
    validation_patterns.clear();
    patterns.clear();
    // ==================== VALIDATION END ==============
    std::sort(single_time.begin(), single_time.end(),
              [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.second > b.second;
              });
    // 打印前十个最大的项
    for (int i = 0; i < std::min(100, static_cast<int>(single_time.size())); ++i) {
        std::cout << "Item " << i + 1 << ": (" << single_time[i].first << ", " << single_time[i].second << ")" << std::endl;
    }
}



void SATEngine::CollectPoByDistance(Gate* gate, std::vector<Gate*>& ordered_pos,
                                    std::vector<int>& fault_cone) {
  int level = 0;
  std::queue<Gate*> que;
  que.emplace(gate);
  fault_cone[gate->GetGId()] = 1;
  while (!que.empty()) {
    auto front = que.front();
    que.pop();

    if (front->GetGType() == GType::G_PO) {
      ordered_pos.emplace_back(front);
      continue;
    }

    for (auto fanout : front->FanoutGates()) {
      if (!fault_cone[fanout->GetGId()]) {
        que.emplace(fanout);
        fault_cone[fanout->GetGId()] = 1;
      }
    }
  }

  std::sort(ordered_pos.begin(), ordered_pos.end(),
            [](Gate* a, Gate* b) { return a->GetDPI() > b->GetDPI(); });
}



/**
 * @brief
 * dump fault list
 * @param report
 */
void SATEngine::DumpSAFList() {
  std::ofstream dump_ofs;
  //  RemoveDir(dump_path_);
  if (access(dump_path_.c_str(), F_OK) == -1) {
    mkdir(dump_path_.c_str(), S_IRWXO | S_IRWXG | S_IRWXU);
  }
  std::string saf_dump_path = dump_path_ + "/ictest_atpg.flist";
  dump_ofs.open(saf_dump_path);
  if (!dump_ofs.is_open()) {
    LOG_ASSERT(false, "open Dump " + saf_dump_path + " fail")
  }

  for (auto fault : saf_list_->GetUncollapsedSAFs()) {
    if (fault->GetSAFType() == SAFType::SA0) {
      dump_ofs << "sa0"
               << "  ";
    } else {
      dump_ofs << "sa1"
               << "  ";
    }

    switch (fault->GetSAFStatus()) {
    case SAFStatus::DS:
      dump_ofs << "DS"
               << "  ";
      break;
    case SAFStatus::UC:
      dump_ofs << "UC"
               << "  ";
      break;
    case SAFStatus::UO:
      dump_ofs << "UO"
               << "  ";
      break;
    case SAFStatus::ATPG_ABORT:
      dump_ofs << "AB"
               << "  ";
      break;
    case SAFStatus::RE:
      dump_ofs << "RE"
               << "  ";
      break;
    case SAFStatus::TI:
      dump_ofs << "TI"
               << "  ";
      break;
    case SAFStatus::BL:
      dump_ofs << "BL"
               << "  ";
      break;
    case SAFStatus::UU:
      dump_ofs << "UU"
               << "  ";
      break;
    case SAFStatus::AU:
      dump_ofs << "AU"
               << "  ";
      break;
    case SAFStatus::DI:
      dump_ofs << "DI"
               << "  ";
      break;
    case SAFStatus::PT:
      dump_ofs << "PT"
               << "  ";
      break;
    case SAFStatus::DT_CHAIN_TEST:
      dump_ofs << "CT"
               << "  ";
      break;
    case SAFStatus::CONFLICT_FREE:
      dump_ofs << "CF"
               << "  ";
      break;
    default:
      LOG_ERROR("no such fault status");
    }

    auto fault_name = fault->GetSAFName();
    if (fault_name.find('=') != -1) {
      fault_name.erase(0, 2);
      int remove_idx = fault_name.find('=');
      fault_name.erase(remove_idx, 2);
      fault_name.insert(remove_idx, 1, '/');
    }
    dump_ofs << fault_name << std::endl;
  }
  dump_ofs.close();
}



}  // namespace ictest