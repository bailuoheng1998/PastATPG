//
// Created by tianpengyu on 23-4-4.
//

#include "atpg/CNFGenerator.h"

namespace ictest {

int Solver::solve()
{
    try {

        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", 0, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", 0, IntRange(0, INT32_MAX));
        BoolOption   strictp("MAIN", "strict", "Validate DIMACS header during parsing.", false);


        Minisat::SATSolver S;
        double initial_time = cpuTime();

        S.verbosity = verb;

        // Try to set resource limits:
        if (cpu_lim != 0) limitTime(cpu_lim);
        if (mem_lim != 0) limitMemory(mem_lim);

        auto cnf_start = std::chrono::steady_clock::now();

        // add variables and clauses
        vector<int> old_clause;
        // record a map from old2new and new2old
        vector<int> var_old2new (max_var_size_, -1);
        vector<int> var_new2old (max_var_size_, -1);
        for (int i = 0; i < cnf_var_intv_.size(); i++) {
            if (cnf_var_intv_[i] == 0) {
                // get a single clause
                vec<Lit> tmp_v;
                for (int j = 0; j < old_clause.size(); ++j) {
                    int abs = 0;
                    bool sign = false;
                    if (old_clause[j] >= 0) {
                        abs = old_clause[j];
                    } else {
                        abs = ::abs(old_clause[j]);
                        sign = true;
                    }
                    if (var_old2new[abs] != -1) {
                        Var new_var = var_old2new[abs];
                        Lit tmp = mkLit(new_var, sign);
                        tmp_v.push(tmp);
                    } else {
                        Var new_var = S.newVar();
                        var_new2old[new_var] = abs;
                        var_old2new[abs] = new_var;
                        Lit tmp = mkLit(new_var, sign);
                        tmp_v.push(tmp);
                    }
                }
                old_clause.clear();
                S.addClause(tmp_v);
                continue;
            } else {
                old_clause.push_back(cnf_var_intv_[i]);
            }
        }

//        if (S.verbosity > 0){
//            printf("============================[ Problem Statistics ]=============================\n");
//            printf("|                                                                             |\n"); }
//
//        if (S.verbosity > 0){
//            printf("|  Number of variables:  %12d                                         |\n", S.nVars());
//            printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
//
//        double parsed_time = cpuTime();
//        if (S.verbosity > 0){
//            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
//            printf("|                                                                             |\n"); }
        auto cnf_end = std::chrono::steady_clock::now();
        cnf_trans_time_ += std::chrono::duration<double>(cnf_end - cnf_start).count();
        auto solve_start = std::chrono::steady_clock::now();
        bool ret = S.solve();
        auto solve_end = std::chrono::steady_clock::now();
        solve_time_ += std::chrono::duration<double>(solve_end - solve_start).count();

        std::string res_now;
        if (S.verbosity > 0){
//            S.printStats();
            printf("\n"); }
        res_now += (ret ? "SATISFIABLE\n" : !ret ? "UNSATISFIABLE\n" : "INDETERMINATE\n");

        if (ret){
            res_now += "SAT\n";
            printf("SATISFIABLE\n");
            for (int i = 0; i < S.nVars(); i++)
                if (S.model[i] != l_Undef)
//                    res_now += ((i==0)?"":" ") + ((S.model[i]==l_True)?"":"-", i+1);
            res_now += " 0\n";
            // record pi bits
            var_vals_ = vector(max_var_size_, -1);
            for (int i = 0; i < var_new2old.size(); ++i) {
                if (var_new2old[i] != -1) {
                    lbool sat = S.modelValue(i);
                    if (sat == l_True) {
                        var_vals_[var_new2old[i]] = 1;
                    } else if (sat == l_False) {
                        var_vals_[var_new2old[i]] = 0;
                    }
                }
            }
            // record X value
            for (int i = 0; i < S.nVars(); i++) {
                lbool val = S.modelValue(i);
                if(val == l_Undef){
                    x_vars_.insert(var_new2old[i]);
                }
            }
            // dump
//            for (int i = 0; i < S.nVars(); i++) {
//                lbool val = S.modelValue(i);
//                if(val == l_Undef){
//                    printf("x%d ",i+1);
//                }
//                else{
//                    if(val == l_False)
//                        printf("-");
//                    printf("%d ",i+1);
//                }
//            }
//            printf("0\n");
        } else if (!ret) {
            printf("UNSATISFIABLE\n");
            res_now += "UNSAT\n";
        } else {
            printf("INDETSATISFIABLE\n");
            res_now += "INDET\n";
        }

        res.push_back(res_now);





#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret ? 10 : !ret ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}

int Solver::solve_help_checkx()
{
    try {
        // helper to check x bits
        vector<vector<int>> cnf_new;

        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", 0, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", 0, IntRange(0, INT32_MAX));
        BoolOption   strictp("MAIN", "strict", "Validate DIMACS header during parsing.", false);


        Minisat::SATSolver S;
        double initial_time = cpuTime();

        S.verbosity = verb;

        // Try to set resource limits:
        if (cpu_lim != 0) limitTime(cpu_lim);
        if (mem_lim != 0) limitMemory(mem_lim);

        auto cnf_start = std::chrono::steady_clock::now();

        // add variables and clauses
        vector<int> old_clause;
        // record a map from old2new and new2old
        vector<int> var_old2new (max_var_size_, -1);
        vector<int> var_new2old (max_var_size_, -1);
        for (int i = 0; i < cnf_var_intv_.size(); i++) {
            if (cnf_var_intv_[i] == 0) {
                // get a single clause
                vec<Lit> tmp_v;
                // get a helper single clause

                for (int j = 0; j < old_clause.size(); ++j) {
                    int abs = 0;
                    bool sign = false;
                    if (old_clause[j] >= 0) {
                        abs = old_clause[j];
                    } else {
                        abs = ::abs(old_clause[j]);
                        sign = true;
                    }
                    if (var_old2new[abs] != -1) {
                        Var new_var = var_old2new[abs];
                        Lit tmp = mkLit(new_var, sign);
                        tmp_v.push(tmp);
                        // if (sign) {
                        //    cla_new.push_back(-tmp.x);
                        // } else {
                        //    cla_new.push_back(tmp.x);
                        // }
                    } else {
                        Var new_var = S.newVar();
                        var_new2old[new_var] = abs;
                        var_old2new[abs] = new_var;
                        Lit tmp = mkLit(new_var, sign);
                        tmp_v.push(tmp);
                        // if (sign) {
                        //    cla_new.push_back(-tmp.x);
                        // } else {
                        //    cla_new.push_back(tmp.x);
                        // }
                    }
                }
                old_clause.clear();
                S.addClause(tmp_v);
                vector<int> cla_new;
                cla_new.clear();
                // printf("z debug * ");
                for(int zi=0; zi< tmp_v.size(); ++zi){
                  Minisat::Lit zl = tmp_v[zi];
                  int zv = zl.x/2 + 1;
                  if(sign(zl))
                    zv = -zv;
                  cla_new.push_back(zv);
                  // printf("%d ",zv);
                }
                // printf("0\n");

                cnf_new.push_back(cla_new);
                continue;
            } else {
                old_clause.push_back(cnf_var_intv_[i]);
            }
        }

//        if (S.verbosity > 0){
//            printf("============================[ Problem Statistics ]=============================\n");
//            printf("|                                                                             |\n"); }
//
//        if (S.verbosity > 0){
//            printf("|  Number of variables:  %12d                                         |\n", S.nVars());
//            printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
//
//        double parsed_time = cpuTime();
//        if (S.verbosity > 0){
//            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
//            printf("|                                                                             |\n"); }
        auto cnf_end = std::chrono::steady_clock::now();
        cnf_trans_time_ += std::chrono::duration<double>(cnf_end - cnf_start).count();
        auto solve_start = std::chrono::steady_clock::now();
        bool ret = S.solve();
        auto solve_end = std::chrono::steady_clock::now();
        solve_time_ += std::chrono::duration<double>(solve_end - solve_start).count();

        std::string res_now;
        if (S.verbosity > 0){
//            S.printStats();
            // ("\n");
             }
        res_now += (ret ? "SATISFIABLE\n" : !ret ? "UNSATISFIABLE\n" : "INDETERMINATE\n");

        if (ret){
            res_now += "SAT\n";
            //printf("SATISFIABLE\n");
            for (int i = 0; i < S.nVars(); i++)
                if (S.model[i] != l_Undef)
//                    res_now += ((i==0)?"":" ") + ((S.model[i]==l_True)?"":"-", i+1);
            res_now += " 0\n";
            // record pi bits
            var_vals_ = vector(max_var_size_, -1);
            for (int i = 0; i < var_new2old.size(); ++i) {
                if (var_new2old[i] != -1) {
                    lbool sat = S.modelValue(i);
                    if (sat == l_True) {
                        var_vals_[var_new2old[i]] = 1;
                    } else if (sat == l_False) {
                        var_vals_[var_new2old[i]] = 0;
                    }
                }
            }
            // record X value
            for (int i = 0; i < S.nVars(); i++) {
                lbool val = S.modelValue(i);
                if(val == l_Undef){
                    x_vars_.insert(var_new2old[i]);
                }
            }
//             dump
            for (int i = 0; i < S.nVars(); i++) {
                lbool val = S.modelValue(i);
                if(val == l_Undef){
                    //printf("x%d ",i+1);
                }
                else{
                    if(val == l_False){
                        //printf("-");
                    }

                    //printf("%d ",i+1);
                }
            }
            //printf("0\n");
        } else if (!ret) {
            //printf("UNSATISFIABLE\n");
            res_now += "UNSAT\n";
        } else {
            //printf("INDETSATISFIABLE\n");
            res_now += "INDET\n";
        }

        res.push_back(res_now);

        // check x bit
        if (ret) {
            int var_new_x = 0;
            bool check = CheckXbits(cnf_new, S);
            if (check) {
                //std::cout << "check satisfied" << std::endl;
            } else {
                //std::cout << "check unsatisfied" << std::endl;
            }
        }
        int var_new_x = 0;
        for (int i = 0; i < S.model.size(); ++i) {
            auto now_val = S.model[i];
            if (S.model[i] != l_Undef) {
                S.model[i] = l_Undef;
                bool check = CheckXbits(cnf_new, S);
                if (!check) {
                    S.model[i] = now_val;
                } else {
                    var_new_x++;
                }
            }
        }
        //std::cout << "var_new_x: " << var_new_x << std::endl;


        /*std::cout<<"var new2old is : ";
        for (int i = 0; i < var_new2old.size(); ++i) {
            std::cout<<" "<<var_new2old[i]<<" ";
        }
        std::cout<<std::endl;*/

        /*std::cout<<"var old2new is : ";
        for (int i = 0; i < var_old2new.size(); ++i) {
            std::cout<<" "<<var_old2new[i]<<" ";
        }
        std::cout<<std::endl;

        std::cout<<"var new2old is : ";
        for (int i = 0; i < var_new2old.size(); ++i) {
            std::cout<<" "<<var_new2old[i]<<" ";
        }
        std::cout<<std::endl;*/


#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret ? 10 : !ret ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}

int Solver::solve(ATPGEngine* d_test_engine, SAF* fault)
{
    try {

        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", 0, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", 0, IntRange(0, INT32_MAX));
        BoolOption   strictp("MAIN", "strict", "Validate DIMACS header during parsing.", false);


        Minisat::SATSolver S;
        double initial_time = cpuTime();

        S.verbosity = verb;

        // Try to set resource limits:
        if (cpu_lim != 0) limitTime(cpu_lim);
        if (mem_lim != 0) limitMemory(mem_lim);

        auto cnf_start = std::chrono::steady_clock::now();

        // add variables and clauses
        vector<int> old_clause;
        // record a map from old2new and new2old
        vector<int> var_old2new (max_var_size_, -1);
        vector<int> var_new2old (max_var_size_, -1);
        for (int i = 0; i < cnf_var_intv_.size(); i++) {
            if (cnf_var_intv_[i] == 0) {
                // get a single clause
                vec<Lit> tmp_v;
                for (int j = 0; j < old_clause.size(); ++j) {
                    int abs = 0;
                    bool sign = false;
                    if (old_clause[j] >= 0) {
                        abs = old_clause[j];
                    } else {
                        abs = ::abs(old_clause[j]);
                        sign = true;
                    }
                    if (var_old2new[abs] != -1) {
                        Var new_var = var_old2new[abs];
                        Lit tmp = mkLit(new_var, sign);
                        tmp_v.push(tmp);
                    } else {
                        Var new_var = S.newVar();
                        var_new2old[new_var] = abs;
                        var_old2new[abs] = new_var;
                        Lit tmp = mkLit(new_var, sign);
                        tmp_v.push(tmp);
                    }
                }
                old_clause.clear();
                S.addClause(tmp_v);
                continue;
            } else {
                old_clause.push_back(cnf_var_intv_[i]);
            }
        }

//        if (S.verbosity > 0){
//            printf("============================[ Problem Statistics ]=============================\n");
//            printf("|                                                                             |\n"); }
//
//        if (S.verbosity > 0){
//            printf("|  Number of variables:  %12d                                         |\n", S.nVars());
//            printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
//
//        double parsed_time = cpuTime();
//        if (S.verbosity > 0){
//            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
//            printf("|                                                                             |\n"); }
        auto cnf_end = std::chrono::steady_clock::now();
        cnf_trans_time_ += std::chrono::duration<double>(cnf_end - cnf_start).count();
        vector<double> initial_activity (S.nVars(), 0);
//        RecordGoodBadDeterminedNumber(var_old2new, var_new2old, S, initial_activity, d_test_engine, fault);
//        RecordDecisionTrace(var_old2new, var_new2old, S, initial_activity, d_test_engine, fault);
//        RecordActPI2PO(var_old2new, var_new2old, S, initial_activity);
//        RecordActFS2POAndPI2FS(var_old2new, var_new2old, S, initial_activity);
//        RecordActFS2POAndPO2PI_inv(var_old2new, var_new2old, S, initial_activity);
//        RecordActFS2POAndPO2PI(var_old2new, var_new2old, S, initial_activity);
        RecordActFS2NEAR(var_old2new, var_new2old, S, initial_activity);
//        RecordActPO2PI(var_old2new, var_new2old, S, initial_activity);
        auto solve_start = std::chrono::steady_clock::now();
        for (int i = 0; i < S.nVars(); ++i) {
            S.activity[i] = initial_activity[i];
        }
//        vector<double> activities_before;
//        for (int i = 0; i < S.nVars(); ++i) {
//            activities_before.push_back(S.activity[i]);
//        }
        bool ret = S.solve();
        auto solve_end = std::chrono::steady_clock::now();
//        vector<double> activities_after;
//        for (int i = 0; i < S.nVars(); ++i) {
//            activities_after.push_back(S.activity[i]);
//        }
//        std::cout << "max: " << max << std::endl;
        solve_time_ += std::chrono::duration<double>(solve_end - solve_start).count();

        std::string res_now;
        if (S.verbosity > 0){
//            S.printStats();
            printf("\n"); }
        res_now += (ret ? "SATISFIABLE\n" : !ret ? "UNSATISFIABLE\n" : "INDETERMINATE\n");

        if (ret){
            res_now += "SAT\n";
            printf("SATISFIABLE\n");
            for (int i = 0; i < S.nVars(); i++)
                if (S.model[i] != l_Undef)
//                    res_now += ((i==0)?"":" ") + ((S.model[i]==l_True)?"":"-", i+1);
            res_now += " 0\n";
            // record pi bits
            var_vals_ = vector(max_var_size_, -1);
            for (int i = 0; i < var_new2old.size(); ++i) {
                if (var_new2old[i] != -1) {
                    lbool sat = S.modelValue(i);
                    if (sat == l_True) {
                        var_vals_[var_new2old[i]] = 1;
                    } else if (sat == l_False) {
                        var_vals_[var_new2old[i]] = 0;
                    }
                }
            }
            // record X value
            for (int i = 0; i < S.nVars(); i++) {
                lbool val = S.modelValue(i);
                if(val == l_Undef){
                    x_vars_.insert(var_new2old[i]);
                }
            }
            // dump
//            for (int i = 0; i < S.nVars(); i++) {
//                lbool val = S.modelValue(i);
//                if(val == l_Undef){
//                    printf("x%d ",i+1);
//                }
//                else{
//                    if(val == l_False)
//                        printf("-");
//                    printf("%d ",i+1);
//                }
//            }
//            printf("0\n");
        } else if (!ret) {
            printf("UNSATISFIABLE\n");
            res_now += "UNSAT\n";
        } else {
            printf("INDETSATISFIABLE\n");
            res_now += "INDET\n";
        }

        res.push_back(res_now);


#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret ? 10 : !ret ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}



void Solver::RecordActFS2NEAR(vector<int> &var_old2new, vector<int> &var_new2old, Minisat::SATSolver& S, vector<double> &initial_activity) {
    int gate_num = max_var_size_ / 3;
    double max_score = 0;
    for (int i = gid_po2pi_gptr_.size() - 1; i > 0; --i) {
        auto gate = gid_po2pi_gptr_[i];
        int gid = gate->GetGId();
        int g_var = gid + 1;
        int b_var = gid + 1 + gate_num;
        int d_var = gid2dchain_vars_[gid];
        assert(max_dpi_ - gate->GetDPI() >= 0);
        double score = double (max_dpi_ - abs(fault_site_->GetDPI() - gate->GetDPI())) / double (max_dpi_) * 10;
        if(var_old2new[g_var] != -1) {
            initial_activity[var_old2new[g_var]] = score;
        }
        if(var_old2new[b_var] != -1) {
            initial_activity[var_old2new[b_var]] = score;
        }
        if(var_old2new[d_var] != -1) {
            initial_activity[var_old2new[d_var]] = score;
        }
        if (max_score < score) {
            max_score = score;
        }
    }
}

//bool Solver::CheckXbits(vector<vector<int>> &cnf_new, Minisat::SATSolver &S) {
//    bool total_res = true;
//    for (int i = 0; i < cnf_new.size(); ++i) {
//        bool cla_res = false;
//        for (int j = 0; j < cnf_new[i].size(); ++j) {
//            int lit = cnf_new[i][j];
//            int var = abs(lit);
//            if (S.model[var] == l_Undef) continue;
//            if ((lit > 0 && S.model[var] == l_True) || (lit < 0 && S.model[var] == l_False)) {
//                cla_res = true;
//                break;
//            }
//        }
//        total_res = total_res & cla_res;
//        if (!total_res) {
//            break;
//        }
//    }
//    return total_res;
//}

bool Solver::CheckXbits(vector<vector<int>> &cnf_new, Minisat::SATSolver &S){
  // printf("z debug * v");
  // for(int i=0; i<S.nVars(); ++i){
  //   if(S.model[i] == l_Undef){
  //     printf("X%d ",i);
  //   }else if(S.model[i] == l_True){
  //     printf(" %d ",i);
  //   }else if(S.model[i] == l_False){
  //     printf("-%d ",i);
  //   }
  // }printf("\n");

    for(auto zcls:cnf_new){
      // printf("z debug * ");
      int sat_ct = 0;
      for(auto zl:zcls){
        int v = abs(zl)-1;
        char c = 'X';
        if(zl>0){
          if(S.model[v]==l_True){
            c = 'T';
          }else if(S.model[v]==l_False){
            c = '-';
          }
        }else if(zl<0){
          if(S.model[v]==l_True){
            c = '-';
          }else if(S.model[v]==l_False){
            c = 'T';
          }
        }
        if(c == 'T')
          sat_ct ++;

        // printf("%d[%c] ", zl, c);
      }
      // printf("(%d)\t\t\n", sat_ct);
      if(sat_ct == 0){
        // printf("z debug * error!\n");
        return false;
      }

    }
    return true;
}



bool CNFGenerator::CollectDChainCircuitClauseSAFL2(SAF* saf_fault, std::vector<int>& fanin_mark) {



  //print the gate
//    for (int i = 0; i < prim_->GetPrimNetlist().size(); ++i) {
//        std::cout<<"gate id :"<<i<<std::endl;
//        std::cout<<"gate type is : "<< int(prim_->GetPrimNetlist()[i]->GetGType())<<std::endl;
//        std::cout<<"fanin number is :"<<prim_->GetPrimNetlist()[i]->FaninGates().size()<<std::endl;
//        std::cout<<"fanin gate_id is :";
//        for (int j = 0; j < prim_->GetPrimNetlist()[i]->FaninSize(); ++j) {
//            std::cout<<prim_->GetPrimNetlist()[i]->FaninGates()[j]->GetGId()<<"   ";
//        }
//        std::cout<<std::endl;
//
//        std::cout<<"fanout number is :"<<prim_->GetPrimNetlist()[i]->FanoutGates().size()<<std::endl;
//        std::cout<<"fanout gate_id is :";
//        for (int j = 0; j < prim_->GetPrimNetlist()[i]->FanoutSize(); ++j) {
//            std::cout<<prim_->GetPrimNetlist()[i]->FanoutGates()[j]->GetGId()<<"   ";
//        }
//        std::cout<<std::endl;
//        std::cout<<std::endl;
//
//    }

    log_info_ = false;


    //print fault gate and type
    //std::cout<<"fault gate is :"<<saf_fault->GetSAFGate()->GetGId()<<",and tpye is "<<int(saf_fault->GetSAFType())<<"id is : "<<saf_fault->GetSAFId()<<std::endl;


  auto faulty_gate = saf_fault->GetSAFGate();
  auto fault_type = saf_fault->GetSAFType();

  std::vector<int> fault_cone(prim_->GetPrimNetlist().size(), 0);
  std::vector<Gate*> reach_prim_outs;
  std::queue<Gate*> que;
  que.emplace(faulty_gate);
  fault_cone[faulty_gate->GetGId()] = 1;

  max_extra_var_id_ = 2 * prim_->GetPrimNetlist().size() + 1;

  if (log_info_) {
    std::cout << "Fault Cone:" << std::endl;
  }

  std::vector<int> d_chain_vars(prim_->GetPrimNetlist().size(), 0);
  std::vector<int> po_xors;

  while (!que.empty()) {
    auto front = que.front();
    que.pop();
    solver_->gid_fsite2po_.emplace_back(front->GetGId());

    if (log_info_) {
      std::cout << front->GetInstName() << std::endl;
    }

    // Collect differ variable clause.
    d_chain_vars[front->GetGId()] = max_extra_var_id_++;
    CollectL2XOR(nullptr, d_chain_vars[front->GetGId()],
                 {GetVarIdOfGoodGate(front), GetVarIdOfBadGate(front)});

    if (front->GetGType() == GType::G_PO) {
      reach_prim_outs.emplace_back(front);
      po_xors.emplace_back(d_chain_vars[front->GetGId()]);
      continue;
    }

    for (auto fanout : front->FanoutGates()) {
      if (!fault_cone[fanout->GetGId()]) {
        que.emplace(fanout);
        fault_cone[fanout->GetGId()] = 1;
      }
    }
  }

  if (reach_prim_outs.empty()) {
    return false;
  }

  if (log_info_) {
    std::cout << "reach po size: " << reach_prim_outs.size() << std::endl;
  }

  for (auto reach_po : reach_prim_outs) {
    que.emplace(reach_po);
    fanin_mark[reach_po->GetGId()] = 1;
  }

  if (log_info_) {
    std::cout << "Fan-in Cone:" << std::endl;
  }

  while (!que.empty()) {
    auto front = que.front();
    que.pop();
    solver_->gid_po2pi_.emplace_back(front->GetGId());
    solver_->gid_po2pi_gptr_.emplace_back(front);

    if (log_info_) {
      std::cout << front->GetInstName() << std::endl;
    }

    if (front == faulty_gate || !fault_cone[front->GetGId()]) {
      std::vector<int> fanin_vars;
      for (auto fanin : front->FaninGates()) {
        // Add branch stem.
        //        if (fanin->GetGType() == GType::G_BRH) {
        //          fanin = fanin->FaninGates()[0];
        //        }
        fanin_vars.emplace_back(GetVarIdOfGoodGate(fanin));
      }
      CollectGateL2(front, GetVarIdOfGoodGate(front), fanin_vars);
    } else {
      std::vector<int> good_fanin_vars;
      std::vector<int> bad_fanin_vars;
      for (auto fanin : front->FaninGates()) {
        // Add branch stem.
        //        if (fanin->GetGType() == GType::G_BRH) {
        //          fanin = fanin->FaninGates()[0];
        //        }
        good_fanin_vars.emplace_back(GetVarIdOfGoodGate(fanin));
        if (fault_cone[fanin->GetGId()]) {
          bad_fanin_vars.emplace_back(GetVarIdOfBadGate(fanin));
        } else {
          bad_fanin_vars.emplace_back(GetVarIdOfGoodGate(fanin));
        }
      }
      CollectGateL2(front, GetVarIdOfGoodGate(front), good_fanin_vars);
      CollectGateL2(front, GetVarIdOfBadGate(front), bad_fanin_vars);
    }

    /*if (front->FaninSize() > 0 && fault_cone[front->GetGId()] &&
        front != faulty_gate) {
      solver_->add(-d_chain_vars[front->GetGId()]);
      for (auto fanin : front->FaninGates()) {
        if (fault_cone[fanin->GetGId()]) {
          solver_->add(d_chain_vars[fanin->GetGId()]);
        }
      }
      solver_->add(0);
    }*/

    for (auto fanin : front->FaninGates()) {
      if (!fanin_mark[fanin->GetGId()]) {
        que.emplace(fanin);
        fanin_mark[fanin->GetGId()] = 1;
      }
    }
  }

  if (log_info_) {
    std::cout << "Fault Activation Clauses:" << std::endl;
  }

  // Fault activation clause.
  if (fault_type == SAFType::SA0) {
    solver_->add(GetVarIdOfGoodGate(faulty_gate));
    solver_->add(0);
    solver_->add(-GetVarIdOfBadGate(faulty_gate));
    solver_->add(0);
  } else {
    solver_->add(-GetVarIdOfGoodGate(faulty_gate));
    solver_->add(0);
    solver_->add(GetVarIdOfBadGate(faulty_gate));
    solver_->add(0);
  }

  if (log_info_) {
    std::cout << "Fault Observation Clauses:" << std::endl;
  }

  if (po_xors.size() > 1) {
    int final_signature = max_extra_var_id_++;
    CollectL2OR(nullptr, final_signature, po_xors);
    solver_->add(final_signature);
    solver_->add(0);
  } else {
    int final_signature = max_extra_var_id_++;
    CollectL2BUF(nullptr, final_signature, po_xors[0]);
    solver_->add(final_signature);
    solver_->add(0);
  }

  // Tied-gate clause.
  for (auto tie : prim_->GetTieGates()) {
    if (tie->GetGType() == GType::G_TIE0) {
      solver_->add(-GetVarIdOfGoodGate(tie));
      solver_->add(0);
    } else if (tie->GetGType() == GType::G_TIE1) {
      solver_->add(GetVarIdOfGoodGate(tie));
      solver_->add(0);
    } else {
      std::cout << "Unknown gate type: " << static_cast<int>(tie->GetGType())
                << std::endl;
      exit(137);
    }
  }

  solver_->gid2dchain_vars_ = d_chain_vars;
  return true;
}

bool CNFGenerator::CollectDChainCircuitClauseSAFL2TPI(
    SAF* saf_fault, std::vector<int>& fanin_mark) {
    log_info_ = false;
    std::set<Gate*> faulty_branch;
  auto faulty_gate = saf_fault->GetSAFGate();
  auto fault_type = saf_fault->GetSAFType();
  if (faulty_gate->GetGType() == GType::G_BRH) {
      assert(faulty_gate->FaninSize() == 1);
      faulty_branch.insert(faulty_gate);
//      for (int i = 0; i < faulty_gate->FaninGates()[0]->FanoutGates().size(); ++i) {
//          assert(faulty_gate->FaninGates()[0]->FanoutGates()[i]->GetGType() == GType::G_BRH);
//          faulty_branch.insert(faulty_gate->FaninGates()[0]->FanoutGates()[i]);
//      }
  }

  std::vector<int> fault_cone(prim_->GetPrimNetlist().size(), 0);
  std::vector<Gate*> reach_prim_outs;
  std::queue<Gate*> que;
  que.emplace(faulty_gate);
  fault_cone[faulty_gate->GetGId()] = 1;

  max_extra_var_id_ = 2 * prim_->GetPrimNetlist().size() + 1;

  if (log_info_) {
    std::cout << "Fault Cone:" << std::endl;
  }

  std::vector<int> d_chain_vars(prim_->GetPrimNetlist().size(), 0);
  std::vector<int> po_xors;

  while (!que.empty()) {
    auto front = que.front();
    que.pop();
    solver_->gid_fsite2po_.emplace_back(front->GetGId());

    if (log_info_) {
      std::cout << front->GetInstName() << std::endl;
    }

    // Collect differ variable clause.
    d_chain_vars[front->GetGId()] = max_extra_var_id_++;
    CollectL2XOR(nullptr, d_chain_vars[front->GetGId()],
                 {GetVarIdOfGoodGate(front), GetVarIdOfBadGate(front)});

    if (front->GetGType() == GType::G_PO) {
      reach_prim_outs.emplace_back(front);
      po_xors.emplace_back(d_chain_vars[front->GetGId()]);
      continue;
    }

    for (auto fanout : front->FanoutGates()) {
      if (!fault_cone[fanout->GetGId()]) {
          if (fanout->GetGType() == GType::G_BRH) {
              if (faulty_branch.count(fanout)) {
                  que.emplace(fanout);
                  fault_cone[fanout->GetGId()] = 1;
              } else {
                  que.emplace(fanout->FanoutGates()[0]);
                  fault_cone[fanout->FanoutGates()[0]->GetGId()] = 1;
              }
          } else {
               que.emplace(fanout);
               fault_cone[fanout->GetGId()] = 1;
          }

      }
    }
  }

  if (reach_prim_outs.empty()) {
    return false;
  }

  if (log_info_) {
    std::cout << "reach po size: " << reach_prim_outs.size() << std::endl;
  }

  for (auto reach_po : reach_prim_outs) {
    que.emplace(reach_po);
    fanin_mark[reach_po->GetGId()] = 1;
  }

  if (log_info_) {
    std::cout << "Fan-in Cone:" << std::endl;
  }

  while (!que.empty()) {
    auto front = que.front();
    que.pop();
    solver_->gid_po2pi_.emplace_back(front->GetGId());
    solver_->gid_po2pi_gptr_.emplace_back(front);

    if (log_info_) {
      std::cout << front->GetInstName() << std::endl;
    }

    if (front == faulty_gate || !fault_cone[front->GetGId()]) {
      std::vector<int> fanin_vars;
      for (auto fanin : front->FaninGates()) {
//         Add branch stem.
        if (fanin->GetGType() == GType::G_BRH && !faulty_branch.count(fanin)) {
          fanin = fanin->FaninGates()[0];
        }
        fanin_vars.emplace_back(GetVarIdOfGoodGate(fanin));
      }
      CollectGateL2(front, GetVarIdOfGoodGate(front), fanin_vars);
    } else {
      std::vector<int> good_fanin_vars;
      std::vector<int> bad_fanin_vars;
      for (auto fanin : front->FaninGates()) {
        // Add branch stem.
        if (fanin->GetGType() == GType::G_BRH && !faulty_branch.count(fanin)) {
          fanin = fanin->FaninGates()[0];
        }
        good_fanin_vars.emplace_back(GetVarIdOfGoodGate(fanin));
        if (fault_cone[fanin->GetGId()]) {
          bad_fanin_vars.emplace_back(GetVarIdOfBadGate(fanin));
        } else {
          bad_fanin_vars.emplace_back(GetVarIdOfGoodGate(fanin));
        }
      }
      CollectGateL2(front, GetVarIdOfGoodGate(front), good_fanin_vars);
      CollectGateL2(front, GetVarIdOfBadGate(front), bad_fanin_vars);
    }

    /*if (front->FaninSize() > 0 && fault_cone[front->GetGId()] &&
        front != faulty_gate) {
      solver_->add(-d_chain_vars[front->GetGId()]);
      for (auto fanin : front->FaninGates()) {
        if (fault_cone[fanin->GetGId()]) {
          solver_->add(d_chain_vars[fanin->GetGId()]);
        }
      }
      solver_->add(0);
    }*/

    for (auto fanin : front->FaninGates()) {
      if (!fanin_mark[fanin->GetGId()]) {
          if (fanin->GetGType() == GType::G_BRH) {
              if (faulty_branch.count(fanin)) {
                 que.emplace(fanin);
                 fanin_mark[fanin->GetGId()] = 1;
              } else {
                  que.emplace(fanin->FaninGates()[0]);
                  fanin_mark[fanin->FaninGates()[0]->GetGId()] = 1;
              }

          } else {
              que.emplace(fanin);
              fanin_mark[fanin->GetGId()] = 1;
          }
      }
    }
  }

  if (log_info_) {
    std::cout << "Fault Activation Clauses:" << std::endl;
  }

  // Fault activation clause.
  if (fault_type == SAFType::SA0) {
    solver_->add(GetVarIdOfGoodGate(faulty_gate));
    solver_->add(0);
    solver_->add(-GetVarIdOfBadGate(faulty_gate));
    solver_->add(0);
  } else {
    solver_->add(-GetVarIdOfGoodGate(faulty_gate));
    solver_->add(0);
    solver_->add(GetVarIdOfBadGate(faulty_gate));
    solver_->add(0);
  }

  if (log_info_) {
    std::cout << "Fault Observation Clauses:" << std::endl;
  }

  if (po_xors.size() > 1) {
    int final_signature = max_extra_var_id_++;
    CollectL2OR(nullptr, final_signature, po_xors);
    solver_->add(final_signature);
    solver_->add(0);
  } else {
    int final_signature = max_extra_var_id_++;
    CollectL2BUF(nullptr, final_signature, po_xors[0]);
    solver_->add(final_signature);
    solver_->add(0);
  }

  // Tied-gate clause.
  for (auto tie : prim_->GetTieGates()) {
    if (tie->GetGType() == GType::G_TIE0) {
      solver_->add(-GetVarIdOfGoodGate(tie));
      solver_->add(0);
    } else if (tie->GetGType() == GType::G_TIE1) {
      solver_->add(GetVarIdOfGoodGate(tie));
      solver_->add(0);
    } else {
      std::cout << "Unknown gate type: " << static_cast<int>(tie->GetGType())
                << std::endl;
      exit(137);
    }
  }

  solver_->gid2dchain_vars_ = d_chain_vars;
  return true;
}


void CNFGenerator::CollectGateL2(Gate* target, int out,
                                 std::vector<int> fanins) {
  switch (target->GetGType()) {
    case GType::G_AND:
      CollectL2AND(target, out, fanins);
      break;
    case GType::G_NAND:
      CollectL2NAND(target, out, fanins);
      break;
    case GType::G_OR:
      CollectL2OR(target, out, fanins);
      break;
    case GType::G_NOR:
      CollectL2NOR(target, out, fanins);
      break;
    case GType::G_MUX:
      assert(target->FaninSize() == 3);
      CollectL2MUX2(target, out, fanins[0], fanins[1], fanins[2]);
      break;
    case GType::G_XOR:
      //      assert(target->FaninSize() == 2);
      //      CollectL2XOR2(out, fanins[0], fanins[1]);
      CollectL2XOR(target, out, fanins);
      break;
    case GType::G_XNOR:
      //      assert(target->FaninSize() == 2);
      //      CollectL2XNOR2(out, fanins[0], fanins[1]);
      CollectL2XNOR(target, out, fanins);
      break;
    case GType::G_BUF:
    case GType::G_ABUF:
    case GType::G_PO:
    case GType::G_BRH:
      assert(target->FaninSize() == 1);
      CollectL2BUF(target, out, fanins[0]);
      break;
    case GType::G_NOT:
      assert(target->FaninSize() == 1);
      CollectL2INV(target, out, fanins[0]);
      break;
    case GType::G_TIE0:
    case GType::G_TIE1:
    case GType::G_PI:
      break;
    default:
      std::cout << "Unknown gate type: " << static_cast<int>(target->GetGType())
                << std::endl;
      exit(137);
  }
}

/**
 * @brief AND Gate CNF: (!z + x)(!z + y)(z + !x + !y)
 */
void CNFGenerator::CollectL2AND(Gate* target, int out,
                                std::vector<int> fanins) {
  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(fanins[i]);
    solver_->add(-out);
    solver_->add(0);
  }

  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(-fanins[i]);
  }
  solver_->add(out);
  solver_->add(0);
}

/**
 * @brief OR Gate CNF: (z + !x)(z + !y)(!z + x + y)
 */
void CNFGenerator::CollectL2OR(Gate* target, int out, std::vector<int> fanins) {
  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(-fanins[i]);
    solver_->add(out);
    solver_->add(0);
  }

  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(fanins[i]);
  }
  solver_->add(-out);
  solver_->add(0);
}

/**
 * @brief NAND Gate CNF: (z + x)(z + y)(!z + !x + !y)
 */
void CNFGenerator::CollectL2NAND(Gate* target, int out,
                                 std::vector<int> fanins) {
  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(fanins[i]);
    solver_->add(out);
    solver_->add(0);
  }

  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(-fanins[i]);
  }
  solver_->add(-out);
  solver_->add(0);
}

/**
 * @brief NOR Gate CNF: (!z + !x)(!z + !y)(z + x + y)
 */
void CNFGenerator::CollectL2NOR(Gate* target, int out,
                                std::vector<int> fanins) {
  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(-fanins[i]);
    solver_->add(-out);
    solver_->add(0);
  }

  for (int i = 0; i < fanins.size(); i++) {
    solver_->add(fanins[i]);
  }
  solver_->add(out);
  solver_->add(0);
}

/**
 * @brief MUX2 Gate CNF: (!i0 + se + y)(i0 + se + !y)(!i1 + !se + y)(i1 + !se +
 * !y)(i0 + i1 + !y)(!i0 + !i1 + y)
 */
void CNFGenerator::CollectL2MUX2(Gate* target, int out, int sel, int in0,
                                 int in1) {
  solver_->add(-in0);
  solver_->add(sel);
  solver_->add(out);
  solver_->add(0);
  solver_->add(in0);
  solver_->add(sel);
  solver_->add(-out);
  solver_->add(0);
  solver_->add(-in1);
  solver_->add(-sel);
  solver_->add(out);
  solver_->add(0);
  solver_->add(in1);
  solver_->add(-sel);
  solver_->add(-out);
  solver_->add(0);
  //  solver_->add(in0);
  //  solver_->add(in1);
  //  solver_->add(-out);
  //  solver_->add(0);
  //  solver_->add(-in0);
  //  solver_->add(-in1);
  //  solver_->add(out);
  //  solver_->add(0);
}

/**
 * @brief XOR Gate CNF: (x + y + !z)(x + !y + z)(!x + y + z)(!x + !y + !z)
 */
void CNFGenerator::CollectL2XOR(Gate* target, int out,
                                std::vector<int> fanins) {
  if (log_info_ && target) {
    std::cout << "Collect Clause at: " << target->GetGId()
              << " Gate type: " << static_cast<int>(target->GetGType())
              << " Inst Name: " << target->GetInstName() << std::endl;

    std::cout << "out: " << out << " fanin: ";
    for (auto fanin : fanins) {
      std::cout << fanin << " ";
    }
    std::cout << std::endl;
  }

  int num_in = fanins.size();
  LOG_ASSERT(pow(2, num_in) <= INT_MAX, "fan-in size overflow");

  for (int i = 0; i < pow(2, num_in); i++) {
    int truth_vals = i;
    int one_cnt = 0;
    for (int j = num_in - 1; j >= 0; j--) {
      if (truth_vals & 1) {
        one_cnt++;
        solver_->add(-fanins[j]);
      } else {
        solver_->add(fanins[j]);
      }
      truth_vals = truth_vals >> 1;
    }
    if (one_cnt & 1) {
      solver_->add(out);
    } else {
      solver_->add(-out);
    }
    solver_->add(0);
  }
}

/**
 * @brief XNOR Gate CNF: (x + y + z)(x + !y + !z)(!x + y + !z)(!x + !y + z)
 */
void CNFGenerator::CollectL2XNOR(Gate* target, int out,
                                 std::vector<int> fanins) {
  if (log_info_ && target) {
    std::cout << "Collect Clause at: " << target->GetGId()
              << " Gate type: " << static_cast<int>(target->GetGType())
              << " Inst Name: " << target->GetInstName() << std::endl;

    std::cout << "out: " << out << " fanin: ";
    for (auto fanin : fanins) {
      std::cout << fanin << " ";
    }
    std::cout << std::endl;
  }

  int num_in = fanins.size();
  LOG_ASSERT(pow(2, num_in) <= INT_MAX, "fan-in size overflow");

  for (int i = 0; i < pow(2, num_in); i++) {
    int truth_vals = i;
    int one_cnt = 0;
    for (int j = num_in - 1; j >= 0; j--) {
      if (truth_vals & 1) {
        one_cnt++;
        solver_->add(-fanins[j]);
      } else {
        solver_->add(fanins[j]);
      }
      truth_vals = truth_vals >> 1;
    }
    if (one_cnt & 1) {
      solver_->add(-out);
    } else {
      solver_->add(out);
    }
    solver_->add(0);
  }
}

/**
 * @brief XOR Gate CNF: (x + y + !z)(x + !y + z)(!x + y + z)(!x + !y + !z)
 */
void CNFGenerator::CollectL2XOR2(Gate* target, int out, int in0, int in1) {
  solver_->add(in0);
  solver_->add(in1);
  solver_->add(-out);
  solver_->add(0);
  solver_->add(in0);
  solver_->add(-in1);
  solver_->add(out);
  solver_->add(0);
  solver_->add(-in0);
  solver_->add(in1);
  solver_->add(out);
  solver_->add(0);
  solver_->add(-in0);
  solver_->add(-in1);
  solver_->add(-out);
  solver_->add(0);
}

/**
 * @brief XNOR Gate CNF: (x + y + z)(x + !y + !z)(!x + y + !z)(!x + !y + z)
 */
void CNFGenerator::CollectL2XNOR2(Gate* target, int out, int in0, int in1) {
  solver_->add(in0);
  solver_->add(in1);
  solver_->add(out);
  solver_->add(0);
  solver_->add(in0);
  solver_->add(-in1);
  solver_->add(-out);
  solver_->add(0);
  solver_->add(-in0);
  solver_->add(in1);
  solver_->add(-out);
  solver_->add(0);
  solver_->add(-in0);
  solver_->add(-in1);
  solver_->add(out);
  solver_->add(0);
}

/**
 * @brief BUF Gate CNF: (!x + y)(x + !y)
 */
void CNFGenerator::CollectL2BUF(Gate* target, int out, int in) {
  solver_->add(-in);
  solver_->add(out);
  solver_->add(0);
  solver_->add(in);
  solver_->add(-out);
  solver_->add(0);
}

/**
 * @brief INV Gate CNF: (x + y)(!x + !y)
 */
void CNFGenerator::CollectL2INV(Gate* target, int out, int in) {
  solver_->add(in);
  solver_->add(out);
  solver_->add(0);
  solver_->add(-in);
  solver_->add(-out);
  solver_->add(0);
}

int CNFGenerator::GetVarIdOfGoodGate(Gate* gate) { return gate->GetGId() + 1; }

int CNFGenerator::GetVarIdOfBadGate(Gate* gate) {
  return gate->GetGId() + 1 + prim_->GetPrimNetlist().size();
}

void CNFGenerator::CollectGateL4(Gate* target, int out, std::vector<int> fanins,
                                 std::vector<int>& x_cone) {
  switch (target->GetGType()) {
    case GType::G_AND:
      CollectL4AND(target, out, fanins, x_cone);
      break;
    case GType::G_NAND:
      CollectL4NAND(target, out, fanins, x_cone);
      //      CollectL4NAND2(target, out, fanins, x_cone);
      break;
    case GType::G_OR:
      CollectL4OR(target, out, fanins, x_cone);
      break;
    case GType::G_NOR:
      CollectL4NOR(target, out, fanins, x_cone);
      break;
    case GType::G_MUX:
      assert(target->FaninSize() == 3);
      CollectL4MUX2(target, out, fanins[0], fanins[1], fanins[2], x_cone);
      break;
    case GType::G_XOR:
      CollectL4XOR(target, out, fanins, x_cone);
      break;
    case GType::G_XNOR:
      CollectL4XNOR(target, out, fanins, x_cone);
      break;
    case GType::G_BUF:
    case GType::G_ABUF:
    case GType::G_PO:
    case GType::G_BRH:
      assert(target->FaninSize() == 1);
      CollectL4BUF(target, out, fanins[0], x_cone);
      break;
    case GType::G_NOT:
      assert(target->FaninSize() == 1);
      CollectL4INV(target, out, fanins[0], x_cone);
      break;
    case GType::G_TIE0:
    case GType::G_TIE1:
    case GType::G_TIEX:
    case GType::G_PI:
    case GType::G_DFF:
      break;
    default:
      std::cout << "Unknown gate type: " << static_cast<int>(target->GetGType())
                << std::endl;
      exit(137);
  }
}

/**
 * @brief Use L4 AND2 gate generate L4 CNF for arbitrary inputs AND gate.
 * L4 AND2 CNF: (Cy + !Cy_s)(Ca + Ca_s + !Cy)(Ca_s + Cb_s + !Cy_s)(Cb + Cb_s +
 * !Cy)(Cy + !Ca + !Cb)(Cy_s + !Ca + !Cb_s)(Cy_s + !Ca_s + !Cb)(Cy_s + !Ca_s +
 * !Cb_s)
 */
void CNFGenerator::CollectL4AND(Gate* target, int out, std::vector<int> fanins,
                                std::vector<int>& x_cone) {
  // If target is real gate, handle hybrid logic system issues.
  if (target) {
    //    for (int i = 0; i < target->FaninSize(); i++) {
    //      auto fanin = target->FaninGates()[i];
    //      if (!x_cone[fanin->GetGId()]) {
    //        solver_->add(-GetVarIdOfL4SecondBit(fanins[i]));
    //        solver_->add(0);
    //      }
    //    }
  }

  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    if (i == fanins.size() - 1) {
      new_out_bit0 = out;
      new_out_bit1 = GetVarIdOfL4SecondBit(out);
    } else {
      new_out_bit0 = max_extra_var_id_++;
      new_out_bit1 = max_extra_var_id_++;
    }

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, a_bit1, -new_out_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, b_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit0, -a_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }
}

/**
 * @brief Use L4 OR2 gate generate L4 CNF for arbitrary inputs OR gate.
 * L4 OR2 CNF: (Cy + !Ca)(Cy + !Cb)(Cy + !Cy_s)(Ca + Cy_s + !Cb_s)
 * (Cb + Cy_s + !Ca_s)(Ca_s + !Ca + !Cy_s)(Cb_s + !Cb + !Cy_s)
 * (Cy_s + !Ca_s + !Cb_s)(Ca + Ca_s + Cb + Cb_s + !Cy)
 */
void CNFGenerator::CollectL4OR(Gate* target, int out, std::vector<int> fanins,
                               std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    if (i == fanins.size() - 1) {
      new_out_bit0 = out;
      new_out_bit1 = GetVarIdOfL4SecondBit(out);
    } else {
      new_out_bit0 = max_extra_var_id_++;
      new_out_bit1 = max_extra_var_id_++;
    }

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -a_bit0});
    AddClauseToSolver({new_out_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, new_out_bit1, -b_bit1});
    AddClauseToSolver({b_bit0, new_out_bit1, -a_bit1});
    AddClauseToSolver({a_bit1, -a_bit0, -new_out_bit1});
    AddClauseToSolver({b_bit1, -b_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, a_bit1, b_bit0, b_bit1, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }
}

/**
 * @brief Use L4 AND2 gate generate L4 CNF for arbitrary inputs NAND gate.
 */
void CNFGenerator::CollectL4NAND(Gate* target, int out, std::vector<int> fanins,
                                 std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    new_out_bit0 = max_extra_var_id_++;
    new_out_bit1 = max_extra_var_id_++;

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, a_bit1, -new_out_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, b_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit0, -a_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }

  // inverter CNF.
  AddClauseToSolver({a_bit0, out});
  AddClauseToSolver({a_bit1, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({out, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({GetVarIdOfL4SecondBit(out), -a_bit1});
  AddClauseToSolver({a_bit1, -a_bit0, -out});
}

void CNFGenerator::CollectL4NAND2(Gate* target, int out,
                                  std::vector<int> fanins,
                                  std::vector<int>& x_cone) {
  //  assert(fanins.size() == 2);
  int a_bit0 = fanins[0];
  int a_bit1 = GetVarIdOfL4SecondBit(a_bit0);
  int b_bit0 = fanins[1];
  int b_bit1 = GetVarIdOfL4SecondBit(b_bit0);
  int y_bit0 = out;
  int y_bit1 = GetVarIdOfL4SecondBit(y_bit0);

  AddClauseToSolver({a_bit0, y_bit0});
  AddClauseToSolver({b_bit0, y_bit0});
  AddClauseToSolver({y_bit0, -y_bit1});
  AddClauseToSolver({-a_bit1, -b_bit1});
  AddClauseToSolver({a_bit0, a_bit1, -y_bit1});
  AddClauseToSolver({b_bit0, b_bit1, -y_bit1});
  AddClauseToSolver({y_bit1, -a_bit0, -b_bit1});
  AddClauseToSolver({y_bit1, -a_bit1, -b_bit0});
  AddClauseToSolver({a_bit1, b_bit1, -a_bit0, -b_bit0, -y_bit0});
}

/**
 * @brief Use L4 OR2 gate generate L4 CNF for arbitrary inputs NOR gate.
 */
void CNFGenerator::CollectL4NOR(Gate* target, int out, std::vector<int> fanins,
                                std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    new_out_bit0 = max_extra_var_id_++;
    new_out_bit1 = max_extra_var_id_++;

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -a_bit0});
    AddClauseToSolver({new_out_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, new_out_bit1, -b_bit1});
    AddClauseToSolver({b_bit0, new_out_bit1, -a_bit1});
    AddClauseToSolver({a_bit1, -a_bit0, -new_out_bit1});
    AddClauseToSolver({b_bit1, -b_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, a_bit1, b_bit0, b_bit1, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }

  // inverter CNF.
  AddClauseToSolver({a_bit0, out});
  AddClauseToSolver({a_bit1, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({out, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({GetVarIdOfL4SecondBit(out), -a_bit1});
  AddClauseToSolver({a_bit1, -a_bit0, -out});
}

/**
 * @brief L4 MUX2 CNF: (Cs + !Cs_s)(cy + !Cy_s)(!Ca_s + !Cs_s)(!Cb_s + !Cs_s)(Cs
 * + Cy + !Ca)(Cs + Cy_s + !Ca_s)(Cs + !Ca + !Cy_s)(Cy + !Cb + !Cs)(Cy_s + !Cb_s
 * + !Cs)(Ca + Ca_s + Cs + !Cy)(!Ca + !Cb + !Cy_s)(Ca + Cy_s + !Cs_s + !Cy)(Cb +
 * Cy_s + !Ca + !Cs_s)(Cs_s + !Cb + !Cs + !Cy_s)(Ca + Cb + Cb_s + !Cs +
 * !Cy_s)(Cb + Cb_s + Cs_s + !Cs + !Cy)
 */
void CNFGenerator::CollectL4MUX2(Gate* target, int out, int sel, int in0,
                                 int in1, std::vector<int>& x_cone) {
  int se_bit0 = sel;
  int se_bit1 = GetVarIdOfL4SecondBit(sel);
  int a_bit0 = in0;
  int a_bit1 = GetVarIdOfL4SecondBit(in0);
  int b_bit0 = in1;
  int b_bit1 = GetVarIdOfL4SecondBit(in1);

  int out_bit0 = out;
  int out_bit1 = GetVarIdOfL4SecondBit(out);

  AddClauseToSolver({se_bit0, -se_bit1});
  AddClauseToSolver({out_bit0, -out_bit1});
  AddClauseToSolver({-a_bit1, -se_bit1});
  AddClauseToSolver({-b_bit1, -se_bit1});
  AddClauseToSolver({se_bit0, out_bit0, -a_bit0});
  AddClauseToSolver({se_bit0, out_bit1, -a_bit1});
  AddClauseToSolver({se_bit0, -a_bit0, -out_bit1});
  AddClauseToSolver({out_bit0, -b_bit0, -se_bit0});
  AddClauseToSolver({out_bit1, -b_bit1, -se_bit0});
  AddClauseToSolver({a_bit0, a_bit1, se_bit0, -out_bit0});
  AddClauseToSolver({-a_bit0, -b_bit0, -out_bit1});
  AddClauseToSolver({a_bit0, out_bit1, -se_bit1, -out_bit0});
  AddClauseToSolver({b_bit0, out_bit1, -a_bit0, -se_bit1});
  AddClauseToSolver({se_bit1, -b_bit0, -se_bit0, -out_bit1});
  AddClauseToSolver({a_bit0, b_bit0, b_bit1, -se_bit0, -out_bit1});
  AddClauseToSolver({b_bit0, b_bit1, se_bit1, -se_bit0, -out_bit0});
}

/**
 * @brief Use L4 XOR2 gate generate L4 CNF for arbitrary inputs XOR gate.
 * L4 XOR2 CNF: (Cy + !Cy_s)(Cy_s + !Ca_s)(Cy_s + !Cb_s)(Ca + Cy + !Cb)(Ca_s +
 * Cb_s + !Cy_s)(Cb + Cy + !Ca)(Ca + Cb + Cy_s + !Cy)(Cy_s + !Ca + !Cb + !Cy)
 */
void CNFGenerator::CollectL4XOR(Gate* target, int out, std::vector<int> fanins,
                                std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    if (i == fanins.size() - 1) {
      new_out_bit0 = out;
      new_out_bit1 = GetVarIdOfL4SecondBit(out);
    } else {
      new_out_bit0 = max_extra_var_id_++;
      new_out_bit1 = max_extra_var_id_++;
    }

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1});
    AddClauseToSolver({new_out_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, new_out_bit0, -b_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, new_out_bit0, -a_bit0});
    AddClauseToSolver({a_bit0, b_bit0, new_out_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit0, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }
}

/**
 * @brief Use L4 XOR2 gate generate L4 CNF for arbitrary inputs XNOR gate.
 */
void CNFGenerator::CollectL4XNOR(Gate* target, int out, std::vector<int> fanins,
                                 std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    new_out_bit0 = max_extra_var_id_++;
    new_out_bit1 = max_extra_var_id_++;

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1});
    AddClauseToSolver({new_out_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, new_out_bit0, -b_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, new_out_bit0, -a_bit0});
    AddClauseToSolver({a_bit0, b_bit0, new_out_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit0, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }

  // inverter CNF.
  AddClauseToSolver({a_bit0, out});
  AddClauseToSolver({a_bit1, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({out, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({GetVarIdOfL4SecondBit(out), -a_bit1});
  AddClauseToSolver({a_bit1, -a_bit0, -out});
}

/**
 * @brief L4 BUF CNF: (Ci_s + !Co_s)(Co + !Ci)(Co + !Co_s)(Co_s + !Ci_s)(Ci +
 * Ci_s + !Co)
 */
void CNFGenerator::CollectL4BUF(Gate* target, int out, int in,
                                std::vector<int>& x_cone) {
  //  if (target) {
  //    assert(x_cone[target->FaninGates()[0]->GetGId()]);
  //  }

  int in_bit0 = in;
  int in_bit1 = GetVarIdOfL4SecondBit(in);
  int out_bit0 = out;
  int out_bit1 = GetVarIdOfL4SecondBit(out);

  // buffer CNF.
  AddClauseToSolver({in_bit1, -out_bit1});
  AddClauseToSolver({out_bit0, -in_bit0});
  AddClauseToSolver({out_bit0, -out_bit1});
  AddClauseToSolver({out_bit1, -in_bit1});
  AddClauseToSolver({in_bit0, in_bit1, -out_bit0});
}

/**
 * @brief L4 INV CNF: (Ci + Co)(Ci_s + !Co_s)(Co + !Co_s)(Co_s + !Ci_s)(Ci_s +
 * !Ci + Co)
 */
void CNFGenerator::CollectL4INV(Gate* target, int out, int in,
                                std::vector<int>& x_cone) {
  //  if (target) {
  //    assert(x_cone[target->FaninGates()[0]->GetGId()]);
  //  }

  int in_bit0 = in;
  int in_bit1 = GetVarIdOfL4SecondBit(in);
  int out_bit0 = out;
  int out_bit1 = GetVarIdOfL4SecondBit(out);

  // inverter CNF.
  AddClauseToSolver({in_bit0, out_bit0});
  AddClauseToSolver({in_bit1, -out_bit1});
  AddClauseToSolver({out_bit0, -out_bit1});
  AddClauseToSolver({out_bit1, -in_bit1});
  AddClauseToSolver({in_bit1, -in_bit0, -out_bit0});
}

void CNFGenerator::CollectGateL3(Gate* target, int out, std::vector<int> fanins,
                                 std::vector<int>& x_cone) {
  switch (target->GetGType()) {
    case GType::G_AND:
      CollectL3AND(target, out, fanins, x_cone);
      break;
    case GType::G_NAND:
      CollectL3NAND(target, out, fanins, x_cone);
      break;
    case GType::G_OR:
      CollectL3OR(target, out, fanins, x_cone);
      break;
    case GType::G_NOR:
      CollectL3NOR(target, out, fanins, x_cone);
      break;
    case GType::G_MUX:
      assert(target->FaninSize() == 3);
      CollectL3MUX2(target, out, fanins[0], fanins[1], fanins[2], x_cone);
      break;
    case GType::G_XOR:
      CollectL3XOR(target, out, fanins, x_cone);
      break;
    case GType::G_XNOR:
      CollectL3XNOR(target, out, fanins, x_cone);
      break;
    case GType::G_BUF:
    case GType::G_ABUF:
    case GType::G_PO:
    case GType::G_BRH:
      assert(target->FaninSize() == 1);
      CollectL3BUF(target, out, fanins[0], x_cone);
      break;
    case GType::G_NOT:
      assert(target->FaninSize() == 1);
      CollectL3INV(target, out, fanins[0], x_cone);
      break;
    case GType::G_TIE0:
    case GType::G_TIE1:
    case GType::G_TIEX:
    case GType::G_PI:
    case GType::G_DFF:
      break;
    default:
      std::cout << "Unknown gate type: " << static_cast<int>(target->GetGType())
                << std::endl;
      exit(137);
  }
}

/**
 * @brief Use L3 AND2 gate generate L3 CNF for arbitrary inputs AND gate.
 * L3 AND2 CNF: (Cy + !Cy_s)(Ca + Ca_s + !Cy)(Ca_s + Cb_s + !Cy_s)(Cb + Cb_s +
 * !Cy)(Cy + !Ca + !Cb)(Cy_s + !Ca + !Cb_s)(Cy_s + !Ca_s + !Cb)(Cy_s + !Ca_s +
 * !Cb_s)
 */
void CNFGenerator::CollectL3AND(Gate* target, int out, std::vector<int> fanins,
                                std::vector<int>& x_cone) {
  // If target is real gate, handle hybrid logic system issues.
  if (target) {
    //    for (int i = 0; i < target->FaninSize(); i++) {
    //      auto fanin = target->FaninGates()[i];
    //      if (!x_cone[fanin->GetGId()]) {
    //        solver_->add(-GetVarIdOfL4SecondBit(fanins[i]));
    //        solver_->add(0);
    //      }
    //    }
  }

  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    if (i == fanins.size() - 1) {
      new_out_bit0 = out;
      new_out_bit1 = GetVarIdOfL4SecondBit(out);
    } else {
      new_out_bit0 = max_extra_var_id_++;
      new_out_bit1 = max_extra_var_id_++;
    }

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, a_bit1, -new_out_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, b_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit0, -a_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }
}

/**
 * @brief Use L3 OR2 gate generate L3 CNF for arbitrary inputs OR gate.
 * L3 OR2 CNF: (Cy + !Ca)(Cy + !Cb)(Cy + !Cy_s)(Ca + Cy_s + !Cb_s)
 * (Cb + Cy_s + !Ca_s)(Ca_s + !Ca + !Cy_s)(Cb_s + !Cb + !Cy_s)
 * (Cy_s + !Ca_s + !Cb_s)(Ca + Ca_s + Cb + Cb_s + !Cy)
 */
void CNFGenerator::CollectL3OR(Gate* target, int out, std::vector<int> fanins,
                               std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    if (i == fanins.size() - 1) {
      new_out_bit0 = out;
      new_out_bit1 = GetVarIdOfL4SecondBit(out);
    } else {
      new_out_bit0 = max_extra_var_id_++;
      new_out_bit1 = max_extra_var_id_++;
    }

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -a_bit0});
    AddClauseToSolver({new_out_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, new_out_bit1, -b_bit1});
    AddClauseToSolver({b_bit0, new_out_bit1, -a_bit1});
    AddClauseToSolver({a_bit1, -a_bit0, -new_out_bit1});
    AddClauseToSolver({b_bit1, -b_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, a_bit1, b_bit0, b_bit1, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }
}

/**
 * @brief Use L3 AND2 gate generate L3 CNF for arbitrary inputs NAND gate.
 */
void CNFGenerator::CollectL3NAND(Gate* target, int out, std::vector<int> fanins,
                                 std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    new_out_bit0 = max_extra_var_id_++;
    new_out_bit1 = max_extra_var_id_++;

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, a_bit1, -new_out_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, b_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit0, -a_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }

  // inverter CNF.
  AddClauseToSolver({a_bit0, out});
  AddClauseToSolver({a_bit1, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({out, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({GetVarIdOfL4SecondBit(out), -a_bit1});
  AddClauseToSolver({a_bit1, -a_bit0, -out});
}

/**
 * @brief Use L3 OR2 gate generate L3 CNF for arbitrary inputs NOR gate.
 */
void CNFGenerator::CollectL3NOR(Gate* target, int out, std::vector<int> fanins,
                                std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    new_out_bit0 = max_extra_var_id_++;
    new_out_bit1 = max_extra_var_id_++;

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -a_bit0});
    AddClauseToSolver({new_out_bit0, -b_bit0});
    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({a_bit0, new_out_bit1, -b_bit1});
    AddClauseToSolver({b_bit0, new_out_bit1, -a_bit1});
    AddClauseToSolver({a_bit1, -a_bit0, -new_out_bit1});
    AddClauseToSolver({b_bit1, -b_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, a_bit1, b_bit0, b_bit1, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }

  // inverter CNF.
  AddClauseToSolver({a_bit0, out});
  AddClauseToSolver({a_bit1, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({out, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({GetVarIdOfL4SecondBit(out), -a_bit1});
  AddClauseToSolver({a_bit1, -a_bit0, -out});
}

/**
 * @brief L3 MUX2 CNF: (Cy + !Cy_s)(!Ca_s + !Cs_s)(!Cb_s + !Cs_s)(Cs + Cy +
 * !Ca)(Cs + Cy_s + !Ca_s)(Cy + !Cb + !Cs)(Cy_s + !Cb_s + !Cs)(Ca + Cs + Cy_s +
 * !Cy)(Ca_s + Cs + Cs_s + !Cy_s)(Ca + Cb + !Cs_s + !Cy_s)(Ca + Cy_s + !Cb +
 * !Cs_s)(Cb + Cy_s + !Ca + !Cs_s)(Cb + Cy_s + !Cs + !Cy)(Cb_s + Cs_s + !Cs +
 * !Cy_s)(!Ca + !Cb + !Cs_s + !Cy_s)
 */
void CNFGenerator::CollectL3MUX2(Gate* target, int out, int sel, int in0,
                                 int in1, std::vector<int>& x_cone) {
  int se_bit0 = sel;
  int se_bit1 = GetVarIdOfL4SecondBit(sel);
  int a_bit0 = in0;
  int a_bit1 = GetVarIdOfL4SecondBit(in0);
  int b_bit0 = in1;
  int b_bit1 = GetVarIdOfL4SecondBit(in1);

  int out_bit0 = out;
  int out_bit1 = GetVarIdOfL4SecondBit(out);

  AddClauseToSolver({out_bit0, -out_bit1});
  AddClauseToSolver({-a_bit1, -se_bit1});
  AddClauseToSolver({-b_bit1, -se_bit1});
  AddClauseToSolver({se_bit0, out_bit0, -a_bit0});
  AddClauseToSolver({se_bit0, out_bit1, -a_bit1});
  AddClauseToSolver({out_bit0, -b_bit0, -se_bit0});
  AddClauseToSolver({out_bit1, -b_bit1, -se_bit0});
  AddClauseToSolver({a_bit0, se_bit0, out_bit1, -out_bit0});
  AddClauseToSolver({a_bit1, se_bit0, se_bit1, -out_bit1});
  AddClauseToSolver({a_bit0, b_bit0, -se_bit1, -out_bit1});
  AddClauseToSolver({a_bit0, out_bit1, -b_bit0, -se_bit1});
  AddClauseToSolver({b_bit0, out_bit1, -a_bit0, -se_bit1});
  AddClauseToSolver({b_bit0, out_bit1, -se_bit0, -out_bit0});
  AddClauseToSolver({b_bit1, se_bit1, -se_bit0, -out_bit1});
  AddClauseToSolver({-a_bit0, -b_bit0, -se_bit1, -out_bit1});
}

/**
 * @brief Use L3 XOR2 gate generate L3 CNF for arbitrary inputs XOR gate.
 * L3 XOR2 CNF: (Cy + !Cy_s)(Cy_s + !Ca_s)(Cy_s + !Cb_s)(Ca + Cy + !Cb)(Ca_s +
 * Cb_s + !Cy_s)(Cb + Cy + !Ca)(Ca + Cb + Cy_s + !Cy)(Cy_s + !Ca + !Cb + !Cy)
 */
void CNFGenerator::CollectL3XOR(Gate* target, int out, std::vector<int> fanins,
                                std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    if (i == fanins.size() - 1) {
      new_out_bit0 = out;
      new_out_bit1 = GetVarIdOfL4SecondBit(out);
    } else {
      new_out_bit0 = max_extra_var_id_++;
      new_out_bit1 = max_extra_var_id_++;
    }

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1});
    AddClauseToSolver({new_out_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, new_out_bit0, -b_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, new_out_bit0, -a_bit0});
    AddClauseToSolver({a_bit0, b_bit0, new_out_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit0, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }
}

/**
 * @brief Use L3 XOR2 gate generate L3 CNF for arbitrary inputs XNOR gate.
 */
void CNFGenerator::CollectL3XNOR(Gate* target, int out, std::vector<int> fanins,
                                 std::vector<int>& x_cone) {
  int a_bit0 = fanins[0], a_bit1 = GetVarIdOfL4SecondBit(fanins[0]);
  int b_bit0, b_bit1;
  int new_out_bit0, new_out_bit1;

  for (int i = 1; i < fanins.size(); i++) {
    new_out_bit0 = max_extra_var_id_++;
    new_out_bit1 = max_extra_var_id_++;

    b_bit0 = fanins[i];
    b_bit1 = GetVarIdOfL4SecondBit(fanins[i]);

    AddClauseToSolver({new_out_bit0, -new_out_bit1});
    AddClauseToSolver({new_out_bit1, -a_bit1});
    AddClauseToSolver({new_out_bit1, -b_bit1});
    AddClauseToSolver({a_bit0, new_out_bit0, -b_bit0});
    AddClauseToSolver({a_bit1, b_bit1, -new_out_bit1});
    AddClauseToSolver({b_bit0, new_out_bit0, -a_bit0});
    AddClauseToSolver({a_bit0, b_bit0, new_out_bit1, -new_out_bit0});
    AddClauseToSolver({new_out_bit1, -a_bit0, -b_bit0, -new_out_bit0});

    a_bit0 = new_out_bit0;
    a_bit1 = new_out_bit1;
  }

  // inverter CNF.
  AddClauseToSolver({a_bit0, out});
  AddClauseToSolver({a_bit1, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({out, -GetVarIdOfL4SecondBit(out)});
  AddClauseToSolver({GetVarIdOfL4SecondBit(out), -a_bit1});
  AddClauseToSolver({a_bit1, -a_bit0, -out});
}

/**
 * @brief L3 BUF CNF: (Ci_s + !Co_s)(Co + !Ci)(Co + !Co_s)(Co_s + !Ci_s)(Ci +
 * Ci_s + !Co)
 */
void CNFGenerator::CollectL3BUF(Gate* target, int out, int in,
                                std::vector<int>& x_cone) {
  //  if (target) {
  //    assert(x_cone[target->FaninGates()[0]->GetGId()]);
  //  }

  int in_bit0 = in;
  int in_bit1 = GetVarIdOfL4SecondBit(in);
  int out_bit0 = out;
  int out_bit1 = GetVarIdOfL4SecondBit(out);

  // buffer CNF.
  AddClauseToSolver({in_bit1, -out_bit1});
  AddClauseToSolver({out_bit0, -in_bit0});
  AddClauseToSolver({out_bit0, -out_bit1});
  AddClauseToSolver({out_bit1, -in_bit1});
  AddClauseToSolver({in_bit0, in_bit1, -out_bit0});
}

/**
 * @brief L3 INV CNF: (Ci + Co)(Ci_s + !Co_s)(Co + !Co_s)(Co_s + !Ci_s)(Ci_s +
 * !Ci + Co)
 */
void CNFGenerator::CollectL3INV(Gate* target, int out, int in,
                                std::vector<int>& x_cone) {
  //  if (target) {
  //    assert(x_cone[target->FaninGates()[0]->GetGId()]);
  //  }

  int in_bit0 = in;
  int in_bit1 = GetVarIdOfL4SecondBit(in);
  int out_bit0 = out;
  int out_bit1 = GetVarIdOfL4SecondBit(out);

  // inverter CNF.
  AddClauseToSolver({in_bit0, out_bit0});
  AddClauseToSolver({in_bit1, -out_bit1});
  AddClauseToSolver({out_bit0, -out_bit1});
  AddClauseToSolver({out_bit1, -in_bit1});
  AddClauseToSolver({in_bit1, -in_bit0, -out_bit0});
}

/**
 * @brief Get Bit1 in L4 logic system of real gate variable.
 */
int CNFGenerator::GetVarIdOfL4SecondBit(int var_id) {
  return var_id + 2 * prim_->GetPrimNetlist().size();
}

int CNFGenerator::EvalGate(Gate* gate) {
  int value;
  int size = gate->FaninSize();
  std::vector<int>& values = model_;

  switch (gate->GetGType()) {
    case GType::G_PI: {
      value = values[gate->GetGId()];
    } break;
    case GType::G_TIE0:
      value = LOGIC_0;
      break;
    case GType::G_TIE1:
      value = LOGIC_1;
      break;
    case GType::G_BUF:
    case GType::G_PO:
    case GType::G_ABUF:
    case GType::G_BRH:
      value = buf_table_nine[values[gate->FaninGates()[0]->GetGId()]];
      break;
    case GType::G_NOT:
      value = not_table_nine[values[gate->FaninGates()[0]->GetGId()]];
      break;
    case GType::G_AND:
      value = values[gate->FaninGates()[0]->GetGId()];

      for (int i = 1; i < size; i++) {
        value = and_table_nine[value][values[gate->FaninGates()[i]->GetGId()]];
      }
      break;
    case GType::G_NAND:
      value = values[gate->FaninGates()[0]->GetGId()];

      for (int i = 1; i < size; i++) {
        value = and_table_nine[value][values[gate->FaninGates()[i]->GetGId()]];
      }
      value = not_table_nine[value];
      break;
    case GType::G_OR:
      value = values[gate->FaninGates()[0]->GetGId()];

      for (int i = 1; i < size; i++) {
        value = or_table_nine[value][values[gate->FaninGates()[i]->GetGId()]];
      }
      break;
    case GType::G_NOR:
      value = values[gate->FaninGates()[0]->GetGId()];

      for (int i = 1; i < size; i++) {
        value = or_table_nine[value][values[gate->FaninGates()[i]->GetGId()]];
      }
      value = not_table_nine[value];
      break;
    case GType::G_XOR:
      value = values[gate->FaninGates()[0]->GetGId()];

      for (int i = 1; i < size; i++) {
        value = xor_table_nine[value][values[gate->FaninGates()[i]->GetGId()]];
      }
      break;
    case GType::G_XNOR:
      value = values[gate->FaninGates()[0]->GetGId()];

      for (int i = 1; i < size; i++) {
        value = xor_table_nine[value][values[gate->FaninGates()[i]->GetGId()]];
      }
      value = not_table_nine[value];
      break;
    case GType::G_MUX: {
      // 0: se, 1: in0, 2: in1.
      uint8_t in0;
      in0 = values[gate->FaninGates()[1]->GetGId()];
      uint8_t in1;
      in1 = values[gate->FaninGates()[2]->GetGId()];
      uint8_t se;
      se = values[gate->FaninGates()[0]->GetGId()];
      value = or_table_nine[and_table_nine[in0][not_table_nine[se]]]
                           [and_table_nine[in1][se]];

      // se is a reconvergent stem, so X can't hold all conditions.
      if (se == LOGIC_x && in0 == in1) {
        value = in0;
      }
    } break;
    default:
      std::cout << "ERROR: this gate type is not supported!" << std::endl;
      std::cout << static_cast<int>(gate->GetGType()) << std::endl;
      assert(false);
  }

  return value;
}

void CNFGenerator::WriteCNF(vector<int>& var_old2new, vector<int>& var_new2old) {
    Minisat::SATSolver S;
    std::ofstream outFile("../../compare_sat/cnf.txt", std::ios::out | std::ios::trunc); // Open in append mode
    if (outFile) {
        // add variables and clauses
        vector<int> old_clause;
        std::string new_clauses;
        int nVar = 0;
        int nClause = 0;
        // record a map from old2new and new2old
        var_old2new = vector<int>(solver_->max_var_size_, -1);
        var_new2old = vector<int>(solver_->max_var_size_, -1);
        for (int i = 0; i < solver_->cnf_var_intv_.size(); i++) {
            if (solver_->cnf_var_intv_[i] == 0) {
                nClause++;
                // get a single clause
                std::string tmp_v;
                for (int j = 0; j < old_clause.size(); ++j) {
                    int abs = 0;
                    bool sign = false;
                    if (old_clause[j] >= 0) {
                        abs = old_clause[j];
                    } else {
                        abs = ::abs(old_clause[j]);
                        sign = true;
                    }
                    if (var_old2new[abs] != -1) {
                        Var new_var = var_old2new[abs];
                        if (j != 0) {
                            if (sign) {
                                tmp_v += (" -" + std::to_string(new_var + 1));
                            } else {
                                tmp_v += (" " +std::to_string(new_var + 1));
                            }
                        } else {
                            if (sign) {
                                tmp_v += ("-" + std::to_string(new_var + 1));
                            } else {
                                tmp_v += std::to_string(new_var + 1);
                            }
                        }
                    } else {
                        Var new_var = S.newVar();
                        nVar++;
                        var_new2old[new_var] = abs;
                        var_old2new[abs] = new_var;
                        if (j != 0) {
                            if (sign) {
                                tmp_v += (" -" + std::to_string(new_var + 1));
                            } else {
                                tmp_v += (" " +std::to_string(new_var + 1));
                            }
                        } else {
                            if (sign) {
                                tmp_v += ("-" + std::to_string(new_var + 1));
                            } else {
                                tmp_v += std::to_string(new_var + 1);
                            }
                        }
                    }
                }
                old_clause.clear();
                new_clauses += (tmp_v + " 0\n");
                continue;
            } else {
                old_clause.push_back(solver_->cnf_var_intv_[i]);
            }
        }
        new_clauses = ("p cnf " + std::to_string(nVar) + " " + std::to_string(nClause) + "\n") + new_clauses;
        outFile << new_clauses << std::endl;
        //std::cout << "Data written to file." << std::endl;
    } else {
        std::cerr << "Unable to open file for writing." << std::endl;
    }
    outFile.close(); // Close the file

}

int CNFGenerator::ReadSATResultMini(vector<bool>& old_vars_signs, vector<int>& var_new2old) {
    std::ifstream inFile("../../compare_sat/output.txt");
    std::string data;
    int re_val = 0;

    if (inFile) {
        std::getline(inFile, data);
        if (data == "SAT") {
            re_val = 10;
        } else if (data == "UNSAT") {
            re_val = 20;
        }
        std::getline(inFile, data);
        std::istringstream vars(data);
        int var;
        while (vars >> var) { // Extract numbers from the stream
            if (var != 0) {
                int var_abs = std::abs(var) - 1;
                if (var_new2old[var_abs] != -1) {
                    if (var > 0) {
                        old_vars_signs[var_new2old[var_abs]] = true;
                    } else {
                        old_vars_signs[var_new2old[var_abs]] = false;
                    }

                }
            }
        }
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
    inFile.close(); // Close the file
    return re_val;
}

int CNFGenerator::ReadSATResultKis(vector<bool>& old_vars_signs, vector<int>& var_new2old) {
    std::ifstream inFile("../../compare_sat/output.txt");
    std::string data;
    int re_val = 0;

    if (inFile) {
        while (std::getline(inFile, data)) {
            if (data == "c ---- [ result ] ------------------------------------------------------------") {
                break;
            }
        }
        std::getline(inFile, data);
        std::getline(inFile, data);

        if (data == "s SATISFIABLE") {
            re_val = 10;
        } else if (data == "s UNSATISFIABLE") {
            re_val = 20;
        }
        if (re_val == 10) {
            while (std::getline(inFile, data)) {
                if (data[0] != 'v') {
                    break;
                }
                std::istringstream vars(data);
                std::string var;
                while (vars >> var) { // Extract numbers from the stream
                    if (var != "0" && var != "v") {
                        int var_abs = std::abs(stoi(var)) - 1;
                        if (var_new2old[var_abs] != -1) {
                            if (stoi(var) > 0) {
                                old_vars_signs[var_new2old[var_abs]] = true;
                            } else {
                                old_vars_signs[var_new2old[var_abs]] = false;
                            }

                        }
                    }
                }
            }

        }
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
    inFile.close(); // Close the file
    return re_val;
}

int CNFGenerator::ReadSATResultCal(vector<bool>& old_vars_signs, vector<int>& var_new2old) {
    std::ifstream inFile("../../compare_sat/output.txt");
    std::string data;
    int re_val = 0;

    if (inFile) {
        while (std::getline(inFile, data)) {
            if (data == "c --- [ result ] -------------------------------------------------------------") {
                break;
            }
        }
        std::getline(inFile, data);
        std::getline(inFile, data);

        if (data == "s SATISFIABLE") {
            re_val = 10;
        } else if (data == "s UNSATISFIABLE") {
            re_val = 20;
        }
        if (re_val == 10) {
            while (std::getline(inFile, data)) {
                if (data[0] != 'v') {
                    break;
                }
                std::istringstream vars(data);
                std::string var;
                while (vars >> var) { // Extract numbers from the stream
                    if (var != "0" && var != "v") {
                        int var_abs = std::abs(stoi(var)) - 1;
                        if (var_new2old[var_abs] != -1) {
                            if (stoi(var) > 0) {
                                old_vars_signs[var_new2old[var_abs]] = true;
                            } else {
                                old_vars_signs[var_new2old[var_abs]] = false;
                            }

                        }
                    }
                }
            }

        }
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
    inFile.close(); // Close the file
    return re_val;
}

int CNFGenerator::ReadSATResultPotopara(vector<bool>& old_vars_signs, vector<int>& var_new2old) {
    std::ifstream inFile("../../compare_sat/output.txt");
    std::string data;
    int re_val = 0;

    if (inFile) {
        std::getline(inFile, data);
        if (data == "s SATISFIABLE") {
            re_val = 10;
        } else if (data == "s UNSATISFIABLE") {
            re_val = 20;
        }
        std::getline(inFile, data);
        if (data[0] != 'v') {
            return re_val;
        }
        std::istringstream vars(data);
        std::string var_s;
        while (vars >> var_s) { // Extract numbers from the stream
            if (var_s != "0" && var_s != "v") {
                int var = std::stoi(var_s);
                int var_abs = std::abs(var) - 1;
                if (var_new2old[var_abs] != -1) {
                    if (var > 0) {
                        old_vars_signs[var_new2old[var_abs]] = true;
                    } else {
                        old_vars_signs[var_new2old[var_abs]] = false;
                    }

                }
            }
        }
    } else {
        std::cerr << "Unable to open file for reading." << std::endl;
    }
    inFile.close(); // Close the file
    return re_val;
}

}  // namespace ictest
