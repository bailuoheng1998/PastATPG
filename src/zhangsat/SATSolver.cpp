/***************************************************************************************[Solver.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <math.h>

#include "zhangsat/Alg.h"
#include "zhangsat/Sort.h"
#include "zhangsat/System.h"
#include "zhangsat/SATSolver.h"
#include <vector>
using std::vector;

using namespace Minisat;

//=================================================================================================
// Options:
system_clock::time_point totoal_time_begin;


double SAT_time = 0.0;      // SAT 调用次数
double UNSAT_time = 0.0;    // 时间
unsigned SAT_ct = 0;        // UNSAT 调用次数
unsigned UNSAT_ct = 0;      // 时间

long long pre_cal_ct = 0;     // IC3 pre层面：计算PA的次数
long long pre_full_sz = 0;    // 状态涉及的变量数目
long long pre_pa_sz = 0;      // PA状态变量数目
long long pre_final_sz = 0;   // 最终core化简变量数目
double pre_cal_time = 0.0;      // 耗时

long long pa_cal_ct = 0;      // SAT层面： 计算PA的次数
long long pa_full_sz = 0;     // 完全赋值变量数目
long long pa_final_sz = 0;    // PA变量数目
double pa_cal_time = 0.0;       // 耗时
const double cutoff_time = 3600.0; 

double get_time_since_begin(){
    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - totoal_time_begin);
    double time_in_sec = double(duration.count()) * microseconds::period::num / microseconds::period::den;
    return time_in_sec;
}

void pa_show_statics(){
    printf("c info SAT_ct = %u SAT_time = %f \n", SAT_ct, SAT_time);
    printf("c info UNSAT_ct = %u UNSAT_time = %f \n", UNSAT_ct, UNSAT_time);
    printf("c info pre_cal_ct = %lld pre_full_sz = %lld pre_pa_sz = %lld pre_final_sz = %lld pre_cal_time = %f \n",
        pre_cal_ct, pre_full_sz, pre_pa_sz, pre_final_sz, pre_cal_time);
    printf("c info pa_cal_ct = %lld pa_full_sz = %lld pa_final_sz = %lld pa_cal_time = %f \n",
        pa_cal_ct, pa_full_sz, pa_final_sz, pa_cal_time);
    printf("c info ALL_time = %f \n", get_time_since_begin());
    fflush(stdout);
}


static const char* _cat = "CORE";

static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.95,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_clause_decay      (_cat, "cla-decay",   "The clause activity decay factor",              0.999,    DoubleRange(0, false, 1, false));
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static DoubleOption  opt_random_seed       (_cat, "rnd-seed",    "Used by the random variable selection",         91648253, DoubleRange(0, false, HUGE_VAL, false));
static IntOption     opt_ccmin_mode        (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);
static BoolOption    opt_luby_restart      (_cat, "luby",        "Use the Luby restart sequence", true);
static IntOption     opt_restart_first     (_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption  opt_restart_inc       (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));
static DoubleOption  opt_garbage_frac      (_cat, "gc-frac",     "The fraction of wasted memory allowed before a garbage collection is triggered",  0.20, DoubleRange(0, false, HUGE_VAL, false));
static IntOption     opt_min_learnts_lim   (_cat, "min-learnts", "Minimum learnt clause limit",  0, IntRange(0, INT32_MAX));


//=================================================================================================
// Constructor/Destructor:


int to_dimacs(Minisat::Lit l){
    int v = l.x/2 + 1;
    if(l.x%2 == 1)
        return -v;
    else
        return v;
}

char to_char(lbool lb){
    if(lb == l_False){
        return 'F';
    }else if(lb == l_True){
        return 'T';
    }else if(lb == l_Undef){
        return 'X';
    }
    return '-';
}


// void SATSolver::check_sat_unsat_ct(){
//     for(int i=0; i<2*nVars(); ++i){
//         vector<CRef> &ws = full_watches[i];
//         for(unsigned j=0; j<ws.size(); ++j){
//             Clause &c = ca[ws[j]];
//             unsigned sat_ct = 0;
//             unsigned unsat_ct = 0;
//             for(int k=0;k<c.size();++k){
//                 if(value(c[k])==l_True)
//                     sat_ct += 1;
//                 else if(value(c[k])==l_False)
//                     unsat_ct += 1;
//             }
//             if(sat_ct != c.sat_ct()){
//                 vis_clause(ws[j], false);
//                 printf("\nerror at level %d \n", decisionLevel());
//                 vis_trail();
//                 exit(0);
//             }
//             if(unsat_ct != c.unsat_ct()){
//                 vis_clause(ws[j], false);
//                 printf("\nerror at level %d \n", decisionLevel());
//                 vis_trail();
//                 exit(0);
//             }
//         }
//     }
// }

void SATSolver::vis_clause(CRef cref, bool check = false){
    Clause &c = ca[cref];
    printf("<%u>[sz%d]{s%d}(u%d)%c| ",(uint32_t)cref, c.size(),c.sat_ct(),c.unsat_ct(),(c.learnt()?'L':'O'));
    unsigned sat_ct = 0;
    unsigned unsat_ct = 0;
    for(int j=0; j<c.size(); ++j){
        printf("%d(%c) ",to_dimacs(c[j]), to_char(value(c[j])));
        if(check){
            if(value(c[j]) == l_True)
                sat_ct ++;
            if(value(c[j]) == l_False)
                unsat_ct ++;
        }
    }
    if(check){
        assert(sat_ct == c.sat_ct());
        assert(unsat_ct == c.unsat_ct());
    }
}

void SATSolver::vis_full_watches(int watch = -1, bool check=true){
    for (int v = 0; v < nVars(); v++){
        if(watch == -1 || watch-1 == v){
            for (int s = 0; s < 2; s++){
                Lit p = mkLit(v, s);
                printf("%d: \t", to_dimacs(p));
                vec<CRef> &ws = full_watches[p]; 
                for(int i=0; i<ws.size(); ++i){
                    vis_clause(ws[i], check);
                    printf("; \n \t");
                }
                printf("\n");
            }
        }
    }
}


void SATSolver::vis_partial_watches(int watch = -1, bool check=true){
    for (int v = 0; v < nVars(); v++){
        if(watch == -1 || watch-1 == v){
            for (int s = 0; s < 2; s++){
                Lit p = mkLit(v, s);
                printf("%d: \t", to_dimacs(p));
                vec<Watcher> &ws = partial_watches[p]; 
                for(int i=0; i<ws.size(); ++i){
                    vis_clause(ws[i].cref, check);
                    printf("; \n \t");
                }
                printf("\n");
            }
        }
    }
}


void SATSolver::vis_trail(){
    printf("[%d]: \t", qhead);
    for(int i=0; i<trail.size(); ++i){
        if(qhead == i+1){
            printf(" --> ");
        }
        printf("%d(%c){L%d} ",to_dimacs(trail[i]), to_char(value(trail[i])), level(var(trail[i])));
    }
    printf("\n");
}

void SATSolver::vis_assumption(){
    for(int i=0; i<assumptions.size(); ++i){
        printf("%d ",to_dimacs(assumptions[i]));
    }puts(" ");
}

void SATSolver::vis_trail_full(){
    printf("[%d](Level %d): \n", qhead, decisionLevel());
    for(int i=0; i<trail.size(); ++i){
        if(qhead == i+1){
            printf(" -->");
        }else{
            printf("    ");
        }
        bool in_ass = false;
        for(int j=0; j<=assumptions.size(); ++j){
            if(assumptions[j] == trail[i]){
                in_ass = true;
                break;
            }
        }
        printf("%d(%c){L%d}<A%c> <---",to_dimacs(trail[i]), to_char(value(trail[i])), level(var(trail[i])), (in_ass?'I':'N'));
        Var v = var(trail[i]);
        CRef cref = reason(v);
        if( cref == CRef_Undef){
            printf("vis: CRef_Undef");
        }else{
            vis_clause(cref, false);
        }
        puts("");
    }
    printf("\n");
}


SATSolver::SATSolver() :

    // Parameters (user settable):
    //
    verbosity        (0)
  , var_decay        (opt_var_decay)
  , clause_decay     (opt_clause_decay)
  , random_var_freq  (opt_random_var_freq)
  , random_seed      (opt_random_seed)
  , luby_restart     (opt_luby_restart)
  , ccmin_mode       (opt_ccmin_mode)
  , phase_saving     (opt_phase_saving)
  , rnd_pol          (false)
  , rnd_init_act     (opt_rnd_init_act)
  , garbage_frac     (opt_garbage_frac)
  , min_learnts_lim  (opt_min_learnts_lim)
  , restart_first    (opt_restart_first)
  , restart_inc      (opt_restart_inc)

    // Parameters (the rest):
    //
  , learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

    // Parameters (experimental):
    //
  , learntsize_adjust_start_confl (100)
  , learntsize_adjust_inc         (1.5)

    // Statistics: (formerly in 'SolverStats')
    //
  , solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0)
  , dec_vars(0), num_clauses(0), num_learnts(0), clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

  , watches            (WatcherDeleted(ca))
  , partial_watches    (WatcherDeleted(ca))
  , full_watches       (CADeleted(ca))
  , order_heap         (VarOrderLt(activity))
  , ok                 (true)
  , cla_inc            (1)
  , var_inc            (1)
  , qhead              (0)
  , simpDB_assigns     (-1)
  , simpDB_props       (0)
  , progress_estimate  (0)
  , remove_satisfied   (true)
  , next_var           (0)

    // Resource constraints:
    //
  , conflict_budget    (-1)
  , propagation_budget (-1)
  , asynch_interrupt   (false)
{}


SATSolver::~SATSolver()
{
}


//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var SATSolver::newVar(lbool upol, bool dvar)
{
    Var v;
    if (free_vars.size() > 0){
        v = free_vars.last();
        free_vars.pop();
    }else
        v = next_var++;

    watches         .init(mkLit(v, false));
    watches         .init(mkLit(v, true ));
    partial_watches .init(mkLit(v, false));
    partial_watches .init(mkLit(v, true ));
    full_watches    .init(mkLit(v, false));
    full_watches    .init(mkLit(v, true ));
    assigns  .insert(v, l_Undef);
    vardata  .insert(v, mkVarData(CRef_Undef, 0));
    activity .insert(v, rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    seen     .insert(v, 0);
    polarity .insert(v, true);
    user_pol .insert(v, upol);
    decision .reserve(v);
    trail    .capacity(v+1);

    // printf("Initialize %d %d\n", v, to_dimacs(mkLit(v, false)));

    setDecisionVar(v, dvar);
    return v;
}


// Note: at the moment, only unassigned variable will be released (this is to avoid duplicate
// releases of the same variable).
void SATSolver::releaseVar(Lit l)
{
    if (value(l) == l_Undef){
        addClause(l);
        released_vars.push(var(l));
    }
}


bool SATSolver::addClause_(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);
    if (!ok) return false;

    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
    Lit p; int i, j;
    for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if (value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);

    if (ps.size() == 0)
        return ok = false;
    else if (ps.size() == 1){
        // printf("Add Unit Clause: %c%d \n", sign(ps[0])?'-':' ', var(ps[0])+1);
        uncheckedEnqueue(ps[0]);
        ok = (propagate() == CRef_Undef);
        return ok; 
    }else{
        CRef cr = ca.alloc(ps, false);
        clauses.push(cr);
        attachClause(cr);
    }

    return true;
}


void SATSolver::attachClause(CRef cr){
    Clause& c = ca[cr];
    assert(c.size() > 1);

    if(in_simp){
        
        watches[~c[0]].push(Watcher(cr, c[1]));
        watches[~c[1]].push(Watcher(cr, c[0]));
        if (c.learnt()) num_learnts++, learnts_literals += c.size();
        else            num_clauses++, clauses_literals += c.size();
        
        
    }else{

        if(c.learnt()){
            watches[~c[0]].push(Watcher(cr, c[1]));
            watches[~c[1]].push(Watcher(cr, c[0]));
            num_learnts++, learnts_literals += c.size();
        }else{
            org_clauses_ct += 1;
            unsigned sat_ct = 0;
            unsigned unsat_ct = 0;
            for(int i=0; i<c.size(); ++i){
                Lit l = c[i];
                if(value(c[i])==l_True) sat_ct += 1;
                else if(value(c[i])==l_False) unsat_ct += 1;
                full_watches[l].push(cr);
            }
            c.set_sat_ct(sat_ct);
            c.set_unsat_ct(unsat_ct);


            // printf("%d %d (v %d) sz = %d  cr = %d\n"
            //     , to_dimacs(~c[0]), to_dimacs(~c[1]), nVars(), c.size(), cr); 
            // for(int i=0; i<c.size(); ++i)
            //     printf("%d ", to_dimacs(c[i]));
            // puts("");
            
            // if(to_dimacs(c[0]) == -1498){
            //     vis_partial_watches(1498);
            //     vis_partial_watches(1661);
            //     partial_watches[~c[0]].push(Watcher(cr, c[1]));
            //     partial_watches[~c[1]].push(Watcher(cr, c[0]));
            //     vis_partial_watches(1498);
            //     vis_partial_watches(1661);
            // }
            // fflush(stdout);

            partial_watches[~c[0]].push(Watcher(cr, c[1]));
            partial_watches[~c[1]].push(Watcher(cr, c[0]));

            if(sat_ct > 0)
                sat_org_clauses += 1;
            num_clauses++, clauses_literals += c.size();
        }
    }
}


void SATSolver::detachClause(CRef cr, bool strict){
    Clause& c = ca[cr];
    assert(c.size() > 1);

    if(in_simp){

        if (strict){
            remove(watches[~c[0]], Watcher(cr, c[1]));
            remove(watches[~c[1]], Watcher(cr, c[0]));
        }else{
            watches.smudge(~c[0]);
            watches.smudge(~c[1]);
        }

        if (c.learnt()) num_learnts--, learnts_literals -= c.size();
        else            num_clauses--, clauses_literals -= c.size();

        
    }else{
    
        if(c.learnt()){
            // Strict or lazy detaching:
            if (strict){
                remove(watches[~c[0]], Watcher(cr, c[1]));
                remove(watches[~c[1]], Watcher(cr, c[0]));
            }else{
                watches.smudge(~c[0]);
                watches.smudge(~c[1]);
            }
            num_learnts--, learnts_literals -= c.size();
        }else{       
            if (strict){
                remove(partial_watches[~c[0]], Watcher(cr, c[1]));
                remove(partial_watches[~c[1]], Watcher(cr, c[0]));
            }else{
                partial_watches.smudge(~c[0]);
                partial_watches.smudge(~c[1]);
            }

            org_clauses_ct -= 1;
            if(c.sat_ct()>0) sat_org_clauses -= 1;
            if(strict){
                for(int i=0; i<c.size(); ++i)
                    remove_swap(full_watches[c[i]], cr);
            }else{
                for(int i=0; i<c.size(); ++i)
                    full_watches.smudge(c[i]);
            }
            num_clauses--, clauses_literals -= c.size();
        }

    }
}


void SATSolver::removeClause(CRef cr) {
    Clause& c = ca[cr];
    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if (locked(c)) vardata[var(c[0])].reason = CRef_Undef;
    c.mark(1); 
    ca.free(cr);
}


bool SATSolver::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; }


// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void SATSolver::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            
            Lit     pos_lit = trail[c];
            Var      x  = var(pos_lit);

            if(!in_simp && assigns [x] != l_Undef){
                Lit     neg_lit = ~pos_lit;
                vec<CRef> &wsp = full_watches[pos_lit];
                vec<CRef> &wsn = full_watches[neg_lit];
                for(int i=0; i<wsp.size(); ++i){
                    Clause& c = ca[wsp[i]];
                    unsigned new_sat_ct = c.sat_ct() - 1;
                    c.set_sat_ct(new_sat_ct);
                    if(new_sat_ct == 0) sat_org_clauses -= 1;
                }
                for(int i=0; i<wsn.size(); ++i){
                    Clause& c = ca[wsn[i]];
                    unsigned new_unsat_ct = c.unsat_ct() - 1;
                    c.set_unsat_ct(new_unsat_ct);
                }
            }
            assigns [x] = l_Undef;
            if (phase_saving > 1 || (phase_saving == 1 && c > trail_lim.last()))
                polarity[x] = sign(pos_lit);
            insertVarOrder(x); }
        qhead = trail_lim[level];
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
    } 
}

void SATSolver::finalPartial() {
    int reduce_num = 0;
    // printf("trail size = %d\n", trail_lim.size());
    if(trail_lim.size() > 0){
        
        int bt_pos = trail_lim[0];
        for (int c = trail.size()-1; c >= bt_pos; c--){
            Lit     pos_lit = trail[c];
            Lit     neg_lit = ~pos_lit;
            vec<CRef> &wsp = full_watches[pos_lit];
            vec<CRef> &wsn = full_watches[neg_lit];

            // printf("Full Lit: %c%d :", sign(pos_lit)?'-':' ', var(pos_lit)+1);
            // for(int i=0; i<wsp.size();i++){
            //     Clause& c = ca[wsp[i]];
            //     printf("[%d,%d]",c.sat_ct(),c.unsat_ct());
            //     for(int j=0; j<c.size(); j++){
            //         printf("%c%d ", sign(c[j])?'-':' ', var(c[j])+1);
            //     }
            //     printf("; ");
            // }printf("\n");
            
            bool can_set_undef = true;
            for(int i=0; i<wsp.size(); ++i){
                Clause& c = ca[wsp[i]];
                if(c.sat_ct() <= 1)
                    can_set_undef = false;
                
            }

            if(can_set_undef){
                reduce_num += 1;
                for(int i=0; i<wsp.size(); ++i){
                    Clause& c = ca[wsp[i]];
                    unsigned new_sat_ct = c.sat_ct() - 1;
                    c.set_sat_ct(new_sat_ct);
                    assert(new_sat_ct > 0);
                }
                for(int i=0; i<wsn.size(); ++i){
                    Clause& c = ca[wsn[i]];
                    unsigned new_unsat_ct = c.unsat_ct() - 1;
                    c.set_unsat_ct(new_unsat_ct);
                }
                Var      x  = var(pos_lit);
                assigns[x] = l_Undef;
            }
        }
    }
    printf("c final reduce number = %d\n",reduce_num);
}



//=================================================================================================
// Major methods:


Lit SATSolver::pickBranchLit()
{
    Var next = var_Undef;

    // Random decision:
    if (drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed,order_heap.size())];
        if (value(next) == l_Undef && decision[next])
            rnd_decisions++; }

    // Activity based decision:
    while (next == var_Undef || value(next) != l_Undef || !decision[next])
        if (order_heap.empty()){
            next = var_Undef;
            break;
        }else
            next = order_heap.removeMin();

    // Choose polarity based on different polarity modes (global or per-variable):
    if (next == var_Undef)
        return lit_Undef;
    else if (user_pol[next] != l_Undef)
        return mkLit(next, user_pol[next] == l_True);
    else if (rnd_pol)
        return mkLit(next, drand(random_seed) < 0.5);
    else
        return mkLit(next, polarity[next]);
}


/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|  
|  Description:
|    Analyze conflict and produce a reason clause.
|  
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|  
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
|        rest of literals. There may be others from the same level though.
|  
|________________________________________________________________________________________________@*/
void SATSolver::analyze(CRef confl, vec<Lit>& out_learnt, int& out_btlevel)
{
    int pathC = 0;
    Lit p     = lit_Undef;

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

    do{
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        if (c.learnt())
            claBumpActivity(c);

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && level(var(q)) > 0){
                varBumpActivity(var(q));
                seen[var(q)] = 1;
                if (level(var(q)) >= decisionLevel())
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }
        
        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i]))
                out_learnt[j++] = out_learnt[i];
        
    }else if (ccmin_mode == 1){
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                for (int k = 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(var(p));
    }

    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}


// Check if 'p' can be removed from a conflict clause.
bool SATSolver::litRedundant(Lit p)
{
    enum { seen_undef = 0, seen_source = 1, seen_removable = 2, seen_failed = 3 };
    assert(seen[var(p)] == seen_undef || seen[var(p)] == seen_source);
    assert(reason(var(p)) != CRef_Undef);

    Clause*               c     = &ca[reason(var(p))];
    vec<ShrinkStackElem>& stack = analyze_stack;
    stack.clear();

    for (uint32_t i = 1; ; i++){
        if (i < (uint32_t)c->size()){
            // Checking 'p'-parents 'l':
            Lit l = (*c)[i];
            
            // Variable at level 0 or previously removable:
            if (level(var(l)) == 0 || seen[var(l)] == seen_source || seen[var(l)] == seen_removable){
                continue; }
            
            // Check variable can not be removed for some local reason:
            if (reason(var(l)) == CRef_Undef || seen[var(l)] == seen_failed){
                stack.push(ShrinkStackElem(0, p));
                for (int i = 0; i < stack.size(); i++)
                    if (seen[var(stack[i].l)] == seen_undef){
                        seen[var(stack[i].l)] = seen_failed;
                        analyze_toclear.push(stack[i].l);
                    }
                    
                return false;
            }

            // Recursively check 'l':
            stack.push(ShrinkStackElem(i, p));
            i  = 0;
            p  = l;
            c  = &ca[reason(var(p))];
        }else{
            // Finished with current element 'p' and reason 'c':
            if (seen[var(p)] == seen_undef){
                seen[var(p)] = seen_removable;
                analyze_toclear.push(p);
            }

            // Terminate with success if stack is empty:
            if (stack.size() == 0) break;
            
            // Continue with top element on stack:
            i  = stack.last().i;
            p  = stack.last().l;
            c  = &ca[reason(var(p))];

            stack.pop();
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|  
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void SATSolver::analyzeFinal(Lit p, LSet& out_conflict)
{
    out_conflict.clear();
    out_conflict.insert(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason(x) == CRef_Undef){
                assert(level(x) > 0);
                out_conflict.insert(~trail[i]);
            }else{
                Clause& c = ca[reason(x)];
                for (int j = 1; j < c.size(); j++)
                    if (level(var(c[j])) > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}


void SATSolver::uncheckedEnqueue(Lit p, CRef from)
{
    if( (!in_simp && use_PA) || decisionLevel() == 0){
        Lit pos_lit = p;
        Lit neg_lit = ~p;
        vec<CRef> &wsp = full_watches.lookup(pos_lit);
        vec<CRef> &wsn = full_watches.lookup(neg_lit);
        for(int i=0; i<wsp.size(); ++i){
            Clause& c = ca[wsp[i]];
            int new_sat_ct = c.sat_ct() + 1;
            c.set_sat_ct(new_sat_ct);
            if(new_sat_ct == 1)
                sat_org_clauses += 1;
        }
        for(int i=0; i<wsn.size(); ++i){
            Clause& c = ca[wsn[i]];
            int new_unsat_ct = c.unsat_ct() + 1;
            c.set_unsat_ct(new_unsat_ct);
        }
    }

    assert(value(p) == l_Undef);
    assigns[var(p)] = lbool(!sign(p));
    vardata[var(p)] = mkVarData(from, decisionLevel());
    trail.push_(p);
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|  
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise CRef_Undef.
|  
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
CRef SATSolver::propagate()
{
    // printf("In PROP\n");
    CRef    confl     = CRef_Undef;
    int     num_props = 0;

    while (qhead < trail.size()){
        Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
        num_props++;

if(!in_simp){

    if(use_PA){
        Lit neg_lit = ~p;
        vec<CRef> &wsn = full_watches[neg_lit];
        for(int z=0; z<wsn.size(); z++){
            CRef cref = wsn[z];
            Clause &cls = ca[cref];
            if((int)cls.unsat_ct() == cls.size() - 1){
                if(cls.sat_ct() == 1)
                    continue;
                int nfalse_idx = -1;
                for(int w=0; w<cls.size(); ++w){
                    if(value(cls[w]) == l_Undef){
                        nfalse_idx = w;
                        break;
                    }
                }
                assert(nfalse_idx != -1);
                Lit ltmp = cls[nfalse_idx];
                cls[nfalse_idx] = cls[0];
                cls[0] = ltmp;
                assert(value(cls[0]) == l_Undef);
                uncheckedEnqueue(cls[0], cref);
            }else if((int)cls.unsat_ct() == cls.size()){
                propagations += num_props;
                simpDB_props -= num_props;
                return cref;
            }
        }
    
    }else{

        assert(confl == CRef_Undef);
        // 
        Watcher        *pi, *pj, *pend;
        vec<Watcher>&  pw  = partial_watches.lookup(p);
        for (pi = pj = (Watcher*)pw, pend = pi + pw.size();  pi != pend;){
            // Try to avoid inspecting the clause:
            Lit blocker = pi->blocker;
            if (value(blocker) == l_True){
                *pj++ = *pi++; continue; }
            
            // Make sure the false literal is data[1]:
            CRef     cr        = pi->cref;
            Clause&  c         = ca[cr];

            // printf("visit CLS: ");
            // for(int o=0; o<c.size(); ++o)
            //     printf("%d ", to_dimacs(c[o]));
            // puts("");

            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            pi++;

            // If 0th watch is true, then clause is already satisfied.
            Lit     first = c[0];
            Watcher w     = Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
                *pj++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    partial_watches[~c[1]].push(w);
                    goto pNextClause; }

            // Did not find watch -- clause is unit under assignment:
            *pj++ = w;
            if (value(first) == l_False){
                confl = cr;
                qhead = trail.size();
                // Copy the remaining watches:
                while (pi < pend)
                    *pj++ = *pi++;
            }else
                uncheckedEnqueue(first, cr);

        pNextClause:;
        }
        pw.shrink(pi - pj);
        assert(pw.size()>0);
        if(confl != CRef_Undef){
            // puts("PROP CONFL");
            return confl;
        }
    }
}

        vec<Watcher>&  ws  = watches.lookup(p);
        Watcher        *i, *j, *end;
        
        for (i = j = (Watcher*)ws, end = i + ws.size();  i != end;){
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (value(blocker) == l_True){
                *j++ = *i++; continue; }

            // Make sure the false literal is data[1]:
            CRef     cr        = i->cref;
            Clause&  c         = ca[cr];
            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            i++;

            // If 0th watch is true, then clause is already satisfied.
            Lit     first = c[0];
            Watcher w     = Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
                *j++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    watches[~c[1]].push(w);
                    goto NextClause; }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (value(first) == l_False){
                confl = cr;
                qhead = trail.size();
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            }else
                uncheckedEnqueue(first, cr);

        NextClause:;
        }
        ws.shrink(i - j);
    }
    propagations += num_props;
    simpDB_props -= num_props;

    return confl;
}


/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|  
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { 
    ClauseAllocator& ca;
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
    bool operator () (CRef x, CRef y) { 
        return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity()); } 
};
void SATSolver::reduceDB()
{
    int     i, j;
    double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity

    sort(learnts, reduceDB_lt(ca));
    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // and clauses with activity smaller than 'extra_lim':
    for (i = j = 0; i < learnts.size(); i++){
        Clause& c = ca[learnts[i]];
        if (c.size() > 2 && !locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim))
            removeClause(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    learnts.shrink(i - j);
    checkGarbage();
}


void SATSolver::removeSatisfied(vec<CRef>& cs)
{
    int i, j;
    for (i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];
        if (satisfied(c))
            removeClause(cs[i]);
        else{
            // Trim clause:
            if(c.learnt()){
                if(c.learnt()) assert(value(c[0]) == l_Undef && value(c[1]) == l_Undef);
                for (int k = 2; k < c.size(); k++)
                    if (value(c[k]) == l_False){
                        c[k--] = c[c.size()-1];
                        c.pop();
                    }
            }else{
                for (int k = 0; k < c.size(); k++)
                    if (value(c[k]) == l_False){
                        c[k--] = c[c.size()-1];
                        c.pop();
                        c.set_unsat_ct(c.unsat_ct()-1);
                    }
            }
            cs[j++] = cs[i];
        }
    }
    cs.shrink(i - j);
}


void SATSolver::rebuildOrderHeap()
{
    vec<Var> vs;
    for (Var v = 0; v < nVars(); v++)
        if (decision[v] && value(v) == l_Undef)
            vs.push(v);
    order_heap.build(vs);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|  
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool SATSolver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != CRef_Undef)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    removeSatisfied(learnts);
    if (remove_satisfied){       // Can be turned off.
        removeSatisfied(clauses);

        // TODO: what todo in if 'remove_satisfied' is false?

        // Remove all released variables from the trail:
        for (int i = 0; i < released_vars.size(); i++){
            assert(seen[released_vars[i]] == 0);
            seen[released_vars[i]] = 1;
        }

        int i, j;
        for (i = j = 0; i < trail.size(); i++)
            if (seen[var(trail[i])] == 0)
                trail[j++] = trail[i];
        trail.shrink(i - j);
        //printf("trail.size()= %d, qhead = %d\n", trail.size(), qhead);
        qhead = trail.size();

        for (int i = 0; i < released_vars.size(); i++)
            seen[released_vars[i]] = 0;

        // Released variables are now ready to be reused:
        append(released_vars, free_vars);
        released_vars.clear();
    }
    checkGarbage();
    rebuildOrderHeap();

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}


/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|  
|  Description:
|    Search for a model the specified number of conflicts. 
|    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
|  
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool SATSolver::search(int nof_conflicts)
{
    // if(get_time_since_begin() > cutoff_time){
    //     printf("c cut by cutoff\n");
    //     pa_show_statics();
    //     exit(0);
    // }
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause;
    starts++;

    for (;;){
        CRef confl = propagate();
        if (confl != CRef_Undef){
            // CONFLICT
            conflicts++; conflictC++;
            if (decisionLevel() == 0) return l_False;

            learnt_clause.clear();
            analyze(confl, learnt_clause, backtrack_level);
            cancelUntil(backtrack_level);

            if (learnt_clause.size() == 1){
                uncheckedEnqueue(learnt_clause[0]);
            }else{
                CRef cr = ca.alloc(learnt_clause, true);
                learnts.push(cr);
                attachClause(cr);
                claBumpActivity(ca[cr]);
                uncheckedEnqueue(learnt_clause[0], cr);
            }

            varDecayActivity();
            claDecayActivity();

            if (--learntsize_adjust_cnt == 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
                max_learnts             *= learntsize_inc;

                if (verbosity >= 1)
                    printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", 
                           (int)conflicts, 
                           (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(), (int)clauses_literals, 
                           (int)max_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progressEstimate()*100);
            }

        }else{

            // NO CONFLICT
            if ((nof_conflicts >= 0 && conflictC >= nof_conflicts) || !withinBudget()){
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                cancelUntil(0);
                return l_Undef; }

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0 && !simplify())
                return l_False;

            if (learnts.size()-nAssigns() >= max_learnts)
                // Reduce the set of learnt clauses:
                reduceDB();

            Lit next = lit_Undef;
            while (decisionLevel() < assumptions.size()){
                // Perform user provided assumption:
                Lit p = assumptions[decisionLevel()];
                if (value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                }else if (value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){

                if (!in_simp && use_PA && sat_org_clauses == org_clauses_ct)
                    return l_True;
                
                // New variable decision:
                decisions++;
                next = pickBranchLit();

                if (next == lit_Undef)
                    // Model found:
                    return l_True;
            }

            // Increase decision level and enqueue 'next'
            newDecisionLevel();
            uncheckedEnqueue(next);
        }
    }
}


double SATSolver::progressEstimate() const
{
    double  progress = 0;
    double  F = 1.0 / nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}

/*
  Finite subsequences of the Luby-sequence:

  0: 1
  1: 1 1 2
  2: 1 1 2 1 1 2 4
  3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
  ...


 */

static double luby(double y, int x){

    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}


#include <chrono>

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool SATSolver::solve_()
{

    // printf("here\n");
    using namespace std::chrono;
    auto start = system_clock::now();
    model.clear();
    conflict.clear();
    if (!ok) {
        UNSAT_ct += 1;
        UNSAT_time += 0;
        return l_False;
    }
    solves++;
    
    // print clauses
    // printf("p cnf %d %d \n", nVars(), nClauses());
    // if(nVars() == 19)
    //     for(int i=0; i<clauses.size(); ++i){
    //         CRef &cref = clauses[i];
    //         Clause &c = ca[cref];
    //         for(int j=0; j<c.size(); j++){
    //             printf("%c%d ", sign(c[j])?'-':' ', var(c[j])+1);
    //         }
    //         printf("\n");
    //     }
    


    max_learnts = nClauses() * learntsize_factor;
    if (max_learnts < min_learnts_lim)
        max_learnts = min_learnts_lim;

    learntsize_adjust_confl   = learntsize_adjust_start_confl;
    learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
    lbool   status            = l_Undef;

//    if (verbosity >= 1){
//        printf("============================[ Search Statistics ]==============================\n");
//        printf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
//        printf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
//        printf("===============================================================================\n");
//    }

    // Search:
    int curr_restarts = 0;
    while (status == l_Undef){
        double rest_base = luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts);
        status = search(rest_base * restart_first);
        if (!withinBudget()) break;
        curr_restarts++;
    }

//    if (verbosity >= 1)
//        printf("===============================================================================\n");
//
    
    int n_undef_ct = 0;
    if (status == l_True){
        finalPartial();

        // Extend & copy model:
        model.growTo(nVars());
//        printf("v ");
        for (int i = 0; i < nVars(); i++) {
            model[i] = value(i);
            if(model[i] == l_Undef){
                n_undef_ct += 1;
//                printf("x%d ",i+1);x
            }
            else{
//                if(model[i] == l_False)
//                    printf("-");
//                printf("%d ",i+1);
            }
        }
        printf("c PA %d in %d (undef_ct = %d )\n", nVars()-n_undef_ct, nVars(), n_undef_ct);
    }else if (status == l_False && conflict.size() == 0)
        ok = false;

    cancelUntil(0);

    

    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    double time_in_sec = double(duration.count()) * microseconds::period::num / microseconds::period::den;
    if(status == l_True){
        SAT_time += time_in_sec;
        SAT_ct += 1;
        pa_cal_ct += 1;
        pa_full_sz += nVars();
        pa_final_sz += (nVars()-n_undef_ct);
        pa_cal_time += time_in_sec;
    }else if(status == l_False){
        UNSAT_time += time_in_sec;
        UNSAT_ct += 1;
    }


    return status;
}


bool SATSolver::implies(const vec<Lit>& assumps, vec<Lit>& out)
{
    trail_lim.push(trail.size());
    for (int i = 0; i < assumps.size(); i++){
        Lit a = assumps[i];

        if (value(a) == l_False){
            cancelUntil(0);
            return false;
        }else if (value(a) == l_Undef)
            uncheckedEnqueue(a);
    }

    unsigned trail_before = trail.size();
    bool     ret          = true;
    if (propagate() == CRef_Undef){
        out.clear();
        for (int j = trail_before; j < trail.size(); j++)
            out.push(trail[j]);
    }else
        ret = false;
    
    cancelUntil(0);
    return ret;
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max)
{
    if (map.size() <= x || map[x] == -1){
        map.growTo(x+1, -1);
        map[x] = max++;
    }
    return map[x];
}


void SATSolver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max)
{
    if (satisfied(c)) return;

    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) != l_False)
            fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max)+1);
    fprintf(f, "0\n");
}


void SATSolver::toDimacs(const char *file, const vec<Lit>& assumps)
{
    FILE* f = fopen(file, "wr");
    if (f == NULL)
        fprintf(stderr, "could not open file %s\n", file), exit(1);
    toDimacs(f, assumps);
    fclose(f);
}


void SATSolver::toDimacs(FILE* f, const vec<Lit>& assumps)
{
    // Handle case when solver is in contradictory state:
    if (!ok){
        fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
        return; }

    vec<Var> map; Var max = 0;

    // Cannot use removeClauses here because it is not safe
    // to deallocate them at this point. Could be improved.
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]]))
            cnt++;
        
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]])){
            Clause& c = ca[clauses[i]];
            for (int j = 0; j < c.size(); j++)
                if (value(c[j]) != l_False)
                    mapVar(var(c[j]), map, max);
        }

    // Assumptions are added as unit clauses:
    cnt += assumps.size();

    fprintf(f, "p cnf %d %d\n", max, cnt);

    for (int i = 0; i < assumps.size(); i++){
        assert(value(assumps[i]) != l_False);
        fprintf(f, "%s%d 0\n", sign(assumps[i]) ? "-" : "", mapVar(var(assumps[i]), map, max)+1);
    }

    for (int i = 0; i < clauses.size(); i++)
        toDimacs(f, ca[clauses[i]], map, max);

    if (verbosity > 0)
        printf("Wrote DIMACS with %d variables and %d clauses.\n", max, cnt);
}


void SATSolver::printStats() const
{
    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    printf("restarts              : %" PRIu64"\n", starts);
    printf("conflicts             : %-12" PRIu64"   (%.0f /sec)\n", conflicts   , conflicts   /cpu_time);
    printf("decisions             : %-12" PRIu64"   (%4.2f %% random) (%.0f /sec)\n", decisions, (float)rnd_decisions*100 / (float)decisions, decisions   /cpu_time);
    printf("propagations          : %-12" PRIu64"   (%.0f /sec)\n", propagations, propagations/cpu_time);
    printf("conflict literals     : %-12" PRIu64"   (%4.2f %% deleted)\n", tot_literals, (max_literals - tot_literals)*100 / (double)max_literals);
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
    printf("CPU time              : %g s\n", cpu_time);
}


//=================================================================================================
// Garbage Collection methods:

void SATSolver::relocAll(ClauseAllocator& to)
{
    // All watchers:
    //
    watches.cleanAll();
    full_watches.cleanAll();
    partial_watches.cleanAll();
    for (int v = 0; v < nVars(); v++)
        for (int s = 0; s < 2; s++){
            Lit p = mkLit(v, s);
            vec<Watcher>& ws = watches[p];
            for (int j = 0; j < ws.size(); j++)
                ca.reloc(ws[j].cref, to);
            vec<Watcher>& pw = partial_watches[p];
            for (int j = 0; j < pw.size(); j++)
                ca.reloc(pw[j].cref, to);
            vec<CRef>& fw = full_watches[p];
            for (int j=0; j<fw.size(); ++j)
                ca.reloc(fw[j], to);
        }

    // All reasons
    //
    for (int i = 0; i < trail.size(); i++){
        Var v = var(trail[i]);

        // Note: it is not safe to call 'locked()' on a relocated clause. This is why we keep
        // 'dangling' reasons here. It is safe and does not hurt.
        if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)]))){
            assert(!isRemoved(reason(v)));
            ca.reloc(vardata[v].reason, to);
        }
    }

    // All learnt:
    //
    int i, j;
    for (i = j = 0; i < learnts.size(); i++)
        if (!isRemoved(learnts[i])){
            ca.reloc(learnts[i], to);
            learnts[j++] = learnts[i];
        }
    learnts.shrink(i - j);

    // All original:
    //
    for (i = j = 0; i < clauses.size(); i++)
        if (!isRemoved(clauses[i])){
            ca.reloc(clauses[i], to);
            clauses[j++] = clauses[i];
        }
    clauses.shrink(i - j);
}


void SATSolver::garbageCollect()
{
    // printf("In GC\n");
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted()); 

    relocAll(to);
    if (verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
               ca.size()*ClauseAllocator::Unit_Size, to.size()*ClauseAllocator::Unit_Size);
    to.moveTo(ca);
}
