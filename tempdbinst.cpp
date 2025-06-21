#include "tempdbinst.h"


/*
    class TempDBInst
*/

TempDBInst::TempDBInst(std::string dir_path,std::string dataset_name, std::string constr_name):DBInst(dir_path, dataset_name, constr_name){
    std::cout << "TempDBInst constructed" << std::endl;
    //fill Ts, Te
    assert(schema_map.find("Ts") != schema_map.end() && schema_map.find("Te") != schema_map.end());
    std::vector<Interval> intervals;
    for(int i=0; i<tuples.size(); i++){
        double t_s = std::stod(tuples[i][schema_map["Ts"]]);
        double t_e = std::stod(tuples[i][schema_map["Te"]]);
        timestamps.insert(t_s);
        timestamps.insert(t_e);
        intervals.push_back(std::pair<double,double>(t_s,t_e));
    }
    interval_len_of_all_tuples = MeasureOfUnion(intervals);
    alpha = 2;

    // Initialize tuple_to_violations
    tuple_to_violations.assign(tuples.size(), {});
    for (int i = 0; i < violations.size(); ++i) {
        for (int tid : violations[i]) {
            tuple_to_violations[tid].push_back(i);
        }
    }

    // 初始化 tuple_to_adjs
    tuple_to_adjs.assign(tuples.size(), {});
    for (const auto& adj_pair : _adj) {
        int u = adj_pair.first;
        int v = adj_pair.second;
        tuple_to_adjs[u].push_back(u);
        tuple_to_adjs[v].push_back(u);
    }
}

TempDBInst::~TempDBInst(){
    std::cout << "TempDBInst destructed" << std::endl;
}


// int TempDBInst::_valid(std::vector<int> tupleIndex){
//     double maxTs = 0, minTe = _WEIGHT_INF;
//     for(int i: tupleIndex){
//         maxTs = std::max(maxTs, get_ts(i));
//         minTe = std::min(minTe, get_te(i));
//     }
//     return maxTs<minTe;
// }

void TempDBInst::Normalize() {
    // 1. 收集所有分割点（所有元组的 Ts 和 Te）
    std::vector<double> split_points(timestamps.begin(), timestamps.end());
    std::sort(split_points.begin(), split_points.end());
    split_points.erase(std::unique(split_points.begin(), split_points.end()), split_points.end());

    // 2. 为每个元组记录其原始索引、Ts、Te
    struct TupleSeg {
        int idx;
        double ts, te;
        TupleSeg(int i, double s, double e) : idx(i), ts(s), te(e) {}
    };
    std::vector<TupleSeg> segs;
    for (int i = 0; i < tuples.size(); ++i) {
        segs.emplace_back(i, get_ts(i), get_te(i));
    }

    // 3. 按分割点扫描，分割区间
    size_t orig_tuple_count = tuples.size();
    for (const auto& sp : split_points) {
        for (size_t i = 0; i < segs.size(); ++i) {
            auto& seg = segs[i];
            // 只处理还未分割且分割点在区间内部的元组
            if (sp > seg.ts && sp < seg.te) {
                // 1. 复制元组
                tuple new_tuple = tuples[seg.idx];
                tuples.push_back(new_tuple);
                weights.push_back(weights[seg.idx] * (seg.te - sp) / (seg.te - seg.ts));
                weights[seg.idx] = weights[seg.idx] * (sp - seg.ts) / (seg.te - seg.ts);
                // 2. 更新原元组终止时间
                set_te(seg.idx, sp);
                // 3. 新元组起止时间
                int new_idx = tuples.size() - 1;
                set_ts(new_idx, sp);
                set_te(new_idx, seg.te);
                // 4. 记录邻接关系
                _adj[seg.idx] = new_idx;
                // _adj_weights[seg.idx] = (weights[seg.idx] + weights[new_idx]);
                // 5. 更新 segs
                seg.te = sp;
                segs.emplace_back(new_idx, sp, seg.te);
            }
        }
    }

    //reweight adj weights
    for(auto adj_pair : _adj){
        _adj_weights[adj_pair.first] =  alpha;
    }
}


double TempDBInst::CalcRepair() {
    double cost = 0.0;
    if (tuples.empty()) return 0;
    // CalcViolations();
    if (violations.empty()) return 0;
    calc_vio_duration();
    std::cout<<"begin repair\n";

    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int V = violations.size();
        int A = _adj.size();
        std::cout<<"N="<<N<<", V="<<V<<", A="<<A<<std::endl;

        // 变量
        IloNumVarArray x(env, N, 0, 1, ILOFLOAT); // relaxed LP
        IloNumVarArray yV(env, V, 0, 1, ILOFLOAT);
        IloNumVarArray yA(env, A, 0, 1, ILOFLOAT);

        // 目标函数
        IloExpr obj(env);
        for (int i = 0; i < N; i++) obj += weights[i] * x[i];
        for (int i = 0; i < V; i++)
            if (vio_weights[i] < _WEIGHT_INF)
                obj += vio_weights[i] * vio_duration[i] * yV[i];
        int adj_pt = 0;
        std::map<int,int> adj_index_to_adj_pt;
        for (auto adj_pair : _adj) {
            obj += _adj_weights[adj_pair.first] * yA[adj_pt];
            adj_index_to_adj_pt[adj_pair.first] = adj_pt;
            adj_pt++;
        }
        model.add(IloMinimize(env, obj));
        obj.end();

        // 约束
        // 1. 违反约束
        for (int i = 0; i < V; i++) {
            IloExpr sumx(env);
            for (int tuple_id : violations[i]) sumx += x[tuple_id];
            if (vio_weights[i] == _WEIGHT_INF) {
                model.add(sumx >= 1);
            } else {
                model.add(sumx + yV[i] >= 1);
            }
            sumx.end();
        }
        // 2. 邻接约束
        adj_pt = 0;
        for (auto adj_pair : _adj) {
            model.add(yA[adj_pt] - x[adj_pair.first] >= 0);
            model.add(yA[adj_pt] - x[adj_pair.second] >= 0);
            adj_pt++;
        }

        // 求解
        IloCplex cplex(model);
        // cplex.setParam(IloCplex::Threads, 8); // 线程数可调
        cplex.setOut(env.getNullStream());    // 不输出日志
        bool ok = cplex.solve();
        std::cout<<"solved\n";
        if (!ok) {
            std::cerr << "Cplex failed!" << std::endl;
            env.end();
            return -1;
        }

        // 取解并舍入
        std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne;
        calculated_repair.clear();

        // x_t
        for (int i = 0; i < N; i++) {
            double val = cplex.getValue(x[i]);
            if (val < 1.0 / (lambda + 1)) {
                calculated_repair.push_back(i);
            } else {
                tupleIndexRoundedToOne.insert(i);
            }
        }
        // y_V
        for (int i = 0; i < V; i++) {
            if (lambda == 2) {
                int x_sum = 0;
                for (int tup_index : violations[i]) {
                    if (tupleIndexRoundedToOne.find(tup_index) != tupleIndexRoundedToOne.end()) x_sum++;
                }
                if (std::max(1 - x_sum, 0) == 1) violationIndexRoundedToOne.insert(i);
            } else {
                double y = cplex.getValue(yV[i]);
                if (y > 1.0 / (lambda + 1)) violationIndexRoundedToOne.insert(i);
            }
        }
        // y_A
        for (auto adj_pair : _adj) {
            if (tupleIndexRoundedToOne.find(adj_pair.first) != tupleIndexRoundedToOne.end() ||
                tupleIndexRoundedToOne.find(adj_pair.second) != tupleIndexRoundedToOne.end()) {
                adjIndexRoundedToOne.insert(adj_pair.first);
                // std::cout<< cplex.getValue(x[adj_pair.first])<<","<< cplex.getValue(x[adj_pair.second])<<" : "<<cplex.getValue(yA[adj_index_to_adj_pt[adj_pair.first]])<<"->1"<<std::endl;
                assert(cplex.getValue(yA[adj_index_to_adj_pt[adj_pair.first]])>= std::max(cplex.getValue(x[adj_pair.first]),cplex.getValue(x[adj_pair.second])));
            }
        }

        // 重新计算cost
        cost = 0.0;
        // x_t
        for (int i : tupleIndexRoundedToOne) {
            cost += weights[i];
        }
        // y_V
        for (int i : violationIndexRoundedToOne) {
            cost += vio_weights[i] * vio_duration[i];
        }
        // y_A
        weighted_sum_of_broken_adj = 0;
        for (int i : adjIndexRoundedToOne) {
            cost += _adj_weights[i];
            weighted_sum_of_broken_adj += 1;
        }
    } catch (IloException& e) {
        std::cerr << "CPLEX exception: " << e.getMessage() << std::endl;
        env.end();
        return -1;
    }
    env.end();
    return cost;
}
double TempDBInst::CalcRepairWOContConstr() {
    double cost = 0.0;
    if (tuples.empty()) return 0;
    // CalcViolations();
    if (violations.empty()) return 0;
    calc_vio_duration();

    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int V = violations.size();

        // 变量
        IloNumVarArray x(env, N, 0, 1, ILOFLOAT);
        IloNumVarArray yV(env, V, 0, 1, ILOFLOAT);

        // 目标函数
        IloExpr obj(env);
        for (int i = 0; i < N; i++) obj += weights[i] * x[i];
        for (int i = 0; i < V; i++)
            if (vio_weights[i] < _WEIGHT_INF)
                obj += vio_weights[i] * vio_duration[i] * yV[i];
        model.add(IloMinimize(env, obj));
        obj.end();

        // 约束
        for (int i = 0; i < V; i++) {
            IloExpr sumx(env);
            for (int tuple_id : violations[i]) sumx += x[tuple_id];
            if (vio_weights[i] == _WEIGHT_INF) {
                model.add(sumx >= 1);
            } else {
                model.add(sumx + yV[i] >= 1);
            }
            sumx.end();
        }

        // 求解
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads, 8);
        cplex.setOut(env.getNullStream());
        bool ok = cplex.solve();
        if (!ok) {
            std::cerr << "Cplex failed!" << std::endl;
            env.end();
            return -1;
        }

        // 取解并舍入
        std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne;
        calculated_repair.clear();

        // x_t
        for (int i = 0; i < N; i++) {
            double val = cplex.getValue(x[i]);
            if (val < 1.0 / (lambda + 1)) {
                calculated_repair.push_back(i);
            } else {
                tupleIndexRoundedToOne.insert(i);
            }
        }
        // y_V
        for (int i = 0; i < V; i++) {
            if (lambda == 2) {
                int x_sum = 0;
                for (int tup_index : violations[i]) {
                    if (tupleIndexRoundedToOne.find(tup_index) != tupleIndexRoundedToOne.end()) x_sum++;
                }
                if (std::max(1 - x_sum, 0) == 1) violationIndexRoundedToOne.insert(i);
            } else {
                double y = cplex.getValue(yV[i]);
                if (y > 1.0 / (lambda + 1)) violationIndexRoundedToOne.insert(i);
            }
        }
        // y_A
        for (auto adj_pair : _adj) {
            if (tupleIndexRoundedToOne.find(adj_pair.first) != tupleIndexRoundedToOne.end() ||
                tupleIndexRoundedToOne.find(adj_pair.second) != tupleIndexRoundedToOne.end()) {
                adjIndexRoundedToOne.insert(adj_pair.first);
            }
        }

        // 重新计算cost
        cost = 0.0;
        // x_t
        for (int i : tupleIndexRoundedToOne) {
            cost += weights[i];
        }
        // y_V
        for (int i : violationIndexRoundedToOne) {
            cost += vio_weights[i] * vio_duration[i];
        }
        // y_A
        weighted_sum_of_broken_adj = 0;
        for (int i : adjIndexRoundedToOne) {
            cost += _adj_weights[i];
            weighted_sum_of_broken_adj += 1;
        }
    } catch (IloException& e) {
        std::cerr << "CPLEX exception: " << e.getMessage() << std::endl;
        env.end();
        return -1;
    }
    env.end();
    return cost;
}
double TempDBInst::ExactRepair() {
    double cost = 0.0;
    if (tuples.empty()) return 0;
    // CalcViolations();
    if (violations.empty()) return 0;
    calc_vio_duration();

    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int V = violations.size();
        int A = _adj.size();

        // 变量（整数0/1）
        IloNumVarArray x(env, N, 0, 1, ILOINT);
        IloNumVarArray yV(env, V, 0, 1, ILOINT);
        IloNumVarArray yA(env, A, 0, 1, ILOINT);

        // 目标函数
        IloExpr obj(env);
        for (int i = 0; i < N; i++) obj += weights[i] * x[i];
        for (int i = 0; i < V; i++)
            if (vio_weights[i] < _WEIGHT_INF)
                obj += vio_weights[i] * vio_duration[i] * yV[i];
        int adj_pt = 0;
        for (auto adj_pair : _adj) {
            obj += _adj_weights[adj_pair.first] * yA[adj_pt];
            adj_pt++;
        }
        model.add(IloMinimize(env, obj));
        obj.end();

        // 约束
        // 1. 违反约束
        for (int i = 0; i < V; i++) {
            IloExpr sumx(env);
            for (int tuple_id : violations[i]) sumx += x[tuple_id];
            if (vio_weights[i] == _WEIGHT_INF) {
                model.add(sumx >= 1);
            } else {
                model.add(sumx + yV[i] >= 1);
            }
            sumx.end();
        }
        // 2. 邻接约束
        adj_pt = 0;
        for (auto adj_pair : _adj) {
            model.add(yA[adj_pt] - x[adj_pair.first] >= 0);
            model.add(yA[adj_pt] - x[adj_pair.second] >= 0);
            adj_pt++;
        }

        // 求解
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads, 8);
        cplex.setOut(env.getNullStream());
        bool ok = cplex.solve();
        if (!ok) {
            std::cerr << "Cplex failed!" << std::endl;
            env.end();
            return -1;
        }

        // 取解并舍入
        std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne;
        calculated_repair.clear();

        // x_t
        for (int i = 0; i < N; i++) {
            int val = static_cast<int>(std::round(cplex.getValue(x[i])));
            if (val == 0) {
                calculated_repair.push_back(i);
            } else {
                tupleIndexRoundedToOne.insert(i);
            }
        }
        // y_V
        for (int i = 0; i < V; i++) {
            int x_sum = 0;
            for (int tup_index : violations[i]) {
                if (tupleIndexRoundedToOne.find(tup_index) != tupleIndexRoundedToOne.end()) x_sum++;
            }
            if (std::max(1 - x_sum, 0) == 1) violationIndexRoundedToOne.insert(i);
        }
        // y_A
        for (auto adj_pair : _adj) {
            int x1 = tupleIndexRoundedToOne.count(adj_pair.first) ? 1 : 0;
            int x2 = tupleIndexRoundedToOne.count(adj_pair.second) ? 1 : 0;
            if (x1 >= 1 || x2 >= 1) {
                adjIndexRoundedToOne.insert(adj_pair.first);
            }
        }

        // 重新计算cost
        cost = 0.0;
        for (int i : tupleIndexRoundedToOne) {
            cost += weights[i];
        }
        for (int i : violationIndexRoundedToOne) {
            cost += vio_weights[i] * vio_duration[i];
        }
        // y_A
        weighted_sum_of_broken_adj = 0;
        for (int i : adjIndexRoundedToOne) {
            cost += _adj_weights[i];
            weighted_sum_of_broken_adj += 1;
        }
    } catch (IloException& e) {
        std::cerr << "CPLEX exception: " << e.getMessage() << std::endl;
        env.end();
        return -1;
    }
    env.end();
    return cost;
}


double TempDBInst::CalcRepair_Postproc() {
    double cost = 0.0;
    if (tuples.empty()) return 0;
    // CalcViolations();
    if (violations.empty()) return 0;
    calc_vio_duration();

    std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne;
    calculated_repair.clear();

    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int V = violations.size();
        int A = _adj.size();

        // 变量
        IloNumVarArray x(env, N, 0, 1, ILOFLOAT); // relaxed LP
        IloNumVarArray yV(env, V, 0, 1, ILOFLOAT);
        IloNumVarArray yA(env, A, 0, 1, ILOFLOAT);

        // 目标函数
        IloExpr obj(env);
        for (int i = 0; i < N; i++) obj += weights[i] * x[i];
        for (int i = 0; i < V; i++)
            if (vio_weights[i] < _WEIGHT_INF)
                obj += vio_weights[i] * vio_duration[i] * yV[i];
        int adj_pt = 0;
        for (auto adj_pair : _adj) {
            obj += _adj_weights[adj_pair.first] * yA[adj_pt];
            adj_pt++;
        }
        model.add(IloMinimize(env, obj));
        obj.end();

        // 约束
        // 1. 违反约束
        for (int i = 0; i < V; i++) {
            IloExpr sumx(env);
            for (int tuple_id : violations[i]) sumx += x[tuple_id];
            if (vio_weights[i] == _WEIGHT_INF) {
                model.add(sumx >= 1);
            } else {
                model.add(sumx + yV[i] >= 1);
            }
            sumx.end();
        }
        // 2. 邻接约束
        adj_pt = 0;
        for (auto adj_pair : _adj) {
            model.add(yA[adj_pt] - x[adj_pair.first] >= 0);
            model.add(yA[adj_pt] - x[adj_pair.second] >= 0);
            adj_pt++;
        }

        // 求解
        IloCplex cplex(model);
        // cplex.setParam(IloCplex::Threads, 8); // 线程数可调
        cplex.setOut(env.getNullStream());    // 不输出日志
        bool ok = cplex.solve();
        if (!ok) {
            std::cerr << "Cplex failed!" << std::endl;
            env.end();
            return -1;
        }

        std::cout<<"cplex finished"<<std::endl;
        
        // 取解并舍入
        // x_t
        for (int i = 0; i < N; i++)
        {
            double val = cplex.getValue(x[i]);
            if (val < 1.0 / (lambda + 1))
            {
                calculated_repair.push_back(i);
            }
            else
            {
                tupleIndexRoundedToOne.insert(i);
            }
        }

        std::cout << "rounding finished" << std::endl;
    } catch (IloException& e) {
        std::cerr << "CPLEX exception: " << e.getMessage() << std::endl;
        env.end();
        return -1;
    }
    env.end();

    

    cost = UpdateCostAndSets(tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne);
    assert(CheckRoundedSolutionValid(tupleIndexRoundedToOne, violationIndexRoundedToOne));

    PostProc(tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne, cost);

    return cost;
}

// 封装cost和相关集合的更新
// tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne, calculated_repair 必须为引用
// 返回cost

double TempDBInst::UpdateCostAndSets(std::set<int>& tupleIndexRoundedToOne,
                                     std::set<int>& violationIndexRoundedToOne,
                                     std::set<int>& adjIndexRoundedToOne) {
    // 重新生成calculated_repair
    calculated_repair.clear();
    for (int i = 0; i < (int)tuples.size(); ++i) {
        if (!tupleIndexRoundedToOne.count(i))
            calculated_repair.push_back(i);
    }
    // y_V
    for (int i = 0; i < violations.size(); i++) {
        int x_sum = 0;
        for (int tup_index : violations[i]) {
            if (tupleIndexRoundedToOne.find(tup_index) != tupleIndexRoundedToOne.end())
                x_sum++;
        }
        if (std::max(1 - x_sum, 0) == 1)
            violationIndexRoundedToOne.insert(i);
    }
    // y_A
    for (auto adj_pair : _adj) {
        if (tupleIndexRoundedToOne.find(adj_pair.first) != tupleIndexRoundedToOne.end() ||
            tupleIndexRoundedToOne.find(adj_pair.second) != tupleIndexRoundedToOne.end()) {
            adjIndexRoundedToOne.insert(adj_pair.first);
        }
    }
    // 重新计算cost（与CalcRepair一致）
    double cost = 0.0;
    for (int i : tupleIndexRoundedToOne) {
        cost += weights[i];
    }
    for (int i : violationIndexRoundedToOne) {
        cost += vio_weights[i] * vio_duration[i];
    }
    weighted_sum_of_broken_adj = 0;
    for (int i : adjIndexRoundedToOne) {
        cost += _adj_weights[i];
        weighted_sum_of_broken_adj += 1;
    }
    return cost;
}

// 增量更新辅助结构：在构造函数或初始化时建立
// tuple_to_violations[tid] = vector of violation indices that contain tid
// 在 TempDBInst 构造函数中添加：
// tuple_to_violations.assign(tuples.size(), {});
// for (int i = 0; i < violations.size(); ++i) {
//     for (int tid : violations[i]) {
//         tuple_to_violations[tid].push_back(i);
//     }
// }
double TempDBInst::UpdateCostAndSetsIncremental(
    std::set<int>& tupleIndexRoundedToOne,
    std::set<int>& violationIndexRoundedToOne,
    std::set<int>& adjIndexRoundedToOne,
    int removed_tid)
{
    // 1. 更新calculated_repair
    calculated_repair.clear();
    for (int i = 0; i < (int)tuples.size(); ++i) {
        if (!tupleIndexRoundedToOne.count(i))
            calculated_repair.push_back(i);
    }

    // 2. 增量更新 violationIndexRoundedToOne
    if (removed_tid >= 0 && removed_tid < (int)tuple_to_violations.size()) {
        for (int vio_idx : tuple_to_violations[removed_tid]) {
            int x_sum = 0;
            for (int tup_index : violations[vio_idx]) {
                if (tupleIndexRoundedToOne.find(tup_index) != tupleIndexRoundedToOne.end())
                    x_sum++;
            }
            if (std::max(1 - x_sum, 0) == 1)
                violationIndexRoundedToOne.insert(vio_idx);
            else
                violationIndexRoundedToOne.erase(vio_idx);
        }
    }

    // 增量更新 adjIndexRoundedToOne
    if (removed_tid >= 0 && removed_tid < (int)tuple_to_adjs.size()) {
        for (int adj_key : tuple_to_adjs[removed_tid]) {
            int u = adj_key;
            int v = _adj[adj_key];
            if (tupleIndexRoundedToOne.count(u) || tupleIndexRoundedToOne.count(v)) {
                adjIndexRoundedToOne.insert(u);
            } else {
                adjIndexRoundedToOne.erase(u);
            }
        }
    }

    // 4. 重新计算cost
    double cost = 0.0;
    for (int i : tupleIndexRoundedToOne) {
        cost += weights[i];
    }
    for (int i : violationIndexRoundedToOne) {
        cost += vio_weights[i] * vio_duration[i];
    }
    weighted_sum_of_broken_adj = 0;
    for (int i : adjIndexRoundedToOne) {
        cost += _adj_weights[i];
        weighted_sum_of_broken_adj += 1;
    }
    return cost;
}

// double TempDBInst::UpdateCostAndSetsIncremental(
//     std::set<int>& tupleIndexRoundedToOne,
//     std::set<int>& violationIndexRoundedToOne,
//     std::set<int>& adjIndexRoundedToOne,
//     int removed_tid)
// {
//     // 1. 更新calculated_repair
//     calculated_repair.clear();
//     for (int i = 0; i < (int)tuples.size(); ++i) {
//         if (!tupleIndexRoundedToOne.count(i))
//             calculated_repair.push_back(i);
//     }

//     // 2. 增量更新 violationIndexRoundedToOne
//     // 只需检查包含 removed_tid 的 violation
//     if (removed_tid >= 0 && removed_tid < (int)tuple_to_violations.size()) {
//         for (int vio_idx : tuple_to_violations[removed_tid]) {
//             int x_sum = 0;
//             for (int tup_index : violations[vio_idx]) {
//                 if (tupleIndexRoundedToOne.find(tup_index) != tupleIndexRoundedToOne.end())
//                     x_sum++;
//             }
//             if (std::max(1 - x_sum, 0) == 1)
//                 violationIndexRoundedToOne.insert(vio_idx);
//             else
//                 violationIndexRoundedToOne.erase(vio_idx);
//         }
//     }

//     // 3. 增量更新 adjIndexRoundedToOne
//     // 只需检查与 removed_tid 相关的邻接约束
//     for (auto& adj_pair : _adj) {
//         if (adj_pair.first != removed_tid && adj_pair.second != removed_tid)
//             continue;
//         if (tupleIndexRoundedToOne.find(adj_pair.first) != tupleIndexRoundedToOne.end() ||
//             tupleIndexRoundedToOne.find(adj_pair.second) != tupleIndexRoundedToOne.end()) {
//             adjIndexRoundedToOne.insert(adj_pair.first);
//         } else {
//             adjIndexRoundedToOne.erase(adj_pair.first);
//         }
//     }

//     // 4. 重新计算cost（与CalcRepair一致）
//     double cost = 0.0;
//     for (int i : tupleIndexRoundedToOne) {
//         cost += weights[i];
//     }
//     for (int i : violationIndexRoundedToOne) {
//         cost += vio_weights[i] * vio_duration[i];
//     }
//     weighted_sum_of_broken_adj = 0;
//     for (int i : adjIndexRoundedToOne) {
//         cost += _adj_weights[i];
//         weighted_sum_of_broken_adj += _adj_weights[i];
//     }
//     return cost;
// }
