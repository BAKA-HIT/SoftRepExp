/*
    dbinst.cpp 

    Created by baka on 3/27

*/

#include <ilcplex/ilocplex.h>
#include "dbinst.h"
#include <algorithm>
#include <iostream>
#include <glpk.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cassert>
#include <queue>
#include <unordered_set>
#include <bitset>
#include <vector>
#include <unordered_map>

// 添加转义函数
auto escape_sql_string = [](const std::string& input) {
    std::string output;
    for (char c : input) {
        if (c == '"' || c == '\'' || c == '/' || c == ';' || c == '\\') {
            c = '_';
        }
        output += c;
    }
    return output;
};

/*
    class DBInst
*/

DBInst::DBInst(std::string dir_path,std::string dataset_name, std::string constr_name) : _dir_path(dir_path), _dataset_name(dataset_name), _constr_name(constr_name){
    //read dataset from dataset_path, the first line being the db schema.
    std::ifstream db_inst(dir_path + dataset_name); 
    std::ifstream dc_inst(dir_path + constr_name);
    
    std::string schema_str;
    std::getline(db_inst, schema_str);
    std::stringstream schema_ss(schema_str);
    std::string s;
    while(std::getline(schema_ss, s, ',')){
        s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
        s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
        schema.push_back(s);
    }

    schema_len = schema.size();

    for(int i=0; i<schema_len; i++){
        schema_map[schema[i]]= i;
    }

    std::string tuple_str;
    while(std::getline(db_inst, tuple_str)){
        tuple t;
        std::stringstream tuple_ss(tuple_str);
        std::string s;
        while(std::getline(tuple_ss, s, ',')){
            s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
            s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
            if(s.empty()) s="#";
            t.push_back(s);
        }

        //enter weight
        std::string weight_str = t.back();
        t.pop_back();
        
        if(t.size() == schema_len){
            tuples.push_back(t);
            
            weights.push_back(std::stod(weight_str));
        }
        else{
            for(auto _t : t)
                std::cout<<_t<<"|";
            std::cout<<std::endl;
            std::cout<<schema_len<<std::endl;
        }

        assert(t.size() == schema_len);
    }

    std::getline(dc_inst, createtable_sql);
    createtable_sql.erase(createtable_sql.find_last_not_of(" \t\n\r\f\v") + 1);
    createtable_sql.erase(0, createtable_sql.find_first_not_of(" \t\n\r\f\v"));

    std::string constr_str;
    while(std::getline(dc_inst, constr_str)){
        std::stringstream constr_ss(constr_str);
        int nTupleInvolved=0;
        double constr_weight = 0;
        std::string s;

        std::getline(constr_ss, s, ',');
        s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
        s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
        nTupleInvolved = std::stoi(s);
        lambda = std::max(lambda, nTupleInvolved);

        std::getline(constr_ss, s, ',');
        s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
        s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
        if(s.compare("inf") == 0) constr_weight = _WEIGHT_INF;
        else constr_weight = std::stod(s);
        constr_weights.push_back(constr_weight);

        std::getline(constr_ss, s);
        s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
        s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
        constraints.push_back(s);
    }

    db_inst.close();
    dc_inst.close();
    
    CalcViolations();
    std::cout<<"csv loaded, tuple number="<<tuples.size()<<", vio number="<<violations.size()<<std::endl;
}


DBInst::~DBInst(){
    if(db) sqlite3_close(db);
}

void DBInst::__calc__violations()
{
    // 清空旧数据
    violations.clear();
    vio_weights.clear();
    // 打开数据库并重建表
    std::string path="./db_temp/";
    #ifndef _EXP_TEST
    path += ".main" + escape_sql_string(_dir_path + _dataset_name + "_" + _constr_name) + ".db";
    #else
    path += ".test" + escape_sql_string(_dir_path + _dataset_name + "_" + _constr_name) + ".db";
    #endif
    std::cout<<path<<std::endl;
    // 如果文件不存在则新建
    std::ifstream f(path);
    if (!f.good()) {
        std::ofstream fnew(path); fnew.close();
    }
    sqlite3_open(path.c_str(), &db);
    
    int rc = 0;
    char *zErrMsg = 0;

    sqlite3_exec(db, "drop table r;", NULL, NULL, NULL);
    rc = sqlite3_exec(db, ("create table r( __id int, " + createtable_sql + ");").c_str(), NULL, NULL, NULL);
    if(rc != SQLITE_OK){
        std::cout << "-----" << std::endl;
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }

    // 插入元组
    // 批量插入优化：用事务包裹所有插入，极大加速
    sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, NULL);
    for (int i = 0; i < tuples.size(); i++)
    {
        std::string sql = "insert into r values(" + std::to_string(i);
        for (int j = 0; j < tuples[i].size(); j++)
        {
            sql += ", \" " + escape_sql_string(tuples[i][j]) + " \"";
        }
        sql += ");";
        rc = sqlite3_exec(db, sql.c_str(), NULL, NULL, &zErrMsg);
        if (rc != SQLITE_OK)
        {
            std::cout << sql << "-----" << std::endl;
            fprintf(stderr, "SQL error: %s\n", zErrMsg);
            sqlite3_free(zErrMsg);
        }
    }
    sqlite3_exec(db, "COMMIT;", NULL, NULL, NULL);

    // 非回调方式执行约束查询
    assert(constraints.size() == constr_weights.size());
    for (int i = 0; i < constraints.size(); i++)
    {
        std::string sql = constraints[i];
        sqlite3_stmt *stmt = nullptr;
        rc = sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL);
        if (rc != SQLITE_OK)
        {
            std::cout << sql << "-----" << std::endl;
            fprintf(stderr, "SQL prepare error: %s\n", sqlite3_errmsg(db));
            continue;
        }
        // 每一行代表一个violation
        while (sqlite3_step(stmt) == SQLITE_ROW)
        {
            std::set<int> vio;
            int ncol = sqlite3_column_count(stmt);
            // assert(ncol == 2);
            for (int col = 0; col < ncol; col++)
            {
                const unsigned char *val = sqlite3_column_text(stmt, col);
                if (val)
                    vio.insert(std::atoi((const char *)val));
            }
            violations.push_back(vio);
            vio_weights.push_back(constr_weights[i]);
        }
        sqlite3_finalize(stmt);
    }
}

void DBInst::CalcViolations() {
    // 优先尝试从文件读取
    std::string viofile = _dir_path + "saved/__vio_saved_" + _constr_name + "_" + _dataset_name;
    std::string weightfile = _dir_path + "saved/__w_saved_" + _constr_name + "_" + _dataset_name;
    std::ifstream vio_ifs(viofile), w_ifs(weightfile);
    bool viofile_ok = vio_ifs.good() && w_ifs.good();
    vio_ifs.close();
    w_ifs.close();

    std::cout << "[INFO] Try loading violations and weights from file...\n";
    if (viofile_ok) {
        
        LoadViolationsAndWeights(viofile, weightfile);
        if (!violations.empty() && violations.size() == vio_weights.size()) {
            std::cout << "[INFO] Loaded violations and weights from file.\n";
        } else {
            std::cout << "[WARN] File exists but data invalid, fallback to calculation.\n";
            __calc__violations();
            SaveViolationsAndWeights(viofile, weightfile);
        }
    }
    else{
        std::cout << "[WARN] File does not exist, fallback to calculation.\n";
        __calc__violations();
        // 计算后保存
        SaveViolationsAndWeights(viofile, weightfile);
    }


}

std::vector<tuple> DBInst::OutputRepair(){
    std::vector<tuple> retTuples;

    for(int i=0; i<calculated_repair.size(); i++){
        retTuples.push_back(tuples[calculated_repair[i]]);
    }

    return retTuples;
}

void DBInst::output(){
    //output tuples
    for(int i=0; i<schema.size(); i++){
        std::cout<<schema[i]<<" \t";
    }
    std::cout<<"weights"<<std::endl;

    for(int i=0; i<tuples.size(); i++){
        tuple curTuple = tuples[i];
        for(int j=0; j<curTuple.size(); j++){
            std::cout<<curTuple[j]<<" \t";
        }
        std::cout<<weights[i]<<std::endl;
    }
}

// CPLEX实现的ExactRepair
// 需在Makefile中加-DUSE_CPLEX并链接CPLEX库

double DBInst::ExactRepair() {
    printf("begin exact repair\n");
    fflush(stdout);
    double cost = 0;
    if (tuples.size() == 0 || violations.size() == 0) return 0;
    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int M = violations.size();
        // 变量
        IloNumVarArray x(env, N, 0, 1, ILOINT); // x_t
        IloNumVarArray y(env, M, 0, 1, ILOINT); // y_V
        // 目标函数
        IloExpr obj(env);
        for (int i = 0; i < N; ++i) obj += weights[i] * x[i];
        for (int i = 0; i < M; ++i) if (vio_weights[i] < _WEIGHT_INF) obj += vio_weights[i] * y[i];
        model.add(IloMinimize(env, obj));
        obj.end();
        // 约束
        for (int i = 0; i < M; ++i) {
            IloExpr sum(env);
            for (int tid : violations[i]) sum += x[tid];
            if (vio_weights[i] == _WEIGHT_INF) {
                model.add(sum >= 1);
            } else {
                model.add(sum + y[i] >= 1);
            }
            sum.end();
        }
        // 求解
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::MIP::Display, 4); 
        cplex.setOut(std::cout);
        bool ok = cplex.solve();
        
        if (!ok) throw std::runtime_error("CPLEX failed to solve ILP");
        // 还原解
        std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne;
        calculated_repair.clear();
        for (int i = 0; i < N; ++i) {
            double xv = cplex.getValue(x[i]);
            if (xv < 0.5) calculated_repair.push_back(i);
            else {
                tupleIndexRoundedToOne.insert(i);
                cost += weights[i];
            }
        }
        for (int i = 0; i < M; ++i) {
            if (vio_weights[i] < _WEIGHT_INF) {
                double yv = cplex.getValue(y[i]);
                if (yv > 0.5) {
                    violationIndexRoundedToOne.insert(i);
                    cost += vio_weights[i];
                }
            }
        }
        assert(CheckRoundedSolutionValid(tupleIndexRoundedToOne, violationIndexRoundedToOne));
    } catch (IloException& e) {
        std::cerr << "CPLEX Exception: " << e << std::endl;
        env.end();
        throw;
    }
    env.end();
    return cost;
}

// MW3SCRepair的CPLEX版本
// min ∑_t w_t x_t + ∑_V w_V y_V
// s.t. ∑_{t in V} x_t - y_V >= 1
// 0 <= x_t, y_V <= 1

double DBInst::MW3SCRepair() {
    double cost = 0.0;
    if (tuples.size() == 0 || violations.size() == 0) return 0;
    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int M = violations.size();
        IloNumVarArray x(env, N, 0, 1, ILOFLOAT);
        IloNumVarArray y(env, M, 0, 1, ILOFLOAT);
        IloExpr obj(env);
        for (int i = 0; i < N; ++i) obj += weights[i] * x[i];
        for (int i = 0; i < M; ++i) obj += vio_weights[i] * y[i];
        model.add(IloMinimize(env, obj));
        obj.end();
        for (int i = 0; i < M; ++i) {
            IloExpr sum(env);
            for (int tid : violations[i]) sum += x[tid];
            sum -= y[i];
            model.add(sum >= 1);
            sum.end();
        }
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads,32);
        cplex.setWarning(env.getNullStream());
        bool ok = cplex.solve();
        if (!ok) throw std::runtime_error("CPLEX failed to solve MW3SCRepair LP");
        calculated_repair.clear();
        for (int i = 0; i < N; ++i) {
            double xv = cplex.getValue(x[i]);
            if (xv < 1.0/(double)(lambda+1)) calculated_repair.push_back(i);
            else cost += weights[i];
        }
        for (int i = 0; i < M; ++i) {
            double yv = cplex.getValue(y[i]);
            if (yv >= 1.0/(double)(lambda+1)) cost += vio_weights[i];
        }
    } catch (IloException& e) {
        std::cerr << "CPLEX Exception: " << e << std::endl;
        env.end();
        throw;
    }
    env.end();
    return cost;
}

// ExactMW3SCRepair的CPLEX版本
// 与MW3SCRepair_Cplex类似，但变量为整数

double DBInst::ExactMW3SCRepair() {
    double cost = 0.0;
    if (tuples.size() == 0 || violations.size() == 0) return 0;
    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int M = violations.size();
        IloNumVarArray x(env, N, 0, 1, ILOINT);
        IloNumVarArray y(env, M, 0, 1, ILOINT);
        IloExpr obj(env);
        for (int i = 0; i < N; ++i) obj += weights[i] * x[i];
        for (int i = 0; i < M; ++i) obj += vio_weights[i] * y[i];
        model.add(IloMinimize(env, obj));
        obj.end();
        for (int i = 0; i < M; ++i) {
            IloExpr sum(env);
            for (int tid : violations[i]) sum += x[tid];
            sum -= y[i];
            model.add(sum >= 1);
            sum.end();
        }
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads,32);
        cplex.setWarning(env.getNullStream());
        bool ok = cplex.solve();
        if (!ok) throw std::runtime_error("CPLEX failed to solve ExactMW3SCRepair ILP");
        calculated_repair.clear();
        for (int i = 0; i < N; ++i) {
            double xv = cplex.getValue(x[i]);
            if (xv < 0.5) calculated_repair.push_back(i);
            else cost += weights[i];
        }
        for (int i = 0; i < M; ++i) {
            double yv = cplex.getValue(y[i]);
            if (yv > 0.5) cost += vio_weights[i];
        }
    } catch (IloException& e) {
        std::cerr << "CPLEX Exception: " << e << std::endl;
        env.end();
        throw;
    }
    env.end();
    return cost;
}

double DBInst::CalcRepair() {
    printf("begin calc repair\n");
    fflush(stdout);
    double cost = 0.0;
    if (tuples.size() == 0 || violations.size() == 0) return 0;
    IloEnv env;
    try {
        IloModel model(env);
        int N = tuples.size();
        int M = violations.size();
        // 变量
        IloNumVarArray x(env, N, 0, 1, ILOFLOAT); // x_t
        IloNumVarArray y(env, M, 0, 1, ILOFLOAT); // y_V
        // 目标函数
        IloExpr obj(env);
        for (int i = 0; i < N; ++i) obj += weights[i] * x[i];
        for (int i = 0; i < M; ++i) obj += vio_weights[i] * y[i];
        model.add(IloMinimize(env, obj));
        obj.end();
        // 约束
        for (int i = 0; i < M; ++i) {
            IloExpr sum(env);
            for (int tid : violations[i]) sum += x[tid];
            if (vio_weights[i] == _WEIGHT_INF) {
                model.add(sum >= 1);
            } else {
                model.add(sum + y[i] >= 1);
            }
            sum.end();
        }
        // 求解
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads,32);
        // cplex.setOut(env.getNullStream());
        // cplex.setWarning(env.getNullStream());
        bool ok = cplex.solve();
        if (!ok) throw std::runtime_error("CPLEX failed to solve LP");
        std::cout << "[CPLEX LP] obj: " << cplex.getObjValue() << std::endl;
        // 还原解
        std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne;
        calculated_repair.clear();
        for (int i = 0; i < N; ++i) {
            double xv = cplex.getValue(x[i]);
            // std::cout<<xv<<" ";
            if (xv < 1.0/(double)(lambda+1)) calculated_repair.push_back(i);
            else {
                tupleIndexRoundedToOne.insert(i);
                cost += weights[i];
            }
        }
        // std::cout<<'\n';
        for (int i = 0; i < M; ++i) {
            if (lambda == 2) {
                int x_sum = 0;
                for (int tup_index : violations[i]) {
                    if (tupleIndexRoundedToOne.count(tup_index)) x_sum++;
                }
                if (x_sum == 0) violationIndexRoundedToOne.insert(i);
            } else {
                double yv = cplex.getValue(y[i]);
                if (yv > 1.0/(double)(lambda+1)) violationIndexRoundedToOne.insert(i);
            }
        }
        for (int i : violationIndexRoundedToOne) cost += vio_weights[i];
        assert(CheckRoundedSolutionValid(tupleIndexRoundedToOne, violationIndexRoundedToOne));
    } catch (IloException& e) {
        std::cerr << "CPLEX Exception: " << e << std::endl;
        env.end();
        throw;
    }
    env.end();
    return cost;
}

double DBInst::GreedyRepair(){

    //Linear-time greedy set cover
    //pseudocode available at https://zhuanlan.zhihu.com/p/682124212
    //calculate on indexes. 0 ~ |I|-1: tuple, |I| ~ |I|+|V|-1: violation
    //Elements: (DC, t1, t2, ...)
    // contained in DC, t1, t2, ...
    
    printf("begin greedy repair\n");
    fflush(stdout);

    // 动态bitset实现（用vector<bool>替代bitset）
    double cost = 0.0;
    if (violations.size() == 0 || tuples.size() == 0) return 0;

    printf("begin calc size\n");
    int n_tuple = tuples.size();
    int n_vio = violations.size();
    int n_set = n_tuple + n_vio;
    printf("ntuple=%d, nvio=%d, nset=%d, end calc size\n",n_tuple, n_vio, n_set);

    // 用unordered_set稀疏存储集合-violation关系
    printf("begin init vec\n");
    std::vector<std::unordered_set<int>> set_to_vios(n_set);
    std::vector<double> set_weight(n_set);
    std::vector<int> set_cover_cnt(n_set); // 当前还能覆盖多少 violation
    std::vector<std::vector<int>> vio_to_sets(n_vio); // 每个 violation 被哪些集合覆盖
    printf("end init vec\n");

    // tuple集合
    printf("begin tuple sets (fast)\n");
    for (int j = 0; j < n_vio; ++j) {
        for (int tid : violations[j]) {
            set_to_vios[tid].insert(j);
            vio_to_sets[j].push_back(tid);
        }
    }
    for (int i = 0; i < n_tuple; ++i) {
        set_weight[i] = weights[i];
        set_cover_cnt[i] = set_to_vios[i].size();
    }
    // violation集合
    printf("begin vio sets\n");
    for (int i = 0; i < n_vio; ++i) {
        set_weight[n_tuple + i] = vio_weights[i];
        set_to_vios[n_tuple + i].insert(i);
        vio_to_sets[i].push_back(n_tuple + i);
        set_cover_cnt[n_tuple + i] = 1;
    }

    // 2. 初始化最大堆
    using HeapNode = std::pair<double, int>; // (profit, set_id)
    auto cmp = [](const HeapNode& a, const HeapNode& b) { return a.first < b.first; };
    std::priority_queue<HeapNode, std::vector<HeapNode>, decltype(cmp)> pq(cmp);
    printf("begin init max heap\n");
    for (int i = 0; i < n_set; ++i) {
        if (set_weight[i] < _WEIGHT_INF && set_cover_cnt[i] > 0) {
            pq.emplace(double(set_cover_cnt[i]), i);
        }
    }

    std::vector<bool> vio_covered(n_vio, false);
    std::vector<bool> set_used(n_set, false);
    int uncovered = n_vio;
    int uncovered_old = uncovered + 1;
    int loopcnt = 0;

    // 3. 主循环（unordered_set加速，保证每次循环uncovered减少）
    printf("begin while loop\n");
    while (uncovered > 0) {
        if(++loopcnt % 10000 == 0){
            printf("%d ",uncovered);
            fflush(stdout);
        }
        
        assert(uncovered < uncovered_old);
        uncovered_old = uncovered;
        HeapNode node;
        std::vector<int> newly_covered;
        do {
            if (pq.empty()) {
                throw std::runtime_error("No valid set left but uncovered > 0");
            }
            node = pq.top(); pq.pop();
            int set_id = node.second;
            if (set_used[set_id] || set_cover_cnt[set_id] == 0) continue;
            double real_profit = set_cover_cnt[set_id] > 0 ? double(set_cover_cnt[set_id]) : 0.0;
            if (fabs(real_profit - node.first) > 1e-8) {
                if (real_profit > 0.0)
                    pq.emplace(real_profit, set_id);
                continue;
            }
            // 计算新覆盖的 violation
            newly_covered.clear();
            for (int v : set_to_vios[set_id]) {
                if (!vio_covered[v]) {
                    newly_covered.push_back(v);
                }
            }
            if (newly_covered.empty()) continue;
            // 找到有效集合，退出do-while
            break;
        } while (true);
        int set_id = node.second;
        set_used[set_id] = true;
        for (int v : newly_covered) {
            vio_covered[v] = true;
            --uncovered;
        }
        // 只对受影响集合做一次批量更新，避免频繁入堆
        for (int v : newly_covered) {
            for (int sid : vio_to_sets[v]) {
                if (!set_used[sid] && set_cover_cnt[sid] > 0) {
                    // O(1)递减更新set_cover_cnt
                    --set_cover_cnt[sid];
                    double new_profit = double(set_cover_cnt[sid]) ;
                    if (set_cover_cnt[sid] > 0 && new_profit > 0.0)
                        pq.emplace(new_profit, sid);
                }
            }
        }
    }

    // 4. 还原解
    std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne;
    calculated_repair.clear();
    for (int i = 0; i < n_tuple; ++i) {
        if (set_used[i]) {
            cost += weights[i];
            tupleIndexRoundedToOne.insert(i);
        } else {
            calculated_repair.push_back(i);
        }
    }
    for (int i = 0; i < n_vio; ++i) {
        if (set_used[n_tuple + i]) {
            cost += vio_weights[i];
            violationIndexRoundedToOne.insert(i);
        }
    }
    assert(CheckRoundedSolutionValid(tupleIndexRoundedToOne, violationIndexRoundedToOne));
    return cost;
}

// LP修复+后处理：先用LP松弛解，记录所有x[i]∈[0.5,1)的元组（RoundedFromHalf），
// 然后贪心地尝试去掉这些元组，只要合法且cost下降就去掉，直到不能再优化。
double DBInst::CalcRepair_Postproc() {
    printf("begin calc repair with pp\n");
    fflush(stdout);
    double cost = 0.0;
    if (tuples.size() == 0 || violations.size() == 0) return 0;
    IloEnv env;
    std::set<int> tupleIndexRoundedToOne, violationIndexRoundedToOne;
    std::vector<int> RoundedFromHalf; // 记录LP解x[i]在[0.5,1)的元组
    try {
        IloModel model(env);
        int N = tuples.size();
        int M = violations.size();
        IloNumVarArray x(env, N, 0, 1, ILOFLOAT);
        IloNumVarArray y(env, M, 0, 1, ILOFLOAT);
        IloExpr obj(env);
        for (int i = 0; i < N; ++i) obj += weights[i] * x[i];
        for (int i = 0; i < M; ++i) obj += vio_weights[i] * y[i];
        model.add(IloMinimize(env, obj));
        obj.end();
        for (int i = 0; i < M; ++i) {
            IloExpr sum(env);
            for (int tid : violations[i]) sum += x[tid];
            if (vio_weights[i] == _WEIGHT_INF) {
                model.add(sum >= 1);
            } else {
                model.add(sum + y[i] >= 1);
            }
            sum.end();
        }
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads,32);
        cplex.setWarning(env.getNullStream());
        bool ok = cplex.solve();
        if (!ok) throw std::runtime_error("CPLEX failed to solve LP");
        calculated_repair.clear();
        tupleIndexRoundedToOne.clear();
        violationIndexRoundedToOne.clear();
        for (int i = 0; i < N; ++i) {
            double xv = cplex.getValue(x[i]);
            if (xv < 1.0/(double)(lambda+1)) {
                calculated_repair.push_back(i);
            } else {
                tupleIndexRoundedToOne.insert(i);
                cost += weights[i];
                if (xv >= 1.0/(double)lambda && xv < 1.0) {
                    RoundedFromHalf.push_back(i);
                }
            }
        }
        for (int i = 0; i < M; ++i) {
            if (lambda == 2) {
                int x_sum = 0;
                for (int tup_index : violations[i]) {
                    if (tupleIndexRoundedToOne.count(tup_index)) x_sum++;
                }
                if (x_sum == 0) violationIndexRoundedToOne.insert(i);
            } else {
                double yv = cplex.getValue(y[i]);
                if (yv > 1.0/(double)(lambda+1)) violationIndexRoundedToOne.insert(i);
            }
        }
        for (int i : violationIndexRoundedToOne) cost += vio_weights[i];
        assert(CheckRoundedSolutionValid(tupleIndexRoundedToOne, violationIndexRoundedToOne));
    } catch (IloException& e) {
        std::cerr << "CPLEX Exception: " << e << std::endl;
        env.end();
        throw;
    }
    env.end();
    // 后处理：堆优化，O(log n) 选取可删元组
    // 1. 初始化
    std::vector<int> vio_covered_cnt(violations.size(), 0);
    std::vector<std::vector<int>> tuple_to_vios(tuples.size());
    for (int i = 0; i < (int)violations.size(); ++i) {
        for (int tid : violations[i]) {
            if (tupleIndexRoundedToOne.count(tid)) vio_covered_cnt[i]++;
            if (tid < (int)tuples.size()) tuple_to_vios[tid].push_back(i);
        }
    }
    std::cout<<"start post process-----\n";
    // 堆节点: (delta, 版本号, tid)
    using HeapNode = std::tuple<double, int, int>;
    auto cmp = [](const HeapNode& a, const HeapNode& b) { return std::get<0>(a) < std::get<0>(b); };
    std::priority_queue<HeapNode, std::vector<HeapNode>, decltype(cmp)> pq(cmp);
    std::unordered_map<int, int> tuple_version; // tid -> version
    // 计算初始delta并入堆
    for (int tid : RoundedFromHalf) {
        bool can_remove = true;
        for (int v : tuple_to_vios[tid]) {
            if (vio_covered_cnt[v] == 1) { can_remove = false; break; }
        }
        if (!can_remove) continue;
        double delta = weights[tid];
        for (int v : tuple_to_vios[tid]) {
            if (vio_covered_cnt[v] == 1) delta -= vio_weights[v];
        }
        tuple_version[tid] = 0;
        pq.emplace(delta, 0, tid);
    }
    // 记录当前RoundedFromHalf集合
    std::unordered_set<int> rounded_set(RoundedFromHalf.begin(), RoundedFromHalf.end());
    while (!pq.empty()) {
        auto [delta, ver, tid] = pq.top(); pq.pop();
        if (!rounded_set.count(tid)) continue; // 已被删
        if (tuple_version[tid] != ver) continue; // 版本过期
        // 检查可删性
        bool can_remove = true;
        for (int v : tuple_to_vios[tid]) {
            if (vio_covered_cnt[v] == 1) { can_remove = false; break; }
        }
        if (!can_remove) continue;
        if (delta <= 1e-8) break; // 没有可删且cost下降的元组
        // 执行删除
        tupleIndexRoundedToOne.erase(tid);
        rounded_set.erase(tid);
        for (int v : tuple_to_vios[tid]) vio_covered_cnt[v]--;
        cost -= delta;
        // 增量更新受影响元组的delta
        for (int v : tuple_to_vios[tid]) {
            if (vio_covered_cnt[v] == 1) {
                // 只有当某violation的covered_cnt变为1时，相关元组的delta才会变
                for (int t2 : violations[v]) {
                    if (t2 == tid) continue;
                    if (!rounded_set.count(t2)) continue;
                    tuple_version[t2]++;
                    // 重新计算delta
                    bool can2 = true;
                    for (int v2 : tuple_to_vios[t2]) {
                        if (vio_covered_cnt[v2] == 1) { can2 = false; break; }
                    }
                    if (!can2) continue;
                    double delta2 = weights[t2];
                    for (int v2 : tuple_to_vios[t2]) {
                        if (vio_covered_cnt[v2] == 1) delta2 -= vio_weights[v2];
                    }
                    pq.emplace(delta2, tuple_version[t2], t2);
                }
            }
        }
    }
    assert(CheckRoundedSolutionValid(tupleIndexRoundedToOne, violationIndexRoundedToOne));
    // 重新生成calculated_repair
    calculated_repair.clear();
    for (int i = 0; i < (int)tuples.size(); ++i) {
        if (!tupleIndexRoundedToOne.count(i)) calculated_repair.push_back(i);
    }
    std::cout << "[CalcRepair_Postproc] Final cost: " << cost << ", repair size: " << calculated_repair.size() << std::endl;
    
    return cost;
}


double DBInst::SumOfTupleWeights(){
    double w = 0;
    for(int i=0; i<tuples.size(); i++){
        w += weights[i];
    }
    return w;
}

// DBInst成员函数：保存violations和vio_weights到文件
void DBInst::SaveViolationsAndWeights(const std::string& viofile, const std::string& weightfile) const {
    std::ofstream ofs_vio(viofile);
    for(const auto& vio : violations) {
        bool first = true;
        for(const auto& idx : vio) {
            if(!first) ofs_vio << ",";
            ofs_vio << idx;
            first = false;
        }
        ofs_vio << "\n";
    }
    ofs_vio.close();
    std::ofstream ofs_weight(weightfile);
    for(const auto& w : vio_weights) {
        ofs_weight << w << "\n";
    }
    ofs_weight.close();
}

// DBInst成员函数：从文件读取violations和vio_weights
void DBInst::LoadViolationsAndWeights(const std::string& viofile, const std::string& weightfile) {
    violations.clear();
    vio_weights.clear();
    std::ifstream ifs_vio(viofile);
    std::string line;
    while(std::getline(ifs_vio, line)) {
        std::set<int> vio;
        std::stringstream ss(line);
        std::string item;
        while(std::getline(ss, item, ',')) {
            if(!item.empty()) vio.insert(std::stoi(item));
        }
        violations.push_back(vio);
    }
    ifs_vio.close();
    std::ifstream ifs_weight(weightfile);
    double w;
    while(ifs_weight >> w) {
        vio_weights.push_back(w);
    }
    ifs_weight.close();
}