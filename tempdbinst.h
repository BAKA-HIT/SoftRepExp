#include <map>
#include <string>
#include "dbinst.h"
#include <cassert>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cassert>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <random>

#ifndef ___SOFTREP__TDB__
#define ___SOFTREP__TDB__

class TempDBInst : public DBInst{

public:
    TempDBInst(std::string dir_path,std::string dataset_name, std::string constr_name = "constr");
    ~TempDBInst();

    void Normalize();
    // void Coalesce();

    double CalcRepair();
    double ExactRepair();
    double CalcRepairWOContConstr();
    double CalcRepair_Postproc();

    void setAlpha(double a){
        alpha = a;
        //reweight adj weights
        // for(auto adj_pair : _adj){
        //     _adj_weights[adj_pair.first] = (weights[adj_pair.first] + weights[adj_pair.second]) * alpha;
        // }
    }

    double WeightedBrokenAdjRatio() {

        std::set<int> tupleIndexRoundedToOne;
        std::set<int> repairs(calculated_repair.begin(), calculated_repair.end());

        for(int i=0; i<tuples.size();i++){
            if(repairs.find(i)==repairs.end())  tupleIndexRoundedToOne.insert(i);
        }

        // 统计所有adj约束的总权重
        double total_weight = 0.0;
        for (const auto& adj_pair : _adj) {
            int adj_key = adj_pair.first;
            total_weight += _adj_weights[adj_key];
        }

        // 统计违反的adj约束的总权重
        double broken_weight = 0.0;
        for (const auto& adj_pair : _adj) {
            int u = adj_pair.first;
            int v = adj_pair.second;
            // 如果u或v在tupleIndexRoundedToOne中，说明违反
            if (tupleIndexRoundedToOne.count(u) || tupleIndexRoundedToOne.count(v)) {
                broken_weight += _adj_weights[u];
            }
        }

        if (total_weight == 0.0) return 0.0;
        return broken_weight / total_weight;
    }

    // 区间结构体
    typedef std::pair<double, double> Interval;

    // 计算一组区间的并集测度（总长度）
    static double MeasureOfUnion(const std::vector<Interval>& intervals) {
        if (intervals.empty()) return 0.0;
        // 拷贝并排序
        std::vector<Interval> segs = intervals;
        std::sort(segs.begin(), segs.end());
        double total = 0.0;
        double cur_l = segs[0].first, cur_r = segs[0].second;
        for (size_t i = 1; i < segs.size(); ++i) {
            if (segs[i].first > cur_r) {
                total += cur_r - cur_l;
                cur_l = segs[i].first;
                cur_r = segs[i].second;
            } else {
                cur_r = std::max(cur_r, segs[i].second);
            }
        }
        total += cur_r - cur_l;
        return total;
    }

    double MeasureOfCoverage(){
        std::vector<Interval> intervals;
        for(auto i : calculated_repair){
            intervals.push_back(std::pair<double,double>(get_ts(i),get_te(i)));
        }

        return MeasureOfUnion(intervals) / interval_len_of_all_tuples;
    }

    // 封装cost和相关集合的更新
    double UpdateCostAndSets(std::set<int>& tupleIndexRoundedToOne,
                             std::set<int>& violationIndexRoundedToOne,
                             std::set<int>& adjIndexRoundedToOne);
    double UpdateCostAndSetsIncremental(std::set<int>&, std::set<int>&, std::set<int>&, int);
    
    void PostProc(
    std::set<int>& tupleIndexRoundedToOne,
    std::set<int>& violationIndexRoundedToOne,
    std::set<int>& adjIndexRoundedToOne,
    double& cost)
    {
        bool improved = true;
        int improve_cnt = 0;
        while (improved && improve_cnt <= 100) {
            improve_cnt++;
            improved = false;
            for (auto it = tupleIndexRoundedToOne.begin(); it != tupleIndexRoundedToOne.end(); ) {
                int tid = *it;
                // 备份仅与tid相关的内容
                auto vio_bak = violationIndexRoundedToOne;
                auto adj_bak = adjIndexRoundedToOne;
                auto repair_bak = calculated_repair;
                double old_cost = cost;

                // 尝试删除tid
                tupleIndexRoundedToOne.erase(it++);
                double new_cost = UpdateCostAndSetsIncremental(tupleIndexRoundedToOne, violationIndexRoundedToOne, adjIndexRoundedToOne, tid);

                if (new_cost < old_cost && CheckRoundedSolutionValid(tupleIndexRoundedToOne, violationIndexRoundedToOne)) {
                    cost = new_cost;
                    improved = true;
                    break; // 每轮只删一个，重新遍历
                } else {
                    // 撤销
                    tupleIndexRoundedToOne.insert(tid);
                    violationIndexRoundedToOne = vio_bak;
                    adjIndexRoundedToOne = adj_bak;
                    calculated_repair = repair_bak;
                }
            }
        }
    }

protected:
    double alpha;

    //adj weight = _adj_weights[first_tuple_id]
    std::map<int,int> _adj;
    std::map<int, double> _adj_weights;

    std::set<double> timestamps;
    std::vector<double> vio_duration;
    std::vector<std::vector<int>> tuple_to_violations, tuple_to_adjs;

    double weighted_sum_of_broken_adj = 0;
    int interval_len_of_all_tuples = 0;

    inline double get_ts(int index){
        assert(schema_map.find("Ts") != schema_map.end());
        return std::stod(tuples[index][schema_map["Ts"]]);
    }

    inline void set_ts(int index, int ts){
        assert(schema_map.find("Ts") != schema_map.end());
        tuples[index][schema_map["Ts"]]=std::to_string(ts);
    }

    inline double get_te(int index){
        assert(schema_map.find("Te") != schema_map.end());
        return std::stod(tuples[index][schema_map["Te"]]);
    }

    inline void set_te(int index, int te){
        assert(schema_map.find("Te") != schema_map.end());
        tuples[index][schema_map["Te"]]=std::to_string(te);
    }

    inline void calc_vio_duration(){
        //vio_duration
        double max_ts = 0.0;
        double min_te = _WEIGHT_INF;
        vio_duration.clear();
        for(int i=0; i<violations.size(); i++){
            double max_ts = 0.0;
            double min_te = _WEIGHT_INF;
            for(int tup_index : violations[i]){
                max_ts = std::max(max_ts, get_ts(tup_index));
                min_te = std::min(min_te, get_te(tup_index));
            }
            // assert(max_ts<min_te);
            vio_duration.push_back(max_ts < min_te ? (min_te - max_ts) : 0);
            assert(i+1 == vio_duration.size());
        }
    }

};

#endif