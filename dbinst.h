/*
    dbinst.h

    Created by baka on 3/27

*/

#include <string>
#include <vector>
#include <map>
#include <limits>
#include <set>
#include "sqlite3.h"
#include <ilcplex/ilocplex.h>

#ifndef ___SOFTREP__DB__
#define ___SOFTREP__DB__

//Relaitonal DB

typedef std::vector<std::string> tuple;
#define _WEIGHT_INF 1.79769e+308

class DBInst{

public:
    DBInst(std::string dir_path,std::string dataset_name, std::string constr_name = "constr");
    ~DBInst();

    // return value: cost, -1 stands for error.
    double CalcRepair();
    double CalcRepair_Postproc();

    double ExactRepair();

    //MW3SC APPROX
    double MW3SCRepair();
    double ExactMW3SCRepair();

    //SC_GREEDY
    double GreedyRepair();

    double SumOfTupleWeights();

    std::vector<tuple> OutputRepair();
    void output();

    void LoadViolationsAndWeights(const std::string& viofile, const std::string& weightfile);
    void SaveViolationsAndWeights(const std::string& viofile, const std::string& weightfile) const;

    // 检查给定的tupleIndexRoundedToOne和violationIndexRoundedToOne是否为合法解
    bool inline CheckRoundedSolutionValid(const std::set<int>& tupleIndexRoundedToOne, const std::set<int>& violationIndexRoundedToOne) const {
        // 对每个violation检查是否被覆盖或豁免
        for (size_t i = 0; i < violations.size(); ++i) {
            bool covered = false;
            // 只要有一个tuple被置为1且属于该violation，则视为被覆盖
            for (int tup : violations[i]) {
                if (tupleIndexRoundedToOne.count(tup)) {
                    covered = true;
                    break;
                }
            }
            // 如果没有被覆盖，则必须被豁免
            if (!covered && !violationIndexRoundedToOne.count(i)) {
                return false;
            }
        }
        return true;
    }

    double getCalcRatio(){
        // return (lambda==2) ? 2 : (lambda+1);
        return 2538.33/1973.26;
    }

    double getCalcPPRatio(){
        return 2010.25/1973.26;
    }

    double getGreedyRatio(){
        // return log(tuples.size() + violations.size());
        return 12714.6/1973.26;
    }

protected:
    sqlite3* db = NULL;

    std::vector<std::string> schema;
    std::map<std::string,int> schema_map;
    int schema_len = 0; //length of schema
    std::vector<tuple> tuples;
    std::vector<double> weights;

    int lambda = 2;

    //storing ids of the tuple
    std::vector<int> calculated_repair;
    
    std::string _dir_path, _dataset_name, _constr_name;
    std::vector<std::string> constraints;
    std::vector<double> constr_weights;
    std::string createtable_sql;

    //Many tuples consist of a violation
    // store indexes only
    std::vector<std::set<int>> violations;
    std::vector<double> vio_weights;

    void __calc__violations();
    void CalcViolations();
};

#endif