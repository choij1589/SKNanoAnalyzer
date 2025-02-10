#ifndef SystematicHelper_h
#define SystematicHelper_h
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include "ROOT/RVec.hxx"
#include "MyCorrection.h"
#include <variant>

using namespace ROOT::VecOps;
using namespace std;

typedef std::unordered_map<std::string, std::variant<std::function<float(MyCorrection::variation, std::string)>, std::function<float()>>> weight_function_map;

class SystematicHelper
{
public:
    struct SYST
    {
        std::string syst;
        std::string source = "total";
        std::string target = "";
        bool evtLoopAgain = false;
        bool oneSided = false;
        bool hasDedicatedSample = false;
        std::string dedicatedSampleKey_up = "";
        std::string dedicatedSampleKey_down = "";
    };

    struct CORRELATION
    {
        std::string name;
        std::string rep_name;                        // the first one in the sources
        unordered_set<std::string> child_syst_names; // the rest of the sources
        RVec<SYST *> sources;
    };

    struct Iter_obj
    {
        std::string iter_name;
        std::string syst_name;
        std::string syst_source;
        MyCorrection::variation variation;
    };


    struct Internal_Iter_obj
    {
        std::string iter_name;
        std::string syst_name;
        std::string syst_source;
        MyCorrection::variation variation;
        void clone(Iter_obj &obj)
        {
            iter_name = obj.iter_name;
            syst_name = obj.syst_name;
            syst_source = obj.syst_source;
            variation = obj.variation;
        }
    };

    SystematicHelper(std::string yaml_path, TString sample);
    ~SystematicHelper();
    SYST* findSystematic(std::string syst_name);
    void assignWeightFunctionMap(const unordered_map < std::string, std::variant<std::function<float(MyCorrection::variation, TString)>, std::function<float()>>> &weight_functions);
    class EventLoopIterator
    {
    public:
        using Iterator = std::vector<Iter_obj>::iterator;

        EventLoopIterator(Iterator it, SystematicHelper* parent) : iter(it), parent(parent) {
            if(iter != parent->systematics_evtLoopAgain.end()){
                parent->current_Iter_obj.clone(*iter);
            }
        }

        EventLoopIterator& operator++()
        {
            ++iter;
            if(iter != parent->systematics_evtLoopAgain.end()){
                parent->current_Iter_obj.clone(*iter);
            }
            return *this;
        }

        bool operator!=(const EventLoopIterator &other) const
        {
            return iter != other.iter;
        }

        const Iter_obj& operator*() const
        {
            return *iter;
        }

    private:
        Iterator iter;
        SystematicHelper* parent;
 
    };
    

    EventLoopIterator begin()
    {
        return EventLoopIterator(systematics_evtLoopAgain.begin(), this);
    }

    // Method to return the end iterator
    EventLoopIterator end()
    {
        return EventLoopIterator(systematics_evtLoopAgain.end(), this);
    }

    inline std::string getCurrentSysName() const { return current_Iter_obj.iter_name; }
    inline std::string getCurrentIterSysTarget() const { return current_Iter_obj.syst_name; }
    inline std::string getCurrentIterSysSource() const { return current_Iter_obj.syst_source; }
    inline MyCorrection::variation getCurrentIterVariation() const { return current_Iter_obj.variation; }

    std::vector<std::string> get_targets_from_name(const std::string &syst_name);
    std::vector<std::string> get_sources_from_name(const std::string &syst_name);
    MyCorrection::variation get_variation_from_name(const std::string &syst_name);
    std::unordered_map<std::string, float> calculateWeight();
    inline float safe_divide(float numerator, float denominator)
    {
        if (abs(denominator) < 1e-8)
        {
            return numerator;
        }
        return numerator / denominator;
    }

private:
    std::string central_name = "Central";
    void checkBadSystematics();
    void make_Iter_obj_EvtLoopAgain();
    void make_map_dedicatedSample();
    bool IsDedicatedSample();
    
    std::unordered_map<std::string, float> calculateWeight_central_case();
    std::unordered_map<std::string, float> calculateWeight_non_central_case();

    std::vector<SYST> systematics;
    std::unordered_map<std::string, CORRELATION> correlations;
    std::unordered_map<std::string, std::function<float(MyCorrection::variation, TString)>> weight_functions;
    std::unordered_map<std::string, std::function<float()>> weight_functions_onesided;

    std::vector<Iter_obj> systematics_evtLoopAgain;
    unordered_map<std::string, std::string> map_dedicatesamplekey_systname;
    
    std::unordered_map<std::string, float> weight_map_nominal;
    std::unordered_map<std::string, float> weight_map_up;
    std::unordered_map<std::string, float> weight_map_down;
    std::unordered_map<MyCorrection::variation, std::string> variation_prefix; 

    bool isDedicatedSample;
    bool weight_functions_assigned;
    std::string sample;
    Internal_Iter_obj current_Iter_obj;
};

#endif
