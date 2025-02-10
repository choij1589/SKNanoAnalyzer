#include "SystematicHelper.h"
#include <yaml-cpp/yaml.h>

SystematicHelper::SystematicHelper(std::string yaml_path, TString sample)
{
    variation_prefix = {
        {MyCorrection::variation::nom, ""},
        {MyCorrection::variation::up, "_Up"},
        {MyCorrection::variation::down, "_Down"}};

    YAML::Node config = YAML::LoadFile(yaml_path);
    for (const auto &node : config["systematics"])
    {
        std::string syst_name = node["syst"].as<std::string>();
        SystematicHelper::SYST syst;
        syst.syst = syst_name;
        if (node["source"].IsDefined())
        {
            syst.source = node["source"].as<std::string>();
        }
        if (node["evtLoopAgain"].IsDefined())
        {
            syst.evtLoopAgain = node["evtLoopAgain"].as<bool>();
        }
        if (node["oneSided"].IsDefined())
        {
            syst.oneSided = node["oneSided"].as<bool>();
        }
        if (node["hasDedicatedSample"].IsDefined())
        {
            syst.hasDedicatedSample = node["hasDedicatedSample"].as<bool>();
        }
        if (node["dedicatedSampleKey_up"].IsDefined())
        {
            syst.dedicatedSampleKey_up = node["dedicatedSampleKey_up"].as<std::string>();
        }
        if (node["dedicatedSampleKey_down"].IsDefined())
        {
            syst.dedicatedSampleKey_down = node["dedicatedSampleKey_down"].as<std::string>();
        }
        if (node["target"].IsDefined())
        {
            syst.target = node["target"].as<std::string>();
        }
        else
        {
            syst.target = syst_name;
        }
        systematics.push_back(syst);
    }

    for (const auto &node : config["correlations"])
    {
        CORRELATION correlation;
        std::string name = node["correlation"].as<std::string>();
        correlation.name = name;
        std::string rep_name;
        for (unsigned int i = 0; i < node["sources"].size(); i++)
        {
            if (i == 0)
            {
                rep_name = node["sources"][i].as<std::string>();
                correlation.rep_name = rep_name;
            }
            else
            {
                correlation.child_syst_names.insert(node["sources"][i].as<std::string>());
            }
            std::string source = node["sources"][i].as<std::string>();
            correlation.sources.push_back(findSystematic(source));
        }
        correlations[name] = correlation;
    }

    make_map_dedicatedSample();
    this->sample = sample.Data();
    isDedicatedSample = IsDedicatedSample();
    checkBadSystematics();
    make_Iter_obj_EvtLoopAgain();
    current_Iter_obj.clone(systematics_evtLoopAgain[0]); // Central
    std::cout << "[SystematicHelper::SystematicHelper] SystematicHelper is created for " <<  sample << std::endl;
    std::cout << "[SystematicHelper::SystematicHelper] Iterator has " << systematics_evtLoopAgain.size() << " elements" << std::endl;
    for (const auto &syst : systematics_evtLoopAgain)
    {
        std::cout << "[SystematicHelper::SystematicHelper] " << syst.iter_name << " " << syst.syst_name << " " << variation_prefix[syst.variation] << std::endl;
    }
    std::cout << "[SystematicHelper::SystematicHelper] All systematics are: " << std::endl;
    for (const auto &syst : systematics)
    {
        std::cout << "[SystematicHelper::SystematicHelper] " << syst.syst << " " << syst.source << " " << syst.evtLoopAgain << " " << syst.oneSided << " " << syst.hasDedicatedSample << " " << syst.dedicatedSampleKey_up << " " << syst.dedicatedSampleKey_down << std::endl;
    }
    std::cout << "[SystematicHelper::SystematicHelper] isDedicatedSample: " << isDedicatedSample << std::endl;
    std::cout << "[SystematicHelper::SystematicHelper] Current systematic: " << current_Iter_obj.iter_name << std::endl;
    std::cout << "[SystematicHelper::SystematicHelper] Current systematic source: " << current_Iter_obj.syst_name << std::endl;
    std::cout << "[SystematicHelper::SystematicHelper] Current systematic variation: " << variation_prefix[current_Iter_obj.variation] << std::endl;
}

SystematicHelper::~SystematicHelper() {}

SystematicHelper::SYST *SystematicHelper::findSystematic(std::string syst_name)
{
    for (auto &syst : systematics)
    {
        if (syst.syst == syst_name)
        {
            return &syst;
        }
    }
    return nullptr;
}

void SystematicHelper::checkBadSystematics()
{
    bool badSysts = false;
    for (const auto &syst : systematics)
    {
        // check 1. if hasDedicatedSample is true, then dedicatedSampleKey_up and dedicatedSampleKey_down must be defined
        if (syst.hasDedicatedSample)
        {
            if (syst.dedicatedSampleKey_up == "" || syst.dedicatedSampleKey_down == "")
            {
                std::cerr << "Systematic " << syst.syst << " has hasDedicatedSample set to true, but dedicatedSampleKey_up and dedicatedSampleKey_down are not defined" << std::endl;
                badSysts = true;
            }
        }

        // check 2. if evtLoopAgain is true, hasDedicatedSample must be false
        if (syst.evtLoopAgain && syst.hasDedicatedSample)
        {
            std::cerr << "Systematic " << syst.syst << " has evtLoopAgain set to true, but hasDedicatedSample is also true" << std::endl;
            badSysts = true;
        }
    }

    for (const auto &correlation : correlations)
    {
        for (const auto &source : correlation.second.sources)
        {
            if (source == nullptr)
            {
                std::cerr << "Correlation " << correlation.first << " has a source that is not defined" << std::endl;
                badSysts = true;
            }
        }
    }

    if (badSysts)
    {
        throw std::runtime_error("Bad systematics found");
        exit(1);
    }
}

void SystematicHelper::make_Iter_obj_EvtLoopAgain()
{

    SystematicHelper::Iter_obj obj_central;
    obj_central.iter_name = central_name;
    obj_central.syst_name = central_name;
    obj_central.syst_source = "total";
    obj_central.variation = MyCorrection::variation::nom;
    systematics_evtLoopAgain.push_back(obj_central);

    // if systematic sample is provided, only loop over central
    if (isDedicatedSample)
        return;

    for (const auto &syst : systematics)
    {
        if (syst.evtLoopAgain)
        {
            SystematicHelper::Iter_obj obj_up;
            SystematicHelper::Iter_obj obj_down;
            obj_up.iter_name = syst.syst + variation_prefix[MyCorrection::variation::up];
            obj_down.iter_name = syst.syst + variation_prefix[MyCorrection::variation::down];
            obj_up.syst_name = syst.syst;
            obj_up.syst_source = syst.source;
            obj_up.variation = MyCorrection::variation::up;
            obj_down.syst_name = syst.syst;
            obj_down.syst_source = syst.source;
            obj_down.variation = MyCorrection::variation::down;

            systematics_evtLoopAgain.push_back(obj_up);
            systematics_evtLoopAgain.push_back(obj_down);
        }
    }
}

void SystematicHelper::make_map_dedicatedSample()
{
    for (const auto &syst : systematics)
    {
        if (syst.hasDedicatedSample)
        {
            map_dedicatesamplekey_systname[syst.dedicatedSampleKey_up] = syst.syst + variation_prefix[MyCorrection::variation::up];
            map_dedicatesamplekey_systname[syst.dedicatedSampleKey_down] = syst.syst + variation_prefix[MyCorrection::variation::down];
        }
    }
}

void SystematicHelper::assignWeightFunctionMap(const unordered_map<std::string, std::variant<std::function<float(MyCorrection::variation, TString)>, std::function<float()>>> &weight_function_map)
{
    std::vector<std::string> syst_no_weight_function;
    // clear weight functions
    weight_functions.clear();
    weight_functions_onesided.clear();
    std::unordered_set<std::string> all_targets;
    for (const auto &syst : systematics)
        all_targets.insert(syst.target);
    for (const auto &syst : systematics)
    {
        if (syst.hasDedicatedSample || syst.evtLoopAgain)
            continue;

        auto it = weight_function_map.find(syst.target);
        if (it == weight_function_map.end())
        {
            syst_no_weight_function.push_back(syst.syst);
            continue;
        }

        const auto &weight_variant = it->second;

        auto assignOneSided = [&]()
        {
            if (weight_variant.index() != 1)
            {
                std::cerr << "Weight function for " << syst.syst << " is not a one-sided weight function" << std::endl;
                exit(1);
            }
            weight_functions_onesided[syst.target] = std::get<std::function<float()>>(weight_variant);
        };

        auto assignTwoSided = [&]()
        {
            if (weight_variant.index() != 0)
            {
                std::cerr << "Weight function for " << syst.syst << " is not a two-sided weight function" << std::endl;
                exit(1);
            }
            weight_functions[syst.target] = std::get<std::function<float(MyCorrection::variation, TString)>>(weight_variant);
        };

        syst.oneSided ? assignOneSided() : assignTwoSided();
    }
    if (syst_no_weight_function.size() > 0)
    {
        std::cerr << "Weight function for " << std::endl;
        for (const auto &syst : syst_no_weight_function)
        {
            std::cerr << syst << std::endl;
        }
        std::cerr << "are not assigned" << std::endl;
        exit(1);
    }
    weight_functions_assigned = true;
}

std::unordered_map<std::string, float> SystematicHelper::calculateWeight()
{
    if(!weight_functions_assigned)
    {
        std::cerr << "Weight functions are not assigned" << std::endl;
        exit(1);
    }
    for (const auto &syst : systematics)
    {
        if (syst.hasDedicatedSample || syst.evtLoopAgain)
        {
            weight_map_nominal[syst.target] = 1.;
            weight_map_up[syst.target] = 1.;
            weight_map_down[syst.target] = 1.;
        }
        else if (syst.oneSided)
        {
            float weight = 1.;
            auto weight_function = weight_functions_onesided[syst.target];
            weight = weight_function();

            weight_map_nominal[syst.target] = weight;
            weight_map_up[syst.syst] = 1.;
            weight_map_down[syst.syst] = weight;
        }
        else
        {
            float weight_up = 1.;
            float weight_down = 1.;
            float weight_nominal = 1.;

            auto weight_function = weight_functions[syst.target];
            weight_nominal = weight_function(MyCorrection::variation::nom, "total");
            weight_up = weight_function(MyCorrection::variation::up, syst.source);
            weight_down = weight_function(MyCorrection::variation::down, syst.source);

            weight_map_nominal[syst.target] = weight_nominal;
            weight_map_up[syst.syst] = weight_up;
            weight_map_down[syst.syst] = weight_down;
        }
    }
    unordered_map<std::string, float> weights;
    if (current_Iter_obj.variation == MyCorrection::variation::nom)
    {
        weights = calculateWeight_central_case();
    }
    else
    {
        weights = calculateWeight_non_central_case();
    }

    return weights;
}

unordered_map<std::string, float> SystematicHelper::calculateWeight_central_case()
{
    std::vector<string> all_weight_systs;
    unordered_map<std::string, float> weights;
    float nominal_weight = 1.;
    unordered_set<std::string> all_weight_targets;

    for (const auto &syst : systematics)
        all_weight_targets.insert(syst.target);

    for (const auto &target : all_weight_targets)
        nominal_weight *= weight_map_nominal[target];

    for (const auto &syst : systematics)
    {
        all_weight_systs.push_back(syst.syst);
    }

    weights[central_name] = nominal_weight;

    for (const auto &correlation : correlations)
    {
        unordered_set<std::string> all_sources_name = correlation.second.child_syst_names;
        all_sources_name.insert(correlation.second.rep_name);
        float weight_up = nominal_weight;
        float weight_down = nominal_weight;
        // If current Iter_obj is Central and correlation include systematic that has dedicated sample of require evtLoopAgain, skip.
        // that correlation set will be calculated in that situation
        bool correlation_should_be_skipped = false;
        for (const auto &sources_in_table : all_sources_name)
        {
            if (findSystematic(sources_in_table)->hasDedicatedSample || findSystematic(sources_in_table)->evtLoopAgain)
            {
                correlation_should_be_skipped = true;
                break;
            }
        }
        
        for (const auto &sources_in_table : all_sources_name)
        {
            std::string this_target = findSystematic(sources_in_table)->target;
            weight_up = safe_divide(weight_up, weight_map_nominal[this_target]);
            weight_down = safe_divide(weight_down, weight_map_nominal[this_target]);
            weight_up *= weight_map_up[sources_in_table];
            weight_down *= weight_map_down[sources_in_table];
            auto it = find(all_weight_systs.begin(), all_weight_systs.end(), sources_in_table);
            if (it != all_weight_systs.end())
            {
                all_weight_systs.erase(it);
            }
        }
        if (correlation_should_be_skipped) continue;
        weights[correlation.second.rep_name + variation_prefix[MyCorrection::variation::up]] = weight_up;
        weights[correlation.second.rep_name + variation_prefix[MyCorrection::variation::down]] = weight_down;
    }

    // systematic that not in correlation table
    for (const auto &syst : all_weight_systs)
    {
        float weight_up = nominal_weight;
        float weight_down = nominal_weight;
        if(findSystematic(syst)->hasDedicatedSample || findSystematic(syst)->evtLoopAgain)
        {
            continue;
        }
        std::string this_target = findSystematic(syst)->target;
        weight_up = safe_divide(weight_up, weight_map_nominal[this_target]);
        weight_down = safe_divide(weight_down, weight_map_nominal[this_target]);
        weight_up *= weight_map_up[syst];
        weight_down *= weight_map_down[syst];
        weights[syst + variation_prefix[MyCorrection::variation::up]] = weight_up;
        weights[syst + variation_prefix[MyCorrection::variation::down]] = weight_down;
    }
    if(isDedicatedSample)
    {
        for (const auto &dedicatedSample : map_dedicatesamplekey_systname)
        {
            std::string this_key = dedicatedSample.first;
            // Check if sample is found as a substring of this_sample
            if (sample.find(this_key) != std::string::npos)
            {
                float central_weight = weights[central_name];
                weights.clear();
                weights[dedicatedSample.second] = central_weight; // Use value from "Central"
                return weights;
            }
        }
        throw std::runtime_error("[SystematicHelper::calculateWeight_central_case] Dedicated sample not found");
    }
    return weights;
}

unordered_map<std::string, float> SystematicHelper::calculateWeight_non_central_case()
{
    // first check current Iter_obj is in correlation table
    bool Iter_obj_in_correlation = false;
    CORRELATION this_correlation;
    unordered_map<std::string, float> weights;
    unordered_set<std::string> all_weight_targets;
    for (const auto &syst : systematics)
        all_weight_targets.insert(syst.target);
    for (const auto &correlation : correlations)
    {
        if (correlation.second.rep_name == current_Iter_obj.syst_name)
        {
            Iter_obj_in_correlation = true;
            this_correlation = correlation.second;
            break;
        }
        for (const auto &child_syst : correlation.second.child_syst_names)
        {
            if (child_syst == current_Iter_obj.syst_name)
            {
                Iter_obj_in_correlation = true;
                this_correlation = correlation.second;
                break;
            }
        }
    }

    weights[current_Iter_obj.iter_name] = 1.;
    for (const auto &target : all_weight_targets)
    {
        weights[current_Iter_obj.iter_name] *= weight_map_nominal[target];
    }
    if (!Iter_obj_in_correlation)
        return weights;
    else
    {

        switch (current_Iter_obj.variation)
        {
        case MyCorrection::variation::up:
            weights[current_Iter_obj.iter_name]  = safe_divide(weights[current_Iter_obj.iter_name], weight_map_nominal[findSystematic(this_correlation.rep_name)->target]);
            weights[current_Iter_obj.iter_name] *= weight_map_up[findSystematic(this_correlation.rep_name)->syst];
            for (const auto &syst : this_correlation.child_syst_names)
            {
                weights[current_Iter_obj.iter_name] = safe_divide(weights[current_Iter_obj.iter_name] ,weight_map_nominal[findSystematic(syst)->target]);
                weights[current_Iter_obj.iter_name] *= weight_map_up[syst];
            }
            break;
        case MyCorrection::variation::down:
            weights[current_Iter_obj.iter_name] = safe_divide(weights[current_Iter_obj.iter_name], weight_map_nominal[findSystematic(this_correlation.rep_name)->target]);
            weights[current_Iter_obj.iter_name] *= weight_map_down[findSystematic(this_correlation.rep_name)->syst];
            for (const auto &syst : this_correlation.child_syst_names)
            {
                weights[current_Iter_obj.iter_name] = safe_divide(weights[current_Iter_obj.iter_name], weight_map_nominal[findSystematic(syst)->target]);
                weights[current_Iter_obj.iter_name] *= weight_map_down[syst];
            }
            break;
        default:
            break;
        }
    }
    return weights;
}

bool SystematicHelper::IsDedicatedSample()
{ 
    for (const auto &dedicatedSample : map_dedicatesamplekey_systname)
    {
        if (sample.find(dedicatedSample.first) != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

std::vector<std::string> SystematicHelper::get_targets_from_name(const std::string &syst_name){
    // remove variation prefix
    std::vector<std::string> targets;
    if(syst_name == central_name){
        targets.push_back(central_name);
        return targets;
    }
    std::string syst_name_no_variation = syst_name;
    for (const auto &prefix : variation_prefix)
    {
        size_t found = syst_name.find(prefix.second);
        if (found != std::string::npos)
        {
            syst_name_no_variation = syst_name.substr(0, found);
            break;
        }
    }
    SYST *syst = findSystematic(syst_name_no_variation);
    
    //check if syst is in correlation table
    for (const auto &correlation : correlations)
    {
        if (correlation.second.rep_name == syst_name_no_variation)
        {
            targets.push_back(findSystematic(correlation.second.rep_name)->target);
            for (const auto &child_syst : correlation.second.child_syst_names)
            {
                targets.push_back(findSystematic(child_syst)->target);
            }
            return targets;
        }
    }
    targets.push_back(syst->target);
    return targets;
}

std::vector<std::string> SystematicHelper::get_sources_from_name(const std::string &syst_name){
    // remove variation prefix
    std::vector<std::string> sources;
    if (syst_name == central_name){
        sources.push_back(central_name);
        return sources;
    }
    std::string syst_name_no_variation = syst_name;
    for (const auto &prefix : variation_prefix)
    {
        size_t found = syst_name.find(prefix.second);
        if (found != std::string::npos)
        {
            syst_name_no_variation = syst_name.substr(0, found);
            break;
        }
    }

    SYST *syst = findSystematic(syst_name_no_variation);
    
    //check if syst is in correlation table
    for (const auto &correlation : correlations)
    {
        if (correlation.second.rep_name == syst_name_no_variation)
        {
            sources.push_back(findSystematic(correlation.second.rep_name)->source);
            for (const auto &child_syst : correlation.second.child_syst_names)
            {
                sources.push_back(findSystematic(child_syst)->source);
            }
            return sources;
        }
    }
    sources.push_back(syst->source);
    return sources;
}

MyCorrection::variation SystematicHelper::get_variation_from_name(const std::string &syst_name){
    //find which prefix is in the syst_name
    if (syst_name == central_name){
        return MyCorrection::variation::nom;
    }
    for (const auto &prefix : variation_prefix)
    {
        size_t found = syst_name.find(prefix.second);
        if (found != std::string::npos)
        {
            return prefix.first;
        }
    }
    throw std::runtime_error("[SystematicHelper::get_variation_from_name] weird syst_name");
}
