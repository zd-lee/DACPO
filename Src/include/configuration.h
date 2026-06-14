/**
 * @file configuration.h
 * @author lizd
 * @version 0.1
 * @date 2024-07-15
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once
#include <boost/property_tree/ptree.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <map>

enum class IPSR_TYPE {
    CLASSIC_IPSR,
    BITREE_IPSR,
    GRAPH_IPSR,
    TEST_GRAPH
};

typedef void (*IPSR_ENTRANCE)(const std::string& input_name, const std::string& output_path, int iters, double pointweight, int depth, int k_neighbors);


void merge_json(nlohmann::json& base, const nlohmann::json& overrideObj);

/**
 * 
 * 
 */
class ConfigManager
{
public:
    static std::string config_floder;

    /**
     * @return nlohmann::json
     */
    static nlohmann::json get_common_config();

    static bool check_config_floder(std::string config_floder);

    static nlohmann::json get_config_in_config_floder(std::string filename, bool if_merge = false); 

};


/**
 */
class JsonConfig
{

private:
    nlohmann::json config;
public:
    JsonConfig (std::string filename);
    JsonConfig (nlohmann::json config):config(config){}
    ~JsonConfig(){}

    nlohmann::json get_config(){
        return config;
    }
};