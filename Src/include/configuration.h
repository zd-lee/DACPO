/**
 * @file configuration.h
 * @author lizd
 * @brief 用于项目的配置文件的读取 
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
 * @brief 用于管理配置文件。
 * 通过修改config_floder来修改配置文件夹的位置，静默模式下默认为default_config_floder
 * 
 * 通过get_config_in_config_floder来获取config_floder下的json文件
 * 
 * 当config_floder不等于default_config_floder时，在每次调用get_config_in_config_floder时
 * 会同时读取default_config_floder下的同名文件，并将两个文件合并：
 *      当两个文件中有相同的key时，以config_floder下的文件为准
 *      其他情况下，以default_config_floder下的文件为准
 */
class ConfigManager
{
public:
    // const static std::string default_config_floder; //默认的配置文件夹
    static std::string config_floder; // ipsr的一些基本配置

    /**
     * @brief 直接获取common.json的内容,其中包含了一些通用的配置
     * 由于这个函数可能需要频繁调用，因此内部有一个静态变量来保存第一次读取common.json的内容
     * @return nlohmann::json
     */
    static nlohmann::json get_common_config();

    static bool check_config_floder(std::string config_floder); //检查配置文件夹下是否有common.json

    // 获取配置文件夹下的json文件
    static nlohmann::json get_config_in_config_floder(std::string filename, bool if_merge = false); 

};


/**
 * @brief 用于读取json文件
 * 默认情况下，会读取default_config_floder下的common.json
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