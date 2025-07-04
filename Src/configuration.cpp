#include "configuration.h"
#include <tools.h>

JsonConfig::JsonConfig(std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cout << "Error: open file " << filename << " failed" << std::endl;
        file.close();
        return;
    }
    file >> config;
    file.close();
}


// const std::string ConfigManager::default_config_floder = "conf/default/";
std::string ConfigManager::config_floder = "";

nlohmann::json ConfigManager::get_common_config()
{
    static std::mutex mtx;
    bool got = false;
    
    static nlohmann::json j;
    mtx.lock();
    if (!got) {
        j = ConfigManager::get_config_in_config_floder("common.json");
        got = true;
    }
    mtx.unlock();
    return j;
}

bool ConfigManager::check_config_floder(std::string config_floder)
{
    std::string json_path = config_floder + "common.json";
    std::ifstream file(json_path);
    if (!file.is_open())
    {
        std::cout << "Error: open file " << json_path << " failed" << std::endl;
        file.close();
        return false;
    }
    file.close();
    return true;
}

nlohmann::json ConfigManager::get_config_in_config_floder(std::string filename,bool if_merge)
{
    JsonConfig tj(ConfigManager::config_floder + filename);
    if (!if_merge) {
        return tj.get_config();
    }
    // nlohmann::json j = tj.get_config();
    // JsonConfig default_tj(ConfigManager::default_config_floder + filename);
    // nlohmann::json default_j = default_tj.get_config();
    // merge_json(default_j, j);
    // return default_j;
    return tj.get_config();

}


void merge_json(nlohmann::json& base, const nlohmann::json& overrideObj)
{
    for (auto& element : overrideObj.items()) {
        if (element.value().is_object()) {
            // 如果要覆盖的对象是嵌套的JSON对象，则递归调用
            if (base.find(element.key()) != base.end() && base[element.key()].is_object()) {
                merge_json(base[element.key()], element.value());
            }
            else {
                if (base.find(element.key()) == base.end()) {
                    if (element.key().find("spconfig") == std::string::npos) {
                        printf("ERROR: key:%s\n not found in basic json!!!", element.key().c_str());
                        assert(false); 
                    }
                    else {
                        printf("WARNING: key:%s\n not found in basic json!!!", element.key().c_str());
                        base[element.key()] = element.value();
                    }
                }else{
                    base[element.key()] = element.value();
                }
            }
        }
        else {
            // 如果要覆盖的对象不是嵌套的JSON对象，则直接覆盖
            base[element.key()] = element.value();
        }
    }
}


