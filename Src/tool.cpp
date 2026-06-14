#include <tools.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>


lzd_tools::Timer::Timer(nlohmann::json& j,std::string func):_j(j),_func_name(func) {
    start = std::chrono::steady_clock::now();
}

lzd_tools::Timer::~Timer(){
    if(!_done){
        Done();
    }
}

void lzd_tools::Timer::Start()
{
    start = std::chrono::steady_clock::now();
}

void lzd_tools::Timer::Done()
{
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) / 1000.0;
    _j[_func_name + "_cost"] = duration.count();
    _done = true;
}

lzd_tools::AcumutableTimer::AcumutableTimer(std::string func_name):_func_name(func_name),_done(false){
    start = std::chrono::steady_clock::now();
}

void lzd_tools::AcumutableTimer::Done(){
    if (_done.get()) {
        return;
    }
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start) / 1000.0;
    _mtx.lock();
    _time_map[_func_name] += duration.count();
    _count_map[_func_name] += 1;
    _mtx.unlock();
    _done.set(true);
    
}

lzd_tools::AcumutableTimer::~AcumutableTimer(){
    if (!_done.get()) {
        Done();
    }
}

double lzd_tools::AcumutableTimer::get_time(std::string func_name){
    _mtx.lock();
    double res = _time_map[func_name];
    _mtx.unlock();
    return res;
}

nlohmann::json lzd_tools::AcumutableTimer::get_time_map(){
    nlohmann::json j;
    _mtx.lock();
    for(auto i:_time_map){
        j[i.first]["time_cost"] = i.second;
    }
    for (auto i : _count_map) {
        j[i.first]["call_count"] = i.second;
    }
    _mtx.unlock();
    return j;
}

std::map<std::string, double> lzd_tools::AcumutableTimer::_time_map;
std::map<std::string, int> lzd_tools::AcumutableTimer::_count_map;
std::mutex lzd_tools::AcumutableTimer::_mtx;




int change_args(std::vector<std::string>& args, std::string key, std::string value){
    for(int i=0;i<args.size()-1;i++){
        if(args[i] == key){
            args[i+1] = value;
            return i;
        }
    }
    printf("\n\n!!chane_args error : key %s not found!!\n\n",key);
    return -1;
}

int find_arg(std::vector<std::string>& args, std::string key){
    for(int i=0;i<args.size()-1;i++){
        if(args[i] == key){
            return i;
        }
    }
    printf("\n\n!!find_arg error : key %s not found!!\n\n",key);
    return -1;
}

std::vector<std::vector<int>> MapLable2Color(std::vector<int> flag) {
    int max = *std::max_element(flag.begin(), flag.end());
    int min = *std::min_element(flag.begin(), flag.end());
    if (min < 0) {
        std::cout << "Error: flag contains a negative value" << std::endl;
        return std::vector<std::vector<int>>();
    }
    if (max > 1e5) {
        std::cout << "Error: flag value is too large: " << max << std::endl;
        return std::vector<std::vector<int>>();
    }

    std::vector<std::vector<int>> color;
    color.resize(max + 1);
    for (int i = 0; i < flag.size(); i++) {
        int r = 0;
        int g = 0;
        int b = 0;
        if (max - min == 0) {
            r = 255;
            g = 255;
            b = 255;
        }
        else {
            r = (flag[i] - min) * 255 / (max - min);
            g = 255 - r;
            b = r / 2;
        }
        color[flag[i]] = { r,g,b };
    }
    return color;
}



std::vector<std::vector<int>> lzd_tools::flag2color(const std::vector<int>& flag){
    std::vector<std::vector<int>> color;
    if(flag.size() == 0){
        return color;
    }

    std::vector<int> flag_unique;
    for(int i=0;i<flag.size();i++){
        if(std::find(flag_unique.begin(), flag_unique.end(), flag[i]) == flag_unique.end()){
            flag_unique.push_back(flag[i]);
        }
    }

    srand((unsigned)time(NULL));
    std::vector<std::vector<int>> color_list;
    for(int i=0;i<flag_unique.size();i++){
        color_list.push_back({rand()%255,rand()%255,rand()%255});
    }
    for(int i=0;i<flag.size();i++){
        for(int j=0;j<flag_unique.size();j++){
            if(flag[i] == flag_unique[j]){
                color.push_back(color_list[j]);
                break;
            }
        }
    }

    return color;
}

int rmkdir(const std::string path){
    if (path.empty()) {
        return 0;
    }
    std::error_code ec;
    std::filesystem::create_directories(path, ec);
    return ec ? -1 : 0;
}

int rmkdir(const std::string base_path,const std::vector<std::string> subpaths){
    std::string cur_path = base_path;
    cur_path += "/";
    for(int i=0;i<subpaths.size();i++){
        cur_path += subpaths[i];
        cur_path += "/";
    }
    return rmkdir(cur_path);
}

lzd_tools::FiniteMap lzd_tools::get_regular_colormap(int min, int max) {
    auto colormap = lzd_tools::FiniteMap();
    int stride = 255 / (max - min);
    assert(stride > 0);
    colormap.insert({ min, {0, 0, 0} });
    for (int i = min+1; i <= max; i++) {
        colormap.insert({ i, {i * stride,255 - i * stride, 0} });
    }
    return colormap;
}

std::vector<int> lzd_tools::squeeze_array(const std::vector<int>& _array, int start) {
    std::vector<int> res(_array.size());
    std::map<int, int> map;
    int count = start;
    for (int i = 0; i < _array.size(); i++) {
        if (map.find(_array[i]) == map.end()) {
            map[_array[i]] = count++;
        }
        res[i] = map[_array[i]];
    }
    return res;
}

std::vector<int> lzd_tools::read_array(const std::string filename) {
    std::vector<int> _label;
    std::ifstream data(filename, std::ios::in);
    if (!data.is_open()) {
        std::cout << "Error: open " << filename << " failed" << std::endl;
        return _label;
    }
    std::string line;
    while (getline(data, line)) {
        std::stringstream ss(line);
        int label;
        ss >> label;
        _label.push_back(label);
    }
    return _label;
}

// enum linked_graph_filp_alg_opt{
//         GREEDY,
//         BRUTE_FORCE,
//         SINGLE_TREE,
//         MUTI_TREE,
//         OPTIMAL_FILP
// };
//    {global_var::linked_graph_filp_alg_opt::GREEDY,"GREEDY"},
//    {global_var::linked_graph_filp_alg_opt::BRUTE_FORCE,"BRUTE_FORCE"},
//    {global_var::linked_graph_filp_alg_opt::SINGLE_TREE,"SINGLE_TREE"},
//    {global_var::linked_graph_filp_alg_opt::MUTI_TREE,"MUTI_TREE"},
//    {global_var::linked_graph_filp_alg_opt::OPTIMAL_FILP,"OPTIMAL_FILP"},
//    {global_var::linked_graph_filp_alg_opt::BEST_TREE,"BEST_TREE"}
//};

std::string global_var::data_output_base = "";


nlohmann::json global_var::config_j;
nlohmann::json global_var::metric_j;
nlohmann::json global_var::res_log_j;

std::vector<std::vector<unsigned int>> global_var::color_map_100 = {
    {217,151,75},
    {149,168,184},
    {80,144,109},
    {250,236,229},
    {116,80,133},
    {248,187,208},
    {100,80,90},
    {245,245,220},
    {0,191,255},
    {75,83,32},
    {64,64,64},
    {224,255,255},
    {255,228,181},
    {255,127,80},
    {154,205,50},
    {218,112,214},
    {205,92,92},
    {255,69,0},
    {240,128,128},
    {152,251,152},
    {135,206,235},
    {255,105,180},
    {102,205,170},
    {72,61,139},
    {25,25,112},
    {245,245,245},
    {128,128,0},
    {127,255,212},
    {47,79,79},
    {139,0,139},
    {173,255,47},
    {0,255,255},
    {70,130,180},
    {255,250,205},
    {85,107,47},
    {147,112,219},
    {244,164,96},
    {50,205,50},
    {255,239,213},
    {255,218,185},
    {105,105,105},
    {176,224,230},
    {255,222,173},
    {189,183,107},
    {60,179,113},
    {143,188,143},
    {255,215,0},
    {199,21,133},
    {139,69,19},
    {0,0,205},
    {255,140,0},
    {34,139,34},
    {255,192,203},
    {173,216,230},
    {210,105,30},
    {139,0,0},
    {220,20,60},
    {255,250,240},
    {46,139,87},
    {255,99,71},
    {255,248,220},
    {32,178,170},
    {255,165,0},
    {0,250,154},
    {255,228,225},
    {144,238,144},
    {0,206,209},
    {255,255,224},
    {160,82,45},
    {135,206,250},
    {255,228,196},
    {138,43,226},
    {165,42,42},
    {255,228,181},
    {127,255,0},
    {72,209,204},
    {255,235,205},
    {255,240,245},
    {244,164,96},
    {222,184,135},
    {245,245,220},
    {75,0,130},
    {255,20,147},
    {147,112,219},
    {128,128,128},
    {124,252,0},
    {176,196,222},
    {106,90,205},
    {0,128,0},
    {100,149,237},
    {255,0,255},
    {0,255,127},
    {233,150,122},
    {240,255,255},
    {192,192,192},
    {255,182,193},
    {0,0,128},
    {255,250,205},
    {0,191,255},
    {169,169,169},
    {255,255,240},
};

//nlohmann::json global_var::Logger::config_j;
//nlohmann::json global_var::Logger::metric_j;
