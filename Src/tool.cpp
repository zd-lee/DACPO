#include <tools.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>


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
    // 遍历args,将值为key的value替换为value
    for(int i=0;i<args.size()-1;i++){
        if(args[i] == key){
            args[i+1] = value;
            return i;
        }
    }
    printf("\n\n!!chane_args error : key %s not found!!\n\n",key);
    return -1;
}

// 返回key在args中的索引。
// 例如:args[find_arg(args,"--bType")+1] = "2",则将bType的值改为2
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
    // 获得flag中的最大值和最小值
    int max = *std::max_element(flag.begin(), flag.end());
    int min = *std::min_element(flag.begin(), flag.end());
    // 判断flag中是否有-1
    if (min < 0) {
        std::cout << "Error: flag中有-1" << std::endl;
        return std::vector<std::vector<int>>();
    }
    if (max > 1e5) {
        std::cout << "Error: flag中有过大的值:" << max << std::endl;
        return std::vector<std::vector<int>>();
    }

    // 生成颜色
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



std::vector<std::vector<int>> flag2color(const std::vector<int>& flag){
    // 将flag转换为颜色
    std::vector<std::vector<int>> color;
    if(flag.size() == 0){
        return color;
    }

    // 方案一:渐变色
    // int _min_flag = *std::min_element(flag.begin(), flag.end());
    // int _max_flag = *std::max_element(flag.begin(), flag.end());
    // for(int i=0;i<flag.size();i++){
    //     int r = 0;
    //     int g = 0;
    //     int b = 0;
    //     if(_max_flag - _min_flag == 0){
    //         r = 255;
    //         g = 255;
    //         b = 255;
    //     }else{
    //         r = (flag[i] - _min_flag) * 255 / (_max_flag - _min_flag);
    //         g = 255 - r;
    //         b = r/2;
    //     }
    //     color.push_back({r,g,b});
    // }

    // 方案二:固定颜色
    // for(int i=0;i<flag.size();i++){
    //     if(flag[i] == 0){
    //         color.push_back({0,0,255});
    //     }else if(flag[i] == 1){
    //         color.push_back({255,0,0});
    //     }else if(flag[i] == 2){
    //         color.push_back({0,255,0});
    //     }else{
    //         color.push_back({255,255,255});
    //     }
    // }

    // 方案三:随机颜色(开销较大)
    // 获得flag中的不同的值的个数
    std::vector<int> flag_unique;
    for(int i=0;i<flag.size();i++){
        if(std::find(flag_unique.begin(), flag_unique.end(), flag[i]) == flag_unique.end()){
            flag_unique.push_back(flag[i]);
        }
    }

    // 每个值对应一个随机颜色
    srand((unsigned)time(NULL));
    // 颜色列表
    std::vector<std::vector<int>> color_list;
    for(int i=0;i<flag_unique.size();i++){
        color_list.push_back({rand()%255,rand()%255,rand()%255});
    }
    // 将flag中的值转换为颜色
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
    // 创建文件夹
    // 返回值:0成功,-1失败
    std::vector<std::string> sub_path;
    std::string tmp = "";
    for(int i=0;i<path.size();i++){
        if(path[i] == '/'){
            sub_path.push_back(tmp);
            tmp = "";
        }else{
            tmp += path[i];
        }
    }
    sub_path.push_back(tmp);
    
    std::string cur_path = "";
    for(int i=0;i<sub_path.size();i++){
        cur_path += sub_path[i];
        if(GetFileAttributesA(cur_path.c_str()) == INVALID_FILE_ATTRIBUTES){
            // 文件夹不存在
            if(CreateDirectoryA(cur_path.c_str(), NULL) == 0){
                // 创建失败
                return -1;
            }else{
                printf("create dir %s\n",cur_path.c_str());
            }
        }
        cur_path += "/";
    }
    return 0;
}

int rmkdir(const std::string base_path,const std::vector<std::string> subpaths){
    // 创建文件夹
    // 返回值:0成功,-1失败
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
    // 得到min到max的渐进色
    int stride = 255 / (max - min);
    assert(stride > 0);
    colormap.insert({ min, {0, 0, 0} });
    for (int i = min+1; i <= max; i++) {
        colormap.insert({ i, {i * stride,255 - i * stride, 0} });
    }
    //// 考虑其他情况
    //colormap.insert({ -1, {255,255,255} });
    return colormap;
}

// 将数组中的值压缩到从start开始的连续整数
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
//std::map<int,std::string> global_var::linked_graph_filp_alg_name = {
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
    {217,151,75},    // 沙漠黄
    {149,168,184},   // 灰蓝
    {80,144,109},    // 海藻绿
    {250,236,229},   // 亮粉
    {116,80,133},    // 紫灰
    {248,187,208},   // 淡粉
    {100,80,90},     // 暗紫
    {245,245,220},   // 米色
    {0,191,255},     // 深天蓝
    {75,83,32},      // 橄榄绿
    {64,64,64},      // 暗灰
    {224,255,255},   // 淡青
    {255,228,181},   // 小麦色
    {255,127,80},    // 珊瑚红
    {154,205,50},    // 黄绿色
    {218,112,214},   // 兰花紫
    {205,92,92},     // 印度红
    {255,69,0},      // 橙红
    {240,128,128},   // 浅珊瑚红
    {152,251,152},   // 苍绿色
    {135,206,235},   // 天蓝色
    {255,105,180},   // 热粉色
    {102,205,170},   // 中绿宝石
    {72,61,139},     // 暗蓝紫
    {25,25,112},     // 午夜蓝
    {245,245,245},   // 亮灰
    {128,128,0},     // 橄榄色
    {127,255,212},   // 绿松石
    {47,79,79},      // 暗石板灰
    {139,0,139},     // 深紫
    {173,255,47},    // 绿黄
    {0,255,255},     // 青色
    {70,130,180},    // 钢蓝
    {255,250,205},   // 柠檬黄
    {85,107,47},     // 橄榄褐
    {147,112,219},   // 中紫
    {244,164,96},    // 沙棕色
    {50,205,50},     // 酸橙绿
    {255,239,213},   // 白杏色
    {255,218,185},   // 桃色
    {105,105,105},   // 暗灰
    {176,224,230},   // 淡钢蓝
    {255,222,173},   // 浅肉色
    {189,183,107},   // 暗卡其色
    {60,179,113},    // 中海绿
    {143,188,143},   // 暗海绿
    {255,215,0},     // 金色
    {199,21,133},    // 红紫色
    {139,69,19},     // 巧克力
    {0,0,205},       // 中蓝色
    {255,140,0},     // 深橙
    {34,139,34},     // 森林绿
    {255,192,203},   // 粉红色
    {173,216,230},   // 淡蓝色
    {210,105,30},    // 朱红
    {139,0,0},       // 暗红
    {220,20,60},     // 猩红
    {255,250,240},   // 象牙白
    {46,139,87},     // 海绿
    {255,99,71},     // 番茄红
    {255,248,220},   // 玉米色
    {32,178,170},    // 轻海绿
    {255,165,0},     // 橙色
    {0,250,154},     // 中春绿
    {255,228,225},   // 薄红
    {144,238,144},   // 淡绿
    {0,206,209},     // 暗青
    {255,255,224},   // 浅黄
    {160,82,45},     // 棕褐色
    {135,206,250},   // 淡天蓝
    {255,228,196},   // 雪色
    {138,43,226},    // 蓝紫色
    {165,42,42},     // 棕色
    {255,228,181},   // 莲雾黄
    {127,255,0},     // 荧光绿
    {72,209,204},    // 深绿松石
    {255,235,205},   // 米黄色
    {255,240,245},   // 薄荷白
    {244,164,96},    // 沙棕
    {222,184,135},   // 黄褐色
    {245,245,220},   // 米色
    {75,0,130},      // 靛色
    {255,20,147},    // 深粉红
    {147,112,219},   // 中蓝紫
    {128,128,128},   // 灰色
    {124,252,0},     // 草绿色
    {176,196,222},   // 浅石板灰
    {106,90,205},    // 板岩蓝
    {0,128,0},       // 绿色
    {100,149,237},   // 矢车菊蓝
    {255,0,255},     // 紫红色
    {0,255,127},     // 春绿
    {233,150,122},   // 暗鲑红
    {240,255,255},   // 苍白蓝
    {192,192,192},   // 银色
    {255,182,193},   // 浅粉色
    {0,0,128},       // 海军蓝
    {255,250,205},   // 柠檬绿
    {0,191,255},     // 深天蓝
    {169,169,169},   // 暗灰色
    {255,255,240},   // 象牙
};

//indicators::ProgressBar bar{
//indicators::option::BarWidth{50},
//indicators::option::Start{"["},
//indicators::option::Fill{"="},
//indicators::option::Lead{">"},
//indicators::option::Remainder{" "},
//indicators::option::End{"]"},
//indicators::option::PostfixText{"Extracting Archive"},
//indicators::option::ForegroundColor{indicators::Color::green},
//indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
//};

//nlohmann::json global_var::Logger::config_j;
//nlohmann::json global_var::Logger::metric_j;