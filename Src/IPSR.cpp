#include <IPSR.h>
//#include <PoissonRecon.h>

std::map<int (*)(Period& p),std::string> update_plane_namelist{
	{interact_update,"interact_update"},
	{diri_ind_ipsr,"diri_ind_ipsr"},
	{diri_dir_ipsr,"diri_dir_ipsr"},
	{neumman_ind_ipsr,"neumman_ind_ipsr"},
	{neumman_dir_ipsr,"neumman_dir_ipsr"},
	{comb_update7,"comb_update7"},
	{comb_update8,"comb_update8"},
    {specify_update,"specify_update"}
};

// std::map<int,std::string> init_type_namelist={
// 	{INIT_RANDOM,"INIT_RANDOM"},
// 	{INIT_PCA,"INIT_PCA"},
// 	{INIT_PCL_EST_NORMAL,"INIT_PCL_EST_NORMAL"},
// 	{INIT_NOTHING,"INIT_NOTHING"},
// };
// 交互式更新
int interact_update(Period& p) {
    if(p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr){
        return -1;
    }
    p.print_mesh = 1;
    static int count = 0;
    if(count  > 0){
        //printf("count = %d\n", count);
        count--;
        return 0;
    }
    p.print_op_freq = 1;
    int type, strategy;
    printf("input Ttype(2/3) strategy(0/1) times\n");
    scanf("%d %d %d", &type, &strategy, &count);
    change_args(*(p.p_cmd),"--bType",std::to_string(type));
    *p.update_strategy = strategy;
    count--;
    return 0;
}

int diri_ind_ipsr(Period& p)
{
    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        return -1;
    }
    *p.update_strategy = INDIRECT_MAPPING;
    change_args(*p.p_cmd, "--bType", "2");
    printf("using Dirichlet and INDIRECT_MAPPING\n");
    return 0;
}

int diri_dir_ipsr(Period& p)
{
    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        assert(false);
        return -1;
    }
    *p.update_strategy = DIRECT_MAPPING;
    change_args(*p.p_cmd, "--bType", "2");
    printf("using Dirichlet and DIRECT_MAPPING\n");
    return 0;
}

int neumman_ind_ipsr(Period& p)
{
    static lzd_tools::thread_safe_bool has_print(false);
    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        return -1;
    }
    *p.update_strategy = INDIRECT_MAPPING;
    if((*p.p_cmd)[find_arg(*p.p_cmd,"--bType") + 1] == "3"){
        return 0;
    }
    change_args(*p.p_cmd, "--bType", "3");
    if (has_print.get()) {
        printf("using Neumman and INDIRECT_MAPPING\n");
        has_print.set(true);
    }
    return 0;
}

int neumman_dir_ipsr(Period& p)
{
    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        return -1;
    }
    *p.update_strategy = DIRECT_MAPPING;
    change_args(*p.p_cmd, "--bType", "3");
    printf("using Neumman and DIRECT_MAPPING\n");
    return 0;
}

int comb_plan5(Period& p) {
    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        return -1;
    }

    if (*(p.p_ep) < 30) {
        change_args(*(p.p_cmd), "--bType", "2");
        *p.update_strategy = INDIRECT_MAPPING;
        printf("using Dirichlet ");
        printf("using INDIRECT_MAPPING\n");
    }
    else {
        change_args(*(p.p_cmd), "--bType", "3");
        printf("using Neumann ");
        int current = *p.update_strategy;
        bool if_change = *p.p_new_var > *p.p_old_var && *p.p_new_var < 1.0;
        if (if_change) {
            if (current == DIRECT_MAPPING) {
                *p.update_strategy = INDIRECT_MAPPING;
                printf("using INDIRECT_MAPPING\n");
            }
            else {
                *p.update_strategy = DIRECT_MAPPING;
                printf("using DIRECT_MAPPING\n");
            }
        }
    }
    return 0;
}

int comb_update6(Period& p) {
    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        return -1;
    }

    int tType_idx = find_arg(*(p.p_cmd), "--bType");

    if (*(p.p_ep) < 20) {
        change_args(*(p.p_cmd), "--bType", "2");
        *p.update_strategy = INDIRECT_MAPPING;
        printf("using Dirichlet ");
        printf("using INDIRECT_MAPPING\n");
        return 0;
    }
    else if (*(p.p_ep) > (p.max_iters - 10)) {
        change_args(*(p.p_cmd), "--bType", "3");
        *p.update_strategy = DIRECT_MAPPING;
        printf("end witg Neumann and DIRECT_MAPPING\n");
        return 0;
    }

    static int count = 0;

    int current = *p.update_strategy;
    bool if_change = (*p.p_new_var > *p.p_old_var && *p.p_new_var < 1.0) || count > 10;
    if (!if_change) {
        count++;
    }
    else {
        count = 0;
        if (current == DIRECT_MAPPING && (*p.p_cmd)[tType_idx + 1] == "3") {
            *p.update_strategy = INDIRECT_MAPPING;
            change_args(*(p.p_cmd), "--bType", "2");
            printf("using Dirichlet ");
            printf("using INDIRECT_MAPPING\n");
        }
        else if (current == INDIRECT_MAPPING && (*p.p_cmd)[tType_idx + 1] == "2") {
            *p.update_strategy = DIRECT_MAPPING;
            change_args(*(p.p_cmd), "--bType", "3");
            printf("using Neumann ");
            printf("using DIRECT_MAPPING\n");
        }
        else {
            printf("\n\nERROR boundary type and strategy not mactch\n\n");
        }

    }
    return 0;
}

/**
 * @brief 
 * 前20次使用Neumann+INDIRECT_MAPPING,最后10次固定使用Dirichlet+DIRECT_MAPPING
 * 中间的情况根据变量的变化来判断
 * @param p 
 * @return int 
 */
int comb_update7(Period& p) {
    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        return -1;
    }

    int tType_idx = find_arg(*(p.p_cmd), "--bType");

    if (*(p.p_ep) < 20) {
        change_args(*(p.p_cmd), "--bType", "3");
        *p.update_strategy = INDIRECT_MAPPING;
        printf("using Neumann ");
        printf("using INDIRECT_MAPPING\n");
        return 0;
    }
    else if (*(p.p_ep) > (p.max_iters - 10)) {
        change_args(*(p.p_cmd), "--bType", "2");
        *p.update_strategy = DIRECT_MAPPING;
        printf("end with Dirichlet and DIRECT_MAPPING\n");
        return 0;
    }

    static int count = 0;
    int current = *p.update_strategy;
    bool if_change = (*p.p_new_var > *p.p_old_var && *p.p_new_var < 1.0) || count > 10;
    if (!if_change) {
        count++;
    }
    else {
        count = 0;
        if (current == DIRECT_MAPPING && (*p.p_cmd)[tType_idx + 1] == "2") {
            *p.update_strategy = INDIRECT_MAPPING;
            change_args(*(p.p_cmd), "--bType", "3");
            printf("using Neumann ");
            printf("using INDIRECT_MAPPING\n");
        }
        else if (current == INDIRECT_MAPPING && (*p.p_cmd)[tType_idx + 1] == "3") {
            *p.update_strategy = DIRECT_MAPPING;
            change_args(*(p.p_cmd), "--bType", "2");
            printf("using Dirichlet ");
            printf("using DIRECT_MAPPING\n");
        }
        else {
            printf("\n\nERROR boundary type and strategy not mactch\n\n");
        }

    }
    return 0;
}

/**
 * @brief 
 * 前10次使用Dirichlet+INDIRECT_MAPPING,以得到一个比较稳定的初值
 * 此后与comb_update7一致
 * @return int 
 */
int comb_update8(Period& p) {
    if(p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr){
        return -1;
    }
    // 注意 这个阶段的ep不能大于comb8第一阶段的ep
    if(*(p.p_ep) < 10){
        change_args(*(p.p_cmd), "--bType", "2");
        *p.update_strategy = INDIRECT_MAPPING;
        printf("using Dirichlet\n");
        printf("using INDIRECT_MAPPING\n");
    }
    else{
        comb_update7(p);
    }
    return 0;
}

int specify_update(Period& p)
{

    if (p.p_ep == nullptr || p.p_old_var == nullptr || p.p_new_var == nullptr || p.p_cmd == nullptr) {
        return -1;
    }
    static nlohmann::json config = ConfigManager::get_config_in_config_floder("specify_update.json",false);


    auto get_iter = [](const int& iter) {
        int lb = 0;
        for (auto& item : config.items()) {
            if (std::stoi(item.key()) <= iter) {
                lb = std::stoi(item.key());
            }
            else {
                return lb;
            }
        }
        return lb;
     };
    nlohmann::json iter_config = config[std::to_string(get_iter(*(p.p_ep)))];
    for (const auto &item : iter_config.items()) {
        assert(!item.key().empty() && !iter_config[item.key()].is_object());
        change_args(*(p.p_cmd), item.key(), iter_config[item.key()].dump());
    }
    if (iter_config.find("update_strategy")!=iter_config.end()) {
        *p.update_strategy = std::stoi(iter_config["update_strategy"].dump());
    }
    else {
        *p.update_strategy = INDIRECT_MAPPING;
    }

    return 0;
}

UpdatePlan get_update_plan_by_name(std::string name)
{
    std::map<std::string, UpdatePlan> update_plan_map;
    for (auto& it : update_plane_namelist) {
        update_plan_map[it.second] = it.first;
    }
    if (update_plan_map.find(name) != update_plan_map.end()) {
        return update_plan_map[name];
    }
    else {
        printf("no such plan %s", name);
        assert(false);
        return nullptr;
    }
}

Period::Period()
{    
}

//lzd_tools::FiniteMap get_mincut_colormap() {
//    auto colormap = lzd_tools::FiniteMap({ {SOURCE_PART_LABEL, {255,0,0}},{SINK_PART_LABEL, {0,255,0}},{SOURCE_LABEL, {0,0,0}},{SINK_LABEL, {0,0,255}} });
//    // 考虑其他情况
//    colormap.insert({ -1, {255,255,255} });
//    return colormap;
//}