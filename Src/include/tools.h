#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <PoissionEntrance.h>
#include "dbscan.h"
#include "nlohmann/json.hpp"
// #include "indicators/indicators.hpp"
// #include <thread>
// #include <chrono>

// 定义一对多的映射，表示两个集合之间的一个一对多的映射，要求两个集合的元素都是自然数。
#define MULTI_MAP std::vector<std::vector<int>>

namespace lzd_tools{
    
    struct thread_safe_bool {
        bool val;
        std::mutex mtx;
        thread_safe_bool(bool val) :val(val) {}
        thread_safe_bool() :val(false) {}
        void set(bool _val) {
            mtx.lock();
            val = _val;
            mtx.unlock();
        }
        bool get() {
            mtx.lock();
            bool res = val;
            mtx.unlock();
            return res;
        }
    };

    struct thread_safe_int {
        int val;
        std::mutex mtx;
        thread_safe_int(int val) :val(val) {}
        thread_safe_int() :val(0) {}
        void operator++() {
            mtx.lock();
            val++;
            mtx.unlock();
        }
        void operator--() {
            mtx.lock();
            val--;
            mtx.unlock();
        }


        int get() {
            mtx.lock();
            int res = val;
            mtx.unlock();
            return res;
        }
        void set(int _val) {
            mtx.lock();
            val = _val;
            mtx.unlock();
        }
    };

    /**
     * @brief 
     * 计时器类,用于计算函数运行时间.
     * 使用方式:
     *     nlohmann::json j;
     *     Timer timer(j,"func_name");
     * 
     *     // timer.Start(); // 可以手动开始计时,也可以不开始,构造函数会自动开始计时
     *     // do something
     *     // timer.Done(); // 可以手动结束计时,也可以不结束,析构函数会自动结束计时  
     */
    class Timer {
    private:
        std::chrono::time_point<std::chrono::steady_clock> start;
        nlohmann::json& _j;
        bool _done = false;
        std::string _func_name;

    public:
        Timer(nlohmann::json& j, std::string func);

        void Start();

        void Done();

        ~Timer();
    };

    class AcumutableTimer {
    private:
        static std::map<std::string, double> _time_map;
        static std::map<std::string, int> _count_map;
        static std::mutex _mtx;
        std::chrono::time_point<std::chrono::steady_clock> start;
        std::string _func_name;
        thread_safe_bool _done;

    public:
        static double get_time(std::string func_name);
        AcumutableTimer(std::string func_name);
        void Done();
        ~AcumutableTimer();
        static nlohmann::json get_time_map();
    };




    
    struct resource_controller {
        thread_safe_int _resource;
        resource_controller(int resource) :_resource(resource) {}

        void P() {
            while (_resource.get() <= 0);
            --_resource;
        }
        void V() {
            ++_resource;
        }


    };
   
    /**************************************颜色映射******************************************/

    // 有限映射 将有限值域映射到RGB空间上
    typedef std::map<int, std::vector<int>> FiniteMap;
  
    /**
     * @brief 
     * 无限映射的虚基类 将一个flag数组中的元素的定义域中的所有值映射到RGB空间上
     * 所有子类需要实现toRGB函数,对于任意的_xval,返回一个RGB颜色
     * 如果flag数组为空,则不标记颜色
     * 使用方式:
     * lzd_tools::InFiniteMap<int>* color_map =new lzd_tools::RandomColorMapper<int>(seg_id);
     * lzd_tools::op2ply(op, path, XForm<REAL, DIM + 1>().Identity(), color_map);
     * ...
     * if (color_map->has_color) {
     *      std::vector<int> colors = (*color_map)[i];
     * }
     * 
     * @tparam REAL 
     */
    template<typename REAL>
    class InFiniteMap {
    public:
        std::vector<REAL> _flag;
        bool has_color;
        
        InFiniteMap(std::vector<REAL> flag) :_flag(flag){
            has_color = flag.size() > 0;
        }

        std::vector<int> operator[](int i) {
            assert(has_color);
            return toRGB(_flag[i]);
        }
        virtual std::vector<int> toRGB(REAL _xval)=0;
    };

    template<typename REAL>
    class RandomColorMapper:public InFiniteMap<REAL> {
    public:
        RandomColorMapper(std::vector<REAL> flag):InFiniteMap<REAL>(flag) {}
        std::vector<int> toRGB(REAL _xval){
            std::vector<int> res;
            static long long prim[3] = { 1980921,34523,9329121 };
            for (int i = 0; i < 3; i++) {
                long long t = _xval * prim[i];
                t %= 255;
                res.push_back(t);
            }
            return res;
        }
    };

    /**
     * @brief 
     * 得到渐变色 将从大到小的值映射为红色-绿色 
     */
    template<typename REAL>
    class GradientMapper :public InFiniteMap<REAL> {
        REAL gap;
    public:
        std::vector<REAL> _flag;
        REAL _max_val;
        REAL _min_val;

        GradientMapper(std::vector<REAL> flag):InFiniteMap<REAL>(flag) {
            if (flag.size() == 0) {
                return;
            }
            _max_val = *std::max_element(flag.begin(), flag.end());
            _min_val = *std::min_element(flag.begin(), flag.end());            
            gap = _max_val - _min_val; // 
        }
        std::vector<int> toRGB(REAL _xval) {
            int t = 254 * (_xval - _min_val) / gap;//防止max变成min..
            t %= 255;
            return { t,255 - t,t / 2 };   
        }
    };

    /**
     * @brief 
     * 提供四种好看的颜色;如果flag中不同的值超过四个,之后出现的值将会被映射到第五种颜色(黑色)上
     * @tparam REAL 
     */
    template<typename REAL>
    class beautifulColorMapper:public InFiniteMap<REAL> {
    private:
        std::vector<std::vector<int>> _color;
        std::map<REAL, int> _map;//_xval -> color_idx

    public:
        std::vector<REAL> _flag;

        enum COLOR_NAME{
            DESERT_YELLOW,
            GRAY_BLUE,
            SEAWEED_GREEN,
            BRIGHT_PINK,
            BLACK
        };

        
        beautifulColorMapper(std::vector<REAL> flag):InFiniteMap<REAL>(flag) {
            _color = {
                {217,151,75},//沙漠黄
                {149,168,184},//灰蓝
                {80,144,109},//海藻绿
                {250,236,229},//亮粉
                {40,30,16}//黑
            };

            //初始化map
            for (int i = 0; i < flag.size(); i++) {
                if (_map.size() >= 4) {
                    break;
                }
                
                if (_map.find(flag[i]) == _map.end()) {
                    _map[flag[i]] = _map.size();
                }

            }
        }

        std::vector<int> toRGB(REAL _xval) {
            if (_map.find(_xval) == _map.end()) {
                return _color[4];
            }
            if (_map[_xval] < 4) {
                return _color[_map[_xval]];
            }
            else {
                return _color[4];
            }
        }
    };
    
    template <typename REAL, int DIM, typename ColorType>
    int op2ply(const ORIENTED_POINTS& op, const std::string path, const XForm<REAL, DIM + 1> xform, lzd_tools::InFiniteMap<ColorType>* color_map)
    {
        std::ofstream out(path);
        out << "ply" << std::endl;
        out << "format ascii 1.0" << std::endl;
        out << "element vertex " << op.size() << std::endl;
        // 如果有flag字段,则将不同的点标记为不同的颜色
        out << "property float x" << std::endl;
        out << "property float y" << std::endl;
        out << "property float z" << std::endl;
        out << "property float nx" << std::endl;
        out << "property float ny" << std::endl;
        out << "property float nz" << std::endl;
        if (color_map->has_color) {
            out << "property uchar red" << std::endl;
            out << "property uchar green" << std::endl;
            out << "property uchar blue" << std::endl;
        }
        out << "end_header" << std::endl;

        for (int i = 0; i < op.size(); i++) {
            Point<REAL, DIM> p = xform * op[i].first;
            auto n = op[i].second.normal;
            out << p[0] << " " << p[1] << " " << p[2];
            out << " " << n[0] << " " << n[1] << " " << n[2];

            if (color_map->has_color) {
                out << " ";
                std::vector<int> colors = (*color_map)[i];
                int r = colors[0], g = colors[1], b = colors[2];
                out << r << " " << g << " " << b << " ";

            }
            out << std::endl;
        }
        return 0;
    }

    /**
     * @brief
     * 绘制有向点云(InFiniteMap版本)
     * @param flagandColorMap 标记,以及标记的颜色列表.
     */
    template <typename REAL, int DIM>
    int op2ply(const ORIENTED_POINTS& op, const std::string path,const XForm<REAL, DIM + 1> xform,lzd_tools::InFiniteMap<int>* color_map)
    {
        std::ofstream out(path);
        out << "ply" << std::endl;
        out << "format ascii 1.0" << std::endl;
        out << "element vertex " << op.size() << std::endl;
        // 如果有flag字段,则将不同的点标记为不同的颜色
        out << "property float x" << std::endl;
        out << "property float y" << std::endl;
        out << "property float z" << std::endl;
        out << "property float nx" << std::endl;
        out << "property float ny" << std::endl;
        out << "property float nz" << std::endl;
        if (color_map->has_color) {
            out << "property uchar red" << std::endl;
            out << "property uchar green" << std::endl;
            out << "property uchar blue" << std::endl;
        }
        out << "end_header" << std::endl;

        for (int i = 0; i < op.size(); i++) {
            Point<REAL, DIM> p = xform * op[i].first;
            auto n = op[i].second.normal;
            out << p[0] << " " << p[1] << " " << p[2];
            out << " " << n[0] << " " << n[1] << " " << n[2];

            if (color_map->has_color) {
                out << " ";
                std::vector<int> colors = (*color_map)[i];
                int r = colors[0],g = colors[1], b = colors[2];
                out << r << " " << g << " " << b << " ";
                
            }
            out << std::endl;
        }
        return 0;
    }

    /**
     * @brief 
     * 绘制有向点云(FiniteMap版本)
     * @param flagandColorMap 标记,以及标记的颜色列表.
     */
    template <typename REAL, int DIM>
    int op2ply(const ORIENTED_POINTS& op, const std::string path, 
        const XForm<REAL, DIM+1> xform = XForm<REAL, DIM+1>().Identity(), 
        std::pair<std::vector<int>, FiniteMap> flagandColorMap= std::make_pair(std::vector<int>(), lzd_tools::FiniteMap()))
    {
        auto flag = flagandColorMap.first;
        auto color_map = flagandColorMap.second;       
        std::ofstream out(path);
        out << "ply" << std::endl;
        out << "format ascii 1.0" << std::endl;
        out << "element vertex " << op.size() << std::endl;
        // 如果有flag字段,则将不同的点标记为不同的颜色
        out << "property float x" << std::endl;
        out << "property float y" << std::endl;
        out << "property float z" << std::endl;
        out << "property float nx" << std::endl;
        out << "property float ny" << std::endl;
        out << "property float nz" << std::endl;
        if(flag.size() == op.size()){
            out << "property uchar red" << std::endl;
            out << "property uchar green" << std::endl;
            out << "property uchar blue" << std::endl;
        }
        out << "end_header" << std::endl;
        auto rand_color = flag2color(flag);

        for(int i=0;i<op.size();i++){
            Point<REAL, DIM> p = xform * op[i].first;
            auto n = op[i].second.normal;
            out << p[0] << " " << p[1] << " " << p[2];
            out << " " << n[0] << " " << n[1] << " " << n[2];

            if(flag.size() == op.size()){
                out<< " ";
                if(color_map.size() == 0){
                    int r = rand_color[flag[i]][0];
                    int g = rand_color[flag[i]][1];
                    int b = rand_color[flag[i]][2];
                    out << r << " " << g << " " << b << " ";
                }else{
                    assert(color_map.find(flag[i])!=color_map.end());
                    int r = color_map[flag[i]][0];
                    int g = color_map[flag[i]][1];
                    int b = color_map[flag[i]][2];
                    out << r << " " << g << " " << b << " ";
                }
            }
            out << std::endl;
        }
        return 0;
    }

    /**
     * @brief 
     * 绘制有色surface(可以指定每个面的颜色)
     */
     template <typename REAL, int DIM>
     int mesh2ply(const MESH& mesh, const std::string path, const XForm<REAL, DIM+1> xform = XForm<REAL, DIM+1>().Identity(), 
        std::pair<std::vector<int>, FiniteMap> flagandColorMap= std::make_pair(std::vector<int>(), lzd_tools::FiniteMap())){
        const int n = mesh.first.size();
        const int m = mesh.second.size();

        bool has_color = flagandColorMap.first.size() == m;

        std::ofstream out(path);
        out << "ply" << std::endl;
        out << "format ascii 1.0" << std::endl;
        out << "element vertex " << n << std::endl;
        out << "property float x" << std::endl;
        out << "property float y" << std::endl;
        out << "property float z" << std::endl;

        out << "element face " << m << std::endl;
        out << "property list uchar int vertex_index" << std::endl;
        if(has_color){
            out << "property uchar red" << std::endl;
            out << "property uchar green" << std::endl;
            out << "property uchar blue" << std::endl;
        }
        out << "end_header" << std::endl;

        for(int i=0;i<n;i++){
            Point<REAL, DIM> p = xform * mesh.first[i];
            out << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
        for(int i=0;i<m;i++){
            out << "3 " << mesh.second[i][0] << " " << mesh.second[i][1] << " " << mesh.second[i][2] << " ";
            if(has_color){
                int r = flagandColorMap.second[flagandColorMap.first[i]][0];
                int g = flagandColorMap.second[flagandColorMap.first[i]][1];
                int b = flagandColorMap.second[flagandColorMap.first[i]][2];
                out << r << " " << g << " " << b << std::endl;
            }
        }
        return 0;
    }
    
    template <typename REAL, int DIM>
    int mesh2ply(const MESH& mesh, const std::string path, const XForm<REAL, DIM+1> xform = XForm<REAL, DIM+1>().Identity(), lzd_tools::InFiniteMap<int>* color_map = new lzd_tools::RandomColorMapper<int>(std::vector<int>())){
        const int n = mesh.first.size();
        const int m = mesh.second.size();
        std::ofstream out(path);
        out << "ply" << std::endl;
        out << "format ascii 1.0" << std::endl;
        out << "element vertex " << n << std::endl;
        // 如果有flag字段,则将不同的点标记为不同的颜色
        out << "property float x" << std::endl;
        out << "property float y" << std::endl;
        out << "property float z" << std::endl;
        if(color_map->has_color){
            out << "property uchar red" << std::endl;
            out << "property uchar green" << std::endl;
            out << "property uchar blue" << std::endl;
        }
        out << "element face " << m << std::endl;
        out << "property list uchar int vertex_index" << std::endl;
        out << "end_header" << std::endl;
        for (int i = 0; i < n; i++) {
            Point<REAL, DIM> p = xform * mesh.first[i];
            out << p[0] << " " << p[1] << " " << p[2];
            if(color_map->has_color){
                out<< " ";
                std::vector<int> colors = (*color_map)[i];
                int r = colors[0],g = colors[1], b = colors[2];
                out << r << " " << g << " " << b << std::endl;
            }else{
                out << std::endl;
            }
        }
        for (int i = 0; i < m; i++) {
            out << "3 " << mesh.second[i][0] << " " << mesh.second[i][1] << " " << mesh.second[i][2] << std::endl;
        }
        out.close();
        return 0;
    }

    template <typename REAL, int DIM, typename ColorType>
    int mesh2ply(const MESH& mesh, const std::string path, const XForm<REAL, DIM+1> xform = XForm<REAL, DIM+1>().Identity(), lzd_tools::InFiniteMap<ColorType>* color_map = new lzd_tools::RandomColorMapper<int>(std::vector<ColorType>())){
        const int n = mesh.first.size();
        const int m = mesh.second.size();
        std::ofstream out(path);
        out << "ply" << std::endl;
        out << "format ascii 1.0" << std::endl;
        out << "element vertex " << n << std::endl;
        // 如果有flag字段,则将不同的点标记为不同的颜色
        out << "property float x" << std::endl;
        out << "property float y" << std::endl;
        out << "property float z" << std::endl;
        if(color_map->has_color){
            out << "property uchar red" << std::endl;
            out << "property uchar green" << std::endl;
            out << "property uchar blue" << std::endl;
        }
        out << "element face " << m << std::endl;
        out << "property list uchar int vertex_index" << std::endl;
        out << "end_header" << std::endl;
        for (int i = 0; i < n; i++) {
            Point<REAL, DIM> p = xform * mesh.first[i];
            out << p[0] << " " << p[1] << " " << p[2];
            if(color_map->has_color){
                out<< " ";
                std::vector<int> colors = (*color_map)[i];
                int r = colors[0],g = colors[1], b = colors[2];
                out << r << " " << g << " " << b << std::endl;
            }else{
                out << std::endl;
            }
        }
        for (int i = 0; i < m; i++) {
            out << "3 " << mesh.second[i][0] << " " << mesh.second[i][1] << " " << mesh.second[i][2] << std::endl;
        }
        out.close();
        return 0;
    } 
    
    ///**
    // * @brief 
    // * 将kidxs对应的点,根据其与P所在切平面的距离排序
    // * @param data 
    // * @param p 查询点
    // * @param kidxs 待排序的点在data中的索引 
    // */
    //template <typename REAL, int DIM>
    //std::vector<int> tangentDistance(const ORIENTED_POINTS& data, const std::pair< Point<REAL,DIM>,Normal<REAL,DIM> >op,const std::vector<int>& kidxs){
    //}

    /**
     * @brief 将min-max映射到0-255 要求min-max不能大于255
     * @param min 
     * @param max 
     * @return FiniteMap 
     */
    FiniteMap get_regular_colormap(int min = 0, int max = 255);

    //lzd_tools::FiniteMap lzd_tools::get_regular_colormap(std::vector<int> op_type);

    /**
     * @brief 
     * 向mesh指定位置添加一个球体
     * @param center 球体的中心
     * @param radius 球体的半径
     * @param n 球体的纬度精细度
     * @param m 球体的经度精细度
     * @return MESH 返回添加球体后的mesh 
     */
    template <typename REAL, int DIM>
    MESH add_sphere(const MESH& mesh, Point<REAL, DIM> center, REAL radius = 0.01, int n = 10, int m = 10){
        int start = mesh.first.size();
        MESH res = mesh;//看看是深拷贝还是浅拷贝
        // 添加顶点
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                REAL theta = 2 * M_PI * i / n;
                REAL phi = M_PI * j / m;
                Point<REAL, DIM> p = center + Point<REAL, DIM>(radius * sin(phi) * cos(theta), radius * sin(phi) * sin(theta), radius * cos(phi));
                res.first.push_back(p);
            }
        }
        // 添加面
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                int a = i * m + j + start;
                int b = i * m + (j + 1) % m + start;
                int c = ((i + 1) % n) * m + j + start;
                int d = ((i + 1) % n) * m + (j + 1) % m + start;
                res.second.push_back({a,b,c});
                res.second.push_back({b,d,c});
            }
        }
        return res;
    }


    template <typename REAL, int DIM>
    void normalize(Point<REAL, DIM>& p){
        REAL sum = 0;
        for(int i=0;i<DIM;i++){
            sum += p[i] * p[i];
        }
        sum = sqrt(sum);
        for(int i=0;i<DIM;i++){
            p[i] /= sum;
        }
    }

    template <typename REAL, int DIM>
    MESH get_sphere(Point<REAL, DIM> center, REAL radius = 0.01, int n = 10, int m = 10){
        MESH res;
        // 添加顶点
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                REAL theta = 2 * M_PI * i / n;
                REAL phi = M_PI * j / m;
                Point<REAL, DIM> p = center + Point<REAL, DIM>(radius * sin(phi) * cos(theta), radius * sin(phi) * sin(theta), radius * cos(phi));
                res.first.push_back(p);
            }
        }
        // 添加面
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                int a = i * m + j;
                int b = i * m + (j + 1) % m;
                int c = ((i + 1) % n) * m + j;
                int d = ((i + 1) % n) * m + (j + 1) % m;
                res.second.push_back({a,b,c});
                res.second.push_back({b,d,c});
            }
        }
        return res;
    }

    /**
     * @brief Get the arrow object
     * 得到一个箭头的Mesh。箭头的起点为start + (end-start)*0.01,终点为end 
     *  箭头为一个圆锥和一个圆柱的组合。圆柱的底面半径为圆柱的1/4
        圆锥的顶点为enod
    * @param radius 圆锥的底面半径
    * @param n 圆锥/圆柱的纬度精细度
     */
    template <typename REAL, int DIM>
    MESH get_arrow(Point<REAL,DIM> start, Point<REAL,DIM> end, REAL radius = 0.01, int n = 10){
        REAL rate = 2;//圆柱半径与圆锥底面半径的比值
        Point<REAL, DIM> dir = end - start;
        Point<REAL, DIM> z = dir;
        lzd_tools::normalize(z);
        Point<REAL, DIM> x = Point<REAL, DIM>(1, 0, 0);
        if (Distance<REAL,DIM>(x,z)<(REAL)0.01) {
            x = Point<REAL, DIM>(0, 1, 0);
        }
        Point<REAL, DIM> y = Point<REAL, DIM>::CrossProduct(z, x);
        lzd_tools::normalize(y);
        x = Point<REAL, DIM>::CrossProduct(y, z);
        lzd_tools::normalize(x);


        end = start + dir * 0.98;
        start = start + dir * 0.02;
        dir = end - start;

        Point<REAL,DIM> cylinder_end = start + dir * 0.9;
        
        MESH res;
        // 圆柱底面顶点
        for(int i=0;i<n;i++){
            REAL theta = 2 * M_PI * i / n;
            Point<REAL, DIM> p = start + radius / rate * (cos(theta) * x + sin(theta) * y);
            res.first.push_back(p);
        }

        // 圆柱顶面顶点
        for(int i=0;i<n;i++){
            REAL theta = 2 * M_PI * i / n;
            Point<REAL, DIM> p = cylinder_end + radius / rate * (cos(theta) * x + sin(theta) * y);
            res.first.push_back(p);
        }
        // 圆柱侧面
        for(int i=0;i<n;i++){
            res.second.push_back({n+i,i,(i+1)%n,});
            res.second.push_back({n+i,(i+1)%n,(i+1)%n+n});
        } 
        int cend = res.first.size();

        // 圆锥底面
        for(int i=0;i<n;i++){
            REAL theta = 2 * M_PI * i / n;
            Point<REAL, DIM> p = cylinder_end + radius * (cos(theta) * x + sin(theta) * y);
            res.first.push_back(p);
        }

        // 圆锥顶点
        res.first.push_back(end);
        int top = res.first.size() - 1;
        
        // 圆锥侧面
        for(int i=0;i<n;i++){
            res.second.push_back({top,i + cend,(i+1)%n+cend});
        }       
        return res;  
    }

    //在op上找到距离op中心最近的点
    template <typename REAL, int DIM>
    Point<REAL, DIM> get_center_on_op(const ORIENTED_POINTS& op){
        Point<REAL, DIM> c,res;
        for(auto& i:op){
            c += i.first;
        }
        c /= op.size();
        // 在op中找到离res最近的点
        REAL min_dis = 1e10;
        for(auto& i:op){
            REAL dist = Distance(i.first, c);
            if(dist < min_dis){
                min_dis = dist;
                res = i.first;
            }
        }        
        return res;
    }

    /**
     * @brief 
     * 翻转mesh 翻转mesh即将mesh中的所有面的顶点顺序翻转
     */
    template <typename REAL, int DIM>
    int revert_mesh(MESH& mesh){
        for(auto& i:mesh.second){
            std::swap(i[0],i[1]);
        }
        return 0;
    }    

    template <typename REAL, int DIM>
    int revert_op(ORIENTED_POINTS& op){
        for(auto& i:op){
            i.second.normal = -i.second.normal;
        }
        return 0;
    }


    template <typename REAL, int DIM>
    int add_topology(MESH& mesh, const MESH& topology){
        int start = mesh.first.size();
        for(auto& i:topology.first){
            mesh.first.push_back(i);
        }
        for(auto i:topology.second){
            std::vector<int> t;
            for(auto j:i){
                t.push_back(j + start);
            }
            mesh.second.push_back(t);
        }
        return 0;
    }

    /**
     * @brief 
     * 删除op中,op.first完全一致的点;其余点的顺序不变
     * @return int 
     */
    template <typename REAL, int DIM>
    int duplicate_filter(const ORIENTED_POINTS& op, ORIENTED_POINTS& res){
        assert(&op != &res);
        std::vector<std::pair<int,Point<REAL, DIM>>> temp;
        for(int i=0;i<op.size();i++){
            temp.push_back({i,op[i].first});
        }
        // 定义比较函数
        auto cmp = [](const std::pair<int,Point<REAL, DIM>>& a, const std::pair<int,Point<REAL, DIM>>& b){
            for(int i=0;i<DIM;i++){
                if(a.second[i] != b.second[i]){
                    return a.second[i] < b.second[i];
                }
            }
            return false;
        };
        auto eq = [](const std::pair<int,Point<REAL, DIM>>& a, const std::pair<int,Point<REAL, DIM>>& b){
            for(int i=0;i<DIM;i++){
                if(a.second[i] != b.second[i]){
                    return false;
                }
            }
            return true;
        };        


        // 根据点的坐标排序
        std::sort(temp.begin(),temp.end(),cmp);
        // 删除重复点
        std::vector<int> flag(op.size(),1);//1表示保留,0表示删除
        int dup_count = 0;
        for(int i=0;i<temp.size();i++){
            if(i+1<temp.size() && eq(temp[i],temp[i+1])){
                flag[temp[i+1].first] = 0;
                dup_count++;
            }
        }
        // 生成res
        res.clear();
        for(int i=0;i<op.size();i++){
            if(flag[i]){
                res.push_back(op[i]);
            }
        }
        return dup_count;
    }    

    template <typename REAL, int DIM>
    struct PointNormalMetric {
        int total_count;
        int angle15_count;
        int angle30_count;
        int angle45_count;
        int angle60_count;
        int angle90_count;
        REAL avg_loss;
        REAL avg_nd_loss;// 保留一下吧与avg_loss区分
        
        // 计算与gt相差大于angle度的点的个数
        // permit_nd: 是否允许一次翻转
        int count_angle_num(const ORIENTED_POINTS& _gt_points_normals, const ORIENTED_POINTS& _points_normals,REAL angle, bool permit_nd = true){
            assert(angle<=90);
            int count = 0, inv_count = 0;
            double inv_angle = 180 - angle;
            for (size_t i = 0; i < _points_normals.size(); i++) {

                double theta = calculate_angle(_gt_points_normals[i].second, _points_normals[i].second) * 180 / M_PI;
                if (theta > angle) {
                    count++;
                }
                if (theta < inv_angle) {
                    inv_count++;
                }
            }
            if (permit_nd) {
                return std::min(count, inv_count);
            }
            return count;
        }

        PointNormalMetric(const ORIENTED_POINTS& op, const ORIENTED_POINTS& gt,bool permit_nd = true) {
            REAL dist = 0;
            assert(op.size() == gt.size());
            total_count = op.size();
            avg_loss = 0;
            for (int i = 0; i < op.size(); i++) {
                for (int d = 0; d < DIM; d++) {
                    dist += op[i].first[d] - gt[i].first[d];
                }
                auto angle = calculate_angle(op[i].second, gt[i].second) / M_PI * 180;
                auto iangle = calculate_angle(op[i].second*-1, gt[i].second) / M_PI * 180;
                avg_loss += angle;
                avg_nd_loss += iangle;
            }
            avg_nd_loss = min(avg_loss, avg_nd_loss) / total_count; 
            avg_loss /= total_count;
            angle15_count = count_angle_num(gt, op, 15, permit_nd);
            angle30_count = count_angle_num(gt, op, 30, permit_nd);
            angle45_count = count_angle_num(gt, op, 45, permit_nd);
            angle60_count = count_angle_num(gt, op, 60, permit_nd);
            angle90_count = count_angle_num(gt, op, 90, permit_nd);
            
            if (dist > 1e5 || dist < -1e5) {
                printf("WARNING:: dist = %f between op and gt\n", dist);
            }
        
        }
        
        nlohmann::json to_json() {
            nlohmann::json res;
            res["total_count"] = total_count;
            res["angle15_count"] = angle15_count;
            res["angle30_count"] = angle30_count;
            res["angle45_count"] = angle45_count;
            res["angle60_count"] = angle60_count;
            res["angle90_count"] = angle90_count;
            res["avg_loss"] = avg_loss;
            res["avg_nd_loss"] = avg_nd_loss;
            return res;
        }

    };
    
    /**
     * @brief 
     * 压缩数组,将数组中的元素压缩为连续的整数
     */
    std::vector<int> squeeze_array(const std::vector<int>& _array, int start = 0);
    
    std::vector<int> read_array(const std::string filename);

 

    // 重新排序op[left,right]，使得op[left:left+k][d] < op[left+k][d] <= op[left+k+1:right][d]
    template <typename REAL, int DIM>
    void quickselect( std::vector<std::pair<Point<REAL, DIM>, int>>& op, int left, int right, int k, int d) {
        if (left < right) {
            // // 在第 d 维度进行划分
            // int pivot_index = partition(op, left, right, d);
            // int rank = pivot_index - left + 1;a
            REAL pivot = op[right].first[d];
            int i = left;
            for (int j = left; j < right; ++j) {
                if (op[j].first[d] < pivot) {
                    std::swap(op[i], op[j]);
                    ++i;
                }
            }
            std::swap(op[i], op[right]);
            int rank = i - left + 1;

            if (rank == k) {
                return;
            }
            else if (rank < k) {
                quickselect(op, i + 1, right, k - rank, d);
            }
            else {
                quickselect(op, left, i - 1, k, d);
            }    
        }
    }

    // 将op中，left到right中的元素，从第d维度开始，递归第二分，直到将op中的元素划分为k个部分
    template <typename REAL, int DIM>
    void _k_select( std::vector<std::pair<Point<REAL, DIM>, int>>& op, int k, int left, int right, int d) {
        if (k == 1) {
            return;
        }
        
        // 递归划分
        int mid = (right - left) / 2;
        quickselect(op, left, right, mid, d);
        int kleft = k / 2, kright = k - kleft;
        d = (d + 1) % DIM;
        _k_select(op, kleft, left, left + mid, d);
        _k_select(op, kright, left + mid + 1, right, d);        
    }

    // 将op中的元素划分为k个部分
    template <typename REAL, int DIM>
    std::vector<int> k_select(const ORIENTED_POINTS& op, int k) {
        std::vector<std::pair<Point<REAL, DIM>, int>> op_with_index;
        for (int i = 0; i < op.size(); ++i) {
            op_with_index.push_back({op[i].first, i});
        }

        _k_select(op_with_index, k, 0, op.size() - 1, 0);
        std::vector<int> res;
        int step = op.size() / k;
        for (int i = 0; i < k; ++i) {
            res.push_back(op_with_index[i * step].second);
        }
        return res;
    }
}

namespace global_var {
    extern std::string data_output_base;

    // enum linked_graph_filp_alg_opt{
	// 		GREEDY,
	// 		BRUTE_FORCE,
	// 		SINGLE_TREE,
	// 		MUTI_TREE,
	// 		OPTIMAL_FILP,
    //         BEST_TREE
	// };
    extern std::map<int,std::string> linked_graph_filp_alg_name;

    extern nlohmann::json config_j;// 控制参数
    extern nlohmann::json metric_j;// 结果
    extern nlohmann::json res_log_j;

    //extern indicators::ProgressBar bar;
    //class Logger {
    //    static nlohmann::json config_j;// 控制参数
    //    static nlohmann::json metric_j;// 效果、日志等
    //};

    extern std::vector<std::vector<unsigned int>> color_map_100;


}


// 求一个一对多的映射的逆映射,给定M0:X->Y,得到映射Mi:Y->X; 需要提供Y的最大值,否则默认Y_max = max({Y|Y->X∈M0}),这可以会导致错误
template <typename REAL, int DIM>
MULTI_MAP inverse_map(MULTI_MAP& map, int MAX_Y = -1){
    // 求B集合中的最大值
    int max = 0;
    if (MAX_Y <= 0) {
        for (auto i : map) {
            for (auto j : i) {
                if (j > max)
                    max = j;
            }
        }
    }
    else {
        max = MAX_Y;
    }
    max += 1;
    // 用于存储逆映射
    MULTI_MAP inverse(max+1);
    for(int i=0;i<map.size();i++){
        for(auto j:map[i]){
            inverse[j].push_back(i);
        }
    }
    return inverse;
}

// 由mesh生成面法向量 overlap_plane 1 10 10 824
template <typename REAL, int DIM>
int gen_face_norm(MESH& mesh, std::vector<Point<REAL,DIM>>& face_norm){
    face_norm.resize(mesh.second.size());
#pragma omp parallel for
    for(int i=0;i<mesh.second.size();i++){
        //if(mesh.second[i].size() != 3){
        //    std::cout << "Error: mesh.second[i].size() != 3" << std::endl;
        //    continue;
        //}
        Point<REAL, DIM> p1 = mesh.first[mesh.second[i][0]];
        Point<REAL, DIM> p2 = mesh.first[mesh.second[i][1]];
        Point<REAL, DIM> p3 = mesh.first[mesh.second[i][2]];
        Point<REAL, DIM> v1 = p2 - p1;
        Point<REAL, DIM> v2 = p3 - p1;
        face_norm[i] = Point<REAL, DIM>::CrossProduct(v1, v2); 
    }
    //if(i<mesh.second.size())
    //    return -1;
    return 0;
}


 /* @brief 
 * @param mesh 待写入的mesh 
 * @param face_norm 面法向量 
 * @param path 输出路径
 * @param flag 用于标记点的颜色,如果flag.size() == 0,则不标记颜色,否则flag.size() == mesh.first.size()
 * @param xform 变换矩阵
 * @return int 
 */
template <typename REAL, int DIM>
int mesh2ply(MESH& mesh, std::vector<Point<REAL,DIM>>& face_norm, const std::string path, const std::vector<int>& flag = std::vector<int>(),const XForm<REAL, DIM+1> xform = XForm<REAL, DIM + 1>().Identity()){
    const int n = mesh.first.size();
    const int m = mesh.second.size();
    std::ofstream out(path);
    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex " << n << std::endl;
    // 如果有flag字段,则将不同的点标记为不同的颜色
    out << "property float x" << std::endl;
    out << "property float y" << std::endl;
    out << "property float z" << std::endl;
    if(flag.size() == n){
        out << "property uchar red" << std::endl;
        out << "property uchar green" << std::endl;
        out << "property uchar blue" << std::endl;
    }


    out << "element face " << m << std::endl;
    out << "property list uchar int vertex_index" << std::endl;
    out << "property float nx" << std::endl;
    out << "property float ny" << std::endl;
    out << "property float nz" << std::endl;
    out << "end_header" << std::endl;
    for (int i = 0; i < n; i++) {
        Point<REAL, DIM> p = xform * mesh.first[i];
        out << p[0] << " " << p[1] << " " << p[2];
        if(flag.size() == n){
            out<< " ";
            if(flag[i] == 0){
                out << "0 0 255" << std::endl;
            }else if(flag[i] == 1){
                out << "255 0 0" << std::endl;
            }else if(flag[i] == 2){
                out << "0 255 0" << std::endl;
            }else{
                out << "255 255 255" << std::endl;
            }
        }else{
            out << std::endl;
        }
    }
    for (int i = 0; i < m; i++) {
        out << "3 " << mesh.second[i][0] << " " << mesh.second[i][1] << " " << mesh.second[i][2] << " ";
        out << face_norm[i][0] << " " << face_norm[i][1] << " " << face_norm[i][2] << std::endl;
    }
    out.close();
    return 0;
}

// 
std::vector<std::vector<int>> flag2color(const std::vector<int>& flag);
/**
 * @brief 
 * 绘制无向点云
 * @param xyz 点云
 * @param flag 用于标记点的颜色,如果flag.size() != xyz.size() 则不标记颜色.本函数会自动根据flag中不同值的个数随机生成颜色
 */
template <typename REAL, int DIM>
int xyz2ply(std::vector<Point<REAL,DIM>>& xyz,const std::string path, const std::vector<int>& flag = std::vector<int>(), const XForm<REAL, DIM+1> xform = XForm<REAL, DIM+1>().Identity()){
    std::ofstream out(path);
    auto colors = flag2color(flag);
    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex " << xyz.size() << std::endl;
    // 如果有flag字段,则将不同的点标记为不同的颜色
    out << "property float x" << std::endl;
    out << "property float y" << std::endl;
    out << "property float z" << std::endl;
    if(flag.size() == xyz.size()){
        out << "property uchar red" << std::endl;
        out << "property uchar green" << std::endl;
        out << "property uchar blue" << std::endl;
    }

    out << "end_header" << std::endl;
    
    for(int i=0;i<xyz.size();i++){
        Point<REAL, DIM> p = xform * xyz[i];
        out << p[0] << " " << p[1] << " " << p[2];
        if(flag.size() == xyz.size()){
            out<< " ";
            int r = colors[i][0];
            int g = colors[i][1];
            int b = colors[i][2];
            out << r << " " << g << " " << b << std::endl;

            // if(flag[i] == 0){
            //     out << "0 0 255" << std::endl;
            // }else if(flag[i] == 1){
            //     out << "255 0 0" << std::endl;
            // }else if(flag[i] == 2){
            //     out << "0 255 0" << std::endl;
            // }else{
            //     out << "255 255 255" << std::endl;
            // }
        }else{
            out << std::endl;
        }
    }

    out.close();
    return 0;
}


int change_args(std::vector<std::string>& args, std::string key, std::string value);

int find_arg(std::vector<std::string>& args, std::string key);

// 计算两个点之间的距离
template <typename REAL, int DIM>
REAL distance(Point<REAL, DIM> p1, Point<REAL, DIM> p2) {
    REAL sum = 0;
    for (int i = 0; i < DIM; i++) {
        sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sqrt(sum);
}

// DBSCAN算法
// 输入:点集,半径,最小点数
// 输出:聚类结果
// 已经弃用 请使用o3d_api中的DBSCAN
template <typename REAL, int DIM>
int DBSCAN(ORIENTED_POINTS& p, std::vector<ORIENTED_POINTS>& res,ORIENTED_POINTS& rest,REAL ep = 0.1, int min_points = 10){
    // 将点集转换为dbscan需要的格式
    std::vector<DBSCAN_joo::Point> joo_points;
    for(auto i:p){
        DBSCAN_joo::Point joo_point = { i.first[0],i.first[1],i.first[2],UNCLASSIFIED };
        joo_points.push_back(joo_point);
    } 
    // 调用dbscan
    DBSCAN_joo::DBSCAN ds(min_points,ep,joo_points);
    ds.run();
    // 将结果转换为需要的格式。对于每一个类别，将其转换为一个点集，然后将这些点集放入res中
    int max = 0;
    for(auto i:joo_points){
        if(i.clusterID > max){
            max = i.clusterID;
        }
    }
    max += 1;
    res.resize(max);
    
    Normal<REAL, DIM> zero_norm;
    for(auto i:joo_points){
        // 将点转换为需要的格式
        Point<REAL, DIM> p(i.x,i.y,i.z);
        Normal<REAL, DIM> n = zero_norm;
        do
        {
            //随机生成一个法向量
            n = Point<REAL, DIM>(rand() % 1001 - 500.0, rand() % 1001 - 500.0, rand() % 1001 - 500.0);   
        } while (n == zero_norm);
        auto pn = std::make_pair(p,n);
        // 将点放入对应的类别中
        if(i.clusterID != UNCLASSIFIED){
            res[i.clusterID].push_back(pn);
        }else{
            rest.push_back(pn);
        }
    }
    return 0;
}

/// <summary>
/// 获取flag值到颜色的映射表。 注意颜色表的长度等于flag中不同值的个数。与flag2color的返回不同
/// </summary>
std::vector<std::vector<int>> MapLable2Color(std::vector<int> flag);


template <typename REAL, int DIM>
Normal<REAL,DIM> get_rand_norm(int seed){
    Normal<REAL, DIM> zero_normal(Point<REAL, DIM>(0, 0, 0));
    do{
        srand(seed);
        Normal<REAL, DIM> n;
        for (int i = 0; i < DIM; i++) {
            n.normal[i] = rand() % 1001 - 500.0;
        }
        if (n == zero_normal) {
            continue;
        }
        else {
            normalize<REAL, DIM>(n);
            return n;
        }
    } while (true);
}

// 递归创建文件夹
// 给定形如"dir1/dir2/dir3"的路径,则创建dir1,dir2,dir3
int rmkdir(const std::string path);

int rmkdir(const std::string base_path,const std::vector<std::string> subpaths);

template <typename REAL>
int linear_regression(const std::vector<REAL> x,const std::vector<REAL> y, REAL& a, REAL& b){
    // 求x,y的均值
    REAL x_mean = 0;
    REAL y_mean = 0;
    for(auto i:x){
        x_mean += i;
    }
    x_mean /= x.size();
    for(auto i:y){
        y_mean += i;
    }
    y_mean /= y.size();
    // 求a,b
    REAL sum1 = 0;
    REAL sum2 = 0;
    for(int i=0;i<x.size();i++){
        sum1 += (x[i] - x_mean) * (y[i] - y_mean);
        sum2 += (x[i] - x_mean) * (x[i] - x_mean);
    }
    a = sum1 / sum2;
    b = y_mean - a * x_mean;
    return 0;
}

// 根据data 辅助判断是否需要改变策略
template<typename REAL>
bool consider_change(std::vector<REAL> data,REAL& a,int consider_count = 5, REAL treshold = 0.0) {
    // 线性拟合
    a = 0;
    if (data.size() < consider_count) {
        return false;
    }
    std::vector<REAL> x, y;
    for (int i = 0; i < consider_count; i++) {
        x.push_back(i);
        y.push_back(data[data.size() - consider_count + i]);
    }
    REAL b;
    linear_regression(x, y, a, b);
    return a < treshold;
}

/**
 * @brief 
 * 由pair<int,int>生成一个唯一的key; 注意，a,b和b,a生成的key是一样的
 * 线程安全
 * @tparam T 
 */
template<typename T>
class UnorderedPairMap {
private:
    std::mutex mtx;
    std::map<long long, T> data;
    static long long hash(int a, int b) {
        if (a > b) {
            std::swap(a, b);
        }
        return (long long)a << 32 | b;
    }
    static std::pair<int, int> unhash(long long key) {
        return std::make_pair(key >> 32, key & 0xffffffff);
    }


public:
    void insert(int a, int b, T value) {
        mtx.lock();
        data[hash(a, b)] = value;
        mtx.unlock();
    }
    T& operator[](std::pair<int, int> p) {
        return data[hash(p.first, p.second)];
    }
    bool has(std::pair<int, int> p) {
        return data.find(hash(p.first, p.second)) != data.end();
    }
    void clear() {
        mtx.lock();
        data.clear();
        mtx.unlock();
    }
    std::vector<std::pair<int, int>> keys() {
        std::vector<std::pair<int, int>> res;
        mtx.lock();
        for (auto i : data) {
            res.push_back(unhash(i.first));
        }
        mtx.unlock();
        return res;
    }
};


