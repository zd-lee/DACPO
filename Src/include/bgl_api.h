/*
 * @Author: ZHUODONG LI 1989218453@qq.com
 * @Date: 2024-04-04 14:43:11
 * @LastEditors: ZHUODONG LI 1989218453@qq.com
 * @LastEditTime: 2025-06-25 00:37:22
 * @FilePath: \ipsr_explore\Src\include\bgl_api.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#pragma once
#include <vector>
#include <algorithm>

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/one_bit_color_map.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/config.hpp>

namespace bgl_api{


    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
        boost::no_property, boost::property< boost::edge_weight_t, int > >
        undirected_graph;
    typedef boost::property_map< undirected_graph, boost::edge_weight_t >::type
        weight_map_type;
    typedef boost::property_traits< weight_map_type >::value_type weight_type;


    typedef std::vector<std::vector<int>> UGRAPH;// 无向图 0表示无边
    typedef int VecType;
    struct edge_t
    {
        VecType first;
        VecType second;
    };

    // 适用于bgl的边集
    template<typename weight_type>
    struct EdgeSet {
        std::vector<edge_t> edge;
        std::vector<weight_type> edge_weight;
        EdgeSet(std::vector<edge_t> edge, std::vector<weight_type> edge_weight) : edge(edge), edge_weight(edge_weight) {}
        EdgeSet() {}
    };
    
    /**
     * @brief 
     * 对input_graph进行切分
     * @param input_graph 无向图，input_graph[i][j] == 0表示ij之间无边
     * @return 切掉的边
     */
    std::vector<std::pair<VecType,VecType>> bgl_stoer_wanger_mincut(const UGRAPH& input_graph); 

    std::vector<std::vector<int>> bgl_connected_components(const UGRAPH& input_graph);

    template<typename weight_type>
    void bgl_stoer_wanger_mincut(
        const EdgeSet<weight_type>& input_graph,
        weight_type& w,
        std::vector<int>& belonging,
        std::vector<std::pair<int,int>>& cut_edges
    ){
        undirected_graph g(input_graph.edge.begin(),input_graph.edge.end(),input_graph.edge_weight.begin(),input_graph.edge.size());
        BOOST_AUTO(parities,
            boost::make_one_bit_color_map(
            num_vertices(g),get(boost::vertex_index,g)));
        
        w = boost::stoer_wagner_min_cut(
            g,get(boost::edge_weight,g),boost::parity_map(parities));

        for(size_t i = 0;i<num_vertices(g);i++){
            belonging.push_back(get(parities,i));
        }
        for(auto ed:boost::make_iterator_range(boost::edges(g))){
            auto s = source(ed,g), t = target(ed,g);
            if(get(parities,s)!=get(parities,t)){
                cut_edges.push_back(std::make_pair<int,int>(s,t));
            }
        }
    }
}