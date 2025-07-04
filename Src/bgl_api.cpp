#include "bgl_api.h"
//            Copyright Daniel Trebbien 2010.
// Distributed under the Boost Software License, Version 1.0.
//   (See accompanying file LICENSE_1_0.txt or the copy at
//         http://www.boost.org/LICENSE_1_0.txt)

using namespace bgl_api;




//// A graphic of the min-cut is available at
//// <http://www.boost.org/doc/libs/release/libs/graph/doc/stoer_wagner_imgs/stoer_wagner.cpp.gif>
//int test_bgl()
//{
//    using namespace std;
//
//    // define the 16 edges of the graph. {3, 4} means an undirected edge between
//    // vertices 3 and 4.
//    edge_t edges[] = { { 3, 4 }, { 3, 6 }, { 3, 5 }, { 0, 4 }, { 0, 1 },
//        { 0, 6 }, { 0, 7 }, { 0, 5 }, { 0, 2 }, { 4, 1 }, { 1, 6 }, { 1, 5 },
//        { 6, 7 }, { 7, 5 }, { 5, 2 }, { 3, 4 } };
//
//    // for each of the 16 edges, define the associated edge weight. ws[i] is the
//    // weight for the edge that is described by edges[i].
//    weight_type ws[] = { 0, 3, 1, 3, 1, 2, 6, 1, 8, 1, 1, 80, 2, 1, 1, 4 };
//
//    // construct the graph object. 8 is the number of vertices, which are
//    // numbered from 0 through 7, and 16 is the number of edges.
//    undirected_graph g(edges, edges + 16, ws, 8, 16);
//
//    // define a property map, `parities`, that will store a boolean value for
//    // each vertex. Vertices that have the same parity after
//    // `stoer_wagner_min_cut` runs are on the same side of the min-cut.
//    BOOST_AUTO(parities,
//        boost::make_one_bit_color_map(
//            num_vertices(g), get(boost::vertex_index, g)));
//
//    // run the Stoer-Wagner algorithm to obtain the min-cut weight. `parities`
//    // is also filled in.
//    int w = boost::stoer_wagner_min_cut(
//        g, get(boost::edge_weight, g), boost::parity_map(parities));
//
//    cout << "The min-cut weight of G is " << w << ".\n" << endl;
//    assert(w == 7);
//
//    cout << "One set of vertices consists of:" << endl;
//    size_t i;
//    for (i = 0; i < num_vertices(g); ++i)
//    {
//        if (get(parities, i))
//            cout << i << endl;
//    }
//    cout << endl;
//
//    cout << "The other set of vertices consists of:" << endl;
//    for (i = 0; i < num_vertices(g); ++i)
//    {
//        if (!get(parities, i))
//            cout << i << endl;
//    }
//    cout << endl;
//
//    return EXIT_SUCCESS;
//}

std::vector<std::pair<int,int>> bgl_api::bgl_stoer_wanger_mincut(const UGRAPH& input_graph)
{
    std::vector<edge_t> _es;
    std::vector<weight_type> ws;
    for(bgl_api::VecType i = 0; i<input_graph.size(); i++) {
        for(bgl_api::VecType j = i+1;j<input_graph[i].size();j++){
            assert(input_graph[i][j] == input_graph[j][i]);
            if(input_graph[i][j] > 0){
                _es.push_back({i,j});
                ws.push_back(input_graph[i][j]);
            }
        }
    }
    undirected_graph g(_es.begin(),_es.end(),ws.begin(),input_graph.size(),_es.size());
    //undirected_graph g = _transfer()
    BOOST_AUTO(parities,
        boost::make_one_bit_color_map(
        num_vertices(g),get(boost::vertex_index,g)));

    int w = boost::stoer_wagner_min_cut(
        g,get(boost::edge_weight,g),boost::parity_map(parities));

    std::vector<int> sw_res;
    // cout << "The min-cut weight of G is " << w << ".\n" << endl;
    for(size_t i = 0;i<num_vertices(g);i++){
        sw_res.push_back(get(parities,i));
    }
    std::vector<std::pair<int,int>> cut_edges;
    for(auto ed:boost::make_iterator_range(boost::edges(g))){
        auto s = source(ed,g), t = target(ed,g);
        if(get(parities,s)!=get(parities,t)){
            cut_edges.push_back(std::make_pair<int,int>(s,t));
        }
    }
    //assert(cut_edges.size() == w);
    return cut_edges;
}

std::vector<std::vector<int>> bgl_api::bgl_connected_components(const UGRAPH &input_graph)
{
    std::vector<edge_t> _es;
    std::vector<weight_type> ws;
    for(bgl_api::VecType i = 0;i<input_graph.size();i++){
        for(bgl_api::VecType j = i+1;j<input_graph[i].size();j++){
            assert(input_graph[i][j] == input_graph[j][i]);
            if(input_graph[i][j]>0){
                _es.push_back({i,j});
                ws.push_back(1);
            }
        }
    }
    undirected_graph g(_es.begin(),_es.end(),ws.begin(),input_graph.size(),_es.size());
    std::vector<int> sp_res(num_vertices(g));
    int num = boost::connected_components(g,&sp_res[0]);
    std::vector<std::vector<int>> component(num);
    for(size_t i = 0;i<sp_res.size();i++){
        component[sp_res[i]].push_back(i);
    }
    return component;
}
