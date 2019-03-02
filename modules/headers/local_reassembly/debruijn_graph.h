//
// Created by Kishwar Shafin on 10/19/18.
//

#ifndef HELEN_DEBRUIJN_GRAPH_H
#define HELEN_DEBRUIJN_GRAPH_H

#include "../dataio/bam_handler.h"
#include <iostream>
#include <set>
#include <stack>
#include <algorithm>
#include <string>
using namespace std;

namespace DeBruijnGraph_options {
    static constexpr int MIN_EDGE_SUPPORT = 2;
    static constexpr int MIN_BASE_QUALITY = 15;
    static constexpr int MIN_MAP_QUALITY = 15;
    static constexpr int MAX_ALLOWED_PATHS = 256;
};

class DeBruijnGraph{
public:
    long long region_start;
    long long region_end;

    // nodes
    map<string, int> node_hash_str_to_int;
    map<int, string> node_hash_int_to_str;
    map<int, bool> good_nodes;
    int current_hash_value;

    // out nodes
    map<int, vector<int> > out_nodes;
    map<int, vector<int> > in_nodes;

    // source and sink nodes
    int _source_node;
    int _sink_node;

    // dfs colors
    map<int, bool> visit_color;
    map<int, bool> stack_color;

    // edge_a -> edge_b (weight, ref_or_not)
    map< pair<int, int>, pair<int, bool> > edges;

    int get_hash(string kmer);
    void add_reference_path(string reference_sequence, int kmer_size);
    void add_read_to_graph(type_read &read, int kmer_size);
    void add_edge(int node_a, int node_b, bool is_ref);
    void remove_edge(int node_a, int node_b);
    void remove_node(int node_a);
    bool detect_cycle(int v);
    bool is_cyclic();
    void prune_graph();
    void graph_trevarsal();
    bool check_if_base_ok(char base);
    vector<string> get_haplotypes();

    DeBruijnGraph(long long region_start, long long region_end);
    static pair<int, int> find_min_k_from_ref(string reference_sequence, int min_kmer, int max_kmer, int step_kmer);
    vector<string> generate_haplotypes(string reference, vector <type_read> reads, int kmer_size);

};
#endif //FRIDAY_DEBRUIJN_GRAPH_H
