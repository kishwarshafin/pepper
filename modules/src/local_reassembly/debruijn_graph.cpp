//
// Created by Kishwar Shafin on 10/19/18.
//

#include "../../headers/local_reassembly/debruijn_graph.h"

pair<int, int> DeBruijnGraph::find_min_k_from_ref(string reference_sequence, int min_kmer, int max_kmer, int step_kmer) {
    int found_min_k = -1;
    int found_max_k = min(max_kmer, (int)(reference_sequence.size()-1));

    for(int k_mer = min_kmer; k_mer <= max_kmer; k_mer += step_kmer) {
        bool has_cycle = false;

        set <string> kmer_set;

        for (int i=0; i < reference_sequence.size() - k_mer + 1; i++) {
            string kmer_str = reference_sequence.substr(i, k_mer);
            if(kmer_set.find(kmer_str) != kmer_set.end()) {
                has_cycle = true;
                break;
            } else{
                kmer_set.insert(kmer_str);
            }
        }
        if(!has_cycle) {
            found_min_k = k_mer;
            break;
        }
    }

    return make_pair(found_min_k, found_max_k);
}

DeBruijnGraph::DeBruijnGraph(long long region_start, long long region_end) {
    current_hash_value = 1;
    this->region_start = region_start;
    this->region_end = region_end;
}

int DeBruijnGraph::get_hash(string kmer){
    if(node_hash_str_to_int.find(kmer) != node_hash_str_to_int.end()) {
        return node_hash_str_to_int[kmer];
    }
    node_hash_str_to_int[kmer] = current_hash_value;
    node_hash_int_to_str[current_hash_value] = kmer;

    current_hash_value += 1;

    return node_hash_str_to_int[kmer];
}

void DeBruijnGraph::add_edge(int node_a, int node_b, bool is_ref) {
    pair<int, int> node_pair = make_pair(node_a, node_b);
    if(edges.find(node_pair) == edges.end()){
        this->out_nodes[node_a].push_back(node_b);
        this->in_nodes[node_b].push_back(node_a);
        pair<int, bool> edge_values = make_pair(is_ref ? 0 : 1, is_ref);
        edges[node_pair] = edge_values;
    } else {
        edges[node_pair].first = edges[node_pair].first + 1;
    }
}


void DeBruijnGraph::remove_node(int node_a) {
    for(auto &v : in_nodes[node_a]) {
        edges.erase(make_pair(v, node_a));
        out_nodes[v].erase(remove(out_nodes[v].begin(), out_nodes[v].end(), node_a), out_nodes[v].end());
    }
    in_nodes.erase(node_a);
}


void DeBruijnGraph::remove_edge(int node_a, int node_b) {
    out_nodes[node_a].erase(remove(out_nodes[node_a].begin(), out_nodes[node_a].end(), node_b), out_nodes[node_a].end());
    in_nodes[node_b].erase(remove(in_nodes[node_b].begin(), in_nodes[node_b].end(), node_b), in_nodes[node_b].end());
    edges.erase(make_pair(node_a, node_b));
}

bool DeBruijnGraph::check_if_base_ok(char base){
    if(base=='A' || base=='C' || base=='G' || base=='T' || base=='a' || base == 'c' || base == 'g' || base == 't')
        return true;
    return false;
}


void DeBruijnGraph::add_read_to_graph(type_read &read, int kmer_size) {
    long long read_start_index = max((long long) 0, (region_start - read.pos - 1));

    int current_position = read_start_index;
    int current_bad_vector_index = 0;

    while(current_position < read.sequence.length() - kmer_size) {
        int next_bad_position = read.bad_indicies[current_bad_vector_index] - 1;
        if(next_bad_position - current_position >= kmer_size) {
            int previous_node = get_hash(read.sequence.substr(current_position, kmer_size));
            for(int i=current_position + 1; i < next_bad_position - kmer_size + 1; i++) {
                int current_node = get_hash(read.sequence.substr(i, kmer_size));
                add_edge(previous_node, current_node, false);
                previous_node = current_node;
            }
        }

        current_position = next_bad_position + 1;
        current_bad_vector_index += 1;
    }
}

void DeBruijnGraph::add_reference_path(string reference_sequence, int kmer_size) {
    int previous_node = get_hash(reference_sequence.substr(0, kmer_size));
    _source_node = previous_node;

    for(int i=1; i < reference_sequence.length() - kmer_size + 1; i++) {
        int current_node = get_hash(reference_sequence.substr(i, kmer_size));
        add_edge(previous_node, current_node, true);
        previous_node = current_node;
    }

    _sink_node = previous_node;
}

bool DeBruijnGraph::detect_cycle(int v) {
    int WHITE = 0;
    int GRAY = 1;
    int BLACK = 2;

    int ENTER = 0;
    int EXIT = 1;

    vector<int> state(current_hash_value + 1, WHITE);
    stack< pair<int, int> > _stack;
    _stack.push(make_pair(ENTER, _source_node));

    while(!_stack.empty()) {
        pair<int, int> current_node = _stack.top();
        _stack.pop();
        if(current_node.first == EXIT)
            state[current_node.second] = BLACK;
        else {
            state[current_node.second] = GRAY;
            _stack.push(make_pair(EXIT, current_node.second));
            for(auto next_v: out_nodes[current_node.second]) {
                if(state[next_v] == GRAY) return true;
                else if (state[next_v] == WHITE) _stack.push(make_pair(ENTER, next_v));
            }
        }
    }
    return false;
}

bool DeBruijnGraph::is_cyclic() {
    visit_color.clear();
    stack_color.clear();

    return detect_cycle(_source_node);
}


void DeBruijnGraph::graph_trevarsal() {
    set<int> fw_visited_nodes;
    set<int> bk_visited_nodes;

    stack<int> st;
    st.push(_source_node);
    visit_color.clear();

    while(!st.empty()) {
        int current_node = st.top();
        st.pop();
        visit_color[current_node] = true;
        fw_visited_nodes.insert(current_node);

        if(current_node == _sink_node)
            continue;

        for(auto &node: out_nodes[current_node]) {
            if(visit_color.find(node) == visit_color.end())
                st.push(node);
        }
    }

    stack<int> st_back;
    st_back.push(_sink_node);
    visit_color.clear();

    while(!st_back.empty()) {
        int current_node = st_back.top();
        st_back.pop();
        visit_color[current_node] = true;
        bk_visited_nodes.insert(current_node);

        if(current_node == _source_node)
            continue;

        for(auto &node: in_nodes[current_node]) {
            if(visit_color.find(node) == visit_color.end()){
                st_back.push(node);
            }
        }
    }

    vector<int> all_legit_nodes;
    set_intersection( fw_visited_nodes.begin(), fw_visited_nodes.end(),
                      bk_visited_nodes.begin(), bk_visited_nodes.end(),
                      back_inserter( all_legit_nodes )  );

    for(auto &node: all_legit_nodes) {
        good_nodes[node] = true;
    }
}

void DeBruijnGraph::prune_graph() {
    for (auto const& edge : edges) {
        int _weight = edge.second.first;
        bool _is_ref = edge.second.second;

        if(_weight < DeBruijnGraph_options::MIN_EDGE_SUPPORT && !_is_ref) {
            remove_edge(edge.first.first, edge.first.second);
        }
    }
    graph_trevarsal();
}

vector<string> DeBruijnGraph::get_haplotypes() {
    vector<string> finished_haplotyes;

    vector< vector<int> > finished_paths;
    stack< vector<int> > running_paths;

    running_paths.push({_source_node});

    while(!running_paths.empty()) {
        if(finished_paths.size() + running_paths.size() > DeBruijnGraph_options::MAX_ALLOWED_PATHS) {
            return {};
        }

        vector<int> current_path = running_paths.top();
        running_paths.pop();

        int running_node = current_path.back();
        for(auto &next_node:out_nodes[running_node]) {
            if(good_nodes.find(next_node) == good_nodes.end()) continue;
            vector<int> new_path(current_path.begin(), current_path.end());
            new_path.push_back(next_node);

            if(next_node == _sink_node) {
                finished_paths.push_back(new_path);
            } else {
                running_paths.push(new_path);
            }
        }
    }

    for(auto &path: finished_paths) {
        string haplotype;
        for(int i=0; i < path.size(); i++) {
            if(i==0) haplotype = node_hash_int_to_str[path[i]];
            else haplotype += node_hash_int_to_str[path[i]].back();
        }
        finished_haplotyes.push_back(haplotype);
    }

    return finished_haplotyes;
}

vector<string> DeBruijnGraph::generate_haplotypes(string reference, vector <type_read> reads, int kmer_size) {
    add_reference_path(reference, kmer_size);
    for(auto &read: reads) {
        if(read.mapping_quality < DeBruijnGraph_options::MIN_MAP_QUALITY) continue;
        if(read.pos_end < region_start) continue;
        if(read.pos > region_end) continue;
        add_read_to_graph(read, kmer_size);
    }
    if(is_cyclic()) {
        return {};
    } else {
        prune_graph();
        return get_haplotypes();
    }
}