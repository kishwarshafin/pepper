//
// Created by Kishwar Shafin on 7/5/18.
//

#include "../../headers/dataio/bed_handler.h"

bed_handler::bed_handler(string path) {
    this->bed_file_path = path;
}

void bed_handler::read_bed_file() {
    string line;
    ifstream bed_file (this->bed_file_path);
    if (bed_file.is_open()) {
        while ( getline (bed_file,line) )
        {
            istringstream iss(line);
            string chromosome_name;
            long long start_pos, end_pos;
            if (!(iss >> chromosome_name >> start_pos >> end_pos)) { break; }
            bed_interval in(start_pos, end_pos);
            bed_interval_map[chromosome_name].push_back(in);
            if(find(chromosome_name_set.begin(), chromosome_name_set.end(), chromosome_name) == chromosome_name_set.end())
                chromosome_name_set.push_back(chromosome_name);
        }
        bed_file.close();
    }
}

bed_handler::~bed_handler() {
}