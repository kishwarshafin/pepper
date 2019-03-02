//
// Created by Kishwar Shafin on 7/5/18.
//
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
using namespace std;

#ifndef FRIDAY_CPP_BED_HANDLER_H
#define FRIDAY_CPP_BED_HANDLER_H

struct bed_interval {
    long long start_pos;
    long long end_pos;
    bed_interval(long long st, long long en) {
        start_pos = st;
        end_pos = en;
    }
};

class bed_handler {
    public:
        bed_handler(string path);
        ~bed_handler();
        void read_bed_file();

        map < string, vector<bed_interval> > bed_interval_map;
        vector <string> chromosome_name_set;
    private:
        string bed_file_path;
};

#endif //FRIDAY_CPP_BED_HANDLER_H
