/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef PERFORMANCE_TIMING_MGT_H
#define PERFORMANCE_TIMING_MGT_H


#include <vector>


#define TIMING_TYPE_COMMUNICATION         1
#define TIMING_TYPE_IO                    2
#define TIMING_TYPE_COMPUTATION           3

#define TIMING_COMMUNICATION_SEND_WAIT   11
#define TIMING_COMMUNICATION_RECV_WAIT   12
#define TIMING_COMMUNICATION_SENDRECV    13
#define TIMING_COMMUNICATION_SEND_QUERRY 14
#define TIMING_COMMUNICATION_RECV_QUERRY 15
#define TIMING_COMMUNICATION_SEND        16
#define TIMING_COMMUNICATION_RECV        17


#define TIMING_IO_INPUT                  21
#define TIMING_IO_OUTPUT                 22
#define TIMING_IO_RESTART                23

#define TIMING_COMPUTATION_ALL           31
#define TIMING_COMPUTATION_H2D_REMAP     32
#define TIMING_COMPUTATION_V1D_REMAP     33
#define TIMING_COMPUTATION_V1D_WEIGHT    34
#define TIMING_COMPUTATION_V1D_COORD     35



class Performance_timing_unit
{
    private:
        int unit_type;
        int unit_behavior;
        int unit_int_keyword;
        char unit_char_keyword[256];
        double previous_time;
        double total_time;
        int comp_id;

        void check_timing_unit(int, int, int, const char*);

    public:
        Performance_timing_unit(int, int, int, int, const char*);
        ~Performance_timing_unit(){}
        void timing_start();
        void timing_stop();
        void timing_output();
        bool match_timing_unit(int, int, int, const char*);
        void timing_add(double time_inc) { total_time += time_inc; }
        void timing_reset() { total_time = 0.0; }
};


class Performance_timing_mgt
{
    private:
        std::vector<Performance_timing_unit*> performance_timing_units;
        int search_timing_unit(int, int, int, const char*);
        int comp_id;

    public: 
        Performance_timing_mgt(int comp_id) { this->comp_id = comp_id; }
        ~Performance_timing_mgt();
        void performance_timing_start(int, int, int, const char*);
        void performance_timing_stop(int, int, int, const char*);
        void performance_timing_add(int, int, int, const char*, double);
        void performance_timing_output();
        void performance_timing_reset();
};

#endif
