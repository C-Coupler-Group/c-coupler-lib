/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_DISTWGT
#define REMAP_OPERATOR_DISTWGT


#include "remap_operator_basis.h"


class Remap_operator_distwgt: public Remap_operator_basis
{
    private:
        double num_power;
        int num_nearest_points;
        double *found_nearest_points_distance;
        long *found_nearest_points_src_indexes;
        double *weigt_values_of_one_dst_cell;
        double threshold_distance;

        void compute_remap_weights_of_one_dst_cell(long);

    public:
        Remap_operator_distwgt(const char*, int, Remap_grid_class **);
        Remap_operator_distwgt();
        ~Remap_operator_distwgt();
        void set_parameter(const char*, const char*);
        int check_parameter(const char*, const char*, char*);
        void calculate_remap_weights();
        void do_remap_values_caculation(double*, double*, int);
        void do_src_decomp_caculation(long*, const long*);
        Remap_operator_basis *duplicate_remap_operator(bool);
        Remap_operator_basis *generate_parallel_remap_operator(Remap_grid_class**, int**);
};


#endif
