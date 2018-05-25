/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_SMOOTH
#define REMAP_OPERATOR_SMOOTH


#include "remap_operator_basis.h"

class Remap_operator_smooth: public Remap_operator_basis
{
    private:
        int num_period;
        void compute_remap_weights_of_one_dst_cell(long);
        
    public:
        Remap_operator_smooth(const char*, int, Remap_grid_class **);
        Remap_operator_smooth() {}
        ~Remap_operator_smooth() {}
        void set_parameter(const char *, const char *);
        int check_parameter(const char *, const char *, char*);
        void calculate_remap_weights();
        void do_remap_values_caculation(double*, double*, int);
        void do_src_decomp_caculation(long*, const long*);
        Remap_operator_basis *duplicate_remap_operator(bool);
        Remap_operator_basis *generate_parallel_remap_operator(Remap_grid_class**, int**);
};


#endif
