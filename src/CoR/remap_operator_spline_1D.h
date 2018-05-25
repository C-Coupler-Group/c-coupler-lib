/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_SPLINE_1D
#define REMAP_OPERATOR_SPLINE_1D


#include "remap_operator_1D_basis.h"


class Remap_operator_spline_1D: public Remap_operator_1D_basis
{
    private:
        bool set_keep_monotonicity;
        bool keep_monotonicity;
        double *data_in_monotonicity_range;
        int *dst_cell_indexes_in_monotonicity_ranges;
        double *array_alpha;
        double *array_mu;
        double *array_lambda;
        double *array_h;
        double *array_d;
        double *temp_array_row;
        double *temp_array_column;
        double *final_factor1;
        double *final_factor2;
        double *final_factor3;
        double *final_factor4;
        double *final_factor5;

        void solve_aperiodic_tridiagonal_system(double*, double*, double*, double*, int);
        void solve_periodic_tridiagonal_system(double*, double*, double*, double*, int);
        void compute_remap_weights_of_one_dst_cell(long);
        void allocate_local_arrays();

    public:
        Remap_operator_spline_1D() {}
        Remap_operator_spline_1D(const char*, int, Remap_grid_class **);
        ~Remap_operator_spline_1D();
        void set_parameter(const char *, const char *);
        int check_parameter(const char *, const char *, char *);
        void calculate_remap_weights();
        void do_remap_values_caculation(double*, double*, int);
        void do_src_decomp_caculation(long*, const long*);
        Remap_operator_basis *duplicate_remap_operator(bool);
        Remap_operator_basis *generate_parallel_remap_operator(Remap_grid_class**, int**);
};


#endif

