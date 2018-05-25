/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_1D_BASIS
#define REMAP_OPERATOR_1D_BASIS


#include "remap_operator_basis.h"


extern double *common_buffer_for_1D_remap_operator;
extern int size_common_buffer_for_1D_remap_operator;


class Remap_operator_1D_basis: public Remap_operator_basis
{
    protected:
        bool periodic;
        bool set_periodic;
        bool set_period;
        bool set_enable_extrapolation;
        bool use_logarithmic_coordinate;
        bool set_use_logarithmic_coordinate;
        bool use_exponent_coordinate;
        bool set_use_exponent_coordinate;
        int extrapolation_approach;
        double period;
        double *coord_values_src;
        double *coord_values_dst;
        double *packed_data_values_src;
        int num_useful_src_cells;
        int array_size_src;
        int *useful_src_cells_global_index;
        int *src_cell_index_left;
        int *src_cell_index_right;

        double base_field_value;

        void initialize_1D_remap_operator();
        void copy_1D_remap_operator_info(Remap_operator_1D_basis*);
        void allocate_1D_remap_operator_common_arrays_space();
        void search_src_cells_around_dst_cell(double, int, int, int&, int&);
        void search_src_cells_around_dst_cell_recursively(double, int, int,  int&, int&);
        void set_common_parameter(const char*, const char*);
        int check_common_parameter(const char*, const char*, char*);
        void calculate_dst_src_mapping_info();
        void preprocess_field_value(double*);
        void postprocess_field_value(double*);

    public:
        Remap_operator_1D_basis() {}
        Remap_operator_1D_basis(const char*, const char*, int, Remap_grid_class **);
        virtual ~Remap_operator_1D_basis();
        virtual void set_parameter(const char*, const char*) = 0;
        virtual int check_parameter(const char*, const char*, char*) = 0;
        virtual void do_remap_values_caculation(double*, double*, int) = 0;
        virtual void do_src_decomp_caculation(long*, const long*) = 0;
        virtual void calculate_remap_weights() = 0;
        virtual Remap_operator_basis *duplicate_remap_operator(bool) = 0;
        virtual Remap_operator_basis *generate_parallel_remap_operator(Remap_grid_class**, int**) = 0;
        virtual void compute_remap_weights_of_one_dst_cell(long) = 0;
};


#endif

