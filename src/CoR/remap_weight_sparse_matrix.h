/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_WEIGHT_SPARSE_MATRIX_H
#define REMAP_WEIGHT_SPARSE_MATRIX_H


#include "remap_grid_class.h"


class Remap_operator_basis;


class Remap_weight_sparse_matrix
{
    private:
        Remap_operator_basis *remap_operator;
        long *cells_indexes_src;
        long *cells_indexes_dst;
        double *weight_values;
        long *remaped_dst_cells_indexes;
        long weight_arrays_size;
        long num_weights;
        long remaped_dst_cells_indexes_array_size;
        long num_remaped_dst_cells_indexes;
        
    public:
        Remap_weight_sparse_matrix(Remap_operator_basis*);
        Remap_weight_sparse_matrix(Remap_operator_basis*, long, long*, long*, double*, long, long*);
        ~Remap_weight_sparse_matrix();
        void clear_weights_info();
        void add_weights(long*, long, double*, int, bool);
        void get_weight(long*, long*, double*, int);
        void remap_values(double*, double*, int);
        void calc_src_decomp(long*, const long*);
        Remap_weight_sparse_matrix *duplicate_remap_weight_of_sparse_matrix();
        Remap_weight_sparse_matrix *generate_parallel_remap_weight_of_sparse_matrix(Remap_grid_class **, int **);
        Remap_operator_basis *get_remap_operator() { return remap_operator; }
        long get_num_weights() { return num_weights; }
        long *get_indexes_src_grid() { return cells_indexes_src; }
        long *get_indexes_dst_grid() { return cells_indexes_dst; }
        long get_num_remaped_dst_cells_indexes() { return num_remaped_dst_cells_indexes; }
        long *get_remaped_dst_cells_indexes() { return remaped_dst_cells_indexes; }
        double *get_weight_values() { return weight_values; }
        void compare_to_another_sparse_matrix(Remap_weight_sparse_matrix*);
        void print();
		Remap_weight_sparse_matrix *gather(int);
};


#endif
