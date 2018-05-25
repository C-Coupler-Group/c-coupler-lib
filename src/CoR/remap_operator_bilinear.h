/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_BILINEAR
#define REMAP_OPERATOR_BILINEAR


#include "remap_operator_basis.h"

class Remap_operator_bilinear: public Remap_operator_basis
{
    private:
        int max_num_found_nearest_points;
        double *found_nearest_points_distance;
        long *found_nearest_points_src_indexes;
        double *weigt_values_of_one_dst_cell;
        int num_nearest_points;
        double num_power;
        double iterative_threshold_distance;

        void compute_remap_weights_of_one_dst_cell(long);
        int search_nearnest_src_points_for_bilinear(double*, long, double&, double&);
        int search_at_least_16_nearnest_src_points_for_bilinear(double*, long, double&, double&);
        int compute_quadrant_of_src_point(double*, double*);
        bool get_near_optimal_bilinear_box(double*, int, long*);
        bool get_near_optimal_bilinear_box_recursively(double**, long**, const int*, int*, const double*, long*, int);
        bool are_three_points_on_the_same_line(double*, double*);
        void solve_two_bilinear_ratios(long*, double*, double&, double&);
        double compute_cross_product_of_counter_lines(double*, double*);
        void bilinear_ratios_solution1(double*, double*, double*, double&, double&);
        void bilinear_ratios_solution2(double*, double*, double*, double&, double&);
        void bilinear_one_ratio_solution_of_quadratic_equation(double*, double*, double*, double&);
        

    public:
        Remap_operator_bilinear(const char*, int, Remap_grid_class **);
        Remap_operator_bilinear();
        ~Remap_operator_bilinear();
        void set_parameter(const char *, const char *);
        int check_parameter(const char*, const char*, char*);
        void calculate_remap_weights();
        void do_remap_values_caculation(double*, double*, int);
        void do_src_decomp_caculation(long*, const long*);
        Remap_operator_basis *duplicate_remap_operator(bool);
        Remap_operator_basis *generate_parallel_remap_operator(Remap_grid_class**, int**);
};


#endif
