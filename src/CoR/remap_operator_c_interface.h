/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_C_INTERFACE
#define REMAP_OPERATOR_C_INTERFACE


#include "grid_cell_search.h"


extern void get_cell_mask_of_src_grid(long, bool*);
extern void get_cell_mask_of_dst_grid(long, bool*);
extern void get_cell_center_coord_values_of_src_grid(long, double*);
extern void get_cell_center_coord_values_of_dst_grid(long, double*);
extern void get_cell_vertex_coord_values_of_src_grid(long, int*, double*, bool);
extern void get_cell_vertex_coord_values_of_dst_grid(long, int*, double*, bool);
extern void search_cell_in_src_grid(double*, long*, bool);
extern void initialize_computing_remap_weights_of_one_cell();
extern void finalize_computing_remap_weights_of_one_cell();
extern void clear_remap_weight_info_in_sparse_matrix();
extern void add_remap_weights_to_sparse_matrix(long*, long, double*, int, int, bool);
extern double compute_difference_of_two_coord_values(double, double, int);
extern void sort_polygon_vertexes(double, double, double*, double*, long*, int);
extern double compute_vector_angle(double*, double*);
extern void clear_src_grid_cell_visiting_info();
extern bool src_cell_and_dst_cell_have_overlap(long, long);
extern bool have_overlapped_src_cells_for_dst_cell(long); 
extern void compute_intersect_points_of_two_great_arcs_of_sphere_grid(double, double, double, double, double, double,
                                                                       double, double, int&, double*, double*);
extern void compute_common_sub_cell_of_src_cell_and_dst_cell_2D(long, long, int&, double*, double*);
extern double compute_area_of_sphere_cell(int, double*, double*);
extern void compute_cell_bounding_box(int, int, double*, double*);
extern void sort_vertexes_of_sphere_cell(int, double*, double*);
extern bool are_the_same_sphere_points(double, double, double, double);

extern long get_size_of_src_grid();
extern long get_size_of_dst_grid();

extern H2D_grid_cell_search_engine *get_current_grid2D_search_engine(bool);


#endif

