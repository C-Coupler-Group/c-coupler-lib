/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_COMMON_UTILS
#define REMAP_COMMON_UTILS


#define DEGREE_TO_RADIAN(data)    (data*PI/180)
#define RADIAN_TO_DEGREE(data)    (data*180/PI)



extern void match_degree_values(double &, double &);
extern bool words_are_the_same(const char*, const char*);
extern bool is_point_in_2D_cell(double, double, double*, double*, int, bool, bool, bool);
extern double compute_three_2D_points_cross_product(double, double, double, double, double, double, bool, bool);
extern double compute_three_3D_points_cross_product(double, double, double, double, double, double, double, double, double);
extern void rotate_sphere_coordinate(double, double, double&, double&);
extern void get_3D_cartesian_coord_of_sphere_coord(double&, double&, double&, double, double);


#endif

