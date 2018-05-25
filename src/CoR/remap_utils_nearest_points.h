/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_UTILS_NEAREST_POINTS_H
#define REMAP_UTILS_NEAREST_POINTS_H


extern void compute_dist_remap_weights_of_one_dst_cell(long, int, double, double*, double*, long*, double*, bool, bool);
extern double calculate_distance_of_two_points_2D(double, double, double, double, bool);


#endif

