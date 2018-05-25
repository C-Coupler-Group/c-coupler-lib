/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef COR_CPL_INTERFACE
#define COR_CPL_INTERFACE

#include "cor_global_data.h"


#define GRID_LATS_GF        "n_grid_lats"
#define GRID_LONS_GF        "n_grid_lons"
#define LAT_GF                        "lat"    
#define LON_GF                    "lon"
#define AREA_GF                    "arear"
#define MASK_GF                    "mask"
#define SURFACE_FIELD_GF                    "surface_field"


extern long cpl_get_grid_size(const char*);
extern double *cpl_get_sphere_grid_center_fields(const char*, const char*);
extern double *cpl_get_sphere_grid_area(const char*);
extern bool *cpl_get_sphere_grid_mask(const char*);
extern long cpl_get_sphere_grid_subgrid_size(const char*, const char*);
extern Remap_grid_data_class *cpl_duplicate_field_with_double_data_type(Remap_grid_data_class*);
extern void cpl_check_remap_weights_format(Remap_weight_of_strategy_class*);

#endif