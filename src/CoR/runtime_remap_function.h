/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef RUNTIME_REMAP_FUNCTION
#define RUNTIME_REMAP_FUNCTION


#include "remap_grid_class.h"
#include "remap_operator_basis.h"
#include "remap_grid_data_class.h"
#include "remap_operator_grid.h"
#include "remap_weight_of_strategy_class.h"


class Runtime_remap_function
{
    private:
        Remap_grid_class *interchanged_grid_src;
        Remap_grid_class *interchanged_grid_dst;
        Remap_grid_class *remap_operator_runtime_grid_src;
        Remap_grid_class *remap_operator_runtime_grid_dst;
        Remap_operator_basis *runtime_remap_operator;
        Remap_grid_data_class *remap_field_data_src;
        Remap_grid_data_class *remap_field_data_dst;
        Remap_grid_data_class *remap_field_data_redundant_mark_field_src;
        Remap_weight_of_strategy_class *remap_weight_of_strategy;
        Remap_weight_of_operator_instance_class *last_remap_weight_of_operator_instance;
        int num_sized_grids_of_interchanged_grid;
        int num_leaf_grids_of_remap_operator_grid_src;
        int num_leaf_grids_of_remap_operator_grid_dst;
        Remap_grid_class *sized_grids_of_interchanged_grid[256];
        Remap_grid_class *leaf_grids_of_remap_operator_grid_src[256];
        Remap_grid_class *leaf_grids_of_remap_operator_grid_dst[256];
        long num_remapping_times;
        long last_remapping_time_iter;
        Remap_operator_grid *runtime_remap_operator_grid_src;
        Remap_operator_grid *runtime_remap_operator_grid_dst;
        bool *last_mask_values_src;
        bool *last_mask_values_dst;
        bool *current_mask_values_src;
        bool *current_mask_values_dst;
        bool *last_redundant_mark_src;
        bool *current_redundant_mark_src;

        void extract_runtime_field(Remap_grid_class *, Remap_grid_data_class*, Remap_grid_data_class*, long);
        bool check_mask_values_status(bool*, bool*, long);
        void check_dimension_order_of_grid_field(Remap_grid_data_class*, Remap_grid_class*);
        
    public:
        Runtime_remap_function(Remap_grid_class*, Remap_grid_class*, Remap_grid_class*, Remap_grid_class*, Remap_operator_basis*, Remap_grid_data_class*, Remap_grid_data_class*, Remap_weight_of_strategy_class*, const char*);
        void calculate_static_remapping_weights(long, const char*, int, bool);
        ~Runtime_remap_function();
};


#endif

