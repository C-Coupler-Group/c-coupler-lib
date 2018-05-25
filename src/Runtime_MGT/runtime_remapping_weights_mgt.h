/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/



#ifndef RUNTIME_REMAPPING_WEIGHTS_MGT_H
#define RUNTIME_REMAPPING_WEIGHTS_MGT_H


#include "remapping_configuration_mgt.h"
#include <vector>
#include "remap_strategy_class.h"
#include "original_grid_mgt.h"
#include "decomp_info_mgt.h"
#include "memory_mgt.h"
#include "remap_weight_of_strategy_class.h"
#include "compset_communicators_info_mgt.h"


class Runtime_remapping_weights
{
    private:
        const char *src_comp_full_name;
        const char *dst_comp_full_name;
        Original_grid_info *src_original_grid;
        Original_grid_info *dst_original_grid;
        Decomp_info *src_decomp_info;
        Decomp_info *dst_decomp_info;
        Remapping_setting *remapping_setting;
        Remap_strategy_class *remapping_strategy;
        Remap_weight_of_strategy_class *sequential_remapping_weights;
        Remap_weight_of_strategy_class *parallel_remapping_weights;
        Field_mem_info *intermediate_V3D_grid_bottom_field;
        Remap_weight_of_operator_class *dynamic_V1D_remap_weight_of_operator;
        Remap_grid_class *runtime_V1D_remap_grid_src;
        Remap_grid_class *runtime_V1D_remap_grid_dst;
        double *src_H2D_grid_area;
        double *dst_H2D_grid_area;
        int size_src_H2D_grid_area;
        int size_dst_H2D_grid_area;

        void generate_parallel_remapping_weights();
        
    public:
        Runtime_remapping_weights(const char*, const char*, Original_grid_info *, Original_grid_info *, Remapping_setting *, Decomp_info*);
        Runtime_remapping_weights();
        ~Runtime_remapping_weights();
        Remap_weight_of_strategy_class *get_sequential_remapping_weights() { return sequential_remapping_weights; }
        Remap_weight_of_strategy_class *get_parallel_remapping_weights() { return parallel_remapping_weights; }
        Original_grid_info *get_src_original_grid() { return src_original_grid; }
        Original_grid_info *get_dst_original_grid() { return dst_original_grid; }
        Decomp_info *get_src_decomp_info() { return src_decomp_info; }
        Decomp_info *get_dst_decomp_info() { return dst_decomp_info; }
        bool match_requirements(const char*, const char*, Original_grid_info *, Original_grid_info *, Remapping_setting *, Decomp_info*);
        Field_mem_info *allocate_intermediate_V3D_grid_bottom_field();
        void renew_dynamic_V1D_remapping_weights();
        void set_H2D_grids_area(const double*, const double*, long, long);
        double *get_src_H2D_grid_area() { return src_H2D_grid_area; }
        double *get_dst_H2D_grid_area() { return dst_H2D_grid_area; }
};


class Runtime_remapping_weights_mgt
{
    private:
        std::vector<Runtime_remapping_weights*> runtime_remapping_weights;

    public:
        Runtime_remapping_weights_mgt() {}
        ~Runtime_remapping_weights_mgt();
        Runtime_remapping_weights *search_or_generate_runtime_remapping_weights(const char *, const char*, Original_grid_info *, Original_grid_info *, Remapping_setting *, Decomp_info*);
        void transfer_runtime_remapping_weights(Runtime_remapping_weights *, Runtime_remapping_weights **, Comp_comm_group_mgt_node *, Comp_comm_group_mgt_node *);
};


#endif

