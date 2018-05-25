/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef RUNTIME_REMAP
#define RUNTIME_REMAP

#include "remap_statement_operand.h"
#include "memory_mgt.h"
#include "remap_weight_of_strategy_class.h"
#include "runtime_remapping_weights_mgt.h"
#include <vector>


class Runtime_transfer_algorithm;
class Runtime_remap_algorithm;


struct Operation_for_dynamic_sigma_grid
{
    Field_mem_info *surface_field_of_sigma_grid_src;
    Field_mem_info *surface_field_of_sigma_grid_dst;
    Remap_grid_class *current_3D_sigma_grid_src;
    Remap_grid_class *current_3D_sigma_grid_dst;
    Runtime_remap_algorithm *runtime_remap_algorithm;
};


class Runtime_remap_algorithm
{
    private:       
        std::vector<Operation_for_dynamic_sigma_grid*> operations_for_dynamic_sigma_grid;
        int comp_id;
        Field_mem_info *specified_src_field_instance;
        Field_mem_info *specified_dst_field_instance;
        Field_mem_info *true_src_field_instance;
        Field_mem_info *true_dst_field_instance;
        Runtime_remapping_weights *runtime_remapping_weights;
        bool transform_data_type;
        
        void do_remap(bool);

    public:    
        Runtime_remap_algorithm(Runtime_remapping_weights *, Field_mem_info *, Field_mem_info *, int);
        Runtime_remapping_weights *get_runtime_remapping_weights() { return runtime_remapping_weights; }
        bool run(bool);
        void allocate_src_dst_fields(bool);
        ~Runtime_remap_algorithm();
};

#endif
