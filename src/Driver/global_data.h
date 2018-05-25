/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef GLOBAL_DATA
#define GLOBAL_DATA

#include "compset_communicators_info_mgt.h"
#include "decomp_info_mgt.h"
#include "field_info_mgt.h"
#include "memory_mgt.h"
#include "timer_mgt.h"
#include "routing_info_mgt.h"
#include "cor_cpl_interface.h"
#include "cor_global_data.h"
#include "fields_gather_scatter_mgt.h"
#include "decomp_grid_mgt.h"
#include "common_utils.h"
#include "execution_report.h"
#include "performance_timing_mgt.h"
#include "ensemble_mgt.h"
#include "object_type_prefix.h"
#include "CCPL_api_mgt.h"
#include "original_grid_mgt.h"
#include "annotation_mgt.h"
#include "inout_interface_mgt.h"
#include "coupling_generator.h"
#include "runtime_trans_algorithm.h"
#include "IO_field_mgt.h"
#include "remapping_configuration_mgt.h"
#include "runtime_remapping_weights_mgt.h"


extern char software_name[];


extern Comp_comm_group_mgt_mgr *comp_comm_group_mgt_mgr;
extern Original_grid_mgt *original_grid_mgr;
extern Routing_info_mgt *routing_info_mgr;
extern Timer_mgt *timer_mgr;
extern Time_mgt *restart_read_timer_mgr;
extern Decomp_info_mgt *decomps_info_mgr;
extern Field_info_mgt *fields_info;
extern Memory_mgt *memory_manager;
extern Remap_mgt *grid_remap_mgr;
extern Fields_gather_scatter_mgt *fields_gather_scatter_mgr;
extern Decomp_grid_mgt *decomp_grids_mgr;
extern Performance_timing_mgt *performance_timing_mgr;
extern Ensemble_mgt *ensemble_mgr;
extern Annotation_mgt *annotation_mgr;
extern Components_time_mgt *components_time_mgrs;
extern Inout_interface_mgt *inout_interface_mgr;
extern IO_field_mgt *IO_fields_mgr;
extern Components_IO_output_procedures_mgt *components_IO_output_procedures_mgr;
extern Remapping_configuration_mgt *remapping_configuration_mgr;
extern Coupling_generator *coupling_generator;
extern Runtime_remapping_weights_mgt *runtime_remapping_weights_mgr;
extern H2D_remapping_wgt_file_container *all_H2D_remapping_wgt_files_info;


#endif
