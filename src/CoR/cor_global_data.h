/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef COR_GLOBAL_DATA
#define COR_GLOBAL_DATA


#include "io_netcdf.h"
#include "io_mgt.h"
#include "remap_strategy_mgt.h"
#include "remap_grid_mgt.h"
#include "remap_operator_mgt.h"
#include "remap_grid_data_mgt.h"
#include "remap_common_utils.h"
#include "remap_weight_of_strategy_mgt.h"
#include "runtime_remap_function.h"
#include "syntax_check.h"
#include "remap_mgt.h"
#include "execution_report.h"
#include <stdio.h>


#define DATA_TYPE_DOUBLE         "real8"
#define DATA_TYPE_FLOAT          "real4"
#define DATA_TYPE_BOOL           "logical"
#define DATA_TYPE_CHAR           "char"
#define DATA_TYPE_LONG           "long"
#define DATA_TYPE_INT            "integer"
#define DATA_TYPE_SHORT          "short"
#define DATA_TYPE_STRING         "string"

#define FILL_VALUE_LABEL                      "_FillValue"
#define MISS_VALUE_LABEL                      "missing_value"
#define DEFAULT_FILL_VALUE                    (1.0e20)


extern IO_mgt *io_manager;
extern Remap_strategy_mgt *remap_strategy_manager;
extern Remap_grid_mgt *remap_grid_manager;
extern Remap_operator_mgt *remap_operator_manager;
extern Remap_field_data_mgt *remap_field_data_manager;
extern Remap_weight_of_strategy_mgt *remap_weights_of_strategy_manager;
extern Remap_weight_of_operator_mgt *sequential_remap_weight_of_operator_manager;
extern Remap_weight_of_operator_mgt *parallel_remap_weight_of_operator_manager;
extern Runtime_remap_function *current_runtime_remap_function;
extern Remap_operator_grid *current_runtime_remap_operator_grid_src;
extern Remap_operator_grid *current_runtime_remap_operator_grid_dst;
extern Remap_operator_basis *current_runtime_remap_operator;
extern int line_number;
extern int execution_phase_number;
extern bool is_coord_unit_degree[];
extern bool is_master_process_in_computing_node;


extern int get_data_type_size(const char*);
extern void check_application_io_datatype_consistency(const char*, const char*, const char*);

#endif
