/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef SYNTAX_CHECK_H
#define SYNTAX_CHECK_H


#include "remap_statement_operand.h"

extern void check_is_parameter_string_type(const char*, int, Remap_statement_operand*, const char*);
extern long get_size_value_from_parameter(const char*, int, Remap_statement_operand*, const char*);
extern long get_int_value_from_parameter(const char*, int, Remap_statement_operand*, const char*);
extern double get_float_value_from_parameter(const char*, int, Remap_statement_operand*, const char*);
extern void check_is_parameter_object_type_grid(const char*, int, Remap_statement_operand*, const char*);
extern void check_is_parameter_object_type_IO(const char*, int, Remap_statement_operand*, const char*);
extern void check_is_parameter_object_type_field_data(const char*, int, Remap_statement_operand*, const char*);
extern void check_is_parameter_object_type_remap_operator(const char*, int, Remap_statement_operand*, const char*);
extern void check_is_parameter_object_type_remap_scheme(const char*, int, Remap_statement_operand*, const char*);
extern void check_is_parameter_object_type_remap_weights(const char*, int, Remap_statement_operand*, const char*);
extern void check_is_parameter_grid_mask_field(const char*, Remap_statement_operand*, const char*);
extern void check_is_parameter_grid_center_field(const char*, Remap_statement_operand*, const char*);
extern void check_is_parameter_grid_boundary(const char*, Remap_statement_operand*, const char*);
extern void check_is_parameter_grid_field(const char*, Remap_statement_operand*, const char*);


#endif

