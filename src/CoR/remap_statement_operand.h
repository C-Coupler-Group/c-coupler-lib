/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_STATEMENT_OPERAND
#define REMAP_STATEMENT_OPERAND


#include "common_utils.h"
#include <vector>


class Remap_grid_class;


struct Remap_statement_object
{
    char object_name[256];
    char object_type[256];
    void *object_pointer;
};


struct Remap_statement_operand
{
    Remap_statement_object *object;
    int num_extension_names;
    char extension_names[16][256];
};


struct Remap_field_attribute
{
    char attribute_name[256];
    char attribute_type[256];
    char attribute_value[8192];
    unsigned long attribute_size;
};


class Remap_data_field
{
    public:
        void *data_buf;
        char field_name_in_application[NAME_STR_SIZE];
        char field_name_in_IO_file[NAME_STR_SIZE];
        char data_type_in_application[NAME_STR_SIZE];
        char data_type_in_IO_file[NAME_STR_SIZE];
        long required_data_size;
        long read_data_size;
        std::vector<Remap_field_attribute> field_attributes;
        bool have_fill_value;
        double fill_value;

        Remap_data_field();
        ~Remap_data_field();
        Remap_data_field *duplicate_remap_data_field(long, bool);
        void interchange_remap_data_field(Remap_data_field*, Remap_grid_class*, Remap_grid_class*);
        void push_back_attribute(Remap_field_attribute field_attribute) { field_attributes.push_back(field_attribute); }
        void read_fill_value();
        void set_fill_value(void*);
        void clean_fill_value();
        void set_field_unit(const char*);
        void set_field_long_name(const char*);
        void initialize_to_fill_value();
        void set_scale_factor_and_add_offset(double, double);
        void clean_scale_factor_and_add_offset_info();
        void read_scale_factor_and_add_offset(double*, double*);
};


#endif
