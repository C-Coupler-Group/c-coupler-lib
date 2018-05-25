/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_FIELD_DATA_CLASS
#define REMAP_FIELD_DATA_CLASS


#include "common_utils.h"
#include "remap_statement_operand.h"


class Remap_grid_class;
class Runtime_remap_function;
class Remap_field_data_mgt;
class Remap_operator_grid;
class IO_basis;


class Remap_grid_data_class
{
    private:
        friend class Remap_grid_class;
        friend class Runtime_remap_function;
        friend class Remap_field_data_mgt;
        friend class Remap_operator_grid;
        
        Remap_grid_class *coord_value_grid;
        std::vector<Remap_grid_class *> sized_grids;
        Remap_data_field *grid_data_field;

        void generate_data_field_info(const char*, Remap_grid_class*);
        void generate_grid_info(Remap_grid_class*);

    public:
        Remap_grid_data_class(Remap_grid_class*, Remap_data_field*);
        Remap_grid_data_class(const char*, Remap_grid_class*, IO_basis*, const char*);
        Remap_grid_data_class(const char*, const char*, const char*, const char*, long, const char*, const char*);
        ~Remap_grid_data_class();
        bool have_data_content() const { return grid_data_field->read_data_size > 0; }
        Remap_grid_class *get_coord_value_grid() { return coord_value_grid; }
        bool is_unit_degree();
        Remap_data_field *get_grid_data_field() const { return grid_data_field; }
        Remap_grid_data_class *duplicate_grid_data_field(Remap_grid_class*, int, bool, bool);
        void interchange_grid_data(Remap_grid_class*);
        void reset_sized_grids(int, Remap_grid_class**);
        bool match_remap_grid_data(const char*);
        void transfer_field_attributes_to_another(Remap_grid_data_class*);
        void generate_analytic_values(const char*);
        void evaluate_error(Remap_grid_data_class*, Remap_grid_data_class*);
        void change_datatype_in_application(const char*);
        void write_grid_data_into_array(char **, long &, long &);
        Remap_grid_data_class(Remap_grid_class *, const char *, long&);
        void set_masked_cell_to_missing_value();
};


#endif
