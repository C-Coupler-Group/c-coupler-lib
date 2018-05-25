/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef IO_FIELD_MGT
#define IO_FIELD_MGT


#include <vector>
#include "common_utils.h"
#include "timer_mgt.h"
#include "inout_interface_mgt.h"


class Coupling_connection;


class IO_field
{
    private:
        int IO_field_id;
        int field_instance_id;
        int comp_id;
        char field_IO_name[NAME_STR_SIZE];
        char field_long_name[NAME_STR_SIZE];
        char field_unit[NAME_STR_SIZE];

    public:
        IO_field(int, int, const char *, const char *);
        IO_field(int, int, int, int, void *, const char *, const char *, const char *, const char *, const char *);
        ~IO_field() {}
        int get_comp_id() { return comp_id; }
        int get_IO_field_id() { return IO_field_id; }
        int get_field_instance_id() { return field_instance_id; }
        const char *get_field_IO_name() { return field_IO_name; }
};


class IO_output_procedure
{
    private:
        int comp_id;
        int procedure_id;
        Inout_interface *export_interface;
        Inout_interface *import_interface;
        Coupling_timer *field_timer;
        Coupling_timer *file_timer;
        Time_mgt *time_mgr;
        std::vector<IO_field*> IO_fields;
        std::vector<Field_mem_info*> data_write_field_insts;
        IO_netcdf *netcdf_file_object;
        bool write_grid_name;
        int inst_or_aver;
        int *field_update_status;

    public:
        IO_output_procedure(int, int, Coupling_timer *, Coupling_timer *, bool);
        void include_all_component_io_fields();
        Coupling_connection *generate_coupling_connection(int);
        void execute();
        ~IO_output_procedure();
};


class IO_field_mgt
{
    private:
        friend class IO_output_procedure;
        std::vector<IO_field*> IO_fields;

    public:
        IO_field_mgt() {}
        ~IO_field_mgt();
        int register_IO_field(int, const char *, const char *);
        int register_IO_fields(int, int, int *, const char *);
        int register_IO_field(int, int, int, void *, const char *, const char *, const char *, const char *, const char *);
        IO_field *search_IO_field(int, const char*);
        void check_for_registering_IO_field(IO_field *, const char *, int);
};


class Component_IO_output_procedures
{
    private:
        std::vector<IO_output_procedure*> IO_output_procedures;
        int comp_id;

    public:
        Component_IO_output_procedures(int, const char*, bool);
        void generate_coupling_connection(std::vector<Coupling_connection*> &, int);
        void execute();
        ~Component_IO_output_procedures();
};


class Components_IO_output_procedures_mgt
{
    private:
        std::vector<Component_IO_output_procedures*> components_IO_output_procedures;

    public:
        Components_IO_output_procedures_mgt() {}
        ~Components_IO_output_procedures_mgt();
        void add_component_IO_output_procedures(int, const char*, bool);
        void add_all_components_IO_output_procedures();
        Component_IO_output_procedures *get_component_IO_output_procedures(int);
};

#endif

