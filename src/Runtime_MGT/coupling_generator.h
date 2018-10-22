/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu and then
  *  modified by Dr. Cheng Zhang and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn or
  *  Dr. Cheng Zhang via zhangc-cess@tsinghua.edu.cn
  ***************************************************************/



#ifndef COUPLING_GENERATOR_H
#define COUPLING_GENERATOR_H


#include <vector>
#include <map>
#include "inout_interface_mgt.h"
#include "dictionary.h"
#include "remapping_configuration_mgt.h"
#include "runtime_remapping_weights_mgt.h"


#define USING_INSTANTANEOUS_VALUE            0
#define USING_AVERAGE_VALUE                  1


class Import_interface_configuration;
class Coupling_generator;
class IO_output_procedure;


struct Interface_field_info
{
    char grid_name[NAME_STR_SIZE];
    char data_type[NAME_STR_SIZE];
    char decomp_name[NAME_STR_SIZE];
    char unit[NAME_STR_SIZE];
    Runtime_remapping_weights *runtime_remapping_weights;
    int bottom_field_indx;
};


struct V3D_grid_bottom_field_coupling_info
{
    Field_mem_info *bottom_field_inst;
    Runtime_remapping_weights *V3D_runtime_remapping_weights;
    Runtime_remapping_weights *H2D_runtime_remapping_weights;
    int field_connection_indx;
    bool is_dynamic_bottom_field;
};


class Coupling_connection
{
    private:
        friend class Coupling_generator;
        friend class Connection_coupling_procedure;
        friend class IO_output_procedure;
        friend class Inout_interface;
        friend class Inout_interface_mgt;
        int connection_id;
        std::vector<const char*> fields_name;
        std::vector<std::pair<const char*, const char*> > src_comp_interfaces;
        char dst_comp_full_name[NAME_STR_SIZE];
        char dst_interface_name[NAME_STR_SIZE];
        std::vector<Interface_field_info*> src_fields_info;
        std::vector<Interface_field_info*> dst_fields_info;
        std::vector<V3D_grid_bottom_field_coupling_info*> src_bottom_fields_coupling_info; 
        std::vector<V3D_grid_bottom_field_coupling_info*> dst_bottom_fields_coupling_info;
        Inout_interface *import_interface;
        Inout_interface *export_interface;
        Coupling_timer *src_timer;
        Coupling_timer *dst_timer;
        int src_time_step_in_second;
        int dst_time_step_in_second;
        int src_inst_or_aver;
        int dst_inst_or_aver;
        int src_current_year, src_current_month, src_current_day, src_current_second;
        int dst_current_year, dst_current_month, dst_current_day, dst_current_second;
        Connection_coupling_procedure *import_procedure;
        Connection_coupling_procedure *export_procedure;
        Comp_comm_group_mgt_node *src_comp_node;
        Comp_comm_group_mgt_node *dst_comp_node;
        int current_proc_id_src_comp;
        int    current_proc_id_dst_comp;
        int src_comp_root_proc_global_id;
        int dst_comp_root_proc_global_id;
        MPI_Comm union_comm;
        int * src_proc_ranks_in_union_comm;
        int * dst_proc_ranks_in_union_comm;

        void write_field_info_into_array(Field_mem_info *, char **, long &, long &);
        void write_connection_fields_info_into_array(Inout_interface *, char **, long &, long &, Coupling_timer**, int&, int&, int&, int&, int&, int&);    
        void read_fields_info_from_array(std::vector<Interface_field_info*> &, const char*, long);
        void read_connection_fields_info_from_array(std::vector<Interface_field_info*>&, const char *, long, int, Coupling_timer **, int&, int&, int&, int&, int&, int&);
        void exchange_connection_fields_info();
        void exchange_bottom_fields_info();
        void generate_interpolation(bool);
        bool exchange_grid(Comp_comm_group_mgt_node *, Comp_comm_group_mgt_node *, const char *);
        void exchange_remapping_setting(int, Remapping_setting &);
        void add_bottom_field_coupling_info(int, Runtime_remapping_weights*, Remapping_setting *);
        void generate_src_bottom_field_coupling_info();

    public:
        Coupling_connection(int);
        ~Coupling_connection();
        void generate_a_coupling_procedure(bool);
        void create_union_comm();
        void generate_data_transfer();
        Field_mem_info *get_bottom_field(bool, int);
        bool get_is_bottom_field_dynamic(int);
        const char *get_dst_comp_full_name() { return dst_comp_full_name; }
        const char *get_dst_interface_name() { return dst_interface_name; }
        const char *get_src_comp_full_name() { return src_comp_interfaces[0].first; }
};


class Import_direction_setting
{
    private:
        char interface_name[NAME_STR_SIZE];
        int status;                      // 0 means off; 1 means on
        int fields_default_setting;       // 0 means off; 1 means all; 2 means remain"
        int components_default_setting;  // 0 means off; 1 means all;
        std::vector<const char*> fields_name;
        std::vector<std::pair<const char*, const char*> > producers_info;
        
        // ... remapping and merge setting;
    
    public:
        Import_direction_setting(int, Import_interface_configuration *, const char*, const char*, TiXmlElement*, const char*, std::vector<const char*>&, int*, bool);
        ~Import_direction_setting();
};


class Import_interface_configuration
{
    private:
        char interface_name[NAME_STR_SIZE];
        std::vector<Import_direction_setting*> import_directions;
        std::vector<const char*> fields_name;
        std::vector<std::vector<std::pair<const char*, const char*> > > fields_src_producers_info;

    public:
        Import_interface_configuration(int, const char*, const char*, TiXmlElement*, const char*, Inout_interface_mgt *, bool);
        ~Import_interface_configuration();
        const char *get_interface_name() { return interface_name; }
        void add_field_src_component(int comp_id, const char*, std::pair<const char*, const char*>);
        void get_field_import_configuration(const char*, std::vector<std::pair<const char*, const char*> >&);
};


class Component_import_interfaces_configuration
{
    private:
        char comp_full_name[NAME_STR_SIZE];
        std::vector<Import_interface_configuration*> import_interfaces_configuration;

    public:    
        Component_import_interfaces_configuration(int, const char *, Inout_interface_mgt *, bool);
        ~Component_import_interfaces_configuration();
        void get_interface_field_import_configuration(const char*, const char*, std::vector<std::pair<const char*,const char*> >&);
};



class Coupling_generator
{
    private:
        std::map<int,std::vector<std::pair<const char*, const char*> > > export_fields_dst_components;
        std::vector<const char*> string_in_export_fields_dst_components;
        std::vector<const char *> all_comp_fullnames_for_coupling_generation;
        std::vector<int> individual_or_family_generation;  // must be 1 (individual) or 2 (family)
        Dictionary<int> *import_field_index_lookup_table;
        Dictionary<int> *export_field_index_lookup_table;
        std::vector<Coupling_connection*> all_coupling_connections;
        std::vector<Coupling_connection*> all_IO_connections;
        int latest_connection_id;
        
        void generate_interface_fields_source_dst(const char*, int);
        void generate_components_connections();        
        void generate_coupling_procedures_common(int, MPI_Comm, bool, bool, const char*);
        void sort_comp_full_names(std::vector<const char*> &, std::vector<int>*);

    public:
        Coupling_generator() { latest_connection_id = 1; import_field_index_lookup_table = NULL; export_field_index_lookup_table = NULL; }
        ~Coupling_generator();
        void clear();
        void generate_coupling_procedures_internal(int, bool, bool, const char*);
        void generate_IO_procedures();
        int apply_connection_id() {  return (++latest_connection_id); }
        int get_latest_connection_id() { return latest_connection_id; }
        void set_latest_connection_id(int connection_id) { latest_connection_id = connection_id; }
        void synchronize_latest_connection_id(MPI_Comm);
        void transfer_interfaces_info_from_one_component_to_another(std::vector<Inout_interface*> &, Comp_comm_group_mgt_node *, Comp_comm_group_mgt_node *);
        void begin_external_coupling_generation();
        void add_comp_for_external_coupling_generation(const char *, int, const char*);
        void do_external_coupling_generation(int, const char*);
        void load_comps_full_names_from_config_file(int, const char *, int, int, int *, const char *);
        void get_one_comp_full_name(int, const char*, int, int, char *, int *, const char *);
        void do_overall_coupling_generation(const char*, const char *);
};


#endif

