/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef COMM_MGT_H
#define COMM_MGT_H


#include <mpi.h>
#include "common_utils.h"
#include "io_netcdf.h"
#include "tinyxml.h"
#include "restart_mgt.h"
#include "performance_timing_mgt.h"
#include <vector>


#define COMP_TYPE_PSEUDO_COUPLED   "pseudo_coupled_system"
#define COMP_TYPE_ACTIVE_COUPLED   "active_coupled_system"
#define COMP_TYPE_CPL              "cpl"
#define COMP_TYPE_ATM              "atm"
#define COMP_TYPE_GLC              "glc"
#define COMP_TYPE_ATM_CHEM         "atm_chem"
#define COMP_TYPE_OCN              "ocn"
#define COMP_TYPE_LND              "lnd"
#define COMP_TYPE_SEA_ICE          "sea_ice"
#define COMP_TYPE_WAVE             "wave"
#define COMP_TYPE_RUNOFF           "roff"
#define COMP_TYPE_ROOT             "ROOT"
#define DATAMODEL_NAME_PREFIX      "DATA_MODEL"
#define DATAINST_NAME_PREFIX       "DATA_INST"
#define ALGMODEL_NAME_PREFIX       "ALGORITHM_MODEL"


class Comp_comm_group_mgt_node
{
    private:
        int comp_id;
        char comp_name[NAME_STR_SIZE];                // The name of component    
        char full_name[NAME_STR_SIZE];
        char comp_type[NAME_STR_SIZE];
        char annotation_start[NAME_STR_SIZE];
        char annotation_end[NAME_STR_SIZE];
        char comp_ccpl_log_file_name[NAME_STR_SIZE];
        char comp_model_log_file_name[NAME_STR_SIZE];
        int comp_model_log_file_device;
        char working_dir[NAME_STR_SIZE];
        Comp_comm_group_mgt_node *parent;
        std::vector<Comp_comm_group_mgt_node*> children;
        MPI_Comm comm_group;
        std::vector<int> local_processes_global_ids;
        long *proc_latest_model_time;
        int current_proc_local_id;
        int current_proc_global_id;
        char *temp_array_buffer;
        long buffer_content_size;
        long buffer_content_iter;
        long buffer_max_size;
        bool definition_finalized;
        bool enabled_in_parent_coupling_generation;
        Restart_mgt *restart_mgr;
        int min_remote_lag_seconds;
        int max_remote_lag_seconds;
        char *log_buffer;
        int log_buffer_content_size;
        Performance_timing_mgt *performance_timing_mgr;

    public:
        Comp_comm_group_mgt_node(const char*, const char*, int, Comp_comm_group_mgt_node*, MPI_Comm&, bool, const char*);
        Comp_comm_group_mgt_node(TiXmlElement *, const char *, const char *);
        ~Comp_comm_group_mgt_node();
        MPI_Comm get_comm_group() const { return comm_group; }
        int get_comp_id() const { return comp_id; }
        int get_current_proc_local_id() const { return current_proc_local_id; }
        void transform_node_into_array();
        int get_buffer_content_size() { return buffer_content_size; }
        int get_buffer_content_iter() { return buffer_content_iter; }
        const char *get_comp_name() const { return comp_name; }
        const char *get_comp_type() const { return comp_type; }
        const char *get_comp_full_name() const { return full_name; }
        int get_num_children() { return children.size(); }
        int get_local_node_id() { return comp_id; }
        Comp_comm_group_mgt_node *get_child(int indx) { return children[indx]; }
        void reset_local_node_id(int);
        void reset_current_proc_local_id(int new_current_proc_local_id) { current_proc_local_id = new_current_proc_local_id; }
        bool is_definition_finalized() { return definition_finalized; }
        const char *get_annotation_start() { return annotation_start; }
        const char *get_annotation_end() { return annotation_end; }
        Comp_comm_group_mgt_node *get_parent() const { return parent; }
        void write_node_into_XML(TiXmlElement *);
        const char *get_working_dir() const { return working_dir; }
        void update_child(const Comp_comm_group_mgt_node*, Comp_comm_group_mgt_node*);
        void transfer_data_buffer(Comp_comm_group_mgt_node*);
        int get_num_procs() const { return local_processes_global_ids.size(); }
        int get_root_proc_global_id() const { return local_processes_global_ids[0]; }
        int get_local_proc_global_id(int);
        void confirm_coupling_configuration_active(int, bool, const char*);
        const char *get_full_name() { return full_name; }
        int open_comp_model_log_file(int*);
        const char *get_comp_ccpl_log_file_name() { return comp_ccpl_log_file_name; } 
        const char *get_comp_model_log_file(int &device_id) { device_id = comp_model_log_file_device; return comp_model_log_file_name; }
        bool is_real_component_model();
        Restart_mgt *get_restart_mgr() { return restart_mgr; }
        bool have_local_process(int);
        void allocate_proc_latest_model_time();
        void set_current_proc_current_time(int, int);
        void set_proc_latest_model_time(int, long);
        long get_proc_latest_model_time(int);
        void get_all_descendant_real_comp_fullnames(int, std::vector<const char*>&, char **, long &, long &);
        int get_min_remote_lag_seconds();
        int get_max_remote_lag_seconds();
        void update_min_max_remote_lag_seconds(int);
        void output_log(const char *, bool);
        Performance_timing_mgt *get_performance_timing_mgr() { return performance_timing_mgr; }        
};


class Comp_comm_group_mgt_mgr
{
    private:
        std::vector<Comp_comm_group_mgt_node*> global_node_array;
        Comp_comm_group_mgt_node *global_node_root;
        std::vector<const char*> root_comps_full_names;
        std::vector<bool> root_comps_enabled_in_parent_coupling_generation;
        bool definition_finalized;
        int current_proc_global_id;
        int num_total_global_procs;
        char executable_name[NAME_STR_SIZE];
        char root_working_dir[NAME_STR_SIZE];
        char internal_H2D_grids_dir[NAME_STR_SIZE];
		char internal_remapping_weights_dir[NAME_STR_SIZE];
        char components_processes_dir[NAME_STR_SIZE];
        char components_exports_dir[NAME_STR_SIZE];
        char active_coupling_connections_dir[NAME_STR_SIZE];
        char comps_ending_config_status_dir[NAME_STR_SIZE];
        char restart_common_dir[NAME_STR_SIZE];
        char runtime_config_root_dir[NAME_STR_SIZE];
        char root_comp_config_dir[NAME_STR_SIZE];
        char CCPL_platform_log_dir[NAME_STR_SIZE];
        char exe_log_file_name[NAME_STR_SIZE];
        int unique_comp_id_indx;
        char *log_buffer;
        int log_buffer_content_size;

    public:
        Comp_comm_group_mgt_mgr(const char*);
        ~Comp_comm_group_mgt_mgr();
        int register_component(const char*, const char*, MPI_Comm&, int, bool, int, const char*);
        bool is_legal_local_comp_id(int, bool);
        void transform_global_node_tree_into_array(Comp_comm_group_mgt_node*, Comp_comm_group_mgt_node**, int&);
        Comp_comm_group_mgt_node *get_global_node_of_local_comp(int, bool, const char*);
        MPI_Comm get_comm_group_of_local_comp(int, const char*);
        const char *get_executable_name() { return executable_name; }
        const char *get_annotation_start() { return global_node_array[0]->get_annotation_start(); }
        Comp_comm_group_mgt_node *get_root_component_model() { return global_node_array[1]; }
        const char *get_comp_ccpl_log_file_name(int);
        const char *get_comp_model_log_file(int, int &);
        void get_output_data_file_header(int, char*);
        Comp_comm_group_mgt_node *search_global_node(int);
        Comp_comm_group_mgt_node *search_global_node(const char*);
        Comp_comm_group_mgt_node *search_comp_with_comp_name(const char*);
        int get_current_proc_id_in_comp(int, const char *);
        int get_num_proc_in_comp(int, const char *);
        int get_current_proc_global_id() { return current_proc_global_id; }
        int get_num_total_global_procs() { return num_total_global_procs; }
        const char *get_root_working_dir() { return root_working_dir; }
        const char *get_internal_H2D_grids_dir() { return internal_H2D_grids_dir; }
		const char *get_internal_remapping_weights_dir() { return internal_remapping_weights_dir; }
        const char *get_components_processes_dir() { return components_processes_dir; }
        const char *get_components_exports_dir() { return components_exports_dir; }
        const char *get_active_coupling_connections_dir() { return active_coupling_connections_dir; }
        const char *get_comps_ending_config_status_dir() { return comps_ending_config_status_dir; }
        const char *get_restart_common_dir() { return restart_common_dir; }
        const char *get_config_root_dir() { return runtime_config_root_dir; }    
        const char *get_root_comp_config_dir() const { return root_comp_config_dir; }
        void get_root_comps_for_overall_coupling_generation(std::vector<const char *> &);
        void confirm_coupling_configuration_active(int, int, bool, const char*);
        const int *get_all_components_ids();
        Comp_comm_group_mgt_node *get_global_node_root() { return global_node_root; }
        bool has_comp_ended_configuration(const char*);
        void push_comp_node(Comp_comm_group_mgt_node *);
        Comp_comm_group_mgt_node *pop_comp_node();
        void check_validation();
        void set_current_proc_current_time(int, int, int);
        Comp_comm_group_mgt_node *load_comp_info_from_XML(int, const char *, MPI_Comm);
        const char *get_CCPL_platform_log_dir() { return CCPL_platform_log_dir; }
        bool is_comp_type_coupled(int, const char *, const char *);
        void output_log(const char *, bool);
        const char *get_exe_log_file_name() { return exe_log_file_name; }        
        void output_performance_timing();
        bool does_comp_name_include_reserved_prefix(const char *);
};


#endif
