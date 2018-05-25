/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef RESTART_MGT_H
#define RESTART_MGT_H


#define RESTART_BUF_TYPE_TIME            "time_restart"
#define RESTART_BUF_TYPE_INTERFACE       "interface"


#include "common_utils.h"
#include "timer_mgt.h"
#include "io_netcdf.h"
#include "memory_mgt.h"
#include <vector>


class Restart_mgt;
class Comp_comm_group_mgt_node;


class Restart_buffer_container
{
    private:
        char comp_full_name[NAME_STR_SIZE];
        char buf_type[NAME_STR_SIZE];
        char keyword[NAME_STR_SIZE];
        char *buffer_content;
        long buffer_content_iter;
        long buffer_content_size;
        long buffer_max_size;
        Restart_mgt *restart_mgr;

    public:
        Restart_buffer_container(const char *, const char *, const char *, Restart_mgt *);
        Restart_buffer_container(const char *, long &, const char *, Restart_mgt *);
        ~Restart_buffer_container() { delete [] buffer_content; }        
        const char *get_buffer_content() { return buffer_content; }
        long get_buffer_content_iter() { return buffer_content_iter; }
        void dump_out(char **, long &, long &);
        bool match(const char *, const char *);
        void dump_in_string(const char *, long);
        void dump_in_data(const void *, long);    
        char **get_buffer_content_ptr() { return &buffer_content; }
        long *get_buffer_content_iter_ptr() { return &buffer_content_iter; }
        long *get_buffer_max_size_ptr() { return &buffer_max_size; }
        void load_restart_data(void *, long);        
        char *load_restart_string(char *, long &, long);
        const char *get_input_restart_mgt_info_file();
        Restart_mgt *get_restart_mgr() { return restart_mgr; }
};


class Restart_mgt
{
    private:
        long last_restart_write_full_time;
        long last_restart_write_elapsed_time;
        std::vector<Restart_buffer_container*> restart_write_buffer_containers;
        std::vector<Restart_buffer_container*> restart_read_buffer_containers;
        std::vector<std::pair<Field_mem_info*, bool> > restarted_field_instances;
        Comp_comm_group_mgt_node *comp_node;
        Time_mgt *time_mgr;
        char *input_restart_mgt_info_file;
        char *restart_read_annotation;
        bool restart_mgt_info_written;
        IO_netcdf *restart_write_data_file;
        IO_netcdf *backup_restart_write_data_file;
        char *restart_read_data_file_name;
        bool restart_normal_fields_enabled;
        bool are_all_restarted_fields_read;
        bool bypass_import_fields_at_read;
        bool bypass_import_fields_at_write;

    public:
        Restart_mgt(Comp_comm_group_mgt_node*);
        ~Restart_mgt();
        void clean(bool);
        void do_restart_write(const char *, bool, bool);
        void write_restart_mgt_into_file();
        void read_restart_mgt_info(bool, const char *, const char *);
        void read_restart_mgt_info(const char *, const char *);
        Restart_buffer_container *search_restart_buffer(const char *, const char*);
        int get_comp_id();
        void get_file_name_in_rpointer_file(char *);
        const char *get_input_restart_mgt_info_file();
        const char *get_restart_read_annotation();
        Restart_buffer_container *apply_restart_buffer(const char *, const char *, const char *);
        bool is_in_restart_write_window(long, bool);
        bool is_in_restart_read_window(long);
        void write_restart_field_data(Field_mem_info *, const char*, const char*, bool);
        void read_restart_field_data(Field_mem_info *, const char *, const char *, bool, const char *, bool, const char*);
        const char *get_restart_read_data_file_name() { return restart_read_data_file_name; }
        void add_restarted_field_instance(Field_mem_info*, bool);
        void get_field_IO_name(char *, Field_mem_info*, const char *, const char*, bool);
        void read_all_restarted_fields(const char*);
        bool check_restart_read_started();
        bool get_are_all_restarted_fields_read() { return are_all_restarted_fields_read; }
        bool get_bypass_import_fields_at_read() { return bypass_import_fields_at_read; }
};


#endif
