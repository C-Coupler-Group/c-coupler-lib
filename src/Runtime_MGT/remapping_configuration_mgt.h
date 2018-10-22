/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAPPING_CONFIGURATION_MGT
#define REMAPPING_CONFIGURATION_MGT


#define REMAP_ALGORITHM_TYPE_H2D     1
#define REMAP_ALGORITHM_TYPE_V1D     2
#define REMAP_ALGORITHM_TYPE_T1D     3


#include <mpi.h>
#include "io_netcdf.h"
#include <vector>
#include "tinyxml.h"
#include "common_utils.h"
#include "original_grid_mgt.h"


class H2D_remapping_wgt_file_info
{
    private:
        const char *wgt_file_name;
        int src_grid_size;
        int dst_grid_size;
        double *src_center_lon;
        double *src_center_lat;
        double *src_area;
        long checksum_src_mask;
        double *dst_center_lon;
        double *dst_center_lat;
        double *dst_area;
        long checksum_dst_mask;
        long num_wgts;
        long *wgts_src_indexes;
        long *wgts_dst_indexes;
        double *wgts_values;
        std::vector<std::pair<Original_grid_info*, Original_grid_info*> > matched_grid_pair;

    public:
        H2D_remapping_wgt_file_info(const char*);
        H2D_remapping_wgt_file_info(const char*, long*);
        ~H2D_remapping_wgt_file_info();
        void clean();
        long get_grid_field_checksum_value(const char *, IO_netcdf *, int, MPI_Comm, bool);
        bool match_H2D_remapping_wgt(Original_grid_info*, Original_grid_info*);
        const char *get_wgt_file_name() { return wgt_file_name; }
        void write_remapping_wgt_file_info_into_array(char **, long &, long &);
        void read_remapping_weights(int);
        long get_num_wgts() { return num_wgts; }
        long *get_wgts_src_indexes () { return wgts_src_indexes; }
        long *get_wgts_dst_indexes () { return wgts_dst_indexes; }
        double *get_wgts_values () { return wgts_values; }        
        void get_checksum_mask(int, const char *, int, long &);
        void read_grid_size(int, const char *, int &);
        void read_weight_grid_data(int, const char *, const char *, void **, int, bool);
        double *get_src_area() { return src_area; }
        double *get_dst_area() { return dst_area; }
};


class H2D_remapping_wgt_file_mgt
{
    private:
        std::vector<H2D_remapping_wgt_file_info*> H2D_remapping_wgt_files;

    public:
        H2D_remapping_wgt_file_mgt() { }
        H2D_remapping_wgt_file_mgt(TiXmlElement*, const char*);
        H2D_remapping_wgt_file_mgt(const char *, long *);
        ~H2D_remapping_wgt_file_mgt() {}        
        void append_remapping_weights(H2D_remapping_wgt_file_mgt *);
        void print();
        void clean() { H2D_remapping_wgt_files.clear(); }
        H2D_remapping_wgt_file_info *search_H2D_remapping_weight(Original_grid_info*, Original_grid_info*, const char *, int);
		H2D_remapping_wgt_file_info *search_H2D_remapping_weight(const char *);
        void write_remapping_wgt_files_info_into_array(char **, long &, long &);
        bool is_the_same_as_another(H2D_remapping_wgt_file_mgt*);
        H2D_remapping_wgt_file_mgt *clone();
        void shrink(Original_grid_info*, Original_grid_info*);
        bool is_empty() { return H2D_remapping_wgt_files.size() == 0; } 
};


class H2D_remapping_wgt_file_container
{
    private:
        std::vector<H2D_remapping_wgt_file_info*> H2D_remapping_wgt_files;

    public:
        H2D_remapping_wgt_file_container() {}
        ~H2D_remapping_wgt_file_container();
        H2D_remapping_wgt_file_info *search_wgt_file_info(const char*);
        void add_wgt_file_info(H2D_remapping_wgt_file_info *);
};



class Remapping_algorithm_specification
{
    private:
        int comp_id;        
        int type_id;
        const char *algorithm_name;
        std::vector<char *> parameters_name;
        std::vector<char *> parameters_value;

    public:
        Remapping_algorithm_specification(const Remapping_algorithm_specification*);
        Remapping_algorithm_specification(const char*, int);
        Remapping_algorithm_specification(const char *, long*);
        Remapping_algorithm_specification(int, TiXmlElement*, const char*, int);
        Remapping_algorithm_specification() { algorithm_name = NULL; }
        ~Remapping_algorithm_specification();
        void print();
        void clean();
        void write_remapping_algorithm_specification_into_array(char **, long &, long &);
        const char *get_algorithm_name() { return algorithm_name; }
        int get_num_parameters() { return parameters_name.size(); }
        void get_parameter(int, char *, char *);
        Remapping_algorithm_specification *clone();
        bool is_the_same_as_another(Remapping_algorithm_specification*);
};


class Remapping_setting
{
    private:
        int comp_id;
        int XML_start_line_number;
        Remapping_algorithm_specification *H2D_remapping_algorithm;
        Remapping_algorithm_specification *V1D_remapping_algorithm;
        Remapping_algorithm_specification *T1D_remapping_algorithm;
        H2D_remapping_wgt_file_mgt *H2D_remapping_wgt_file_mgr;
        int field_specification_manner;    // 0 means default of a component; 1 means type; 2 means name
        std::vector<const char *> fields_specification;

    public:
        Remapping_setting();
        Remapping_setting(const char*, const char*);
        Remapping_setting(int, TiXmlElement*, const char*);
        ~Remapping_setting();
        void detect_conflict(Remapping_setting*, const char*);
        void reset_remapping_setting();
        int get_field_specification_manner() { return field_specification_manner; }
        void get_field_remapping_setting(Remapping_setting&, const char*);
        Remapping_algorithm_specification *get_H2D_remapping_algorithm() { return H2D_remapping_algorithm; }
        Remapping_algorithm_specification *get_V1D_remapping_algorithm() { return V1D_remapping_algorithm; }
        Remapping_algorithm_specification *get_T1D_remapping_algorithm() { return T1D_remapping_algorithm; }
        void write_remapping_setting_into_array(char **, long &, long &);
        void read_remapping_setting_from_array(const char *, long &);
        void print();
        Remapping_setting *clone();
        bool is_the_same_as_another(Remapping_setting*);
        void append_H2D_remapping_weights(Remapping_setting *);
        H2D_remapping_wgt_file_info *search_H2D_remapping_weight(Original_grid_info *, Original_grid_info*, int);
        long calculate_checksum();
        void shrink(Original_grid_info*, Original_grid_info*);
};


class Remapping_configuration
{
    private:
        int comp_id;
        std::vector<Remapping_setting*> remapping_settings;

    public:
        Remapping_configuration();
        Remapping_configuration(int, const char*, TiXmlDocument *XML_file);
        ~Remapping_configuration();
        int get_comp_id() { return comp_id; }
        bool get_field_remapping_setting(Remapping_setting&, const char*);
};


class Remapping_configuration_mgt
{
    private:
        std::vector<Remapping_configuration*> remapping_configurations;

    public:
        Remapping_configuration_mgt() {}
        ~Remapping_configuration_mgt();
        void add_remapping_configuration(int);
        Remapping_configuration *search_remapping_configuration(int);
        void get_field_remapping_setting(Remapping_setting &, int, const char*);
};


#endif

