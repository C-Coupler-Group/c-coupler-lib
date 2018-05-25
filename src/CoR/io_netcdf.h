/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/



#ifndef IO_NETCDF
#define IO_NETCDF


#include "io_basis.h"
#include <netcdf.h>
#include "remap_weight_of_strategy_class.h"


#define SCRIP_CENTER_LON_LABEL          "grid_center_lon"
#define SCRIP_CENTER_LAT_LABEL          "grid_center_lat"
#define SCRIP_VERTEX_LON_LABEL          "grid_corner_lon"
#define SCRIP_VERTEX_LAT_LABEL          "grid_corner_lat"
#define SCRIP_MASK_LABEL                "grid_imask"


class IO_netcdf: public IO_basis
{
    private:
        int ncfile_id;
        int rcode;
        bool io_with_time_info;
        int time_dim_id;
        int time_count;
        bool is_external_file;
        
        void write_field_data(Remap_grid_data_class*, Remap_grid_class*, bool, const char*, int, bool, bool);
        void datatype_from_netcdf_to_application(nc_type, char*, const char*);
        void datatype_from_application_to_netcdf(const char*, nc_type*);
        void report_nc_error();
        bool get_file_field_attribute(const char *, const char *, char *, char *);

    public:
        IO_netcdf(int);
        IO_netcdf(const char*, const char*, const char*, bool);
        ~IO_netcdf();
        bool read_data(Remap_data_field*, int, bool);
        void write_grided_data(Remap_grid_data_class*, bool, int, int, bool);
        void write_remap_weights(Remap_weight_of_strategy_class*);
        long get_dimension_size(const char*, MPI_Comm, bool);
        void read_remap_weights(Remap_weight_of_strategy_class*, Remap_strategy_class*, bool);
        void put_global_attr(const char*, const void*, const char *, const char *, int);
        void read_file_field(const char*, void**, int*, char*, MPI_Comm, bool);
        bool get_file_field_string_attribute(const char*, const char *, char*, char *, MPI_Comm, bool);
        void write_grid(Remap_grid_class*, bool, bool);
};


#endif 
