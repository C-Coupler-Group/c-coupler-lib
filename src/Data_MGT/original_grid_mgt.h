/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef ORIGINAL_GRID_MGT
#define ORIGINAL_GRID_MGT


#include "remap_grid_class.h"
#include "remap_mgt.h"
#include <vector>


#define BOTTOM_FIELD_VARIATION_UNSET             ((int)-1)
#define BOTTOM_FIELD_VARIATION_STATIC            ((int) 0)
#define BOTTOM_FIELD_VARIATION_DYNAMIC           ((int) 1)
#define BOTTOM_FIELD_VARIATION_EXTERNAL          ((int) 2)


class Original_grid_info
{
    private: 
        int comp_id;
        int grid_id;
        const char *grid_name;
        const char *bottom_field_name;
        const char *comp_full_name;
        Remap_grid_class *original_CoR_grid;
        Remap_grid_class *H2D_sub_CoR_grid;
        Remap_grid_class *V1D_sub_CoR_grid;
        Remap_grid_class *T1D_sub_CoR_grid;
        int H2D_sub_grid_order;
        int V1D_sub_grid_order;
        int T1D_sub_grid_order;
		int V3D_lev_field_id;
		int V3D_lev_field_variation_type;   // 0: static; 1: dynamic
        int bottom_field_id;
        int bottom_field_variation_type;   // 0: static; 1: dynamic; 2: external
        Original_grid_info *mid_point_grid;
        Original_grid_info *interface_level_grid;
        long checksum_H2D_mask;

        void generate_remapping_grids();
        
    public:
        Original_grid_info(int, int, const char*, const char*, Remap_grid_class*, bool);
        const char *get_annotation();
        const char *get_grid_name() const { return grid_name; }
        int get_local_grid_id() const { return grid_id; }
        ~Original_grid_info();
        int get_comp_id() const { return comp_id; }
        int get_grid_id() const { return grid_id; }
        int get_bottom_field_variation_type() const { return bottom_field_variation_type; }
		int get_V3D_lev_field_variation_type() const { return V3D_lev_field_variation_type; }
        const char *get_bottom_field_name() const { return bottom_field_name; }
        void set_bottom_field_variation_type(int type) { bottom_field_variation_type = type; }
		void set_V3D_lev_field_variation_type(int type) { V3D_lev_field_variation_type = type; }
        void set_grid_checksum(long);
        void set_unique_bottom_field(int, int, const char*);
		void set_unique_3D_lev_field(int, const char *, const char *);
        Remap_grid_class *get_original_CoR_grid() const { return original_CoR_grid; }
        Remap_grid_class *get_H2D_sub_CoR_grid() { return H2D_sub_CoR_grid; }
        Remap_grid_class *get_V1D_sub_CoR_grid() { return V1D_sub_CoR_grid; }
        Remap_grid_class *get_T1D_sub_CoR_grid() { return T1D_sub_CoR_grid; }
        bool is_V1D_sub_grid_after_H2D_sub_grid();
        bool is_3D_grid() { return H2D_sub_CoR_grid != NULL && V1D_sub_CoR_grid != NULL && T1D_sub_CoR_grid == NULL; }
        bool is_H2D_grid() { return H2D_sub_CoR_grid != NULL && V1D_sub_CoR_grid == NULL && T1D_sub_CoR_grid == NULL; } 
        void write_grid_into_array(char **, long &, long &);
        int get_bottom_field_id() { return bottom_field_id; }        
		int get_V3D_lev_field_id() { return V3D_lev_field_id; }
        void get_grid_data(int, const char*, const char*, int, char*, const char*, const char*);
        Original_grid_info *get_interface_level_grid() { return interface_level_grid; }
        Original_grid_info *get_mid_point_grid() { return mid_point_grid; }
        void set_mid_point_grid(Original_grid_info*);
        long get_checksum_H2D_mask() { return checksum_H2D_mask; }
        bool is_H2D_grid_and_the_same_as_another_grid(Original_grid_info *);
        double *get_center_lon_values();
        double *get_center_lat_values();
        void reset_grid_data();
        const char *get_comp_full_name() { return comp_full_name; }
};


class Original_grid_mgt
{
    private:
        std::vector<Original_grid_info*> original_grids;
        char CoR_script_name[NAME_STR_SIZE];
        Remap_mgt *CoR_grids;

    public:
        Original_grid_mgt();
        ~Original_grid_mgt();
        void initialize_CoR_grids();
        int get_CoR_defined_grid(int, const char*, const char*, const char*);
        Original_grid_info *search_grid_info(const char*, int);
        Original_grid_info *search_grid_info(int);
        Remap_grid_class *get_original_CoR_grid(int) const;
        Original_grid_info *get_original_grid(int) const;
        bool is_grid_id_legal(int) const;        
        int get_comp_id_of_grid(int) const;
        const char *get_name_of_grid(int) const;
        int get_grid_size(int, const char*) const;
        int get_grid_id(int, const char*, const char*);
        int add_original_grid(int, const char*, Remap_grid_class*);
        int get_num_grid_levels(int);
        bool is_V1D_sub_grid_after_H2D_sub_grid(int);
        void common_checking_for_H2D_registration_via_data(int, const char *, const char *, const char *, char *, const char *, int, int, int, int, int, int *, char *, char *, char *, char *, char *, char *, char *, char *, const char *, int);
        int create_H2D_grid_from_global_data(int, const char *, const char *, const char *, const char *, int, int, int, int, int, int, int, int, int, char *, char *, char *, char *, char *, char *, int *, char *, char *, char *, const char *);
        int register_H2D_grid_via_global_data(int, const char *, const char *, const char *, char *, const char *, int, int, int, int, 
                                               int, int, int, int, char *, char *, char *, char *, char *, char *, int *, char *, char *, char *, const char *, int);
        int register_H2D_grid_via_local_data(int, const char *, const char *, const char *, char *, const char *, int, int, int, int, int, 
                                             int, int, int, int, int *, char *, char *, char *, char *, char *, char *, int *, char *, char *, char *, const char *, int *, const char *, int);
        int register_H2D_grid_via_file(int, const char *, const char *, const char *);
        int register_H2D_grid_via_comp(int, const char *, const char *);
        int register_V1D_grid_via_data(int, int, const char *, int, const char *, int, double, const double *, const double *, const char *);
        int register_md_grid_via_multi_grids(int, const char*, int, int, int, int, int *, const char*);
        void set_3d_grid_bottom_field(int, int, int, int, int, const char*, const char*);		
		void set_3D_grid_3D_vertical_coord_field_inst(int, int, const char *, const char *);
        void register_mid_point_grid(int, int*, int*, int, const int*, const char*, const char *);
        void delete_external_original_grids();
        void calculate_min_max_H2D_coord_value(int, char *, char *, int, int, const char *, double &, double &);
		Original_grid_info *search_or_register_internal_grid(int, Remap_grid_class *);
};




#endif

