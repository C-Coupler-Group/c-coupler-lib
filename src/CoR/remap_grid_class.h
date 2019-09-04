/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_GRID_CLASS
#define REMAP_GRID_CLASS


#include "common_utils.h"
#include "remap_statement_operand.h"
#include <vector>
#include <stdio.h>
#include <string.h>


#define GRID_BOUNDARY_LABEL                   "boundary"
#define GRID_CENTER_LABEL                     "center"
#define GRID_VERTEX_LABEL                     "vertex"
#define GRID_NUM_VERTEX_LABEL                 "num_vertex"
#define GRID_MASK_LABEL                       "mask"
#define GRID_AREA_LABEL                       "area"
#define PI                                    ((double) 3.14159265358979323846)
#define PARTIAL_AREA_BOUND_TYPE_INDEX         "index"
#define PARTIAL_AREA_BOUND_TYPE_VALUE         "value"
#define NULL_COORD_VALUE                      ((double) -999999.)
#define TOLERABLE_ERROR                       ((double) 0.0000001)
#define COORD_LABEL_LAT                       "lat"
#define COORD_LABEL_LON                       "lon"
#define COORD_LABEL_LEV                       "lev"
#define COORD_LABEL_TIME                      "time"
#define COORD_LABEL_TRACER                    "tracer"
#define COORD_UNIT_RADIANS                    "radians"
#define COORD_UNIT_DEGREES                    "degrees"
#define GRID_FIELD_ATTRIBUTE_UNIT             "unit"
#define SPHERE_GRID_ROTATION_LAT_THRESHOLD    ((double) 70.0)

#define V3D_GRID_3D_LEVEL_FIELD_NAME          "V3D_grid_3D_level_field"

class Runtime_remap_function;
class Remap_grid_data_class;
class Remap_operator_basis;
class Remap_grid_class;
class Remap_operator_grid;


struct Partial_area_bound
{
    Remap_grid_class *area_bound_grid;
    char bound_type[256];
    double bound_start;
    double bound_end;
};


struct Partial_area
{
    char area_name[256];
    std::vector<Partial_area_bound> area_bounds;
};


class Remap_grid_class
{
    private:
        friend class Runtime_remap_function;
        friend class Remap_operator_grid;
        char decomp_name[NAME_STR_SIZE];
        char grid_name[NAME_STR_SIZE];
        char coord_label[NAME_STR_SIZE];
        char coord_unit[NAME_STR_SIZE];
        long grid_size;
        int num_dimensions;
        int num_vertexes;
        Remap_grid_class *whole_grid;
        Remap_grid_class *mid_point_grid;
        Remap_grid_class *interface_level_grid;
        std::vector<Partial_area *> partial_areas;
        Remap_grid_class *super_grid_of_setting_coord_values;
        char name_super_grid_of_setting_coord_values[NAME_STR_SIZE];
        Remap_grid_class *first_super_grid_of_enable_setting_coord_value;
        char name_first_super_grid_of_enable_setting_coord_value[NAME_STR_SIZE];
        Remap_grid_class *duplicated_grid;
        Remap_grid_class *original_grid;
        bool generated_from_duplication;
        std::vector<Remap_grid_class *> sub_grids;
        std::vector<Remap_grid_data_class *> grid_center_fields;
        std::vector<Remap_grid_data_class *> grid_vertex_fields;
        Remap_grid_data_class *grid_mask_field;
        Remap_grid_data_class *imported_area;
        std::vector<Remap_grid_class *> super_grids_of_setting_mask_value; 
        Remap_grid_data_class *original_grid_mask_field;
        bool masks_are_known;        
        bool cyclic;
        bool enable_to_set_coord_values;
        bool are_vertex_values_set_in_default;
        Remap_grid_data_class *redundant_cell_mark_field;
        bool *redundant_cell_mark;
        double *area_or_volumn;
        double boundary_min_lon;
        double boundary_max_lon;
        double boundary_max_lat;
        double boundary_min_lat;
		bool using_V3D_level_coord;
        Remap_grid_data_class *hybrid_grid_coefficient_field;
        Remap_grid_data_class *sigma_grid_sigma_value_field;
        Remap_grid_data_class *level_V3D_coord_trigger_field;
        double sigma_grid_top_value;
        double sigma_grid_scale_factor;
        bool level_V3D_coord_trigger_field_specified;
        Remap_grid_data_class *level_V3D_coord_dynamic_trigger_field;
        

        /* Functions of checking the coordinate values of grid */
        void initialize_grid_class_data();
        void check_grid_field_can_be_set();
        void check_center_coord_value_can_be_set(const char*);
        void check_vertex_coord_value_can_be_set(const char*);
        void check_mask_value_can_be_set();
        bool check_vertex_fields_completeness(Remap_operator_basis*);
        void check_coord_values_range() const;
        void check_center_vertex_values_consistency_2D();
        void check_center_fields_sorting_order();
        void check_and_set_first_super_grid_of_enable_setting_coord_value(Remap_grid_class*);
        void formalize_sphere_grid();

        /* Functions of setting default coordinate values of grid */
        void generate_voronoi_grid();
        void set_1D_coord_vertex_values_in_default(const double*, double*, long, bool, bool);
        void set_2D_coord_vertex_values_in_default(const double*, double*, long, bool, bool, const double*, double*, long, bool, bool);

        /* Functions of formalizing the coordinate values of grid */
        void transform_coord_unit_from_radian_to_degrees();
        void remove_redundant_vertex_values();
        void detect_redundant_cells();
        void formalize_cyclic_coord_values(Remap_grid_data_class*);
        void transform_coord_values_from_radian_to_degrees(Remap_grid_data_class*);

        /* Functions of getting grid properties */
        bool has_partial_sub_grid() const;
        long get_lower_dimension_size_of_leaf_grid(Remap_grid_data_class*, Remap_grid_class*);
        long get_higher_dimension_size_of_leaf_grid(Remap_grid_data_class*, Remap_grid_class*);
        Remap_grid_class *get_similar_grids_setting_coord_values();

        /* Functions of getting grids relation */
        bool is_superset_of_grid(Remap_grid_class *);

        void calculate_area_of_sphere_grid();
        void calculate_area_or_volumn();
        void generate_default_grid_mask();
        void reset_super_grids_of_setting_mask_value(Remap_grid_class *);

    public:

        /* Functions of constructor and destructor */
        Remap_grid_class() {}
        Remap_grid_class(const char*, const char*, const char*, const char*, long);
        Remap_grid_class(const char*, int, Remap_grid_class**, long);
        Remap_grid_class(Remap_grid_class*, Remap_grid_class*, Remap_grid_class*, bool);
        Remap_grid_class(const char*, const char*);
        Remap_grid_class(const char*, const char*, const char*);
        ~Remap_grid_class();

        /* Functions of getting grid properties or data */
        long get_grid_size() const { return grid_size; }
        int get_num_dimensions() const { return num_dimensions; }
        int get_num_vertexes() const { return num_vertexes; }
        const char *get_grid_name() const { return grid_name; }
        const char *get_coord_label() const { return coord_label; }
        const char *get_coord_unit() const { return coord_unit; }
        const char *get_decomp_name() const { return decomp_name; }
        void set_decomp_name(const char*);
        bool get_grid_cyclic() const { return cyclic; }
		void allocate_default_center_field();
        Remap_grid_data_class *get_grid_mask_field() const { return grid_mask_field; }
        Remap_grid_data_class *get_grid_imported_area() const { return imported_area; }
        bool get_are_vertex_values_set_in_default() const { return are_vertex_values_set_in_default; }
        const Remap_grid_class *get_whole_grid() const { return whole_grid; }
        Remap_grid_class *get_super_grid_of_setting_coord_values() const { return super_grid_of_setting_coord_values; }
		void set_super_grid_of_setting_coord_values(Remap_grid_class *super_grid) { super_grid_of_setting_coord_values = super_grid; }
        Remap_grid_class *get_first_super_grid_of_enable_setting_coord_value() { return first_super_grid_of_enable_setting_coord_value; }
        void get_leaf_grids(int *, Remap_grid_class**, const Remap_grid_class*) const;
        Remap_grid_class *get_a_leaf_grid_of_sigma_or_hybrid();
        Remap_grid_class *get_a_leaf_grid(const char*);
        void get_sized_sub_grids(int*, Remap_grid_class**);
        void get_masked_sub_grids(int*, Remap_grid_class**);
        void get_partial_grid_mask_fields(int*, Remap_grid_data_class**, Remap_grid_class*);
        void get_grid_index_interchange_table(Remap_grid_class*, int*);
        bool is_partial_grid() const;
        Remap_grid_class *get_original_grid() { return original_grid; }
        Remap_grid_data_class *get_grid_center_field(const char*) const;
        Remap_grid_data_class *get_grid_center_field() const;
        Remap_grid_data_class *get_grid_vertex_field() const;
        bool has_grid_coord_label(const char*) const;
        bool get_is_sphere_grid() const;
        bool are_all_vertex_fields_specified_by_user();
        double *get_area_or_volumn() { return area_or_volumn; }

        /* Functions of getting grids relation */
        bool match_grid(const char*) const;
        bool match_grid(int, Remap_grid_class**);
        void read_grid_data_from_IO(char extension_names[16][256], const char*, const char*, int);
        void read_grid_data_from_array(const char*, const char*, const char*, const char*, int);
        void read_grid_data_through_span(char[16][256], const char*, const char*, long, const char*);
        void extract_mask(const char*, const char*, const char*);
        void compute_ocn_mask(const char*, double);
        bool have_overlap_with_grid(Remap_grid_class *);
        bool is_subset_of_grid(Remap_grid_class *);
        bool is_the_same_grid_with(Remap_grid_class *);
        bool is_similar_grid_with(Remap_grid_class *);

        /* Functions of checking grid data */
        void end_grid_definition_stage(Remap_operator_basis*);

        /* Functions of inputing or computing grid data */
        void add_partial_grid_area(const char*);
        void add_partial_grid_area_bounds(const char*, const char*, const char*, const char*, const char*);
        void compute_partial_grid_mask();
        Remap_grid_data_class *expand_to_generate_full_coord_value(Remap_grid_data_class*);
        void compute_remap_field_data_runtime_mask(Remap_grid_class*, Remap_grid_class **, int*, Remap_grid_data_class **);

        /* Functions for remaping process */
        void generate_interchange_grids(Remap_grid_class*, Remap_grid_class**, Remap_grid_class**, int);
        Remap_grid_class *generate_remap_operator_runtime_grid(Remap_grid_class*, Remap_operator_basis*, Remap_grid_data_class *);
        void interchange_grid_fields_for_remapping(Remap_grid_class*, Remap_grid_class*, Remap_grid_data_class *);
        Remap_grid_class *duplicate_grid(Remap_grid_class*);
        Remap_grid_class *generate_decomp_grid(const int*, int, const char*);
        void generate_3D_grid_decomp_sigma_values(Remap_grid_class*, Remap_grid_class*, const int*, int);
        void gen_lev_coord_from_sigma_or_hybrid(char extension_names[16][256], const char*, const char*, const char*, const char*, double);
        void calculate_lev_sigma_values();
		void update_grid_center_3D_level_field_from_external();
        bool is_sigma_grid();
		bool does_use_V3D_level_coord();
		void set_using_V3D_level_coord();
        Remap_grid_data_class *get_sigma_grid_sigma_value_field();
        Remap_grid_data_class *get_hybrid_grid_coefficient_field() { return hybrid_grid_coefficient_field; }
        Remap_grid_data_class *get_level_V3D_coord_trigger_field() { return level_V3D_coord_trigger_field; }
        bool is_level_V3D_coord_trigger_field_specified() { return level_V3D_coord_trigger_field_specified; }
        void allocate_sigma_grid_specific_fields(Remap_grid_data_class*, Remap_grid_data_class*, Remap_grid_data_class*, double, double);
        void set_level_V3D_coord_dynamic_trigger_field(Remap_grid_data_class *); 
        Remap_grid_data_class *get_level_V3D_coord_dynamic_trigger_field() { return level_V3D_coord_dynamic_trigger_field; }
        bool is_level_V3D_coord_trigger_field_updated();
        void set_lev_grid_sigma_info(const char*, double, double, const char*);
        void set_lev_grid_sigma_info(double, const double *, const double *, double);
        double get_sigma_grid_top_value() { return sigma_grid_top_value; }
        void renew_lev_grid_coord_values(double*, double*);
        bool has_super_grids_of_setting_mask_value() { return super_grids_of_setting_mask_value.size() > 0; }
        void set_original_grid(Remap_grid_class *grid) { original_grid = grid; }

        /* Function for checking coordinate values consistency with coupler */
        bool check_coord_values_consistency(const char*, const char*, const void*);
        bool check_mask_values_consitency(const char*, const void*);

        void set_grid_boundary(double, double, double, double);

        Remap_grid_data_class *get_unique_center_field();
        Remap_grid_data_class *get_unique_vertex_field();

        void set_coord_vertex_values_in_default();        
        void write_grid_field_into_array(Remap_grid_data_class *, char **, long&, long&);
        void read_grid_field_from_array(Remap_grid_data_class **, const char *, long &);
        Remap_grid_class *get_linked_grid_from_array(Remap_grid_class *, const char *, char *);
        void link_grids(Remap_grid_class *, const char *);
        void write_grid_name_into_array(Remap_grid_class *, char **, long &, long &);
        void write_grid_into_array(char **, long &, long &);
        bool format_sub_grids(Remap_grid_class *);
        bool is_sub_grid_of_grid(Remap_grid_class *);
        Remap_grid_class *search_sub_grid(const char*);
        Remap_grid_class *get_sphere_sub_grid();
        Remap_grid_class(Remap_grid_class*, const char *, const char *, long &);
        Remap_grid_data_class *generate_mid_point_grid_field(Remap_grid_data_class *);
        Remap_grid_class *generate_mid_point_grid();
		Remap_grid_class *get_mid_point_grid() { return mid_point_grid; }
		void set_mid_point_grid(Remap_grid_class *mid_point_grid) { this->mid_point_grid = mid_point_grid; }
        double get_boundary_min_lon() { return boundary_min_lon; }
        double get_boundary_max_lon() { return boundary_max_lon; }
        double get_boundary_min_lat() { return boundary_min_lat; }
        double get_boundary_max_lat() { return boundary_max_lat; }
};

#endif
