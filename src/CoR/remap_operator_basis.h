/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_OPERATOR_BASIS
#define REMAP_OPERATOR_BASIS

#include "remap_grid_class.h"
#include "remap_operator_c_interface.h"
#include "remap_weight_sparse_matrix.h"
#include "remap_weight_of_strategy_class.h"
#include <vector>


#define REMAP_OPERATOR_NAME_BILINEAR               "bilinear"
#define REMAP_OPERATOR_NAME_CONSERV_2D             "conserv_2D"
#define REMAP_OPERATOR_NAME_SPLINE_1D              "spline_1D"
#define REMAP_OPERATOR_NAME_DISTWGT                "distwgt"
#define REMAP_OPERATOR_NAME_LINEAR                 "linear"
#define REMAP_OPERATOR_NAME_SMOOTH                 "smooth"
#define REMAP_OPERATOR_NAME_REGRID                 "regrid"


class Remap_operator_basis
{
    private:
        void register_remap_grids(int, Remap_grid_class **);
        
    protected:
        friend class Remap_weight_of_operator_instance_class;
        friend class Remap_weight_of_strategy_class;
        char object_name[256];
        char operator_name[256];
        Remap_grid_class *src_grid;
        Remap_grid_class *dst_grid;
        int num_dimensions;
        bool enable_to_set_parameters;
        bool is_operator_regriding;
        bool require_grid_vertex_values;
        bool require_grid_cell_neighborhood;
        std::vector<Remap_weight_sparse_matrix*> remap_weights_groups;
        long *displ_src_cells_overlap_with_dst_cells;
        long *index_src_cells_overlap_with_dst_cells;
        long size_index_src_cells_overlap_with_dst_cells;
        bool enable_extrapolate;

    public:
        Remap_operator_basis(const char*, const char*, int, bool, bool, bool, int, Remap_grid_class **);
        Remap_operator_basis();
        virtual ~Remap_operator_basis();
        virtual void set_parameter(const char*, const char*) = 0;
        virtual int check_parameter(const char*, const char*, char*) = 0;
        virtual void do_remap_values_caculation(double*, double*, int) = 0;
        virtual void do_src_decomp_caculation(long*, const long*) = 0;
        virtual void calculate_remap_weights() = 0;
        virtual Remap_operator_basis *duplicate_remap_operator(bool) = 0;
        virtual Remap_operator_basis *generate_parallel_remap_operator(Remap_grid_class**, int**) = 0;
        virtual void compute_remap_weights_of_one_dst_cell(long) = 0;
        bool match_remap_operator(const char*);
        bool match_remap_operator(Remap_grid_class*, Remap_grid_class*, const char*);
        void calculate_grids_overlaping();
        void copy_remap_operator_basic_data(Remap_operator_basis*, bool);
        void generate_parallel_remap_weights(Remap_operator_basis*, Remap_grid_class**, int**);
        Remap_grid_class *get_src_grid() { return src_grid; }
        Remap_grid_class *get_dst_grid() { return dst_grid; }
        char *get_object_name() { return object_name; }
        void disable_to_set_parameters() { enable_to_set_parameters = false; }
        bool get_is_operator_regridding() { return is_operator_regriding; }
        bool does_require_grid_vertex_values() { return require_grid_vertex_values; }
        bool does_require_grid_cell_neighborhood() { return require_grid_cell_neighborhood; }
        Remap_weight_sparse_matrix *get_remap_weights_group(int index) { return remap_weights_groups[index]; }
        int get_num_remap_weights_groups() { return remap_weights_groups.size(); }
        const char *get_operator_name() { return operator_name; }
        int get_num_dimensions() { return num_dimensions; }
        bool get_is_sphere_grid() { return src_grid->get_is_sphere_grid(); }
        void add_weight_sparse_matrix(Remap_weight_sparse_matrix *sparse_matrix) { remap_weights_groups.push_back(sparse_matrix); }
        void update_unique_weight_sparse_matrix(Remap_weight_sparse_matrix *);
        void change_remap_operator_info(const char*, Remap_grid_class*, Remap_grid_class*);
        void set_src_grid(Remap_grid_class *new_src_grid) { src_grid = new_src_grid; }
        void set_dst_grid(Remap_grid_class *new_dst_grid) { dst_grid = new_dst_grid; }
		Remap_operator_basis *gather(int);
};


extern bool *H2D_grid_decomp_mask;


#endif
