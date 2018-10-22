/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "global_data.h"
#include "remap_operator_basis.h"
#include "quick_sort.h"
#include <string.h>


bool *H2D_grid_decomp_mask = NULL;


Remap_operator_basis::Remap_operator_basis()
{
    this->displ_src_cells_overlap_with_dst_cells = NULL;
    this->index_src_cells_overlap_with_dst_cells = NULL;
}


Remap_operator_basis::Remap_operator_basis(const char *object_name, 
                                          const char *operator_name, 
                                          int num_dimensions,
                                          bool is_operator_regriding,
                                          bool require_grid_vertex_values,
                                          bool require_grid_cell_neighborhood,
                                          int num_remap_grids, 
                                          Remap_grid_class **remap_grids)
{
    this->is_operator_regriding = is_operator_regriding;
    this->enable_to_set_parameters = true;
    this->require_grid_vertex_values = require_grid_vertex_values;
    this->require_grid_cell_neighborhood = require_grid_cell_neighborhood;
    this->num_dimensions = num_dimensions;
    this->displ_src_cells_overlap_with_dst_cells = NULL;
    this->index_src_cells_overlap_with_dst_cells = NULL;
    strcpy(this->object_name, object_name);
    strcpy(this->operator_name, operator_name);

    register_remap_grids(num_remap_grids, remap_grids);
}


Remap_operator_basis::~Remap_operator_basis()
{
    for (int i = 0; i < remap_weights_groups.size(); i ++)
        delete remap_weights_groups[i];

    if (displ_src_cells_overlap_with_dst_cells != NULL) {
        delete [] displ_src_cells_overlap_with_dst_cells;
        delete [] index_src_cells_overlap_with_dst_cells;
    }    
}


void Remap_operator_basis::register_remap_grids(int num_remap_grids, Remap_grid_class **remap_grids)
{
    int i, j, num_leaf_grids;
    Remap_grid_class *leaf_grids_src[256], *leaf_grids_dst[256];

    if (is_operator_regriding)
        EXECUTION_REPORT(REPORT_ERROR, -1, num_remap_grids == 2, 
                     "when generating remap operator %s object \"%s\", must input two grid objects where one is source grid and the other is destination grid\n",
                     operator_name, object_name);
    else 
        EXECUTION_REPORT(REPORT_ERROR, -1, num_remap_grids == 1, 
                     "when generating operator %s object \"%s\", there must be only one input grid object as %s does not regrid\n",
                     operator_name, object_name, operator_name);

    src_grid = remap_grids[0];
    if (num_remap_grids == 2)
        dst_grid = remap_grids[1];
    else dst_grid = remap_grids[0];
 
    for (i = 0; i < num_remap_grids; i ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_grids[i]->get_num_dimensions() == num_dimensions, 
                     "the dimension number of grid object \"%s\" is different from the implicit dimension number of remap operator %s\n",
                     remap_grids[i]->get_grid_name(), operator_name);    
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_grids[i]->get_grid_size()> 0, 
                     "the size of grid object involved in remaping must be specified before building remap operators. However, the size of grid object \"%s\" have not been specified explicitly or implicitly\n",
                     remap_grids[i]->get_grid_name());        
        EXECUTION_REPORT(REPORT_ERROR, -1, !remap_grids[i]->is_partial_grid(),
                     "\"%s\" is a partial grid. It can not be used as source or destination grid of remap operator\n",
                     remap_grids[i]->get_grid_name());
    }

    remap_grids[0]->get_leaf_grids(&num_leaf_grids, leaf_grids_src, remap_grids[0]);
    if (num_remap_grids == 2) {
        remap_grids[1]->get_leaf_grids(&num_leaf_grids, leaf_grids_dst, remap_grids[1]);
        for (i = 0; i < num_leaf_grids; i ++)
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(leaf_grids_src[i]->get_coord_label(), leaf_grids_dst[i]->get_coord_label()),
                         "the arrangement of coordinate labels of grid \"%s\" and \"%s\" must be the same\n",
                         remap_grids[0]->get_grid_name(),
                         remap_grids[1]->get_grid_name());            
    }

    if (get_is_sphere_grid() &&
        words_are_the_same(leaf_grids_src[0]->get_coord_label(), COORD_LABEL_LAT) && 
        words_are_the_same(leaf_grids_src[1]->get_coord_label(), COORD_LABEL_LON)) {
        if (num_remap_grids == 2)
            EXECUTION_REPORT(REPORT_ERROR, -1, false,
                         "for the 2D remapping of sphere grid, the coordinate labels of grid \"%s\" and \"%s\" must be in the order of \"lon\" - \"lat\"\n",
                         remap_grids[0]->get_grid_name(),
                         remap_grids[1]->get_grid_name());
        else 
            EXECUTION_REPORT(REPORT_ERROR, -1, false,
                         "for the 2D remapping of sphere grid, the coordinate labels of grid \"%s\" must be in the order of \"lon\" - \"lat\"\n",
                         remap_grids[0]->get_grid_name(),
                         remap_grids[1]->get_grid_name());
    }
}


bool Remap_operator_basis::match_remap_operator(const char *object_name)
{
    return words_are_the_same(this->object_name, object_name);
}


bool Remap_operator_basis::match_remap_operator(Remap_grid_class *grid_src, Remap_grid_class *grid_dst, const char *operator_name)
{
    return words_are_the_same(this->operator_name, operator_name) && this->src_grid == grid_src && this->dst_grid == grid_dst;
}


void Remap_operator_basis::calculate_grids_overlaping()
{
    double center_coord_values_dst[2];
    int num_vertexes_dst, num_grid_dimensions_dst, i;
    long cell_index_src, overlapping_src_cells_indexes[409600];
    int num_overlapping_src_cells;
    long temp_array_iter, *temp_array, cell_index_dst;
    bool dst_cell_mask;


    if (displ_src_cells_overlap_with_dst_cells == NULL) {
        displ_src_cells_overlap_with_dst_cells = new long [dst_grid->get_grid_size()+1];
        index_src_cells_overlap_with_dst_cells = new long [dst_grid->get_grid_size()];
        size_index_src_cells_overlap_with_dst_cells = dst_grid->get_grid_size();
    }
    temp_array_iter = 0;
    num_grid_dimensions_dst = current_runtime_remap_operator_grid_src->get_num_grid_dimensions();

    for (i = 0; i < dst_grid->get_grid_size(); i ++) 
        displ_src_cells_overlap_with_dst_cells[i] = 0;

    for (cell_index_dst = 0; cell_index_dst < dst_grid->get_grid_size(); cell_index_dst ++) {
        finalize_computing_remap_weights_of_one_cell();
        initialize_computing_remap_weights_of_one_cell();
        get_cell_mask_of_dst_grid(cell_index_dst, &dst_cell_mask);
        if (!dst_cell_mask)
            continue;
        get_cell_center_coord_values_of_dst_grid(cell_index_dst, center_coord_values_dst);  
        get_current_grid2D_search_engine(true)->search_overlapping_cells(num_overlapping_src_cells, overlapping_src_cells_indexes, get_current_grid2D_search_engine(false)->get_cell(cell_index_dst), true, false);
        displ_src_cells_overlap_with_dst_cells[cell_index_dst] = num_overlapping_src_cells;
        if (num_overlapping_src_cells+temp_array_iter >= size_index_src_cells_overlap_with_dst_cells) {
            size_index_src_cells_overlap_with_dst_cells *= 2;
            temp_array = new long [size_index_src_cells_overlap_with_dst_cells];
            for (i = 0; i < temp_array_iter+1; i ++)
                temp_array[i] = index_src_cells_overlap_with_dst_cells[i];
            delete [] index_src_cells_overlap_with_dst_cells;
            index_src_cells_overlap_with_dst_cells = temp_array;
        }
        for (i = 0; i < num_overlapping_src_cells; i ++)
            index_src_cells_overlap_with_dst_cells[temp_array_iter++] = overlapping_src_cells_indexes[i];
        do_quick_sort(overlapping_src_cells_indexes, (long*)NULL, 0, num_overlapping_src_cells-1);
    }
    finalize_computing_remap_weights_of_one_cell();

    displ_src_cells_overlap_with_dst_cells[dst_grid->get_grid_size()] = temp_array_iter;
    for (i = dst_grid->get_grid_size()-1; i >= 0; i --) {
        temp_array_iter -= displ_src_cells_overlap_with_dst_cells[i];
        displ_src_cells_overlap_with_dst_cells[i] = temp_array_iter;
    }
}


void Remap_operator_basis::copy_remap_operator_basic_data(Remap_operator_basis *another_remap_operator, bool fully_copy)
{
    long i;


    strcpy(another_remap_operator->object_name, this->object_name);
    strcpy(another_remap_operator->operator_name, this->operator_name);
    another_remap_operator->src_grid = this->src_grid;
    another_remap_operator->dst_grid = this->dst_grid;
    another_remap_operator->num_dimensions = this->num_dimensions;
    another_remap_operator->enable_to_set_parameters = this->enable_to_set_parameters;
    another_remap_operator->is_operator_regriding = this->is_operator_regriding;
    another_remap_operator->require_grid_vertex_values = this->require_grid_vertex_values;
    another_remap_operator->require_grid_cell_neighborhood = this->require_grid_cell_neighborhood;
    another_remap_operator->displ_src_cells_overlap_with_dst_cells = NULL;
    another_remap_operator->index_src_cells_overlap_with_dst_cells = NULL;

    if (fully_copy)
        for (i = 0; i < this->remap_weights_groups.size(); i ++)
            another_remap_operator->remap_weights_groups.push_back(this->remap_weights_groups[i]->duplicate_remap_weight_of_sparse_matrix());
}


void Remap_operator_basis::generate_parallel_remap_weights(Remap_operator_basis *another_remap_operator, Remap_grid_class **decomp_original_grids, int **global_cells_local_indexes_in_decomps)
{
    for (int i = 0; i < this->remap_weights_groups.size(); i ++)
        another_remap_operator->remap_weights_groups.push_back(this->remap_weights_groups[i]->generate_parallel_remap_weight_of_sparse_matrix(decomp_original_grids, global_cells_local_indexes_in_decomps));
}


void Remap_operator_basis::change_remap_operator_info(const char *operator_name, Remap_grid_class *grid_src, Remap_grid_class *grid_dst)
{
    strcpy(this->operator_name, operator_name);
    this->src_grid = grid_src;
    this->dst_grid = grid_dst;
}


void Remap_operator_basis::update_unique_weight_sparse_matrix(Remap_weight_sparse_matrix *new_sparse_matrix)
{
    if (remap_weights_groups.size() > 0) {
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights_groups.size() == 1 && remap_weights_groups[0]->get_num_weights() == 0, "Software error in Remap_operator_basis::update_unique_weight_sparse_matrix");
        delete remap_weights_groups[0];
        remap_weights_groups.clear();
    }

    remap_weights_groups.push_back(new_sparse_matrix);
}


Remap_operator_basis *Remap_operator_basis::gather(int comp_id)
{
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	Remap_operator_basis *overall_remap_operator = NULL;

	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_weights_groups.size() == 1, "software error in Remap_operator_basis::gather");
	if (comp_node->get_current_proc_local_id() == 0)
		overall_remap_operator = duplicate_remap_operator(false);

	for (int i = 0; i < remap_weights_groups.size(); i ++) {
		Remap_weight_sparse_matrix *overall_remap_weight_sparse_matrix = remap_weights_groups[i]->gather(comp_id);
		if (comp_node->get_current_proc_local_id() == 0)
			overall_remap_operator->remap_weights_groups.push_back(overall_remap_weight_sparse_matrix);
	}
	
	return overall_remap_operator;
}

