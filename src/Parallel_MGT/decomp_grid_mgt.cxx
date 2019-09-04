/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "decomp_grid_mgt.h"
#include "cor_cpl_interface.h"


Decomp_grid_info::Decomp_grid_info(int decomp_id, Remap_grid_class *original_grid)
{
    Decomp_info *decomp;
    Remap_grid_class *decomp_info_grid, *decomp_2D_grid;
    Remap_grid_class *leaf_grids[256], *sub_grids[256];
    int num_leaf_grids, num_sub_grids, i, j, decomp_leaf_grid_indexes[2];
    char decomp_grid_name[NAME_STR_SIZE];
    int comp_id;


    EXECUTION_REPORT(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(decomp_id), "Software error in Decomp_grid_info::Decomp_grid_info(int decomp_id)");
    decomp = decomps_info_mgr->get_decomp_info(decomp_id);
    decomp_info_grid = original_grid_mgr->get_original_CoR_grid(decomp->get_grid_id());
    comp_id = decomp->get_comp_id();
    this->decomp_id = decomp_id;
    strcpy(this->decomp_name, decomp->get_decomp_name());
    this->original_grid = original_grid;
    EXECUTION_REPORT_LOG(REPORT_LOG, decomp->get_host_comp_id(), true, 
                     "Generate decomposition grid for the grid \"%s\" on the parallel decomposition \"%s\"", 
                     original_grid->get_grid_name(), decomp_name);

    if (decomp->get_num_local_cells() == 0) {
        this->decomp_grid = NULL;
        return;
    }

    if (this->original_grid->get_is_sphere_grid()) {
        EXECUTION_REPORT(REPORT_ERROR, -1, decomp_info_grid == original_grid, "Software error in Decomp_grid_info::Decomp_grid_info: inconsistent H2D grid: %s  %s  %s", original_grid->get_grid_name(), decomp->get_grid_name(), decomp_name);
        EXECUTION_REPORT_LOG(REPORT_LOG, decomp->get_host_comp_id(), true, 
                         "generate decomposition sphere grid (%s %s) with size %d", 
                         decomp_name, original_grid->get_grid_name(), decomp->get_num_local_cells());
        this->decomp_grid = this->original_grid->generate_decomp_grid(decomp->get_local_cell_global_indx(), decomp->get_num_local_cells(), decomp_name);
    }
    else {
        if (original_grid->has_grid_coord_label(COORD_LABEL_LON) || original_grid->has_grid_coord_label(COORD_LABEL_LAT))
            EXECUTION_REPORT(REPORT_ERROR, -1, decomp_info_grid->is_subset_of_grid(original_grid), "Software error in Decomp_grid_info::Decomp_grid_info: inconsistent H2D and MD grid",
                             original_grid->get_grid_name(), decomp_info_grid->get_grid_name(), decomp_name);
        original_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, original_grid);
        for (i = 0, j = 0; i < num_leaf_grids; i ++)
            if (leaf_grids[i]->is_subset_of_grid(decomp_info_grid))
                decomp_leaf_grid_indexes[j++] = i;
        EXECUTION_REPORT(REPORT_ERROR,-1, j == 2, "Software error in Decomp_grid_info::Decomp_grid_info: wrong j"); 
        EXECUTION_REPORT(REPORT_ERROR,-1, decomp_leaf_grid_indexes[1]-decomp_leaf_grid_indexes[0] == 1, "C-Coupler error in Decomp_grid_info::Decomp_grid_info\n");
        for (i = 0, num_sub_grids = 0; i < num_leaf_grids; i ++) {
            if (leaf_grids[i] == NULL)
                continue;
            if (leaf_grids[i]->is_subset_of_grid(decomp_info_grid)) {
                for (j = i+1; j < num_leaf_grids; j ++)
                    if (leaf_grids[j]->is_subset_of_grid(decomp_info_grid))
                        leaf_grids[j] = NULL;
                decomp_2D_grid = decomp_grids_mgr->search_decomp_grid_info(decomp_id, decomp_info_grid, false)->get_decomp_grid();
                sub_grids[num_sub_grids++] = decomp_2D_grid;
            }
            else {
                sub_grids[num_sub_grids++] = leaf_grids[i]->duplicate_grid(leaf_grids[i]); 
                remap_grid_manager->add_temp_grid(sub_grids[num_sub_grids-1]);
            }
        }
        sprintf(decomp_grid_name, "%s_at_DECOMP_%s", original_grid->get_grid_name(), decomp->get_decomp_name());
        EXECUTION_REPORT_LOG(REPORT_LOG, decomp->get_host_comp_id(), true, "Generate decomp grid \"%s\" based on H2D grid \"%s\" %lx", decomp_grid_name, decomp_2D_grid->get_grid_name(), decomp_2D_grid);
        this->decomp_grid = new Remap_grid_class(decomp_grid_name, num_sub_grids, sub_grids, 0);
        this->decomp_grid->set_decomp_name(decomp_name);
        this->decomp_grid->set_original_grid(original_grid);
        if (original_grid->has_grid_coord_label(COORD_LABEL_LEV))
            this->decomp_grid->generate_3D_grid_decomp_sigma_values(original_grid, decomp_2D_grid, decomp->get_local_cell_global_indx(), decomp->get_num_local_cells());
		if (original_grid->does_use_V3D_level_coord())
			this->decomp_grid->set_using_V3D_level_coord();
    }
}


Decomp_grid_info::~Decomp_grid_info()
{
    if (decomp_grid != NULL && decomp_grid != original_grid) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "remove decomp grid %s\n", decomp_grid->get_grid_name());
        delete decomp_grid;
    }
}


Decomp_grid_info *Decomp_grid_mgt::search_decomp_grid_info(int decomp_id, Remap_grid_class *original_grid, bool diag)
{
    for (int i = 0; i < decomp_grids_info.size(); i ++)
        if (decomp_grids_info[i]->get_decomp_id() == decomp_id && decomp_grids_info[i]->get_original_grid() == original_grid)
            return decomp_grids_info[i];

    if (diag)
        EXECUTION_REPORT(REPORT_ERROR, -1, true, "C-Coupler error in Decomp_grid_mgt::search_decomp_grid_info");
        
    decomp_grids_info.push_back(new Decomp_grid_info(decomp_id, original_grid));
    return decomp_grids_info[decomp_grids_info.size()-1];
}


void Decomp_grid_mgt::set_decomp_grids_using_3D_level_coord(Remap_grid_class *original_V3D_grid)
{
	for (int i = 0; i < decomp_grids_info.size(); i ++)
		if (decomp_grids_info[i]->get_original_grid() == original_V3D_grid)
			decomp_grids_info[i]->get_decomp_grid()->set_using_V3D_level_coord();
}


Decomp_grid_mgt::~Decomp_grid_mgt()
{
    for (int i = 0; i < decomp_grids_info.size(); i ++)
        delete decomp_grids_info[i];
}


Remap_grid_class *Decomp_grid_mgt::search_decomp_grid_original_grid(int decomp_id, Remap_grid_class *decomp_grid)
{
	for (int i = 0; i < decomp_grids_info.size(); i ++)
        if (decomp_grids_info[i]->get_decomp_id() == decomp_id && decomp_grids_info[i]->get_decomp_grid() == decomp_grid)
            return decomp_grids_info[i]->get_original_grid();	
		
	EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Decomp_grid_mgt::search_decomp_grid_original_grid");
	return NULL;
}


