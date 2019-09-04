/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_grid_class.h"
#include "remap_grid_data_class.h"
#include "remap_common_utils.h"
#include "remap_operator_grid.h"
#include "parse_special_words.h"
#include "remap_operator_regrid.h"
#include "radix_sort.h"
#include "execution_report.h"
#include "delaunay_voronoi.h"
#include <math.h>
#include <string.h>


void Remap_grid_class::initialize_grid_class_data()
{
    strcpy(this->grid_name, "\0");
    strcpy(this->coord_label, "\0");
    strcpy(this->coord_unit, "\0");
    strcpy(this->decomp_name, "\0");
    this->grid_size = 0;
    this->num_dimensions = 0;
    this->num_vertexes = 0;
    this->whole_grid = NULL;
    this->partial_areas.clear();
    this->super_grid_of_setting_coord_values = NULL;
    this->first_super_grid_of_enable_setting_coord_value = NULL;
    this->duplicated_grid = NULL;
    this->generated_from_duplication = false;
    this->original_grid = NULL;
    this->sub_grids.clear();
    this->grid_center_fields.clear();
    this->grid_vertex_fields.clear();
    this->grid_mask_field = NULL;
    this->imported_area = NULL;
    this->super_grids_of_setting_mask_value.clear();
    this->original_grid_mask_field = NULL;
    this->masks_are_known = false;    
    this->cyclic = false;
    this->enable_to_set_coord_values = true;    
    this->are_vertex_values_set_in_default = false;
    this->redundant_cell_mark = NULL;
    this->redundant_cell_mark_field = NULL;
    this->area_or_volumn = NULL;
    this->sigma_grid_scale_factor = 0;
    this->sigma_grid_top_value = 0;
	this->using_V3D_level_coord = false;
    this->boundary_min_lon = NULL_COORD_VALUE;
    this->boundary_max_lon = NULL_COORD_VALUE;
    this->boundary_min_lat = NULL_COORD_VALUE;
    this->boundary_max_lat = NULL_COORD_VALUE;
    this->hybrid_grid_coefficient_field = NULL;
    this->sigma_grid_sigma_value_field = NULL;
    this->level_V3D_coord_trigger_field = NULL;
    this->level_V3D_coord_trigger_field_specified = false;
    this->level_V3D_coord_dynamic_trigger_field = NULL;
    this->mid_point_grid = NULL;
    this->interface_level_grid = NULL;
}


/* Constructor for 1D grid */
Remap_grid_class::Remap_grid_class(const char *grid_name, 
                                   const char *coord_label, 
                                   const char *coord_unit, 
                                   const char *cyclic_or_acyclic,
                                   long grid_size)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(coord_label, COORD_LABEL_LON) || words_are_the_same(coord_label, COORD_LABEL_LAT) ||
                 words_are_the_same(coord_label, COORD_LABEL_LEV) || words_are_the_same(coord_label, COORD_LABEL_TIME) ||
                 words_are_the_same(coord_label, COORD_LABEL_TRACER),
                 "the coordinate name of 1D grid must be one of lon, lat, lev, time and tracer\n");
    if (words_are_the_same(coord_unit, COORD_UNIT_DEGREES) || words_are_the_same(coord_unit, COORD_UNIT_RADIANS))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(coord_label, COORD_LABEL_LAT) || words_are_the_same(coord_label, COORD_LABEL_LON), 
                     "the coordinate label of grid \"%s\" does not match the unit. The coordinate label must be \"lon\" or \"lat\" when coordinate unit is \"degrees\" or \"radians\"\n",
                     grid_name);
    if (words_are_the_same(coord_label, COORD_LABEL_LON) || words_are_the_same(coord_label, COORD_LABEL_LAT))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(coord_unit, COORD_UNIT_DEGREES) || words_are_the_same(coord_unit, COORD_UNIT_RADIANS), 
                     "the coordinate label of grid \"%s\" does not match the unit. The coordinate unit must be \"degrees\" or \"radians\" when coordinate label is \"lon\" or \"lat\"\n",
                     grid_name);
    if (words_are_the_same(coord_label, COORD_LABEL_LON))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(cyclic_or_acyclic, COORD_BOUND_CYCLIC) || words_are_the_same(cyclic_or_acyclic, COORD_BOUND_ACYCLIC),
                     "the cyclic label must be \"cyclic\" or \"acyclic\"\n");
    if (words_are_the_same(coord_label, COORD_LABEL_LEV) || words_are_the_same(coord_label, COORD_LABEL_TIME) || words_are_the_same(coord_label, COORD_LABEL_TRACER))
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_size > 0, "to define the 1D grid of %s, the grid size must be known\n", coord_label);

    initialize_grid_class_data();
    strcpy(this->grid_name, grid_name);
    strcpy(this->coord_label, coord_label);
    strcpy(this->coord_unit, coord_unit); 
    this->grid_size = grid_size;
    this->num_dimensions = 1;

    if (grid_size > 0)
        this->first_super_grid_of_enable_setting_coord_value = this;
    if (words_are_the_same(coord_label, COORD_LABEL_LON))
        if (words_are_the_same(cyclic_or_acyclic, COORD_BOUND_CYCLIC))
            this->cyclic = true;
        else if (words_are_the_same(cyclic_or_acyclic, COORD_BOUND_ACYCLIC))
            this->cyclic = false;
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "the cyclic label for coordinate lon must be \"cyclic\" or \"acyclic\"\n");

	if (get_is_sphere_grid())
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Allocate H2D grid \"%s\" with size %ld", grid_name, grid_size);
}


/* Constructor for multi-dimension grid */
Remap_grid_class::Remap_grid_class(const char *grid_name, 
                                   int num_sub_grids, 
                                   Remap_grid_class **sub_grids, 
                                   long grid_size_by_user)
{
    long grid_size_determined_by_sub_grid = 1;
    int i, j;
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];


    initialize_grid_class_data();
    strcpy(this->grid_name, grid_name);

    for (i = 0; i < num_sub_grids; i ++) {
        if (i > 0)
            EXECUTION_REPORT(REPORT_ERROR, -1, sub_grids[i-1]->get_grid_size() == 0 && sub_grids[i]->get_grid_size() == 0 ||
                         sub_grids[i-1]->get_grid_size() > 0 && sub_grids[i]->get_grid_size() > 0, 
                         "the size of sub grids must be known or known at the same time\n");
        grid_size_determined_by_sub_grid *= sub_grids[i]->get_grid_size();
        this->sub_grids.push_back(sub_grids[i]);
        this->num_dimensions += sub_grids[i]->get_num_dimensions();
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, grid_size_determined_by_sub_grid > 0 || grid_size_by_user > 0,
                 "the size of multi dimensional grid \"%s\" must be known: specified by user or computed from sub grids\n",
                 grid_name);
    if (grid_size_determined_by_sub_grid > 0 && grid_size_by_user > 0)
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_size_by_user == grid_size_determined_by_sub_grid,
                     "the size of combined grid \"%s\" is specified by input parameter, but it is different from the multiple of sub grid sizes, please check\n",
                     grid_name);
    this->grid_size = grid_size_by_user > 0? grid_size_by_user : grid_size_determined_by_sub_grid;

    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++)
        for (j = i+1; j < num_leaf_grids; j ++)
            EXECUTION_REPORT(REPORT_ERROR, -1, !words_are_the_same(leaf_grids[i]->get_coord_label(), leaf_grids[j]->get_coord_label()), 
                         "coordinate label conflict: sub grids \"%s\" and \"%s\" have the same coordinate label \"%s\"\n",
                         leaf_grids[i]->get_grid_name(), leaf_grids[j]->get_grid_name(), leaf_grids[i]->get_coord_label()); 

    this->masks_are_known = true;
    for (i = 0; i < this->sub_grids.size(); i ++)
        if (!this->sub_grids[i]->masks_are_known)
            this->masks_are_known = false;

    for (i = 0; i < num_leaf_grids; i ++)
        leaf_grids[i]->check_and_set_first_super_grid_of_enable_setting_coord_value(this);

	if (get_is_sphere_grid())
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Allocate H2D grid \"%s\" with size %ld", grid_name, grid_size);
}


/* Constructor to generate the grid of destination field data of remap operator */
Remap_grid_class::Remap_grid_class(Remap_grid_class *field_data_grid, 
                                       Remap_grid_class *remap_src_grid, 
                                       Remap_grid_class *remap_dst_grid, 
                                       bool is_remap_operator_regridding)
{
    int num_leaf_grids_src, num_leaf_grids_dst, num_leaf_grids_field_data, num_sized_grids_field_data;
    Remap_grid_class *leaf_grids_src[256], *leaf_grids_dst[256], *leaf_grids_field_data[256], *sized_grids_field_data[256];
    int i, j;


    remap_src_grid->get_leaf_grids(&num_leaf_grids_src, leaf_grids_src, remap_src_grid);
    remap_dst_grid->get_leaf_grids(&num_leaf_grids_dst, leaf_grids_dst, remap_dst_grid);
    field_data_grid->get_leaf_grids(&num_leaf_grids_field_data, leaf_grids_field_data, field_data_grid);

    for (i = 0; i < num_leaf_grids_src; i ++)
        for (j = 0; j < num_leaf_grids_field_data; j ++) 
            if (leaf_grids_field_data[j] == leaf_grids_src[i]) {
                leaf_grids_field_data[j] = leaf_grids_dst[i];
                break;
            }

    initialize_grid_class_data();
    this->num_dimensions = num_leaf_grids_field_data;
    this->grid_size = field_data_grid->grid_size / remap_src_grid->grid_size * remap_dst_grid->grid_size;
    this->enable_to_set_coord_values = false;
    sprintf(this->grid_name, "REMAP_DST %s from %s to %s\0", field_data_grid->grid_name, remap_src_grid->grid_name, remap_dst_grid->grid_name);
    for (i = 0; i < num_leaf_grids_field_data; i ++)
        sub_grids.push_back(leaf_grids_field_data[i]);
    get_sized_sub_grids(&num_sized_grids_field_data, sized_grids_field_data);
    sub_grids.clear();
    for (i = 0; i < num_sized_grids_field_data; i ++)
        sub_grids.push_back(sized_grids_field_data[i]);

    if (is_sigma_grid())
        allocate_sigma_grid_specific_fields(NULL, NULL, NULL, 0, 0);

	if (get_is_sphere_grid())
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Allocate H2D grid \"%s\" with size %ld", grid_name, grid_size);
}


/* constructor for partial grid */
Remap_grid_class::Remap_grid_class(const char *grid_name, const char *whole_grid_name)
{
    Remap_grid_class *whole_grid;


    whole_grid = remap_grid_manager->search_remap_grid_with_grid_name(whole_grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, !whole_grid->is_partial_grid(), "\"%s\" is a partial grid or has partial sub grids. It can not be used to generate another partial grid\n", whole_grid_name);
    
    initialize_grid_class_data();
    strcpy(this->grid_name, grid_name);
    this->whole_grid = whole_grid;
    this->enable_to_set_coord_values = false;
    this->num_dimensions = whole_grid->num_dimensions;
    this->grid_size = whole_grid->grid_size;    
    whole_grid->end_grid_definition_stage(NULL);
}


Remap_grid_class::~Remap_grid_class()
{
	if (get_is_sphere_grid())
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Deallocate H2D grid \"%s\" with size %ld", grid_name, grid_size);

    for (int i = 0; i < grid_center_fields.size(); i ++)
        delete grid_center_fields[i];
    for (int i = 0; i < grid_vertex_fields.size(); i ++)
        delete grid_vertex_fields[i];
    if (grid_mask_field != NULL)
        delete grid_mask_field;

    if (generated_from_duplication)
        for (int i = 0; i < sub_grids.size(); i ++)
            delete sub_grids[i];

    for (int i = 0; i < partial_areas.size(); i ++)
        delete partial_areas[i];

    if (redundant_cell_mark_field != NULL)
        delete redundant_cell_mark_field;

    if (area_or_volumn != NULL)
        delete [] area_or_volumn;

    if (sigma_grid_sigma_value_field != NULL)
        delete sigma_grid_sigma_value_field;

    if (hybrid_grid_coefficient_field != NULL)
        delete hybrid_grid_coefficient_field;

    if (level_V3D_coord_trigger_field != NULL)
        delete level_V3D_coord_trigger_field;
}


void Remap_grid_class::check_and_set_first_super_grid_of_enable_setting_coord_value(Remap_grid_class *super_grid)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 1 && super_grid->get_grid_size() > 0, 
        "remap software error: only one-dimension grid can check and set super grid of setting coordinate value size\n");

    if (first_super_grid_of_enable_setting_coord_value == NULL)
        first_super_grid_of_enable_setting_coord_value = super_grid;
    else EXECUTION_REPORT(REPORT_ERROR, -1, first_super_grid_of_enable_setting_coord_value->is_subset_of_grid(super_grid), 
                      "%s is the first super grid with known grid size of grid %s. However, it is not a sub grid of %s\n",
                      first_super_grid_of_enable_setting_coord_value->get_grid_name(), grid_name, super_grid->get_grid_name());
}


void Remap_grid_class::generate_interchange_grids(Remap_grid_class *remap_grid, 
                                                Remap_grid_class **interchange_grid, 
                                                Remap_grid_class **runtime_mask_sub_grids,
                                                int num_runtime_mask_sub_grids)
{
    int num_leaf_grids_remap;
    Remap_grid_class *leaf_grids_remap[256];    
    int num_size_grids_this, num_size_grids_remap, num_size_grids_interchange;
    Remap_grid_class *size_grids_this[256], *size_grids_remap[256], *size_grids_interchange[256];   
    Remap_grid_class *existing_grid;
    char new_grid_name[256];
    int i, j, k;


    this->get_sized_sub_grids(&num_size_grids_this, size_grids_this);
    for (i = 0; i < num_size_grids_this; i ++)
        if (size_grids_this[i]->have_overlap_with_grid(remap_grid))
            EXECUTION_REPORT(REPORT_ERROR, -1, size_grids_this[i]->is_subset_of_grid(remap_grid), "remap software error1 when generating interchange grids\n");

    remap_grid->get_sized_sub_grids(&num_size_grids_remap, size_grids_remap);
    remap_grid->get_leaf_grids(&num_leaf_grids_remap, leaf_grids_remap, remap_grid);

    /* Order the sized grid in the interchange grid: sized grids of remap operator are in the lowest dimensions, 
           the sized grids which the sized grids of remap operator depend on are in the middle dimensions, and the
           remained sized grids are in the highest dimensions */
    for (i = 0, num_size_grids_interchange = 0; i < num_size_grids_remap; i ++)
        size_grids_interchange[num_size_grids_interchange++] = size_grids_remap[i];
    for (i = 0; i < num_size_grids_this; i ++) {
        if (size_grids_this[i]->is_subset_of_grid(remap_grid))
            continue;
        for (j = 0; j < num_runtime_mask_sub_grids; j ++)
            if (size_grids_this[i]->is_subset_of_grid(runtime_mask_sub_grids[j]) &&
                remap_grid->have_overlap_with_grid(runtime_mask_sub_grids[j])) 
                    break;
        for (k = 0; k < num_leaf_grids_remap; k ++)
            if (leaf_grids_remap[k]->get_super_grid_of_setting_coord_values() != NULL 
                && size_grids_this[i]->is_subset_of_grid(leaf_grids_remap[k]->get_super_grid_of_setting_coord_values()))
                break;
        if (j < num_runtime_mask_sub_grids || k < num_leaf_grids_remap) 
            size_grids_interchange[num_size_grids_interchange++] = size_grids_this[i];
    }    
    for (i = 0; i < num_size_grids_this; i ++) {
        for (j = 0; j < num_size_grids_interchange; j ++)
            if (size_grids_this[i] == size_grids_interchange[j])
                break;
        if (j == num_size_grids_interchange)
            size_grids_interchange[num_size_grids_interchange++] = size_grids_this[i];
    }

    *interchange_grid = new Remap_grid_class(grid_name, num_size_grids_interchange, size_grids_interchange, 0);
    existing_grid = remap_grid_manager->search_same_remap_grid(*interchange_grid);
    if (existing_grid != NULL) {
        delete *interchange_grid;
        *interchange_grid = existing_grid;
    }
    else {
		remap_grid_manager->add_temp_grid(*interchange_grid);
		sprintf(new_grid_name, "%s_INTERCHANGE", grid_name);
		(*interchange_grid)->get_leaf_grids(&num_leaf_grids_remap, leaf_grids_remap, *interchange_grid);
		for (i = 0; i < num_leaf_grids_remap; i ++)
			sprintf(new_grid_name, "%s_%s", new_grid_name, leaf_grids_remap[i]->get_coord_label());
		strcpy((*interchange_grid)->grid_name, new_grid_name);
    }
    
    EXECUTION_REPORT(REPORT_ERROR, -1, this->is_similar_grid_with(*interchange_grid), "remap software error2 when generating interchange grids\n");
}


void Remap_grid_class::interchange_grid_fields_for_remapping(Remap_grid_class *interchange_grid, 
                                                             Remap_grid_class *remap_grid, 
                                                             Remap_grid_data_class *runtime_mask_field)
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];
    Remap_grid_class *grid_of_setting_coord_values;
    int i, j;


    /* Interchange center and vertex values according to the super_grid_of_setting_coord_values of leaf grids */
    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++) {
        grid_of_setting_coord_values = leaf_grids[i]->get_super_grid_of_setting_coord_values();
        if (grid_of_setting_coord_values == NULL || (remap_grid != NULL && !leaf_grids[i]->is_subset_of_grid(remap_grid)))
            continue;
        for (j = 0; j < grid_of_setting_coord_values->grid_center_fields.size(); j ++)
            grid_of_setting_coord_values->grid_center_fields[j]->interchange_grid_data(interchange_grid);
        for (j = 0; j < grid_of_setting_coord_values->grid_vertex_fields.size(); j ++)
            grid_of_setting_coord_values->grid_vertex_fields[j]->interchange_grid_data(interchange_grid);
    }

    /* Interchange mask value according to masked sub grids */
    if (runtime_mask_field != NULL)
        runtime_mask_field->interchange_grid_data(interchange_grid);
}


/* Function duplicate_grid recursively duplicates a grid for runtime remap operator in preorder */
Remap_grid_class *Remap_grid_class::duplicate_grid(Remap_grid_class *top_grid)
{
    int i;
    Remap_grid_class *duplicated_grid;
    Remap_grid_class *similar_grid;


    EXECUTION_REPORT(REPORT_ERROR, -1, this->duplicated_grid == NULL && !this->is_partial_grid(), "remap software error1 in duplicate_grid\n");

    duplicated_grid = new Remap_grid_class;
    duplicated_grid->initialize_grid_class_data();
    duplicated_grid->original_grid = this;
    duplicated_grid->generated_from_duplication = true;
    this->duplicated_grid = duplicated_grid;
    for (i = 0; i < this->sub_grids.size(); i ++)
        duplicated_grid->sub_grids.push_back(this->sub_grids[i]->duplicate_grid(top_grid));
    
    strcpy(duplicated_grid->grid_name, this->grid_name);
    strcpy(duplicated_grid->coord_unit, this->coord_unit);
    strcpy(duplicated_grid->coord_label, this->coord_label);
    duplicated_grid->grid_size = this->grid_size;
    duplicated_grid->num_dimensions = this->num_dimensions;
    duplicated_grid->num_vertexes = this->num_vertexes;
    duplicated_grid->first_super_grid_of_enable_setting_coord_value = this->first_super_grid_of_enable_setting_coord_value;
    duplicated_grid->cyclic = this->cyclic;
    duplicated_grid->super_grid_of_setting_coord_values = this->super_grid_of_setting_coord_values;
    duplicated_grid->enable_to_set_coord_values = false;
    duplicated_grid->are_vertex_values_set_in_default = this->are_vertex_values_set_in_default;
    duplicated_grid->sigma_grid_scale_factor = this->sigma_grid_scale_factor;
    duplicated_grid->sigma_grid_top_value = this->sigma_grid_top_value;
    duplicated_grid->boundary_min_lon = this->boundary_min_lon;
    duplicated_grid->boundary_max_lon = this->boundary_max_lon;
    duplicated_grid->boundary_min_lat = this->boundary_min_lat;
    duplicated_grid->boundary_max_lat = this->boundary_max_lat;
    duplicated_grid->whole_grid = this->whole_grid;
	duplicated_grid->using_V3D_level_coord = this->using_V3D_level_coord;
    if (this->sigma_grid_sigma_value_field != NULL)
        duplicated_grid->sigma_grid_sigma_value_field = this->sigma_grid_sigma_value_field->duplicate_grid_data_field(duplicated_grid, 1, true, true);
    if (this->hybrid_grid_coefficient_field != NULL)
        duplicated_grid->hybrid_grid_coefficient_field = this->hybrid_grid_coefficient_field->duplicate_grid_data_field(duplicated_grid, 1, true, true);
    if (this->redundant_cell_mark_field != NULL) {
        duplicated_grid->redundant_cell_mark_field = this->redundant_cell_mark_field->duplicate_grid_data_field(this, 1, true, true);
        duplicated_grid->redundant_cell_mark = (bool *) duplicated_grid->redundant_cell_mark_field->grid_data_field->data_buf;
    }

    /* super_grid_of_setting_coord_values will be redirected when it is a subset of top grid */
    if (this->num_dimensions == 1 && this->super_grid_of_setting_coord_values != NULL && 
        this->super_grid_of_setting_coord_values->is_subset_of_grid(top_grid)) {
        similar_grid = NULL;
        for (i = 0; i < remap_grid_manager->remap_grids.size(); i ++) 
            if (remap_grid_manager->remap_grids[i]->is_similar_grid_with(this->super_grid_of_setting_coord_values) && remap_grid_manager->remap_grids[i]->duplicated_grid != NULL)  {
                similar_grid = remap_grid_manager->remap_grids[i];
                break;
            }
        EXECUTION_REPORT(REPORT_ERROR, -1, similar_grid != NULL, "remap software error2 in duplicate_grid\n");
        duplicated_grid->super_grid_of_setting_coord_values = similar_grid->duplicated_grid;
    }

    /* duplicate grid center and vertex data fields */
    for (int i = 0; i < grid_center_fields.size(); i ++)
        duplicated_grid->grid_center_fields.push_back(grid_center_fields[i]->duplicate_grid_data_field(duplicated_grid, 1, true, true));
    for (int i = 0; i < grid_vertex_fields.size(); i ++)
        duplicated_grid->grid_vertex_fields.push_back(grid_vertex_fields[i]->duplicate_grid_data_field(duplicated_grid, this->num_vertexes, true, true));
    
    this->duplicated_grid = NULL;

	if (get_is_sphere_grid())
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Allocate H2D grid \"%s\" with size %ld from the grid \"%s\"", duplicated_grid->grid_name, grid_size, grid_name);

    return duplicated_grid;
}


Remap_grid_class *Remap_grid_class::generate_remap_operator_runtime_grid(Remap_grid_class *remap_grid, 
                                                                        Remap_operator_basis *remap_operator, 
                                                                        Remap_grid_data_class *runtime_mask)

{
    long i;
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256], *super_grid, *runtime_remap_grid;
    Remap_grid_data_class *vertex_value_field, *duplicated_grid_data;
    Remap_data_field *mask_data_field;


    EXECUTION_REPORT(REPORT_ERROR, -1, remap_grid->is_subset_of_grid(this), "remap software error1 in generate_remap_operator_runtime_grid\n");

    runtime_remap_grid = remap_grid->duplicate_grid(remap_grid);
    
    if (!remap_operator->get_is_operator_regridding())
        return runtime_remap_grid;

    runtime_remap_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, runtime_remap_grid);
    for (i = 0; i < num_leaf_grids; i ++) {
        if (leaf_grids[i]->has_grid_coord_label(COORD_LABEL_LEV) && (leaf_grids[i]->get_sigma_grid_sigma_value_field() != NULL || leaf_grids[i]->does_use_V3D_level_coord())) {
			if (leaf_grids[i]->grid_center_fields.size() > 0)
				continue;
			leaf_grids[i]->allocate_default_center_field();				
            leaf_grids[i]->num_vertexes = 2;
			for (int j = 0; j < leaf_grids[i]->grid_size; j ++)
				 ((double*)leaf_grids[i]->grid_center_fields[0]->get_grid_data_field()->data_buf)[j] = j % 2 + 1;
            leaf_grids[i]->grid_vertex_fields.push_back(leaf_grids[i]->grid_center_fields[0]->duplicate_grid_data_field(remap_grid, leaf_grids[i]->num_vertexes, false, false));
        }
        else {
            super_grid = leaf_grids[i]->super_grid_of_setting_coord_values;
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, super_grid != NULL && super_grid->is_subset_of_grid(runtime_remap_grid), "remap software error2 in generate_remap_operator_runtime_grid \"%s\"   \"%s\" %lx  \"%s\"\n", grid_name, leaf_grids[i]->get_grid_name(), leaf_grids[i], super_grid->get_grid_name());
        }
    }

    if (runtime_mask != NULL) {
        if (runtime_mask->get_coord_value_grid()->is_similar_grid_with(remap_grid)) {
            runtime_remap_grid->grid_mask_field = runtime_mask->duplicate_grid_data_field(remap_grid, 1, true, true);
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "The runtime remap grid \"%s\" directly uses the runtime_mask when remapping the field with grid \"%s\"", runtime_remap_grid->get_grid_name(), grid_name);
        }
        else if (runtime_mask->get_coord_value_grid()->have_overlap_with_grid(remap_grid)) {
			EXECUTION_REPORT(REPORT_ERROR, -1, remap_grid->is_subset_of_grid(runtime_mask->get_coord_value_grid()), "software error in generate_remap_operator_runtime_grid");
            runtime_remap_grid->original_grid_mask_field = runtime_mask;
            runtime_mask->interchange_grid_data(this);
            runtime_remap_grid->grid_mask_field = runtime_mask->duplicate_grid_data_field(remap_grid, 1, false, false);
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "The runtime remap grid \"%s\" uses a super runtime_mask when remapping the field with grid \"%s\"", runtime_remap_grid->get_grid_name(), grid_name);
        }
		else {			
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "The runtime remap grid \"%s\" uses a outer runtime_mask when remapping the field with grid \"%s\"", runtime_remap_grid->get_grid_name(), grid_name);
			runtime_mask->interchange_grid_data(this);
		}
    }

    return runtime_remap_grid;
}


void Remap_grid_class::allocate_default_center_field()
{
	Remap_data_field *remap_data_field;

	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_dimensions == 1 && grid_center_fields.size() == 0, "Software error in Remap_grid_class::allocate_default_center_field");

	remap_data_field = new Remap_data_field;	
	strcpy(remap_data_field->field_name_in_application, coord_label);
	strcpy(remap_data_field->field_name_in_IO_file, coord_label);
	strcpy(remap_data_field->field_name_in_IO_file, coord_label);
	strcpy(remap_data_field->data_type_in_application, DATA_TYPE_DOUBLE);
	strcpy(remap_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE);
	remap_data_field->required_data_size = remap_data_field->read_data_size = grid_size;
	remap_data_field->data_buf = new double [grid_size];
	remap_data_field->set_field_unit(coord_unit);
	grid_center_fields.push_back(new Remap_grid_data_class(this, remap_data_field));
}


void Remap_grid_class::read_grid_data_through_span(char extension_names[16][256], 
                                                  const char *bound_start, 
                                                  const char *bound_end, 
                                                  long num_elements,
                                                  const char *span_function)
{
    Remap_field_attribute grid_field_attribute;    


    EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 1 && grid_size > 0, "only the coordinate center values of 1D grid whose size is known can be generated through spanning\n");
    check_center_coord_value_can_be_set(this->coord_label);

    if (words_are_the_same(span_function, FUNCTION_WORD_ISPAN))
        grid_center_fields.push_back(new Remap_grid_data_class(this->coord_label,
                                                               this->grid_name,
                                                               bound_start,
                                                               bound_end,
                                                               num_elements,
                                                               DATA_TYPE_LONG,
                                                               DATA_TYPE_DOUBLE));
    else if (words_are_the_same(span_function, FUNCTION_WORD_FSPAN))
        grid_center_fields.push_back(new Remap_grid_data_class(this->coord_label,
                                                               this->grid_name,
                                                               bound_start,
                                                               bound_end,
                                                               num_elements,
                                                               DATA_TYPE_DOUBLE,
                                                               DATA_TYPE_DOUBLE));
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in read_grid_data_through_span\n");
    if (words_are_the_same(this->coord_label, COORD_LABEL_LON)) {
        grid_center_fields[0]->get_grid_data_field()->set_field_long_name("longitude");
        grid_center_fields[0]->get_grid_data_field()->set_field_unit(COORD_UNIT_DEGREES);
    }
    if (words_are_the_same(this->coord_label, COORD_LABEL_LAT)) {
        grid_center_fields[0]->get_grid_data_field()->set_field_long_name("latitude");
        grid_center_fields[0]->get_grid_data_field()->set_field_unit(COORD_UNIT_DEGREES);
    }
}


void Remap_grid_class::generate_default_grid_mask()
{
    Remap_data_field *mask_data_field;


    mask_data_field = new Remap_data_field;
    strcpy(mask_data_field->field_name_in_application, GRID_MASK_LABEL);
    strcpy(mask_data_field->data_type_in_application, DATA_TYPE_BOOL);
    strcpy(mask_data_field->data_type_in_IO_file, DATA_TYPE_BOOL);
    strcpy(mask_data_field->field_name_in_IO_file, GRID_MASK_LABEL);
    mask_data_field->required_data_size = grid_size;    
    mask_data_field->read_data_size = grid_size;
    mask_data_field->data_buf = new char [mask_data_field->required_data_size*get_data_type_size(mask_data_field->data_type_in_application)];
    for (int i = 0; i < grid_size; i ++)
        ((bool*)(mask_data_field->data_buf))[i] = true;

    grid_mask_field = new Remap_grid_data_class(this, mask_data_field);
}


void Remap_grid_class::extract_mask(const char *field_data_name, const char *fill_value_lower_bound_str, const char *fill_value_higher_bound_str)
{
    Remap_grid_data_class *field_data;
    double fill_value_lower_bound, fill_value_higher_bound;
    double *field_data_values;
    bool *mask_values;
    long i;


    check_mask_value_can_be_set();
    field_data = remap_field_data_manager->search_remap_field_data(field_data_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, field_data->have_data_content(), "the data values of field %s have not been set\n", field_data->get_grid_data_field()->field_name_in_application);
    EXECUTION_REPORT(REPORT_ERROR, -1, this->is_subset_of_grid(field_data->get_coord_value_grid()), 
                 "grid %s should be a subset of grid %s which is the associative grid of field data \"%s\"\n",
                 this->grid_name, field_data->get_coord_value_grid()->get_grid_name(), field_data->get_grid_data_field()->field_name_in_application);
    sscanf(fill_value_lower_bound_str, "%lf", &fill_value_lower_bound);
    sscanf(fill_value_higher_bound_str, "%lf", &fill_value_higher_bound);
    EXECUTION_REPORT(REPORT_ERROR, -1, fill_value_lower_bound <= fill_value_higher_bound, "the last input parameter must be no smaller than the second input parameter\n");

    generate_default_grid_mask();
    
    field_data->interchange_grid_data(this);
    field_data_values = (double*) field_data->get_grid_data_field()->data_buf;
    mask_values = (bool*) grid_mask_field->get_grid_data_field()->data_buf;
    for (i = 0; i < grid_size; i ++)
        mask_values[i] = !(field_data_values[i] >= fill_value_lower_bound && field_data_values[i] <= fill_value_higher_bound);
}


void Remap_grid_class::compute_ocn_mask(const char *topo_field_name, double unit_trans)
{
    Remap_grid_class *topo_grid, *lonlat_sub_grid_this, *lev_sub_grid;
    Remap_grid_data_class *topo_field, *lon_field_topo_grid, *lat_field_topo_grid;
    double *lev_vertex_values;
    double *lon_values_topo_grid, *lat_values_topo_grid, *topo_values, topo_coords[2];
    Remap_grid_class *leaf_grids[256], *sized_sub_grids[256], *sub_grids[256], *remap_grids[256];
    Remap_operator_regrid *remap_operator;
    Remap_operator_grid *remap_operator_grid_this, *current_operator_grid;
    double *ocn_counts, *total_counts;
    int i, num_leaf_grids, num_sized_sub_grids;
    long j, k, lonlat_src_cell_index, tmp_grid_size;
    bool *mask_values;
    double last_time, current_time;
    char tmp_str[256];


    /* Check for this grid */
    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_dimensions == 2 || this->num_dimensions == 3,
                 "grid %s which is used to set ocn mask should be a 2D or 3D grid\n", this->grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, this->has_grid_coord_label(COORD_LABEL_LAT) && this->has_grid_coord_label(COORD_LABEL_LON),
                 "grid %s must have longitude and latitude coordinates", this->grid_name);
    if (this->num_dimensions == 3)
        EXECUTION_REPORT(REPORT_ERROR, -1, this->has_grid_coord_label(COORD_LABEL_LEV), 
                     "grid %s is a 3D grid which must have coordinate level\n", this->grid_name);
    this->check_mask_value_can_be_set();
    this->end_grid_definition_stage(NULL);
    this->get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->get_grid_center_field() != NULL,
                     "the center value of coordinate %s of %s has not been set\n",
                     leaf_grids[i]->grid_name, this->grid_name);
        if (words_are_the_same(leaf_grids[i]->coord_label, COORD_LABEL_LEV))
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->super_grid_of_setting_coord_values == leaf_grids[i], 
                         "Can not generate ocean mask of grid %s because it is a SIGMA grid.\n",
                         this->grid_name);
    }
    for (i = 0; i < num_leaf_grids; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->get_grid_vertex_field() != NULL,
                     "the vertex value of coordinate %s of %s has not been set by user or in default\n",
                     leaf_grids[i]->grid_name, this->grid_name);

    /* Check for the topo field as well its grid */
    topo_field = remap_field_data_manager->search_remap_field_data(topo_field_name);
    topo_grid = topo_field->get_coord_value_grid();
    EXECUTION_REPORT(REPORT_ERROR, -1, topo_grid != NULL, "the topo field is not a gridded data\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, topo_grid->get_is_sphere_grid(),
                 "the associative grid of field %s is %s, it must be a 2D sphere grid with longitude and latitude coordinates\n",
                 topo_field_name, topo_grid->grid_name);
    topo_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->get_grid_center_field() != NULL,
                     "the center value of coordinate %s of %s has not been set\n",
                     leaf_grids[i]->grid_name, topo_grid->grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, topo_field->have_data_content(), "the data values of %s have not been set\n", topo_field_name);

    this->generate_default_grid_mask();
    
    /* Generate interchange grid for this grid and topo grid */
    this->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
    tmp_grid_size = 1;
    lev_sub_grid = NULL;
    for (i = 0, j = 0; i < num_sized_sub_grids; i ++)
        if (sized_sub_grids[i]->has_grid_coord_label(COORD_LABEL_LON) || sized_sub_grids[i]->has_grid_coord_label(COORD_LABEL_LAT)) {
            sub_grids[j++] = sized_sub_grids[i];
            tmp_grid_size *= sized_sub_grids[i]->grid_size;
        }
        else lev_sub_grid = sized_sub_grids[i];
    lonlat_sub_grid_this = new Remap_grid_class("TMP_GRID", j, sub_grids, tmp_grid_size);
    this->interchange_grid_fields_for_remapping(lonlat_sub_grid_this, NULL, NULL);
    topo_grid->interchange_grid_fields_for_remapping(topo_grid, NULL, NULL);
    topo_field->interchange_grid_data(topo_grid);

    /* Get level coordinate values */
    lev_vertex_values = NULL;
    if (this->num_dimensions == 3) {
        lev_vertex_values = (double*) lev_sub_grid->get_grid_vertex_field()->get_grid_data_field()->data_buf;
        EXECUTION_REPORT(REPORT_ERROR, -1, lev_sub_grid->num_vertexes == 2, "remap software error3 in compute_ocn_mask\n");
        for (i = 1; i < lev_sub_grid->grid_size*2; i ++)
            if (fabs(lev_vertex_values[0]) < fabs(lev_vertex_values[lev_sub_grid->grid_size*2-1]))
                EXECUTION_REPORT(REPORT_ERROR, -1, fabs(lev_vertex_values[i]) >= fabs(lev_vertex_values[i-1]),
                             "the coordinate values of ocn level grid %s are not sorted\n",
                             lev_sub_grid->grid_name);
            else EXECUTION_REPORT(REPORT_ERROR, -1, fabs(lev_vertex_values[i]) <= fabs(lev_vertex_values[i-1]),
                              "the coordinate values of ocn level grid %s are not sorted\n",
                              lev_sub_grid->grid_name);
    }

    if (execution_phase_number == 0)
        return;

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "generate ocn mask of grid %s from terrain data", this->grid_name);

    /* Prepare remap operator grid and lon/lat values of topo grid, etc */
    remap_grids[0] = lonlat_sub_grid_this;
    remap_grids[1] = topo_grid;
    remap_operator = new Remap_operator_regrid("TMP_remap_operator", 2, remap_grids, 2);
    remap_operator_grid_this = new Remap_operator_grid(lonlat_sub_grid_this, remap_operator, true, false);
    topo_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, topo_grid);
    remap_operator_grid_this->update_operator_grid_data();
    lon_field_topo_grid = topo_grid->expand_to_generate_full_coord_value(leaf_grids[0]->get_grid_center_field());
    lat_field_topo_grid = topo_grid->expand_to_generate_full_coord_value(leaf_grids[1]->get_grid_center_field());
    lon_values_topo_grid = (double*) lon_field_topo_grid->get_grid_data_field()->data_buf;
    lat_values_topo_grid = (double*) lat_field_topo_grid->get_grid_data_field()->data_buf;
    topo_values = (double*) topo_field->get_grid_data_field()->data_buf;
    ocn_counts = new double [this->grid_size];
    total_counts = new double [this->grid_size];
    for (j = 0; j < this->grid_size; j ++) {
        ocn_counts[j] = 0;
        total_counts[j] = 0;
    }

    for (j = 0; j < topo_grid->grid_size; j ++) {
        topo_coords[0] = lon_values_topo_grid[j];
        topo_coords[1] = lat_values_topo_grid[j];
        current_operator_grid = remap_operator_grid_this;
        if (fabs(lat_values_topo_grid[j]) > SPHERE_GRID_ROTATION_LAT_THRESHOLD) {
            rotate_sphere_coordinate(lon_values_topo_grid[j], lat_values_topo_grid[j], topo_coords[0], topo_coords[1]);
            current_operator_grid = remap_operator_grid_this->get_rotated_remap_operator_grid();
        }
        lonlat_src_cell_index = current_operator_grid->search_cell_of_locating_point(topo_coords, true);
        if (lonlat_src_cell_index == -1)
            lonlat_src_cell_index = current_operator_grid->search_cell_of_locating_point(topo_coords, false);
        if (lonlat_src_cell_index == -1)
            continue;
        if (this->num_dimensions == 2) {
            if (topo_values[j] <= 0)
                ocn_counts[lonlat_src_cell_index] ++;
            total_counts[lonlat_src_cell_index] ++;
        }
        else {
            for (k = 0; k < lev_sub_grid->grid_size; k ++)
                total_counts[lonlat_src_cell_index*lev_sub_grid->grid_size+k] ++;
            if (topo_values[j] <= 0)
                for (k = 0; k < lev_sub_grid->grid_size; k ++)
                    if (fabs(topo_values[j]) >= fabs(lev_vertex_values[k*2]) || fabs(topo_values[j]) >= fabs(lev_vertex_values[k*2+1]))
                        ocn_counts[lonlat_src_cell_index*lev_sub_grid->grid_size+k] ++;
        }
    }

    mask_values = (bool*) grid_mask_field->get_grid_data_field()->data_buf;
    if (this->num_dimensions == 2)
        this->grid_mask_field->interchange_grid_data(lonlat_sub_grid_this);
    else this->grid_mask_field->interchange_grid_data(lev_sub_grid);
    for (j = 0; j < this->grid_size; j ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, total_counts[j] > 0, "the terrain grid is too coarse\n");
        if (ocn_counts[j]/total_counts[j] >= 0.5)
            mask_values[j] = true;
        else mask_values[j] = false;
    }

    this->grid_mask_field->interchange_grid_data(this);
    
    delete [] ocn_counts;
    delete [] total_counts;
    delete remap_operator_grid_this;
    delete lon_values_topo_grid;
    delete lat_values_topo_grid;
    delete remap_operator;
    delete lonlat_sub_grid_this;
}


void Remap_grid_class::gen_lev_coord_from_sigma_or_hybrid(char extension_names[16][256], 
                                                        const char *lev_coord_bot_str, 
                                                        const char *lev_coord_top_str,
                                                        const char *lev_sigma_coord_str, 
                                                        const char *hybrid_grid_coefficient_str,
                                                        double scale_factor)
{
    Remap_grid_data_class *field_lev_coord_bot;
    double value_lev_coord_top;
    Remap_grid_data_class *lev_sigma_coord;
    Remap_grid_class *similar_grid;
    Remap_grid_class *leaf_grids[256], *lev_leaf_grid, *sized_grids[256];
    long i, j;
    int num_leaf_grids, num_sized_grids;
    Remap_grid_data_class *remap_grid_data_field, *hybrid_grid_coefficient = NULL;
    double *data_grid_field, *data_array_bot, *data_array_ptr, *data_sigma, data_bot, data_top, *lev_sigma_coord_values;
    Remap_grid_class *interchanged_grid;
    double full_ratio;
    bool *horizontal_grid_mask_values, ascending_order;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_dimensions == 3 && this->has_grid_coord_label(COORD_LABEL_LON) &&
                 this->has_grid_coord_label(COORD_LABEL_LAT) && this->has_grid_coord_label(COORD_LABEL_LEV), 
                 "grid %s must be a 3D spacial grid with lon, lat and lev coordinate\n", this->get_grid_name());
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(extension_names[1], COORD_LABEL_LEV), "this function can only generate the center value of level coordinate values\n");
    check_center_coord_value_can_be_set(COORD_LABEL_LEV);
    this->get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->coord_label, COORD_LABEL_LEV)) {
            lev_leaf_grid = leaf_grids[i];
            break;
        }
    EXECUTION_REPORT(REPORT_ERROR, -1, lev_leaf_grid->grid_size > 0, "remap software error0 in en_lev_coord_from_sigma\n");
    field_lev_coord_bot = remap_field_data_manager->search_remap_field_data(lev_coord_bot_str);
    lev_sigma_coord = remap_field_data_manager->search_remap_field_data(lev_sigma_coord_str);
    EXECUTION_REPORT(REPORT_ERROR, -1, field_lev_coord_bot->get_coord_value_grid() != NULL, "%s must be a field data with associative grid\n", lev_coord_bot_str);
    EXECUTION_REPORT(REPORT_ERROR, -1, field_lev_coord_bot->get_coord_value_grid()->is_subset_of_grid(this) && field_lev_coord_bot->get_coord_value_grid()->get_is_sphere_grid(),
                 "the grid of field \"%s\" must be a subgrid of \"%s\" only without coordinate level\n",
                 lev_coord_bot_str, this->grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, lev_sigma_coord->get_grid_data_field()->required_data_size == lev_leaf_grid->grid_size,
                 "the size of field %s must be the same with 1D level subgrid of grid %s\n", lev_sigma_coord_str, this->grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, lev_leaf_grid->grid_size >= 3, "the level subgrid %s of grid %s must have at least 3 levels. Please check.", lev_leaf_grid->get_grid_name(), this->get_grid_name());
    EXECUTION_REPORT(REPORT_ERROR, -1, field_lev_coord_bot->have_data_content(), "the data value of field %s has not been set\n", lev_coord_bot_str);

    lev_sigma_coord_values = (double*)lev_sigma_coord->get_grid_data_field()->data_buf;
    for (i = 1; i < lev_leaf_grid->grid_size; i ++) {
        if (lev_sigma_coord_values[i-1] < lev_sigma_coord_values[i]) {
            ascending_order = true;
            break;
        }
        if (lev_sigma_coord_values[i-1] > lev_sigma_coord_values[i]) {
            ascending_order = false;
            break;
        }
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, i < lev_leaf_grid->grid_size, "All sigma values in field %s for grid %s are the same. Please check.", lev_sigma_coord_str, grid_name);
    for (i = 1; i < lev_leaf_grid->grid_size; i ++)
        if (ascending_order)
            EXECUTION_REPORT(REPORT_ERROR, -1, (ascending_order && (lev_sigma_coord_values[i-1] <= lev_sigma_coord_values[i])) || (!ascending_order && (lev_sigma_coord_values[i-1] >= lev_sigma_coord_values[i])), 
            "sigma values in field %s for grid %s must be in ascending order or descending order. Please check.", lev_sigma_coord_str, grid_name);

    if (hybrid_grid_coefficient_str != NULL) {
        hybrid_grid_coefficient = remap_field_data_manager->search_remap_field_data(hybrid_grid_coefficient_str);
        EXECUTION_REPORT(REPORT_ERROR, -1, hybrid_grid_coefficient->get_grid_data_field()->required_data_size == lev_leaf_grid->get_grid_size(), 
                         "the size of field %s must be the same with 1D level subgrid of grid %s\n", hybrid_grid_coefficient_str, this->grid_name);
    }
    sscanf(lev_coord_top_str, "%lf", &data_top);
    this->level_V3D_coord_trigger_field_specified = true;
    allocate_sigma_grid_specific_fields(hybrid_grid_coefficient, lev_sigma_coord, field_lev_coord_bot, data_top, scale_factor);
    calculate_lev_sigma_values();
}


void Remap_grid_class::set_lev_grid_sigma_info(double top_value, const double *sigma_values, const double *hybrid_grid_coefficients, double scale_factor)
{
    Remap_grid_data_class *sigma_value_field = NULL, *hybrid_grid_coefficient_field = NULL;
    Remap_data_field *remap_data_field;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 1 && words_are_the_same(coord_label, COORD_LABEL_LEV), "Software error in Remap_grid_class::set_lev_grid_sigma_info: only lev grid can be set sigma information", grid_name);
    remap_data_field = new Remap_data_field;    
    strcpy(remap_data_field->field_name_in_application, "lev_sigma");
    strcpy(remap_data_field->data_type_in_application, DATA_TYPE_DOUBLE);
    strcpy(remap_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE);
    remap_data_field->required_data_size = remap_data_field->read_data_size = grid_size;
    remap_data_field->data_buf = new double [grid_size];
    memcpy(remap_data_field->data_buf, sigma_values, grid_size*sizeof(double));
    sigma_value_field = new Remap_grid_data_class(this, remap_data_field);
    if (hybrid_grid_coefficients != NULL) {
        strcpy(remap_data_field->field_name_in_application, "lev_hybrid_coefB");
        remap_data_field = new Remap_data_field;    
        strcpy(remap_data_field->field_name_in_application, "lev_hybrid_coefA");
        strcpy(remap_data_field->data_type_in_application, DATA_TYPE_DOUBLE);
        strcpy(remap_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE);
        remap_data_field->required_data_size = remap_data_field->read_data_size = grid_size;
        remap_data_field->data_buf = new double [grid_size];
        memcpy(remap_data_field->data_buf, hybrid_grid_coefficients, grid_size*sizeof(double));
        hybrid_grid_coefficient_field = new Remap_grid_data_class(this, remap_data_field);
    }
    allocate_sigma_grid_specific_fields(hybrid_grid_coefficient_field, sigma_value_field, NULL, top_value, scale_factor);
}


void Remap_grid_class::set_lev_grid_sigma_info(const char *sigma_value_field_name, double top_value, double scale_factor, const char *hybrid_grid_coefficient_field_name)
{
    Remap_grid_data_class *sigma_value_field, *hybrid_grid_coefficient_field;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 1 && words_are_the_same(coord_label, COORD_LABEL_LEV), 
                     "Grid %s must be a 1D grid at vertical direction (the label of grid is \"lev\": only lev grid can be set sigma information", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_coord_values, "Grid %s has been used (for example, interpolation) before. It cannot be set with sigma information", grid_name);
    sigma_value_field = remap_field_data_manager->search_remap_field_data(sigma_value_field_name);
    if (hybrid_grid_coefficient_field_name != NULL) {
        hybrid_grid_coefficient_field = remap_field_data_manager->search_remap_field_data(hybrid_grid_coefficient_field_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, hybrid_grid_coefficient_field != NULL, "C-Coupler error1 in set_lev_grid_sigma_info");
    }
    else hybrid_grid_coefficient_field = NULL;
    EXECUTION_REPORT(REPORT_ERROR, -1, sigma_value_field != NULL, "C-Coupler error2 in set_lev_grid_sigma_info");
    allocate_sigma_grid_specific_fields(hybrid_grid_coefficient_field, sigma_value_field, NULL, top_value, scale_factor);
}


bool Remap_grid_class::does_use_V3D_level_coord()
{
	if (!has_grid_coord_label(COORD_LABEL_LEV))
		return false;
	
	return get_a_leaf_grid(COORD_LABEL_LEV)->using_V3D_level_coord;
}


void Remap_grid_class::set_using_V3D_level_coord()
{
	Remap_grid_class *leaf_grid = get_a_leaf_grid(COORD_LABEL_LEV);	
	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, leaf_grid != NULL, "Software error in Remap_grid_class::set_using_V3D_level_coord");

	leaf_grid->using_V3D_level_coord = true;
}


Remap_grid_class *Remap_grid_class::get_a_leaf_grid_of_sigma_or_hybrid()
{
    Remap_grid_class *leaf_grid;

    
    leaf_grid = get_a_leaf_grid(COORD_LABEL_LEV);
    if (leaf_grid->sigma_grid_sigma_value_field != NULL)
        return leaf_grid;

    return NULL;
}


bool Remap_grid_class::is_sigma_grid()
{
    Remap_grid_class *leaf_grid;

    
    if (num_dimensions != 3)
        return false;

    if (!has_grid_coord_label(COORD_LABEL_LON) || !has_grid_coord_label(COORD_LABEL_LAT) || !has_grid_coord_label(COORD_LABEL_LEV))
        return false;

    if (get_a_leaf_grid_of_sigma_or_hybrid() != NULL)
        return true;

    return false;
}


bool Remap_grid_class::is_level_V3D_coord_trigger_field_updated()
{
    bool result = false;
    

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (this->is_sigma_grid() || this->does_use_V3D_level_coord()) && level_V3D_coord_dynamic_trigger_field != NULL && level_V3D_coord_trigger_field != NULL, "C-Coupler error1 in Remap_grid_class::is_level_V3D_coord_trigger_field_updated: \"%s\" at %lx", this->grid_name, (char*)this); 
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->required_data_size == level_V3D_coord_trigger_field->get_grid_data_field()->required_data_size, "C-Coupler error1 in Remap_grid_class::is_level_V3D_coord_trigger_field_updated");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(level_V3D_coord_trigger_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "C-Coupler error in Remap_grid_class::is_level_V3D_coord_trigger_field_updated: wrong data type");
    
    if (words_are_the_same(level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT)) {
        for (int i = 0; i < level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->required_data_size; i ++) {
            if (!level_V3D_coord_trigger_field_specified || ((float*) level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->data_buf)[i] != ((double*) level_V3D_coord_trigger_field->get_grid_data_field()->data_buf)[i])
                result = true;
            ((double*) level_V3D_coord_trigger_field->get_grid_data_field()->data_buf)[i] = ((float*) level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->data_buf)[i];
        }
    }
    else {
        for (int i = 0; i < level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->required_data_size; i ++) {
            if (!level_V3D_coord_trigger_field_specified || ((double*) level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->data_buf)[i] != ((double*) level_V3D_coord_trigger_field->get_grid_data_field()->data_buf)[i])
                result = true;
            ((double*) level_V3D_coord_trigger_field->get_grid_data_field()->data_buf)[i] = ((double*) level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->data_buf)[i];
        }
    }

    level_V3D_coord_trigger_field_specified = true;

    if (result)
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "surface field for sigma grid %s has been updated", grid_name);
    
    return result;
}


void Remap_grid_class::set_level_V3D_coord_dynamic_trigger_field(Remap_grid_data_class *value_field)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, this->is_sigma_grid() || this->does_use_V3D_level_coord(), "C-Coupler error1 in Remap_grid_class::set_level_V3D_coord_dynamic_trigger_field");
    EXECUTION_REPORT(REPORT_ERROR, -1, value_field->get_coord_value_grid()->is_subset_of_grid(this), "C-Coupler error2 in Remap_grid_class::set_level_V3D_coord_dynamic_trigger_field");
    EXECUTION_REPORT(REPORT_ERROR, -1, level_V3D_coord_dynamic_trigger_field == NULL, "The surface field of grid %s has been specified by the models before. Please check.", grid_name);

    level_V3D_coord_dynamic_trigger_field = value_field;
    if (level_V3D_coord_trigger_field == NULL)
        level_V3D_coord_trigger_field = value_field->duplicate_grid_data_field(value_field->get_coord_value_grid(), 1, false, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, level_V3D_coord_trigger_field->get_grid_data_field()->required_data_size == level_V3D_coord_dynamic_trigger_field->get_grid_data_field()->required_data_size, "C-Coupler error in Remap_grid_class::set_level_V3D_coord_dynamic_trigger_field: wrong field size");
    if (words_are_the_same(level_V3D_coord_trigger_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT)) {
        strcpy(level_V3D_coord_trigger_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE);
        strcpy(level_V3D_coord_trigger_field->get_grid_data_field()->data_type_in_IO_file, DATA_TYPE_DOUBLE);
        delete [] level_V3D_coord_trigger_field->get_grid_data_field()->data_buf;
        level_V3D_coord_trigger_field->get_grid_data_field()->data_buf = new double [level_V3D_coord_trigger_field->get_grid_data_field()->required_data_size];
    }
    if (grid_center_fields.size() == 0) {
        grid_center_fields.push_back(level_V3D_coord_trigger_field->duplicate_grid_data_field(this, 1, false, false));
        strcpy(grid_center_fields[0]->get_grid_data_field()->field_name_in_application, COORD_LABEL_LEV);
        strcpy(grid_center_fields[0]->get_grid_data_field()->field_name_in_IO_file, COORD_LABEL_LEV);
    }
}


void Remap_grid_class::allocate_sigma_grid_specific_fields(Remap_grid_data_class *hybrid_grid_coefficient_field, Remap_grid_data_class *sigma_grid_sigma_value_field, Remap_grid_data_class *level_V3D_coord_trigger_field, double data_top, double scale_factor)
{
    int num_sized_sub_grids, i, j;
    Remap_grid_class *sized_sub_grids[256];
    Remap_grid_class *sphere_grid, *existing_grid;


    if (this->level_V3D_coord_trigger_field != NULL)
        return;

    if (sigma_grid_sigma_value_field != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, enable_to_set_coord_values && this->get_a_leaf_grid(COORD_LABEL_LEV)->enable_to_set_coord_values, "Grid %s cannot be set to a sigma grid. It may be because the grid has been used (for example, interpolation) before or its level subgrid is a level middle grid.", grid_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, this->get_a_leaf_grid(COORD_LABEL_LEV)->hybrid_grid_coefficient_field == NULL, "Hybrid coefficient of 1D grid %s has been set before", this->get_a_leaf_grid(COORD_LABEL_LEV)->get_grid_name());
        EXECUTION_REPORT(REPORT_ERROR, -1, this->get_a_leaf_grid(COORD_LABEL_LEV)->sigma_grid_sigma_value_field == NULL, "Sigma information of 1D grid %s has been set before", this->get_a_leaf_grid(COORD_LABEL_LEV)->get_grid_name());
        EXECUTION_REPORT(REPORT_ERROR, -1, this->get_a_leaf_grid(COORD_LABEL_LEV)->grid_center_fields.size() == 0, "1D grid %s has been set to non-sigma grid before. Therefore it cannot be set sigma information", this->get_a_leaf_grid(COORD_LABEL_LEV)->get_grid_name());
        EXECUTION_REPORT(REPORT_ERROR, -1, this->grid_center_fields.size() == 0, "C-Coupler error5 in allocate_sigma_grid_specific_fields");
        EXECUTION_REPORT(REPORT_ERROR, -1, sigma_grid_sigma_value_field->get_grid_data_field()->required_data_size == this->get_a_leaf_grid(COORD_LABEL_LEV)->get_grid_size(), 
                         "The size of field %s is different from the size of grid %s. Field %s cannot be used for sigma or hybrid grid %s",
                         sigma_grid_sigma_value_field->get_grid_data_field()->field_name_in_application, this->get_a_leaf_grid(COORD_LABEL_LEV)->get_grid_name(), sigma_grid_sigma_value_field->get_grid_data_field()->field_name_in_application);
        this->get_a_leaf_grid(COORD_LABEL_LEV)->sigma_grid_sigma_value_field = sigma_grid_sigma_value_field->duplicate_grid_data_field(get_a_leaf_grid(COORD_LABEL_LEV), 1, true, true);
        if (hybrid_grid_coefficient_field != NULL) {
            EXECUTION_REPORT(REPORT_ERROR, -1, hybrid_grid_coefficient_field->get_grid_data_field()->required_data_size == this->get_a_leaf_grid(COORD_LABEL_LEV)->get_grid_size(), 
                             "The size of field %s is different from the size of grid %s. Field %s cannot be used for sigma or hybrid grid %s",
                             hybrid_grid_coefficient_field->get_grid_data_field()->field_name_in_application, this->get_a_leaf_grid(COORD_LABEL_LEV)->get_grid_name(), hybrid_grid_coefficient_field->get_grid_data_field()->field_name_in_application);
            this->get_a_leaf_grid(COORD_LABEL_LEV)->hybrid_grid_coefficient_field = hybrid_grid_coefficient_field->duplicate_grid_data_field(get_a_leaf_grid(COORD_LABEL_LEV), 1, true, true);
        }
        if (level_V3D_coord_trigger_field != NULL) {
            this->level_V3D_coord_trigger_field = level_V3D_coord_trigger_field->duplicate_grid_data_field(level_V3D_coord_trigger_field->get_coord_value_grid(), 1, true, true);
            grid_center_fields.push_back(sigma_grid_sigma_value_field->duplicate_grid_data_field(get_a_leaf_grid(COORD_LABEL_LEV), 1, false, false));  // temp level coord value field for global grid that should not be further used
            strcpy(grid_center_fields[0]->get_grid_data_field()->field_name_in_application, COORD_LABEL_LEV);
            strcpy(grid_center_fields[0]->get_grid_data_field()->field_name_in_IO_file, COORD_LABEL_LEV);
        }
        this->get_a_leaf_grid(COORD_LABEL_LEV)->sigma_grid_top_value = data_top;
        this->get_a_leaf_grid(COORD_LABEL_LEV)->sigma_grid_scale_factor = scale_factor;        
        this->get_a_leaf_grid(COORD_LABEL_LEV)->end_grid_definition_stage(NULL);
        return;
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Allocate sigma grid specific fields for grid %s", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, is_sigma_grid(), "C-Coupler error0 in allocate_sigma_grid_specific_fields");
        
    get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
    for (i = 0, j = 0; i < num_sized_sub_grids; i ++)
        if (!sized_sub_grids[i]->has_grid_coord_label(COORD_LABEL_LEV))
            sized_sub_grids[j++] = sized_sub_grids[i];
    sphere_grid = new Remap_grid_class("sphere_grid_for_level_V3D_coord_trigger_field", j, sized_sub_grids, 0);
    existing_grid = remap_grid_manager->search_same_remap_grid(sphere_grid);
    if (existing_grid != NULL) {
        delete sphere_grid;
        sphere_grid = existing_grid;
    }
    else remap_grid_manager->add_temp_grid(sphere_grid);

    this->level_V3D_coord_trigger_field = this->get_a_leaf_grid_of_sigma_or_hybrid()->sigma_grid_sigma_value_field->duplicate_grid_data_field(sphere_grid, 1, false, false);
    strcpy(this->level_V3D_coord_trigger_field->get_grid_data_field()->field_name_in_application, COORD_LABEL_LEV);
    strcpy(this->level_V3D_coord_trigger_field->get_grid_data_field()->field_name_in_IO_file, COORD_LABEL_LEV);

    if (this->grid_center_fields.size() != 0)
        EXECUTION_REPORT(REPORT_ERROR, -1, this->grid_center_fields.size() == 1, "C-Coupler error5 in allocate_sigma_grid_specific_fields");
    else {
        grid_center_fields.push_back(this->get_a_leaf_grid_of_sigma_or_hybrid()->get_sigma_grid_sigma_value_field()->duplicate_grid_data_field(get_a_leaf_grid(COORD_LABEL_LEV), 1, false, false));  // temp level coord value field for global grid that should not be further used
        strcpy(grid_center_fields[0]->get_grid_data_field()->field_name_in_application, COORD_LABEL_LEV);
        strcpy(grid_center_fields[0]->get_grid_data_field()->field_name_in_IO_file, COORD_LABEL_LEV);
    }
}


Remap_grid_data_class *Remap_grid_class::get_sigma_grid_sigma_value_field()
{
    return sigma_grid_sigma_value_field;
}


void Remap_grid_class::update_grid_center_3D_level_field_from_external()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->num_dimensions == 3 && this->does_use_V3D_level_coord() && this->grid_center_fields.size() == 1, "Software error in Remap_grid_class::update_grid_center_3D_level_field_from_external");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, get_data_type_size(level_V3D_coord_trigger_field->get_grid_data_field()->data_type_in_application) == 8 && get_data_type_size(grid_center_fields[0]->get_grid_data_field()->data_type_in_application) == 8, "Software error in Remap_grid_class::update_grid_center_3D_level_field_from_external");
	this->grid_center_fields[0]->interchange_grid_data(level_V3D_coord_trigger_field->get_coord_value_grid());
	memcpy(this->grid_center_fields[0]->get_grid_data_field()->data_buf, level_V3D_coord_trigger_field->get_grid_data_field()->data_buf, this->grid_center_fields[0]->get_grid_data_field()->required_data_size*get_data_type_size(level_V3D_coord_trigger_field->get_grid_data_field()->data_type_in_application));	
}


void Remap_grid_class::calculate_lev_sigma_values()
{
    long i, j, lev_grid_size;
    Remap_grid_class *leaf_grids[256], *lev_leaf_grid_of_sigma_or_hybrid, *lev_leaf_grid, *interchanged_grid;
    double *data_grid_field, *data_array_bot, *data_array_ptr, *data_sigma, data_bot, full_ratio, *hybrid_grid_coefficient_values, *tmp_vertical_coord_values;
    bool *horizontal_grid_mask_values;
    double local_sigma_grid_top_value, local_sigma_grid_scale_factor;
    

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->is_sigma_grid() && this->level_V3D_coord_trigger_field != NULL, "C-Coupler error2 in calculate_lev_sigma_values %s", grid_name);    
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(this->level_V3D_coord_trigger_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE),  "C-Coupler error in Remap_grid_class::calculate_lev_sigma_values: wrong data type of surface field");

    lev_leaf_grid_of_sigma_or_hybrid = get_a_leaf_grid_of_sigma_or_hybrid();
    lev_leaf_grid = get_a_leaf_grid(COORD_LABEL_LEV);

    hybrid_grid_coefficient_values = NULL;
    if (lev_leaf_grid_of_sigma_or_hybrid->hybrid_grid_coefficient_field != NULL) {
        hybrid_grid_coefficient_values = (double*) lev_leaf_grid_of_sigma_or_hybrid->hybrid_grid_coefficient_field->get_grid_data_field()->data_buf;
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(lev_leaf_grid_of_sigma_or_hybrid->hybrid_grid_coefficient_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE),  "C-Coupler error in Remap_grid_class::calculate_lev_sigma_values: wrong data type of hybrid coefficient");
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Sigma grid \"%s\" is a hybrid grid", grid_name);
    }

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(lev_leaf_grid_of_sigma_or_hybrid->get_sigma_grid_sigma_value_field()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE),    "C-Coupler error in Remap_grid_class::calculate_lev_sigma_values: wrong data type of sigma field");

    full_ratio = 0;
    data_sigma = (double*) lev_leaf_grid_of_sigma_or_hybrid->get_sigma_grid_sigma_value_field()->get_grid_data_field()->data_buf;
    lev_grid_size = lev_leaf_grid_of_sigma_or_hybrid->get_sigma_grid_sigma_value_field()->get_coord_value_grid()->grid_size;
    local_sigma_grid_top_value = lev_leaf_grid_of_sigma_or_hybrid->sigma_grid_top_value;
    local_sigma_grid_scale_factor = lev_leaf_grid_of_sigma_or_hybrid->sigma_grid_scale_factor;
    for (i = 0; i < lev_grid_size; i ++) {
        EXECUTION_REPORT(REPORT_WARNING, -1, fabs(data_sigma[i]) <= 1.0, "the sigma value in grid %s should be between 0.0 and 1.0\n", lev_leaf_grid_of_sigma_or_hybrid->get_sigma_grid_sigma_value_field()->get_coord_value_grid()->get_grid_name());
        if (fabs(data_sigma[i]) > 1.0)
            if (data_sigma[i] > 0)
                data_sigma[i] = 1.0;
            else data_sigma[i] = -1.0;
        if (full_ratio == 1.0)
            EXECUTION_REPORT(REPORT_ERROR, -1, data_sigma[i] >= 0, "the sigma value in grid %s must be all positive or negative", lev_leaf_grid_of_sigma_or_hybrid->get_sigma_grid_sigma_value_field()->get_coord_value_grid()->get_grid_name());
        if (full_ratio == -1.0)
            EXECUTION_REPORT(REPORT_ERROR, -1, data_sigma[i] <= 0, "the sigma value in grid %s must be all positive or negative", lev_leaf_grid_of_sigma_or_hybrid->get_sigma_grid_sigma_value_field()->get_coord_value_grid()->get_grid_name());
        if (data_sigma[i] > 0)
            full_ratio = 1.0;
        if (data_sigma[i] < 0)
            full_ratio = -1.0;
    }
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, full_ratio == 1.0 || full_ratio == -1.0, "the sigma value in grid %s must be all positive or negative", lev_leaf_grid_of_sigma_or_hybrid->get_sigma_grid_sigma_value_field()->get_coord_value_grid()->get_grid_name());
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_center_fields.size() == 1 && grid_center_fields[0]->grid_data_field->required_data_size == grid_size, "Software error in ");

    leaf_grids[0] = lev_leaf_grid;
    leaf_grids[1] = level_V3D_coord_trigger_field->get_coord_value_grid();
    interchanged_grid = new Remap_grid_class("TEMP_GRID_INTERCHANGE", 2, leaf_grids, 0);
    grid_center_fields[0]->interchange_grid_data(interchanged_grid);
    delete interchanged_grid;
    data_grid_field = (double*) grid_center_fields[0]->get_grid_data_field()->data_buf;
    level_V3D_coord_trigger_field->interchange_grid_data(level_V3D_coord_trigger_field->get_coord_value_grid());
    data_array_bot = (double*) level_V3D_coord_trigger_field->get_grid_data_field()->data_buf;
    if (level_V3D_coord_trigger_field->get_coord_value_grid()->get_grid_mask_field() != NULL)
        horizontal_grid_mask_values = (bool*) (level_V3D_coord_trigger_field->get_coord_value_grid()->get_grid_mask_field()->get_grid_data_field()->data_buf);
    else horizontal_grid_mask_values = NULL;
    tmp_vertical_coord_values = new double [lev_leaf_grid_of_sigma_or_hybrid->grid_size];
    for (i = 0; i < level_V3D_coord_trigger_field->get_coord_value_grid()->get_grid_size(); i ++) {
        data_bot = data_array_bot[i];
        if (horizontal_grid_mask_values != NULL && !horizontal_grid_mask_values[i]) {
            for (j = 0; j < lev_leaf_grid_of_sigma_or_hybrid->grid_size; j ++)
                tmp_vertical_coord_values[j] = data_array_bot[i];
        }
        else {
            if (hybrid_grid_coefficient_values == NULL) {
                for (j = 0; j < lev_leaf_grid_of_sigma_or_hybrid->grid_size; j ++)
                    tmp_vertical_coord_values[j] = (fabs(data_sigma[j])*(data_bot-local_sigma_grid_top_value) + local_sigma_grid_top_value)*full_ratio*local_sigma_grid_scale_factor;
            }
            else {
                for (j = 0; j < lev_leaf_grid_of_sigma_or_hybrid->grid_size; j ++)
                    tmp_vertical_coord_values[j] = (local_sigma_grid_top_value*hybrid_grid_coefficient_values[j]+fabs(data_sigma[j])*data_bot)*full_ratio*local_sigma_grid_scale_factor;        
            }
        }
        data_array_ptr = data_grid_field + i*lev_leaf_grid->grid_size;
        if (lev_leaf_grid == lev_leaf_grid_of_sigma_or_hybrid) {
            for (j = 0; j < lev_leaf_grid->grid_size; j ++)
                data_array_ptr[j] = tmp_vertical_coord_values[j];
        }
        else {
            for (j = 0; j < lev_leaf_grid->grid_size; j ++)
                data_array_ptr[j] = (tmp_vertical_coord_values[j]+tmp_vertical_coord_values[j+1])/2;
        }
    }
    delete [] tmp_vertical_coord_values;    
}


Remap_grid_data_class *Remap_grid_class::get_unique_center_field()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_center_fields.size() == 1 && (is_sigma_grid() || does_use_V3D_level_coord()), "C-Coupler error in Remap_grid_class::get_unique_center_field");
    return grid_center_fields[0];
}


Remap_grid_data_class *Remap_grid_class::get_unique_vertex_field()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_vertex_fields.size() == 1 && is_sigma_grid() && num_vertexes == 2, "C-Coupler error in Remap_grid_class::get_unique_vertex_field");
    return grid_vertex_fields[0];
}


void Remap_grid_class::read_grid_data_from_array(const char *coord_label, const char *coord_name, const char *data_type, const char *array_data, int num_vertexes)
{
    Remap_data_field *remap_data_field;
    Remap_grid_data_class *remap_grid_data_field;
    Remap_grid_class *similar_grid;


    if (num_vertexes != 0)
        this->num_vertexes = num_vertexes;

    EXECUTION_REPORT(REPORT_ERROR, -1, grid_size > 0, "Software error in Remap_grid_class::read_grid_data_from_array: wrong grid size");

    remap_data_field = new Remap_data_field;    
    if (words_are_the_same(coord_label, GRID_MASK_LABEL)) {
        strcpy(remap_data_field->field_name_in_application, GRID_MASK_LABEL);
        strcpy(remap_data_field->data_type_in_application, DATA_TYPE_BOOL);
        check_mask_value_can_be_set();
    }
    else {
        if (words_are_the_same(coord_label, GRID_VERTEX_LABEL))
            check_vertex_coord_value_can_be_set(coord_name);
        else if (words_are_the_same(coord_label, GRID_CENTER_LABEL))
            check_center_coord_value_can_be_set(coord_name);
        strcpy(remap_data_field->field_name_in_application, coord_name);
        strcpy(remap_data_field->data_type_in_application, DATA_TYPE_DOUBLE);
    }
    remap_data_field->required_data_size = grid_size;
    if (words_are_the_same(coord_label, GRID_VERTEX_LABEL))  
        remap_data_field->required_data_size *= num_vertexes;
    remap_data_field->read_data_size = remap_data_field->required_data_size;
    remap_data_field->data_buf = new char [remap_data_field->required_data_size*get_data_type_size(remap_data_field->data_type_in_application)];
    if (words_are_the_same(coord_label, GRID_VERTEX_LABEL))
        sprintf(remap_data_field->field_name_in_IO_file, "%s_%s", GRID_VERTEX_LABEL, coord_name);
    else strcpy(remap_data_field->field_name_in_IO_file, coord_name);
    if (words_are_the_same(coord_label, GRID_MASK_LABEL)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(data_type, DATA_TYPE_INT), "Software error in Remap_grid_class::read_grid_data_from_array: wrong data type for mask");
        for (int i = 0; i < grid_size; i ++)
            ((bool*)remap_data_field->data_buf)[i] = ((int*) array_data)[i] == 1;
    }
    else {
        int array_size = grid_size; 
        if (num_vertexes > 0)
            array_size *= num_vertexes;
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(data_type, DATA_TYPE_FLOAT) || words_are_the_same(data_type, DATA_TYPE_DOUBLE), "Software error in Remap_grid_class::read_grid_data_from_array: wrong data type for others");
        if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
            for (int i = 0; i < array_size; i ++)
                ((double*)remap_data_field->data_buf)[i] = ((float*) array_data)[i];
        }
        else {
            for (int i = 0; i < array_size; i ++)
                ((double*)remap_data_field->data_buf)[i] = ((double*) array_data)[i];            
        }
    }
    remap_grid_data_field = new Remap_grid_data_class(this, remap_data_field);
    if (words_are_the_same(coord_label, GRID_MASK_LABEL))
        grid_mask_field = remap_grid_data_field;
    else if (words_are_the_same(coord_label, GRID_CENTER_LABEL) || words_are_the_same(coord_label, GRID_VERTEX_LABEL)) {
        transform_coord_values_from_radian_to_degrees(remap_grid_data_field);
        formalize_cyclic_coord_values(remap_grid_data_field);
        if (words_are_the_same(coord_label, GRID_CENTER_LABEL))
            grid_center_fields.push_back(remap_grid_data_field);
        else grid_vertex_fields.push_back(remap_grid_data_field);
    }
    else imported_area = remap_grid_data_field;

    if (words_are_the_same(coord_name, COORD_LABEL_LON) || words_are_the_same(coord_name, COORD_LABEL_LAT) || words_are_the_same(coord_name, COORD_LABEL_LAT)) {
        Remap_grid_class *leaf_grids[256]; 
        int num_leaf_grids;
        get_leaf_grids(&num_leaf_grids, leaf_grids, this);
        for (int i = 0; i < num_leaf_grids; i ++)
            if (words_are_the_same(leaf_grids[i]->get_coord_label(), coord_name)) {
                if (words_are_the_same(coord_name, COORD_LABEL_LON))
                    remap_grid_data_field->get_grid_data_field()->set_field_long_name("longitude");
                else if (words_are_the_same(coord_name, COORD_LABEL_LAT))
                    remap_grid_data_field->get_grid_data_field()->set_field_long_name("latitude");
                remap_grid_data_field->get_grid_data_field()->set_field_unit(leaf_grids[i]->get_coord_unit());
            }
    }

    strcpy(remap_data_field->data_type_in_IO_file, data_type);
}


void Remap_grid_class::read_grid_data_from_IO(char extension_names[16][256], 
                                              const char *IO_object_name, 
                                              const char *variable_name_in_IO_file,
                                              int num_vertexes)
{
    Remap_data_field *remap_data_field;
    Remap_grid_data_class *remap_grid_data_field;
    Remap_grid_class *similar_grid;


    if (words_are_the_same(extension_names[1], COORD_LABEL_TRACER))
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "can not set the coordinate values of grid %s\n", extension_names[1]);
    if ((words_are_the_same(extension_names[0], GRID_CENTER_LABEL) || words_are_the_same(extension_names[0], GRID_VERTEX_LABEL))) {
        if (words_are_the_same(extension_names[1], COORD_LABEL_LON) || words_are_the_same(extension_names[1], COORD_LABEL_LAT)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions <= 2, "the coordinate values of %s can only be set in 1D or 2D grid\n", extension_names[1]);
            if (num_dimensions == 2)
                EXECUTION_REPORT(REPORT_ERROR, -1, this->get_is_sphere_grid(), "the 2D grid %s to be set lon or lat values must be a sphere grid with coordinates lon, and lat\n",
                             this->grid_name);
        }
        if (words_are_the_same(extension_names[1], COORD_LABEL_TIME))
            EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 1, "the coordinate values of %s can only be set in 1D grid\n", extension_names[1]);
        if (words_are_the_same(extension_names[1], COORD_LABEL_LEV)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 1 || num_dimensions == 3, "the coordinate values of %s can only be set in 1D or 3D grid\n", extension_names[1]);
            if (num_dimensions == 3)
                EXECUTION_REPORT(REPORT_ERROR, -1, this->has_grid_coord_label(COORD_LABEL_LON) && this->has_grid_coord_label(COORD_LABEL_LAT),
                             "the 3D grid %s to be set lev values must be a spacial grid with coordinates lon, lat and lev\n",
                             this->grid_name);
        }
    }
    if (words_are_the_same(extension_names[0], GRID_VERTEX_LABEL) && this->num_vertexes != 0) {
        EXECUTION_REPORT(REPORT_ERROR, -1, this->num_vertexes == num_vertexes, 
                     "the number of vertexes of grid \"%s\" specified by input parameter is different from the numbe of vertexes specified by the previous statements\n",
                     grid_name);
    }

    if (num_vertexes != 0)
        this->num_vertexes = num_vertexes;

    EXECUTION_REPORT(REPORT_ERROR, -1, grid_size > 0, "the grid size of grid object \"%s\" is unknown, and so can not read grid field data from IO file\n", grid_name);

    remap_data_field = new Remap_data_field;    
    if (words_are_the_same(extension_names[0], GRID_MASK_LABEL)) {
        strcpy(remap_data_field->field_name_in_application, GRID_MASK_LABEL);
        strcpy(remap_data_field->data_type_in_application, DATA_TYPE_BOOL);
        check_mask_value_can_be_set();
    }
    else {
        if (words_are_the_same(extension_names[0], GRID_VERTEX_LABEL))
            check_vertex_coord_value_can_be_set(extension_names[1]);
        else check_center_coord_value_can_be_set(extension_names[1]);
        strcpy(remap_data_field->field_name_in_application, extension_names[1]);
        strcpy(remap_data_field->data_type_in_application, DATA_TYPE_DOUBLE);
    }
    remap_data_field->required_data_size = grid_size;
    if (words_are_the_same(extension_names[0], GRID_VERTEX_LABEL))  
        remap_data_field->required_data_size *= num_vertexes;
    remap_data_field->data_buf = new char [remap_data_field->required_data_size*get_data_type_size(remap_data_field->data_type_in_application)];
    strcpy(remap_data_field->field_name_in_IO_file, variable_name_in_IO_file);
    io_manager->read_data(IO_object_name, remap_data_field, true);
    remap_grid_data_field = new Remap_grid_data_class(this, remap_data_field);
    if (words_are_the_same(extension_names[0], GRID_MASK_LABEL))
        grid_mask_field = remap_grid_data_field;
    else {
        transform_coord_values_from_radian_to_degrees(remap_grid_data_field);
        formalize_cyclic_coord_values(remap_grid_data_field);
        if (words_are_the_same(extension_names[0], GRID_CENTER_LABEL))
            grid_center_fields.push_back(remap_grid_data_field);
        else grid_vertex_fields.push_back(remap_grid_data_field);
    }
}


bool Remap_grid_class::has_grid_coord_label(const char *coord_label) const
{
    int num_leaf_grids, i;
    Remap_grid_class *leaf_grids[256];


    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++) 
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), coord_label))
            return true;

    return false;
}


bool Remap_grid_class::match_grid(const char *grid_name) const
{
    return words_are_the_same(this->grid_name, grid_name);
}


bool Remap_grid_class::match_grid(int num_sized_sub_grids, Remap_grid_class **sized_sub_grids)
{
    int this_num_sized_sub_grids, i;
    Remap_grid_class *this_sized_sub_grids[256];


    if (this->grid_size == 0)
        return false;

    this->get_sized_sub_grids(&this_num_sized_sub_grids, this_sized_sub_grids);
    if (num_sized_sub_grids != this_num_sized_sub_grids)
        return false;
    for (i = 0; i < num_sized_sub_grids; i ++)
        if (sized_sub_grids[i] != this_sized_sub_grids[i])
            break;
    return i == num_sized_sub_grids;
}


bool Remap_grid_class::have_overlap_with_grid(Remap_grid_class *another_grid) 
{
    int num_leaf_grids_this, num_leaf_grids_another, i, j;
    Remap_grid_class *leaf_grids_this[128], *leaf_grids_another[128];


    get_leaf_grids(&num_leaf_grids_this, leaf_grids_this, this);
    another_grid->get_leaf_grids(&num_leaf_grids_another, leaf_grids_another, another_grid);
    for (i = 0; i < num_leaf_grids_this; i ++)
        for (j = 0; j < num_leaf_grids_another; j ++)
            if (leaf_grids_this[i] == leaf_grids_another[j])
                return true;

    return false;
}


bool Remap_grid_class::is_subset_of_grid(Remap_grid_class *another_grid) 
{
    int num_leaf_grids_this, num_leaf_grids_another, i, j;
    Remap_grid_class *leaf_grids_this[128], *leaf_grids_another[128];


    get_leaf_grids(&num_leaf_grids_this, leaf_grids_this, this);
    another_grid->get_leaf_grids(&num_leaf_grids_another, leaf_grids_another, another_grid);
    for (i = 0; i < num_leaf_grids_this; i ++) {
        for (j = 0; j < num_leaf_grids_another; j ++)
            if (words_are_the_same(leaf_grids_this[i]->grid_name, leaf_grids_another[j]->grid_name))
                break;
        if (j == num_leaf_grids_another)
            return false;
    }

    if (this->get_num_dimensions() == another_grid->get_num_dimensions())
        return this->grid_size == another_grid->grid_size;

    return true;
}


bool Remap_grid_class::is_the_same_grid_with(Remap_grid_class *another_grid) 
{
    int num_leaf_grids_this, num_leaf_grids_another, i;
    Remap_grid_class *leaf_grids_this[128], *leaf_grids_another[128];


    get_leaf_grids(&num_leaf_grids_this, leaf_grids_this, this);
    another_grid->get_leaf_grids(&num_leaf_grids_another, leaf_grids_another, another_grid);

    if (num_leaf_grids_this != num_leaf_grids_another)
        return false;

    for (i = 0; i < num_leaf_grids_this; i ++)
        if (!words_are_the_same(leaf_grids_this[i]->grid_name, leaf_grids_another[i]->grid_name))
            return false;

    return true;
}


bool Remap_grid_class::is_similar_grid_with(Remap_grid_class *another_grid) 
{
    return is_subset_of_grid(another_grid) && another_grid->is_subset_of_grid(this);
}


bool Remap_grid_class::is_superset_of_grid(Remap_grid_class *another_grid) 
{
    return another_grid->is_subset_of_grid(this);
}


Remap_grid_class *Remap_grid_class::get_a_leaf_grid(const char *coord_label)
{
    Remap_grid_class *leaf_grids[256];
    int i, num_leaf_grids;


    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->coord_label, coord_label))
            return leaf_grids[i];

    EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error in get_a_leaf_grid");

    return NULL;
}


void Remap_grid_class::get_leaf_grids(int *num_leaf_grids, Remap_grid_class **leaf_grids, const Remap_grid_class *top_grid) const
{
    if (this == top_grid)
        *num_leaf_grids = 0;

    if (whole_grid != NULL) {
        whole_grid->get_leaf_grids(num_leaf_grids, leaf_grids, top_grid);
        return;
    }
    
    if (sub_grids.size() == 0)
        leaf_grids[(*num_leaf_grids)++] = (Remap_grid_class*) this;

    for (int i = 0; i < sub_grids.size(); i ++)
        sub_grids[i]->get_leaf_grids(num_leaf_grids, leaf_grids, top_grid);
}


void Remap_grid_class::get_sized_sub_grids(int *num_sized_sub_grids, Remap_grid_class **sized_sub_grids) 
{
    int num_leaf_grids, tmp_num_dimensions, i, j;
    Remap_grid_class *sized_super_grid;
    Remap_grid_class *leaf_grids[256];
    long tmp_grid_size;


    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    *num_sized_sub_grids = 0;
    for (i = 0; i < num_leaf_grids; i ++) {
        if (leaf_grids[i] == NULL)
            continue;
        sized_super_grid = leaf_grids[i]->get_first_super_grid_of_enable_setting_coord_value();
        EXECUTION_REPORT(REPORT_ERROR, -1, sized_super_grid != NULL, "remap software error1 in function get_sized_sub_grids\n");
        sized_sub_grids[(*num_sized_sub_grids)++] = sized_super_grid;
        if (sized_super_grid != leaf_grids[i])
            for (j = i+1; j < num_leaf_grids; j ++) {
                if (leaf_grids[j] == NULL)
                    continue;
                if (leaf_grids[j]->is_subset_of_grid(sized_super_grid)) {
                    EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[j]->get_grid_size() == 0, "remap software error2 in function get_sized_sub_grids\n");
                    leaf_grids[j] = NULL;
                }
            }
    }

    tmp_grid_size = 1;
    tmp_num_dimensions = 0;
    for (i = 0; i < *num_sized_sub_grids; i ++) {
        tmp_grid_size *= sized_sub_grids[i]->get_grid_size();
        tmp_num_dimensions += sized_sub_grids[i]->get_num_dimensions();
        if (original_grid == NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, is_superset_of_grid(sized_sub_grids[i]), 
                         "remap software error3 in function get_sized_sub_grids\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, !sized_sub_grids[i]->is_partial_grid(), 
                     "remap software error4 in function get_sized_sub_grids\n");
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, tmp_grid_size == grid_size && tmp_num_dimensions == num_dimensions, 
                 "remap software error5 in function get_sized_sub_grids\n");
}


void Remap_grid_class::get_masked_sub_grids(int *num_masked_sub_grids, Remap_grid_class **masked_sub_grids) 
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];  
    Remap_grid_class *super_grid;  
    int i, j;


    *num_masked_sub_grids = 0;
    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++) {
        if (leaf_grids[i] == NULL)
            continue;
        if (leaf_grids[i]->super_grids_of_setting_mask_value.size() == 0)
            continue;
        super_grid = NULL;
        for (j = 0; j < leaf_grids[i]->super_grids_of_setting_mask_value.size(); j ++) {
            if (leaf_grids[i]->super_grids_of_setting_mask_value[j]->is_subset_of_grid(this)) {
                if (super_grid == NULL || super_grid->is_subset_of_grid(leaf_grids[i]->super_grids_of_setting_mask_value[j]))
                    super_grid = leaf_grids[i]->super_grids_of_setting_mask_value[j];
            }
        }
        if (super_grid == NULL)
            continue;
        masked_sub_grids[(*num_masked_sub_grids)++] = super_grid;
        for (j = i+1; j < num_leaf_grids; j ++) {
            if (leaf_grids[j] == NULL)
                continue;
            if (leaf_grids[j]->is_subset_of_grid(super_grid))
                leaf_grids[j] = NULL;
        }
    }

    for (i = 0; i < *num_masked_sub_grids; i ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, masked_sub_grids[i]->grid_mask_field != NULL, "remap software error1 in get_masked_sub_grids\n");
        for (j = i+1; j < *num_masked_sub_grids; j ++)
            EXECUTION_REPORT(REPORT_ERROR, -1, !masked_sub_grids[i]->have_overlap_with_grid(masked_sub_grids[j]), "remap software error2 in get_masked_sub_grids\n");
    }
}


void Remap_grid_class::get_partial_grid_mask_fields(int *num_mask_fields_partial, 
                                                    Remap_grid_data_class **mask_fields_partial, 
                                                    Remap_grid_class *top_grid)
{
    if (this == top_grid)
        *num_mask_fields_partial = 0;

    for (int i = 0; i < sub_grids.size(); i ++)
        sub_grids[i]->get_partial_grid_mask_fields(num_mask_fields_partial, mask_fields_partial, top_grid);

    if (whole_grid != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_mask_field != NULL, "remap software error in get_partial_grid_mask_fields\n");
        mask_fields_partial[(*num_mask_fields_partial)++] = grid_mask_field;
    }
}


Remap_grid_data_class *Remap_grid_class::get_grid_center_field(const char *label) const
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];
    
    
    get_leaf_grids(&num_leaf_grids, leaf_grids, this);

    for (int i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), label))
            return leaf_grids[i]->get_grid_center_field();

    EXECUTION_REPORT(REPORT_ERROR, -1, "Software error in Remap_grid_class::get_grid_center_field: fail to find a leaf grid with coordinate \"%s\"", label);

    return NULL;
}


Remap_grid_data_class *Remap_grid_class::get_grid_center_field() const
{
    Remap_grid_class *super_grid;
    int i;
    Remap_grid_data_class *existing_data_field = NULL;
    

    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_dimensions == 1, "remap software error1 in get_grid_center_field\n");    
    
    super_grid = this->get_super_grid_of_setting_coord_values();
    if (super_grid == NULL)
        return NULL;

    for (i = 0; i < super_grid->grid_center_fields.size(); i ++)
        if (words_are_the_same(super_grid->grid_center_fields[i]->grid_data_field->field_name_in_application, this->coord_label)) {
            existing_data_field = super_grid->grid_center_fields[i];
            break;
        }

    EXECUTION_REPORT(REPORT_ERROR, -1, existing_data_field != NULL, "remap software error2 in get_grid_center_field\n");

    return existing_data_field;
}


Remap_grid_data_class *Remap_grid_class::get_grid_vertex_field() const
{
    Remap_grid_class *super_grid;
    int i;
    Remap_grid_data_class *existing_data_field = NULL;
    

    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_dimensions == 1,
                 "remap software error1 in get_grid_vertex_field\n");    
    
    super_grid = this->get_super_grid_of_setting_coord_values();
    if (super_grid == NULL)
        return NULL;

    for (i = 0; i < super_grid->grid_vertex_fields.size(); i ++)
        if (words_are_the_same(super_grid->grid_vertex_fields[i]->grid_data_field->field_name_in_application, this->coord_label)) {
            existing_data_field = super_grid->grid_vertex_fields[i];
            break;
        }

    return existing_data_field;
}


void Remap_grid_class::check_grid_field_can_be_set()
{
    Remap_grid_class *similar_grid;


    EXECUTION_REPORT(REPORT_ERROR, -1, this->whole_grid == NULL, "\"%s\" is a partial grid, whose coordinate values can not be set\n", this->grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, this->enable_to_set_coord_values, "can not set the coordinate values of grid \"%s\", because it or its super grids have been used for remapping\n", this->grid_name);
    similar_grid = get_similar_grids_setting_coord_values();
    if (similar_grid != NULL)
        EXECUTION_REPORT(REPORT_ERROR, -1, false,
                     "can not set coordinate values in grid \"%s\" because its similar grid \"%s\" have been used to set coordinate values: all coordinated values can only be set in one of similar grids\n",
                     grid_name, similar_grid->grid_name);
}


void Remap_grid_class::check_vertex_coord_value_can_be_set(const char *coord_label)
{
    int i;


    check_grid_field_can_be_set();
    
    for (i = 0; i < grid_center_fields.size(); i ++)
        if (words_are_the_same(grid_center_fields[i]->grid_data_field->field_name_in_application, coord_label))
            break;

    EXECUTION_REPORT(REPORT_ERROR, -1, i < grid_center_fields.size(), 
                 "the vertex values of coordinate \"%s\" can not be set in grid object \"%s\" because the corresonding center values of coordinate \"%s\" have not been set\n",
                 coord_label, grid_name, coord_label);

    for (i = 0; i < grid_vertex_fields.size(); i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, !words_are_the_same(grid_vertex_fields[i]->grid_data_field->field_name_in_application, coord_label), 
                     "the vertex values of coordinate \"%s\" can not be set in grid object \"%s\" because these values have been set before\n",
                     coord_label, grid_name);
}


void Remap_grid_class::check_center_coord_value_can_be_set(const char *coord_label)
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];
    int i, j;


    check_grid_field_can_be_set();

    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(coord_label, leaf_grids[i]->coord_label)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->get_grid_center_field() == NULL,
                         "grid %s and %s have the same leaf grid object %s, they can not both set the coordinate values corresponding to this leaf grid\n", 
                         this->get_grid_name(), leaf_grids[i]->super_grid_of_setting_coord_values->get_grid_name(), grid_name);            
            leaf_grids[i]->super_grid_of_setting_coord_values = this;
            for (j = 0; j < remap_grid_manager->remap_grids.size(); j ++)
                if (leaf_grids[i]->is_subset_of_grid(remap_grid_manager->remap_grids[j]))
                    EXECUTION_REPORT(REPORT_ERROR, -1, this->is_subset_of_grid(remap_grid_manager->remap_grids[j]) ||
                                 remap_grid_manager->remap_grids[j]->is_subset_of_grid(this),
                                 "grid %s have coordinate %s. However, current grid %s which sets the coordinate values is not a subgrid or super grid of this grid\n",
                                 remap_grid_manager->remap_grids[j]->get_grid_name(), leaf_grids[i]->get_coord_label(), this->get_grid_name());
            break;
        }
        
    EXECUTION_REPORT(REPORT_ERROR, -1, i < num_leaf_grids, "remap software error in check_center_coord_value_can_be_set\n");
}


void Remap_grid_class::check_mask_value_can_be_set()
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];
    int i;


    check_grid_field_can_be_set();

    EXECUTION_REPORT(REPORT_ERROR, -1, grid_mask_field == NULL, 
                 "the mask values of grid object \"%s\" can not be set again because these mask values have been set before\n",
                 grid_name);

    EXECUTION_REPORT(REPORT_ERROR, -1, !masks_are_known, 
                 "the mask values of grid \"%s\" can not be set because these mask values can be computed according to the sub grids\n",
                 grid_name);

    masks_are_known = true;

    for (i = 0; i < remap_grid_manager->remap_grids.size(); i ++)
        if (this != remap_grid_manager->remap_grids[i] && 
            this->is_similar_grid_with(remap_grid_manager->remap_grids[i]))
            EXECUTION_REPORT(REPORT_ERROR, -1, !remap_grid_manager->remap_grids[i]->masks_are_known, 
                         "the mask values of grid object \"%s\" can not be set because the mask values of its similar grid \"%s\" mask values are known\n",
                         grid_name, remap_grid_manager->remap_grids[i]->grid_name);

    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++) 
        leaf_grids[i]->super_grids_of_setting_mask_value.push_back(this);
}


void Remap_grid_class::set_1D_coord_vertex_values_in_default(const double *center_values, 
                                                            double *vertex_values,
                                                            long grid_size, 
                                                            bool cyclic,
                                                            bool is_unit_degree)
{
    long i, prev_i, next_i;
    double data_value1, data_value2, data_value3, data_value4;


    for (i = 0; i < grid_size; i ++) {
        if (cyclic) {
            prev_i = (i-1+grid_size) % grid_size;
            while (redundant_cell_mark != NULL && redundant_cell_mark[prev_i])
                prev_i = (prev_i-1+grid_size) % grid_size;
            next_i = (i+1) % grid_size;
            while (redundant_cell_mark != NULL && redundant_cell_mark[next_i])
                next_i = (next_i+1) % grid_size;
        }
        else {
            prev_i = i - 1;
            while (prev_i >= 0 && (redundant_cell_mark != NULL && redundant_cell_mark[prev_i]))
                prev_i --;
            if (prev_i < 0)
                prev_i = i;
            next_i = i + 1;
            while (next_i < grid_size && (redundant_cell_mark != NULL && redundant_cell_mark[next_i]))
                next_i ++;
            if (next_i >= grid_size)
                next_i = i;
        }
        data_value1 = center_values[prev_i];
        data_value4 = center_values[next_i];
        data_value2 = center_values[i];
        data_value3 = data_value2;
        if (is_unit_degree) {
            match_degree_values(data_value1, data_value2);
            match_degree_values(data_value4, data_value3);        
        }
        vertex_values[2*i] = (data_value1+data_value2)/2;
        vertex_values[2*i+1] = (data_value3+data_value4)/2;
    }
}


void Remap_grid_class::set_2D_coord_vertex_values_in_default(const double *center_values_dim1, 
                                                             double *vertex_values_dim1,
                                                             long grid_size_dim1, 
                                                             bool cyclic_dim1,
                                                             bool is_unit_degree_dim1,
                                                             const double *center_values_dim2, 
                                                             double *vertex_values_dim2,
                                                             long grid_size_dim2, 
                                                             bool cyclic_dim2,
                                                             bool is_unit_degree_dim2)
{
    double **temp_vertex_values_dim1, **temp_vertex_values_dim2;
    long i, j, k, l, indx1, indx2, left_i, right_i, top_j, bot_j, cell_index;
    long box_vertex_start_dim1, box_vertex_end_dim1, box_vertex_start_dim2, box_vertex_end_dim2;
    double max_coord_value, min_coord_value, sum_value_dim1, sum_value_dim2;
    int num_non_null_cells;
    bool cross_360_degree_dim1;


    EXECUTION_REPORT(REPORT_ERROR, -1, get_is_sphere_grid(), "grid \"%s\" does not match the constraint that the grid which can be set 2D vertex values must be 2D sphere grid", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, !cyclic_dim2, "remap software error1 in set_2D_coord_vertex_values_in_default\n");

    temp_vertex_values_dim1 = new double* [grid_size_dim2+1];
    temp_vertex_values_dim2 = new double* [grid_size_dim2+1];
    for (i = 0; i < grid_size_dim2+1; i ++) {
        temp_vertex_values_dim1[i] = new double [grid_size_dim1+1];
        temp_vertex_values_dim2[i] = new double [grid_size_dim1+1];
    }

    for (j = 0; j < grid_size_dim2+1; j ++) {
        for (i = 0; i < grid_size_dim1+1; i ++) {
            box_vertex_start_dim1 = i - 1;
            box_vertex_end_dim1 = i;
            box_vertex_start_dim2 = j - 1;
            box_vertex_end_dim2 = j;
            if (box_vertex_start_dim2 < 0) 
                box_vertex_start_dim2 = 0;
            if (box_vertex_end_dim2 >= grid_size_dim2)
                box_vertex_end_dim2 = grid_size_dim2 - 1;
            if (!cyclic_dim1) {
				if (redundant_cell_mark != NULL)
	                while ((redundant_cell_mark[box_vertex_start_dim2*grid_size_dim1+box_vertex_start_dim1] ||
    	                   redundant_cell_mark[box_vertex_end_dim2*grid_size_dim1+box_vertex_start_dim1]) &&
        	               box_vertex_start_dim1 >= 0)
            	        box_vertex_start_dim1 --;
				if (redundant_cell_mark != NULL)
	                while ((redundant_cell_mark[box_vertex_start_dim2*grid_size_dim1+box_vertex_end_dim1] ||
    	                   redundant_cell_mark[box_vertex_end_dim2*grid_size_dim1+box_vertex_end_dim1]) &&
        	               box_vertex_end_dim1 < grid_size_dim1)
            	        box_vertex_end_dim1 ++;
                if (box_vertex_start_dim1 < 0)
                    box_vertex_start_dim1 = box_vertex_end_dim1;
                if (box_vertex_end_dim1 >= grid_size_dim1) 
                    box_vertex_end_dim1 = box_vertex_start_dim1;
                EXECUTION_REPORT(REPORT_ERROR, -1, box_vertex_start_dim1 >= 0 && box_vertex_end_dim1 < grid_size_dim1, "remap software error2 in set_2D_coord_vertex_values_in_default\n");
            }
            else {
                box_vertex_start_dim1 = (box_vertex_start_dim1+grid_size_dim1) % grid_size_dim1;
				if (redundant_cell_mark != NULL)
	                while (redundant_cell_mark[box_vertex_start_dim2*grid_size_dim1+box_vertex_start_dim1] ||
    	                   redundant_cell_mark[box_vertex_end_dim2*grid_size_dim1+box_vertex_start_dim1])
        	            box_vertex_start_dim1 = (box_vertex_start_dim1+grid_size_dim1-1) % grid_size_dim1;
                box_vertex_end_dim1 = (box_vertex_end_dim1+grid_size_dim1) % grid_size_dim1;
				if (redundant_cell_mark != NULL)
	                while (redundant_cell_mark[box_vertex_start_dim2*grid_size_dim1+box_vertex_end_dim1] ||
    	                   redundant_cell_mark[box_vertex_end_dim2*grid_size_dim1+box_vertex_end_dim1])
        	            box_vertex_end_dim1 = (box_vertex_end_dim1+grid_size_dim1+1) % grid_size_dim1;
            }
            num_non_null_cells = 0;
            sum_value_dim1 = 0;
            sum_value_dim2 = 0;
            if (box_vertex_start_dim1 > box_vertex_end_dim1)
                box_vertex_end_dim1 += grid_size_dim1;
            cross_360_degree_dim1 = false;
            if (is_unit_degree_dim1) {
                min_coord_value = (double)DEFAULT_FILL_VALUE;
                max_coord_value = -min_coord_value;
                for (k = box_vertex_start_dim2; k <= box_vertex_end_dim2; k ++)
                    for (l = box_vertex_start_dim1; l <= box_vertex_end_dim1; l ++) {
                        indx1 = l%grid_size_dim1;
                        indx2 = k%grid_size_dim2;
                        if (redundant_cell_mark == NULL || !redundant_cell_mark[indx2*grid_size_dim1+indx1]) {
                            if (max_coord_value < center_values_dim1[indx2*grid_size_dim1+indx1])
                                max_coord_value = center_values_dim1[indx2*grid_size_dim1+indx1];
                            if (min_coord_value > center_values_dim1[indx2*grid_size_dim1+indx1])
                                min_coord_value = center_values_dim1[indx2*grid_size_dim1+indx1];
                        }
                    }
                if (max_coord_value - min_coord_value >= 180)
                    cross_360_degree_dim1 = true;
            }

            for (k = box_vertex_start_dim2; k <= box_vertex_end_dim2; k ++)
                for (l = box_vertex_start_dim1; l <= box_vertex_end_dim1; l ++) {
                    indx1 = l%grid_size_dim1;
                    indx2 = k%grid_size_dim2;
                    if (redundant_cell_mark == NULL || !redundant_cell_mark[indx2*grid_size_dim1+indx1]) {
                        sum_value_dim1 += center_values_dim1[indx2*grid_size_dim1+indx1];
                        sum_value_dim2 += center_values_dim2[indx2*grid_size_dim1+indx1];
                        num_non_null_cells ++;
                        if (cross_360_degree_dim1 && center_values_dim1[indx2*grid_size_dim1+indx1] > 180)
                            sum_value_dim1 -= 360;
                    }
                }
            temp_vertex_values_dim1[j][i] = sum_value_dim1/num_non_null_cells;
            temp_vertex_values_dim2[j][i] = sum_value_dim2/num_non_null_cells;
            if (cross_360_degree_dim1)
                while (temp_vertex_values_dim1[j][i] < 0)
                    temp_vertex_values_dim1[j][i] += 360;
        }
    }

    for (j = 0; j < grid_size_dim2; j ++) {
        for (i = 0; i < grid_size_dim1; i ++) {
            left_i = i; 
            right_i = i+1;
            top_j = j;
            bot_j = j+1;
            cell_index = j*grid_size_dim1 + i;
            vertex_values_dim1[cell_index*4+0] = temp_vertex_values_dim1[top_j][left_i];
            vertex_values_dim1[cell_index*4+1] = temp_vertex_values_dim1[bot_j][left_i];
            vertex_values_dim1[cell_index*4+2] = temp_vertex_values_dim1[bot_j][right_i];
            vertex_values_dim1[cell_index*4+3] = temp_vertex_values_dim1[top_j][right_i];
            vertex_values_dim2[cell_index*4+0] = temp_vertex_values_dim2[top_j][left_i];
            vertex_values_dim2[cell_index*4+1] = temp_vertex_values_dim2[bot_j][left_i];
            vertex_values_dim2[cell_index*4+2] = temp_vertex_values_dim2[bot_j][right_i];
            vertex_values_dim2[cell_index*4+3] = temp_vertex_values_dim2[top_j][right_i];
        }
    }
    
    for (i = 0; i < grid_size_dim2+1; i ++) {
        delete [] temp_vertex_values_dim1[i];
        delete [] temp_vertex_values_dim2[i];
    }
    delete [] temp_vertex_values_dim1;
    delete [] temp_vertex_values_dim2;
}


Remap_grid_class *Remap_grid_class::get_similar_grids_setting_coord_values()
{
    for (int i = 0; i < remap_grid_manager->remap_grids.size(); i ++) {
        if (remap_grid_manager->remap_grids[i]->is_partial_grid())
            continue;
        if (remap_grid_manager->remap_grids[i] != this && remap_grid_manager->remap_grids[i]->is_similar_grid_with(this)) 
            if (remap_grid_manager->remap_grids[i]->grid_center_fields.size() > 0 || 
                remap_grid_manager->remap_grids[i]->grid_mask_field != NULL) {
                return remap_grid_manager->remap_grids[i];
            }
    }

    return NULL;
}


void Remap_grid_class::set_grid_boundary(double min_lon, double max_lon, double min_lat, double max_lat)
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256], *lon_sub_grid;
    double eps = 1.0000001;


    EXECUTION_REPORT(REPORT_ERROR, -1, this->get_is_sphere_grid(), "%s is not a sphere grid, while only sphere grid can be set boundary", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, boundary_min_lon == NULL_COORD_VALUE, "the boundary of grid %s has been set before", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, min_lon >= -360*eps && min_lon <= 360*eps, "the minimum longitude (%lf) of the boundary of grid %s must be between 0 and 360 degrees", min_lon, grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, max_lon >= -360*eps && max_lon <= 360*eps, "the maximum longitude (%lf) of the boundary of grid %s must be between 0 and 360 degrees", max_lon, grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, min_lat >= -90*eps && min_lat <= 90*eps, "the minimum latitude (%lf) of the boundary of grid %s must be between -90 and +90 degrees", min_lat, grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, max_lat >= -90*eps && max_lat <= 90*eps, "the maximum latitude (%lf) of the boundary of grid %s must be between -90 and +90 degrees", max_lat, grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, min_lat < max_lat, "the minimum latitude (%lf) of the boundary must be smaller than the maximum latitude (%lf) of the boundary", min_lat, max_lat, grid_name);

    if (are_floating_values_equal(min_lat, (double)-90))
        min_lat = -90;
    if (are_floating_values_equal(max_lat, (double)90))
        max_lat = 90;

    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    if (words_are_the_same(leaf_grids[0]->get_coord_label(), COORD_LABEL_LON))
        lon_sub_grid = leaf_grids[0];
    else lon_sub_grid = leaf_grids[1];

    boundary_min_lon = min_lon;
    boundary_max_lon = max_lon;
    boundary_min_lat = min_lat;
    boundary_max_lat = max_lat;
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Set boundary of grid \"%s\" to %lf %lf %lf %lf", grid_name, boundary_min_lon, boundary_max_lon, boundary_min_lat, boundary_max_lat);
}


void Remap_grid_class::generate_voronoi_grid()
{
    bool is_global_grid;
    double *center_lon_values, *center_lat_values;
    double *vertex_lon_values, *vertex_lat_values;
    long i;
    Delaunay_Voronoi *delaunay_triangularization;

  
    EXECUTION_REPORT(REPORT_ERROR, -1, this->get_is_sphere_grid(), "remap software error1 in generate_voronoi_grid\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, !are_floating_values_equal(NULL_COORD_VALUE, boundary_min_lon), "remap software error2 in generate_voronoi_grid\n");
    EXECUTION_REPORT(REPORT_WARNING, -1, boundary_min_lon != NULL_COORD_VALUE, "the boundary of area of grid %s (%lx) has not been set by user. Default boundary area (global area) will be used to generate the voronoi grid\n", grid_name, this);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Generate voronoi grid for \"%s\": %lf  %lf  %lf  %lf", grid_name, boundary_min_lon, boundary_max_lon, boundary_min_lat, boundary_max_lat);

    are_vertex_values_set_in_default = true;

    is_global_grid = boundary_min_lat == -90 && boundary_max_lat == 90 && fabs(boundary_min_lon-boundary_max_lon) == 360;
    center_lon_values = (double*) get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf;
    center_lat_values = (double*) get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf;

    if (!is_global_grid)
        for (i = 0; i < grid_size; i ++) {
            EXECUTION_REPORT(REPORT_ERROR, -1, center_lat_values[i] >= boundary_min_lat && center_lat_values[i] <= boundary_max_lat, 
                             "points in sphere grid %s are out of the latitude boundary. Please enlarge the region of latitude", this->grid_name);
            if (boundary_min_lon < boundary_max_lon)
                EXECUTION_REPORT(REPORT_ERROR, -1, center_lon_values[i] >= boundary_min_lon && center_lon_values[i] <= boundary_max_lon, 
                                 "points in sphere grid %s are out of the longitude boundary. Please enlarge the region of longitude", this->grid_name);
            else EXECUTION_REPORT(REPORT_ERROR, -1, center_lon_values[i] <= boundary_min_lon || center_lon_values[i] >= boundary_max_lon, 
                                  "points in sphere grid %s are out of the longitude boundary. Please enlarge the region of longitude", this->grid_name);
        }

    delaunay_triangularization = new Delaunay_Voronoi(grid_size, center_lat_values, center_lon_values, is_global_grid,
                                               boundary_min_lon, boundary_max_lon, boundary_min_lat, boundary_max_lat,
                                               redundant_cell_mark, &vertex_lon_values, &vertex_lat_values, &num_vertexes);
    delete delaunay_triangularization;

    EXECUTION_REPORT(REPORT_ERROR, -1, grid_vertex_fields.size() == 0, "remap software error2 in generate_voronoi_grid\n");
    grid_vertex_fields.push_back(get_grid_center_field(COORD_LABEL_LON)->duplicate_grid_data_field(this, num_vertexes, false, false));
    grid_vertex_fields.push_back(get_grid_center_field(COORD_LABEL_LAT)->duplicate_grid_data_field(this, num_vertexes, false, false));
    delete [] grid_vertex_fields[0]->grid_data_field->data_buf;
    delete [] grid_vertex_fields[1]->grid_data_field->data_buf;
    grid_vertex_fields[0]->grid_data_field->data_buf = vertex_lon_values;
    grid_vertex_fields[1]->grid_data_field->data_buf = vertex_lat_values;
}


void Remap_grid_class::set_coord_vertex_values_in_default()
{
    bool is_unit_degree_dim1, is_unit_degree_dim2;
    double *center_values, *vertex_values;
    long i, j, k;
    long cell_id, group_grid_cell_index[256];
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256], *current_leaf_grid, *leaf_grid_dim1, *leaf_grid_dim2;


    if (grid_center_fields.size() == 0)
        return;

    if (grid_vertex_fields.size() > 0)
        return;

    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++) {
        if (leaf_grids[i]->super_grid_of_setting_coord_values == this && leaf_grids[i]->grid_size == 0)
            return;
    }

    are_vertex_values_set_in_default = true;
    num_vertexes = 1;
    for (i = 0; i < grid_center_fields.size(); i ++)
        num_vertexes *= 2;
    for (i = 0; i < grid_center_fields.size(); i ++) 
        grid_vertex_fields.push_back(grid_center_fields[i]->duplicate_grid_data_field(this, num_vertexes, false, false));

    if (grid_center_fields.size() == 1) {
        for (j = 0; j < num_leaf_grids; j ++)
            if (is_sigma_grid()) {
                if (leaf_grids[j]->has_grid_coord_label(COORD_LABEL_LEV)){
                    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(grid_center_fields[0]->get_grid_data_field()->field_name_in_application, COORD_LABEL_LEV), "remap software error0 in set_coord_vertex_values_in_default\n");
                    break;
                }
            }
            else {
                if (leaf_grids[j]->get_grid_center_field() == grid_center_fields[0])
                    break;
            }
        EXECUTION_REPORT(REPORT_ERROR, -1, j < num_leaf_grids, "remap software error1 in set_coord_vertex_values_in_default\n");
        current_leaf_grid = leaf_grids[j];
        EXECUTION_REPORT(REPORT_ERROR, -1, current_leaf_grid->grid_size > 0, "remap software error2 in set_coord_vertex_values_in_default\n");
        grid_center_fields[0]->interchange_grid_data(current_leaf_grid);
        grid_vertex_fields[0]->interchange_grid_data(current_leaf_grid);
        is_unit_degree_dim1 = words_are_the_same(current_leaf_grid->coord_unit, COORD_UNIT_DEGREES);
        for (j = 0; j < grid_size/current_leaf_grid->grid_size; j ++) {
            center_values = (double*)grid_center_fields[0]->grid_data_field->data_buf + j*current_leaf_grid->grid_size;
            vertex_values = (double*)grid_vertex_fields[0]->grid_data_field->data_buf + 2*j*current_leaf_grid->grid_size;
            set_1D_coord_vertex_values_in_default(center_values,
                                                  vertex_values,
                                                  current_leaf_grid->grid_size,
                                                  current_leaf_grid->cyclic,
                                                  is_unit_degree_dim1);
        }        
        grid_vertex_fields[0]->interchange_grid_data(this);
        grid_center_fields[0]->interchange_grid_data(this);
    }
    else if (grid_center_fields.size() == 2) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 2, "detect unsupported grid in set_coord_vertex_values_in_default, please contact the authors of CoR\n");
        if (leaf_grids[1]->cyclic) {
            leaf_grid_dim1 = leaf_grids[1]; 
            leaf_grid_dim2 = leaf_grids[0];
        }
        else {
            leaf_grid_dim1 = leaf_grids[0];
            leaf_grid_dim2 = leaf_grids[1];
        }
        is_unit_degree_dim1 = words_are_the_same(leaf_grid_dim1->coord_unit, COORD_UNIT_DEGREES);
        is_unit_degree_dim2 = words_are_the_same(leaf_grid_dim2->coord_unit, COORD_UNIT_DEGREES);
        set_2D_coord_vertex_values_in_default((double*)leaf_grid_dim1->get_grid_center_field()->grid_data_field->data_buf,
                                              (double*)leaf_grid_dim1->get_grid_vertex_field()->grid_data_field->data_buf,
                                              leaf_grid_dim1->grid_size, 
                                              leaf_grid_dim1->cyclic,
                                              is_unit_degree_dim1,
                                              (double*)leaf_grid_dim2->get_grid_center_field()->grid_data_field->data_buf,
                                              (double*)leaf_grid_dim2->get_grid_vertex_field()->grid_data_field->data_buf,
                                              leaf_grid_dim2->grid_size, 
                                              leaf_grid_dim2->cyclic,
                                              is_unit_degree_dim2);
        return;
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "detect unsupported grid2 in set_coord_vertex_values_in_default, please contact the authors of CoR\n");
}


bool Remap_grid_class::check_vertex_fields_completeness(Remap_operator_basis *remap_operator)
{
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];    
    int i;


    if (this->grid_vertex_fields.size() == 1)
        EXECUTION_REPORT(REPORT_ERROR, -1, this->num_vertexes == 2, 
                     "as grid \"%s\" sets the vertex values of only one coordinate \"%s\", the number of vertexes of each cell must be 2\n", 
                     this->grid_name, this->grid_vertex_fields[0]->grid_data_field->field_name_in_application);    
    if (this->grid_vertex_fields.size() > 0)
        EXECUTION_REPORT(REPORT_ERROR, -1, this->grid_center_fields.size() == this->grid_vertex_fields.size(), 
                     "when vertex values of some coordinates are set in grid \"%s\", the vertex values of all coordinates whose center values are set must be given by users\n",
                     this->grid_name, this->grid_vertex_fields[0]->grid_data_field->field_name_in_application);
    
    this->get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    
    for (i = 0; i < num_leaf_grids; i ++)
         if (leaf_grids[i]->get_grid_vertex_field() == NULL)
             return false;

    return true;
}


long Remap_grid_class::get_lower_dimension_size_of_leaf_grid(Remap_grid_data_class *remap_grid_data_field, Remap_grid_class *leaf_grid) 
{
    long lower_dimension_size = 1;
    int i;


    for (i = 0; i < remap_grid_data_field->sized_grids.size(); i ++) {
        if (leaf_grid->is_subset_of_grid(remap_grid_data_field->sized_grids[i])) {
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grid == remap_grid_data_field->sized_grids[i], 
                         "remap software error1 in get_lower_dimension_size_of_leaf_grid\n");            
            return lower_dimension_size;
        }
        lower_dimension_size *= remap_grid_data_field->sized_grids[i]->get_grid_size();
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in get_lower_dimension_size_of_leaf_grid\n");    

    return 0;
}


long Remap_grid_class::get_higher_dimension_size_of_leaf_grid(Remap_grid_data_class *remap_grid_data_field, Remap_grid_class *leaf_grid) 
{
    long higher_dimension_size = 1;
    int i;


    for (i = remap_grid_data_field->sized_grids.size()-1; i >= 0; i --) {
        if (leaf_grid->is_subset_of_grid(remap_grid_data_field->sized_grids[i])) {
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grid == remap_grid_data_field->sized_grids[i], 
                         "remap software error1 in get_higher_dimension_size_of_leaf_grid\n");            
            return higher_dimension_size;
        }
        higher_dimension_size *= remap_grid_data_field->sized_grids[i]->get_grid_size();
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in get_higher_dimension_size_of_leaf_grid\n");    

    return 0;  
}


void Remap_grid_class::get_grid_index_interchange_table(Remap_grid_class *interchange_grid, int *index_interchange_table) 
{
    int i, j;
    int num_sized_sub_grids_this, num_sized_sub_grids_interchange;
    Remap_grid_class *sized_sub_grids_this[256], *sized_sub_grids_interchange[256];


    EXECUTION_REPORT(REPORT_ERROR, -1, this->is_similar_grid_with(interchange_grid),
                 "remap software error in get_grid_index_interchange_table\n");        

    this->get_sized_sub_grids(&num_sized_sub_grids_this, sized_sub_grids_this);
    interchange_grid->get_sized_sub_grids(&num_sized_sub_grids_interchange, sized_sub_grids_interchange);

    for (i = 0; i < num_sized_sub_grids_this; i ++) 
        for (j = 0; j < num_sized_sub_grids_interchange; j ++)
            if (sized_sub_grids_this[i] == sized_sub_grids_interchange[j]) {
                index_interchange_table[i] = j;
                break;
            }
}


void Remap_grid_class::add_partial_grid_area(const char *area_name)
{
    int i;
    Partial_area *partial_area;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, this->whole_grid != NULL, "%s is not a partial grid. Can not add partial grid area to it\n", this->grid_name);    
    EXECUTION_REPORT(REPORT_ERROR, -1, this->enable_to_set_coord_values, "can not set the partial area of grid \"%s\", because it has been used for remapping\n", this->grid_name);
    for (i = 0; i < partial_areas.size(); i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, !words_are_the_same(partial_areas[i]->area_name, area_name),
                     "partial area \"%s\" has been added in grid \"%s\". Can not add the same partial area more than once\n", 
                     area_name, this->grid_name);    

    partial_area = new Partial_area;
    strcpy(partial_area->area_name, area_name);
    partial_areas.push_back(partial_area);
}


void Remap_grid_class::add_partial_grid_area_bounds(const char *area_name, 
                                                    const char *area_bound_grid_name, 
                                                    const char *bound_type, 
                                                    const char *bound_start_string, 
                                                    const char *bound_end_string)
{
    int i;
    Remap_grid_class *area_bound_grid;
    Partial_area *partial_area;
    Partial_area_bound partial_area_bound;
    double bound_start, bound_end;


    sscanf(bound_start_string, "%lf", &bound_start);
    sscanf(bound_end_string, "%lf", &bound_end);
    
    partial_area = NULL;
    for (i = 0; i < partial_areas.size(); i ++)
        if (words_are_the_same(partial_areas[i]->area_name, area_name)) {
            partial_area = partial_areas[i];
            break;
        }
        
    EXECUTION_REPORT(REPORT_ERROR, -1, partial_area != NULL, "\"%s\" has not been declared as a partial area of grid \"%s\"\n", area_name, this->grid_name);

    area_bound_grid = remap_grid_manager->search_remap_grid_with_grid_name(area_bound_grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, area_bound_grid->is_subset_of_grid(this), "grid %s can not locate the area bounds becase it is not a sub grid of %s\n", area_bound_grid_name, this->grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, !area_bound_grid->is_partial_grid(), "grid %s can not locate the area bounds of %s becase it is not a whole grid\n", area_bound_grid_name, this->grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(bound_type, PARTIAL_AREA_BOUND_TYPE_INDEX) || words_are_the_same(bound_type, PARTIAL_AREA_BOUND_TYPE_VALUE),
                 "the type of partial grid area bound must be \"index\" or \"value\"\n");

    if (words_are_the_same(bound_type, PARTIAL_AREA_BOUND_TYPE_INDEX)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, area_bound_grid->grid_size > 0, "the area bound type can not be \"index\" because the size of located grid \"%s\" is unknown\n", area_bound_grid->grid_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, bound_start <= bound_end, 
                     "the two area bounds \"%s\" and \"%s\" must be in ascending order when area bound type is \"index\"\n",
                     bound_start_string, bound_end_string, area_bound_grid->grid_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, bound_start >= 0 && bound_end < area_bound_grid->grid_size, 
                     "the two area bounds \"%s\" and \"%s\" are out of grid size of \"%s\"\n",
                     bound_start_string, bound_end_string, area_bound_grid->grid_name);
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, area_bound_grid->num_dimensions == 1, 
                     "\"%s\" can not be the grid of locating the partial area when area bound type is \"value\", because it is not a 1D grid\n",
                     area_bound_grid->grid_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, area_bound_grid->super_grid_of_setting_coord_values != NULL,
                     "\"%s\" can not be the grid of locating the partial area when area bound type is \"value\", because the corresponding center coordinate values are not set\n",
                     area_bound_grid->grid_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, area_bound_grid->super_grid_of_setting_coord_values->is_subset_of_grid(this), 
                     "\"%s\" can not be the grid of locating the partial area when area bound type is \"value\", because grid \"%s\" which sets the corresponding center coordinate values is not a subset of grid \"%s\"\n",
                     area_bound_grid->grid_name, area_bound_grid->super_grid_of_setting_coord_values->grid_name, this->grid_name);
        if (!area_bound_grid->cyclic)
            EXECUTION_REPORT(REPORT_ERROR, -1, bound_start <= bound_end, 
                         "the two area bounds \"%s\" and \"%s\" must be in ascending order because area bound type is \"value\" and the located grid \"%s\" is acyclic\n",
                         bound_start_string, bound_end_string, area_bound_grid->grid_name);
    }

    for (i = 0; i < partial_area->area_bounds.size(); i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, !area_bound_grid->have_overlap_with_grid(partial_area->area_bounds[i].area_bound_grid), 
                     "\"%s\" can not be the grid of locating the partial area because it has overlap with another grid \"%s\" of locating the bounds of the same partial area\n",
                     area_bound_grid->grid_name, partial_area->area_bounds[i].area_bound_grid->grid_name);

    partial_area_bound.area_bound_grid = area_bound_grid;
    strcpy(partial_area_bound.bound_type, bound_type);
    partial_area_bound.bound_start = bound_start;
    partial_area_bound.bound_end = bound_end;
    partial_area->area_bounds.push_back(partial_area_bound);
}


void Remap_grid_class::compute_partial_grid_mask()
{
    Remap_data_field *remap_data_field;
    long i, j, k, l;
    bool *mask_values, *tmp_mask_values, *tmp_mask_values_segment;
    int num_leaf_grids, num_sized_grids;
    Remap_grid_class *leaf_grids[256], *sized_grids[256];
    Remap_grid_class *interchanged_grid;
    Remap_grid_data_class *tmp_grid_mask_field;
    Remap_grid_data_class *expanded_center_value_field;
    double *expanded_center_values;
    double bound_end, bound_start;


    if (whole_grid == NULL || grid_mask_field != NULL)
        return;
    
    EXECUTION_REPORT(REPORT_ERROR, -1, partial_areas.size() > 0, "%s is a partial grid. However, its partial grid areas are not set before using it\n", this->grid_name);
    for (i = 0; i < partial_areas.size(); i ++) 
        EXECUTION_REPORT(REPORT_ERROR, -1, partial_areas[i]->area_bounds.size() > 0, "no area bound has been set in the area %s of partial grid %s\n",
                     partial_areas[i]->area_name, this->grid_name);

    remap_data_field = new Remap_data_field;
    remap_data_field->required_data_size = grid_size;
    remap_data_field->read_data_size = grid_size;
    remap_data_field->data_buf = new bool [grid_size];
    strcpy(remap_data_field->field_name_in_application, GRID_MASK_LABEL);
    strcpy(remap_data_field->data_type_in_application, DATA_TYPE_BOOL);
    grid_mask_field = new Remap_grid_data_class(whole_grid, remap_data_field);
    tmp_grid_mask_field = grid_mask_field->duplicate_grid_data_field(whole_grid, 1, false, false);
    mask_values = (bool*) remap_data_field->data_buf;
    tmp_mask_values = (bool*) tmp_grid_mask_field->grid_data_field->data_buf;
    for (i = 0; i < whole_grid->grid_size; i ++)
        mask_values[i] = false;

    whole_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, whole_grid);

    for (i = 0; i < partial_areas.size(); i ++) {
        for (j = 0; j < grid_size; j ++)
            tmp_mask_values[j] = true;
        
        /* compute mask according to index boundary */
        for (j = 0; j < partial_areas[i]->area_bounds.size(); j ++) 
            if (words_are_the_same(partial_areas[i]->area_bounds[j].bound_type, PARTIAL_AREA_BOUND_TYPE_INDEX))    {
                bound_start = partial_areas[i]->area_bounds[j].bound_start;
                bound_end = partial_areas[i]->area_bounds[j].bound_end;
                whole_grid->get_sized_sub_grids(&num_sized_grids, sized_grids);
                for (k = 0, l = 0; k < num_sized_grids; k ++)
                    if (!sized_grids[k]->is_subset_of_grid(partial_areas[i]->area_bounds[j].area_bound_grid))
                        sized_grids[l++] = sized_grids[k];
                sized_grids[l++] = partial_areas[i]->area_bounds[j].area_bound_grid;
                interchanged_grid = new Remap_grid_class("TEMP_GRID", l, sized_grids, 0);
                tmp_grid_mask_field->interchange_grid_data(interchanged_grid);
                for (k = 0; k < partial_areas[i]->area_bounds[j].area_bound_grid->grid_size; k ++) 
                    if (k < bound_start || k > bound_end) {
                        tmp_mask_values_segment = tmp_mask_values + k*grid_size/partial_areas[i]->area_bounds[j].area_bound_grid->grid_size;
                        for (l = 0; l < grid_size/partial_areas[i]->area_bounds[j].area_bound_grid->grid_size; l ++)
                            tmp_mask_values_segment[l] = false;
                    }
                delete interchanged_grid;
            }
        tmp_grid_mask_field->interchange_grid_data(whole_grid);

        /* compute mask according to value boundary */
        for (j = 0; j < partial_areas[i]->area_bounds.size(); j ++) 
            if (words_are_the_same(partial_areas[i]->area_bounds[j].bound_type, PARTIAL_AREA_BOUND_TYPE_VALUE)) {
                bound_start = partial_areas[i]->area_bounds[j].bound_start;
                bound_end = partial_areas[i]->area_bounds[j].bound_end;
                for (k = 0; k < num_leaf_grids; k ++)
                    if (leaf_grids[k] == partial_areas[i]->area_bounds[j].area_bound_grid) {
                        expanded_center_value_field = whole_grid->expand_to_generate_full_coord_value(leaf_grids[k]->get_grid_center_field());
                        break;
                    }
                EXECUTION_REPORT(REPORT_ERROR, -1, k < num_leaf_grids, "remap software error2 in compute_partial_grid_mask\n");
                expanded_center_values = (double*) expanded_center_value_field->grid_data_field->data_buf;
                if (bound_start <= bound_end) {
                    for (k = 0; k < grid_size; k ++)
                        if (expanded_center_values[k] < bound_start || expanded_center_values[k] > bound_end)
                            tmp_mask_values[k] = false;
                }
                else {
                    for (k = 0; k < grid_size; k ++)
                        if (expanded_center_values[k] < bound_start && expanded_center_values[k] > bound_end)
                            tmp_mask_values[k] = false;
                }
                delete expanded_center_value_field;
            }

        for (j = 0; j < grid_size; j ++)
            mask_values[j] = mask_values[j] || tmp_mask_values[j];
    }

    delete tmp_grid_mask_field;
}


bool Remap_grid_class::has_partial_sub_grid() const
{
    for (int i = 0; i < sub_grids.size(); i ++) {
        if (sub_grids[i]->whole_grid != NULL)
            return true;
        if (sub_grids[i]->has_partial_sub_grid())
            return true;
    }

    return false;
}


Remap_grid_data_class *Remap_grid_class::expand_to_generate_full_coord_value(Remap_grid_data_class *existing_data_field)
{
    Remap_grid_data_class *expanded_data_field;
    Remap_grid_class *sized_grids[256];
    int num_sized_grids, i, j, num_points_per_cell;
    char *existing_data, *expanded_data;
    long existing_data_size, expanded_data_size;


    EXECUTION_REPORT(REPORT_ERROR, -1, existing_data_field->get_coord_value_grid()->is_subset_of_grid(this),
                 "remap software error1 in expand_to_generate_full_coord_value %s %s\n", 
                 existing_data_field->get_coord_value_grid()->get_grid_name(), grid_name);

    num_points_per_cell = existing_data_field->get_grid_data_field()->read_data_size / existing_data_field->get_coord_value_grid()->get_grid_size();
    
    if (existing_data_field->get_coord_value_grid()->is_similar_grid_with(this)) {
        expanded_data_field = existing_data_field->duplicate_grid_data_field(this, num_points_per_cell, true, true);
        return expanded_data_field;
    }

    expanded_data_field = existing_data_field->duplicate_grid_data_field(this, num_points_per_cell, false, false);
    expanded_data_size = expanded_data_field->grid_data_field->required_data_size * get_data_type_size(existing_data_field->grid_data_field->data_type_in_application);
    expanded_data = (char*) expanded_data_field->grid_data_field->data_buf;
    existing_data_size = existing_data_field->grid_data_field->required_data_size * get_data_type_size(existing_data_field->grid_data_field->data_type_in_application);
    existing_data = (char*) existing_data_field->grid_data_field->data_buf;
    for (i = 0; i < expanded_data_size/existing_data_size; i ++)
        memcpy(expanded_data+existing_data_size*i, existing_data, existing_data_size);

    this->get_sized_sub_grids(&num_sized_grids, sized_grids);
    expanded_data_field->sized_grids.clear();
    for (i = 0; i < existing_data_field->sized_grids.size(); i ++) {
        expanded_data_field->sized_grids.push_back(existing_data_field->sized_grids[i]);
        for (j = 0; j < num_sized_grids; j ++)
            if (existing_data_field->sized_grids[i] == sized_grids[j]) {
                sized_grids[j] = NULL;
                break;
            }
        EXECUTION_REPORT(REPORT_ERROR, -1, j < num_sized_grids, "remap software error2 in expand_to_generate_full_coord_value\n");
    }
    for (i = 0; i < num_sized_grids; i ++)
        if (sized_grids[i] != NULL)
            expanded_data_field->sized_grids.push_back(sized_grids[i]);

    expanded_data_field->interchange_grid_data(this);
    return expanded_data_field;
}


void Remap_grid_class::transform_coord_unit_from_radian_to_degrees()
{
    if (this->num_dimensions > 1 ||
        !words_are_the_same(this->coord_unit, COORD_UNIT_RADIANS) ||
        super_grid_of_setting_coord_values == NULL ||
        whole_grid != NULL)
        return;

    strcpy(coord_unit, COORD_UNIT_DEGREES);
}


void Remap_grid_class::check_coord_values_range() const
{
    double range_start, range_end;
    Remap_grid_data_class *center_field, *vertex_field;
    long i;
    char tmp_string[256];

    
    if (this->num_dimensions > 1 || !words_are_the_same(this->coord_unit, COORD_UNIT_DEGREES) ||
        super_grid_of_setting_coord_values == NULL || whole_grid != NULL)
        return;

    if (words_are_the_same(coord_label, COORD_LABEL_LON)) {
        range_start = 0;
        range_end = 360;
    }
    else if (words_are_the_same(coord_label, COORD_LABEL_LAT)){
        range_start = -90;
        range_end = 90;
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error in check_coord_values_range\n");
    
    center_field = get_grid_center_field();
    vertex_field = get_grid_vertex_field();

    sprintf(tmp_string, "%lf and %lf", range_start, range_end);
    for (i = 0; i < center_field->grid_data_field->required_data_size; i ++) {
        if (((double*) center_field->grid_data_field->data_buf)[i] == NULL_COORD_VALUE)
            continue;
        EXECUTION_REPORT(REPORT_WARNING, -1, ((double*) center_field->grid_data_field->data_buf)[i] >= range_start &&
                       ((double*) center_field->grid_data_field->data_buf)[i] <= range_end,
                       "the center values of coordinate \"%s\" of grid \"%s\" must be between %s\n",
                       coord_label, super_grid_of_setting_coord_values->grid_name, tmp_string);
        if (((double*) center_field->grid_data_field->data_buf)[i] < range_start)
            ((double*) center_field->grid_data_field->data_buf)[i] = range_start;
        if (((double*) center_field->grid_data_field->data_buf)[i] > range_end)
            ((double*) center_field->grid_data_field->data_buf)[i] = range_end;
    }
    if (vertex_field != NULL)
        for (i = 0; i < vertex_field->grid_data_field->required_data_size; i ++) 
            if (((double*) vertex_field->grid_data_field->data_buf)[i] != NULL_COORD_VALUE) {
                EXECUTION_REPORT(REPORT_WARNING, -1, ((double*) vertex_field->grid_data_field->data_buf)[i] >= range_start &&
                               ((double*) vertex_field->grid_data_field->data_buf)[i] <= range_end,
                               "the vertex values of coordinate \"%s\" of grid \"%s\" must be between %s\n",
                               coord_label, super_grid_of_setting_coord_values->grid_name, tmp_string);
                if (((double*) vertex_field->grid_data_field->data_buf)[i] < range_start)
                    ((double*) vertex_field->grid_data_field->data_buf)[i] = range_start;
                if (((double*) vertex_field->grid_data_field->data_buf)[i] > range_end)
                    ((double*) vertex_field->grid_data_field->data_buf)[i] = range_end;
            }
}


void Remap_grid_class::check_center_fields_sorting_order()
{
    long i, j, k;
    int num_sized_sub_grids;
    Remap_grid_class *sized_sub_grids[256], *sub_grids[256], *interchange_grid;
    double *current_coord_values;
    bool have_set_sorting_order, ascending_order;

    
    if (grid_center_fields.size() == 0)
        return;
    for (i = 0; i < grid_center_fields.size(); i ++) {
        if (words_are_the_same(grid_center_fields[i]->get_grid_data_field()->field_name_in_application, COORD_LABEL_LAT) || 
            words_are_the_same(grid_center_fields[i]->get_grid_data_field()->field_name_in_application, COORD_LABEL_LON))
            continue;
		if (grid_center_fields[i]->get_grid_data_field()->required_data_size != grid_size)  // pseudo 3-D vertical center field
			continue;
        get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
        for (j = 0, k = 0; j < num_sized_sub_grids; j ++)
            if (sized_sub_grids[j]->has_grid_coord_label(grid_center_fields[i]->get_grid_data_field()->field_name_in_application)) {
                EXECUTION_REPORT(REPORT_ERROR, -1, sized_sub_grids[j]->num_dimensions == 1 && sized_sub_grids[j]->grid_size > 0, "remap software error1 in check_center_fields_sorting_order");
                sub_grids[k++] = sized_sub_grids[j];
            }
        if (sub_grids[0]->grid_size == 1)
            continue;
        for (j = 0; j < num_sized_sub_grids; j ++)
            if (sized_sub_grids[j] != sub_grids[0])
                sub_grids[k++] = sized_sub_grids[j];
        interchange_grid = new Remap_grid_class("TMP_GRID", k, sub_grids, this->grid_size);
        grid_center_fields[i]->interchange_grid_data(interchange_grid);
        have_set_sorting_order = false;
        for (j = 0; j < this->grid_size/sub_grids[0]->grid_size; j ++) {
            current_coord_values = (double*)(grid_center_fields[i]->get_grid_data_field()->data_buf) + j*sub_grids[0]->grid_size;
            for (k = 1; k < sub_grids[0]->grid_size; k ++) {
                if (!have_set_sorting_order && current_coord_values[k-1] != current_coord_values[k]) {
                    have_set_sorting_order = true;
                    ascending_order = current_coord_values[k-1] < current_coord_values[k];
                }
                EXECUTION_REPORT(REPORT_ERROR, -1, current_coord_values[k-1] == current_coord_values[k] || ascending_order == current_coord_values[k-1] < current_coord_values[k], "the center coordinate values of subgrid %s of grid %s are not sorted in ascending or descending order\n", sub_grids[0]->get_grid_name(), this->grid_name);
            }
        }
        delete interchange_grid;
    }
}


void Remap_grid_class::check_center_vertex_values_consistency_2D()
{
    Remap_grid_data_class *center_field1, *center_field2, *vertex_field1, *vertex_field2;
    long i, j;
    double lat_rotated, lon_rotated;
    double center_coord1_value, center_coord2_value, vertex_coord1_values[65536], vertex_coord2_values[65536];


    EXECUTION_REPORT(REPORT_ERROR, -1, num_vertexes <= 65536, "Software error in Remap_grid_class::check_center_vertex_values_consistency_2D: too big num_vertexes of the grid", num_vertexes);

    if (are_vertex_values_set_in_default || grid_center_fields.size() != 2 || grid_vertex_fields.size() == 0)
        return;
    
    EXECUTION_REPORT(REPORT_ERROR, -1, grid_center_fields.size() == grid_vertex_fields.size(), "remap software error1 in check_center_vertex_values_relation_2D\n");
    for (i = 0; i < grid_center_fields.size(); i ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_center_fields[i]->grid_data_field->required_data_size == grid_size, "remap software error2 in check_center_vertex_values_relation_2D\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_vertex_fields[i]->grid_data_field->required_data_size == grid_size*num_vertexes, "remap software error3 in check_center_vertex_values_relation_2D\n");
    }

    center_field1 = grid_center_fields[0];
    center_field2 = grid_center_fields[1];
    vertex_field1 = grid_vertex_fields[0];
    vertex_field2 = grid_vertex_fields[1];

    if (get_is_sphere_grid())
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(center_field1->get_grid_data_field()->field_name_in_application, COORD_LABEL_LON) &&
                     words_are_the_same(center_field2->get_grid_data_field()->field_name_in_application, COORD_LABEL_LAT) &&
                     words_are_the_same(vertex_field1->get_grid_data_field()->field_name_in_application, COORD_LABEL_LON) &&
                     words_are_the_same(vertex_field2->get_grid_data_field()->field_name_in_application, COORD_LABEL_LAT), 
                     "remap software error4 in check_center_vertex_values_relation_2D\n");

    for (i = 0; i < grid_size; i ++) {
        center_coord1_value = ((double*) center_field1->grid_data_field->data_buf)[i];
        center_coord2_value = ((double*) center_field2->grid_data_field->data_buf)[i];
        for (j = 0; j < num_vertexes; j ++) {
            vertex_coord1_values[j] = (((double*) vertex_field1->grid_data_field->data_buf)+i*num_vertexes)[j];
            vertex_coord2_values[j] = (((double*) vertex_field2->grid_data_field->data_buf)+i*num_vertexes)[j];
        }
        if (get_is_sphere_grid() && fabs(center_coord2_value) > SPHERE_GRID_ROTATION_LAT_THRESHOLD) {
            rotate_sphere_coordinate(center_coord1_value, center_coord2_value, center_coord1_value, center_coord2_value);
            for (j = 0; j < num_vertexes; j ++)
                if (vertex_coord1_values[j] != NULL_COORD_VALUE)
                    rotate_sphere_coordinate(vertex_coord1_values[j], vertex_coord2_values[j], vertex_coord1_values[j], vertex_coord2_values[j]);
        }
        if (!is_point_in_2D_cell(center_coord1_value,
                                 center_coord2_value,
                                 vertex_coord1_values,
                                 vertex_coord2_values,
                                 num_vertexes,
                                 center_field1->is_unit_degree(),
                                 center_field2->is_unit_degree(),
                                 get_is_sphere_grid())) {
            char error_string1[65536], error_string2[1024];
			sprintf(error_string1, "( ");
			for (j = 0; j < num_vertexes; j ++) {
				sprintf(error_string2, "(lon=%lf,lat=%lf) ", vertex_coord1_values[j], vertex_coord2_values[j]); 
				strcat(error_string1, error_string2);
			}
			strcat(error_string1,")");
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "At the %dth point of the grid \"%s\", the center coordinate values (lon=%lf,lat=%lf) are not consistent with vertex coordinate values (%s): each center point must be surrounded by the polygon determined by the corresponding vertex points.", i+1, grid_name, center_coord1_value, center_coord2_value, error_string1);
        }
    }
}


void Remap_grid_class::formalize_sphere_grid()
{
    int i, j;
    Remap_grid_data_class *tmp_field1, *tmp_field2;


    if (!get_is_sphere_grid() || grid_center_fields.size() < 2)
        return;

    /* Order the center and vertex fields of sphere grid: lon before lat */
    for (i = 0; i < 2; i ++)
        if (words_are_the_same(grid_center_fields[i]->grid_data_field->field_name_in_application, COORD_LABEL_LON))
            tmp_field1 = grid_center_fields[i];
        else tmp_field2 = grid_center_fields[i];
    grid_center_fields.clear();
    grid_center_fields.push_back(tmp_field1);
    grid_center_fields.push_back(tmp_field2);

    if (grid_vertex_fields.size() > 0) {
        for (i = 0; i < 2; i ++)
            if (words_are_the_same(grid_vertex_fields[i]->grid_data_field->field_name_in_application, COORD_LABEL_LON))
                tmp_field1 = grid_vertex_fields[i];
            else tmp_field2 = grid_vertex_fields[i];
        grid_vertex_fields.clear();
        grid_vertex_fields.push_back(tmp_field1);
        grid_vertex_fields.push_back(tmp_field2);
    }
}


void Remap_grid_class::end_grid_definition_stage(Remap_operator_basis *remap_operator)
{
    Remap_grid_class *similar_grid;


    EXECUTION_REPORT(REPORT_ERROR, -1, grid_size > 0, "the grid size of \"%s\" is unknown. This grid can not be used to generator field data or remap operator\n", grid_name);

    for (int i = 0; i < remap_grid_manager->remap_grids.size(); i ++)
        if (remap_grid_manager->remap_grids[i]->is_subset_of_grid(this)) {
            if (!remap_grid_manager->remap_grids[i]->enable_to_set_coord_values)
                continue;
            remap_grid_manager->remap_grids[i]->formalize_sphere_grid();
            remap_grid_manager->remap_grids[i]->compute_partial_grid_mask();
            remap_grid_manager->remap_grids[i]->transform_coord_unit_from_radian_to_degrees();
            remap_grid_manager->remap_grids[i]->detect_redundant_cells();
            remap_grid_manager->remap_grids[i]->set_coord_vertex_values_in_default();
            remap_grid_manager->remap_grids[i]->remove_redundant_vertex_values();
            remap_grid_manager->remap_grids[i]->check_center_vertex_values_consistency_2D();
            remap_grid_manager->remap_grids[i]->check_coord_values_range();
            remap_grid_manager->remap_grids[i]->check_center_fields_sorting_order();
            remap_grid_manager->remap_grids[i]->enable_to_set_coord_values = false;
        }

    similar_grid = get_similar_grids_setting_coord_values();
    if (similar_grid != NULL)
        this->num_vertexes = similar_grid->num_vertexes;

    this->calculate_area_or_volumn();
}


bool Remap_grid_class::is_partial_grid() const
{
    return whole_grid != NULL || has_partial_sub_grid();
}


void Remap_grid_class::compute_remap_field_data_runtime_mask(Remap_grid_class *field_data_partial_grid, 
                                                             Remap_grid_class **mask_sub_grids,
                                                             int *num_mask_sub_grids,
                                                             Remap_grid_data_class **runtime_mask)
{
    int temp_num_mask_sub_grids;
    Remap_grid_data_class *sub_mask_fields[256], *mask_field_expanded;
    bool *runtime_mask_values, *mask_values_expanded;
    long i, j;


    this->get_masked_sub_grids(num_mask_sub_grids, mask_sub_grids);
    for (i = 0; i < *num_mask_sub_grids; i ++)
        sub_mask_fields[i] = mask_sub_grids[i]->grid_mask_field;
    
    if (field_data_partial_grid != NULL) {
        field_data_partial_grid->get_partial_grid_mask_fields(&temp_num_mask_sub_grids, sub_mask_fields+*num_mask_sub_grids, field_data_partial_grid);
        for (i = 0; i < temp_num_mask_sub_grids; i ++)
            mask_sub_grids[i+(*num_mask_sub_grids)] = (Remap_grid_class*) sub_mask_fields[i+(*num_mask_sub_grids)]->get_coord_value_grid();
        (*num_mask_sub_grids) += temp_num_mask_sub_grids;
    }

    if ((*num_mask_sub_grids) == 0) {
        *runtime_mask = NULL;
        return;
    }

    if ((*num_mask_sub_grids) == 1) {
        *runtime_mask = sub_mask_fields[0]->duplicate_grid_data_field(sub_mask_fields[0]->get_coord_value_grid(), 1, true, true);
        return;
    }

    *runtime_mask = sub_mask_fields[0]->duplicate_grid_data_field(this, 1, false, false);
    runtime_mask_values = (bool*) (*runtime_mask)->grid_data_field->data_buf;
    for (i = 0; i < this->grid_size; i ++)
        runtime_mask_values[i] = true;

    for (i = 0; i < *num_mask_sub_grids; i ++) {
        mask_field_expanded = this->expand_to_generate_full_coord_value(sub_mask_fields[i]);
        mask_values_expanded = (bool*) mask_field_expanded->grid_data_field->data_buf;
        for (j = 0; j < this->grid_size; j ++)
            runtime_mask_values[j] = runtime_mask_values[j] && mask_values_expanded[j];
        delete mask_field_expanded;
    }
}


void Remap_grid_class::detect_redundant_cells()
{
    double *full_center_coord_values[256];
    Remap_grid_data_class *full_center_fields[256];
    int num_defined_center_fields, num_leaf_grids;
    Remap_grid_class *leaf_grids[256];
    long *cell_index, i, next_i, j;
    Radix_sort<double, long> *radix_sort;
    Remap_data_field *mask_data_field;
	bool has_redundant_cells = false;


    if (grid_center_fields.size() == 0)
        return;

    mask_data_field = new Remap_data_field;
    strcpy(mask_data_field->data_type_in_application, DATA_TYPE_BOOL);
    mask_data_field->required_data_size = grid_size;    
    mask_data_field->read_data_size = grid_size;
    mask_data_field->data_buf = new char [mask_data_field->required_data_size];
    redundant_cell_mark_field = new Remap_grid_data_class(this, mask_data_field);
    redundant_cell_mark = (bool*) mask_data_field->data_buf;
    for (i = 0; i < grid_size; i ++)
        redundant_cell_mark[i] = false;

    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0, num_defined_center_fields = 0; i < num_leaf_grids; i ++) {
        full_center_fields[i] = leaf_grids[i]->get_grid_center_field();
        if (full_center_fields[i] != NULL) 
            num_defined_center_fields ++;
    }

    if (num_defined_center_fields != num_dimensions) 
        return;

    for (i = 0; i < num_dimensions; i ++) {
        full_center_fields[i] = this->expand_to_generate_full_coord_value(full_center_fields[i]);
        full_center_coord_values[i] = (double*) full_center_fields[i]->grid_data_field->data_buf;
    }

    cell_index = new long [grid_size];
    for (i = 0; i < grid_size; i ++)
        cell_index[i] = i;
    
    radix_sort = new Radix_sort<double, long>(full_center_coord_values, 
                                              num_dimensions, 
                                              cell_index,
                                              grid_size,
                                              TOLERABLE_ERROR);
    radix_sort->do_radix_sort();    

    for (i = 0; i < grid_size; i ++) {
        next_i = (i+1) % grid_size;
        for (j = 0; j < num_dimensions; j ++)
            if (fabs(radix_sort->radix_values[j][i]-radix_sort->radix_values[j][next_i]) > TOLERABLE_ERROR)
                break;
        if (j == num_dimensions) {
            redundant_cell_mark[radix_sort->content[next_i]] = true;
			has_redundant_cells = true;
        }
    }

    delete radix_sort;
    delete [] cell_index;    

    if (get_is_sphere_grid()) {
        int north_pole_cell_index = -1, south_pole_cell_index = -1;
        for (i = 0; i < grid_size; i ++) {
            if (((float)full_center_coord_values[1][i]) == (float)90.0)
                if (north_pole_cell_index == -1)
                    north_pole_cell_index = i;
                else {
					redundant_cell_mark[i] = true;
					has_redundant_cells = true;
                }
            if (((float)full_center_coord_values[1][i]) == (float)(-90.0))
                if (south_pole_cell_index == -1)
                    south_pole_cell_index = i;
                else {
					redundant_cell_mark[i] = true;
					has_redundant_cells = true;
                }
        }
    }

    for (i = 0; i < num_dimensions; i ++)
        delete full_center_fields[i];

	if (!has_redundant_cells) {
		delete redundant_cell_mark_field;
		redundant_cell_mark_field = NULL;
		redundant_cell_mark = NULL;
	}
}


void Remap_grid_class::remove_redundant_vertex_values()
{
    long i, j, k, m, j_1;
    double *vertex_values[256], value1, value2;
    bool is_lon_vertex[256];


    if (grid_vertex_fields.size() > 0) {
        for (i = 0; i < grid_vertex_fields.size(); i ++) {
            vertex_values[i] = (double *) grid_vertex_fields[i]->grid_data_field->data_buf;
            is_lon_vertex[i] = words_are_the_same(grid_vertex_fields[i]->grid_data_field->field_name_in_application, COORD_LABEL_LON);
        }
        for (i = 0; i < grid_size; i ++) {
            for (j = 0; j < grid_vertex_fields.size(); j ++)
                if (((double*)grid_center_fields[j]->grid_data_field->data_buf)[i] == NULL_COORD_VALUE)
                    break;
            if (j < grid_vertex_fields.size())
                for (j = 0; j < grid_vertex_fields.size(); j ++)
                    for (k = 0; k < num_vertexes; k ++)
                        vertex_values[j][i*num_vertexes+k] = NULL_COORD_VALUE;
        }
        for (i = 0; i < grid_size; i ++) {
            for (j = 0; j < num_vertexes; j ++) {
                j_1 = (j+1) % num_vertexes;
                for (k = 0; k < grid_vertex_fields.size(); k ++) {
                    value1 = vertex_values[k][i*num_vertexes+j];
                    value2 = vertex_values[k][i*num_vertexes+j_1];
                    if (is_lon_vertex[k])
                        match_degree_values(value1, value2);
                    if (fabs(value1-value2) > TOLERABLE_ERROR)
                        break;
                }
                if (k == grid_vertex_fields.size())
                    for (k = 0; k < grid_vertex_fields.size(); k ++)
                        vertex_values[k][i*num_vertexes+j] = NULL_COORD_VALUE;
            }
            for (j = 0, m = 0; j < num_vertexes; j ++) {
                if (vertex_values[0][i*num_vertexes+j] != NULL_COORD_VALUE) {
                    for (k = 0; k < grid_vertex_fields.size(); k ++)
                        vertex_values[k][i*num_vertexes+m] = vertex_values[k][i*num_vertexes+j];
                    m ++;
                }
            }
            if (num_vertexes == 2 && m == 1 || num_vertexes > 2 && m <= 2)
                m = 0;
            for (; m < num_vertexes; m ++)
                for (k = 0; k < grid_vertex_fields.size(); k ++)
                    vertex_values[k][i*num_vertexes+m] = NULL_COORD_VALUE;
        }        
    }
}


void Remap_grid_class::transform_coord_values_from_radian_to_degrees(Remap_grid_data_class *remap_grid_data)
{
    long i, array_size;
    int num_leaf_grids;
    Remap_grid_class *leaf_grids[256];
    double *coord_value_array;

    
    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->coord_label, remap_grid_data->grid_data_field->field_name_in_application))
            break;
    if (!words_are_the_same(leaf_grids[i]->coord_unit, COORD_UNIT_RADIANS))
        return;

    strcpy(leaf_grids[i]->coord_unit, COORD_UNIT_DEGREES);
    coord_value_array = (double*) remap_grid_data->grid_data_field->data_buf;
    array_size = remap_grid_data->grid_data_field->required_data_size;    
    for (i = 0; i < array_size; i ++)
        if (coord_value_array[i] != NULL_COORD_VALUE) {
            coord_value_array[i] = RADIAN_TO_DEGREE(coord_value_array[i]);
            if (are_floating_values_equal((double)90, coord_value_array[i]))
                coord_value_array[i] = 90;
            if (are_floating_values_equal((double)-90, coord_value_array[i]))
                coord_value_array[i] = -90;
        }
    for (i = 0; i < remap_grid_data->grid_data_field->field_attributes.size(); i ++)
        if (words_are_the_same(remap_grid_data->grid_data_field->field_attributes[i].attribute_name, GRID_FIELD_ATTRIBUTE_UNIT) ||
            words_are_the_same(remap_grid_data->grid_data_field->field_attributes[i].attribute_value, COORD_UNIT_RADIANS))
            strcpy(remap_grid_data->grid_data_field->field_attributes[i].attribute_value, COORD_UNIT_DEGREES);
}


void Remap_grid_class::formalize_cyclic_coord_values(Remap_grid_data_class *remap_grid_data)
{
    double coord_value_limit;
    long i, array_size;
    int num_leaf_grids;
    double *coord_value_array;
    Remap_grid_class *leaf_grids[256];


    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->coord_label, remap_grid_data->grid_data_field->field_name_in_application))
            break;
    if (!leaf_grids[i]->cyclic)
        return;

    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(leaf_grids[i]->coord_unit, COORD_UNIT_DEGREES), 
                 "remap software error1 in formalize_cyclic_coord_values\n");

    coord_value_limit = 360;
    coord_value_array = (double*) remap_grid_data->grid_data_field->data_buf;
    array_size = remap_grid_data->grid_data_field->required_data_size;

    for (i = 0; i < array_size; i ++) {
        if (coord_value_array[i] == NULL_COORD_VALUE)
            continue;
        while(coord_value_array[i] < 0)
            coord_value_array[i] += coord_value_limit;
        while(coord_value_array[i] >= coord_value_limit)
            coord_value_array[i] -= coord_value_limit;
        EXECUTION_REPORT(REPORT_ERROR, -1, coord_value_array[i] >= 0 && coord_value_array[i] < coord_value_limit,
                     "remap software error2 in formalize_cyclic_coord_values\n");
    }
}


bool Remap_grid_class::are_all_vertex_fields_specified_by_user()
{
    Remap_grid_class *leaf_grids[256];
    int num_leaf_grids, i;

    this->get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < num_leaf_grids; i ++) {
        if (leaf_grids[i]->get_grid_vertex_field() == NULL)
            return false;
        if (leaf_grids[i]->get_super_grid_of_setting_coord_values()->are_vertex_values_set_in_default)
            return false;
    }

    return true;
}


void Remap_grid_class::calculate_area_of_sphere_grid() 
{
    Remap_grid_data_class *lon_vertex_field, *lat_vertex_field;
    Remap_grid_class *leaf_grids[256];
    int tmp_num_vertex, current_num_vertex, num_leaf_grids;
    double *lon_vertex_values, *lat_vertex_values;
    double current_lon_vertex_values[65536], current_lat_vertex_values[65536], value1, value2;
    long i, j, k;
    bool should_rotate;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_dimensions == 2 && this->has_grid_coord_label(COORD_LABEL_LON) && this->has_grid_coord_label(COORD_LABEL_LAT), "remap software error1 in calculate_area_of_sphere_grid\n");
    
    if (!this->are_all_vertex_fields_specified_by_user())
        return;

    this->get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (i = 0; i < this->num_dimensions; i ++)
        if (leaf_grids[i]->has_grid_coord_label(COORD_LABEL_LON))
            lon_vertex_field = this->expand_to_generate_full_coord_value(leaf_grids[i]->get_grid_vertex_field());
        else lat_vertex_field = this->expand_to_generate_full_coord_value(leaf_grids[i]->get_grid_vertex_field());
    lon_vertex_values = (double*) lon_vertex_field->get_grid_data_field()->data_buf;
    lat_vertex_values = (double*) lat_vertex_field->get_grid_data_field()->data_buf;
    
    tmp_num_vertex = lon_vertex_field->get_grid_data_field()->read_data_size / this->grid_size;
    EXECUTION_REPORT(REPORT_ERROR, -1, tmp_num_vertex == lat_vertex_field->get_grid_data_field()->read_data_size / this->grid_size && tmp_num_vertex >= 2 && tmp_num_vertex <= 65536,
                 "remap_software error2 in calculate_area_of_sphere_grid\n");

    area_or_volumn = new double [this->grid_size];
    for (i = 0; i < this->grid_size; i ++)
        area_or_volumn[i] = 0.0;
    
    for (i = 0; i < this->grid_size; i ++) {
        if (tmp_num_vertex > 2) {
            should_rotate = false;
            current_num_vertex = 0;
            for (j = 0; j < tmp_num_vertex; j ++) {
                if (lon_vertex_values[i*tmp_num_vertex+j] == NULL_COORD_VALUE)
                    continue;
                current_lon_vertex_values[current_num_vertex] = lon_vertex_values[i*tmp_num_vertex+j];
                current_lat_vertex_values[current_num_vertex] = lat_vertex_values[i*tmp_num_vertex+j];
                if (fabs(current_lat_vertex_values[current_num_vertex]) > SPHERE_GRID_ROTATION_LAT_THRESHOLD)
                    should_rotate = true;
                current_num_vertex ++;
            }
        }
        else {    
            EXECUTION_REPORT(REPORT_ERROR, -1, lon_vertex_values[i*2] != NULL_COORD_VALUE && lon_vertex_values[i*2+1] != NULL_COORD_VALUE,
                         "remap software error3 in calculate_area_of_sphere_grid\n");
            EXECUTION_REPORT(REPORT_ERROR, -1, lat_vertex_values[i*2] != NULL_COORD_VALUE && lat_vertex_values[i*2+1] != NULL_COORD_VALUE,
                         "remap software error4 in calculate_area_of_sphere_grid : %d\n", i);
            current_num_vertex = 4;
            current_lon_vertex_values[0] = lon_vertex_values[i*2];
            current_lat_vertex_values[0] = lat_vertex_values[i*2];
            current_lon_vertex_values[1] = lon_vertex_values[i*2];
            current_lat_vertex_values[1] = lat_vertex_values[i*2+1];
            current_lon_vertex_values[2] = lon_vertex_values[i*2+1];
            current_lat_vertex_values[2] = lat_vertex_values[i*2+1];
            current_lon_vertex_values[3] = lon_vertex_values[i*2+1];
            current_lat_vertex_values[3] = lat_vertex_values[i*2];
            should_rotate = fabs(lat_vertex_values[i*2]) > SPHERE_GRID_ROTATION_LAT_THRESHOLD || fabs(lat_vertex_values[i*2+1]) > SPHERE_GRID_ROTATION_LAT_THRESHOLD;
        }
        if (should_rotate)
            for (j = 0; j < current_num_vertex; j ++)
                rotate_sphere_coordinate(current_lon_vertex_values[j], current_lat_vertex_values[j], current_lon_vertex_values[j], current_lat_vertex_values[j]);
        sort_vertexes_of_sphere_cell(current_num_vertex, current_lon_vertex_values, current_lat_vertex_values);
        for (j = current_num_vertex - 1; j >= 0 ; j --) {
            value1 = current_lon_vertex_values[j];
            value2 = current_lon_vertex_values[(j+current_num_vertex-1)%current_num_vertex];
            match_degree_values(value1, value2);
            if (fabs(value1-value2) <= TOLERABLE_ERROR && fabs(current_lat_vertex_values[j]-current_lat_vertex_values[(j+current_num_vertex-1)%current_num_vertex]) <= TOLERABLE_ERROR) {
                current_lon_vertex_values[j] = NULL_COORD_VALUE;
                current_lat_vertex_values[j] = NULL_COORD_VALUE;
            }
        }
        k = current_num_vertex;
        current_num_vertex = 0;
        for (j = 0; j < k; j ++) 
            if (current_lon_vertex_values[j] != NULL_COORD_VALUE) {
                current_lon_vertex_values[current_num_vertex] = current_lon_vertex_values[j];
                current_lat_vertex_values[current_num_vertex] = current_lat_vertex_values[j];
                current_num_vertex ++;
            }
        if (current_num_vertex > 0)
            area_or_volumn[i] = compute_area_of_sphere_cell(current_num_vertex, current_lon_vertex_values, current_lat_vertex_values);
    }

    delete lon_vertex_field;
    delete lat_vertex_field;
}


void Remap_grid_class::calculate_area_or_volumn()
{
    if (this->area_or_volumn != NULL)
        return;
    
    if (this->has_grid_coord_label(COORD_LABEL_LON) && this->has_grid_coord_label(COORD_LABEL_LAT))
        if (this->num_dimensions == 2)
            calculate_area_of_sphere_grid();
        else if (this->num_dimensions == 3 && this->has_grid_coord_label(COORD_LABEL_LEV)) {
        }
}


bool Remap_grid_class::get_is_sphere_grid() const
{
    return this->num_dimensions == 2 && this->has_grid_coord_label(COORD_LABEL_LON) && this->has_grid_coord_label(COORD_LABEL_LAT);
}


void Remap_grid_class::set_decomp_name(const char *decomp_name)
{
    strcpy(this->decomp_name, decomp_name);
}


Remap_grid_class *Remap_grid_class::generate_decomp_grid(const int *local_cell_indexes, int num_local_cells, const char *decomp_name)
{
    Remap_grid_class *leaf_grids[256];
    Remap_grid_class *decomp_leaf_grids[256];
    Remap_grid_class *decomp_grid;
    Remap_grid_data_class *original_center_field_lon, *original_center_field_lat;
    char decomp_grid_name[NAME_STR_SIZE*2];
    int num_leaf_grids;
    double *decomp_lon_center_values, *decomp_lat_center_values;
    double *this_lon_center_values, *this_lat_center_values;
    bool *decomp_mask_values, *this_mask_values;
    int i, j;


    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == 2 && leaf_grids[0]->has_grid_coord_label(COORD_LABEL_LON) && leaf_grids[1]->has_grid_coord_label(COORD_LABEL_LAT),
                 "the grid %s to decomp grid must be a lat-lon grid\n", grid_name);
    end_grid_definition_stage(NULL);

    if (leaf_grids[0]->cyclic)
        decomp_leaf_grids[0] = new Remap_grid_class(leaf_grids[0]->grid_name, leaf_grids[0]->coord_label, leaf_grids[0]->coord_unit, COORD_BOUND_CYCLIC, 0);
    else decomp_leaf_grids[0] = new Remap_grid_class(leaf_grids[0]->grid_name, leaf_grids[0]->coord_label, leaf_grids[0]->coord_unit, COORD_BOUND_ACYCLIC, 0);
    decomp_leaf_grids[1] = new Remap_grid_class(leaf_grids[1]->grid_name, leaf_grids[1]->coord_label, leaf_grids[1]->coord_unit, COORD_BOUND_ACYCLIC, 0);
    sprintf(decomp_grid_name, "%s_at_DECOMP_%s", grid_name, decomp_name);
    decomp_grid = new Remap_grid_class(decomp_grid_name, 2, decomp_leaf_grids, num_local_cells);
    decomp_leaf_grids[0]->super_grid_of_setting_coord_values = decomp_grid;
    decomp_leaf_grids[1]->super_grid_of_setting_coord_values = decomp_grid;
    decomp_leaf_grids[0]->super_grids_of_setting_mask_value.push_back(decomp_grid);
    decomp_leaf_grids[1]->super_grids_of_setting_mask_value.push_back(decomp_grid);

    EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[0]->get_grid_center_field() != NULL && leaf_grids[1]->get_grid_center_field() != NULL,
                 "user must specify the center coordinate of the horizontal 2D grid %s\n", this->grid_name);

    decomp_grid->grid_center_fields.push_back(leaf_grids[0]->get_grid_center_field()->duplicate_grid_data_field(decomp_grid, 1, false, false));
    decomp_grid->grid_center_fields.push_back(leaf_grids[1]->get_grid_center_field()->duplicate_grid_data_field(decomp_grid, 1, false, false));

    if (this->grid_mask_field == NULL)
        this->generate_default_grid_mask();
    decomp_grid->grid_mask_field = this->grid_mask_field->duplicate_grid_data_field(decomp_grid, 1, false, false);

    if (area_or_volumn != NULL)
        decomp_grid->area_or_volumn = new double [num_local_cells];
    if (leaf_grids[0]->get_grid_center_field()->get_coord_value_grid()->get_num_dimensions() == 2) {
        original_center_field_lon = leaf_grids[0]->get_grid_center_field();
        original_center_field_lat = leaf_grids[1]->get_grid_center_field();
    }
    else {
        original_center_field_lon = expand_to_generate_full_coord_value(leaf_grids[0]->get_grid_center_field());
        original_center_field_lat = expand_to_generate_full_coord_value(leaf_grids[1]->get_grid_center_field());
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, original_center_field_lon->get_coord_value_grid()->get_num_dimensions() == 2, "C-Coupler error in generate_decomp_grid\n");
    decomp_lon_center_values = (double*) decomp_grid->grid_center_fields[0]->get_grid_data_field()->data_buf;
    decomp_lat_center_values = (double*) decomp_grid->grid_center_fields[1]->get_grid_data_field()->data_buf;
    decomp_mask_values = (bool*) decomp_grid->grid_mask_field->get_grid_data_field()->data_buf;
    this_lon_center_values = (double*) original_center_field_lon->get_grid_data_field()->data_buf;
    this_lat_center_values = (double*) original_center_field_lat->get_grid_data_field()->data_buf;
    this_mask_values = (bool*) this->grid_mask_field->get_grid_data_field()->data_buf;
    for (i = 0; i < num_local_cells; i ++) {
        if (local_cell_indexes[i] != CCPL_NULL_INT) {
            decomp_lon_center_values[i] = this_lon_center_values[local_cell_indexes[i]];
            decomp_lat_center_values[i] = this_lat_center_values[local_cell_indexes[i]];
            decomp_mask_values[i] = this_mask_values[local_cell_indexes[i]];
            if (area_or_volumn != NULL)
                decomp_grid->area_or_volumn[i] = this->area_or_volumn[local_cell_indexes[i]];
        }
        else {
            decomp_lon_center_values[i] = NULL_COORD_VALUE;
            decomp_lat_center_values[i] = NULL_COORD_VALUE;
            decomp_mask_values[i] = false;
            if (area_or_volumn != NULL)
                decomp_grid->area_or_volumn[i] = NULL_COORD_VALUE;
        }
    }

    if (leaf_grids[0]->get_grid_center_field()->get_coord_value_grid()->get_num_dimensions() == 1) {
        delete original_center_field_lon;
        delete original_center_field_lat;
    }

    strcpy(decomp_grid->decomp_name, decomp_name);
    decomp_grid->original_grid = this;

    remap_grid_manager->add_temp_grid(decomp_leaf_grids[0]);
    remap_grid_manager->add_temp_grid(decomp_leaf_grids[1]);
    
    return decomp_grid;
}


void Remap_grid_class::generate_3D_grid_decomp_sigma_values(Remap_grid_class *original_3D_grid,Remap_grid_class *decomp_grid, const int *local_cell_indexes, int num_local_cells)
{
    Remap_grid_class *leaf_grids[256];
    int num_leaf_grids;
    double *local_sigma_grid_surface_values;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, this->has_grid_coord_label(COORD_LABEL_LON) && this->has_grid_coord_label(COORD_LABEL_LAT) && this->has_grid_coord_label(COORD_LABEL_LEV),
                     "C-Coupler error1 in generate_3D_grid_decomp_sigma_values\n");

    if (!original_3D_grid->is_sigma_grid())
        return;

    if (original_3D_grid->grid_center_fields.size() == 0)
        return;

    EXECUTION_REPORT(REPORT_ERROR, -1, original_3D_grid->grid_center_fields.size() == 1 && this->grid_center_fields.size() == 0, "C-Coupler error2 in generate_3D_grid_decomp_sigma_values\n");
    
    grid_center_fields.push_back(original_3D_grid->grid_center_fields[0]->duplicate_grid_data_field(this, 1, false, false));
    level_V3D_coord_trigger_field = original_3D_grid->get_a_leaf_grid_of_sigma_or_hybrid()->get_sigma_grid_sigma_value_field()->duplicate_grid_data_field(decomp_grid, 1, false, false);
    local_sigma_grid_surface_values = (double*) level_V3D_coord_trigger_field->get_grid_data_field()->data_buf;
    for (int i = 0; i < num_local_cells; i ++)
        local_sigma_grid_surface_values[i] = 0.0;

    level_V3D_coord_trigger_field_specified = original_3D_grid->level_V3D_coord_trigger_field_specified;
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "generate decomp sigma values for 3D grid %s: %ld %ld %ld\n", grid_name, get_a_leaf_grid_of_sigma_or_hybrid()->get_sigma_grid_sigma_value_field()->get_coord_value_grid()->get_grid_size(),
                     level_V3D_coord_trigger_field->get_coord_value_grid()->get_grid_size(), this->get_grid_size());

    this->get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (int i = 0; i < num_leaf_grids; i ++)
        if (leaf_grids[i]->has_grid_coord_label(COORD_LABEL_LEV)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->grid_size > 0 && leaf_grids[i]->grid_center_fields.size() == 0, "C-Coupler error3 in generate_3D_grid_decomp_sigma_values\n");            
            leaf_grids[i]->super_grid_of_setting_coord_values = this;
        }
        
    for (int i = 0; i < num_local_cells; i ++)
        if (local_cell_indexes[i] != CCPL_NULL_INT)
            local_sigma_grid_surface_values[i] = ((double*)(original_3D_grid->level_V3D_coord_trigger_field->get_grid_data_field()->data_buf))[local_cell_indexes[i]];
        else local_sigma_grid_surface_values[i] = NULL_COORD_VALUE;
}


void Remap_grid_class::renew_lev_grid_coord_values(double *new_center_coord_values, double *new_vertex_coord_values)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, grid_center_fields.size() == 1, "C-Coupler error in Remap_grid_class::renew_lev_grid_coord_value: NULL grid_center_fields");
    for (int i = 0; i < grid_size; i ++)
        ((double*)grid_center_fields[0]->get_grid_data_field()->data_buf)[i] = new_center_coord_values[i];
    if (new_vertex_coord_values != NULL)
        for (int i = 0; i < grid_size*2; i ++)
            ((double*)grid_vertex_fields[0]->get_grid_data_field()->data_buf)[i] = new_vertex_coord_values[i];
}


bool Remap_grid_class::check_coord_values_consistency(const char *coord_label, const char *data_type, const void *coord_values_model)
{
    Remap_grid_class *leaf_grids[256];
    int num_leaf_grids;
    double *coord_values_this;


    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(coord_label, COORD_LABEL_LON) || words_are_the_same(coord_label, COORD_LABEL_LAT) || words_are_the_same(coord_label, COORD_LABEL_LEV),
                 "C-Coupler error1 in check_coord_values_consistency\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, this->has_grid_coord_label(coord_label), "grid %s does not have coordinate %s\n", grid_name, coord_label);

    get_leaf_grids(&num_leaf_grids, leaf_grids, this);
    for (int i = 0; i < num_leaf_grids; i ++)
        if (words_are_the_same(leaf_grids[i]->coord_label, coord_label)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, leaf_grids[i]->super_grid_of_setting_coord_values == this,
                         "grid %s does not set coordinate %s, please use another grid to check\n", grid_name, coord_label);
            coord_values_this = (double*) leaf_grids[i]->get_grid_center_field()->get_grid_data_field()->data_buf;
        }

    if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
        for (long j = 0; j < grid_size; j ++)
            if (coord_values_this[j] != ((float*)coord_values_model)[j]) {
                printf("different coord value at %ld: %lf %f\n", j, coord_values_this[j], ((float*)coord_values_model)[j]);
                return false;
            }
    }
    else {
        for (long j = 0; j < grid_size; j ++)
            if (coord_values_this[j] != ((double*)coord_values_model)[j]) {
                printf("different coord value at %ld: %0.18lf %0.18lf %0.18lf\n", j, coord_values_this[j], ((double*)coord_values_model)[j], (double)((float)(((double*)coord_values_model)[j])));
                return false;
            }
    }

    return true;
}


bool Remap_grid_class::check_mask_values_consitency(const char *data_type, const void*mask_values_model)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, grid_mask_field != NULL, "the masks of grid %s has not been set, so as that we can not check the consistency of masks\n");

    if (words_are_the_same(data_type, DATA_TYPE_INT)) {
        for (long i = 0; i < grid_size; i ++)
            if (((bool*)grid_mask_field->get_grid_data_field()->data_buf)[i] != (((int*)mask_values_model)[i]>0)) {
                printf("different integer mask value at %ld: %d, %d\n", i, (bool)(((int*)mask_values_model)[i]), ((bool*)grid_mask_field->get_grid_data_field()->data_buf)[i]);
                return false;
            }
    }
    else if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
        for (long i = 0; i < grid_size; i ++)
            if (((bool*)grid_mask_field->get_grid_data_field()->data_buf)[i] != (((float*)mask_values_model)[i]>0)) {
                printf("different float mask value at %ld: %f, %d\n", ((float*)mask_values_model)[i], ((bool*)grid_mask_field->get_grid_data_field()->data_buf)[i]);
                return false;
            }
    }
    else {
        for (long i = 0; i < grid_size; i ++)
            if (((bool*)grid_mask_field->get_grid_data_field()->data_buf)[i] != (((double*)mask_values_model)[i]>0)) {
                printf("different double mask value at %ld: %lf, %d\n", ((double*)mask_values_model)[i], ((bool*)grid_mask_field->get_grid_data_field()->data_buf)[i]);
                return false;
            }
    }
    return true;
}


void Remap_grid_class::write_grid_field_into_array(Remap_grid_data_class *grid_field, char **array, long &buffer_max_size, long &buffer_content_size)
{
    int temp_int;

    
    if (grid_field != NULL) {
        temp_int = grid_size;
        grid_field->write_grid_data_into_array(array, buffer_max_size, buffer_content_size);
    }
    else temp_int = 0;
    write_data_into_array_buffer(&temp_int, sizeof(int), array, buffer_max_size, buffer_content_size);
}


void Remap_grid_class::read_grid_field_from_array(Remap_grid_data_class **grid_field, const char *array, long &buffer_content_iter)
{
    int temp_int;


    read_data_from_array_buffer(&temp_int, sizeof(int), array, buffer_content_iter, true);
    if (temp_int != 0)
        *grid_field = new Remap_grid_data_class(this, array, buffer_content_iter);
}


void Remap_grid_class::write_grid_name_into_array(Remap_grid_class *grid, char **array, long &buffer_max_size, long &buffer_content_size)
{
    char temp_grid_name[NAME_STR_SIZE];


    if (grid == NULL)
        sprintf(temp_grid_name, "NULL");
    else strcpy(temp_grid_name, grid->get_grid_name());
    write_data_into_array_buffer(temp_grid_name, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
}


Remap_grid_class *Remap_grid_class::get_linked_grid_from_array(Remap_grid_class *top_grid, const char *grid_name_suffix, char *temp_grid_name)
{
    if (!words_are_the_same(temp_grid_name, "NULL")) {
		strcat(temp_grid_name, "_FROM_");
        strcat(temp_grid_name, grid_name_suffix);
        if (words_are_the_same(this->grid_name, temp_grid_name))
            return this;
        Remap_grid_class *linked_grid = top_grid->search_sub_grid(temp_grid_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, linked_grid != NULL, "Software error in Remap_grid_class::Remap_grid_class: wrong linked grid");
        return linked_grid;
    }
    
    return NULL;
}



void Remap_grid_class::write_grid_into_array(char **array, long &buffer_max_size, long &buffer_content_size)
{
    int temp_int;

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, " write grid \"%s\"", grid_name);
    
    for (int i = sub_grids.size()-1; i >= 0 ; i --)
        sub_grids[i]->write_grid_into_array(array, buffer_max_size, buffer_content_size);
    temp_int = sub_grids.size();
    write_data_into_array_buffer(&temp_int, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_grid_name_into_array(first_super_grid_of_enable_setting_coord_value, array, buffer_max_size, buffer_content_size);
    write_grid_name_into_array(super_grid_of_setting_coord_values, array, buffer_max_size, buffer_content_size);
    write_grid_field_into_array(level_V3D_coord_dynamic_trigger_field, array, buffer_max_size, buffer_content_size);
    write_grid_field_into_array(hybrid_grid_coefficient_field, array, buffer_max_size, buffer_content_size);
    write_grid_field_into_array(sigma_grid_sigma_value_field, array, buffer_max_size, buffer_content_size);
    write_grid_field_into_array(level_V3D_coord_trigger_field, array, buffer_max_size, buffer_content_size);
    write_grid_field_into_array(grid_mask_field, array, buffer_max_size, buffer_content_size);
    write_grid_field_into_array(redundant_cell_mark_field, array, buffer_max_size, buffer_content_size);
    for (int i = grid_vertex_fields.size()-1; i >= 0; i --)
        grid_vertex_fields[i]->write_grid_data_into_array(array, buffer_max_size, buffer_content_size);
    temp_int = grid_vertex_fields.size();
    write_data_into_array_buffer(&temp_int, sizeof(int), array, buffer_max_size, buffer_content_size);
    for (int i = grid_center_fields.size()-1; i >=0 ; i --)
        grid_center_fields[i]->write_grid_data_into_array(array, buffer_max_size, buffer_content_size);
    temp_int = grid_center_fields.size();
    write_data_into_array_buffer(&temp_int, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&sigma_grid_scale_factor, sizeof(double), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&sigma_grid_top_value, sizeof(double), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&level_V3D_coord_trigger_field_specified, sizeof(bool), array, buffer_max_size, buffer_content_size);    
	write_data_into_array_buffer(&using_V3D_level_coord, sizeof(bool), array, buffer_max_size, buffer_content_size); 
    write_data_into_array_buffer(&are_vertex_values_set_in_default, sizeof(bool), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&cyclic, sizeof(bool), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&masks_are_known, sizeof(bool), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&num_vertexes, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&num_dimensions, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&boundary_min_lon, sizeof(double), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&boundary_max_lon, sizeof(double), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&boundary_min_lat, sizeof(double), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&boundary_max_lat, sizeof(double), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&grid_size, sizeof(long), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(decomp_name, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(coord_unit, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(coord_label, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(grid_name, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
}


void Remap_grid_class::reset_super_grids_of_setting_mask_value(Remap_grid_class *super_grid)
{
    for (int i = 0; i < super_grid->super_grids_of_setting_mask_value.size(); i ++)
        this->super_grids_of_setting_mask_value.push_back(super_grid->super_grids_of_setting_mask_value[i]);
    if (this->grid_mask_field != NULL)
        this->super_grids_of_setting_mask_value.push_back(this);
    for (int i = 0; i < sub_grids.size(); i ++)
        sub_grids[i]->reset_super_grids_of_setting_mask_value(this);
}


Remap_grid_class::Remap_grid_class(Remap_grid_class *top_grid, const char *grid_name_suffix, const char *array, long &buffer_content_iter)
{
    int temp_int;
    char temp_grid_name[NAME_STR_SIZE];

    
    if (top_grid == NULL)
        top_grid = this;
    
    initialize_grid_class_data();

    read_data_from_array_buffer(grid_name, NAME_STR_SIZE, array, buffer_content_iter, true);
	strcat(grid_name, "_FROM_");
    strcat(grid_name, grid_name_suffix);
    if (remap_grid_manager->search_remap_grid_with_grid_name(grid_name) == NULL)
        remap_grid_manager->add_remap_grid(this);
    else remap_grid_manager->add_temp_grid(this);
    read_data_from_array_buffer(coord_label, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(coord_unit, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(decomp_name, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(&grid_size, sizeof(long), array, buffer_content_iter, true);
    read_data_from_array_buffer(&boundary_max_lat, sizeof(double), array, buffer_content_iter, true);
    read_data_from_array_buffer(&boundary_min_lat, sizeof(double), array, buffer_content_iter, true);
    read_data_from_array_buffer(&boundary_max_lon, sizeof(double), array, buffer_content_iter, true);
    read_data_from_array_buffer(&boundary_min_lon, sizeof(double), array, buffer_content_iter, true);
    read_data_from_array_buffer(&num_dimensions, sizeof(int), array, buffer_content_iter, true);
    read_data_from_array_buffer(&num_vertexes, sizeof(int), array, buffer_content_iter, true);
    read_data_from_array_buffer(&masks_are_known, sizeof(bool), array, buffer_content_iter, true);
    read_data_from_array_buffer(&cyclic, sizeof(bool), array, buffer_content_iter, true);
    read_data_from_array_buffer(&are_vertex_values_set_in_default, sizeof(bool), array, buffer_content_iter, true);
	read_data_from_array_buffer(&using_V3D_level_coord, sizeof(bool), array, buffer_content_iter, true);
    read_data_from_array_buffer(&level_V3D_coord_trigger_field_specified, sizeof(bool), array, buffer_content_iter, true);
    read_data_from_array_buffer(&sigma_grid_top_value, sizeof(double), array, buffer_content_iter, true);
    read_data_from_array_buffer(&sigma_grid_scale_factor, sizeof(double), array, buffer_content_iter, true);
    read_data_from_array_buffer(&temp_int, sizeof(int), array, buffer_content_iter, true);
    for (int i = 0; i < temp_int; i ++)
        grid_center_fields.push_back(new Remap_grid_data_class(this, array, buffer_content_iter));
    read_data_from_array_buffer(&temp_int, sizeof(int), array, buffer_content_iter, true);
    for (int i = 0; i < temp_int; i ++)
        grid_vertex_fields.push_back(new Remap_grid_data_class(this, array, buffer_content_iter));
    read_grid_field_from_array(&redundant_cell_mark_field, array, buffer_content_iter);
    read_grid_field_from_array(&grid_mask_field, array, buffer_content_iter);
    read_grid_field_from_array(&level_V3D_coord_trigger_field, array, buffer_content_iter);
    read_grid_field_from_array(&sigma_grid_sigma_value_field, array, buffer_content_iter);
    read_grid_field_from_array(&hybrid_grid_coefficient_field, array, buffer_content_iter);
    read_grid_field_from_array(&level_V3D_coord_dynamic_trigger_field, array, buffer_content_iter);
    read_data_from_array_buffer(name_super_grid_of_setting_coord_values, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(name_first_super_grid_of_enable_setting_coord_value, NAME_STR_SIZE, array, buffer_content_iter, true);
    read_data_from_array_buffer(&temp_int, sizeof(int), array, buffer_content_iter, true);
    for (int i = 0; i < temp_int; i ++) {
        Remap_grid_class *child_grid = new Remap_grid_class(top_grid, grid_name_suffix, array, buffer_content_iter);
        Remap_grid_class *existing_grid = remap_grid_manager->search_remap_grid_with_grid_name(child_grid->get_grid_name());
        sub_grids.push_back(existing_grid);
    }
    
    if (this == top_grid) {
        link_grids(top_grid, grid_name_suffix);
        reset_super_grids_of_setting_mask_value(top_grid);
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "read boundary of grid \"%s\" (%lx) to %lf %lf %lf %lf", grid_name, this, boundary_min_lon, boundary_max_lon, boundary_min_lat, boundary_max_lat);
}


void Remap_grid_class::link_grids(Remap_grid_class *top_grid, const char *grid_name_suffix)
{
    if (super_grid_of_setting_coord_values == NULL)
        super_grid_of_setting_coord_values = get_linked_grid_from_array(top_grid, grid_name_suffix, name_super_grid_of_setting_coord_values);
    if (first_super_grid_of_enable_setting_coord_value == NULL)
        first_super_grid_of_enable_setting_coord_value = get_linked_grid_from_array(top_grid, grid_name_suffix, name_first_super_grid_of_enable_setting_coord_value);

    if (num_dimensions == 1)
        EXECUTION_REPORT(REPORT_ERROR, -1, first_super_grid_of_enable_setting_coord_value != NULL, "Software error in Remap_grid_class::link_grids");

    for (int i = 0; i < sub_grids.size(); i ++)
        sub_grids[i]->link_grids(top_grid, grid_name_suffix);    
    for (int i = 0; i < grid_center_fields.size(); i ++)
        grid_center_fields[i]->generate_grid_info(grid_center_fields[i]->coord_value_grid);
    for (int i = 0; i < grid_vertex_fields.size(); i ++)
        grid_vertex_fields[i]->generate_grid_info(grid_vertex_fields[i]->coord_value_grid);
    if (grid_mask_field != NULL)
        grid_mask_field->generate_grid_info(grid_mask_field->coord_value_grid);
    if (imported_area != NULL)
        imported_area->generate_grid_info(imported_area->coord_value_grid);
    if (original_grid_mask_field != NULL)
        original_grid_mask_field->generate_grid_info(original_grid_mask_field->coord_value_grid);
    if (redundant_cell_mark_field != NULL)
        redundant_cell_mark_field->generate_grid_info(redundant_cell_mark_field->coord_value_grid);
    if (hybrid_grid_coefficient_field != NULL)
        hybrid_grid_coefficient_field->generate_grid_info(hybrid_grid_coefficient_field->coord_value_grid);
    if (sigma_grid_sigma_value_field != NULL)
        sigma_grid_sigma_value_field->generate_grid_info(sigma_grid_sigma_value_field->coord_value_grid);
    if (level_V3D_coord_trigger_field != NULL)
        level_V3D_coord_trigger_field->generate_grid_info(level_V3D_coord_trigger_field->coord_value_grid);
    if (level_V3D_coord_dynamic_trigger_field != NULL)
        level_V3D_coord_dynamic_trigger_field->generate_grid_info(level_V3D_coord_dynamic_trigger_field->coord_value_grid);
}


bool Remap_grid_class::is_sub_grid_of_grid(Remap_grid_class *another_grid)
{
    if (this == another_grid)
        return true;

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "compare %s vs %s %d %d", this->get_grid_name(), another_grid->get_grid_name(), another_grid->sub_grids.size(), another_grid->get_num_dimensions());

    for (int i = 0; i < another_grid->sub_grids.size(); i ++)
        if (is_subset_of_grid(another_grid->sub_grids[i]))
            return true;
    return false;
}


bool Remap_grid_class::format_sub_grids(Remap_grid_class *top_grid)
{
    std::vector<Remap_grid_class *> formated_sub_grids;    
    std::vector<int> last_grid_ids;
    int i, j;
    

    formated_sub_grids.clear();
    last_grid_ids.clear();
    for (i = 0; i < sub_grids.size(); i ++) {
        for (j = 0; j < formated_sub_grids.size(); j ++)
            if (sub_grids[i]->is_subset_of_grid(formated_sub_grids[j]))
                break;
        if (j < formated_sub_grids.size() && last_grid_ids[j]+1 != i)
            return false;
        if (j == formated_sub_grids.size()) {
            if (sub_grids[i]->super_grid_of_setting_coord_values == NULL || sub_grids[i]->super_grid_of_setting_coord_values->is_sub_grid_of_grid(top_grid))
                formated_sub_grids.push_back(sub_grids[i]);
            else formated_sub_grids.push_back(sub_grids[i]->super_grid_of_setting_coord_values);
            last_grid_ids.push_back(i);
        }
    }

    sub_grids.clear();
    for (j = 0; j < formated_sub_grids.size(); j ++)
        sub_grids.push_back(formated_sub_grids[j]);

    for (int i = 0; i < sub_grids.size(); i ++)
        if (!sub_grids[i]->format_sub_grids(top_grid))
            return false;

    return true;
}


Remap_grid_class *Remap_grid_class::search_sub_grid(const char *grid_name)
{
    if (words_are_the_same(this->grid_name, grid_name))
        return this;

    for (int i = 0; i < sub_grids.size(); i ++) {
        Remap_grid_class *required_grid = sub_grids[i]->search_sub_grid(grid_name);
        if (required_grid != NULL)
            return required_grid;
    }

    return NULL;
}


Remap_grid_class *Remap_grid_class::get_sphere_sub_grid()
{
    if (get_is_sphere_grid())
        return this;

    for (int i = 0; i < sub_grids.size(); i ++)
        if (sub_grids[i]->get_sphere_sub_grid() != NULL)
            return sub_grids[i]->get_sphere_sub_grid();

    return NULL;
}


Remap_grid_data_class *Remap_grid_class::generate_mid_point_grid_field(Remap_grid_data_class *interface_level_grid_field)
{
    Remap_grid_data_class *mid_point_grid_field = interface_level_grid_field->duplicate_grid_data_field(mid_point_grid, 1, false, false);
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(mid_point_grid_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "Software error in Remap_grid_class::generate_mid_point_grid_field: wrong data type");
    for (int i = 0; i < mid_point_grid->grid_size; i ++)
        ((double*)mid_point_grid_field->get_grid_data_field()->data_buf)[i] = (((double*)interface_level_grid_field->get_grid_data_field()->data_buf)[i]+((double*)interface_level_grid_field->get_grid_data_field()->data_buf)[i+1]) / 2; 

    return mid_point_grid_field;
}


Remap_grid_class *Remap_grid_class::generate_mid_point_grid()
{
    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_dimensions == 1 && has_grid_coord_label(COORD_LABEL_LEV), "Software error in Remap_grid_class::generate_mid_point_grid: not a V1D grid");
    EXECUTION_REPORT(REPORT_ERROR, -1, this->interface_level_grid == NULL, "Software error in Remap_grid_class::generate_mid_point_grid: for mid-point grid of a mid-point grid");
    if (this->mid_point_grid != NULL)
        return this->mid_point_grid;

    mid_point_grid = duplicate_grid(this);
    mid_point_grid->interface_level_grid = this;
    mid_point_grid->grid_size = this->grid_size - 1;
    mid_point_grid->first_super_grid_of_enable_setting_coord_value = mid_point_grid;
    sprintf(mid_point_grid->grid_name, "mid_point_grid_for_%s", this->grid_name);
    if (this->grid_center_fields.size() == 1) {
        delete mid_point_grid->grid_center_fields[0];
        mid_point_grid->grid_center_fields[0] = this->generate_mid_point_grid_field(this->grid_center_fields[0]);
    }
    if (this->sigma_grid_sigma_value_field != NULL) {
        delete mid_point_grid->sigma_grid_sigma_value_field;
        mid_point_grid->sigma_grid_sigma_value_field = this->generate_mid_point_grid_field(this->sigma_grid_sigma_value_field);
    }
    if (this->hybrid_grid_coefficient_field != NULL) {
        delete mid_point_grid->hybrid_grid_coefficient_field;
        mid_point_grid->hybrid_grid_coefficient_field = this->generate_mid_point_grid_field(this->hybrid_grid_coefficient_field);
    }
    if (mid_point_grid->grid_vertex_fields.size() == 1) {
        delete mid_point_grid->grid_vertex_fields[0];
        mid_point_grid->grid_vertex_fields.clear();
    }

    return mid_point_grid;
}

