/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "original_grid_mgt.h"
#include <netcdf.h>
#include <unistd.h>
#include <math.h>


Original_grid_info::Original_grid_info(int comp_id, int grid_id, const char *grid_name, const char *annotation, Remap_grid_class *original_CoR_grid, bool model_registration)
{
    this->comp_id = comp_id;
    this->grid_id = grid_id;
    this->original_CoR_grid = original_CoR_grid;
    this->bottom_field_variation_type = BOTTOM_FIELD_VARIATION_UNSET;
	this->V3D_lev_field_variation_type = BOTTOM_FIELD_VARIATION_UNSET;
    this->bottom_field_name = NULL;
	this->V3D_lev_field_id = -1;
    this->bottom_field_id = -1;
    this->checksum_H2D_mask = -1;
    this->mid_point_grid = NULL;
    this->interface_level_grid = NULL;
    this->grid_name = strdup(grid_name);
    comp_full_name = strdup(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "Original_grid_info")->get_full_name());
    annotation_mgr->add_annotation(this->grid_id, "grid_registration", annotation);
    generate_remapping_grids();

    if (H2D_sub_CoR_grid != NULL && model_registration) {
        char *grid_data = (char *) (new double [H2D_sub_CoR_grid->get_grid_size()]);
        get_grid_data(-1, "mask", DATA_TYPE_INT, H2D_sub_CoR_grid->get_grid_size(), grid_data, "internal", "internal");
        checksum_H2D_mask = calculate_checksum_of_array(grid_data, H2D_sub_CoR_grid->get_grid_size(), sizeof(int), NULL, NULL);
        delete [] grid_data;
    }

    if (model_registration && H2D_sub_CoR_grid != NULL && V1D_sub_CoR_grid == NULL && T1D_sub_CoR_grid == NULL) {
        char nc_file_name[NAME_STR_SIZE];
        sprintf(nc_file_name, "%s/%s@%s.nc", comp_comm_group_mgt_mgr->get_internal_H2D_grids_dir(), grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "Original_grid_info")->get_full_name());
        if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in register_h2d_grid_with_data") == 0) {
            IO_netcdf *netcdf_file_object = new IO_netcdf("H2D_grid_data", nc_file_name, "w", false);
            netcdf_file_object->write_grid(H2D_sub_CoR_grid, false, true);
            netcdf_file_object->put_global_attr("edge_type", "LON_LAT", DATA_TYPE_STRING, DATA_TYPE_STRING, -1);   // to be modified
            Remap_grid_class *leaf_grids[256];
            int num_leaf_grids;
            H2D_sub_CoR_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, H2D_sub_CoR_grid);
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(leaf_grids[0]->get_coord_label(), COORD_LABEL_LON), "software error in Original_grid_info::Original_grid_info");
            double temp_double = H2D_sub_CoR_grid->get_boundary_min_lon();
            netcdf_file_object->put_global_attr("min_lon", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            temp_double = H2D_sub_CoR_grid->get_boundary_max_lon();
            netcdf_file_object->put_global_attr("max_lon", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            temp_double = H2D_sub_CoR_grid->get_boundary_min_lat();
            netcdf_file_object->put_global_attr("min_lat", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            temp_double = H2D_sub_CoR_grid->get_boundary_max_lat();
            netcdf_file_object->put_global_attr("max_lat", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            if (leaf_grids[0]->get_grid_cyclic())                    
                netcdf_file_object->put_global_attr("cyclic_or_acyclic", "cyclic", DATA_TYPE_STRING, DATA_TYPE_STRING, -1);
            else netcdf_file_object->put_global_attr("cyclic_or_acyclic", "acyclic", DATA_TYPE_STRING, DATA_TYPE_STRING, -1);            
            netcdf_file_object->put_global_attr("title", grid_name, DATA_TYPE_STRING, DATA_TYPE_STRING, -1);
            delete netcdf_file_object;
            char status_file_name[NAME_STR_SIZE];
            sprintf(status_file_name, "%s/%s@%s.nc.end", comp_comm_group_mgt_mgr->get_comps_ending_config_status_dir(), grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "Original_grid_info")->get_full_name());
            FILE *status_file = fopen(status_file_name, "w+");
            EXECUTION_REPORT(REPORT_ERROR, -1, status_file != NULL, "Software error in Comp_comm_group_mgt_node::merge_comp_comm_info: configuration ending status file cannot be created: %s", status_file_name);
            fclose(status_file);
        }
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Original_grid_info::Original_grid_info"));
    }
}


Original_grid_info::~Original_grid_info()
{
    delete [] grid_name;
    delete [] comp_full_name;
    if (bottom_field_name != NULL)
        delete [] bottom_field_name;
}


void Original_grid_info::reset_grid_data()
{
    if (original_CoR_grid != NULL) {
        delete original_CoR_grid;
        original_CoR_grid = NULL;
    }  
}



const char *Original_grid_info::get_annotation()
{
    return annotation_mgr->get_annotation(this->grid_id, "grid_registration");
}


void Original_grid_info::generate_remapping_grids()
{
    Remap_grid_class *lon_sub_grid = NULL, *lat_sub_grid = NULL, *leaf_grids[256];
    int num_leaf_grids;

    
    H2D_sub_CoR_grid = NULL;
    V1D_sub_CoR_grid = NULL;
    T1D_sub_CoR_grid = NULL;
    H2D_sub_grid_order = -1;
    V1D_sub_grid_order = -1;
    T1D_sub_grid_order = -1;
    
    if (original_CoR_grid->get_num_dimensions() == 2)
        EXECUTION_REPORT(REPORT_ERROR, -1, original_CoR_grid->get_is_sphere_grid(), "Software error in Original_grid_info::generate_remapping_grids: not a sphere grid");
    if (original_CoR_grid->get_num_dimensions() == 1)
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(original_CoR_grid->get_coord_label(),COORD_LABEL_LEV) || words_are_the_same(original_CoR_grid->get_coord_label(),COORD_LABEL_TIME), "Software error in Original_grid_info::generate_remapping_grids: not a vertical grid or not a time grid");
    H2D_sub_CoR_grid = original_CoR_grid->get_sphere_sub_grid();
    
    original_CoR_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, original_CoR_grid);
    for (int i = 0; i < num_leaf_grids; i ++) {
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LON)) {
            lon_sub_grid = leaf_grids[i];
            H2D_sub_grid_order = i;
        }
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LAT))
            lat_sub_grid = leaf_grids[i];
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LEV)) {
            V1D_sub_CoR_grid = leaf_grids[i];
            V1D_sub_grid_order = i;
        }
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_TIME)) {
            T1D_sub_CoR_grid = leaf_grids[i];
            T1D_sub_grid_order = i;
        }
    }

    if (H2D_sub_CoR_grid == NULL && (lon_sub_grid != NULL || lat_sub_grid != NULL)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, lon_sub_grid != NULL && lat_sub_grid != NULL, "Software error in Original_grid_info::generate_remapping_grids: wrong sub grids for a sphere grid");
        leaf_grids[0] = lon_sub_grid;
        leaf_grids[1] = lat_sub_grid;
        char H2D_grid_name[NAME_STR_SIZE];
        sprintf(H2D_grid_name, "H2D_grid_for_%s", grid_name);
        H2D_sub_CoR_grid = new Remap_grid_class(H2D_grid_name, 2, leaf_grids, 0);
    }    

    EXECUTION_REPORT(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL || V1D_sub_CoR_grid != NULL || T1D_sub_CoR_grid != NULL, "Software error in Original_grid_info::generate_remapping_grids: empty grid");
}


bool Original_grid_info::is_V1D_sub_grid_after_H2D_sub_grid()
{
    if (H2D_sub_CoR_grid == NULL || V1D_sub_CoR_grid == NULL)
        return true;

    return H2D_sub_grid_order < V1D_sub_grid_order;
}


void Original_grid_info::set_unique_3D_lev_field(int field_id, const char *type, const char *annotation) 
{ 
    EXECUTION_REPORT(REPORT_ERROR, -1, bottom_field_variation_type == BOTTOM_FIELD_VARIATION_UNSET && V3D_lev_field_variation_type == BOTTOM_FIELD_VARIATION_UNSET, "Software error in Original_grid_info::set_unique_3D_lev_field");
    annotation_mgr->add_annotation(grid_id, "set 3-D level field", annotation);
    V3D_lev_field_variation_type = words_are_the_same(type, "constant")? 0 : 1;
	V3D_lev_field_id = field_id;
} 


void Original_grid_info::set_unique_bottom_field(int field_id, int type, const char *annotation) 
{ 
    EXECUTION_REPORT(REPORT_ERROR, -1, bottom_field_variation_type == BOTTOM_FIELD_VARIATION_UNSET, "Software error in Original_grid_info::set_unique_bottom_field");
    annotation_mgr->add_annotation(grid_id, "set surface field", annotation);    
    bottom_field_variation_type = type;
    if (type != BOTTOM_FIELD_VARIATION_EXTERNAL) {
        Field_mem_info *field_inst = memory_manager->get_field_instance(field_id);
        bottom_field_name = strdup(field_inst->get_field_name());
        bottom_field_id = field_id;
        decomp_grids_mgr->search_decomp_grid_info(field_inst->get_decomp_id(), get_original_CoR_grid(), false)->get_decomp_grid()->set_level_V3D_coord_dynamic_trigger_field(field_inst->get_field_data());
    }
} 


void Original_grid_info::write_grid_into_array(char **temp_array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    get_original_CoR_grid()->write_grid_into_array(temp_array_buffer, buffer_max_size, buffer_content_size);
	write_data_into_array_buffer(&V3D_lev_field_variation_type, sizeof(int), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&bottom_field_variation_type, sizeof(int), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&checksum_H2D_mask, sizeof(long), temp_array_buffer, buffer_max_size, buffer_content_size);
}


void Original_grid_info::get_grid_data(int decomp_id, const char *label, const char *data_type, int array_size, char *grid_data, const char *annotation, const char *API_label)
{
    Remap_grid_data_class *grid_field;
    Remap_grid_class *field_CoR_grid;
    
        
    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(label, COORD_LABEL_LON) || words_are_the_same(label, COORD_LABEL_LAT) || words_are_the_same(label, GRID_MASK_LABEL), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the label (currently is \"%s\") is wrong (must be \"lon\", \"lat\" or \"mask\"). Please verify the model code with the annotation \"%s\".", API_label, grid_name, label, annotation);    
    if (words_are_the_same(label, GRID_MASK_LABEL))
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type, DATA_TYPE_INT), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the data type (currently is %s) corresponding to the parameter of \"grid_data\" is wrong when the label is \"mask\" (the right data type must be integer). Please verify the model code with the annotation \"%s\".", API_label, grid_name, data_type, annotation);
    else EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type, DATA_TYPE_FLOAT) || words_are_the_same(data_type, DATA_TYPE_DOUBLE), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the data type (currently is %s) corresponding to the parameter of \"grid_data\" is wrong when the label (currently is \"%s\") is not \"mask\" (the right data type must be floating-point). Please verify the model code with the annotation \"%s\".", API_label, grid_name, data_type, label, annotation);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, array_size > 0, "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the array of \"grid_data\" has not been allocated. Please verify the model code with the annotation \"%s\".", API_label, grid_name, annotation);
    if (decomp_id != -1) {
        field_CoR_grid = decomp_grids_mgr->search_decomp_grid_info(decomp_id, original_CoR_grid, false)->get_decomp_grid();
        EXECUTION_REPORT(REPORT_ERROR, comp_id, array_size == decomps_info_mgr->get_decomp_info(decomp_id)->get_num_local_cells(), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the array size (currently is %d) of the parameter of \"grid_data\" is different from the size corresponding to the parallel decomposition (currently is %d). Please verify the model code with the annotation \"%s\".", API_label, grid_name, array_size, decomps_info_mgr->get_decomp_info(decomp_id)->get_num_local_cells(), annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_original_CoR_grid(decomps_info_mgr->get_decomp_info(decomp_id)->get_grid_id()) == original_CoR_grid, "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the grid_id and decomp_id do not correspond to the same H2D grid. Please verify the model code with the annotation \"%s\".", API_label, grid_name, annotation);
    }
    else {
        field_CoR_grid = H2D_sub_CoR_grid;
        EXECUTION_REPORT(REPORT_ERROR, comp_id, array_size == field_CoR_grid->get_grid_size(), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the array size (currently is %d) of the parameter of \"grid_data\" is different from the size corresponding to the parallel decomposition (currently is %d). Please verify the model code with the annotation \"%s\".", API_label, grid_name, array_size, field_CoR_grid->get_grid_size(), annotation);
    }

    if (!words_are_the_same(label, GRID_MASK_LABEL)) {
        grid_field = field_CoR_grid->get_grid_center_field(label);
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_field != NULL, "Software error in Original_grid_info::get_grid_data: NULL grid center field");
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(grid_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "Software error in Original_grid_info::get_grid_data: wrong center field data type");
        if (array_size != grid_field->get_grid_data_field()->required_data_size)
            grid_field = field_CoR_grid->expand_to_generate_full_coord_value(grid_field);
        if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
            transform_datatype_of_arrays((double*)grid_field->get_grid_data_field()->data_buf, (float*)grid_data, array_size);
        else transform_datatype_of_arrays((double*)grid_field->get_grid_data_field()->data_buf, (double*)grid_data, array_size);
        if (grid_field != field_CoR_grid->get_grid_center_field(label))
            delete grid_field;
    }
    else {
        grid_field = field_CoR_grid->get_grid_mask_field();
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_field != NULL, "Software error in Original_grid_info::get_grid_data: NULL grid mask field");            
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(grid_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_BOOL), "Software error in Original_grid_info::get_grid_data: wrong mask field data type");
        EXECUTION_REPORT(REPORT_ERROR, -1, array_size == grid_field->get_grid_data_field()->required_data_size, "Software error in Original_grid_info::get_grid_data: inconsistent field size");
        transform_datatype_of_arrays((bool*)grid_field->get_grid_data_field()->data_buf, (int*)grid_data, array_size);
    }
}


void Original_grid_info::set_mid_point_grid(Original_grid_info *new_grid)
{
    this->mid_point_grid = new_grid;
    new_grid->interface_level_grid = this;
    if (this->bottom_field_variation_type != BOTTOM_FIELD_VARIATION_UNSET)
        new_grid->set_unique_bottom_field(this->bottom_field_id, this->bottom_field_variation_type, annotation_mgr->get_annotation(this->grid_id,"set surface field"));
}


void Original_grid_info::set_grid_checksum(long checksum_mask)
{
    this->checksum_H2D_mask = checksum_mask;
}


double *Original_grid_info::get_center_lon_values()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL, "Software error in Original_grid_info::get_center_lon_values");
	double *center_lon_values = new double [H2D_sub_CoR_grid->get_grid_size()];
	get_grid_data(-1, "lon", DATA_TYPE_DOUBLE, H2D_sub_CoR_grid->get_grid_size(), (char*) center_lon_values, "internal", "internal");
	return center_lon_values;
}


double *Original_grid_info::get_center_lat_values()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL, "Software error in Original_grid_info::get_center_lat_values");
	double *center_lat_values = new double [H2D_sub_CoR_grid->get_grid_size()];
	get_grid_data(-1, "lat", DATA_TYPE_DOUBLE, H2D_sub_CoR_grid->get_grid_size(), (char*) center_lat_values, "internal", "internal");
	return center_lat_values;
}


bool Original_grid_info::is_H2D_grid_and_the_same_as_another_grid(Original_grid_info *another_grid)
{
	bool check_result;
	double *this_center_lon, *this_center_lat, *another_center_lon, *another_center_lat;

	
    if (!this->is_H2D_grid() || !another_grid->is_H2D_grid())
        return false;
	
	this_center_lon = this->get_center_lon_values();
	another_center_lon = another_grid->get_center_lon_values();
	check_result = are_two_coord_arrays_same(this_center_lon, another_center_lon, this->get_H2D_sub_CoR_grid()->get_grid_size(), another_grid->get_H2D_sub_CoR_grid()->get_grid_size());
	delete [] this_center_lon;
	delete [] another_center_lon;	
	if (!check_result)
        return false;

	this_center_lat = this->get_center_lat_values();
	another_center_lat = another_grid->get_center_lat_values();
	check_result = are_two_coord_arrays_same(this_center_lat, another_center_lat, this->get_H2D_sub_CoR_grid()->get_grid_size(), another_grid->get_H2D_sub_CoR_grid()->get_grid_size());
	delete [] this_center_lat;
	delete [] another_center_lat;
	return check_result;
}


Original_grid_mgt::Original_grid_mgt()
{    
    original_grids.clear();
    CoR_script_name[0] = '\0';
    CoR_grids = NULL;
}


void Original_grid_mgt::initialize_CoR_grids()
{
    if (CoR_grids != NULL)
        return;
    
    sprintf(CoR_script_name, "%s/CCPL_grid.cor", comp_comm_group_mgt_mgr->get_root_comp_config_dir());
    FILE *fp = fopen(CoR_script_name, "r");
    if (fp == NULL)
        CoR_script_name[0] = '\0';
    else fclose(fp);
    if (strlen(CoR_script_name) != 0) {
        char current_dir[NAME_STR_SIZE], grids_dir[NAME_STR_SIZE];
        EXECUTION_REPORT(REPORT_ERROR, -1, getcwd(current_dir,NAME_STR_SIZE) != NULL, "Cannot get the current working directory for running the model");
        sprintf(grids_dir, "%s/grids_weights", comp_comm_group_mgt_mgr->get_root_comp_config_dir());
        EXECUTION_REPORT(REPORT_ERROR, -1, chdir(grids_dir) == 0, "Fail to access the directory of the CoR grid data files: \"%s\". Please verify.", grids_dir);
        CoR_grids = new Remap_mgt(CoR_script_name);
        chdir(current_dir);
    }
    else CoR_grids = new Remap_mgt(NULL);
}


Original_grid_mgt::~Original_grid_mgt()
{
    for (int i = 0; i < original_grids.size(); i ++)
        if (original_grids[i] != NULL)
            delete original_grids[i];
    
    delete CoR_grids;
}


void Original_grid_mgt::delete_external_original_grids()
{
    for (int i = 0; i < original_grids.size(); i ++) {
        if (original_grids[i] == NULL)
            continue;
        Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grids[i]->get_comp_full_name());
        if (comp_node == NULL || comp_node->get_current_proc_local_id() == -1) {
            original_grids[i]->reset_grid_data();
            delete original_grids[i];
            original_grids[i] = NULL;
        }
    }
}


Original_grid_info *Original_grid_mgt::search_grid_info(const char *grid_name, int comp_id)
{
    for (int i = 0; i < original_grids.size(); i ++)
        if (original_grids[i] != NULL)
            if (words_are_the_same(original_grids[i]->get_grid_name(), grid_name) && original_grids[i]->get_comp_id() == comp_id)
                return original_grids[i];

    return NULL;
}


Original_grid_info *Original_grid_mgt::search_grid_info(int grid_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_grid_id_legal(grid_id), "software error in Original_grid_mgt::search_grid_info based on grid_id");    
    return original_grids[grid_id&TYPE_ID_SUFFIX_MASK];
}


int Original_grid_mgt::register_H2D_grid_via_comp(int comp_id, const char *grid_name, const char *annotation)
{
    char XML_file_name[NAME_STR_SIZE], nc_file_name[NAME_STR_SIZE], status_file_name[NAME_STR_SIZE];
    const char *another_comp_full_name = NULL, *another_comp_grid_name = NULL;
    int line_number;
    

    sprintf(XML_file_name, "%s/all/coupling_connections/%s.coupling_connections.xml", comp_comm_group_mgt_mgr->get_config_root_dir(), comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in register_H2D_grid_via_comp")->get_full_name());
    TiXmlDocument *XML_file = open_XML_file_to_read(comp_id, XML_file_name, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_H2D_grid_via_comp"), false);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, XML_file != NULL, "Error happens when calling the API \"CCPL_register_H2D_grid_from_another_component\" to register an H2D grid \"%s\": the coupling connection configuration file (\"%s\") does not exist. The API call is at the model code with the annotation \"%s\". ", grid_name, XML_file_name, annotation);

    TiXmlElement *root_XML_element;
    TiXmlNode *root_XML_element_node = get_XML_first_child_of_unique_root(comp_id, XML_file_name, XML_file);
    for (; root_XML_element_node != NULL; root_XML_element_node = root_XML_element_node->NextSibling()) {
        if (root_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        root_XML_element = root_XML_element_node->ToElement();
        if (words_are_the_same(root_XML_element->Value(),"local_grids"))
            break;
    }
    if (root_XML_element_node != NULL) {
        for (TiXmlNode *grid_XML_element_node = root_XML_element->FirstChild(); grid_XML_element_node != NULL; grid_XML_element_node = grid_XML_element_node->NextSibling()) {
            if (grid_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                continue;
            TiXmlElement *grid_XML_element = grid_XML_element_node->ToElement();
            const char *xml_grid_name = get_XML_attribute(comp_id, CCPL_NAME_STR_LEN, grid_XML_element, "local_grid_name", XML_file_name, line_number, "grid name of the current component", "the coupling connection configuration file", true);
            if (words_are_the_same(xml_grid_name, grid_name)) {
                another_comp_full_name = get_XML_attribute(comp_id, 512, grid_XML_element, "another_comp_full_name", XML_file_name, line_number, "the full name of the another component", "the coupling connection configuration file", true);
                another_comp_grid_name = get_XML_attribute(comp_id, CCPL_NAME_STR_LEN, grid_XML_element, "another_comp_grid_name", XML_file_name, line_number, "the grid name of the another component", "the coupling connection configuration file", true);
                EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(another_comp_grid_name) > 0, "Error happens when calling the API \"CCPL_register_H2D_grid_from_another_component\" to register an H2D grid \"%s\": the coupling connection configuration file (\"%s\") specifies an empty name of the remote grid. Please check the XML file around line number %d", grid_name, XML_file_name, line_number);
                break;
            }
        }    
    }

    EXECUTION_REPORT(REPORT_ERROR, comp_id, another_comp_full_name != NULL && another_comp_grid_name != NULL, "Error happens when calling the API \"CCPL_register_H2D_grid_from_another_component\" to register an H2D grid \"%s\": the coupling connection configuration file (\"%s\") does not contain the information for this grid. The API call is at the model code with the annotation \"%s\". ", grid_name, XML_file_name, annotation);
    sprintf(nc_file_name, "%s/%s@%s.nc", comp_comm_group_mgt_mgr->get_internal_H2D_grids_dir(), another_comp_grid_name, another_comp_full_name);
    sprintf(status_file_name, "%s/%s@%s.nc.end", comp_comm_group_mgt_mgr->get_comps_ending_config_status_dir(), another_comp_grid_name, another_comp_full_name);
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Wait to read NetCDF file \"%s\" to register H2D grid \"%s\" based on the grid \"%s\" of remote component \"%s\". Dead wait will be encounted if the full name of the remote component is wrong. So please make sure the full name of the remote component is correct in the the coupling connection configuration file (\"%s\")", nc_file_name, grid_name, another_comp_grid_name, another_comp_full_name, XML_file_name);
    if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in register_H2D_grid_via_comp") == 0) {
        while (true) {
            FILE *status_file = fopen(status_file_name, "r");
            if (status_file == NULL) {
                if (comp_comm_group_mgt_mgr->has_comp_ended_configuration(another_comp_full_name))
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Fail to read NetCDF file \"%s\" to register H2D grid \"%s\": the remote component model \"%s\" has ended its coupling configuration stage without registering the required grid (\"%s\") before. Please check the configuration file (\"%s\") or the corresponding model code.", nc_file_name, grid_name, another_comp_full_name, another_comp_grid_name, XML_file_name);
                continue;
            }
            fclose(status_file);
            break;
        }
    }
    
    delete XML_file;
    synchronize_comp_processes_for_API(comp_id, API_ID_GRID_MGT_REG_H2D_GRID_VIA_COMP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_H2D_grid_via_comp"), "register_H2D_grid_via_comp", annotation);
    if (!report_error_enabled)
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_H2D_grid_via_comp"));

    return register_H2D_grid_via_file(comp_id, grid_name, nc_file_name, annotation);
}


void Original_grid_mgt::calculate_min_max_H2D_coord_value(int comp_id, char *center_coord, char *vertex_coord, int size_center_coord, int size_vertex_coord, const char *data_type, double &min_coord_value, double &max_coord_value)
{
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
        float float_min_value = (float)NULL_COORD_VALUE, float_max_value = (float)NULL_COORD_VALUE;
        get_min_value_in_array((float*)center_coord, size_center_coord, true, (float)NULL_COORD_VALUE, float_min_value);
        get_min_value_in_array((float*)vertex_coord, size_vertex_coord, true, (float)NULL_COORD_VALUE, float_min_value);
        get_max_value_in_array((float*)center_coord, size_center_coord, true, (float)NULL_COORD_VALUE, float_max_value);
        get_max_value_in_array((float*)vertex_coord, size_vertex_coord, true, (float)NULL_COORD_VALUE, float_max_value);
        min_coord_value = float_min_value;
        max_coord_value = float_max_value;
    }
    else {
        min_coord_value = NULL_COORD_VALUE, max_coord_value = NULL_COORD_VALUE;
        get_min_value_in_array((double*)center_coord, size_center_coord, true, NULL_COORD_VALUE, min_coord_value);
        get_min_value_in_array((double*)vertex_coord, size_vertex_coord, true, NULL_COORD_VALUE, min_coord_value);
        get_max_value_in_array((double*)center_coord, size_center_coord, true, NULL_COORD_VALUE, max_coord_value);
        get_max_value_in_array((double*)vertex_coord, size_vertex_coord, true, NULL_COORD_VALUE, max_coord_value);
    }
    double *all_min_coord_values = new double [comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_num_procs()];
    double *all_max_coord_values = new double [comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_num_procs()];
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The min and max values calculated by C-Coupler is %lf and %lf", min_coord_value, max_coord_value);
    MPI_Gather(&min_coord_value, 1, MPI_DOUBLE, all_min_coord_values, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,""));
    MPI_Gather(&max_coord_value, 1, MPI_DOUBLE, all_max_coord_values, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,""));
    if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id,"") == 0)
        for (int i = 0; i < comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_num_procs(); i ++) {
            if (are_floating_values_equal(NULL_COORD_VALUE, min_coord_value) || !are_floating_values_equal(NULL_COORD_VALUE, all_min_coord_values[i]) && min_coord_value > all_min_coord_values[i])
                min_coord_value = all_min_coord_values[i];
            if (are_floating_values_equal(NULL_COORD_VALUE, max_coord_value) || !are_floating_values_equal(NULL_COORD_VALUE, all_max_coord_values[i]) && max_coord_value < all_max_coord_values[i])
                max_coord_value = all_max_coord_values[i];
        }
    MPI_Bcast(&min_coord_value, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,""));
    MPI_Bcast(&max_coord_value, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,""));
    delete [] all_min_coord_values;
    delete [] all_max_coord_values;
}


void Original_grid_mgt::common_checking_for_H2D_registration_via_data(int comp_id, const char *grid_name, const char *edge_type, const char *coord_unit, char *cyclic_or_acyclic, const char *data_type, int size_mask, int size_center_lon, 
                                                                      int size_center_lat, int size_vertex_lon, int size_vertex_lat, int *mask, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, char *vertex_lon, char *vertex_lat, const char *annotation, int API_id)
{
    char API_label[NAME_STR_SIZE], hint[NAME_STR_SIZE];
    double eps = 1.0000001;
    double min_lon_value, max_lon_value, min_lat_value, max_lat_value;

    
    get_API_hint(comp_id, API_id, API_label);
    sprintf(hint, "registering an H2D grid \"%s\"", grid_name);
    
    transform_datatype_of_arrays(min_lon, (char*)(&min_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(max_lon, (char*)(&max_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(min_lat, (char*)(&min_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(max_lat, (char*)(&max_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, edge_type, "edge_type", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, cyclic_or_acyclic, "cyclic_or_acyclic", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, data_type, "implicit data type", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, min_lon_value, "\"min_lon\"", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, max_lon_value, "\"max_lon\"", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, min_lat_value, "\"min_lat\"", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, max_lat_value, "\"max_lat\"", annotation);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(edge_type, "LON_LAT") || words_are_the_same(edge_type, "GREAT_ARC") || words_are_the_same(edge_type, "XY") || words_are_the_same(edge_type, "TriPolar"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is \"%s\") of parameter \"edge_type\" is not \"LON_LAT\", \"GREAT_ARC\", \"TriPolar\" or \"XY\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, edge_type, annotation);
    if (words_are_the_same(edge_type, "LON_LAT") || words_are_the_same(edge_type, "GREAT_ARC") || words_are_the_same(edge_type, "TriPolar"))
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(coord_unit, COORD_UNIT_DEGREES) || words_are_the_same(coord_unit, COORD_UNIT_RADIANS), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the parameter \"coord_unit\" (currently is \"%s\") is not \"degrees\" or \"radians\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, coord_unit, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(cyclic_or_acyclic, "cyclic") || words_are_the_same(cyclic_or_acyclic, "acyclic"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is \"%s\") of parameter \"cyclic_or_acyclic\" is not \"cyclic\" or \"acyclic\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, cyclic_or_acyclic, annotation);
    if (words_are_the_same(edge_type, "XY"))
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(cyclic_or_acyclic, "acyclic"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is \"%s\") of parameter \"cyclic_or_acyclic\" is not \"acyclic\" when the parameter \"edge_type\" has been set to \"XY\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, cyclic_or_acyclic, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type, DATA_TYPE_FLOAT) || words_are_the_same(data_type, DATA_TYPE_DOUBLE), "software error in register_h2d_grid_with_data: wrong implicit data type");

    EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries("integer", mask, size_mask, 0, 1, 0, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"mask\" are wrong (not 0 and 1). Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, are_floating_values_equal(NULL_COORD_VALUE, min_lon_value) && are_floating_values_equal(NULL_COORD_VALUE, max_lon_value) || !are_floating_values_equal(NULL_COORD_VALUE, min_lon_value) && !are_floating_values_equal(NULL_COORD_VALUE, max_lon_value), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": \"min_lon\" (currently is %lf) and \"max_lon\" (currently is %lf) must be CCPL_NULL_COORD_VALUE (%lf) or not at the same time. Please check the model code related to the annotation \"%s\"", grid_name, API_label, min_lon_value, max_lon_value, NULL_COORD_VALUE, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, are_floating_values_equal(NULL_COORD_VALUE, min_lat_value) && are_floating_values_equal(NULL_COORD_VALUE, max_lat_value) || !are_floating_values_equal(NULL_COORD_VALUE, min_lat_value) && !are_floating_values_equal(NULL_COORD_VALUE, max_lat_value), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": \"min_lat\" (currently is %lf) and \"max_lat\" (currently is %lf) must be CCPL_NULL_COORD_VALUE (%lf) or not at the same time. Please check the model code related to the annotation \"%s\"", grid_name, API_label, min_lat_value, max_lat_value, NULL_COORD_VALUE, annotation);

    if (!are_floating_values_equal(NULL_COORD_VALUE, min_lat_value))
        EXECUTION_REPORT(REPORT_ERROR, comp_id, max_lat_value > min_lat_value, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": \"min_lat\" (currently is %lf) is not smaller than \"max_lat\" (currently is %lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, max_lat_value, annotation);
    else {
        calculate_min_max_H2D_coord_value(comp_id, center_lat, vertex_lat, size_center_lat, size_vertex_lat, data_type, min_lat_value, max_lat_value);
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The min and max latitude values of the grid \"%s\" calculated by C-Coupler is %lf and %lf", grid_name, min_lat_value, max_lat_value);
        transform_datatype_of_arrays((char*)(&min_lat_value), min_lat, DATA_TYPE_DOUBLE, data_type, 1);
        transform_datatype_of_arrays((char*)(&max_lat_value), max_lat, DATA_TYPE_DOUBLE, data_type, 1);
    }    
    if (words_are_the_same(coord_unit, COORD_UNIT_DEGREES)) {
        if (!are_floating_values_equal(NULL_COORD_VALUE, min_lat_value)) {
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lat_value, 1, (double) -90.0*eps, (double) 90.0*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lat\" is wrong (not between -90 and 90). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, annotation);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lat_value, 1, (double) -90.0*eps, (double) 90.0*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lat\" is wrong (not between -90 and 90). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lat_value, annotation);
        }
        if (are_floating_values_equal((double)-90.0, min_lat_value) || are_floating_values_equal((double)90.0, max_lat_value))
            strcpy(cyclic_or_acyclic, "cyclic");    
        if (words_are_the_same(cyclic_or_acyclic, "cyclic")) {
            min_lon_value = -360;
            max_lon_value = 360;
        }
        if (!are_floating_values_equal(NULL_COORD_VALUE, min_lon_value)) {
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lon_value, 1, (double) -360.0*eps, (double) 360.0*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lon\" is wrong (not between -360 and 360). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, annotation);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lon_value, 1, (double) -360.0*eps, (double) 360.0*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lon\" is wrong (not between -360 and 360). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lon_value, annotation);
            if (words_are_the_same(cyclic_or_acyclic, "acyclic"))
                EXECUTION_REPORT(REPORT_ERROR, comp_id, fabs(max_lon_value-min_lon_value) <= ((double)360.0)*eps, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the difference between \"min_lon\" and \"max_lon\" (%lf) is wrong (not between -360 and 360). Please check the model code related to the annotation \"%s\".", grid_name, API_label, fabs(max_lon_value-min_lon_value), annotation);            
        }
    }
    else if (words_are_the_same(coord_unit, COORD_UNIT_RADIANS)) {
        if (!are_floating_values_equal(NULL_COORD_VALUE, min_lat_value)) {
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lat_value, 1, -((double)3.1415927)/2*eps, ((double)3.1415927)/2*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lat\" are wrong (not between -PI/2 and PI/2). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, annotation);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lat_value, 1, -((double)3.1415927)/2*eps, ((double)3.1415927)/2*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lat\" are wrong (not between -PI/2 and PI/2). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lat_value, annotation);
        }
        if (are_floating_values_equal((double)-PI/2, min_lat_value) || are_floating_values_equal((double)PI/2, max_lat_value))
            strcpy(cyclic_or_acyclic, "cyclic");
        if (words_are_the_same(cyclic_or_acyclic, "cyclic")) {
            min_lon_value = -((double)3.1415927)*2;
            max_lon_value = ((double)3.1415927)*2;
        }
        if (!are_floating_values_equal(NULL_COORD_VALUE, min_lon_value)) {
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lon_value, 1, -((double)3.1415927)*2*eps, ((double)3.1415927)*2*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lon\" are wrong (not between -2PI and 2PI). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, annotation);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lon_value, 1, -((double)3.1415927)*2*eps, ((double)3.1415927)*2*eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lon\" are wrong (not between -2PI and 2PI). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lon_value, annotation);
            if (words_are_the_same(cyclic_or_acyclic, "acyclic"))
                EXECUTION_REPORT(REPORT_ERROR, comp_id, fabs(max_lon_value-min_lon_value) <= ((double)3.1415927)*2*eps, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the difference between \"min_lon\" and \"max_lon\" is wrong (not between -2PI and 2PI). Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);    
        }
    }
    if (are_floating_values_equal(NULL_COORD_VALUE, min_lon_value)) {
        calculate_min_max_H2D_coord_value(comp_id, center_lon, vertex_lon, size_center_lon, size_vertex_lon, data_type, min_lon_value, max_lon_value);
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The min and max longitude values of the grid \"%s\" calculated by C-Coupler is %lf and %lf", grid_name, min_lon_value, max_lon_value);
        transform_datatype_of_arrays((char*)(&min_lon_value), min_lon, DATA_TYPE_DOUBLE, data_type, 1);
        transform_datatype_of_arrays((char*)(&max_lon_value), max_lon, DATA_TYPE_DOUBLE, data_type, 1);
    }
    if (!are_floating_values_equal(NULL_COORD_VALUE, min_lon_value)) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) center_lon, size_center_lon, min_lon_value, max_lon_value, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"center_lon\" are not between \"min_lon\" (%lf) and \"max_lon\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, max_lon_value, annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) vertex_lon, size_vertex_lon, min_lon_value, max_lon_value, (double) NULL_COORD_VALUE, true), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"vertex_lon\" are not between \"min_lon\" (%lf) and \"max_lon\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, max_lon_value, annotation);
    }
    if (!are_floating_values_equal(NULL_COORD_VALUE, min_lat_value)) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) center_lat, size_center_lat, min_lat_value, max_lat_value, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"center_lat\" are not between \"min_lat\" (%lf) and \"max_lat\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, max_lat_value, annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) vertex_lat, size_vertex_lat, min_lat_value, max_lat_value, (double) NULL_COORD_VALUE, true), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"vertex_lat\" are not between \"min_lat\" (%lf) and \"max_lat\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, max_lat_value, annotation);
    }
}


int Original_grid_mgt::register_H2D_grid_via_local_data(int comp_id, const char *grid_name, const char *edge_type, const char *coord_unit, char *cyclic_or_acyclic, const char *data_type, int grid_size, int num_local_cells, int size_local_cells_global_index, int size_center_lon, int size_center_lat, 
                                                   int size_mask, int size_area, int size_vertex_lon, int size_vertex_lat, int *local_cells_global_index, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *decomp_name, int *decomp_id, const char *annotation, int API_id)
{
    int data_type_size;
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data");
    char API_label[NAME_STR_SIZE], hint[NAME_STR_SIZE];
    int global_center_lon_size, global_center_lat_size, global_area_size, global_mask_size, global_vertex_lon_size, global_vertex_lat_size;
    


    get_API_hint(comp_id, API_id, API_label);
    sprintf(hint, "registering an H2D grid \"%s\"", grid_name);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_local_cells >= 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the specified number (%d) of local grid cells is wrong (smaller than 0). Please check the model code related to \"%s\"", API_label, grid_name, num_local_cells, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_local_cells <= size_local_cells_global_index, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array size (currently is %d) of \"local_cells_global_index\" is smaller than \"num_local_cells\" (currently is %d). Please check the model code related to \"%s\"", API_label, grid_name, num_local_cells, size_local_cells_global_index, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, *decomp_id == 0 && strlen(decomp_name) > 0 || *decomp_id == -1 && strlen(decomp_name) == 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": parameters \"decomp_name\" and \"decomp_id\" are not speicified/unspecified at the same time. Please check the model code related to the annotation \"%s\"", API_label, grid_name, annotation);
    
    common_checking_for_H2D_registration_via_data(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, data_type, size_mask, size_center_lon, size_center_lat, 
                                                  size_vertex_lon, size_vertex_lat, mask, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, vertex_lon, vertex_lat, annotation, API_id);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, grid_size > 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the \"grid_size\" (currently is %d) is too small (smaller than 1). Please check the model code related to \"%s\"", API_label, grid_name, grid_size, annotation);
    check_API_parameter_int(comp_id, API_id, comm, NULL, grid_size, "\"grid_size\"", annotation);
    for (int i = 0; i < num_local_cells; i ++)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, local_cells_global_index[i] >= 1 && local_cells_global_index[i] <= grid_size, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": some values in \"local_cells_global_index\" (for example, %d) are out of bound (smaller than 1 or larger than \"grid_size\" (%d)). Please check the model code related to the annotation \"%s\".", API_label, grid_name, local_cells_global_index[i], grid_size, annotation);

    data_type_size = get_data_type_size(data_type);
    char *global_center_lon = check_and_aggregate_local_grid_data(comp_id, API_id, comm, hint, grid_size, size_center_lon, data_type_size, (char*) center_lon, "center_lon", num_local_cells, local_cells_global_index, global_center_lon_size, annotation);
    char *global_center_lat = check_and_aggregate_local_grid_data(comp_id, API_id, comm, hint, grid_size, size_center_lat, data_type_size, (char*) center_lat, "center_lat", num_local_cells, local_cells_global_index, global_center_lat_size, annotation);
    char *global_area = check_and_aggregate_local_grid_data(comp_id, API_id, comm, hint, grid_size, size_area, data_type_size, (char*) area, "area", num_local_cells, local_cells_global_index, global_area_size, annotation);
    char *global_vertex_lon = check_and_aggregate_local_grid_data(comp_id, API_id, comm, hint, grid_size, size_vertex_lon, data_type_size, (char*) vertex_lon, "vertex_lon", num_local_cells, local_cells_global_index, global_vertex_lon_size, annotation);
    char *global_vertex_lat = check_and_aggregate_local_grid_data(comp_id, API_id, comm, hint, grid_size, size_vertex_lat, data_type_size, (char*) vertex_lat, "vertex_lat", num_local_cells, local_cells_global_index, global_vertex_lat_size, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, global_vertex_lon_size == global_vertex_lat_size, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": \"vertex_lon\" or \"vertex_lat\" are not consistent with each other (not specified/unspecified at the same time or have different array sizes (%d and %d respectively)). Please check the model code related to the annotation \"%s\".", API_label, grid_name, global_vertex_lon_size, global_vertex_lat_size, annotation);

    int *global_mask = (int*) check_and_aggregate_local_grid_data(comp_id, API_id, comm, hint, grid_size, size_mask, sizeof(int), (char*) mask, "mask", num_local_cells, local_cells_global_index, global_mask_size, annotation);
    if (global_mask == NULL) {
        int *active_mask = NULL;
        if (num_local_cells > 0) {
            active_mask = new int [num_local_cells];
            for (int i = 0; i < num_local_cells; i ++)
                active_mask[i] = 1;
        }
        global_mask = (int*) check_and_aggregate_local_grid_data(comp_id, API_id, comm, hint, grid_size, num_local_cells, sizeof(int), (char*) active_mask, "mask", num_local_cells, local_cells_global_index, global_mask_size, annotation);
        if (active_mask != NULL)
            delete [] active_mask;
    }
    
    int grid_id = create_H2D_grid_from_global_data(comp_id, grid_name, coord_unit, cyclic_or_acyclic, data_type, grid_size, 0, global_vertex_lon_size/grid_size, global_center_lon_size, global_center_lat_size, 
                                                   global_mask_size, global_area_size, global_vertex_lon_size, global_vertex_lat_size, min_lon, max_lon, min_lat, max_lat, global_center_lon, global_center_lat, global_mask, global_area, global_vertex_lon, global_vertex_lat, annotation);
    
    if (strlen(decomp_name) != 0)
        *decomp_id = decomps_info_mgr->register_H2D_parallel_decomposition(decomp_name, grid_id, num_local_cells, local_cells_global_index, annotation);

    return grid_id;
}


int Original_grid_mgt::create_H2D_grid_from_global_data(int comp_id, const char *grid_name, const char *coord_unit, const char *cyclic_or_acyclic, const char *data_type, int dim_size1, int dim_size2, int num_vertex, int size_center_lon,
                                                   int size_center_lat, int size_mask, int size_area, int size_vertex_lon, int size_vertex_lat, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *annotation)
{
    char true_H2D_grid_name[NAME_STR_SIZE], true_lon_grid_name[NAME_STR_SIZE], true_lat_grid_name[NAME_STR_SIZE];
    Remap_grid_class *CoR_H2D_grid, *CoR_lon_grid, *CoR_lat_grid, *sub_grids[256];
    int grid_size = dim_size2 == 0? dim_size1 : dim_size1*dim_size2;
    double min_lon_value, max_lon_value, min_lat_value, max_lat_value;


    sprintf(true_H2D_grid_name, "%s@%s", grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, annotation)->get_full_name());
    sprintf(true_lon_grid_name, "lon_%s", true_H2D_grid_name);
    sprintf(true_lat_grid_name, "lat_%s", true_H2D_grid_name);
    if (dim_size2 == 0) {
        CoR_lon_grid = new Remap_grid_class(true_lon_grid_name, "lon", coord_unit, cyclic_or_acyclic, 0);
        CoR_lat_grid = new Remap_grid_class(true_lat_grid_name, "lat", coord_unit, cyclic_or_acyclic, 0);
    }
    else {
        CoR_lon_grid = new Remap_grid_class(true_lon_grid_name, "lon", coord_unit, cyclic_or_acyclic, dim_size1);
        CoR_lat_grid = new Remap_grid_class(true_lat_grid_name, "lat", coord_unit, cyclic_or_acyclic, dim_size2);        
    }
    sub_grids[0] = CoR_lon_grid;
    sub_grids[1] = CoR_lat_grid;
    CoR_H2D_grid = new Remap_grid_class(true_H2D_grid_name, 2, sub_grids, grid_size);
    if (size_mask > 0)
        CoR_H2D_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)mask, 0);
    else {
        int *tmp_mask = new int [CoR_H2D_grid->get_grid_size()];
        for (int i = 0; i < CoR_H2D_grid->get_grid_size(); i ++)
            tmp_mask[i] = 1;
        CoR_H2D_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)tmp_mask, 0);
        delete [] tmp_mask;
    }
    if (size_area > 0)
        CoR_H2D_grid->read_grid_data_from_array("area", "area", data_type, (const char*)area, 0);    
    if (size_center_lon == grid_size) {
        CoR_H2D_grid->read_grid_data_from_array("center", COORD_LABEL_LON, data_type, (const char*)center_lon, 0);
        CoR_H2D_grid->read_grid_data_from_array("center", COORD_LABEL_LAT, data_type, (const char*)center_lat, 0);
        if (size_vertex_lon >0) {
            CoR_H2D_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, data_type, (const char*)vertex_lon, num_vertex);
            CoR_H2D_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, data_type, (const char*)vertex_lat, num_vertex);        
        }
    }
    else {
        CoR_lon_grid->read_grid_data_from_array("center", COORD_LABEL_LON, data_type, (const char*)center_lon, 0);
        CoR_lat_grid->read_grid_data_from_array("center", COORD_LABEL_LAT, data_type, (const char*)center_lat, 0);        
        if (size_vertex_lon > 0) {
            CoR_lon_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, data_type, (const char*)vertex_lon, num_vertex);
            CoR_lat_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, data_type, (const char*)vertex_lat, num_vertex);        
        }
    }

    transform_datatype_of_arrays(min_lon, (char*)(&min_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(max_lon, (char*)(&max_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(min_lat, (char*)(&min_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(max_lat, (char*)(&max_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
    if (words_are_the_same(coord_unit,COORD_UNIT_RADIANS)) {
        min_lon_value = min_lon_value*180/PI;
        max_lon_value = max_lon_value*180/PI;
        min_lat_value = min_lat_value*180/PI;
        max_lat_value = max_lat_value*180/PI;
    }
    if (words_are_the_same(cyclic_or_acyclic, "cyclic")) {
        min_lon_value = 0.0;
        max_lon_value = 360.0;
    }
    
    remap_grid_manager->add_remap_grid(CoR_lon_grid);
    remap_grid_manager->add_remap_grid(CoR_lat_grid);
    remap_grid_manager->add_remap_grid(CoR_H2D_grid);
    CoR_H2D_grid->set_grid_boundary(min_lon_value, max_lon_value, min_lat_value, max_lat_value);
    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, CoR_H2D_grid, true));
    CoR_H2D_grid->end_grid_definition_stage(NULL);
	
    return original_grids[original_grids.size()-1]->get_grid_id();
}


int Original_grid_mgt::register_H2D_grid_via_global_data(int comp_id, const char *grid_name, const char *edge_type, const char *coord_unit, char *cyclic_or_acyclic, const char *data_type, int dim_size1, int dim_size2, int size_center_lon, int size_center_lat, 
                                                   int size_mask, int size_area, int size_vertex_lon, int size_vertex_lat, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *annotation, int API_id)
{
    int data_type_size, grid_size, num_vertex;
    char true_H2D_grid_name[NAME_STR_SIZE], true_lon_grid_name[NAME_STR_SIZE], true_lat_grid_name[NAME_STR_SIZE];
    char coord_label[NAME_STR_SIZE], coord_name[NAME_STR_SIZE], API_label[NAME_STR_SIZE], hint[NAME_STR_SIZE];
    Remap_grid_class *CoR_H2D_grid, *CoR_lon_grid, *CoR_lat_grid, *sub_grids[256];
    

    get_API_hint(comp_id, API_id, API_label);
    sprintf(hint, "registering an H2D grid \"%s\"", grid_name);
    
    common_checking_for_H2D_registration_via_data(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, data_type, size_mask, size_center_lon, size_center_lat, 
                                                  size_vertex_lon, size_vertex_lat, mask, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, vertex_lon, vertex_lat, annotation, API_id);
    
    data_type_size = get_data_type_size(data_type);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), hint, size_center_lon, data_type_size, (const char*)center_lon, "center_lon", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), hint, size_center_lat, data_type_size, (const char*)center_lat, "center_lat", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), hint, size_mask, sizeof(int), (const char*)mask, "mask", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), hint, size_area, data_type_size, (const char*)area, "area", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), hint, size_vertex_lon, data_type_size, (const char*)vertex_lon, "vertex_lon", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), hint, size_vertex_lat, data_type_size, (const char*)vertex_lat, "vertex_lat", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), NULL, dim_size1, "dim_size1", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), NULL, dim_size2, "dim_size2", annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, dim_size1 > 3, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is %d) of the parameter \"dim_size1\" is wrong. It must be larger than 3. Please check the model code related to the annotation \"%s\".", grid_name, API_label, dim_size1, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, dim_size2 == 0 || dim_size2 > 3, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is %d) of the parameter \"dim_size2\" is wrong. It must be 0 or a positive value larger than 3. Please check the model code related to the annotation \"%s\".", grid_name, API_label, dim_size2, annotation);
    if (dim_size2 == 0) {
        grid_size = dim_size1;
        EXECUTION_REPORT(REPORT_ERROR, comp_id, size_center_lon == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lon\" is different from the grid size that is determined by \"dim_size1\" (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lon, grid_size, annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, size_center_lat == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lat\" is different from the grid size that is determined by \"dim_size1\" (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lat, grid_size, annotation);
    }
    else {
        grid_size = (dim_size1)*(dim_size2);
        if (size_center_lon == dim_size1)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, size_center_lat == dim_size2, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lat\" is different from \"dim_size2\" (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lat, dim_size2, annotation);
        else {            
            EXECUTION_REPORT(REPORT_ERROR, comp_id, size_center_lon == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lon\" is different from \"dim_size1\" (currently is %d) and the grid size (currently is %d) that is determined by the product of \"dim_size1\" and \"dim_size2\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lon, dim_size1, grid_size, annotation);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, size_center_lat == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lat\" is different from the grid size (currently is %d) that is determined by the product of \"dim_size1\" and \"dim_size2\" (the array size of \"center_lon\" is the same as the grid size). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lat, grid_size, annotation);
        }    
    }
    if (size_mask > 0)        
        EXECUTION_REPORT(REPORT_ERROR, comp_id, size_mask == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"mask\" is different from the grid size (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_mask, grid_size, annotation);
    if (size_area > 0)        
        EXECUTION_REPORT(REPORT_ERROR, comp_id, size_area == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"area\" is different from the grid size (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_area, grid_size, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, (size_vertex_lon == -1 && size_vertex_lat == -1) || (size_vertex_lon > 0 && size_vertex_lat > 0), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the optional parameters \"vertex_lon\" and \"vertex_lat\" are not set/unset at the same time. Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);
    if (size_vertex_lon > 0) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, (size_vertex_lon) % (size_center_lon) == 0, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size of \"vertex_lon\" (currently is %d) is not an integral multiple of the array size (currently is %d) of \"center_lon\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_vertex_lon, size_center_lon, annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, (size_vertex_lat) % (size_center_lat) == 0, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size of \"vertex_lat\" (currently is %d) is not an integral multiple of the array size (currently is %d) of \"center_lat\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_vertex_lat, size_center_lat, annotation);
        int num_vertex_lon = (size_vertex_lon) / (size_center_lon);
        int num_vertex_lat = (size_vertex_lat) / (size_center_lat);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, num_vertex_lon == num_vertex_lat, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the numbers of vertexes determined by \"vertex_lon\" (currently is %d) and \"vertex_lat\" (currently is %d) respectively are not the same. Please check the model code related to the annotation \"%s\".", grid_name, API_label, num_vertex_lon, num_vertex_lat, annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, size_vertex_lon > grid_size && size_vertex_lat > grid_size || size_vertex_lon < grid_size && size_vertex_lat < grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array sizes of \"vertex_lon\" and \"vertex_lat\" (currently are %d and %d respectively) are not bigger/smaller than the grid size at the same time. Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_vertex_lon, size_vertex_lat, annotation);
        if (size_center_lon == grid_size) 
            EXECUTION_REPORT(REPORT_ERROR, comp_id, num_vertex_lon >= 3, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the number of vertexes (currently is %d) is wrong as it is smaller than 3. Please check the model code related to the annotation \"%s\".", grid_name, API_label, num_vertex_lon, annotation);
        else EXECUTION_REPORT(REPORT_ERROR, comp_id ,num_vertex_lon == 2, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the number of vertexes (currently is %d) is wrong as it is not 2. Please check the model code related to the annotation \"%s\".", grid_name, API_label, num_vertex_lon, annotation);
        num_vertex = num_vertex_lon;
    }

    return create_H2D_grid_from_global_data(comp_id, grid_name, coord_unit, cyclic_or_acyclic, data_type, dim_size1, dim_size2, num_vertex, size_center_lon, size_center_lat, 
                                            size_mask, size_area, size_vertex_lon, size_vertex_lat, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation);
}


int Original_grid_mgt::register_H2D_grid_via_file(int comp_id, const char *grid_name, const char *data_file_name, const char *annotation)
{
    int rcode, ncfile_id, grid_id;
    int size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat;
    long dim_lon_size, dim_lat_size, dim_H2D_size, dim_size1, dim_size2;
    char *center_lon, *center_lat, *vertex_lon, *vertex_lat, *area;
    char min_lon[NAME_STR_SIZE], max_lon[NAME_STR_SIZE], min_lat[NAME_STR_SIZE], max_lat[NAME_STR_SIZE];
    int *mask;
    char data_type_for_center_lat[NAME_STR_SIZE], data_type_for_center_lon[NAME_STR_SIZE], data_type_for_vertex_lon[NAME_STR_SIZE], data_type_for_vertex_lat[NAME_STR_SIZE], data_type_for_mask[NAME_STR_SIZE], data_type_for_area[NAME_STR_SIZE];
    char data_type_temp[NAME_STR_SIZE];
    char edge_type[NAME_STR_SIZE], cyclic_or_acyclic[NAME_STR_SIZE], unit_center_lon[NAME_STR_SIZE], unit_center_lat[NAME_STR_SIZE], unit_vertex_lon[NAME_STR_SIZE], unit_vertex_lat[NAME_STR_SIZE];
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_H2D_grid_via_file");
    bool is_root_proc = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in register_H2D_grid_via_file") == 0;
    

    check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE, comm, "registering an H2D grid", data_file_name, "data_file_name", annotation);

    IO_netcdf *netcdf_file_object = new IO_netcdf("H2D_grid_data", data_file_name, "r", false);
    dim_lon_size = netcdf_file_object->get_dimension_size(COORD_LABEL_LON, comm, is_root_proc);
    dim_lat_size = netcdf_file_object->get_dimension_size(COORD_LABEL_LAT, comm, is_root_proc);
    dim_H2D_size = netcdf_file_object->get_dimension_size("grid_size", comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_CENTER_LON_LABEL, (void**)(&center_lon), &size_center_lon, data_type_for_center_lon, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_CENTER_LAT_LABEL, (void**)(&center_lat), &size_center_lat, data_type_for_center_lat, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_VERTEX_LON_LABEL, (void**)(&vertex_lon), &size_vertex_lon, data_type_for_vertex_lon, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_VERTEX_LAT_LABEL, (void**)(&vertex_lat), &size_vertex_lat, data_type_for_vertex_lat, comm, is_root_proc);
    netcdf_file_object->read_file_field("area", (void**)(&area), &size_area, data_type_for_area, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_MASK_LABEL, (void**)(&mask), &size_mask, data_type_for_mask, comm, is_root_proc);
    if (dim_lon_size > 0 && dim_lat_size > 0 && dim_H2D_size > 0)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, dim_H2D_size == dim_lon_size*dim_lat_size, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the size (currently is %d) of dimension \"grid_size\" is different from the multiple (currently is %d) of sizes of dimensions \"lon\" and \"lat\"", grid_name, annotation, dim_H2D_size, dim_lon_size*dim_lat_size, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, dim_H2D_size > 0 || (dim_lon_size > 0 && dim_lat_size > 0), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: the dimension size (dimensions \"lon\" and \"lat\" in the file) or the grid size (dimension \"grid_size\" in the file) is not correctly specified in the file \"%s\". Please verify.", grid_name, annotation, data_file_name);
    if (dim_lon_size > 0 && dim_lat_size > 0) {
        dim_size1 = dim_lon_size;
        dim_size2 = dim_lat_size;
    }
    else {
        dim_size1 = dim_H2D_size;
        dim_size2 = 0;
    }

    EXECUTION_REPORT(REPORT_ERROR, comp_id, center_lon != NULL, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: the longitude value for the center of each grid point (variable \"lon\" in the file) is not specified in the data file \"%s\". ", 
                     grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, center_lat != NULL, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: the latitude value for the center of each grid point (variable \"lat\" in the file) is not specified in the data file \"%s\". ", 
                     grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, vertex_lon != NULL && vertex_lat != NULL || vertex_lon == NULL && vertex_lat == NULL, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the longitude and latitude values for each vertex (variables \"vertex_lon\" and \"vertex_lat\" in the file) of each grid point must be specified/unspecified at the same time", 
                     grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon,data_type_for_center_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variables \"lon\" and \"lat\" are not the same", grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon,DATA_TYPE_FLOAT) || words_are_the_same(data_type_for_center_lon,DATA_TYPE_DOUBLE), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variables \"lon\" is not floating-point", grid_name, annotation, data_file_name);
    if (vertex_lon != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon,data_type_for_vertex_lon), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"vertex_lon\" is different from the data type of variable \"lon\".", grid_name, annotation, data_file_name);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon,data_type_for_vertex_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"vertex_lat\" is different from the data type of variable \"lat\".", grid_name, annotation, data_file_name);
    }
    if (area != NULL)        
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon,data_type_for_area), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"area\" is different from the data type of variable \"lon\".", grid_name, annotation, data_file_name);
    if (mask != NULL)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_mask, DATA_TYPE_INT), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"mask\" is not \"integer\".", grid_name, annotation, data_file_name);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "edge_type", edge_type, data_type_temp, comm, is_root_proc) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the global attribute \"edge_type\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "cyclic_or_acyclic", cyclic_or_acyclic, data_type_temp, comm, is_root_proc) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the global attribute \"cyclic_or_acyclic\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(edge_type) > 0, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", \"edge_type\" is not specified as a global attribute", grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(cyclic_or_acyclic) > 0, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", \"cyclic_or_acyclic\" is not specified as a global attribute", grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, (netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LON_LABEL, "unit", unit_center_lon, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LON_LABEL, "units", unit_center_lon, data_type_temp, comm, is_root_proc)) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"%s\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, SCRIP_CENTER_LON_LABEL, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, (netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LAT_LABEL, "unit", unit_center_lat, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LAT_LABEL, "units", unit_center_lat, data_type_temp, comm, is_root_proc)) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"%s\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, SCRIP_CENTER_LAT_LABEL, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "min_lon", min_lon, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"min_lon\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);
    transform_datatype_of_arrays(min_lon, min_lon, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "max_lon", max_lon, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"max_lon\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);
    transform_datatype_of_arrays(max_lon, max_lon, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "min_lat", min_lat, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"min_lat\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);    
    transform_datatype_of_arrays(min_lat, min_lat, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "max_lat", max_lat, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"max_lat\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);
    transform_datatype_of_arrays(max_lat, max_lat, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(unit_center_lon,unit_center_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the units of variables \"lon\" and \"lat\" are different", grid_name, annotation, data_file_name);
    if (vertex_lon != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LON_LABEL, "unit", unit_vertex_lon, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LON_LABEL, "units", unit_vertex_lon, data_type_temp, comm, is_root_proc), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"vertex_lon\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LAT_LABEL, "unit", unit_vertex_lat, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LAT_LABEL, "units", unit_vertex_lat, data_type_temp, comm, is_root_proc), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"vertex_lat\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(unit_center_lon,unit_vertex_lon), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the units of variables \"lon\" and \"vertex_lon\" are different", grid_name, annotation, data_file_name);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(unit_center_lat,unit_vertex_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the units of variables \"lat\" and \"vertex_lat\" are different", grid_name, annotation, data_file_name);
    }

    grid_id = register_H2D_grid_via_global_data(comp_id, grid_name, edge_type, unit_center_lon, cyclic_or_acyclic, data_type_for_center_lon, dim_size1, dim_size2,size_center_lon, size_center_lat,
                                         size_mask, size_area, size_vertex_lon, size_vertex_lat, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation, API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE);

    delete [] center_lon;
    delete [] center_lat;
    if (vertex_lon != NULL) {
        delete [] vertex_lon;
        delete [] vertex_lat;
    }
    if (mask != NULL)
        delete [] mask;
    if (area != NULL)
        delete [] area;

    return grid_id;
}


int Original_grid_mgt::register_V1D_grid_via_data(int API_id, int comp_id, const char *grid_name, int grid_type, const char *coord_unit, int grid_size, 
                                                  double value1, const double *value2, const double *value3, const char *annotation)
{
    char full_grid_name[NAME_STR_SIZE];
    Remap_grid_class *CoR_V1D_grid;
    
    
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_V1D_grid_via_data"), NULL, grid_size, "implicit grid size", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_V1D_grid_via_data"), "registering a V1D grid", 1, sizeof(double), (const char*)(&value1), "floating-point parameters", annotation);
	if (value2 != NULL)
	    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_V1D_grid_via_data"), "registering a V1D grid", grid_size, sizeof(double), (const char*)(value2), "floating-point parameters", annotation);
	if (value3 != NULL)
	    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_V1D_grid_via_data"), "registering a V1D grid", grid_size, sizeof(double), (const char*)(value3), "floating-point parameters", annotation);

    sprintf(full_grid_name, "%s@%s", grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"in register_V1D_grid_via_data")->get_full_name());
	if (grid_type == 0)
		CoR_V1D_grid = new Remap_grid_class(full_grid_name, COORD_LABEL_TRACER, "", NULL, grid_size);
	else if (grid_type <= 4)
		CoR_V1D_grid = new Remap_grid_class(full_grid_name, COORD_LABEL_LEV, coord_unit, NULL, grid_size);
	switch(grid_type) {
		case 1:
			CoR_V1D_grid->read_grid_data_from_array("center", COORD_LABEL_LEV, DATA_TYPE_DOUBLE, (const char*)value2, 0);
			break;
		case 2:
			CoR_V1D_grid->set_lev_grid_sigma_info(value1, value2, NULL, 1.0);
			break;
		case 3:
			CoR_V1D_grid->set_lev_grid_sigma_info(value1, value2, value3, 1.0);
			break;
		default:
			break;
	}

    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, CoR_V1D_grid, true));

    remap_grid_manager->add_remap_grid(CoR_V1D_grid);
    
    return original_grids[original_grids.size()-1]->get_grid_id();    
}


int Original_grid_mgt::register_md_grid_via_multi_grids(int comp_id, const char *grid_name, int sub_grid1_id, int sub_grid2_id, int sub_grid3_id, int size_mask, int *mask, const char *annotation)
{
    int num_sub_grids = 2, num_H2D_grid = 0, num_V1D_grid = 0, num_T1D_grid = 0;
    char full_grid_name[NAME_STR_SIZE];
    Remap_grid_class *CoR_MD_grid, *sub_grids[3];

    
    MPI_Comm local_comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Original_grid_mgt::register_md_grid_via_multi_grids");
    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(sub_grid1_id), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": \"sub_grid1_id\" (currently 0x%x) is wrong. Please check the model code related to the annotation \"%s\".", grid_name, sub_grid1_id, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_comp_id_of_grid(sub_grid1_id) == comp_id, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the component model corresponding to the grid with id of \"sub_grid1_id\" is different from the current component with the id of comp_id. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, "for registering a multi-dimension grid", original_grid_mgr->get_name_of_grid(sub_grid1_id), "\"sub_grid1_id\" (the corresponding grid name)", annotation);        
    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(sub_grid2_id), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": \"sub_grid2_id\" is wrong. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_comp_id_of_grid(sub_grid2_id) == comp_id, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the component corresponding to the grid with id of \"sub_grid2_id\" is different from the current component with the id of comp_id. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    check_API_parameter_data_array(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in register_h2d_grid_with_data"), "registering a multi-dimension grid", size_mask, sizeof(int), (const char*)mask, "mask", annotation);    
    check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, "for registering a multi-dimension grid", original_grid_mgr->get_name_of_grid(sub_grid2_id), "\"sub_grid2_id\" (the corresponding grid name)", annotation);
    sub_grids[0] = original_grid_mgr->get_original_CoR_grid(sub_grid1_id);
    sub_grids[1] = original_grid_mgr->get_original_CoR_grid(sub_grid2_id);
    int have_sub_grid3 = sub_grid3_id != -1? 1: 0;
    check_API_parameter_int(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, NULL, have_sub_grid3, "\"sub_grid3_id\" (specified or not)", annotation);
    if (sub_grid3_id != -1) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(sub_grid3_id), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": \"sub_grid3_id\" is wrong. Please check the model code related to the annotation \"%s\".", grid_name, annotation);        
        check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, "for registering a multi-dimension grid", original_grid_mgr->get_name_of_grid(sub_grid3_id), "\"sub_grid3_id\" (the corresponding grid name)", annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_comp_id_of_grid(sub_grid3_id) == comp_id, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the component corresponding to the grid with id of \"sub_grid3_id\" is different from the current component with the id of comp_id. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
        sub_grids[num_sub_grids++] = original_grid_mgr->get_original_CoR_grid(sub_grid3_id);
    }    

    for (int i = 0; i < num_sub_grids; i ++)
        if (sub_grids[i]->get_is_sphere_grid())
            num_H2D_grid ++;
        else if (sub_grids[i]->has_grid_coord_label(COORD_LABEL_LEV))
            num_V1D_grid ++;
        else if (sub_grids[i]->has_grid_coord_label(COORD_LABEL_TIME))
            num_T1D_grid ++;
        else EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": grid \"%s\" cannot be used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, sub_grids[i]->get_grid_name(), annotation);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_H2D_grid <= 1, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": more than one H2D grid are used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_V1D_grid <= 1, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": more than one V1D grid are used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_T1D_grid <= 1, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": more than one T1D grid are used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    
    sprintf(full_grid_name, "%s@%s", grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"register_md_grid_via_multi_grids")->get_full_name());
    CoR_MD_grid = new Remap_grid_class(full_grid_name, num_sub_grids, sub_grids, 0);
    if (size_mask > 0)    {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, size_mask == CoR_MD_grid->get_grid_size(), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the array size (currently is %d) of \"mask\" is different from the grid size (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, size_mask, CoR_MD_grid->get_grid_size(), annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, are_array_values_between_boundaries("integer", mask, size_mask, 0, 1, 0, false), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": some values of the parameter \"mask\" are wrong (not 0 and 1). Please check the model code related to the annotation \"%s\".", grid_name, annotation);
        CoR_MD_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)mask, 0);
    }    

    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, CoR_MD_grid, true));    
    remap_grid_manager->add_remap_grid(CoR_MD_grid);
    
    return original_grids[original_grids.size()-1]->get_grid_id();
}


void Original_grid_mgt::set_3D_grid_3D_vertical_coord_field_inst(int grid_id, int field_id, const char *static_or_dynamic, const char *annotation)
{
	char API_label[NAME_STR_SIZE];
	Original_grid_info *original_grid;
	Field_mem_info *field_inst;
	int comp_id;
	MPI_Comm local_comm;


	get_API_hint(-1, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, API_label); 
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id), "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of a 3-D grid: the parameter of \"grid_id\" is wrong. Please verify the model code with the annotation \"%s\".", annotation);
	original_grid = original_grid_mgr->get_original_grid(grid_id);
	comp_id = original_grid->get_comp_id();
    local_comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Original_grid_mgt::set_3D_grid_3D_vertical_coord_field_inst");
    check_API_parameter_string(comp_id, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, local_comm, "setting the 3-D vertical coordinate field of a 3-D grid", original_grid->get_grid_name(), "\"grid_id\" (the name of the 3-D grid)", annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid->is_3D_grid(), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": this grid is not a 3-D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !original_grid->get_original_CoR_grid()->is_sigma_grid(), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": cannot set a 3-D vertical field to this grid because its V1D sub grid is a SIGMA or HYBRID grid. Please verify the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !original_grid->get_original_CoR_grid()->does_use_V3D_level_coord(), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": cannot set a 3-D vertical field to this grid because its vertical coordinate values have already been set. Please verify the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, memory_manager->check_is_legal_field_instance_id(field_id), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": the parameter of \"field_id\" is wrong. Please verify the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
	check_API_parameter_field_instance(comp_id, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, local_comm, "setting the 3-D vertical coordinate field of a 3-D grid", field_id, "\"field_id\" (the corresponding 3-D field)", annotation);
	field_inst = memory_manager->get_field_instance(field_id);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, field_inst->get_grid_id() == grid_id, "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": the given field instance is not on the given 3-D grid. Please verify the model code with the annotation \"%s\".", original_grid->get_grid_name(), annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(original_grid->get_original_CoR_grid()->get_a_leaf_grid(COORD_LABEL_LEV)->get_coord_unit(), field_inst->get_unit()), "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": The unit (\"%s\") of the grid coordinate is different from the unit (\"%s\") of the given field instance (\"%s\"). Please verify the model code with the annotation \"%s\".", original_grid->get_grid_name(), original_grid->get_original_CoR_grid()->get_a_leaf_grid(COORD_LABEL_LEV)->get_coord_unit(), field_inst->get_unit(), field_inst->get_field_name(), annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(static_or_dynamic, "constant") || words_are_the_same(static_or_dynamic, "variable"), "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": the input parameter of \"label\" (the current value is \"%s\") must be \"constant\" or \"variable\". Please verify the model code with the annotation \"%s\".", original_grid->get_grid_name(), static_or_dynamic, annotation);
    check_API_parameter_string(comp_id, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, local_comm, "setting the 3-D vertical coordinate field of a 3-D grid", original_grid->get_grid_name(), "\"label\" (specification for constant or variable of 3-D field)", annotation);
    if (original_grid->get_V3D_lev_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": the 3-D vertical coordinate field has been set before at the model code with the annotation \"%s\" and cannot be set again at the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation_mgr->get_annotation(grid_id, "set 3-D level field"), annotation);
	original_grid->set_unique_3D_lev_field(field_id, static_or_dynamic, annotation);	
	original_grid->get_original_CoR_grid()->set_using_V3D_level_coord();
	decomp_grids_mgr->set_decomp_grids_using_3D_level_coord(original_grid->get_original_CoR_grid());
}


void Original_grid_mgt::set_3d_grid_bottom_field(int comp_id, int grid_id, int field_id, int static_or_dynamic_or_external, int API_id, const char *API_label, const char *annotation)
{
    Original_grid_info *original_grid = original_grid_mgr->get_original_grid(grid_id);
    Original_grid_info *field_original_grid;
    MPI_Comm local_comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Original_grid_mgt::set_3d_grid_bottom_field");


    synchronize_comp_processes_for_API(comp_id, API_id, local_comm, "setting the surface field of a 3-D grid", annotation);
    check_API_parameter_string(comp_id, API_id, local_comm, "setting the surface field of a 3-D grid", original_grid->get_grid_name(), "\"grid_id\" (the name of the 3-D grid)", annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid->is_3D_grid(), "Error happens when calling the API \"%s\" to set the surface field of the grid \"%s\": this grid is not a 3-D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid->get_original_CoR_grid()->is_sigma_grid(), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": cannot set the surface field to this grid because its V1D sub grid is not a SIGMA or HYBRID grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    if (original_grid->get_bottom_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the surface field has been set before at the model code with the annotation \"%s\" and cannot be set again at the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation_mgr->get_annotation(grid_id, "set surface field"), annotation);
    if (static_or_dynamic_or_external != BOTTOM_FIELD_VARIATION_EXTERNAL) {
        check_API_parameter_field_instance(comp_id, API_id, local_comm, "setting the surface field of a 3-D grid", field_id, "\"field_id\" (the corresponding surface field)", annotation);
        Field_mem_info *field_inst = memory_manager->get_field_instance(field_id);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, field_inst->get_grid_id() != -1, "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the surface field \"%s\" is not on a grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
        field_original_grid = original_grid_mgr->get_original_grid(field_inst->get_grid_id());
        EXECUTION_REPORT(REPORT_ERROR, comp_id, field_original_grid->is_H2D_grid(), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the grid \"%s\" corresponding to the surface field \"%s\" is not a H2D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, field_original_grid->get_original_CoR_grid()->is_subset_of_grid(original_grid->get_original_CoR_grid()), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the grid \"%s\" corresponding to the surface field \"%s\" is not a sub grid of the 3-D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(field_inst->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT) || words_are_the_same(field_inst->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the data type of the surface field \"%s\" is not floating-point. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
    }
    original_grid->set_unique_bottom_field(field_id, static_or_dynamic_or_external, annotation);
    if (original_grid->get_mid_point_grid() != NULL)
        original_grid->get_mid_point_grid()->set_unique_bottom_field(field_id, static_or_dynamic_or_external, annotation);
    if (original_grid->get_interface_level_grid() != NULL)
        original_grid->get_interface_level_grid()->set_unique_bottom_field(field_id, static_or_dynamic_or_external, annotation);
}


int Original_grid_mgt::get_CoR_defined_grid(int comp_id, const char *grid_name, const char *CoR_grid_name, const char *annotation)
{
    Original_grid_info *original_grid;
    Remap_grid_class *original_CoR_grid;
    char CoR_script_file_name[NAME_STR_SIZE];

    
    sprintf(CoR_script_file_name, "%s/CCPL_grid.cor", comp_comm_group_mgt_mgr->get_root_comp_config_dir());
    original_CoR_grid = remap_grid_manager->search_remap_grid_with_grid_name(CoR_grid_name);
    if (original_CoR_grid == NULL)
        if (strlen(CoR_script_name) > 0)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not define the corresponding CoR grid. Please check the CoR script or the model code corresponding to the annotation \"%s\"",
                             grid_name, CoR_grid_name, CoR_script_name, annotation);    
        else EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not exist. Please verify", 
                              grid_name, CoR_grid_name, CoR_script_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, original_CoR_grid->format_sub_grids(original_CoR_grid), "Please modify the definition of grid \"%s\" in the CoR script \"%s\". We propose to order the dimensions of the grid into the order such as lon, lat, level and time");
    original_CoR_grid->end_grid_definition_stage(NULL);
    if (original_CoR_grid->get_is_sphere_grid() && are_floating_values_equal(NULL_COORD_VALUE, original_CoR_grid->get_boundary_min_lon()))
        original_CoR_grid->set_grid_boundary(0, 360.0, -90.0, 90.0);
    original_grid = new Original_grid_info(comp_id, original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, original_CoR_grid, true);
    original_grids.push_back(original_grid);
    if (original_grid->get_H2D_sub_CoR_grid() != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_CoR_grid->get_grid_center_field(COORD_LABEL_LON) != NULL, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not specify the center values of longitude (X) of the CoR grid. Please verify.", grid_name, CoR_grid_name, CoR_script_file_name);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_CoR_grid->get_grid_center_field(COORD_LABEL_LAT) != NULL, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not specify the center values of latitude (Y) of the CoR grid. Please verify.", grid_name, CoR_grid_name, CoR_script_file_name);
    }
    if (original_grid->get_V1D_sub_CoR_grid() != NULL) {
        Remap_grid_data_class *level_field = original_CoR_grid->get_grid_center_field(COORD_LABEL_LEV);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, level_field != NULL && level_field->get_coord_value_grid()->get_grid_size() == original_grid->get_V1D_sub_CoR_grid()->get_grid_size(),
                         "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not specify the center values of the vertical coordinate of the CoR grid or the vertical coordinate is not a Z grid. Please verify.", grid_name, CoR_grid_name, CoR_script_file_name);
    }
    return original_grid->get_local_grid_id();
}


bool Original_grid_mgt::is_grid_id_legal(int grid_id) const
{
    int true_grid_id = grid_id & TYPE_ID_SUFFIX_MASK;


    if ((grid_id & TYPE_ID_PREFIX_MASK) != TYPE_GRID_LOCAL_ID_PREFIX)
        return false;

    if (true_grid_id < 0 || true_grid_id >= original_grids.size())
        return false;

    return true;
}


int Original_grid_mgt::get_comp_id_of_grid(int grid_id) const
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id&TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_comp_id_of_grid");
    return original_grids[grid_id&TYPE_ID_SUFFIX_MASK]->get_comp_id();
}


const char *Original_grid_mgt::get_name_of_grid(int grid_id) const
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id&TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_name_of_grid");
    return original_grids[grid_id&TYPE_ID_SUFFIX_MASK]->get_grid_name();
}


Remap_grid_class *Original_grid_mgt::get_original_CoR_grid(int grid_id) const
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id&TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_original_CoR_grid");
    return original_grids[grid_id&TYPE_ID_SUFFIX_MASK]->get_original_CoR_grid();
}


Original_grid_info *Original_grid_mgt::get_original_grid(int grid_id) const
{
    if (!is_grid_id_legal(grid_id))
        return NULL;
    
    return original_grids[grid_id&TYPE_ID_SUFFIX_MASK];
}

int Original_grid_mgt::get_grid_size(int grid_id, const char *annotation) const
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_grid_id_legal(grid_id), "Error happens when getting the size of a grid: the parameter of grid ID is wrong. Please verify the model code corresponding to the annotation \"%s\"", annotation);        
    return get_original_CoR_grid(grid_id)->get_grid_size();
}


int Original_grid_mgt::get_grid_id(int comp_id, const char *grid_name, const char *annotation)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id,false), "Error happens when getting the ID of a grid: the parameter of component ID is wrong. Please verify the model code corresponding to the annotation \"%s\"", annotation);
    Original_grid_info *original_grid = search_grid_info(grid_name, comp_id);
    if (original_grid == NULL)
        return -1;
    else return original_grid->get_grid_id();
}


int Original_grid_mgt::add_original_grid(int comp_id, const char *grid_name, Remap_grid_class *original_CoR_grid)
{
    Original_grid_info *existing_grid = search_grid_info(grid_name, comp_id);
    if (existing_grid != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_CoR_grid == existing_grid->get_original_CoR_grid(), "Software error in Original_grid_mgt::add_original_grid");
        return existing_grid->get_local_grid_id();
    }

    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, grid_name, "Original_grid_mgt::add_original_grid", original_CoR_grid, false));
    return original_grids[original_grids.size()-1]->get_local_grid_id();
}


int Original_grid_mgt::get_num_grid_levels(int grid_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id&TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_num_grid_levels: wrong grid id");        
    if (original_grids[grid_id&TYPE_ID_SUFFIX_MASK]->get_V1D_sub_CoR_grid() == NULL)
        return 1;

    EXECUTION_REPORT(REPORT_ERROR, -1, original_grids[grid_id&TYPE_ID_SUFFIX_MASK]->get_V1D_sub_CoR_grid()->get_grid_size() > 0, "Software error in Original_grid_mgt::get_num_grid_levels: wrong size of vertical grid");
    return original_grids[grid_id&TYPE_ID_SUFFIX_MASK]->get_V1D_sub_CoR_grid()->get_grid_size();
}


bool Original_grid_mgt::is_V1D_sub_grid_after_H2D_sub_grid(int grid_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id&TYPE_ID_SUFFIX_MASK] != NULL, "Original_grid_info::is_V1D_sub_grid_after_H2D_sub_grid: wrong grid id");        
    return original_grids[grid_id&TYPE_ID_SUFFIX_MASK]->is_V1D_sub_grid_after_H2D_sub_grid();
}


void Original_grid_mgt::register_mid_point_grid(int level_3D_grid_id, int *mid_3D_grid_id, int *mid_1D_grid_id, int size_mask, const int *mask, const char *annotation, const char *API_label)
{
    Original_grid_info *level_3D_grid, *mid_3D_grid, *mid_1D_grid;
    Remap_grid_class *level_1D_CoR_grid, *level_H2D_CoR_grid, *mid_1D_CoR_grid, *mid_3D_CoR_grid, *sub_grids[256], *existing_level_1D_CoR_grid;
    char grid_name[NAME_STR_SIZE];


    level_3D_grid = get_original_grid(level_3D_grid_id);
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(level_3D_grid->get_comp_id(), "in Original_grid_mgt::register_mid_point_grid");
    synchronize_comp_processes_for_API(level_3D_grid->get_comp_id(), API_ID_GRID_MGT_REG_MID_POINT_GRID, comm, "registering a mid-point grid", annotation);
    check_API_parameter_string(level_3D_grid->get_comp_id(), API_ID_GRID_MGT_REG_MID_POINT_GRID, comm, "registering a mid-point grid", level_3D_grid->get_grid_name(), "\"level_3D_grid_id\" (the name of the interface-level grid)", annotation);
    EXECUTION_REPORT(REPORT_ERROR, level_3D_grid->get_comp_id(), level_3D_grid->is_3D_grid(), "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the grid specified by the parameter \"level_3D_grid_id\" is not a 3-D grid. Please verify the model code with the annotation \"%s\".", API_label, annotation);
    if (level_3D_grid->get_mid_point_grid() != NULL)    
        EXECUTION_REPORT(REPORT_ERROR, level_3D_grid->get_comp_id(), false, "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the mid-point grid of the grid \"%s\" has been registered before (at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\".", API_label, level_3D_grid->get_grid_name(), level_3D_grid->get_mid_point_grid()->get_annotation(), annotation);
    if (level_3D_grid->get_interface_level_grid() != NULL)
        EXECUTION_REPORT(REPORT_ERROR, level_3D_grid->get_comp_id(), false, "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the specified interface-level grid itself is a mid-point grid (registered at the model code with the annotation \"%s\") of the grid \"%s\". It cannot be used to register another mid-point grid. Please verify the model code with the annotation \"%s\".", API_label, level_3D_grid->get_annotation(), level_3D_grid->get_interface_level_grid()->get_grid_name(), annotation);
    check_API_parameter_data_array(level_3D_grid->get_comp_id(), API_ID_GRID_MGT_REG_MID_POINT_GRID, comm, "registering a mid-point grid", size_mask, sizeof(int), (const char*)mask, "\"mask\"", annotation);
    level_1D_CoR_grid = level_3D_grid->get_V1D_sub_CoR_grid();
    level_H2D_CoR_grid = level_3D_grid->get_H2D_sub_CoR_grid(); 
    EXECUTION_REPORT(REPORT_ERROR, -1, level_1D_CoR_grid != NULL && level_H2D_CoR_grid != NULL, "Software error in Original_grid_mgt::register_mid_point_grid: NULL level_1D_CoR_grid or NULL level_H2D_CoR_grid");
	existing_level_1D_CoR_grid = level_1D_CoR_grid->get_mid_point_grid();
    mid_1D_CoR_grid = level_1D_CoR_grid->generate_mid_point_grid();
    sub_grids[0] = level_H2D_CoR_grid;
    sub_grids[1] = mid_1D_CoR_grid;
    sprintf(grid_name, "mid_grid_for_%s", level_3D_grid->get_original_CoR_grid()->get_grid_name());
	mid_3D_CoR_grid = level_3D_grid->get_original_CoR_grid()->get_mid_point_grid();
	if (mid_3D_CoR_grid == NULL) {
	    mid_3D_CoR_grid = new Remap_grid_class(grid_name, 2, sub_grids, 0);
    	remap_grid_manager->add_remap_grid(mid_3D_CoR_grid);
	}
    mid_1D_CoR_grid->end_grid_definition_stage(NULL);
    mid_3D_CoR_grid->end_grid_definition_stage(NULL);
    mid_1D_grid = search_grid_info(mid_1D_CoR_grid->get_grid_name(), level_3D_grid->get_comp_id());
    if (mid_1D_grid == NULL) {
        mid_1D_grid = new Original_grid_info(level_3D_grid->get_comp_id(), original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, mid_1D_CoR_grid->get_grid_name(), annotation, mid_1D_CoR_grid, true);
        original_grids.push_back(mid_1D_grid);
		if (existing_level_1D_CoR_grid == NULL)
	        remap_grid_manager->add_remap_grid(mid_1D_grid->get_original_CoR_grid());
    }
    *mid_1D_grid_id = mid_1D_grid->get_grid_id();
    mid_3D_grid = search_grid_info(mid_3D_CoR_grid->get_grid_name(), level_3D_grid->get_comp_id());
    EXECUTION_REPORT(REPORT_ERROR, -1, mid_3D_grid == NULL, "Software error in Original_grid_mgt::register_mid_point_grid: mid_3D_grid exists");
    mid_3D_grid = new Original_grid_info(level_3D_grid->get_comp_id(), original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, mid_3D_CoR_grid->get_grid_name(), annotation, mid_3D_CoR_grid, true);
    original_grids.push_back(mid_3D_grid);
    *mid_3D_grid_id = mid_3D_grid->get_grid_id();
    level_3D_grid->set_mid_point_grid(mid_3D_grid);
    if (size_mask > 0) {
        EXECUTION_REPORT(REPORT_ERROR, level_3D_grid->get_comp_id(), size_mask == mid_3D_CoR_grid->get_grid_size(), "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the array size of \"mask\" is different from the size of the mid-point grid. Please verify the model code with the annotation \"%s\".", API_label, annotation);
        EXECUTION_REPORT(REPORT_ERROR, level_3D_grid->get_comp_id(), are_array_values_between_boundaries("integer", mask, size_mask, 0, 1, 0, false), "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: some values of the parameter \"mask\" are wrong (not 0 and 1). Please check the model code related to the annotation \"%s\".", API_label, annotation);
        mid_3D_CoR_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)mask, 0);
    }
    if (level_1D_CoR_grid->get_super_grid_of_setting_coord_values() == level_1D_CoR_grid)
        mid_1D_CoR_grid->set_super_grid_of_setting_coord_values(mid_1D_CoR_grid);
	else if (level_1D_CoR_grid->get_super_grid_of_setting_coord_values() != NULL && level_1D_CoR_grid->get_super_grid_of_setting_coord_values()->get_num_dimensions() == 3)
		mid_1D_CoR_grid->set_super_grid_of_setting_coord_values(mid_3D_CoR_grid);
	level_3D_grid->get_original_CoR_grid()->set_mid_point_grid(mid_3D_CoR_grid);
}


Original_grid_info *Original_grid_mgt::search_or_register_internal_grid(int comp_id, Remap_grid_class *CoR_grid)
{
	for (int i = 0; i < original_grids.size(); i ++)
		if (original_grids[i]->get_comp_id() == comp_id && CoR_grid == original_grids[i]->get_original_CoR_grid())
			return original_grids[i];
	original_grids.push_back(new Original_grid_info(comp_id, original_grids.size()|TYPE_GRID_LOCAL_ID_PREFIX, CoR_grid->get_grid_name(), "Original_grid_mgt::add_original_grid", CoR_grid, false));
	return original_grids[original_grids.size()-1];
}


