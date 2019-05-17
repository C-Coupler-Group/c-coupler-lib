/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "remapping_configuration_mgt.h"
#include "CCPL_api_mgt.h"


H2D_remapping_wgt_file_info::H2D_remapping_wgt_file_info(const char *wgt_file_name)
{    
    this->wgt_file_name = strdup(wgt_file_name);
    num_wgts = 0;
    wgts_src_indexes = NULL;
    wgts_dst_indexes = NULL;
    wgts_values = NULL;
    src_center_lon = NULL;
    src_center_lat = NULL;
    src_area = NULL;
    dst_center_lon = NULL;
    dst_center_lat = NULL;
    dst_area = NULL;
    checksum_src_mask = -1;
    checksum_dst_mask = -1;
    src_grid_size = -1;
    dst_grid_size = -1;
}


H2D_remapping_wgt_file_info::H2D_remapping_wgt_file_info(const char *array, long *buffer_content_iter)
{    
    long str_size; 


    wgt_file_name = load_string(NULL, str_size, 4096, array, *buffer_content_iter, NULL);
    read_data_from_array_buffer(&checksum_dst_mask, sizeof(long), array, *buffer_content_iter, true);
    read_data_from_array_buffer(&dst_grid_size, sizeof(int), array, *buffer_content_iter, true);
    read_data_from_array_buffer(&checksum_src_mask, sizeof(long), array, *buffer_content_iter, true);
    read_data_from_array_buffer(&src_grid_size, sizeof(int), array, *buffer_content_iter, true);
    
    num_wgts = 0;
    wgts_src_indexes = NULL;
    wgts_dst_indexes = NULL;
    wgts_values = NULL;    
    src_center_lon = NULL;
    src_center_lat = NULL;
    src_area = NULL;
    dst_center_lon = NULL;
    dst_center_lat = NULL;
    dst_area = NULL;
    checksum_src_mask = -1;
    checksum_dst_mask = -1;
    src_grid_size = -1;
    dst_grid_size = -1;
}


void H2D_remapping_wgt_file_info::write_remapping_wgt_file_info_into_array(char **array, long &buffer_max_size, long &buffer_content_size)
{
    write_data_into_array_buffer(&src_grid_size, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&checksum_src_mask, sizeof(long), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&dst_grid_size, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&checksum_dst_mask, sizeof(long), array, buffer_max_size, buffer_content_size);
    dump_string(wgt_file_name, -1, array, buffer_max_size, buffer_content_size);
}


long H2D_remapping_wgt_file_info::get_grid_field_checksum_value(const char *field_name, IO_netcdf *netcdf_file_object, int grid_size, MPI_Comm comm, bool is_root_proc)
{
    char *data_buffer, data_type[NAME_STR_SIZE];
    int field_size;
    long checksum;


    netcdf_file_object->read_file_field(field_name, (void**)(&data_buffer), &field_size, data_type, comm, is_root_proc);
    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), field_size > 0, "Error happens when reading the remapping weights file \"%s\": variable \"%s\" does not exist in the file. Please verify.", wgt_file_name, field_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), field_size == grid_size, "Error happens when reading the remapping weights file \"%s\": the array size of the variable \"%s\" is not the size of the corresponding grid. Please verify.", wgt_file_name, field_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(data_type, DATA_TYPE_INT), "Software error in H2D_remapping_wgt_file_info::get_grid_field_checksum_value");
    checksum = calculate_checksum_of_array(data_buffer, field_size, get_data_type_size(data_type), NULL, NULL);

    delete [] data_buffer;

    return checksum;
}


bool H2D_remapping_wgt_file_info::match_H2D_remapping_wgt(Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid)
{
	bool check_result;
	double *original_grid_center_lon, *original_grid_center_lat;

	
    for (int i = 0; i < matched_grid_pair.size(); i ++)
        if (matched_grid_pair[i].first == src_original_grid && matched_grid_pair[i].second == dst_original_grid)
            return true;

    read_grid_size(dst_original_grid->get_comp_id(), "n_a", src_grid_size);
    if (src_grid_size != src_original_grid->get_H2D_sub_CoR_grid()->get_grid_size()) {
        clean();
        return false;
    }
    read_grid_size(dst_original_grid->get_comp_id(), "n_b", dst_grid_size);
    if (dst_grid_size != dst_original_grid->get_H2D_sub_CoR_grid()->get_grid_size()) {
        clean();
        return false;
    }
    get_checksum_mask(dst_original_grid->get_comp_id(), "mask_a", src_grid_size, checksum_src_mask);
    if (src_original_grid->get_checksum_H2D_mask() != this->checksum_src_mask) {
        clean();
        return false;
    }
    get_checksum_mask(dst_original_grid->get_comp_id(), "mask_b", dst_grid_size, checksum_dst_mask);
    if (dst_original_grid->get_checksum_H2D_mask() != this->checksum_dst_mask) {
        clean();
        return false;
    }
    read_weight_grid_data(dst_original_grid->get_comp_id(), "xc_a", DATA_TYPE_DOUBLE, (void**)(&src_center_lon), src_grid_size, true);
	original_grid_center_lon = src_original_grid->get_center_lon_values();
    check_result = are_two_coord_arrays_same(original_grid_center_lon, this->src_center_lon, src_original_grid->get_H2D_sub_CoR_grid()->get_grid_size(), this->src_grid_size);
	delete [] original_grid_center_lon;
    if (!check_result) {
        clean();
        return false;
    }
    read_weight_grid_data(dst_original_grid->get_comp_id(), "xc_b", DATA_TYPE_DOUBLE, (void**)(&dst_center_lon), dst_grid_size, true);
	original_grid_center_lon = dst_original_grid->get_center_lon_values();
    check_result = are_two_coord_arrays_same(original_grid_center_lon, this->dst_center_lon, dst_original_grid->get_H2D_sub_CoR_grid()->get_grid_size(), this->dst_grid_size);
	delete [] original_grid_center_lon;
    if (!check_result) {
        clean();
        return false;
    }
    read_weight_grid_data(dst_original_grid->get_comp_id(), "yc_a", DATA_TYPE_DOUBLE, (void**)(&src_center_lat), src_grid_size, true);
	original_grid_center_lat = src_original_grid->get_center_lat_values();
    check_result = are_two_coord_arrays_same(original_grid_center_lat, this->src_center_lat, src_original_grid->get_H2D_sub_CoR_grid()->get_grid_size(), this->src_grid_size);
	delete [] original_grid_center_lat;
    if (!check_result) {
        clean();
        return false;
    }
    read_weight_grid_data(dst_original_grid->get_comp_id(), "yc_b", DATA_TYPE_DOUBLE, (void**)(&dst_center_lat), dst_grid_size, true);
	original_grid_center_lat = dst_original_grid->get_center_lat_values();
    check_result = are_two_coord_arrays_same(original_grid_center_lat, this->dst_center_lat, dst_original_grid->get_H2D_sub_CoR_grid()->get_grid_size(), this->dst_grid_size);
	delete [] original_grid_center_lat;
    if (!check_result) {
        clean();
        return false;
    }
    read_weight_grid_data(dst_original_grid->get_comp_id(), "area_a", DATA_TYPE_DOUBLE, (void**)(&src_area), src_grid_size, false);
    read_weight_grid_data(dst_original_grid->get_comp_id(), "area_b", DATA_TYPE_DOUBLE, (void**)(&dst_area), dst_grid_size, false);
    read_remapping_weights(dst_original_grid->get_comp_id());

    matched_grid_pair.push_back(std::make_pair(src_original_grid, dst_original_grid));
    return true;
}


void H2D_remapping_wgt_file_info::get_checksum_mask(int comp_id, const char *mask_label, int grid_size, long &checksum_mask)
{
    int local_proc_id_in_file_read_comm, num_procs_in_file_read_comm;
    int wgts_status_tag = checksum_mask != -1? 1 : 0;
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in H2D_remapping_wgt_file_info::read_remapping_weights");
    MPI_Comm file_read_comm;


    MPI_Comm_split(comp_node->get_comm_group(), wgts_status_tag, 0, &file_read_comm);
    if (wgts_status_tag == 1) {
        MPI_Comm_free(&file_read_comm);
        return;
    }
    MPI_Comm_rank(file_read_comm, &local_proc_id_in_file_read_comm);
    MPI_Comm_size(file_read_comm, &num_procs_in_file_read_comm);
    
    IO_netcdf *netcdf_file_object = new IO_netcdf("remapping weights file for H2D interpolation", wgt_file_name, "r", false);    
    checksum_mask = get_grid_field_checksum_value(mask_label, netcdf_file_object, grid_size, file_read_comm, local_proc_id_in_file_read_comm == 0);    
    delete netcdf_file_object;
    
    MPI_Comm_free(&file_read_comm);
}


void H2D_remapping_wgt_file_info::read_grid_size(int comp_id, const char *label, int &grid_size)
{
    int field_size, i, local_proc_id_in_file_read_comm, num_procs_in_file_read_comm;
    char data_type[NAME_STR_SIZE];
    int wgts_status_tag = grid_size != -1? 1 : 0;
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in H2D_remapping_wgt_file_info::read_remapping_weights");
    MPI_Comm file_read_comm;


    MPI_Comm_split(comp_node->get_comm_group(), wgts_status_tag, 0, &file_read_comm);
    if (wgts_status_tag == 1) {
        MPI_Comm_free(&file_read_comm);
        return;
    }
    MPI_Comm_rank(file_read_comm, &local_proc_id_in_file_read_comm);
    MPI_Comm_size(file_read_comm, &num_procs_in_file_read_comm);

    IO_netcdf *netcdf_file_object = new IO_netcdf("remapping weights file for H2D interpolation", wgt_file_name, "r", false);
    grid_size = netcdf_file_object->get_dimension_size(label, file_read_comm, local_proc_id_in_file_read_comm == 0);
    delete netcdf_file_object;

    MPI_Comm_free(&file_read_comm);
}
            

void H2D_remapping_wgt_file_info::read_weight_grid_data(int comp_id, const char *label, const char *required_data_type, void **buffer_ptr, int buffer_size, bool necessary)
{
    int field_size, i, local_proc_id_in_file_read_comm, num_procs_in_file_read_comm;
    char data_type[NAME_STR_SIZE];
    int wgts_status_tag = (*buffer_ptr) != NULL? 1 : 0;
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in H2D_remapping_wgt_file_info::read_remapping_weights");
    MPI_Comm file_read_comm;


    MPI_Comm_split(comp_node->get_comm_group(), wgts_status_tag, 0, &file_read_comm);
    if (wgts_status_tag == 1) {
        MPI_Comm_free(&file_read_comm);
        return;
    }
    MPI_Comm_rank(file_read_comm, &local_proc_id_in_file_read_comm);
    MPI_Comm_size(file_read_comm, &num_procs_in_file_read_comm);

    IO_netcdf *netcdf_file_object = new IO_netcdf("remapping weights file for H2D interpolation", wgt_file_name, "r", false);
    netcdf_file_object->read_file_field(label, buffer_ptr, &field_size, data_type, file_read_comm, local_proc_id_in_file_read_comm == 0);
	if (necessary)
	    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), field_size == buffer_size && words_are_the_same(data_type, required_data_type), "Error happens when reading the remapping weights file \"%s\": fail to read the variable \"%s\" because of wrong array size or wrong data type (should be %s). Please verify.", wgt_file_name, label, required_data_type);
    delete netcdf_file_object;

    MPI_Comm_free(&file_read_comm);
}

    
void H2D_remapping_wgt_file_info::read_remapping_weights(int comp_id)
{
    int field_size, i, local_proc_id_in_file_read_comm, num_procs_in_file_read_comm;
    char data_type[NAME_STR_SIZE];
    int *temp_wgts_src_indexes, *temp_wgts_dst_indexes;
    int wgts_status_tag = wgts_src_indexes != NULL? 1 : 0;
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in H2D_remapping_wgt_file_info::read_remapping_weights");
    MPI_Comm file_read_comm;


    MPI_Comm_split(comp_node->get_comm_group(), wgts_status_tag, 0, &file_read_comm);

    if (wgts_status_tag == 1) {
        MPI_Comm_free(&file_read_comm);
        return;
    }

    MPI_Comm_rank(file_read_comm, &local_proc_id_in_file_read_comm);
    MPI_Comm_size(file_read_comm, &num_procs_in_file_read_comm);

    if (num_procs_in_file_read_comm < comp_node->get_num_procs()) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Partially load remapping weight file %s", wgt_file_name);
    }
    else EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Load remapping weight file %s", wgt_file_name);

    IO_netcdf *netcdf_file_object = new IO_netcdf("remapping weights file for H2D interpolation", wgt_file_name, "r", false);
    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), src_grid_size > 0, "Error happens when reading the remapping weights file \"%s\": fail to read the size of the source grid (dimension \"n_a\" in the file). Please verify.", wgt_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), dst_grid_size > 0, "Error happens when reading the remapping weights file \"%s\": fail to read the size of the target grid (dimension \"n_a\" in the file). Please verify.", wgt_file_name);
    num_wgts = netcdf_file_object->get_dimension_size("n_s", file_read_comm, local_proc_id_in_file_read_comm == 0);    
    netcdf_file_object->read_file_field("col", (void**)(&temp_wgts_src_indexes), &field_size, data_type, file_read_comm, local_proc_id_in_file_read_comm == 0);
    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), field_size == num_wgts && words_are_the_same(data_type, DATA_TYPE_INT), "Error happens when reading the remapping weights file \"%s\": fail to read the variable \"col\" because of wrong array size (should be dimension of \"n_s\") or wrong data type (should be integer). Please verify.", wgt_file_name);
    netcdf_file_object->read_file_field("row", (void**)(&temp_wgts_dst_indexes), &field_size, data_type, file_read_comm, local_proc_id_in_file_read_comm == 0);
    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), field_size == num_wgts && words_are_the_same(data_type, DATA_TYPE_INT), "Error happens when reading the remapping weights file \"%s\": fail to read the variable \"row\" because of wrong array size (should be dimension of \"n_s\") or wrong data type (should be integer). Please verify.", wgt_file_name);
    netcdf_file_object->read_file_field("S", (void**)(&wgts_values), &field_size, data_type, file_read_comm, local_proc_id_in_file_read_comm == 0);
    EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), field_size == num_wgts && words_are_the_same(data_type, DATA_TYPE_DOUBLE), "Error happens when reading the remapping weights file \"%s\": fail to read the variable \"S\" because of wrong array size (should be dimension of \"n_s\") or wrong data type (should be double). Please verify.", wgt_file_name);

    wgts_src_indexes = new long [num_wgts];
    wgts_dst_indexes = new long [num_wgts];
    for (int i = 0; i < num_wgts; i ++) {
        wgts_src_indexes[i] = temp_wgts_src_indexes[i] - 1;
        wgts_dst_indexes[i] = temp_wgts_dst_indexes[i] - 1;
        EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), wgts_src_indexes[i] >= 0 && wgts_src_indexes[i] < src_grid_size, "Error happens when reading the remapping weights file \"%s\": some values in the variable \"col\" are out of the bound of source grid size. Please verify.", wgt_file_name);
        EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), wgts_dst_indexes[i] >= 0 && wgts_dst_indexes[i] < dst_grid_size, "Error happens when reading the remapping weights file \"%s\": some values in the variable \"col\" are out of the bound of target grid size. Please verify.", wgt_file_name);
    }

    delete netcdf_file_object;
    delete [] temp_wgts_dst_indexes;
    delete [] temp_wgts_src_indexes;
    MPI_Comm_free(&file_read_comm);
}


void H2D_remapping_wgt_file_info::clean()
{
    if (wgts_src_indexes != NULL) {
        delete [] wgts_src_indexes;
        delete [] wgts_dst_indexes;
        delete [] wgts_values;
        wgts_src_indexes = NULL;
        wgts_dst_indexes = NULL;
        wgts_values = NULL;
    }
    if (src_center_lon != NULL) {
        delete [] src_center_lon;
        src_center_lon = NULL;
    }
    if (src_center_lat != NULL) {
        delete [] src_center_lat;
        src_center_lat = NULL;
    }
    if (dst_center_lon != NULL) {
        delete [] dst_center_lon;
        dst_center_lon = NULL;
    }
    if (dst_center_lat != NULL) {
        delete [] dst_center_lat;
        dst_center_lat = NULL;        
    }
    if (src_area != NULL) {
        delete [] src_area;
        src_area = NULL;
    }
    if (dst_area != NULL) {
        delete [] dst_area;
        dst_area = NULL;
    }
}


H2D_remapping_wgt_file_info::~H2D_remapping_wgt_file_info()
{
    clean();
    delete [] wgt_file_name;
}


H2D_remapping_wgt_file_mgt::H2D_remapping_wgt_file_mgt(TiXmlElement *XML_element, const char *XML_file_name)
{
    char overall_XML_file_name[NAME_STR_SIZE];


    sprintf(overall_XML_file_name, "%s/all/overall_remapping_configuration.xml", comp_comm_group_mgt_mgr->get_config_root_dir());
    
    for (TiXmlNode *detailed_element_node = XML_element->FirstChild(); detailed_element_node != NULL; detailed_element_node = detailed_element_node->NextSibling()) {
        if (detailed_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *detailed_element = detailed_element_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), words_are_the_same(detailed_element->Value(), "file"), "When setting a remapping weights file in the remapping configuration in the XML file \"%s\", \"%s\" is not a legal attribute. Please verify the XML file arround the line number %d.", XML_file_name, detailed_element->Value(), detailed_element->Row());
        const char *short_file_name = get_XML_attribute(comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id(), 1000, detailed_element, "name", XML_file_name, line_number, "the name of a remapping weights file",  "remapping configuration", true);
        char full_file_name[NAME_STR_SIZE];
        if (words_are_the_same(overall_XML_file_name, XML_file_name))
            sprintf(full_file_name, "%s/all/grids_weights/%s", comp_comm_group_mgt_mgr->get_config_root_dir(), short_file_name);
        else sprintf(full_file_name, "%s/grids_weights/%s", comp_comm_group_mgt_mgr->get_root_comp_config_dir(), short_file_name);
        H2D_remapping_wgt_file_info *existing_wgt_file_info = all_H2D_remapping_wgt_files_info->search_wgt_file_info(full_file_name);
        if (existing_wgt_file_info == NULL) {
            H2D_remapping_wgt_files.push_back(new H2D_remapping_wgt_file_info(full_file_name));
            all_H2D_remapping_wgt_files_info->add_wgt_file_info(H2D_remapping_wgt_files[H2D_remapping_wgt_files.size()-1]);
        }
        else H2D_remapping_wgt_files.push_back(existing_wgt_file_info);
    }
}


H2D_remapping_wgt_file_mgt::H2D_remapping_wgt_file_mgt(const char *array, long *buffer_content_iter)
{
    int temp_int;
    read_data_from_array_buffer(&temp_int, sizeof(int), array, *buffer_content_iter, true);
    for (int i = 0; i < temp_int; i ++) {
        H2D_remapping_wgt_file_info *new_wgt_file_info = new H2D_remapping_wgt_file_info(array, buffer_content_iter);
        H2D_remapping_wgt_file_info *existing_wgt_file_info = all_H2D_remapping_wgt_files_info->search_wgt_file_info(new_wgt_file_info->get_wgt_file_name());
        if (existing_wgt_file_info != NULL) {
            delete new_wgt_file_info;
            H2D_remapping_wgt_files.push_back(existing_wgt_file_info);
        }
        else {
            H2D_remapping_wgt_files.push_back(new_wgt_file_info);
            all_H2D_remapping_wgt_files_info->add_wgt_file_info(new_wgt_file_info);
        }
    }
}


void H2D_remapping_wgt_file_mgt::write_remapping_wgt_files_info_into_array(char **array, long &buffer_max_size, long &buffer_content_size)
{
    for (int i = H2D_remapping_wgt_files.size()-1; i >= 0; i --)
        H2D_remapping_wgt_files[i]->write_remapping_wgt_file_info_into_array(array, buffer_max_size, buffer_content_size);
    int temp_int = H2D_remapping_wgt_files.size();    
    write_data_into_array_buffer(&temp_int, sizeof(int), array, buffer_max_size, buffer_content_size);    
}


void H2D_remapping_wgt_file_mgt::append_remapping_weights(H2D_remapping_wgt_file_mgt *another_mgr)
{
    for (int i = 0; i < another_mgr->H2D_remapping_wgt_files.size(); i ++)
        this->H2D_remapping_wgt_files.push_back(another_mgr->H2D_remapping_wgt_files[i]);
}


H2D_remapping_wgt_file_info *H2D_remapping_wgt_file_mgt::search_H2D_remapping_weight(Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid, const char *algorithm_name, int comp_id)
{
	char full_default_wgt_file_name[NAME_STR_SIZE*2];

	
    if (src_original_grid->get_original_CoR_grid() != src_original_grid->get_H2D_sub_CoR_grid() && 
        (src_original_grid->get_original_CoR_grid()->get_grid_mask_field() != NULL || dst_original_grid->get_original_CoR_grid()->get_grid_mask_field() != NULL))
        return NULL;

    for (int i = 0; i < H2D_remapping_wgt_files.size(); i ++)
        if (H2D_remapping_wgt_files[i]->match_H2D_remapping_wgt(src_original_grid, dst_original_grid))
            return H2D_remapping_wgt_files[i];

	if (src_original_grid->get_H2D_sub_CoR_grid() != NULL && algorithm_name != NULL) {	
		sprintf(full_default_wgt_file_name, "%s/DEFAULT_WGT_of___%s___FROM___%s___TO___%s___AT___%s.nc", comp_comm_group_mgt_mgr->get_internal_remapping_weights_dir(), algorithm_name, src_original_grid->get_H2D_sub_CoR_grid()->get_grid_name(), dst_original_grid->get_H2D_sub_CoR_grid()->get_grid_name(), comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_full_name());
		EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Default remapping weight file \"%s\" will be checked if avaiable", full_default_wgt_file_name);
		H2D_remapping_wgt_file_info *existing_wgt_file = all_H2D_remapping_wgt_files_info->search_wgt_file_info(full_default_wgt_file_name);
		if (existing_wgt_file != NULL)
			return existing_wgt_file;
		if (does_file_exist(full_default_wgt_file_name)) {
			EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Default remapping weight file \"%s\" exists", full_default_wgt_file_name);
			H2D_remapping_wgt_file_info *new_wgt_file = new H2D_remapping_wgt_file_info(full_default_wgt_file_name);
			if (new_wgt_file->match_H2D_remapping_wgt(src_original_grid, dst_original_grid)) {
				EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Default remapping weight file \"%s\" will be used", full_default_wgt_file_name);
				all_H2D_remapping_wgt_files_info->add_wgt_file_info(new_wgt_file);
				return new_wgt_file;
			}
			else delete new_wgt_file;
		}
	}
        
    return NULL;    
}


H2D_remapping_wgt_file_info *H2D_remapping_wgt_file_mgt::search_H2D_remapping_weight(const char *file_name)
{
	for (int i = 0; i < H2D_remapping_wgt_files.size(); i ++)
		if (words_are_the_same(H2D_remapping_wgt_files[i]->get_wgt_file_name(), file_name))
			return H2D_remapping_wgt_files[i];
		
	return NULL;
}


void H2D_remapping_wgt_file_mgt::print()
{
    printf("   H2D remapping wgts files\n");
    for(int i = 0; i < H2D_remapping_wgt_files.size(); i ++)
        printf("        %s\n", H2D_remapping_wgt_files[i]->get_wgt_file_name());
}


bool H2D_remapping_wgt_file_mgt::is_the_same_as_another(H2D_remapping_wgt_file_mgt *another_mgr)
{
    if (this->H2D_remapping_wgt_files.size() != another_mgr->H2D_remapping_wgt_files.size())
        return false;

    for (int i = 0; i < this->H2D_remapping_wgt_files.size(); i ++)
        if (this->H2D_remapping_wgt_files[i] != another_mgr->H2D_remapping_wgt_files[i])
            return false;

    return true;
}


H2D_remapping_wgt_file_mgt *H2D_remapping_wgt_file_mgt::clone()
{
    H2D_remapping_wgt_file_mgt *cloned_mgr = new H2D_remapping_wgt_file_mgt();
    for (int i = 0; i < this->H2D_remapping_wgt_files.size(); i ++)
        cloned_mgr->H2D_remapping_wgt_files.push_back(this->H2D_remapping_wgt_files[i]);

    return cloned_mgr;
}


void H2D_remapping_wgt_file_mgt::shrink(Original_grid_info *src_grid, Original_grid_info *dst_grid)
{
    H2D_remapping_wgt_file_info *remapping_file = search_H2D_remapping_weight(src_grid, dst_grid, NULL, -1);
    H2D_remapping_wgt_files.clear();
    if (remapping_file != NULL)
        H2D_remapping_wgt_files.push_back(remapping_file);
}


H2D_remapping_wgt_file_container::~H2D_remapping_wgt_file_container()
{
    for (int i = 0; i < H2D_remapping_wgt_files.size(); i ++)
        delete H2D_remapping_wgt_files[i];
}


H2D_remapping_wgt_file_info *H2D_remapping_wgt_file_container::search_wgt_file_info(const char *wgt_file_name)
{
    for (int i = 0; i < H2D_remapping_wgt_files.size(); i ++)
        if (words_are_the_same(H2D_remapping_wgt_files[i]->get_wgt_file_name(), wgt_file_name))
            return H2D_remapping_wgt_files[i];

    return NULL;
}


void H2D_remapping_wgt_file_container::add_wgt_file_info(H2D_remapping_wgt_file_info *wgt_file_info)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, search_wgt_file_info(wgt_file_info->get_wgt_file_name()) == NULL, "Software error in H2D_remapping_wgt_file_container::add_wgt_file_info");
    H2D_remapping_wgt_files.push_back(wgt_file_info);
}


Remapping_algorithm_specification::Remapping_algorithm_specification(const Remapping_algorithm_specification *src_specification)
{
    this->type_id = src_specification->type_id;
    this->algorithm_name = strdup(src_specification->algorithm_name);
    for (int i = 0; i < src_specification->parameters_name.size(); i ++) {
        this->parameters_name.push_back(new char [NAME_STR_SIZE]);
        strcpy(this->parameters_name[this->parameters_name.size()-1], src_specification->parameters_name[i]);
        this->parameters_value.push_back(new char [NAME_STR_SIZE]);
        strcpy(this->parameters_value[this->parameters_value.size()-1], src_specification->parameters_value[i]);
    }
}


Remapping_algorithm_specification::Remapping_algorithm_specification(const char *algorithm_name, int algorithm_type)
{
    this->algorithm_name = strdup(algorithm_name);
    this->type_id = algorithm_type;
}


Remapping_algorithm_specification::Remapping_algorithm_specification(int comp_id, TiXmlElement *XML_element, const char *XML_file_name, int algorithm_type)
{
    int line_number;

    
    this->type_id = algorithm_type;
    this->comp_id = comp_id;
    const char *algorithm_name = get_XML_attribute(comp_id, -1, XML_element, "name", XML_file_name, line_number, "the name of a remapping algorithm", "remapping configuration", true);
    if (algorithm_type == REMAP_ALGORITHM_TYPE_H2D)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, remap_operator_manager->get_remap_operator_num_dim(algorithm_name) == 2, "\"%s\" is not a legal remapping operator or not a 2D remapping operator. Please verify the XML file \"%s\" around the line number %d", algorithm_name, XML_file_name, line_number);
    else if (algorithm_type == REMAP_ALGORITHM_TYPE_V1D || algorithm_type == REMAP_ALGORITHM_TYPE_T1D)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, remap_operator_manager->get_remap_operator_num_dim(algorithm_name) == 1, "\"%s\" is not a legal remapping operator or not a 1D remapping operator. Please verify the XML file \"%s\" around the line number %d", algorithm_name, XML_file_name, line_number);
    else EXECUTION_REPORT(REPORT_ERROR, -1, "Software error in Remapping_algorithm_specification::Remapping_algorithm_specification");

    this->algorithm_name = strdup(algorithm_name);
    for (TiXmlNode *detailed_element_node = XML_element->FirstChild(); detailed_element_node != NULL; detailed_element_node = detailed_element_node->NextSibling()) {
        if (detailed_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *detailed_element = detailed_element_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(detailed_element->Value(), "parameter"), "When setting the remapping configuration in the XML file \"%s\", \"%s\" is not a legal attribute. Please verify the XML file arround the line number %d.", XML_file_name, detailed_element->Value(), detailed_element->Row());
        const char *parameter_name = get_XML_attribute(comp_id, -1, detailed_element, "name", XML_file_name, line_number, "the name of a parameter of the corresponding remapping algorithm",  "remapping configuration", true);
        const char *parameter_value = get_XML_attribute(comp_id, 1000, detailed_element, "value", XML_file_name, line_number, "the value of a parameter of the corresponding remapping algorithm",  "remapping configuration", true);
        char error_string[NAME_STR_SIZE];
        int parameter_check_result = remap_operator_manager->check_operator_parameter(algorithm_name, parameter_name, parameter_value, error_string);
        if (parameter_check_result == 0)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "The remapping algorithm \"%s\" does not have a parameter named \"%s\". Please verify the XML file \"%s\" around the line number %d", algorithm_name, parameter_name, XML_file_name, line_number);
        else if (parameter_check_result == 1)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "The value of the parameter \"%s\" of the remapping algorithm \"%s\" is wrong. %s. Please verify the XML file \"%s\" around the line number %d", parameter_name, algorithm_name, error_string, XML_file_name, line_number);
        parameters_name.push_back(strdup(parameter_name));
        parameters_value.push_back(strdup(parameter_value));
    }
}


Remapping_algorithm_specification::Remapping_algorithm_specification(const char *array, long *buffer_content_iter)
{
    long str_size;

    
    read_data_from_array_buffer(&type_id, sizeof(int), array, *buffer_content_iter, true);
    algorithm_name = load_string(NULL, str_size, NAME_STR_SIZE, array, *buffer_content_iter, NULL);
    int temp_int;
    read_data_from_array_buffer(&temp_int, sizeof(int), array, *buffer_content_iter, true);
    for (int i = 0; i < temp_int; i ++) {
        this->parameters_name.push_back(new char [NAME_STR_SIZE]);
        read_data_from_array_buffer(this->parameters_name[this->parameters_name.size()-1], NAME_STR_SIZE, array, *buffer_content_iter, true);
        this->parameters_value.push_back(new char [NAME_STR_SIZE]);
        read_data_from_array_buffer(this->parameters_value[this->parameters_value.size()-1], NAME_STR_SIZE, array, *buffer_content_iter, true);        
    }
}


Remapping_algorithm_specification::~Remapping_algorithm_specification()
{
    for (int i = 0; i < parameters_name.size(); i ++) {
        delete [] parameters_name[i];
        delete [] parameters_value[i];
    }
    delete [] algorithm_name;
}


void Remapping_algorithm_specification::print()
{
    if (type_id == REMAP_ALGORITHM_TYPE_H2D)
        printf("   H2D remapping algorithm \"%s\" ", algorithm_name);
    if (type_id == REMAP_ALGORITHM_TYPE_V1D)
        printf("   V1D remapping algorithm \"%s\" ", algorithm_name);
    if (type_id == REMAP_ALGORITHM_TYPE_T1D)
        printf("   T1D remapping algorithm \"%s\" ", algorithm_name);
    for (int i = 0; i < parameters_name.size(); i ++)
        printf(": \"%s\"(\"%s\") ", parameters_name[i], parameters_value[i]);
    printf("\n");
}


void Remapping_algorithm_specification::clean()
{
    for (int i = 0; i < parameters_name.size(); i ++) {
        delete [] parameters_value[i];
        delete [] parameters_name[i];
    }
    parameters_name.clear();
    parameters_value.clear();
    delete [] algorithm_name;
    algorithm_name = NULL;
    type_id = -1;
}


void Remapping_algorithm_specification::write_remapping_algorithm_specification_into_array(char **array, long &buffer_max_size, long &buffer_content_size)
{
    for (int i = parameters_name.size()-1; i >= 0; i --) {
        write_string_into_array_buffer(parameters_value[i], NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
        write_string_into_array_buffer(parameters_name[i], NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    }
    int temp_int = parameters_name.size();
    write_data_into_array_buffer(&temp_int, sizeof(int), array, buffer_max_size, buffer_content_size);
    dump_string(algorithm_name, -1, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&type_id, sizeof(int), array, buffer_max_size, buffer_content_size);
}


void Remapping_algorithm_specification::get_parameter(int i, char *parameter_name, char *parameter_value)
{
    EXECUTION_REPORT(REPORT_ERROR, comp_id, i >= 0 && i < parameters_name.size(), "Software error in Remapping_algorithm_specification::get_parameter");
    strcpy(parameter_name, parameters_name[i]);
    strcpy(parameter_value, parameters_value[i]);
}


Remapping_algorithm_specification *Remapping_algorithm_specification::clone()
{
    Remapping_algorithm_specification *cloned_specification = new Remapping_algorithm_specification;
    cloned_specification->comp_id = this->comp_id;
    cloned_specification->type_id = this->type_id;
    if (this->algorithm_name != NULL)
        cloned_specification->algorithm_name = strdup(this->algorithm_name);
    for (int i = 0; i < this->parameters_name.size(); i ++) {
        cloned_specification->parameters_name.push_back(strdup(this->parameters_name[i]));
        cloned_specification->parameters_value.push_back(strdup(this->parameters_value[i]));
    }

    return cloned_specification;
}


bool Remapping_algorithm_specification::is_the_same_as_another(Remapping_algorithm_specification *another)
{
    if (another->type_id != this->type_id)
        return false;
    if (!words_are_the_same(another->algorithm_name,this->algorithm_name))
        return false;
    if (another->parameters_name.size() != this->parameters_name.size())
        return false;
    for (int i = 0; i < this->parameters_name.size(); i ++) {
        if (!words_are_the_same(another->parameters_name[i],this->parameters_name[i]) || !words_are_the_same(another->parameters_value[i],this->parameters_value[i]))
            return false;
    }

    return true;    
}


Remapping_setting::Remapping_setting()
{
    H2D_remapping_algorithm = NULL;
    V1D_remapping_algorithm = NULL;
    T1D_remapping_algorithm = NULL;
    H2D_remapping_wgt_file_mgr = NULL;
}


Remapping_setting::Remapping_setting(const char *H2D_remapping_algorithm_name, const char *field_type)
{
    this->comp_id = comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id();
    H2D_remapping_algorithm = new Remapping_algorithm_specification(H2D_remapping_algorithm_name, REMAP_ALGORITHM_TYPE_H2D);
    V1D_remapping_algorithm = new Remapping_algorithm_specification(REMAP_OPERATOR_NAME_LINEAR, REMAP_ALGORITHM_TYPE_V1D);
    T1D_remapping_algorithm = new Remapping_algorithm_specification(REMAP_OPERATOR_NAME_LINEAR, REMAP_ALGORITHM_TYPE_T1D);
    H2D_remapping_wgt_file_mgr = NULL;
    field_specification_manner = 1;
    fields_specification.push_back(strdup(field_type));
}


Remapping_setting::Remapping_setting(int comp_id, TiXmlElement *XML_element, const char *XML_file_name)
{
    int i, line_number;
    int num_remapping_algorithm_section = 0, num_fields_section = 0, num_remapping_weights_section = 0;


    this->comp_id = comp_id;
    this->XML_start_line_number = XML_element->Row();
    H2D_remapping_algorithm = NULL;
    V1D_remapping_algorithm = NULL;
    T1D_remapping_algorithm = NULL;
    H2D_remapping_wgt_file_mgr = NULL;
    
    for (TiXmlNode *detailed_element_node = XML_element->FirstChild(); detailed_element_node != NULL; detailed_element_node = detailed_element_node->NextSibling()) {
        if (detailed_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *detailed_element = detailed_element_node->ToElement();
        if (words_are_the_same(detailed_element->Value(), "remapping_algorithms")) {        
            if (!is_XML_setting_on(comp_id, detailed_element, XML_file_name, "the status of the configuration of a section about remapping algorithms", "remapping configuration"))
                continue;
            EXECUTION_REPORT(REPORT_ERROR, comp_id, num_remapping_algorithm_section == 0, "When setting the remapping configuration in the XML file \"%s\", there are more than one section for specifying remapping algorithms. That is not allowed. Please verify the XML file arround the line number %d", XML_file_name, detailed_element->Row());        
            for (TiXmlNode *algorithm_element_node = detailed_element_node->FirstChild(); algorithm_element_node != NULL; algorithm_element_node = algorithm_element_node->NextSibling()) {
                if (algorithm_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                    continue;
                TiXmlElement *algorithm_element = algorithm_element_node->ToElement();
                if (words_are_the_same(algorithm_element->Value(), "H2D_algorithm")) {
                    if (!is_XML_setting_on(comp_id, algorithm_element, XML_file_name, "the status of the configuration of a section about H2D_algorithm", "remapping configuration"))
                        continue;
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, H2D_remapping_algorithm == NULL, "When setting the remapping configuration in the XML file \"%s\", H2D_algorithm has been set more than once. Please verify the XML file arround the line number %d", XML_file_name, detailed_element->Row());
                    H2D_remapping_algorithm = new Remapping_algorithm_specification(comp_id, algorithm_element, XML_file_name, REMAP_ALGORITHM_TYPE_H2D);
                }
                else if (words_are_the_same(algorithm_element->Value(), "V1D_algorithm")) {
                    if (!is_XML_setting_on(comp_id, algorithm_element, XML_file_name, "the status of the configuration of a section about V1D_algorithm", "remapping configuration"))
                        continue;
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, V1D_remapping_algorithm == NULL, "When setting the remapping configuration in the XML file \"%s\", V1D_algorithm has been set more than once. Please verify the XML file arround the line number %d", XML_file_name, detailed_element->Row());
                    V1D_remapping_algorithm = new Remapping_algorithm_specification(comp_id, algorithm_element, XML_file_name, REMAP_ALGORITHM_TYPE_V1D);
                }
                else if (words_are_the_same(algorithm_element->Value(), "H2D_weights")) {
                    if (!is_XML_setting_on(comp_id, algorithm_element, XML_file_name, "the status of the configuration of a section about remapping weights files", "remapping configuration"))
                        continue;
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_remapping_weights_section == 0, "When setting the remapping configuration in the XML file \"%s\", there are more than one active section for specifying remapping weights files. That is not allowed. Please verify the XML file arround the line number %d", XML_file_name, detailed_element->Row());        
                    H2D_remapping_wgt_file_mgr = new H2D_remapping_wgt_file_mgt(algorithm_element, XML_file_name);
                    num_remapping_weights_section ++;
                }
                else EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "When setting the remapping configuration in the XML file \"%s\", \"%s\" is not a legal attribute of remapping algorithm. Please verify the XML file arround the line number %d", XML_file_name, algorithm_element->Value(), algorithm_element->Row());
            }
            num_remapping_algorithm_section ++;
        }
        else if (words_are_the_same(detailed_element->Value(),"fields")) {        
            if (!is_XML_setting_on(comp_id, detailed_element, XML_file_name, "the status of the configuration of a section about fields", "remapping configuration"))
                continue;
            EXECUTION_REPORT(REPORT_ERROR, comp_id, num_fields_section == 0, "When setting the remapping configuration in the XML file \"%s\", there are more than one section for specifying fields. That is not allowed. Please verify the XML file arround the line number %d", XML_file_name, detailed_element->Row());        
            const char *specification_type = get_XML_attribute(comp_id, -1, detailed_element, "specification", XML_file_name, line_number, "how to specify the fields corresponding to a remapping setting", "remapping configuration", true);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(specification_type, "type") || words_are_the_same(specification_type, "default") || words_are_the_same(specification_type, "name"), "In the XML file \"%s\", the manner for how to specify fields must be \"type\", \"default\" or \"name\". Please verify the XML file arround the line number %d.", XML_file_name, line_number);
            if (words_are_the_same(specification_type, "type")) {
                field_specification_manner = 1;
                for (TiXmlNode *type_element_node = detailed_element_node->FirstChild(); type_element_node != NULL; type_element_node = type_element_node->NextSibling()) {
                    if (type_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                        continue;
                    TiXmlElement *type_element = type_element_node->ToElement();
                    const char *field_type = get_XML_attribute(comp_id, -1, type_element, "value", XML_file_name, line_number, "the field type corresponding to a remapping setting", "remapping configuration", true);
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(field_type, "state") || words_are_the_same(field_type, "flux"), "In the XML file \"%s\" for remapping configuration, the field type \"%s\" is wrong. C-Coupler only supports field types \"state\" and \"flux\" at this time. Please verify the XML file arround the line number %d.", XML_file_name, field_type, line_number);
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, fields_specification.size() == 0, "In the XML file \"%s\" for remapping configuration, there are more than one field type specified while only one field type can be set for a remapping setting. Please verify the XML file arround the line number %d.", XML_file_name, type_element->Row());
                    fields_specification.push_back(strdup(field_type));
                }
                EXECUTION_REPORT(REPORT_ERROR, comp_id, fields_specification.size() > 0, "In the XML file \"%s\" for remapping configuration, no field type has been specified for a remapping setting. Please verify the XML file arround the line number %d.", XML_file_name, detailed_element->Row());
            }
            else if (words_are_the_same(specification_type, "default"))
                field_specification_manner = 0;
            else {
                field_specification_manner = 2;
                for (TiXmlNode *field_element_node = detailed_element_node->FirstChild(); field_element_node != NULL; field_element_node = field_element_node->NextSibling()) {
                    if (field_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                        continue;
                    TiXmlElement *field_element = field_element_node->ToElement();
                    const char *field_name = get_XML_attribute(comp_id, CCPL_NAME_STR_LEN, field_element, "value", XML_file_name, line_number, "the field name corresponding to a remapping setting", "remapping configuration", true);
                    for (i = 0; i < fields_specification.size(); i ++)
                        if (words_are_the_same(fields_specification[i], field_name))
                            break;
                    if (i == fields_specification.size())
                        fields_specification.push_back(strdup(field_name));
                }
                EXECUTION_REPORT(REPORT_ERROR, comp_id, fields_specification.size() > 0, "In the XML file \"%s\" for remapping configuration, no field name has been specified for a remapping setting. Please verify the XML file arround the line number %d.", XML_file_name, detailed_element->Row());
            }
            num_fields_section ++;
        }
        else EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "When setting the remapping configuration in the XML file \"%s\", \"%s\" is not a legal attribute. Please verify the XML file arround the line number %d.", XML_file_name, detailed_element->Value(), detailed_element->Row());
    }    
    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_remapping_algorithm_section == 1 && (H2D_remapping_algorithm != NULL || V1D_remapping_algorithm != NULL || H2D_remapping_wgt_file_mgr != NULL), "For the XML file \"%s\" that is for remapping configuration, no remapping algorithms or remapping weights files is specified for the remapping setting starting from the line number %d. Please verify.", XML_file_name, XML_element->Row());
    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_fields_section == 1, "For the XML file \"%s\" that is for remapping configuration, no fields is specified for the remapping setting starting from the line number %d. Please verify.", XML_file_name, XML_element->Row());
}


Remapping_setting::~Remapping_setting()
{
    for (int i = 0; i < fields_specification.size(); i ++)
        delete [] fields_specification[i];

    if (H2D_remapping_algorithm != NULL)
        delete H2D_remapping_algorithm;
    if (V1D_remapping_algorithm != NULL)
        delete V1D_remapping_algorithm;
    if (T1D_remapping_algorithm != NULL)
        delete T1D_remapping_algorithm;
    if (H2D_remapping_wgt_file_mgr != NULL)
        delete H2D_remapping_wgt_file_mgr;
}


void Remapping_setting::reset_remapping_setting()
{
    for (int i = 0; i < fields_specification.size(); i ++)
        delete [] fields_specification[i];
    fields_specification.clear();
    field_specification_manner = -1;
    H2D_remapping_algorithm = NULL;
    V1D_remapping_algorithm = NULL;
    T1D_remapping_algorithm = NULL;
}


void Remapping_setting::detect_conflict(Remapping_setting *another_setting, const char *XML_file_name)
{
    if (this->field_specification_manner == 0 && another_setting->field_specification_manner == 0) 
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "In the XML file \%s\" that is for remapping configuration, there is conflict between the remapping settings starting from line %d and %d respectively: both settings specify \"default\". Please verify. ", XML_file_name, another_setting->XML_start_line_number, this->XML_start_line_number);

    if (this->field_specification_manner == 1 && another_setting->field_specification_manner == 1)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, !words_are_the_same(this->fields_specification[0], another_setting->fields_specification[0]), "In the XML file \%s\" that is for remapping configuration, there is conflict (the same type of fields: \"%s\") between the remapping settings starting from line %d and %d respectively. Please verify. ", XML_file_name, this->fields_specification[0], another_setting->XML_start_line_number, this->XML_start_line_number);

    if (this->field_specification_manner == 2 && another_setting->field_specification_manner == 2) {
        for (int i = 0; i < this->fields_specification.size(); i ++)
            for (int j = 0; j < another_setting->fields_specification.size(); j ++)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(this->fields_specification[i], another_setting->fields_specification[j]), "In the XML file \%s\" that is for remapping configuration, there is conflict (the same field name: \"%s\") between the remapping settings starting from line %d and %d respectively. Please verify. ", XML_file_name, this->fields_specification[i], another_setting->XML_start_line_number, this->XML_start_line_number); 
    }
}


void Remapping_setting::append_H2D_remapping_weights(Remapping_setting *comp_remapping_setting)
{
    if (comp_remapping_setting->H2D_remapping_wgt_file_mgr == NULL)
        return;

    if (this->H2D_remapping_wgt_file_mgr == NULL)
        this->H2D_remapping_wgt_file_mgr = new H2D_remapping_wgt_file_mgt();

    this->H2D_remapping_wgt_file_mgr->append_remapping_weights(comp_remapping_setting->H2D_remapping_wgt_file_mgr);
}


void Remapping_setting::get_field_remapping_setting(Remapping_setting &field_remapping_configuration, const char *field_name)
{
    bool transfer_remapping_algorithms = false;

        
    if (!(field_remapping_configuration.H2D_remapping_algorithm == NULL && this->H2D_remapping_algorithm != NULL ||
          field_remapping_configuration.H2D_remapping_algorithm == NULL && this->H2D_remapping_wgt_file_mgr != NULL ||
          field_remapping_configuration.V1D_remapping_algorithm == NULL && this->V1D_remapping_algorithm != NULL ||
          field_remapping_configuration.T1D_remapping_algorithm == NULL && this->T1D_remapping_algorithm != NULL))
        return;

    if (field_specification_manner == 0)
        transfer_remapping_algorithms = true;
    else if (field_specification_manner == 1) {
        if (fields_info->search_field_info(field_name) == NULL)
            transfer_remapping_algorithms = true;
        else if (words_are_the_same(fields_info->search_field_info(field_name)->field_type, fields_specification[0]))
            transfer_remapping_algorithms = true;
    }
    else {
        for (int i = 0; i < fields_specification.size(); i ++)
            if (words_are_the_same(fields_specification[i], field_name)) {
                transfer_remapping_algorithms = true;
                break;
            }
    }

    if (transfer_remapping_algorithms) {
        if (field_remapping_configuration.H2D_remapping_algorithm == NULL) {
            if (this->H2D_remapping_algorithm != NULL)
                field_remapping_configuration.H2D_remapping_algorithm = new Remapping_algorithm_specification(this->H2D_remapping_algorithm);
            if (this->H2D_remapping_wgt_file_mgr != NULL)
                field_remapping_configuration.append_H2D_remapping_weights(this);
        }
        if (field_remapping_configuration.V1D_remapping_algorithm == NULL && this->V1D_remapping_algorithm != NULL)
            field_remapping_configuration.V1D_remapping_algorithm = new Remapping_algorithm_specification(this->V1D_remapping_algorithm);
        if (field_remapping_configuration.T1D_remapping_algorithm == NULL && this->T1D_remapping_algorithm != NULL)
            field_remapping_configuration.T1D_remapping_algorithm = new Remapping_algorithm_specification(this->T1D_remapping_algorithm);
    }
}


void Remapping_setting::print()
{
    printf("\n\nprint a remapping setting for \"%s\":\n", comp_comm_group_mgt_mgr->get_global_node_root()->get_full_name());
    if (H2D_remapping_algorithm != NULL)
        H2D_remapping_algorithm->print();
    if (V1D_remapping_algorithm != NULL)
        V1D_remapping_algorithm->print();
    if (T1D_remapping_algorithm != NULL)
        T1D_remapping_algorithm->print();
    if (H2D_remapping_wgt_file_mgr != NULL)
        H2D_remapping_wgt_file_mgr->print();
    printf("\n\n");
}


void Remapping_setting::write_remapping_setting_into_array(char **array, long &buffer_max_size, long &buffer_content_size)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, H2D_remapping_algorithm != NULL && V1D_remapping_algorithm != NULL && T1D_remapping_algorithm != NULL, "software error in Remapping_setting::write_remapping_setting_into_array");
    T1D_remapping_algorithm->write_remapping_algorithm_specification_into_array(array, buffer_max_size, buffer_content_size);
    V1D_remapping_algorithm->write_remapping_algorithm_specification_into_array(array, buffer_max_size, buffer_content_size);
    H2D_remapping_algorithm->write_remapping_algorithm_specification_into_array(array, buffer_max_size, buffer_content_size);
    if (H2D_remapping_wgt_file_mgr == NULL)
        H2D_remapping_wgt_file_mgr = new H2D_remapping_wgt_file_mgt();
    H2D_remapping_wgt_file_mgr->write_remapping_wgt_files_info_into_array(array, buffer_max_size, buffer_content_size);
}


void Remapping_setting::read_remapping_setting_from_array(const char *array, long &buffer_content_iter)
{    
    H2D_remapping_wgt_file_mgr = new H2D_remapping_wgt_file_mgt(array, &buffer_content_iter);
    H2D_remapping_algorithm = new Remapping_algorithm_specification(array, &buffer_content_iter);
    V1D_remapping_algorithm = new Remapping_algorithm_specification(array, &buffer_content_iter);
    T1D_remapping_algorithm = new Remapping_algorithm_specification(array, &buffer_content_iter);
    EXECUTION_REPORT(REPORT_ERROR, -1, buffer_content_iter == 0, "Software error in Remapping_setting::read_remapping_setting_from_array");
}


Remapping_setting *Remapping_setting::clone()
{
    Remapping_setting *cloned_setting = new Remapping_setting;
    cloned_setting->comp_id = this->comp_id;
    cloned_setting->XML_start_line_number = this->XML_start_line_number;
    cloned_setting->field_specification_manner = this->field_specification_manner;
    if (this->H2D_remapping_algorithm != NULL)
        cloned_setting->H2D_remapping_algorithm = this->H2D_remapping_algorithm->clone();
    if (this->V1D_remapping_algorithm != NULL)
        cloned_setting->V1D_remapping_algorithm = this->V1D_remapping_algorithm->clone();
    if (this->T1D_remapping_algorithm != NULL)
        cloned_setting->T1D_remapping_algorithm = this->T1D_remapping_algorithm->clone();
    if (this->H2D_remapping_wgt_file_mgr != NULL)
        cloned_setting->H2D_remapping_wgt_file_mgr = this->H2D_remapping_wgt_file_mgr->clone();
    for(int i = 0; i < this->fields_specification.size(); i ++)
        cloned_setting->fields_specification.push_back(strdup(this->fields_specification[i]));

    return cloned_setting;
}


bool Remapping_setting::is_the_same_as_another(Remapping_setting *another)
{
    if (another->field_specification_manner != this->field_specification_manner)
        return false;
    if (another->H2D_remapping_algorithm == NULL && this->H2D_remapping_algorithm != NULL || another->H2D_remapping_algorithm != NULL && this->H2D_remapping_algorithm == NULL)
        return false;
    if (another->V1D_remapping_algorithm == NULL && this->V1D_remapping_algorithm != NULL || another->V1D_remapping_algorithm != NULL && this->V1D_remapping_algorithm == NULL)
        return false;
    if (another->T1D_remapping_algorithm == NULL && this->T1D_remapping_algorithm != NULL || another->T1D_remapping_algorithm != NULL && this->T1D_remapping_algorithm == NULL)
        return false;

    EXECUTION_REPORT(REPORT_ERROR, -1, this->H2D_remapping_wgt_file_mgr != NULL && another->H2D_remapping_wgt_file_mgr != NULL, "Software error in Remapping_setting::is_the_same_as_another: empty remapping file managers");

    if (another->H2D_remapping_algorithm != NULL && !another->H2D_remapping_algorithm->is_the_same_as_another(this->H2D_remapping_algorithm))
        return false;
    if (another->V1D_remapping_algorithm != NULL && !another->V1D_remapping_algorithm->is_the_same_as_another(this->V1D_remapping_algorithm))
        return false;
    if (another->T1D_remapping_algorithm != NULL && !another->T1D_remapping_algorithm->is_the_same_as_another(this->T1D_remapping_algorithm))
        return false;
    if (!another->H2D_remapping_wgt_file_mgr->is_the_same_as_another(this->H2D_remapping_wgt_file_mgr))
        return false;

    if (another->fields_specification.size() != this->fields_specification.size())
        return false;
    for (int i = 0; i < this->fields_specification.size(); i ++)
        if (!words_are_the_same(another->fields_specification[i], this->fields_specification[i]))
            return false;

    return true;
}


H2D_remapping_wgt_file_info *Remapping_setting::search_H2D_remapping_weight(Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid, int remapping_host_comp_id)
{
    if (H2D_remapping_wgt_file_mgr == NULL)
        return NULL;

	if (H2D_remapping_algorithm != NULL)
	    return H2D_remapping_wgt_file_mgr->search_H2D_remapping_weight(src_original_grid, dst_original_grid, H2D_remapping_algorithm->get_algorithm_name(), remapping_host_comp_id);
	else return H2D_remapping_wgt_file_mgr->search_H2D_remapping_weight(src_original_grid, dst_original_grid, NULL, -1);
}


long Remapping_setting::calculate_checksum()
{
    char *temp_array = NULL;
    long buffer_max_size, buffer_content_size;


    write_remapping_setting_into_array(&temp_array, buffer_max_size, buffer_content_size);
    long checksum = calculate_checksum_of_array(temp_array, buffer_content_size, 1, NULL, NULL);
    delete [] temp_array;

    return checksum;
}


void Remapping_setting::shrink(Original_grid_info *src_grid, Original_grid_info *dst_grid)
{
    if (src_grid->get_V1D_sub_CoR_grid() == NULL)
        V1D_remapping_algorithm->clean();
    if (src_grid->get_T1D_sub_CoR_grid() == NULL)
        T1D_remapping_algorithm->clean();
    if (src_grid->get_H2D_sub_CoR_grid() == NULL)
        H2D_remapping_algorithm->clean();
    if (src_grid->get_H2D_sub_CoR_grid() == NULL && H2D_remapping_wgt_file_mgr != NULL)
        H2D_remapping_wgt_file_mgr->clean();
    if (H2D_remapping_wgt_file_mgr != NULL)
        H2D_remapping_wgt_file_mgr->shrink(src_grid, dst_grid);
    if (H2D_remapping_wgt_file_mgr != NULL && !H2D_remapping_wgt_file_mgr->is_empty()) {
        delete H2D_remapping_algorithm;
        H2D_remapping_algorithm = new Remapping_algorithm_specification(REMAP_OPERATOR_NAME_BILINEAR, REMAP_ALGORITHM_TYPE_H2D);
    }
}


Remapping_configuration::Remapping_configuration()
{
    comp_id = -1;
    remapping_settings.push_back(new Remapping_setting(REMAP_OPERATOR_NAME_BILINEAR, "state"));
    remapping_settings.push_back(new Remapping_setting(REMAP_OPERATOR_NAME_CONSERV_2D, "flux"));
}


Remapping_configuration::Remapping_configuration(int comp_id, const char *XML_file_name, TiXmlDocument *XML_file)
{
    this->comp_id = comp_id;
    for (TiXmlNode *XML_element_node = get_XML_first_child_of_unique_root(comp_id,XML_file_name,XML_file); XML_element_node != NULL; XML_element_node = XML_element_node->NextSibling()) {
        if (XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *XML_element = XML_element_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(XML_element->Value(), "remapping_setting"), "\"%s\" is not a legal attribute (the legal is \"remapping_setting\") for defining a remapping setting. Please verify the XML file arround the line number %d.", XML_element->Value(), XML_element->Row());
        if (!is_XML_setting_on(comp_id, XML_element, XML_file_name, "the status of a remapping setting", "remapping configuration"))
            continue;
        remapping_settings.push_back(new Remapping_setting(comp_id, XML_element, XML_file_name));
    }

    for (int i = 0; i < remapping_settings.size(); i ++)
        for (int j = i+1; j < remapping_settings.size(); j ++)
            remapping_settings[i]->detect_conflict(remapping_settings[j], XML_file_name);

    std::vector<Remapping_setting*> temp_remapping_settings;
    for (int i = 0; i < remapping_settings.size(); i ++)
        if (remapping_settings[i]->get_field_specification_manner() == 2)
            temp_remapping_settings.push_back(remapping_settings[i]);
    for (int i = 0; i < remapping_settings.size(); i ++)
        if (remapping_settings[i]->get_field_specification_manner() == 1)
            temp_remapping_settings.push_back(remapping_settings[i]);
    for (int i = 0; i < remapping_settings.size(); i ++)
        if (remapping_settings[i]->get_field_specification_manner() == 0)
            temp_remapping_settings.push_back(remapping_settings[i]);

    EXECUTION_REPORT(REPORT_ERROR, -1, temp_remapping_settings.size() == remapping_settings.size(), "Software error in Remapping_configuration::Remapping_configuration");

    remapping_settings.clear();
    for (int i = 0; i < temp_remapping_settings.size(); i ++)
        remapping_settings.push_back(temp_remapping_settings[i]);
}


Remapping_configuration::~Remapping_configuration()
{
    for (int i = 0; i < remapping_settings.size(); i ++)
        delete remapping_settings[i];
}


bool Remapping_configuration::get_field_remapping_setting(Remapping_setting &field_remapping_configuration, const char *field_name)
{
    for (int i = 0; i < remapping_settings.size(); i ++) {
        remapping_settings[i]->get_field_remapping_setting(field_remapping_configuration, field_name);
        if (field_remapping_configuration.get_H2D_remapping_algorithm() != NULL && field_remapping_configuration.get_V1D_remapping_algorithm() != NULL && field_remapping_configuration.get_T1D_remapping_algorithm() != NULL)
            return true;
    }
    return false;
}


Remapping_configuration_mgt::~Remapping_configuration_mgt()
{
    for (int i = 0; i < remapping_configurations.size(); i ++)
        delete remapping_configurations[i];
}


void Remapping_configuration_mgt::add_remapping_configuration(int comp_id)
{
    if (comp_id == comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id())
        remapping_configurations.push_back(new Remapping_configuration());

    char XML_file_name[NAME_STR_SIZE];
    Comp_comm_group_mgt_node *current_comp_node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false,"in Remapping_configuration_mgt::add_remapping_configuration");
    if (comp_id == comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id())
        sprintf(XML_file_name, "%s/all/overall_remapping_configuration.xml", comp_comm_group_mgt_mgr->get_config_root_dir());
    else sprintf(XML_file_name, "%s/remapping_configs/%s.remapping_configuration.xml", comp_comm_group_mgt_mgr->get_root_comp_config_dir(), comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false,"in Remapping_configuration_mgt::add_remapping_configuration")->get_full_name());
    TiXmlDocument *XML_file = open_XML_file_to_read(comp_id, XML_file_name, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in Remapping_configuration::Remapping_configuration"), false);
    if (XML_file == NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The remapping configuration file \"%s\" for the current component does not exist.", XML_file_name);
        return;
    }
    remapping_configurations.push_back(new Remapping_configuration(comp_id, XML_file_name, XML_file));
    delete XML_file;
}


Remapping_configuration *Remapping_configuration_mgt::search_remapping_configuration(int comp_id)
{
    for (int i = 0; i < remapping_configurations.size(); i ++)
        if (remapping_configurations[i]->get_comp_id() == comp_id)
            return remapping_configurations[i];

    return NULL;
}


void Remapping_configuration_mgt::get_field_remapping_setting(Remapping_setting &field_remapping_setting, int comp_id, const char *field_name)
{
    field_remapping_setting.reset_remapping_setting();
    Comp_comm_group_mgt_node *current_comp_node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in Remapping_configuration_mgt::get_current_remapping_setting");
    EXECUTION_REPORT(REPORT_ERROR, -1, current_comp_node != NULL, "Software error in Remapping_configuration_mgt::get_field_remapping_setting1");
    for (; current_comp_node != NULL; current_comp_node = current_comp_node->get_parent()) {
        Remapping_configuration *current_remapping_configuration = search_remapping_configuration(current_comp_node->get_comp_id());
        if (current_remapping_configuration != NULL)
            if (current_remapping_configuration->get_field_remapping_setting(field_remapping_setting, field_name)) {
//                field_remapping_setting.print();
                return;
            }
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, remapping_configurations[0]->get_field_remapping_setting(field_remapping_setting, field_name), "Software error in Remapping_configuration_mgt::get_field_remapping_setting for field %s", field_name);
//    field_remapping_setting.print();
}

