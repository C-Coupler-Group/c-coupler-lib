/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu and then
  *  modified by Dr. Cheng Zhang and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn or
  *  Dr. Cheng Zhang via zhangc-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "coupling_generator.h"



void clean_string_pair_vector(std::vector< std::pair<const char*, const char *> > &string_pair_vector)
{
    for (int i = 0; i < string_pair_vector.size(); i ++) {
        delete [] string_pair_vector[i].first;
        delete [] string_pair_vector[i].second;
    }
    string_pair_vector.clear();
}


MPI_Comm create_union_comm_common(MPI_Comm comp1, MPI_Comm comm2, int current_proc_id1, int current_proc_id2, std::vector<int> &procs_global_ids1, std::vector<int> &procs_global_ids2, int connection_id, int *proc_ranks_in_union_comm1, int *proc_ranks_in_union_comm2)
{
    MPI_Group common_group, comm_group1, comm_group2, union_group;
    MPI_Group intersection_group;
    MPI_Comm inter_comm, exclusive_comm, union_comm;
    int *proc_ranks1, *proc_ranks2;
    int intersection_size;
    int comp_num_procs1 = procs_global_ids1.size();
    int comp_num_procs2 = procs_global_ids2.size();
    

    proc_ranks1 = new int[comp_num_procs1];
    proc_ranks2 = new int[comp_num_procs2];

    for (int i = 0; i < comp_num_procs1; i ++)
        proc_ranks1[i] = procs_global_ids1[i];
    for (int i = 0; i < comp_num_procs2; i ++)
        proc_ranks2[i] = procs_global_ids2[i];

    if (current_proc_id1 != -1)
        MPI_Barrier(comp1);
    if (current_proc_id2 != -1)
        MPI_Barrier(comm2);

    MPI_Comm_group(MPI_COMM_WORLD, &common_group);
    MPI_Group_incl(common_group, comp_num_procs1, proc_ranks1, &comm_group1);
    MPI_Group_incl(common_group, comp_num_procs2, proc_ranks2, &comm_group2);
    MPI_Group_intersection(comm_group1, comm_group2, &intersection_group);
    MPI_Group_size(intersection_group, &intersection_size);

    if (intersection_size == comp_num_procs1) 
        union_comm = comm2;
    else if (intersection_size == comp_num_procs2)
        union_comm = comp1;
    else if (intersection_size == 0) {
        if (current_proc_id1 != -1) {
            MPI_Intercomm_create(comp1, 0, MPI_COMM_WORLD, procs_global_ids2[0], connection_id, &inter_comm);
            MPI_Intercomm_merge(inter_comm, true, &union_comm);
        }
        else if (current_proc_id2 != -1) {
            MPI_Intercomm_create(comm2, 0, MPI_COMM_WORLD, procs_global_ids1[0], connection_id, &inter_comm);
            MPI_Intercomm_merge(inter_comm, true, &union_comm);
        }
    }
    else {
        int * translate_ranks = new int[comp_num_procs2];
        for (int i = 0; i < comp_num_procs2; i ++) proc_ranks2[i] = i;
        MPI_Group_translate_ranks(comm_group2, comp_num_procs2, proc_ranks2, comm_group1, translate_ranks);
        if (current_proc_id2 != -1) {
            int color = 1;
            if (current_proc_id1 != -1) 
                color = 0;
            MPI_Comm_split(comm2, color, 0, &exclusive_comm);
        }
        if (current_proc_id1 != -1) {
            int root_indx_in_exclusive_comm = -1;
            for (int i = 0; i < comp_num_procs2; i ++)
                if (translate_ranks[i] < 0) {
                    root_indx_in_exclusive_comm = i;
                    break;
                }
            MPI_Intercomm_create(comp1, 0, MPI_COMM_WORLD, procs_global_ids2[root_indx_in_exclusive_comm], connection_id, &inter_comm);
            MPI_Intercomm_merge(inter_comm, true, &union_comm);
        }
        else {
            MPI_Intercomm_create(exclusive_comm, 0, MPI_COMM_WORLD, procs_global_ids1[0], connection_id, &inter_comm);
            MPI_Intercomm_merge(inter_comm, true, &union_comm);
        }
        delete [] translate_ranks;
    }

    for (int i = 0; i < comp_num_procs1; i ++) 
        proc_ranks1[i] = i;
    for (int i = 0; i < comp_num_procs2; i ++) 
        proc_ranks2[i] = i;

    MPI_Comm_group(union_comm, &union_group);
    if (proc_ranks_in_union_comm1 != NULL)
        MPI_Group_translate_ranks(comm_group1, comp_num_procs1, proc_ranks1, union_group, proc_ranks_in_union_comm1);
    if (proc_ranks_in_union_comm2 != NULL)
        MPI_Group_translate_ranks(comm_group2, comp_num_procs2, proc_ranks2, union_group, proc_ranks_in_union_comm2);

    delete [] proc_ranks1;
    delete [] proc_ranks2;

    return union_comm;
}


Coupling_connection::Coupling_connection(int id)
{
    import_interface = NULL;
    export_interface = NULL;
    import_procedure = NULL;
    export_procedure = NULL;
    connection_id = id;
    union_comm = MPI_COMM_NULL;
    src_proc_ranks_in_union_comm = NULL;
    dst_proc_ranks_in_union_comm = NULL;
    if (connection_id > coupling_generator->get_latest_connection_id())
        coupling_generator->set_latest_connection_id(connection_id);
}


Coupling_connection::~Coupling_connection()
{
    clean_string_pair_vector(src_comp_interfaces);
}


void Coupling_connection::generate_a_coupling_procedure(bool has_frac_remapping)
{
    src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_comp_interfaces[0].first);
    dst_comp_node =comp_comm_group_mgt_mgr->search_global_node(dst_comp_full_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, src_comp_node != NULL, "Software error in Coupling_connection::generate_a_coupling_procedure: Cannot find the src comp_node for \"%s\"", src_comp_interfaces[0].first);
    EXECUTION_REPORT(REPORT_ERROR, -1, dst_comp_node != NULL, "Software error in Coupling_connection::generate_a_coupling_procedure: Cannot find the dst comp_node for \"%s\"", dst_comp_full_name);
    current_proc_id_src_comp = src_comp_node->get_current_proc_local_id();
    current_proc_id_dst_comp = dst_comp_node->get_current_proc_local_id();
    src_comp_root_proc_global_id = src_comp_node->get_root_proc_global_id();
    dst_comp_root_proc_global_id = dst_comp_node->get_root_proc_global_id();
    
    if (current_proc_id_src_comp == -1 && current_proc_id_dst_comp == -1)
        return;

    if (current_proc_id_src_comp != -1) {
        export_interface = inout_interface_mgr->get_interface(src_comp_interfaces[0].first, src_comp_interfaces[0].second);
        EXECUTION_REPORT(REPORT_ERROR, -1, export_interface != NULL, "Software error in Coupling_connection::generate_a_coupling_procedure: NULL export interface");
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "start to generate a coupling connection from \"%s\" (current component) to \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    }
    if (current_proc_id_dst_comp != -1) {
        import_interface = inout_interface_mgr->get_interface(dst_comp_full_name, dst_interface_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, import_interface != NULL, "Software error in Coupling_connection::generate_a_coupling_procedure: NULL import interface");
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "start to generate a coupling connection from \"%s\" to \"%s\" (current component). The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    }

    create_union_comm();
    exchange_connection_fields_info();
    generate_interpolation(has_frac_remapping);

    if (current_proc_id_src_comp != -1) {
        export_procedure = new Connection_coupling_procedure(export_interface, this);
        export_interface->add_coupling_procedure(export_procedure);
    }
    if (current_proc_id_dst_comp != -1) {
        import_procedure = new Connection_coupling_procedure(import_interface, this);
        import_interface->add_coupling_procedure(import_procedure);
    }

    generate_data_transfer();

    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "Finish generating a coupling connection from \"%s\" (current component) to \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "Finish generating a coupling connection from \"%s\" to \"%s\" (current component). The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
}


void Coupling_connection::create_union_comm()
{
    int src_comp_num_procs = src_comp_node->get_num_procs();
    int dst_comp_num_procs = dst_comp_node->get_num_procs();
    std::vector<int> src_procs_global_ids, dst_procs_global_ids;


    if (current_proc_id_src_comp != -1) {
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "start to create union comm between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
        MPI_Barrier(src_comp_node->get_comm_group());
    }
    if (current_proc_id_dst_comp != -1) {
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "start to create union comm between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
        MPI_Barrier(dst_comp_node->get_comm_group());
    }

    for (int i = 0; i < src_comp_num_procs; i ++)
        src_procs_global_ids.push_back(src_comp_node->get_local_proc_global_id(i));
    for (int i = 0; i < dst_comp_num_procs; i ++)
        dst_procs_global_ids.push_back(dst_comp_node->get_local_proc_global_id(i));
    if (src_proc_ranks_in_union_comm == NULL) 
        src_proc_ranks_in_union_comm = new int[src_comp_num_procs]; 
    if (dst_proc_ranks_in_union_comm == NULL) 
        dst_proc_ranks_in_union_comm = new int[dst_comp_num_procs]; 

    union_comm = create_union_comm_common(src_comp_node->get_comm_group(), dst_comp_node->get_comm_group(), current_proc_id_src_comp, current_proc_id_dst_comp, src_procs_global_ids, dst_procs_global_ids, connection_id, src_proc_ranks_in_union_comm, dst_proc_ranks_in_union_comm);
    
    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "Finish creating union comm between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "Finish creating union comm between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
}


void Coupling_connection::generate_data_transfer()
{
    Field_mem_info **src_fields_mem = new Field_mem_info *[src_fields_info.size()];
    Field_mem_info **dst_fields_mem = new Field_mem_info *[src_fields_info.size()];
    Routing_info **fields_router = new Routing_info *[src_fields_info.size()];
    Runtime_trans_algorithm * send_algorithm_object = NULL;
    Runtime_trans_algorithm * recv_algorithm_object = NULL;
    MPI_Win data_win, tag_win;
    int dst_comp_id;
    long content_size = NAME_STR_SIZE;
    char *temp_dst_decomp_name = new char [NAME_STR_SIZE];


    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "Start to generate runtime data transfer algorithm from component \"%s\" (current component) to \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "Start to generate runtime data transfer algorithm from component \"%s\" to \"%s\" (current component). The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);

    for (int i = 0; i < src_fields_info.size(); i ++) { 
        if (dst_fields_info[i]->runtime_remapping_weights != NULL && dst_fields_info[i]->runtime_remapping_weights->get_src_decomp_info() != NULL) {
            strcpy(temp_dst_decomp_name, dst_fields_info[i]->runtime_remapping_weights->get_src_decomp_info()->get_decomp_name());
            dst_comp_id = dst_fields_info[i]->runtime_remapping_weights->get_src_decomp_info()->get_comp_id();
        }
        else {
            strcpy(temp_dst_decomp_name, dst_fields_info[i]->decomp_name);
            dst_comp_id = dst_comp_node->get_comp_id();
        }
        transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), &temp_dst_decomp_name, content_size);
        fields_router[i] = routing_info_mgr->search_or_add_router(src_comp_node->get_comp_id(), dst_comp_id, src_fields_info[i]->decomp_name, temp_dst_decomp_name);
        if (current_proc_id_src_comp != -1)
            src_fields_mem[i] = export_procedure->get_data_transfer_field_instance(i);
        if (current_proc_id_dst_comp != -1)
            dst_fields_mem[i] = import_procedure->get_data_transfer_field_instance(i);
    }

    if (current_proc_id_dst_comp != -1) {
        recv_algorithm_object = new Runtime_trans_algorithm(false, src_fields_info.size(), dst_fields_mem, fields_router, union_comm, src_proc_ranks_in_union_comm, connection_id);
        import_procedure->add_data_transfer_algorithm(recv_algorithm_object);
        inout_interface_mgr->add_runtime_receive_algorithm(recv_algorithm_object);
    }
    if (current_proc_id_src_comp != -1) {
        send_algorithm_object = new Runtime_trans_algorithm(true, src_fields_info.size(), src_fields_mem, fields_router, union_comm, dst_proc_ranks_in_union_comm, connection_id);
        export_procedure->add_data_transfer_algorithm(send_algorithm_object);
    }
    if (current_proc_id_dst_comp != -1) {
        MPI_Win_create(recv_algorithm_object->get_total_buf(), recv_algorithm_object->get_total_buf_size()*sizeof(char), sizeof(char), MPI_INFO_NULL, union_comm, &data_win);
        //MPI_Win_create(recv_algorithm_object->get_data_buf(), recv_algorithm_object->get_data_buf_size()*sizeof(char), sizeof(char), MPI_INFO_NULL, union_comm, &data_win);
        //MPI_Win_create(recv_algorithm_object->get_tag_buf(), recv_algorithm_object->get_tag_buf_size()*sizeof(long), sizeof(long), MPI_INFO_NULL, union_comm, &tag_win);
        recv_algorithm_object->set_data_win(data_win);
        //recv_algorithm_object->set_tag_win(tag_win);
        inout_interface_mgr->add_MPI_win(data_win);
        //inout_interface_mgr->add_MPI_win(tag_win);
    }
    if (current_proc_id_src_comp != -1) {
        if (current_proc_id_dst_comp == -1) {
            MPI_Win_create(NULL, 0, sizeof(char), MPI_INFO_NULL, union_comm, &data_win);
            //MPI_Win_create(NULL, 0, sizeof(long), MPI_INFO_NULL, union_comm, &tag_win);
            inout_interface_mgr->add_MPI_win(data_win);
            //inout_interface_mgr->add_MPI_win(tag_win);
        }
        send_algorithm_object->set_data_win(data_win);
        //send_algorithm_object->set_tag_win(tag_win);
    }

    delete [] src_fields_mem;
    delete [] dst_fields_mem;
    delete [] fields_router;
    delete [] temp_dst_decomp_name;

    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "Finish generating runtime data transfer algorithm from component \"%s\" (current component) to \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "Finish generating runtime data transfer algorithm from component \"%s\" to \"%s\" (current component). The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
}


bool Coupling_connection::exchange_grid(Comp_comm_group_mgt_node *sender_comp_node, Comp_comm_group_mgt_node *receiver_comp_node, const char *grid_name)
{
    char *temp_array_buffer = NULL;
    long buffer_max_size, buffer_content_size;
    int original_grid_status, *all_original_grid_status, num_processes, bottom_field_variation_type, V3D_lev_field_variation_type;
    long checksum_lon, checksum_lat, checksum_mask;
    bool should_exchange_grid = false;


    Original_grid_info *sender_original_grid = original_grid_mgr->search_grid_info(grid_name, sender_comp_node->get_comp_id());
    Original_grid_info *receiver_original_grid = original_grid_mgr->search_grid_info(grid_name, receiver_comp_node->get_comp_id());
    
    original_grid_status = 0;
    if (sender_original_grid != NULL && receiver_original_grid != NULL && sender_original_grid->get_original_CoR_grid() == receiver_original_grid->get_original_CoR_grid())
        original_grid_status = 1;
    MPI_Comm_size(union_comm, &num_processes);
    all_original_grid_status = new int [num_processes];
    MPI_Allgather(&original_grid_status, 1, MPI_INT, all_original_grid_status, 1, MPI_INT, union_comm);
    for (int i = 0; i < num_processes; i ++)
        if (all_original_grid_status[i] == 1) {
            delete [] all_original_grid_status;
            if (sender_comp_node->get_current_proc_local_id() != -1)
                EXECUTION_REPORT_LOG(REPORT_LOG, sender_comp_node->get_comp_id(), true, "Does not exchange grid \"%s\" from \"%s\" to \"%s\" because the CoR grid is the same", grid_name, sender_comp_node->get_comp_full_name(), receiver_comp_node->get_comp_full_name());
            if (receiver_comp_node->get_current_proc_local_id() != -1)
                EXECUTION_REPORT_LOG(REPORT_LOG, receiver_comp_node->get_comp_id(), true, "Does not exchange grid \"%s\" from \"%s\" to \"%s\" because the CoR grid is the same", grid_name, sender_comp_node->get_comp_full_name(), receiver_comp_node->get_comp_full_name());
            return false;
        }
    
    original_grid_status = sender_original_grid == NULL? 0 : 1;
    MPI_Allgather(&original_grid_status, 1, MPI_INT, all_original_grid_status, 1, MPI_INT, union_comm);
    for (int i = 0; i < num_processes; i ++)
        if (all_original_grid_status[i] == 0) {
            should_exchange_grid = true;
            break;
        }
    delete [] all_original_grid_status;
    if (!should_exchange_grid) {
        if (sender_comp_node->get_current_proc_local_id() != -1)
            EXECUTION_REPORT_LOG(REPORT_LOG, sender_comp_node->get_comp_id(), true, "Does not exchange grid \"%s\" from \"%s\" to \"%s\" again", grid_name, sender_comp_node->get_comp_full_name(), receiver_comp_node->get_comp_full_name());
        if (receiver_comp_node->get_current_proc_local_id() != -1)
            EXECUTION_REPORT_LOG(REPORT_LOG, receiver_comp_node->get_comp_id(), true, "Does not exchange grid \"%s\" from \"%s\" to \"%s\" again", grid_name, sender_comp_node->get_comp_full_name(), receiver_comp_node->get_comp_full_name());
        return true;
    }

    if (sender_comp_node->get_current_proc_local_id() != -1) 
        EXECUTION_REPORT_LOG(REPORT_LOG, sender_comp_node->get_comp_id(), true, "Send grid %s to component \"%s\"", grid_name, receiver_comp_node->get_full_name());
    if (receiver_comp_node->get_current_proc_local_id() != -1) 
        EXECUTION_REPORT_LOG(REPORT_LOG, receiver_comp_node->get_comp_id(), true, "Receive grid %s from component \"%s\"", grid_name, sender_comp_node->get_full_name());

    if (sender_original_grid != NULL)
        sender_original_grid->write_grid_into_array(&temp_array_buffer, buffer_max_size, buffer_content_size);
    transfer_array_from_one_comp_to_another(sender_comp_node->get_current_proc_local_id(), sender_comp_node->get_root_proc_global_id(), receiver_comp_node->get_current_proc_local_id(), receiver_comp_node->get_root_proc_global_id(), receiver_comp_node->get_comm_group(), &temp_array_buffer, buffer_content_size);

    if (original_grid_status == 0) {
        read_data_from_array_buffer(&checksum_mask, sizeof(long), temp_array_buffer, buffer_content_size, true);
        read_data_from_array_buffer(&bottom_field_variation_type, sizeof(int), temp_array_buffer, buffer_content_size, true);
        read_data_from_array_buffer(&V3D_lev_field_variation_type, sizeof(int), temp_array_buffer, buffer_content_size, true);
        Remap_grid_class *mirror_grid = new Remap_grid_class(NULL, sender_comp_node->get_full_name(), temp_array_buffer, buffer_content_size);
        mirror_grid = remap_grid_manager->search_remap_grid_with_grid_name(mirror_grid->get_grid_name());
        EXECUTION_REPORT(REPORT_ERROR, -1, buffer_content_size == 0, "software error in Coupling_connection::exchange_grid: wrong buffer_content_size");
        receiver_original_grid = original_grid_mgr->get_original_grid(original_grid_mgr->add_original_grid(sender_comp_node->get_comp_id(), grid_name, mirror_grid));
        if (receiver_original_grid->get_bottom_field_variation_type() != bottom_field_variation_type)
            EXECUTION_REPORT(REPORT_ERROR, -1, receiver_original_grid->get_original_CoR_grid()->is_sigma_grid(), "Software error in Coupling_connection::exchange_grid regarding bottom_field_variation_type");
        receiver_original_grid->set_bottom_field_variation_type(bottom_field_variation_type);
		receiver_original_grid->set_V3D_lev_field_variation_type(V3D_lev_field_variation_type);
        receiver_original_grid->set_grid_checksum(checksum_mask);
    }

    if (temp_array_buffer != NULL)
        delete [] temp_array_buffer;

    return true;
}


void Coupling_connection::exchange_remapping_setting(int i, Remapping_setting &field_remapping_setting)
{
    char *array = NULL;
    long buffer_max_size, buffer_content_size;
    
    if (current_proc_id_src_comp != -1) {
        remapping_configuration_mgr->get_field_remapping_setting(field_remapping_setting, src_comp_node->get_comp_id(), fields_name[i]);
        field_remapping_setting.write_remapping_setting_into_array(&array, buffer_max_size, buffer_content_size);
    }
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), &array, buffer_content_size);    
    if (current_proc_id_dst_comp != -1)
        field_remapping_setting.read_remapping_setting_from_array(array, buffer_content_size);

    delete [] array;
}


void Coupling_connection::add_bottom_field_coupling_info(int field_connection_indx, Runtime_remapping_weights *V3D_remapping_weights, Remapping_setting *remapping_setting)
{
    for (int i = 0; i < dst_bottom_fields_coupling_info.size(); i ++)
        if (dst_bottom_fields_coupling_info[i]->V3D_runtime_remapping_weights == V3D_remapping_weights)
            return;

    V3D_grid_bottom_field_coupling_info *bottom_field_coupling_info = new V3D_grid_bottom_field_coupling_info;
    bottom_field_coupling_info->V3D_runtime_remapping_weights = V3D_remapping_weights;
    bottom_field_coupling_info->field_connection_indx = field_connection_indx;
	if (V3D_remapping_weights->get_src_original_grid()->get_original_CoR_grid()->is_sigma_grid())
	    bottom_field_coupling_info->is_dynamic_bottom_field = V3D_remapping_weights->get_src_original_grid()->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_DYNAMIC;
	else bottom_field_coupling_info->is_dynamic_bottom_field = V3D_remapping_weights->get_src_original_grid()->get_V3D_lev_field_variation_type() == BOTTOM_FIELD_VARIATION_DYNAMIC;
    bottom_field_coupling_info->bottom_field_inst = V3D_remapping_weights->allocate_intermediate_V3D_grid_bottom_field();
	if (V3D_remapping_weights->get_src_original_grid()->get_original_CoR_grid()->is_sigma_grid())
	    bottom_field_coupling_info->H2D_runtime_remapping_weights = runtime_remapping_weights_mgr->search_or_generate_runtime_remapping_weights(src_comp_node->get_comp_full_name(), dst_comp_node->get_comp_full_name(), 
	        original_grid_mgr->get_original_grid(V3D_remapping_weights->get_src_decomp_info()->get_grid_id()), original_grid_mgr->get_original_grid(V3D_remapping_weights->get_dst_decomp_info()->get_grid_id()), 
    	    remapping_setting, V3D_remapping_weights->get_dst_decomp_info());
	else bottom_field_coupling_info->H2D_runtime_remapping_weights = V3D_remapping_weights;

    if (V3D_remapping_weights->get_parallel_remapping_weights() != NULL) {
        Remap_weight_of_operator_class *dynamic_V1D_remap_weight_of_operator = V3D_remapping_weights->get_parallel_remapping_weights()->get_dynamic_V1D_remap_weight_of_operator();
        EXECUTION_REPORT(REPORT_ERROR, -1, dynamic_V1D_remap_weight_of_operator != NULL, "Software error in Coupling_connection::add_bottom_field_coupling_info: do not have dynamic_V1D_remap_weight_of_operator");
        EXECUTION_REPORT(REPORT_ERROR, -1, !bottom_field_coupling_info->bottom_field_inst->get_field_data()->get_coord_value_grid()->is_sigma_grid(), "Software error in Coupling_connection::add_bottom_field_coupling_info: wrong surface field grid or wrong 2-D+1-D order");
        if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->get_level_V3D_coord_dynamic_trigger_field() != NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->get_level_V3D_coord_dynamic_trigger_field() == bottom_field_coupling_info->bottom_field_inst->get_field_data(), "Software error in Coupling_connection::add_bottom_field_coupling_info: the surface field of the same grid has been set to different data fields");
        else dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->set_level_V3D_coord_dynamic_trigger_field(bottom_field_coupling_info->bottom_field_inst->get_field_data());

        if (V3D_remapping_weights->get_dst_original_grid()->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_EXTERNAL) {
            if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_level_V3D_coord_dynamic_trigger_field() != NULL)
                EXECUTION_REPORT(REPORT_ERROR, -1, dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_level_V3D_coord_dynamic_trigger_field() == bottom_field_coupling_info->bottom_field_inst->get_field_data(), "Software error in Coupling_connection::add_bottom_field_coupling_info: the surface field of the same grid has been set to different data fields");
            else dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->set_level_V3D_coord_dynamic_trigger_field(bottom_field_coupling_info->bottom_field_inst->get_field_data());
        }
    }

    dst_bottom_fields_coupling_info.push_back(bottom_field_coupling_info);
}


void Coupling_connection::generate_src_bottom_field_coupling_info()
{
    int *bottom_fields_indx = NULL;
    long buffer_max_size = 0, buffer_content_size = 0;


    if (current_proc_id_dst_comp != -1)
        for (int i = dst_bottom_fields_coupling_info.size()-1; i >= 0; i --)
            write_data_into_array_buffer(&(dst_bottom_fields_coupling_info[i]->field_connection_indx), sizeof(int), (char**)(&bottom_fields_indx), buffer_max_size, buffer_content_size);
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), (char**)(&bottom_fields_indx), buffer_content_size);
    if (current_proc_id_src_comp != -1) {
        for (int i = 0; i < buffer_content_size / 4; i ++) {
            V3D_grid_bottom_field_coupling_info *bottom_field_coupling_info = new V3D_grid_bottom_field_coupling_info;
            bottom_field_coupling_info->V3D_runtime_remapping_weights = NULL;
            bottom_field_coupling_info->H2D_runtime_remapping_weights = NULL; 
            bottom_field_coupling_info->field_connection_indx = bottom_fields_indx[i];
            Original_grid_info *src_original_grid = original_grid_mgr->search_grid_info(src_fields_info[bottom_fields_indx[i]]->grid_name, comp_comm_group_mgt_mgr->search_global_node(src_comp_interfaces[0].first)->get_comp_id());
			if (src_original_grid->get_bottom_field_id() != -1) {				
				EXECUTION_REPORT(REPORT_ERROR, -1, src_original_grid->is_3D_grid() && src_original_grid->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_DYNAMIC || src_original_grid->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_STATIC, "Software error in Coupling_connection::generate_src_bottom_field_coupling_info: wrong grid");
				bottom_field_coupling_info->is_dynamic_bottom_field = src_original_grid->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_DYNAMIC;
	            bottom_field_coupling_info->bottom_field_inst = memory_manager->get_field_instance(src_original_grid->get_bottom_field_id());
			}
			else if (src_original_grid->get_V3D_lev_field_id() != -1) {
				EXECUTION_REPORT(REPORT_ERROR, -1, src_original_grid->get_V3D_lev_field_variation_type() == BOTTOM_FIELD_VARIATION_DYNAMIC || src_original_grid->get_V3D_lev_field_variation_type() == BOTTOM_FIELD_VARIATION_STATIC, "Software error in Coupling_connection::generate_src_bottom_field_coupling_info: wrong grid");
				bottom_field_coupling_info->is_dynamic_bottom_field = src_original_grid->get_V3D_lev_field_variation_type() == BOTTOM_FIELD_VARIATION_DYNAMIC;
	            bottom_field_coupling_info->bottom_field_inst = memory_manager->get_field_instance(src_original_grid->get_V3D_lev_field_id());
			}
			else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Coupling_connection::generate_src_bottom_field_coupling_info");
            EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->get_original_CoR_grid(bottom_field_coupling_info->bottom_field_inst->get_grid_id())->is_subset_of_grid(src_original_grid->get_original_CoR_grid()), "Software error in Coupling_connection::generate_src_bottom_field_coupling_info: wrong grid relation");
            src_bottom_fields_coupling_info.push_back(bottom_field_coupling_info);
        }        
    }

    if (bottom_fields_indx != NULL)
        delete [] bottom_fields_indx;
}


void Coupling_connection::generate_interpolation(bool has_frac_remapping)
{
    Original_grid_info *src_original_grid = NULL, *dst_original_grid = NULL;
    int i, j;

    
    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "start to generate interpolation between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "start to generate interpolation between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);

    for (int i = 0; i < fields_name.size(); i ++) {
        src_fields_info[i]->runtime_remapping_weights = NULL;
        dst_fields_info[i]->runtime_remapping_weights = NULL;
        if (words_are_the_same(dst_fields_info[i]->grid_name, "NULL")) {
			if (current_proc_id_dst_comp != -1)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, dst_comp_node->get_comp_id(), words_are_the_same(src_fields_info[i]->grid_name, "NULL"), "Error happens when trying to generate data interpolation for the field \"%s\" in the coupling from the component model \"%s\" (interface is \"%s\") to the component model \"%s\" (interface is \"%s\"): the field is scalar (not on a grid) on the target component model while is on a grid \"%s\" on the source component model. Please verify.", fields_name[i], src_comp_interfaces[0].first, src_comp_interfaces[0].second, dst_comp_full_name, dst_interface_name, dst_fields_info[i]->grid_name);
            continue;
        }
        if (words_are_the_same(src_fields_info[i]->grid_name, "NULL")) {
			if (current_proc_id_dst_comp != -1)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, dst_comp_node->get_comp_id(), words_are_the_same(dst_fields_info[i]->grid_name, "NULL"), "Error happens when trying to generate data interpolation for the field \"%s\" in the coupling from the component model \"%s\" (interface is \"%s\") to the component model \"%s\" (interface is \"%s\"): the field is scalar (not on a grid) on the source component model while is on a grid \"%s\" on the target component model. Please verify.", fields_name[i], src_comp_interfaces[0].first, src_comp_interfaces[0].second, dst_comp_full_name, dst_interface_name, dst_fields_info[i]->grid_name);
            continue;
        }
        if (src_comp_node == dst_comp_node && words_are_the_same(src_fields_info[i]->grid_name, dst_fields_info[i]->grid_name))
            continue;
        Remapping_setting field_remapping_setting;    
        if (words_are_the_same(fields_name[i], "remap_frac")) {
            EXECUTION_REPORT(REPORT_ERROR, -1, fields_name.size() > 1 && i == fields_name.size()-1, "Software error in Coupling_connection::generate_interpolation");
            exchange_remapping_setting(0, field_remapping_setting);
        }
        else exchange_remapping_setting(i, field_remapping_setting);
//        exchange_grid(dst_comp_node, src_comp_node, dst_fields_info[i]->grid_name);
        exchange_grid(src_comp_node, dst_comp_node, src_fields_info[i]->grid_name);
        if (current_proc_id_dst_comp != -1) {
            dst_original_grid = original_grid_mgr->search_grid_info(dst_fields_info[i]->grid_name, comp_comm_group_mgt_mgr->search_global_node(dst_comp_full_name)->get_comp_id());
            src_original_grid = original_grid_mgr->search_grid_info(src_fields_info[i]->grid_name, comp_comm_group_mgt_mgr->search_global_node(src_comp_interfaces[0].first)->get_comp_id());
            if (src_original_grid->is_H2D_grid_and_the_same_as_another_grid(dst_original_grid)) {
                EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "The data interpolation from grid \"%s\" to \"%s\" is bypassed as these too grids are the same", src_original_grid->get_grid_name(), dst_original_grid->get_grid_name());
                continue;
            }
            else EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "The data interpolation from grid \"%s\" to \"%s\" is necessary", src_original_grid->get_grid_name(), dst_original_grid->get_grid_name());
            EXECUTION_REPORT(REPORT_ERROR, dst_comp_node->get_comp_id(), words_are_the_same(src_fields_info[i]->data_type, DATA_TYPE_FLOAT) || words_are_the_same(src_fields_info[i]->data_type, DATA_TYPE_DOUBLE),  "Error happens when trying to generate data interpolation for the field \%s\" in the coupling from the component model \"%s\" (interface is \"%s\") to the component model \"%s\" (interface is \"%s\"): the data type in the source component model is not floating-point but %s. Please verify.", fields_name[i], src_comp_interfaces[0].first, src_comp_interfaces[0].second, dst_comp_full_name, dst_interface_name, src_fields_info[i]->data_type);
            EXECUTION_REPORT(REPORT_ERROR, dst_comp_node->get_comp_id(), words_are_the_same(dst_fields_info[i]->data_type, DATA_TYPE_FLOAT) || words_are_the_same(dst_fields_info[i]->data_type, DATA_TYPE_DOUBLE),  "Error happens when trying to generate data interpolation for the field \%s\" in the coupling from the component model \"%s\" (interface is \"%s\") to the component model \"%s\" (interface is \"%s\"): the data type in the source component model is not floating-point but %s. Please verify.", fields_name[i], src_comp_interfaces[0].first, src_comp_interfaces[0].second, dst_comp_full_name, dst_interface_name, dst_fields_info[i]->data_type);    
            if (src_original_grid->get_original_CoR_grid()->is_sigma_grid()) {
                EXECUTION_REPORT(REPORT_ERROR, dst_comp_node->get_comp_id(), src_original_grid->get_bottom_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET, "Fail to generate an interpolation from component model \"%s\" to \"%s\": the surface field of the source 3-D grid \"%s\" with SIGMA or HYBRID vertical coordinate has not been specified. Please verify.", src_comp_interfaces[0].first, dst_comp_full_name, src_original_grid->get_grid_name());
                EXECUTION_REPORT(REPORT_ERROR, dst_comp_node->get_comp_id(), src_original_grid->get_bottom_field_variation_type() != BOTTOM_FIELD_VARIATION_EXTERNAL, "Fail to generate an interpolation from component model \"%s\" to \"%s\": the source 3-D grid \"%s\" with SIGMA or HYBRID vertical coordinate has an external surface field while a 3-D grid with external surface field cannot be referred as a source grid in interpolation. Please verify.", src_comp_interfaces[0].first, dst_comp_full_name, src_original_grid->get_grid_name());
            }
            if (dst_original_grid->get_original_CoR_grid()->is_sigma_grid()) {
                EXECUTION_REPORT(REPORT_ERROR, dst_comp_node->get_comp_id(), dst_original_grid->get_bottom_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET, "Fail to generate an interpolation from component model \"%s\" to \"%s\": the surface field of the target 3-D grid \"%s\" that includes SIGMA or HYBRID vertical coordinate has not been specified. Please verify.", src_comp_interfaces[0].first, dst_comp_full_name, dst_original_grid->get_grid_name());
                if (dst_original_grid->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_EXTERNAL)
                    EXECUTION_REPORT(REPORT_ERROR, dst_comp_node->get_comp_id(), src_original_grid->get_original_CoR_grid()->is_sigma_grid(), "Fail to generate an interpolation from component model \"%s\" to \"%s\": the target 3-D grid \"%s\" has an external surface field; at this time, the source 3-D grid \"%s\" should include SIGMA or HYBRID vertical coordinate but actually not. Please verify. ", src_comp_interfaces[0].first, dst_comp_full_name, dst_original_grid->get_grid_name(), src_original_grid->get_grid_name());
            }    
            dst_fields_info[i]->runtime_remapping_weights = runtime_remapping_weights_mgr->search_or_generate_runtime_remapping_weights(src_comp_node->get_comp_full_name(), dst_comp_node->get_comp_full_name(), src_original_grid, dst_original_grid, &field_remapping_setting, decomps_info_mgr->search_decomp_info(dst_fields_info[i]->decomp_name, dst_comp_node->get_comp_id()));
            if (src_original_grid->get_original_CoR_grid()->is_sigma_grid() || src_original_grid->get_original_CoR_grid()->does_use_V3D_level_coord())
                add_bottom_field_coupling_info(i, dst_fields_info[i]->runtime_remapping_weights, &field_remapping_setting);
        }
    }
    
    int *remapping_weights_index_table = new int [fields_name.size()];
    long array_size = fields_name.size()*sizeof(int);
    if (current_proc_id_dst_comp != -1) 
        for (int i = 0; i < fields_name.size(); i ++)
            if (dst_fields_info[i]->runtime_remapping_weights == NULL || !dst_fields_info[i]->runtime_remapping_weights->get_src_original_grid()->is_H2D_grid())
                remapping_weights_index_table[i] = -1;
            else {
                for (j = 0; j < i; j ++)
                    if (dst_fields_info[i]->runtime_remapping_weights == dst_fields_info[j]->runtime_remapping_weights)
                        break;
                remapping_weights_index_table[i] = j;
            }
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_node->get_local_proc_global_id(0), current_proc_id_src_comp, src_comp_node->get_local_proc_global_id(0), src_comp_node->get_comm_group(), (char**)(&remapping_weights_index_table), array_size);
    for (int i = 0; i < fields_name.size(); i ++) {
        if (remapping_weights_index_table[i] == -1)
            continue;
        if (remapping_weights_index_table[i] != i) {
            if (current_proc_id_src_comp != -1)
                src_fields_info[i]->runtime_remapping_weights = src_fields_info[remapping_weights_index_table[i]]->runtime_remapping_weights;    
        }
        else runtime_remapping_weights_mgr->transfer_runtime_remapping_weights(dst_fields_info[i]->runtime_remapping_weights, &(src_fields_info[i]->runtime_remapping_weights), dst_comp_node, src_comp_node);
    }
    delete [] remapping_weights_index_table;

    generate_src_bottom_field_coupling_info();
    exchange_bottom_fields_info();

    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "finish generating interpolation between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "finish generating interpolation between components \"%s\" and \"%s\". The connection id is %d", src_comp_interfaces[0].first, dst_comp_full_name, connection_id);
}


void Coupling_connection::exchange_connection_fields_info()
{
    char *src_fields_info_array = NULL, *dst_fields_info_array = NULL;
    long src_fields_info_array_size, dst_fields_info_array_size, buffer_max_size, comp_id;


    if (current_proc_id_dst_comp != -1)
        write_connection_fields_info_into_array(import_interface, &dst_fields_info_array, buffer_max_size, dst_fields_info_array_size, &dst_timer, dst_inst_or_aver, dst_time_step_in_second, dst_current_year, dst_current_month, dst_current_day, dst_current_second);
    if (current_proc_id_src_comp != -1)
        write_connection_fields_info_into_array(export_interface, &src_fields_info_array, buffer_max_size, src_fields_info_array_size, &src_timer, src_inst_or_aver, src_time_step_in_second, src_current_year, src_current_month, src_current_day, src_current_second);
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), &src_fields_info_array, src_fields_info_array_size);
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), &dst_fields_info_array, dst_fields_info_array_size);
    comp_id = export_interface != NULL? export_interface->get_comp_id() : import_interface->get_comp_id();
    read_connection_fields_info_from_array(src_fields_info, src_fields_info_array, src_fields_info_array_size, comp_id, &src_timer, src_inst_or_aver, src_time_step_in_second, src_current_year, src_current_month, src_current_day, src_current_second);
    comp_id = import_interface != NULL? import_interface->get_comp_id() : export_interface->get_comp_id();
    read_connection_fields_info_from_array(dst_fields_info, dst_fields_info_array, dst_fields_info_array_size, comp_id, &dst_timer, dst_inst_or_aver, dst_time_step_in_second, dst_current_year, dst_current_month, dst_current_day, dst_current_second);
    EXECUTION_REPORT(REPORT_ERROR, -1, fields_name.size() == src_fields_info.size() && fields_name.size() == dst_fields_info.size(), "Software error in Coupling_connection::exchange_connection_fields_info");

    src_timer->reset_remote_lag_count();
    
    delete [] src_fields_info_array;
    delete [] dst_fields_info_array;
}


void Coupling_connection::read_fields_info_from_array(std::vector<Interface_field_info*> &fields_info, const char *array_buffer, long buffer_content_iter)
{
    while (buffer_content_iter > 0) {
        Interface_field_info *field_info = new Interface_field_info;
        field_info->bottom_field_indx = -1;
        field_info->runtime_remapping_weights = NULL;
        read_data_from_array_buffer(field_info->decomp_name, NAME_STR_SIZE, array_buffer, buffer_content_iter, true);
        read_data_from_array_buffer(field_info->grid_name, NAME_STR_SIZE, array_buffer, buffer_content_iter, true);
        read_data_from_array_buffer(field_info->unit, NAME_STR_SIZE, array_buffer, buffer_content_iter, true);
        read_data_from_array_buffer(field_info->data_type, NAME_STR_SIZE, array_buffer, buffer_content_iter, true);
        fields_info.push_back(field_info);
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, buffer_content_iter == 0, "Software error in Coupling_connection::read_fields_info_from_array: wrong buffer_content_iter");
}


void Coupling_connection::read_connection_fields_info_from_array(std::vector<Interface_field_info*> &fields_info, const char *array_buffer, long buffer_content_iter, int comp_id, Coupling_timer **timer, int &inst_or_aver, int &time_step_in_second,
                                                                                 int &current_year, int &current_month, int &current_day, int &current_second)
{
    bool successful;

    read_data_from_array_buffer(&current_second, sizeof(int), array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(&current_day, sizeof(int), array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(&current_month, sizeof(int), array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(&current_year, sizeof(int), array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(&time_step_in_second, sizeof(int), array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(&inst_or_aver, sizeof(int), array_buffer, buffer_content_iter, true);
    *timer = new Coupling_timer(array_buffer, buffer_content_iter, comp_id, true, successful);

    read_fields_info_from_array(fields_info, array_buffer, buffer_content_iter);
    EXECUTION_REPORT(REPORT_ERROR, -1, fields_info.size() == fields_name.size(), "Software error in Coupling_connection::read_connection_fields_info_from_array: wrong size of fields_info");
}


void Coupling_connection::write_field_info_into_array(Field_mem_info *field, char **array, long &buffer_max_size, long &buffer_content_size)
{
    char tmp_string[NAME_STR_SIZE];

    
    EXECUTION_REPORT(REPORT_ERROR, -1, field != NULL, "Software error in Coupling_connection::write_field_info_into_array");
    write_data_into_array_buffer(field->get_field_data()->get_grid_data_field()->data_type_in_application, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(field->get_unit(), NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    const char *grid_name = field->get_grid_name();
    if (grid_name == NULL) {
        strcpy(tmp_string, "NULL");
        write_data_into_array_buffer(tmp_string, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    }
    else write_data_into_array_buffer(grid_name, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    const char *decomp_name = field->get_decomp_name();
    if (decomp_name == NULL) {
        strcpy(tmp_string, "NULL");
        write_data_into_array_buffer(tmp_string, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
    }
    else write_data_into_array_buffer(decomp_name, NAME_STR_SIZE, array, buffer_max_size, buffer_content_size);
}
 

void Coupling_connection::write_connection_fields_info_into_array(Inout_interface *inout_interface, char **array, long &buffer_max_size, long &buffer_content_size, Coupling_timer **timer, int &inst_or_aver, int &time_step_in_second, 
                                                                                   int &current_year, int &current_month, int &current_day, int &current_second)
{
    char tmp_string[NAME_STR_SIZE];
    int field_local_index;

    
    for (int i = fields_name.size() - 1; i >= 0; i --) {
        Field_mem_info *field = inout_interface->search_registered_field_instance(fields_name[i], field_local_index);
        if (field == NULL)
            EXECUTION_REPORT(REPORT_ERROR, inout_interface->get_comp_id(), field != NULL, "Software error in Coupling_connection::write_connection_fields_info_into_array: %s %s", inout_interface->get_interface_name(), fields_name[i]);
        write_field_info_into_array(field, array, buffer_max_size, buffer_content_size);
    }
    *timer = inout_interface->get_timer();
    inst_or_aver = inout_interface->get_inst_or_aver();
    time_step_in_second = components_time_mgrs->get_time_mgr(inout_interface->get_comp_id())->get_time_step_in_second();
    components_time_mgrs->get_time_mgr(inout_interface->get_comp_id())->get_current_time(current_year, current_month, current_day, current_second, 0, "in Coupling_connection::write_connection_fields_info_into_array");
    (*timer)->write_timer_into_array(array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&inst_or_aver, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&time_step_in_second, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&current_year, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&current_month, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&current_day, sizeof(int), array, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&current_second, sizeof(int), array, buffer_max_size, buffer_content_size);
}


void Coupling_connection::exchange_bottom_fields_info()
{
    char *src_fields_info_array = NULL, *dst_fields_info_array = NULL;
    long src_fields_info_array_size = 0, dst_fields_info_array_size = 0, buffer_max_size;


    if (current_proc_id_dst_comp != -1)
        for (int i = dst_bottom_fields_coupling_info.size() - 1; i >= 0; i --)
            write_field_info_into_array(dst_bottom_fields_coupling_info[i]->bottom_field_inst, &dst_fields_info_array, buffer_max_size, dst_fields_info_array_size);
    if (current_proc_id_src_comp != -1)
        for (int i = src_bottom_fields_coupling_info.size() - 1; i >= 0; i --)
            write_field_info_into_array(src_bottom_fields_coupling_info[i]->bottom_field_inst, &src_fields_info_array, buffer_max_size, src_fields_info_array_size);

    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), &src_fields_info_array, src_fields_info_array_size);
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), &dst_fields_info_array, dst_fields_info_array_size);

    read_fields_info_from_array(src_fields_info, src_fields_info_array, src_fields_info_array_size);
    read_fields_info_from_array(dst_fields_info, dst_fields_info_array, dst_fields_info_array_size);

    EXECUTION_REPORT(REPORT_ERROR, -1, dst_fields_info.size() == src_fields_info.size(), "Software error in Coupling_connection::exchange_bottom_fields_info");
    
    for (int i = fields_name.size(); i < dst_fields_info.size(); i ++) {
        dst_fields_info[i]->bottom_field_indx = i - fields_name.size();
        src_fields_info[i]->bottom_field_indx = dst_fields_info[i]->bottom_field_indx;
        if (current_proc_id_dst_comp != -1)
            dst_fields_info[i]->runtime_remapping_weights = dst_bottom_fields_coupling_info[dst_fields_info[i]->bottom_field_indx]->H2D_runtime_remapping_weights;
    }

    if (src_fields_info_array != NULL)
        delete [] src_fields_info_array;
    if (dst_fields_info_array != NULL)
        delete [] dst_fields_info_array;
}


Field_mem_info *Coupling_connection::get_bottom_field(bool export_or_import, int bottom_field_indx)
{
    if (export_or_import) {
        EXECUTION_REPORT(REPORT_ERROR, -1, bottom_field_indx < src_bottom_fields_coupling_info.size(), "Software error in Coupling_connection::get_bottom_field: wrong bottom_field_indx");
        return src_bottom_fields_coupling_info[bottom_field_indx]->bottom_field_inst;
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, bottom_field_indx < dst_bottom_fields_coupling_info.size(), "Software error in Coupling_connection::get_bottom_field: wrong bottom_field_indx");

        return dst_bottom_fields_coupling_info[bottom_field_indx]->bottom_field_inst;
    }
}


bool Coupling_connection::get_is_bottom_field_dynamic(int field_indx)
{
    if (current_proc_id_src_comp != -1)
        return src_bottom_fields_coupling_info[src_fields_info[field_indx]->bottom_field_indx]->is_dynamic_bottom_field;

    return dst_bottom_fields_coupling_info[dst_fields_info[field_indx]->bottom_field_indx]->is_dynamic_bottom_field;
}


Import_direction_setting::Import_direction_setting(int host_comp_id, Import_interface_configuration *interface_configuration, const char *comp_full_name, const char *interface_name, TiXmlElement *redirection_element, const char *XML_file_name, std::vector<const char*> &interface_fields_name, int *fields_count, bool check_comp_existence)
{
    TiXmlElement *fields_element = NULL, *components_element = NULL, *remapping_element = NULL, *merge_element = NULL;
    int i, line_number;
    std::pair<const char*, const char*> producer_info;


    strcpy(this->interface_name, interface_name);
    for (TiXmlNode *detailed_element_node = redirection_element->FirstChild(); detailed_element_node != NULL; detailed_element_node = detailed_element_node->NextSibling()) {        
        if (detailed_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *detailed_element = detailed_element_node->ToElement();
        if (words_are_the_same(detailed_element->Value(), "fields")) {
            if (!is_XML_setting_on(host_comp_id, detailed_element, XML_file_name, "the status of \"fields\"", "import interface configuration file"))
                continue;
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, fields_element == NULL, "When setting a coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the attribute of \"fields\" can be set only once but currently has been set at least twice. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, detailed_element->Row());
            fields_element = detailed_element;
            const char *default_str = get_XML_attribute(host_comp_id, -1, detailed_element, "default", XML_file_name, line_number, "default setting of \"fields\"", "import interface configuration file", true);
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, words_are_the_same(default_str, "off") || words_are_the_same(default_str, "all") || words_are_the_same(default_str, "remain"), "When setting a coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the value of \"default\" for the attribute of \"fields\" is wrong (legal values are \"off\", \"all\" and \"remain\"). Please verify the XML file arround the line number %d.", interface_name, XML_file_name, line_number);
            if (words_are_the_same(default_str, "off")) {
                fields_default_setting = 0;
                for (TiXmlNode *field_element_node = detailed_element->FirstChild(); field_element_node != NULL; field_element_node = field_element_node->NextSibling()) {
                    if (field_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                        continue;                    
                    TiXmlElement *field_element = field_element_node->ToElement();
                    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, words_are_the_same(field_element->Value(),"field"), "When setting the attribute \"fields\" for an coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", please use the keyword \"field\" for the name of a field. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, field_element->Row());
                    const char *field_name = get_XML_attribute(host_comp_id, -1, field_element, "name", XML_file_name, line_number, "the name of a field", "import interface configuration file", true);    
                    check_and_verify_name_format_of_string_for_XML(host_comp_id, field_name, "the field", XML_file_name, line_number);
                    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, fields_info->search_field_info(field_name) != NULL, "When setting the attribute \"fields\" for an coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", an illegal field name (\"%s\") is detected. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, field_name, field_element->Row());
                    for (i = 0; i < interface_fields_name.size(); i ++)
                        if (words_are_the_same(interface_fields_name[i], field_name))
                            break;
                    if (i < interface_fields_name.size())
                        fields_name.push_back(strdup(interface_fields_name[i]));
                    else EXECUTION_REPORT(REPORT_WARNING, host_comp_id, false, "When setting the attribute \"fields\" for the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the interface does not contain a field with the name of \"%s\"", interface_name, XML_file_name, field_name);
                }
                EXECUTION_REPORT(REPORT_WARNING, host_comp_id, fields_name.size() > 0, "When setting a coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", there are no fields specified. Please note the XML file arround the line number %d.", interface_name, XML_file_name, detailed_element->Row());
            }
            else if (words_are_the_same(default_str, "all")) {
                fields_default_setting = 1;
                for (i = 0; i < interface_fields_name.size(); i ++) {
                    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, fields_count[i] == 0, "When setting the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the configuration information of field \"%s\" has been set more than once. This is not allowed. Please note that the default value \"all\" means all fields. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, interface_fields_name[i], detailed_element->Row()); 
                    fields_count[i] ++;
                    fields_name.push_back(strdup(interface_fields_name[i]));                    
                }
            }
            else {
                fields_default_setting = 2;
                for (i = 0; i < interface_fields_name.size(); i ++) {
                    if (fields_count[i] != 0)
                        continue;
                    fields_count[i] ++;
                    fields_name.push_back(strdup(interface_fields_name[i]));                    
                }
            }
        }
        else if (words_are_the_same(detailed_element->Value(), "components")) {
            if (!is_XML_setting_on(host_comp_id, detailed_element, XML_file_name, "the status of \"components\"", "import interface configuration file"))
                continue;
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, components_element == NULL, "When setting the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the attribute of \"components\" has been set at least twice in a coupling connection specification. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, redirection_element->Row());
            components_element = detailed_element;
            const char *default_str = get_XML_attribute(host_comp_id, -1, detailed_element, "default", XML_file_name, line_number, "default setting for components", "import interface configuration file", true);
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, words_are_the_same(default_str, "off") || words_are_the_same(default_str, "all"), "When setting the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the value of \"default\" for the attribute of \"componets\" is wrong (legal values are \"off\" and \"all\"). Please verify the XML file arround the line number %d.", interface_name, XML_file_name, line_number);
            if (words_are_the_same(default_str, "off")) {
                components_default_setting = 0;
                for (TiXmlNode *component_element_node = detailed_element->FirstChild(); component_element_node != NULL; component_element_node = component_element_node->NextSibling()) {
                    if (component_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                        continue;
                    TiXmlElement *component_element = component_element_node->ToElement();
                    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, words_are_the_same(component_element->Value(),"component") || words_are_the_same(component_element->Value(),"interface") || words_are_the_same(component_element->Value(),"entry"), "When setting the attribute \"components\" for the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", please use the keyword \"component\", \"interface\" or \"entry\" for a specification of source component models or export interfaces (current wrong keyword is \"%s\"). Please verify the XML file arround the line number %d.", interface_name, XML_file_name, component_element->Value(), component_element->Row());
                    const char *full_name = get_XML_attribute(host_comp_id, 512, component_element, "comp_full_name", XML_file_name, line_number, "the full name of a component", "import interface configuration file", false);
//                    if (check_comp_existence)
//                        EXECUTION_REPORT(REPORT_ERROR, host_comp_id, comp_comm_group_mgt_mgr->search_global_node(full_name) != NULL, "When setting the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the full component name (\"%s\") is wrong. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, full_name, line_number);
                    const char *interface_name = component_element->Attribute("interface_name", &line_number);
                    if (full_name == NULL && interface_name == NULL)
                        continue;
                    if (full_name == NULL)
                        producer_info.first = strdup("\0");
                    else producer_info.first = strdup(full_name);
                    if (interface_name == NULL)
                        producer_info.second = strdup("\0");
                    else producer_info.second = strdup(interface_name);
                    for (i = 0; i < producers_info.size(); i ++)
                        if (words_are_the_same(producers_info[i].first, producer_info.first) && words_are_the_same(producers_info[i].second, producer_info.second))
                            break;
                    if (i == producers_info.size())
                        producers_info.push_back(producer_info);                
                }
            }
            else {
                components_default_setting = 1;
                const int *all_components_ids = comp_comm_group_mgt_mgr->get_all_components_ids();
                for (i = 1; i < all_components_ids[0]; i ++) {
                    producer_info.first = strdup(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(all_components_ids[i],false,"in Import_direction_setting()")->get_comp_full_name());
                    producer_info.second = strdup("\0");
                    producers_info.push_back(producer_info);
                }
            }
        }
        else if (words_are_the_same(detailed_element->Value(), "merge_setting")) {
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "When setting the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the attribute of \"merge_setting\" is not supported currently", interface_name, XML_file_name);
        }
        else EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "When setting the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", \"%s\" is not a legal attribute. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, detailed_element->Value(), detailed_element->Row());
    }        
    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, fields_element != NULL, "For a coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the information about fields is not set. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, redirection_element->Row());
    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, components_element != NULL, "For a coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the information about source component models is not set. Please verify the XML file arround the line number %d.", interface_name, XML_file_name, redirection_element->Row());
    for (i = 0; i < fields_name.size(); i ++)
        for (int j = 0; j < producers_info.size(); j ++)
            interface_configuration->add_field_src_component(host_comp_id, fields_name[i], producers_info[j]);
}


Import_direction_setting::~Import_direction_setting()
{
    clean_string_pair_vector(producers_info);
}


Import_interface_configuration::Import_interface_configuration(int host_comp_id, const char *comp_full_name, const char *interface_name, TiXmlElement *interface_element, const char *XML_file_name, Inout_interface_mgt *all_interfaces_mgr, bool check_comp_existence)
{
    int *fields_count, line_number;
    Inout_interface *interface_ptr = all_interfaces_mgr->get_interface(comp_full_name, interface_name);
    std::vector<std::pair<const char*, const char*> > producers_info;
    
    
    strcpy(this->interface_name, interface_name);
    interface_ptr->get_fields_name(&fields_name);
    fields_count = new int [fields_name.size()];
    for (int i = 0; i < fields_name.size(); i ++)
        fields_count[i] = 0;

    for (int i = 0; i < fields_name.size(); i ++)
        fields_src_producers_info.push_back(producers_info);

    for (TiXmlNode *redirection_element_node = interface_element->FirstChild(); redirection_element_node != NULL; redirection_element_node = redirection_element_node->NextSibling()) {
        if (redirection_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *redirection_element = redirection_element_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, interface_ptr->get_comp_id(), words_are_the_same(redirection_element->Value(),"import_connection"), "When setting the coupling connection configuration of the import interface \"%s\" in the XML file \"%s\", the XML element for specifying the coupling connection configuration should be named \"import_connection\" (current wrong name is \"%s\"). Please verify the XML file arround the line number %d.", interface_name, XML_file_name, redirection_element->Value(), redirection_element->Row());
        if (!is_XML_setting_on(interface_ptr->get_comp_id(), redirection_element, XML_file_name, "the status of some coupling connection configurations for an import interface", "import interface configuration file"))
            continue;
        import_directions.push_back(new Import_direction_setting(host_comp_id, this, comp_full_name, interface_name, redirection_element, XML_file_name, fields_name, fields_count, check_comp_existence));
    }

    delete [] fields_count;
}


void Import_interface_configuration::add_field_src_component(int comp_id, const char *field_name, std::pair<const char*, const char*> producer_info)
{
    int i;
    
    for (i = 0; i < fields_name.size(); i ++)
        if (words_are_the_same(field_name, fields_name[i]))
            break;
    EXECUTION_REPORT(REPORT_ERROR, comp_id, i < fields_name.size(), "Software error in Import_interface_configuration::add_field_src_component");    
    fields_src_producers_info[i].push_back(std::make_pair(strdup(producer_info.first), strdup(producer_info.second)));
}


void Import_interface_configuration::get_field_import_configuration(const char *field_name, std::vector<std::pair<const char*, const char*> > &producers_info)
{
    int i;


    for (i = 0; i < fields_name.size(); i ++)
        if (words_are_the_same(fields_name[i], field_name)) {
            for (int j = 0; j < fields_src_producers_info[i].size(); j ++)
                producers_info.push_back(std::make_pair(strdup(fields_src_producers_info[i][j].first), strdup(fields_src_producers_info[i][j].second)));
            break;
        }

    EXECUTION_REPORT(REPORT_ERROR, -1, i < fields_name.size(), "Software error in Import_interface_configuration::get_field_import_configuration");
}


Import_interface_configuration::~Import_interface_configuration()
{
    for (int i = 0; i < fields_src_producers_info.size(); i ++)
        clean_string_pair_vector(fields_src_producers_info[i]);
}


Component_import_interfaces_configuration::Component_import_interfaces_configuration(int host_comp_id, const char *comp_full_name, Inout_interface_mgt *interface_mgr, bool check_comp_existence)
{
    char XML_file_name[NAME_STR_SIZE];
    int line_number;


    strcpy(this->comp_full_name, comp_full_name);
    sprintf(XML_file_name, "%s/all/coupling_connections/%s.coupling_connections.xml", comp_comm_group_mgt_mgr->get_config_root_dir(), comp_full_name);
    TiXmlDocument *XML_file = open_XML_file_to_read(host_comp_id, XML_file_name, MPI_COMM_NULL, false);    
    if (XML_file == NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "As there is no import interface configuration file (the file name should be \"%s.coupling_connections.xml\") specified for the component \"%s\", the coupling procedures of the import/export interfaces of this component will be generated automatically", 
                         comp_full_name, comp_full_name);
        return;
    }
    TiXmlElement *root_XML_element = XML_file->FirstChildElement();
    TiXmlNode *root_XML_element_node = get_XML_first_child_of_unique_root(host_comp_id, XML_file_name, XML_file);
    for (; root_XML_element_node != NULL; root_XML_element_node = root_XML_element_node->NextSibling()) {    
        if (root_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        root_XML_element = root_XML_element_node->ToElement();
        if (words_are_the_same(root_XML_element->Value(),"local_import_interfaces"))
            break;
    }
    if (root_XML_element_node == NULL) {
        delete XML_file;
        return;
    }

    for (TiXmlNode *interface_XML_element_node = root_XML_element->FirstChild(); interface_XML_element_node != NULL; interface_XML_element_node = interface_XML_element_node->NextSibling()) {
        if (interface_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *interface_XML_element = interface_XML_element_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(interface_XML_element->Value(),"import_interface"), "Detect a wrong name (\"%s\") of an XML element in the XML configuration file \"%s\": the XML element for specifying the coupling connection configuration of an import interface should be named \"import_interface\". Please verify the XML file arround the line number %d.", interface_XML_element->Value(), XML_file_name, interface_XML_element->Row());
        const char *interface_name = get_XML_attribute(host_comp_id, CCPL_NAME_STR_LEN, interface_XML_element, "name", XML_file_name, line_number, "the \"name\" of an import interface", "import interface configuration file", true);
        if (!is_XML_setting_on(host_comp_id, interface_XML_element, XML_file_name, "the \"status\" of the coupling connection configurations for an import interface", "import interface configuration file"))
            continue;
        check_and_verify_name_format_of_string_for_XML(-1, interface_name, "the import interface", XML_file_name, line_number);
        Inout_interface *import_interface = interface_mgr->get_interface(comp_full_name, interface_name);
        if (import_interface == NULL) {
            EXECUTION_REPORT(REPORT_WARNING, -1, false, "The coupling connection configuration of the import interface named \"%s\" has been specified in the XML configuration file \"%s\", while the component \"%s\" does not register the corresponding import interface. So this coupling connection configuration information is negleted.\"", 
                             interface_name, XML_file_name, comp_full_name);
            continue;
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, import_interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT, "Detect the coupling connection configuration of a non-import interface \"%s\" of the component model \"%s\" in the XML configuration file \"%s\": please note that it is meaningless and not allowed to specify coupling connections of a non-export interface. Please verify the model code or the XML file arround the line number %d",
                         interface_name, comp_full_name, XML_file_name, interface_XML_element->Row());
        for (int i = 0; i < import_interfaces_configuration.size(); i ++)
            EXECUTION_REPORT(REPORT_ERROR, -1, !words_are_the_same(import_interfaces_configuration[i]->get_interface_name(), import_interface->get_interface_name()), "The coupling connection configuration of the import interface named \"%s\" has been set more than once in the XML file \"%s\", which is not allowed (only once for an interface). Please verify.", import_interface->get_interface_name(), XML_file_name);
        import_interfaces_configuration.push_back(new Import_interface_configuration(host_comp_id, comp_full_name, import_interface->get_interface_name(), interface_XML_element, XML_file_name, interface_mgr, check_comp_existence));
    }

    delete XML_file;
    EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Finish loading the configuration of import interfaces from the XML file %s", XML_file_name);
}


Component_import_interfaces_configuration::~Component_import_interfaces_configuration()
{
}


void Component_import_interfaces_configuration::get_interface_field_import_configuration(const char *interface_name, const char *field_name, std::vector<std::pair<const char*, const char*> > &producers_info)
{
    clean_string_pair_vector(producers_info);
    for (int i = 0; i < import_interfaces_configuration.size(); i ++)
        if (words_are_the_same(import_interfaces_configuration[i]->get_interface_name(), interface_name)) 
            import_interfaces_configuration[i]->get_field_import_configuration(field_name, producers_info);
}


Coupling_generator::~Coupling_generator()
{
    clear();
}


void Coupling_generator::clear()
{
    if (import_field_index_lookup_table != NULL)
        delete import_field_index_lookup_table;
    if (export_field_index_lookup_table != NULL)
        delete export_field_index_lookup_table;
    import_field_index_lookup_table = NULL;
    export_field_index_lookup_table = NULL;

    for (int i = 0; i < string_in_export_fields_dst_components.size(); i ++)
        delete [] string_in_export_fields_dst_components[i];
    string_in_export_fields_dst_components.clear();
    export_fields_dst_components.clear();

    all_coupling_connections.clear();
    all_IO_connections.clear();

    for (int i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        delete [] all_comp_fullnames_for_coupling_generation[i];
    all_comp_fullnames_for_coupling_generation.clear();
    individual_or_family_generation.clear();
}


void Coupling_generator::synchronize_latest_connection_id(MPI_Comm comm)
{
    int overall_latest_connection_id;
    MPI_Allreduce(&latest_connection_id, &overall_latest_connection_id, 1, MPI_INT, MPI_MAX, comm);
    latest_connection_id = overall_latest_connection_id;
}


void Coupling_generator::generate_coupling_procedures_common(int API_id, MPI_Comm comm, bool is_overall_generation, bool is_internal_generation, const char *annotation)
{
    char *temp_array_buffer = NULL, field_name[NAME_STR_SIZE];
    long current_array_buffer_size, max_array_buffer_size;
    int temp_int, num_pushed_comp_node = 0;    
    Coupling_connection *coupling_connection;
    std::pair<const char *,const char*> src_comp_interface;
    int current_proc_local_id;
    Comp_comm_group_mgt_node *local_comp_node = NULL, *temp_comp_node, *existing_comp_node;
    char API_label[NAME_STR_SIZE];
    

    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_rank(comm, &current_proc_local_id) == MPI_SUCCESS);    
    get_API_hint(-1, API_id, API_label);
    if (current_proc_local_id == 0)
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to generate coupling procedures commonly");
    synchronize_latest_connection_id(comm);
    inout_interface_mgr->get_all_unconnected_inout_interface_fields_info(all_comp_fullnames_for_coupling_generation, &temp_array_buffer, current_array_buffer_size, comm);
    bcast_array_in_one_comp(current_proc_local_id, &temp_array_buffer, current_array_buffer_size, comm);

    for (int i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        local_comp_node = comp_comm_group_mgt_mgr->search_global_node(all_comp_fullnames_for_coupling_generation[i]);
        if (local_comp_node != NULL && local_comp_node->get_current_proc_local_id() != -1)
            break;
        local_comp_node = NULL;
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, local_comp_node != NULL, "Software error in Coupling_generator::generate_coupling_procedures: wrong local_comp_node");

    for (int i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        temp_comp_node = comp_comm_group_mgt_mgr->load_comp_info_from_XML(local_comp_node->get_comp_id(), all_comp_fullnames_for_coupling_generation[i], comm);
        existing_comp_node = comp_comm_group_mgt_mgr->search_global_node(all_comp_fullnames_for_coupling_generation[i]);
        if (existing_comp_node == NULL) {
            comp_comm_group_mgt_mgr->push_comp_node(temp_comp_node);
            num_pushed_comp_node ++;
            existing_comp_node = temp_comp_node;
        }
        else delete temp_comp_node;
        if (existing_comp_node->get_current_proc_local_id() != -1)
            EXECUTION_REPORT(REPORT_ERROR, existing_comp_node->get_comp_id(), !components_time_mgrs->get_time_mgr(existing_comp_node->get_comp_id())->get_time_has_been_advanced(), "Error happens when calling the API \"%s\": fail to generate coupling procedures because the time of the component model \"%s\" has been advanced before this API call. Please check the model code related to the annotation \"%s\"", API_label, all_comp_fullnames_for_coupling_generation[i], annotation);
    }

    if (current_proc_local_id == 0) {
        Inout_interface_mgt *all_interfaces_mgr = new Inout_interface_mgt(temp_array_buffer, current_array_buffer_size);
        generate_interface_fields_source_dst(temp_array_buffer, current_array_buffer_size);
        for (int i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
            std::vector<Inout_interface*> import_interfaces_of_a_component;
            all_interfaces_mgr->get_all_import_interfaces_of_a_component(import_interfaces_of_a_component, all_comp_fullnames_for_coupling_generation[i]);
            Component_import_interfaces_configuration *comp_import_interfaces_config = new Component_import_interfaces_configuration(local_comp_node->get_comp_id(), all_comp_fullnames_for_coupling_generation[i], all_interfaces_mgr, is_overall_generation);
            for (int j = 0; j < import_interfaces_of_a_component.size(); j ++) {
                std::vector<const char*> import_fields_name;
                import_interfaces_of_a_component[j]->get_fields_name(&import_fields_name);
                for (int k = 0; k < import_fields_name.size(); k ++) {
                    std::vector<std::pair<const char*, const char*> > configuration_export_producer_info;
                    coupling_connection = new Coupling_connection(coupling_generator->apply_connection_id());
                    comp_import_interfaces_config->get_interface_field_import_configuration(import_interfaces_of_a_component[j]->get_interface_name(), import_fields_name[k], configuration_export_producer_info);
                    strcpy(coupling_connection->dst_comp_full_name, all_comp_fullnames_for_coupling_generation[i]);
                    strcpy(coupling_connection->dst_interface_name, import_interfaces_of_a_component[j]->get_interface_name());
                    coupling_connection->fields_name.push_back(strdup(import_fields_name[k]));                    
                    int field_index = export_field_index_lookup_table->search(import_fields_name[k],false);
                    if (field_index != 0) {
                        if (configuration_export_producer_info.size() == 0) {                        
                            for (int l = 0; l < export_fields_dst_components[field_index].size(); l ++) {
                                if (!is_internal_generation && words_are_the_same(export_fields_dst_components[field_index][l].first, all_comp_fullnames_for_coupling_generation[i]))
                                    continue;
                                src_comp_interface.first = strdup(export_fields_dst_components[field_index][l].first);
                                src_comp_interface.second = strdup(export_fields_dst_components[field_index][l].second);
                                coupling_connection->src_comp_interfaces.push_back(src_comp_interface);
                            }
                        }
                        else {
                            for (int l = 0; l < configuration_export_producer_info.size(); l ++) {
                                for (int m = 0; m < export_fields_dst_components[field_index].size(); m ++) {
                                    bool matching = false;
                                    if (strlen(configuration_export_producer_info[l].first) == 0)
                                        matching = words_are_the_same(configuration_export_producer_info[l].second, export_fields_dst_components[field_index][m].second);
                                    else if (words_are_the_same(configuration_export_producer_info[l].first, export_fields_dst_components[field_index][m].first))
                                        matching = strlen(configuration_export_producer_info[l].second) == 0 || words_are_the_same(configuration_export_producer_info[l].second, export_fields_dst_components[field_index][m].second);
                                    if (matching) {
                                        src_comp_interface.first = strdup(export_fields_dst_components[field_index][m].first);
                                        src_comp_interface.second = strdup(export_fields_dst_components[field_index][m].second);
                                        coupling_connection->src_comp_interfaces.push_back(src_comp_interface);
                                    }
                                }
                            }
                        }
                    }
                    char *report_string = NULL;
                    if (coupling_connection->src_comp_interfaces.size() > 0) {
                        report_string = new char [NAME_STR_SIZE*coupling_connection->src_comp_interfaces.size()];
                        report_string[0] = '\0';
                        for (int j = 0; j < coupling_connection->src_comp_interfaces.size(); j ++)
                            sprintf(report_string, "%s             %d) Component model is \"%s\", export interface is \"%s\"\n", report_string, j+1, coupling_connection->src_comp_interfaces[j].first, coupling_connection->src_comp_interfaces[j].second);                    
                    }
                    if (coupling_connection->src_comp_interfaces.size() == 1) {
                        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Field \"%s\" of the import interface \"%s\" in the component model \"%s\" have one source as follows. %s\n", coupling_connection->fields_name[0], coupling_connection->dst_interface_name, coupling_connection->dst_comp_full_name, report_string);
                        all_coupling_connections.push_back(coupling_connection);
                    }    
                    else if (coupling_connection->src_comp_interfaces.size() > 1)
                        EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when calling the API \"%s\" for coupling generation: Field \"%s\" of the import interface \"%s\" in the component model \"%s\" have more than one source as follows. Please verify.\n%s\n", API_label, coupling_connection->fields_name[0], coupling_connection->dst_interface_name, coupling_connection->dst_comp_full_name, report_string);
                    else delete coupling_connection;
                    if (report_string != NULL)
                        delete [] report_string;
                }
            }
            delete comp_import_interfaces_config;
        }
        
        for (int j, i = all_coupling_connections.size() - 1; i >= 0; i --) {
            for (j = 0; j < i; j ++)
                if (words_are_the_same(all_coupling_connections[i]->src_comp_interfaces[0].first, all_coupling_connections[j]->src_comp_interfaces[0].first) &&
                    words_are_the_same(all_coupling_connections[i]->src_comp_interfaces[0].second, all_coupling_connections[j]->src_comp_interfaces[0].second) &&
                    words_are_the_same(all_coupling_connections[i]->dst_comp_full_name, all_coupling_connections[j]->dst_comp_full_name) &&
                    words_are_the_same(all_coupling_connections[i]->dst_interface_name, all_coupling_connections[j]->dst_interface_name))
                    break;
            if (j < i) {
                EXECUTION_REPORT(REPORT_ERROR, -1, all_coupling_connections[i]->fields_name.size() == 1,  "software error in Coupling_generator::generate_coupling_procedures: %d", all_coupling_connections[i]->fields_name.size());
                all_coupling_connections[j]->fields_name.push_back(all_coupling_connections[i]->fields_name[0]);
                all_coupling_connections.erase(all_coupling_connections.begin()+i);
            }
        }

        delete all_interfaces_mgr;

        if (temp_array_buffer != NULL)
            delete [] temp_array_buffer;
        temp_array_buffer = NULL;
        current_array_buffer_size = 0;

        for (int i = all_coupling_connections.size() - 1; i >= 0; i --) {
            // all_coupling_connections[i]->src_comp_interfaces.size() is 1
            for (int j = all_coupling_connections[i]->src_comp_interfaces.size()-1; j >= 0; j --) {
                write_data_into_array_buffer(all_coupling_connections[i]->src_comp_interfaces[j].second, NAME_STR_SIZE, &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
                write_data_into_array_buffer(all_coupling_connections[i]->src_comp_interfaces[j].first, NAME_STR_SIZE, &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
            }
            temp_int = all_coupling_connections[i]->src_comp_interfaces.size();
            write_data_into_array_buffer(&temp_int, sizeof(int), &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
            for (int j = all_coupling_connections[i]->fields_name.size() - 1; j >= 0; j --)
                write_data_into_array_buffer(all_coupling_connections[i]->fields_name[j], NAME_STR_SIZE, &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
            temp_int = all_coupling_connections[i]->fields_name.size();
            write_data_into_array_buffer(&temp_int, sizeof(int), &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);            
            write_data_into_array_buffer(all_coupling_connections[i]->dst_interface_name, NAME_STR_SIZE, &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
            write_data_into_array_buffer(all_coupling_connections[i]->dst_comp_full_name, NAME_STR_SIZE, &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
            write_data_into_array_buffer(&(all_coupling_connections[i]->connection_id), sizeof(int), &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
        }
        temp_int = all_coupling_connections.size();
        write_data_into_array_buffer(&temp_int, sizeof(int), &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
    }
    
    bcast_array_in_one_comp(current_proc_local_id, &temp_array_buffer, current_array_buffer_size, comm);
    if (current_proc_local_id != 0) {
        int num_connections, num_fields, num_sources;
        long buffer_content_iter = current_array_buffer_size;
        read_data_from_array_buffer(&num_connections, sizeof(int), temp_array_buffer, buffer_content_iter, true);
        for (int i = 0; i < num_connections; i ++) {
            int connection_id;
            read_data_from_array_buffer(&connection_id, sizeof(int), temp_array_buffer, buffer_content_iter, true);
            coupling_connection = new Coupling_connection(connection_id);
            read_data_from_array_buffer(coupling_connection->dst_comp_full_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
            read_data_from_array_buffer(coupling_connection->dst_interface_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
            read_data_from_array_buffer(&num_fields, sizeof(int), temp_array_buffer, buffer_content_iter, true);
            for (int j = 0; j < num_fields; j ++) {
                read_data_from_array_buffer(field_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
                coupling_connection->fields_name.push_back(strdup(field_name));    
                EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "connection field name %s for interface %s", coupling_connection->dst_interface_name, field_name);
            }        
            read_data_from_array_buffer(&num_sources, sizeof(int), temp_array_buffer, buffer_content_iter, true);
            for (int j = 0; j < num_sources; j ++) {
                char tmp_str[NAME_STR_SIZE];
                read_data_from_array_buffer(tmp_str, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
                src_comp_interface.first = strdup(tmp_str);
                read_data_from_array_buffer(tmp_str, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
                src_comp_interface.second = strdup(tmp_str);
                coupling_connection->src_comp_interfaces.push_back(src_comp_interface);
            }
            all_coupling_connections.push_back(coupling_connection);
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, buffer_content_iter == 0, "Software error in Coupling_generator::generate_coupling_procedures: %d", buffer_content_iter);
    }

    delete [] temp_array_buffer;
    for (int i = 0; i < all_coupling_connections.size(); i ++) {
        all_coupling_connections[i]->generate_a_coupling_procedure(false);
    }
    for (int i = 0; i < num_pushed_comp_node; i ++)
        delete comp_comm_group_mgt_mgr->pop_comp_node();

    clear();

    //original_grid_mgr->delete_external_original_grids();
    
    if (current_proc_local_id == 0)
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish generating coupling procedure");
}



void Coupling_generator::generate_coupling_procedures_internal(int comp_id, bool family_generation, bool is_internal_generation, const char *annotation)
{
    char *temp_array_buffer = NULL;
    long current_array_buffer_size, max_array_buffer_size;


    if (family_generation) {
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, true, "in Coupling_generator::generate_coupling_procedures")->get_all_descendant_real_comp_fullnames(comp_id, all_comp_fullnames_for_coupling_generation, &temp_array_buffer, max_array_buffer_size, current_array_buffer_size);
        generate_coupling_procedures_common(API_ID_COUPLING_GEN_FAMILY, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Coupling_generator::generate_coupling_procedures"), (comp_id & TYPE_ID_SUFFIX_MASK)==0, is_internal_generation, annotation);
    }
    else {
        all_comp_fullnames_for_coupling_generation.push_back(strdup(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, true, "in Coupling_generator::generate_coupling_procedures")->get_full_name()));
        generate_coupling_procedures_common(API_ID_COUPLING_GEN_INDIVIDUAL, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Coupling_generator::generate_coupling_procedures"), (comp_id & TYPE_ID_SUFFIX_MASK)==0, is_internal_generation, annotation);
    }
}


void Coupling_generator::generate_IO_procedures()
{
    return;   ////// to be modified
/*
    const int *sorted_comp_ids = comp_comm_group_mgt_mgr->get_sorted_comp_ids();

    components_IO_output_procedures_mgr->add_all_components_IO_output_procedures();
    for (int i = 1; i < sorted_comp_ids[0]; i ++) {
        if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(sorted_comp_ids[i], "in Coupling_generator::generate_IO_procedures") == -1)
            continue;
        components_IO_output_procedures_mgr->get_component_IO_output_procedures(sorted_comp_ids[i])->generate_coupling_connection(all_IO_connections, all_coupling_connections.size());
    }

    for (int i = 0; i < all_IO_connections.size(); i ++) {
        all_IO_connections[i]->generate_a_coupling_procedure(false);
    }
*/
}


void Coupling_generator::generate_interface_fields_source_dst(const char *temp_array_buffer, int buffer_content_size)
{
    char comp_full_name[NAME_STR_SIZE], interface_name[NAME_STR_SIZE], field_name[NAME_STR_SIZE];
    std::vector<const char*> distinct_import_fields_name;
    std::vector<const char*> distinct_export_fields_name;


    import_field_index_lookup_table = new Dictionary<int>(1024);
    export_field_index_lookup_table = new Dictionary<int>(1024);

    long buffer_content_iter = buffer_content_size;
    int import_or_export, field_id_iter = 100, field_index, num_fields;
    while (buffer_content_iter > 0) {
        read_data_from_array_buffer(interface_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
        read_data_from_array_buffer(comp_full_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
        read_data_from_array_buffer(&import_or_export, sizeof(int), temp_array_buffer, buffer_content_iter, true);
        read_data_from_array_buffer(&num_fields, sizeof(int), temp_array_buffer, buffer_content_iter, true);
        for (int i = 0; i < num_fields; i ++) {
            read_data_from_array_buffer(field_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
            if (import_or_export == 0) {
                if (import_field_index_lookup_table->search(field_name, false) == 0) {
                    import_field_index_lookup_table->insert(field_name, field_id_iter++);
                    distinct_import_fields_name.push_back(strdup(field_name));
                }
                field_index = import_field_index_lookup_table->search(field_name, false);                }
            else {
                if (export_field_index_lookup_table->search(field_name, false) == 0) {
                    export_field_index_lookup_table->insert(field_name, field_id_iter++);
                    distinct_export_fields_name.push_back(strdup(field_name));
                }
                field_index = export_field_index_lookup_table->search(field_name, false);
                char *str1 = strdup(comp_full_name);
                char *str2 = strdup(interface_name);
                export_fields_dst_components[field_index].push_back(std::pair<const char*,const char*>(str1,str2));
                string_in_export_fields_dst_components.push_back(str1);
                string_in_export_fields_dst_components.push_back(str2);
            }
        }
    }
}


void Coupling_generator::transfer_interfaces_info_from_one_component_to_another(std::vector<Inout_interface*> &interfaces, Comp_comm_group_mgt_node *comp_node_src, Comp_comm_group_mgt_node *comp_node_dst)
{
    long buffer_max_size = 0, buffer_content_size = 0;
    char *temp_array_buffer = NULL;

    
    for (int i = 0; i < interfaces.size(); i ++)
        interfaces[i]->transform_interface_into_array(&temp_array_buffer, buffer_max_size, buffer_content_size);
    transfer_array_from_one_comp_to_another(comp_node_src->get_current_proc_local_id(), comp_node_src->get_root_proc_global_id(), comp_node_dst->get_current_proc_local_id(), comp_node_dst->get_root_proc_global_id(), comp_node_dst->get_comm_group(), &temp_array_buffer, buffer_content_size);

    if (comp_node_src->get_current_proc_local_id() < 0 && comp_node_dst->get_current_proc_local_id() >= 0) {
        EXECUTION_REPORT(REPORT_ERROR, -1, interfaces.size() == 0, "software error in Coupling_generator::transfer_interfaces_info_from_one_component_to_another");
        while (buffer_content_size > 0)
            interfaces.push_back(new Inout_interface(temp_array_buffer, buffer_content_size));
    }

    if (temp_array_buffer != NULL)
        delete [] temp_array_buffer;
}


void Coupling_generator::begin_external_coupling_generation()
{
    EXECUTION_REPORT(REPORT_ERROR, -1, all_comp_fullnames_for_coupling_generation.size() == 0 && individual_or_family_generation.size() == 0, "Software error in Coupling_generator::begin_external_coupling_generation");
}


void Coupling_generator::add_comp_for_external_coupling_generation(const char *comp_full_name, int individual_or_family, const char *annotation)
{
    int i;

    
    EXECUTION_REPORT(REPORT_ERROR, -1, strlen(comp_full_name) > 0, "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the full name of a given component model cannot be an empty string. Please verify the model code with the annotation \"%s\"", annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, individual_or_family == 1 || individual_or_family == 2, "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the value of the parameter \"individual_or_family\" must be 1 (individual) or 2 (family). Please verify the model code with the annotation \"%s\"", annotation);
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        if (words_are_the_same(comp_full_name, all_comp_fullnames_for_coupling_generation[i])) {
            if (individual_or_family_generation[i] == 1)
                individual_or_family_generation[i] = individual_or_family;
            break;
        }
    if (i == all_comp_fullnames_for_coupling_generation.size()) {
        all_comp_fullnames_for_coupling_generation.push_back(strdup(comp_full_name));
        individual_or_family_generation.push_back(individual_or_family);
    }
}


void Coupling_generator::sort_comp_full_names(std::vector<const char*> &comp_full_names, std::vector<int> *comp_index)
{
    std::vector<const char *> temp_comp_full_names;
    std::vector<int> temp_comp_index;
    int i, j, k, num_comps = comp_full_names.size();


    for (i = 0; i < num_comps; i ++) {
        k = 0;
        for (j = 1; j < comp_full_names.size(); j ++)
            if (strcmp(comp_full_names[k], comp_full_names[j]) > 0)
                k = j;
        temp_comp_full_names.push_back(comp_full_names[k]);
        comp_full_names.erase(comp_full_names.begin()+k);
        if (comp_index != NULL) {
            temp_comp_index.push_back((*comp_index)[k]);
            comp_index->erase(comp_index->begin()+k);
        }
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, num_comps == temp_comp_full_names.size() && comp_full_names.size() == 0, "Software error in Coupling_generator::do_external_coupling_generation");

    for (i = 0; i < temp_comp_full_names.size(); i ++) {
        comp_full_names.push_back(temp_comp_full_names[i]);
        if (comp_index != NULL)
            comp_index->push_back(temp_comp_index[i]);
    }    
}


void Coupling_generator::do_overall_coupling_generation(const char *local_root_comp_full_name, const char *annotation)
{
    int i; 

    
    comp_comm_group_mgt_mgr->get_root_comps_for_overall_coupling_generation(all_comp_fullnames_for_coupling_generation);
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        if (words_are_the_same(local_root_comp_full_name, all_comp_fullnames_for_coupling_generation[i]))
            break;
    if (i == all_comp_fullnames_for_coupling_generation.size()) {
        clear();
        return;
    }

    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        individual_or_family_generation.push_back(2);    
    do_external_coupling_generation(API_ID_COMP_MGT_END_COMP_REG, annotation);
}


void Coupling_generator::do_external_coupling_generation(int API_id, const char *annotation)
{
    int i, j, k, num_comps = all_comp_fullnames_for_coupling_generation.size(), *temp_int_array;
    std::vector<Comp_comm_group_mgt_node*> all_comp_nodes, loaded_comp_nodes;
    Comp_comm_group_mgt_node *temp_comp_node, *local_comp_node, *existing_comp_node;
    char *local_temp_array_buffer = NULL, *remote_temp_array_buffer = NULL;
    long local_current_array_buffer_size, local_max_array_buffer_size, remote_current_array_buffer_size, remote_max_array_buffer_size, temp_array_size, str_size;
    int current_connection_id, max_connection_id, remote_connection_id;
    MPI_Request request_send, request_recv;
    MPI_Status status;
    MPI_Comm union_comm = MPI_COMM_NULL;
    bool is_current_proc_in_union_comm = false;
    int current_proc_id_in_union_comm, size_union_comm;
    std::vector<int> proc_global_ids_in_union_comm, proc_global_ids_in_current_comm;
    std::vector<const char*> temp_all_comp_fullnames_for_coupling_generation;


    sort_comp_full_names(all_comp_fullnames_for_coupling_generation, &individual_or_family_generation);
    temp_int_array = new int [individual_or_family_generation.size()];
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        temp_int_array[i] = individual_or_family_generation[i];
    }

    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        if (comp_comm_group_mgt_mgr->search_global_node(all_comp_fullnames_for_coupling_generation[i]) != NULL)
            break;
    EXECUTION_REPORT(REPORT_ERROR, -1, i < all_comp_fullnames_for_coupling_generation.size(), "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the current process is not in any component model specified in the parameter \"comps_full_names\". Please verify the model code with the annotation \"%s\"", annotation);
        
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        local_comp_node = comp_comm_group_mgt_mgr->search_global_node(all_comp_fullnames_for_coupling_generation[i]);
        if (local_comp_node == NULL)
            continue;
        for (j = 0; j < all_comp_fullnames_for_coupling_generation.size(); j ++) {
            EXECUTION_REPORT_LOG(REPORT_LOG, local_comp_node->get_comp_id(), true, "The API call of \"CCPL_do_external_coupling_generation\" at the model code with the annotation \"%s\" try to access the component model \"%s\". Deadlock will be happen if the full name of the component model is wrong.", annotation, all_comp_fullnames_for_coupling_generation[j]);
            temp_comp_node = comp_comm_group_mgt_mgr->load_comp_info_from_XML(local_comp_node->get_comp_id(), all_comp_fullnames_for_coupling_generation[j], local_comp_node->get_comm_group());
            existing_comp_node = comp_comm_group_mgt_mgr->search_global_node(all_comp_fullnames_for_coupling_generation[j]);
            if (existing_comp_node != NULL) {
                delete temp_comp_node;
                temp_comp_node = existing_comp_node;
            }
            for (k = 0; k < all_comp_nodes.size(); k ++)
                if (words_are_the_same(all_comp_nodes[k]->get_full_name(), all_comp_fullnames_for_coupling_generation[j]))
                    break;
            if (k == all_comp_nodes.size()) {
                all_comp_nodes.push_back(temp_comp_node);
                if (existing_comp_node == NULL)
                    loaded_comp_nodes.push_back(temp_comp_node);
            }
            else if (existing_comp_node == NULL)
                delete temp_comp_node;
        }
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, all_comp_nodes.size() == all_comp_fullnames_for_coupling_generation.size(), "Software error in Coupling_generator::do_external_coupling_generation: wrong all_comp_nodes.size()");

    if (all_comp_fullnames_for_coupling_generation.size() == 1) {
        generate_coupling_procedures_internal(all_comp_nodes[0]->get_comp_id(), individual_or_family_generation[0] == 2, false, annotation);
        delete [] temp_int_array;
        return;
    }    

    for (int i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        dump_string(all_comp_fullnames_for_coupling_generation[i], -1, &local_temp_array_buffer, local_max_array_buffer_size, local_current_array_buffer_size);
    write_data_into_array_buffer(local_temp_array_buffer, local_current_array_buffer_size, &remote_temp_array_buffer, remote_max_array_buffer_size, remote_current_array_buffer_size);
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        if (all_comp_nodes[i]->get_current_proc_local_id() >= 0) {
            EXECUTION_REPORT_LOG(REPORT_LOG, all_comp_nodes[i]->get_comp_id(), true, "In the flowchart of executing the API \"CCPL_do_external_coupling_generation\" at the model code with the annotation \"%s\": before checking consistency of the parameters between all processes of this component model. Deadlock may happen if not all processes of this component model calls this API at the same time", annotation);
            check_API_parameter_data_array(all_comp_nodes[i]->get_comp_id(), API_id, all_comp_nodes[i]->get_comm_group(), "coupling generation", local_current_array_buffer_size, sizeof(char), local_temp_array_buffer, "comps_full_names", annotation);
            check_API_parameter_data_array(all_comp_nodes[i]->get_comp_id(), API_id, all_comp_nodes[i]->get_comm_group(), "coupling generation", individual_or_family_generation.size(), sizeof(int), (const char*) temp_int_array, "individual_or_familiy", annotation);
            coupling_generator->synchronize_latest_connection_id(all_comp_nodes[i]->get_comm_group());
            EXECUTION_REPORT_LOG(REPORT_LOG, all_comp_nodes[i]->get_comp_id(), true, "In the flowchart of executing the API \"CCPL_do_external_coupling_generation\" at the model code with the annotation \"%s\": after checking consistency of the parameters between all processes of this component model. Deadlock may happen if not all processes of this component model calls this API at the same time", annotation);            
        }
    max_connection_id = coupling_generator->get_latest_connection_id();
    for (i = 1; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        if (all_comp_nodes[0]->get_current_proc_local_id() == 0) {
            EXECUTION_REPORT_LOG(REPORT_LOG, all_comp_nodes[0]->get_comp_id(), true, "In the flowchart of executing the API \"CCPL_do_external_coupling_generation\" at the model code with the annotation \"%s\": check consistency of the parameters between the compoment models \"%s\" and \"%s\". Deadlock may happen if these two component models do not call this API at the same time", annotation, all_comp_fullnames_for_coupling_generation[0], all_comp_fullnames_for_coupling_generation[i]);
            MPI_Irecv(&remote_connection_id, 1, MPI_INT, all_comp_nodes[i]->get_local_proc_global_id(0), 10101, MPI_COMM_WORLD, &request_recv);
        }
        if (all_comp_nodes[i]->get_current_proc_local_id() == 0) {
            EXECUTION_REPORT_LOG(REPORT_LOG, all_comp_nodes[i]->get_comp_id(), true, "In the flowchart of executing the API \"CCPL_do_external_coupling_generation\" at the model code with the annotation \"%s\": check consistency of the parameters between the compoment models \"%s\" and \"%s\". Deadlock may happen if these two component models do not call this API at the same time", annotation, all_comp_fullnames_for_coupling_generation[i], all_comp_fullnames_for_coupling_generation[0]); 
            MPI_Isend(&max_connection_id, 1, MPI_INT, all_comp_nodes[0]->get_local_proc_global_id(0), 10101, MPI_COMM_WORLD, &request_send);
        }
        if (all_comp_nodes[0]->get_current_proc_local_id() == 0) {
            MPI_Wait(&request_recv, &status);
            max_connection_id = remote_connection_id > max_connection_id? remote_connection_id : max_connection_id;
        }
        if (all_comp_nodes[i]->get_current_proc_local_id() == 0)
            MPI_Wait(&request_send, &status);
    }
    for (i = 1; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        if (all_comp_nodes[0]->get_current_proc_local_id() == 0)
            MPI_Isend(&max_connection_id, 1, MPI_INT, all_comp_nodes[i]->get_local_proc_global_id(0), 10101, MPI_COMM_WORLD, &request_send);
        if (all_comp_nodes[i]->get_current_proc_local_id() == 0)
            MPI_Irecv(&max_connection_id, 1, MPI_INT, all_comp_nodes[0]->get_local_proc_global_id(0), 10101, MPI_COMM_WORLD, &request_recv);
        if (all_comp_nodes[0]->get_current_proc_local_id() == 0)
            MPI_Wait(&request_send, &status);
        if (all_comp_nodes[i]->get_current_proc_local_id() == 0)
            MPI_Wait(&request_recv, &status);
    }
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        if (all_comp_nodes[i]->get_current_proc_local_id() != -1)
            MPI_Bcast(&max_connection_id, 1, MPI_INT, 0, all_comp_nodes[i]->get_comm_group());
    coupling_generator->set_latest_connection_id(max_connection_id);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "The coupling connection id in Coupling_generator::do_external_coupling_generation is %d", coupling_generator->get_latest_connection_id());

    for (i = 1; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        if (all_comp_nodes[0]->get_current_proc_local_id() != -1 || all_comp_nodes[i]->get_current_proc_local_id() != -1) {
            transfer_array_from_one_comp_to_another(all_comp_nodes[0]->get_current_proc_local_id(), all_comp_nodes[0]->get_local_proc_global_id(0), all_comp_nodes[i]->get_current_proc_local_id(), all_comp_nodes[i]->get_local_proc_global_id(0), all_comp_nodes[i]->get_comm_group(), &remote_temp_array_buffer, remote_current_array_buffer_size);
            temp_array_size = individual_or_family_generation.size() * sizeof(int);
            transfer_array_from_one_comp_to_another(all_comp_nodes[0]->get_current_proc_local_id(), all_comp_nodes[0]->get_local_proc_global_id(0), all_comp_nodes[i]->get_current_proc_local_id(), all_comp_nodes[i]->get_local_proc_global_id(0), all_comp_nodes[i]->get_comm_group(), (char**)(&temp_int_array), temp_array_size);
            if (all_comp_nodes[i]->get_current_proc_local_id() != -1) {
                EXECUTION_REPORT(REPORT_ERROR, all_comp_nodes[i]->get_comp_id(), remote_current_array_buffer_size == local_current_array_buffer_size, "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the full names of component models specified through the input parameters are not consistency between the component models \"%s\" and \"%s\". Please check the model code with the annotation \"%s\"", all_comp_fullnames_for_coupling_generation[i], all_comp_fullnames_for_coupling_generation[0], annotation);
                for (j = 0; j < local_current_array_buffer_size; j ++)
                    if (local_temp_array_buffer[j] != remote_temp_array_buffer[j])
                        break;
                EXECUTION_REPORT(REPORT_ERROR, all_comp_nodes[i]->get_comp_id(), j == local_current_array_buffer_size, "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the full names of component models specified through the input parameters are not consistency between the component models \"%s\" and \"%s\". Please check the model code with the annotation \"%s\"", all_comp_fullnames_for_coupling_generation[i], all_comp_fullnames_for_coupling_generation[0], annotation);
                for (j = 0; j < individual_or_family_generation.size(); j ++)    
                    if (individual_or_family_generation[j] != temp_int_array[j])
                        break;
                EXECUTION_REPORT(REPORT_ERROR, all_comp_nodes[i]->get_comp_id(), j == individual_or_family_generation.size(), "ERROR happens when calling the API \"CCPL_do_external_coupling_generation\": the values of the parameter \"individual_or_family\" are not consistency between the component models \"%s\" and \"%s\". Please check the model code with the annotation \"%s\"", all_comp_fullnames_for_coupling_generation[i], all_comp_fullnames_for_coupling_generation[0], annotation);
            }
        }
    if (local_temp_array_buffer != NULL) {
        delete [] local_temp_array_buffer;
        local_temp_array_buffer = NULL;
    }
    if (remote_temp_array_buffer != NULL) {
        delete [] remote_temp_array_buffer;
        remote_temp_array_buffer = NULL;
    }    

    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        current_connection_id = coupling_generator->apply_connection_id();
        if (i == 0) {
            if (all_comp_nodes[i]->get_current_proc_local_id() != -1)
                union_comm = all_comp_nodes[i]->get_comm_group();
        }
        else if (is_current_proc_in_union_comm || all_comp_nodes[i]->get_current_proc_local_id() != -1) {
            int *temp_proc_global_ids = NULL, proc_global_id;
            long array_size;
            MPI_Comm_rank(MPI_COMM_WORLD, &proc_global_id);
            if (is_current_proc_in_union_comm) {
                MPI_Comm_rank(union_comm, &current_proc_id_in_union_comm);
                MPI_Comm_size(union_comm, &size_union_comm);
                array_size = size_union_comm * sizeof(int);
                temp_proc_global_ids = new int [size_union_comm];
                MPI_Allgather(&proc_global_id, 1, MPI_INT, temp_proc_global_ids, 1, MPI_INT, union_comm);
            }
            else current_proc_id_in_union_comm = -1;
            if (proc_global_id == all_comp_nodes[0]->get_local_proc_global_id(0) || all_comp_nodes[i]->get_current_proc_local_id() != -1) {
                if (proc_global_id == all_comp_nodes[0]->get_local_proc_global_id(0))
                    transfer_array_from_one_comp_to_another(0, all_comp_nodes[0]->get_local_proc_global_id(0), all_comp_nodes[i]->get_current_proc_local_id(), all_comp_nodes[i]->get_local_proc_global_id(0), all_comp_nodes[i]->get_comm_group(), (char**) (&temp_proc_global_ids), array_size);
                else transfer_array_from_one_comp_to_another(-1, all_comp_nodes[0]->get_local_proc_global_id(0), all_comp_nodes[i]->get_current_proc_local_id(), all_comp_nodes[i]->get_local_proc_global_id(0), all_comp_nodes[i]->get_comm_group(), (char**) (&temp_proc_global_ids), array_size);
            }
            proc_global_ids_in_union_comm.clear();
            for (int i = 0; i < array_size/sizeof(int); i ++)
                proc_global_ids_in_union_comm.push_back(temp_proc_global_ids[i]);
            if (temp_proc_global_ids != NULL)
                delete [] temp_proc_global_ids;
            proc_global_ids_in_current_comm.clear();
            for (j = 0; j < all_comp_nodes[i]->get_num_procs(); j ++)
                proc_global_ids_in_current_comm.push_back(all_comp_nodes[i]->get_local_proc_global_id(j));
            union_comm = create_union_comm_common(union_comm, all_comp_nodes[i]->get_comm_group(), current_proc_id_in_union_comm, all_comp_nodes[i]->get_current_proc_local_id(), proc_global_ids_in_union_comm, proc_global_ids_in_current_comm, current_connection_id, NULL, NULL);
        }
        if (all_comp_nodes[i]->get_current_proc_local_id() != -1)
            is_current_proc_in_union_comm = true;
    }

    MPI_Comm_size(union_comm, &size_union_comm);
    MPI_Comm_rank(union_comm, &current_proc_id_in_union_comm);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish generating the union comm (%d) for external coupling generation", size_union_comm);

    remote_max_array_buffer_size = 0; 
    remote_current_array_buffer_size = 0;
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++) {
        if (all_comp_nodes[i]->get_current_proc_local_id() != -1) {
            if (individual_or_family_generation[i] == 2)
                all_comp_nodes[i]->get_all_descendant_real_comp_fullnames(all_comp_nodes[i]->get_comp_id(), temp_all_comp_fullnames_for_coupling_generation, &local_temp_array_buffer, local_max_array_buffer_size, local_current_array_buffer_size);
            else temp_all_comp_fullnames_for_coupling_generation.push_back(strdup(all_comp_fullnames_for_coupling_generation[i]));
            if (all_comp_nodes[i]->get_current_proc_local_id() == 0)
                for (j = 0; j < temp_all_comp_fullnames_for_coupling_generation.size(); j ++)
                    dump_string(temp_all_comp_fullnames_for_coupling_generation[j], -1, &remote_temp_array_buffer, remote_max_array_buffer_size, remote_current_array_buffer_size);
            for (j = 0; j < temp_all_comp_fullnames_for_coupling_generation.size(); j ++)
                delete [] temp_all_comp_fullnames_for_coupling_generation[j];
            temp_all_comp_fullnames_for_coupling_generation.clear();
        }
    }
    local_temp_array_buffer = NULL;
    local_current_array_buffer_size = 0;
    local_max_array_buffer_size = 0;
    gather_array_in_one_comp(size_union_comm, current_proc_id_in_union_comm, remote_temp_array_buffer, remote_current_array_buffer_size, sizeof(char), NULL, (void**)(&local_temp_array_buffer), local_current_array_buffer_size, union_comm);
    bcast_array_in_one_comp(current_proc_id_in_union_comm, &local_temp_array_buffer, local_current_array_buffer_size, union_comm);
    for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
        delete [] all_comp_fullnames_for_coupling_generation[i];
    all_comp_fullnames_for_coupling_generation.clear();
    individual_or_family_generation.clear();
    while (local_current_array_buffer_size > 0) {
        const char *full_name = load_string(NULL, str_size, -1, local_temp_array_buffer, local_current_array_buffer_size, NULL);
        add_comp_for_external_coupling_generation(full_name, 1, annotation);
    }    
    EXECUTION_REPORT(REPORT_ERROR, -1, local_current_array_buffer_size == 0, "Software error in Coupling_generator::do_external_coupling_generation");
    if (local_temp_array_buffer != NULL)
        delete [] local_temp_array_buffer;
    if (remote_temp_array_buffer != NULL)
        delete [] remote_temp_array_buffer;

    if (current_proc_id_in_union_comm == 0) {        
        for (int i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
            EXECUTION_REPORT_LOG(REPORT_ERROR, -1, true, "comp for external generation %s at proc %d\n", all_comp_fullnames_for_coupling_generation[i], current_proc_id_in_union_comm);
    }
    generate_coupling_procedures_common(API_id, union_comm, false, false, annotation);

    clear();
    delete [] temp_int_array;
}


void Coupling_generator::load_comps_full_names_from_config_file(int comp_id, const char *keyword, int size_comps_full_names, int size_individual_or_family, int *num_comps, const char *annotation)
{
    char XML_file_name[NAME_STR_SIZE];
    const char *current_comp_full_name = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, true, "in Coupling_generator::load_comps_full_names_from_config_file")->get_full_name();
    const char *temp_full_name;
    const char *XML_keyword;
    int line_number, i;


    *num_comps = 0;
    
    sprintf(XML_file_name, "%s/all/coupling_connections/%s.coupling_connections.xml", comp_comm_group_mgt_mgr->get_config_root_dir(), current_comp_full_name);
    TiXmlDocument *XML_file = open_XML_file_to_read(comp_id, XML_file_name, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Coupling_generator::load_comps_full_names_from_config_file"), false);
    if (XML_file == NULL)
        return;
    
    TiXmlElement *root_XML_element, *XML_element, *detailed_XML_element;
    TiXmlNode *root_XML_element_node = get_XML_first_child_of_unique_root(comp_id, XML_file_name, XML_file), *XML_element_node = NULL, *detailed_XML_element_node = NULL;
    for (; root_XML_element_node != NULL; root_XML_element_node = root_XML_element_node->NextSibling()) {
        if (root_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        root_XML_element = root_XML_element_node->ToElement();
        if (words_are_the_same(root_XML_element->Value(),"component_full_names_sets"))
            break;
    }
    if (root_XML_element_node == NULL)
        return;

    for (XML_element_node = root_XML_element->FirstChild(); XML_element_node != NULL; XML_element_node = XML_element_node->NextSibling()) {
        XML_element = XML_element_node->ToElement();
        if (!is_XML_setting_on(comp_id, XML_element, XML_file_name, "the status for the full names of a set of component models", "the coupling connection configuration file"))
            continue;
        XML_keyword = get_XML_attribute(comp_id, -1, XML_element, "keyword", XML_file_name, line_number, "the keyword of the set of component models' full names", "the coupling connection configuration file", true);
        if (words_are_the_same(keyword, XML_keyword))
            break;
    }    
    if (XML_element_node == NULL)
        return;

    for (detailed_XML_element_node = XML_element->FirstChild(); detailed_XML_element_node != NULL; detailed_XML_element_node = detailed_XML_element_node->NextSibling()) {
        detailed_XML_element = detailed_XML_element_node->ToElement();
        temp_full_name = get_XML_attribute(comp_id, -1, detailed_XML_element, "comp_full_name", XML_file_name, line_number, "the full name of a component model", "the coupling connection configuration file", true);
        if (temp_full_name != NULL) {
            for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
                if (words_are_the_same(temp_full_name, all_comp_fullnames_for_coupling_generation[i]))
                    break;            
            if (i == all_comp_fullnames_for_coupling_generation.size())
                all_comp_fullnames_for_coupling_generation.push_back(strdup(temp_full_name));
            else EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "ERROR happens when calling the API \"CCPL_get_configurable_comps_full_names\": when loading the set of component models corresponding to the keyword \"%s\" from the XML file \"%s\", there are more than one entry of the component model \"%s\". Please check the XML file around the line number %d", keyword, XML_file_name, temp_full_name, line_number);
            const char *individual_or_family = detailed_XML_element->Attribute("individual_or_family", &line_number);
            int value_individual_or_family = 1;
            if (individual_or_family != NULL) {
                EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(individual_or_family, "individual") || words_are_the_same(individual_or_family, "family"), "ERROR happens when calling the API \"CCPL_get_configurable_comps_full_names\": the value of \"individual_or_family\" in the XML file \"%s\" must be \"individual\" or \"family\" (currently is \"%s\"). Please verify the XML file at line %d", XML_file_name, individual_or_family, line_number);
                if (words_are_the_same(individual_or_family, "family"))
                    value_individual_or_family = 2;
            }
            individual_or_family_generation.push_back(value_individual_or_family);
        }
    }

    if (all_comp_fullnames_for_coupling_generation.size() > 0) {
        for (i = 0; i < all_comp_fullnames_for_coupling_generation.size(); i ++)
            if (words_are_the_same(current_comp_full_name, all_comp_fullnames_for_coupling_generation[i]))
                break;
        if (i == all_comp_fullnames_for_coupling_generation.size())
            all_comp_fullnames_for_coupling_generation.push_back(strdup(current_comp_full_name));
    }
    
    *num_comps = all_comp_fullnames_for_coupling_generation.size();

    EXECUTION_REPORT(REPORT_ERROR, comp_id, size_comps_full_names >= *num_comps, "ERROR happens when calling the API \"CCPL_get_configurable_comps_full_names\" with the keyword \"%s\": the array size (currently is %d) of the input parameter \"comps_full_names\" is smaller than the number of component models (currently is %d) specified in the XML file \"%s\". Please verify the model code with the annotation \"%s\" or the configuration file.", keyword, size_comps_full_names, XML_file_name, *num_comps, annotation);
    if (size_individual_or_family != -1)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, size_individual_or_family >= *num_comps, "ERROR happens when calling the API \"CCPL_get_configurable_comps_full_names\" with the keyword \"%s\": the array size (currently is %d) of the input parameter \"individual_or_family\" is smaller than the number of component models (currently is %d) specified in the XML file \"%s\". Please verify the model code with the annotation \"%s\" or the configuration file.", keyword, size_comps_full_names, XML_file_name, *num_comps, annotation);
    
    delete XML_file;
}


void Coupling_generator::get_one_comp_full_name(int comp_id, const char *keyword, int str_size, int index, char *comp_full_name, int *local_individual_or_family, const char *annotation)
{
    char XML_file_name[NAME_STR_SIZE];


    sprintf(XML_file_name, "%s/all/coupling_connections/%s.coupling_connections.xml", comp_comm_group_mgt_mgr->get_config_root_dir(), comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, true, "in Coupling_generator::load_comps_full_names_from_config_file")->get_full_name());
    EXECUTION_REPORT(REPORT_ERROR, comp_id, str_size >= strlen(all_comp_fullnames_for_coupling_generation[index]), "Error happens when calling the API \"CCPL_get_configurable_comps_full_names\" with the keyword \"%s\": the string length (currently is %d) of the input parameter \"comps_full_names\" is smaller than the length (currently is %d) of the full name of a component model (\"%s\") that is loaded from the XML file \"%s\". Please verify the model code with the annotation \"%s\".", keyword, str_size, strlen(all_comp_fullnames_for_coupling_generation[index]), all_comp_fullnames_for_coupling_generation[index], XML_file_name, annotation);
    strncpy(comp_full_name, all_comp_fullnames_for_coupling_generation[index], strlen(all_comp_fullnames_for_coupling_generation[index]));
    *local_individual_or_family = individual_or_family_generation[index];
    for (int i = strlen(all_comp_fullnames_for_coupling_generation[index]); i < str_size; i ++)
        comp_full_name[i] = ' ';
}

