/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu and then
  *  modified by Dr. Cheng Zhang and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn or
  *  Dr. Cheng Zhang via zhangc-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "routing_info_mgt.h"
#include "global_data.h"
#include "cor_global_data.h"
#include "CCPL_api_mgt.h"
#include <stdio.h>
#include <string.h>

Routing_info *Routing_info_mgt::search_or_add_router(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
    Routing_info *router;

    router = search_router(src_comp_id, dst_comp_id, src_decomp_name, dst_decomp_name);

    if (router != NULL)
        return router;

    router = new Routing_info(src_comp_id, dst_comp_id, src_decomp_name, dst_decomp_name);
    routers.push_back(router);

    return router;
}


Routing_info_mgt::~Routing_info_mgt()
{
    for (int i = 0; i < routers.size(); i ++)
        delete routers[i];
}


Routing_info *Routing_info_mgt::search_router(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
    for (int i = 0; i < routers.size(); i ++)
        if (routers[i]->match_router(src_comp_id, dst_comp_id, src_decomp_name, dst_decomp_name))
            return routers[i];

    return NULL;
}


Routing_info::Routing_info(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
    src_decomp_info = decomps_info_mgr->search_decomp_info(src_decomp_name, src_comp_id);
    dst_decomp_info = decomps_info_mgr->search_decomp_info(dst_decomp_name, dst_comp_id);
    this->src_comp_id = src_comp_id;
    this->dst_comp_id = dst_comp_id;
    src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_comp_id);
    dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_comp_id);
    strcpy(this->src_comp_full_name, src_comp_node->get_comp_full_name());
    strcpy(this->index_dst_comp_full_name, dst_comp_node->get_comp_full_name());
    if (dst_decomp_info != NULL)
        dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_decomp_info->get_host_comp_id());
    strcpy(this->true_dst_comp_full_name, dst_comp_node->get_comp_full_name());
    src_comp_node_id = src_comp_node->get_comp_id();
    dst_comp_node_id = dst_comp_node->get_comp_id();
    strcpy(this->src_decomp_name, src_decomp_name);
    strcpy(this->dst_decomp_name, dst_decomp_name);
    src_decomp_size = 0;
    dst_decomp_size = 0;
    current_proc_id_src_comp = src_comp_node->get_current_proc_local_id();
    current_proc_id_dst_comp = dst_comp_node->get_current_proc_local_id();

    if (current_proc_id_src_comp != 0)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_id, true, "Start to generate router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);
    if (current_proc_id_dst_comp != 0)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_id, true, "Start to generate router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);

    if (words_are_the_same(src_decomp_name, "NULL")) {
        EXECUTION_REPORT(REPORT_ERROR,-1, words_are_the_same(dst_decomp_name, "NULL"), "for router of scalar variables, the local and remote decompositions must be \"NULL\"\n");
        num_dimensions = 0;
        if (current_proc_id_src_comp != -1) 
            src_decomp_size = 1;
        if (current_proc_id_dst_comp != -1) 
            dst_decomp_size = 1;
    }
    else {
        num_dimensions = 2;
        build_2D_router();
        if (current_proc_id_src_comp != -1) 
            src_decomp_size = src_decomp_info->get_num_local_cells();
        if (current_proc_id_dst_comp != -1) 
            dst_decomp_size = dst_decomp_info->get_num_local_cells();
    }

    if (current_proc_id_src_comp != 0)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_id, true, "Finish generating router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);
    if (current_proc_id_dst_comp != 0)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_id, true, "Finish generating router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);
}


Routing_info::~Routing_info()
{
    for (int i = 0; i < recv_from_remote_procs_routing_info.size(); i ++) {
        if (recv_from_remote_procs_routing_info[i]->num_elements_transferred > 0) {
            delete [] recv_from_remote_procs_routing_info[i]->local_indx_segment_starts;
            delete [] recv_from_remote_procs_routing_info[i]->local_indx_segment_lengths;
        }
    }
    for (int i = 0; i < send_to_remote_procs_routing_info.size(); i ++) {
        if (send_to_remote_procs_routing_info[i]->num_elements_transferred > 0) {
            delete [] send_to_remote_procs_routing_info[i]->local_indx_segment_starts;
            delete [] send_to_remote_procs_routing_info[i]->local_indx_segment_lengths;
        }
    }
}


bool Routing_info::match_router(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
    return (words_are_the_same(src_comp_full_name, comp_comm_group_mgt_mgr->search_global_node(src_comp_id)->get_full_name()) && words_are_the_same(index_dst_comp_full_name, comp_comm_group_mgt_mgr->search_global_node(dst_comp_id)->get_full_name()) &&
            words_are_the_same(this->src_decomp_name, src_decomp_name) && words_are_the_same(this->dst_decomp_name, dst_decomp_name));
}


void Routing_info::build_2D_router()
{
    int num_src_procs = src_comp_node->get_num_procs();
    int *num_cells_each_src_proc = new int [num_src_procs];
    int num_dst_procs = dst_comp_node->get_num_procs();
    int * num_cells_each_dst_proc = new int [num_dst_procs];
    int num_local_src_cells, num_local_dst_cells, *num_global_src_cells = new int [1], *num_global_dst_cells = new int [1];
    int *cells_indx_each_src_proc = NULL;
    int *cells_indx_each_dst_proc = NULL;
    int src_comp_root_proc_global_id = src_comp_node->get_root_proc_global_id();
    int dst_comp_root_proc_global_id = dst_comp_node->get_root_proc_global_id();
    Routing_info_with_one_process *routing_info;
    long total_src_cells, total_dst_cells;


    if (current_proc_id_src_comp != -1) {
        EXECUTION_REPORT(REPORT_ERROR, -1, src_decomp_info != NULL, "Software error in Routing_info::build_2D_router: NULL src decomp info");
        num_local_src_cells = src_decomp_info->get_num_local_cells();
        *num_global_src_cells = src_decomp_info->get_num_global_cells();
        gather_array_in_one_comp(num_src_procs, current_proc_id_src_comp, (void*)src_decomp_info->get_local_cell_global_indx(), num_local_src_cells, 
                                 sizeof(int), num_cells_each_src_proc, (void**)(&cells_indx_each_src_proc), total_src_cells, src_comp_node->get_comm_group());
    }
    if (current_proc_id_dst_comp != -1) {
        EXECUTION_REPORT(REPORT_ERROR, -1, dst_decomp_info != NULL, "Software error in Routing_info::build_2D_router: NULL dst decomp info");
        num_local_dst_cells = dst_decomp_info->get_num_local_cells();
        *num_global_dst_cells = dst_decomp_info->get_num_global_cells();
        gather_array_in_one_comp(num_dst_procs, current_proc_id_dst_comp, (void*)dst_decomp_info->get_local_cell_global_indx(), num_local_dst_cells, 
                                 sizeof(int), num_cells_each_dst_proc, (void**)(&cells_indx_each_dst_proc), total_dst_cells, dst_comp_node->get_comm_group());
    }

    long temp_size = num_src_procs*sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), (char**)(&num_cells_each_src_proc), temp_size);
    temp_size = num_dst_procs*sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), (char**)(&num_cells_each_dst_proc), temp_size);
    temp_size = sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), (char**)(&num_global_src_cells), temp_size);    
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT(REPORT_ERROR, -1, *num_global_src_cells == *num_global_dst_cells, "Software error in Routing_info::build_2D_router: different global decomp grid size: %d vs %d", *num_global_src_cells, *num_global_dst_cells);
    total_src_cells = 0;
    for (int i = 0; i < num_src_procs; i ++) 
        total_src_cells += num_cells_each_src_proc[i] * sizeof(int);
    total_dst_cells = 0;
    for (int i = 0; i < num_dst_procs; i ++) 
        total_dst_cells += num_cells_each_dst_proc[i] * sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), (char**)(&cells_indx_each_src_proc), total_src_cells);
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), (char**)(&cells_indx_each_dst_proc), total_dst_cells);
    
    if (current_proc_id_src_comp != -1) {
        int tmp_displs = 0;
        for (int i = 0; i < num_dst_procs; i ++) {
            routing_info = compute_routing_info_between_decomps(num_local_src_cells, src_decomp_info->get_local_cell_global_indx(), num_cells_each_dst_proc[i], cells_indx_each_dst_proc+tmp_displs, 
                                                                src_decomp_info->get_num_global_cells(), comp_comm_group_mgt_mgr->get_current_proc_global_id(), dst_comp_node->get_local_proc_global_id(i), true);
            tmp_displs += num_cells_each_dst_proc[i];
            send_to_remote_procs_routing_info.push_back(routing_info);
        }
    }

    if (current_proc_id_dst_comp != -1) {
        int tmp_displs = 0;
        for (int i = 0; i < num_src_procs; i ++) {
            routing_info = compute_routing_info_between_decomps(num_local_dst_cells, dst_decomp_info->get_local_cell_global_indx(), num_cells_each_src_proc[i], cells_indx_each_src_proc+tmp_displs, 
                                                                dst_decomp_info->get_num_global_cells(), comp_comm_group_mgt_mgr->get_current_proc_global_id(), src_comp_node->get_local_proc_global_id(i), false);
            tmp_displs += num_cells_each_src_proc[i];
            recv_from_remote_procs_routing_info.push_back(routing_info);
        }
    }
    
    if (cells_indx_each_src_proc != NULL) 
        delete [] cells_indx_each_src_proc;
    if (cells_indx_each_dst_proc != NULL) 
        delete [] cells_indx_each_dst_proc;
    delete [] num_cells_each_src_proc; 
    delete [] num_cells_each_dst_proc;
    delete [] num_global_src_cells;
    delete [] num_global_dst_cells;
}


Routing_info_with_one_process *Routing_info::compute_routing_info_between_decomps(int num_local_cells_local, const int *local_cells_global_indexes_local, 
                                                  int num_local_cells_remote, const int *local_cells_global_indexes_remote, 
                                                  int num_global_cells, int local_proc_id, int remote_proc_id, bool is_src)
{
    Routing_info_with_one_process *routing_info;
    const int *reference_cell_indx;
    int *logical_indx_lookup_table_local, *logical_indx_lookup_table_remote; 
    int num_reference_cells;
    int last_local_logical_indx;
    int j;


    routing_info = new Routing_info_with_one_process;
    routing_info->num_elements_transferred = 0;
    routing_info->num_local_indx_segments = 0;
    routing_info->remote_proc_global_id = remote_proc_id;

    if (num_local_cells_local == 0)
        return routing_info;
    
    /* Determine the reference cell index table according to the table size */
    if (is_src) {
        reference_cell_indx = local_cells_global_indexes_remote;
        num_reference_cells = num_local_cells_remote;  
    }
    else {
        reference_cell_indx = local_cells_global_indexes_local;
        num_reference_cells = num_local_cells_local; 
    }

    logical_indx_lookup_table_remote = new int [num_global_cells];
    logical_indx_lookup_table_local = new int [num_global_cells];
    for (j = 0; j < num_global_cells; j ++) {
        logical_indx_lookup_table_local[j] = -1;
        logical_indx_lookup_table_remote[j] = -1;
    }
    
    for (j = 0; j < num_local_cells_local; j ++)
        if (local_cells_global_indexes_local[j] >= 0)
            if (local_cells_global_indexes_local[j] != CCPL_NULL_INT)
                logical_indx_lookup_table_local[local_cells_global_indexes_local[j]] = j;
    for (j = 0; j < num_local_cells_remote; j ++)
        if (local_cells_global_indexes_remote[j] >= 0)
            if (local_cells_global_indexes_remote[j] != CCPL_NULL_INT)
                logical_indx_lookup_table_remote[local_cells_global_indexes_remote[j]] = j;

    /* Compute the number of common cells and the number of segments of common cells */
    last_local_logical_indx = -100;
    for (j = 0; j < num_reference_cells; j ++) 
        if (reference_cell_indx[j] != CCPL_NULL_INT && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
			if (reference_cell_indx == local_cells_global_indexes_local) {
	            if (last_local_logical_indx + 1 != j) 
	                routing_info->num_local_indx_segments ++;
	            last_local_logical_indx = j;
	            routing_info->num_elements_transferred ++;				
			}
			else {
	            if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]]) 
	                routing_info->num_local_indx_segments ++;
	            last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
	            routing_info->num_elements_transferred ++;
			}
        }

    /* Compute the info of segments when there are common cells */
    last_local_logical_indx = -100;
    if (routing_info->num_elements_transferred > 0) {
        routing_info->local_indx_segment_starts = new int [routing_info->num_local_indx_segments];
        routing_info->local_indx_segment_lengths = new int [routing_info->num_local_indx_segments];
        routing_info->num_local_indx_segments = 0;
        for (j = 0; j < num_reference_cells; j ++) 
            if (reference_cell_indx[j] != CCPL_NULL_INT && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
				if (reference_cell_indx == local_cells_global_indexes_local) {
	                if (last_local_logical_indx + 1 != j) {
	                    routing_info->local_indx_segment_starts[routing_info->num_local_indx_segments] = j;
	                    routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments] = 1;
	                    routing_info->num_local_indx_segments ++;
	                }
	                else routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments - 1] ++;
	                last_local_logical_indx = j;
				}
				else {
	                if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]]) {
	                    routing_info->local_indx_segment_starts[routing_info->num_local_indx_segments] = logical_indx_lookup_table_local[reference_cell_indx[j]];
	                    routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments] = 1;
	                    routing_info->num_local_indx_segments ++;
	                }
	                else routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments - 1] ++;
	                last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
				}
            }
    }

    delete [] logical_indx_lookup_table_remote;
    delete [] logical_indx_lookup_table_local;

    return routing_info;
}


Routing_info_with_one_process *Routing_info::get_routing_info(bool is_send, int i)
{
    if (is_send) {
        EXECUTION_REPORT(REPORT_ERROR, -1, i >= 0 && i < send_to_remote_procs_routing_info.size(), "Software error in Routing_info::get_num_elements_transferred_with_remote_proc: wrong i at sender");
        return send_to_remote_procs_routing_info[i];
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, i >= 0 && i < recv_from_remote_procs_routing_info.size(), "Software error in Routing_info::get_num_elements_transferred_with_remote_proc: wrong i at receiver");
        return recv_from_remote_procs_routing_info[i];
    }
}


Comp_comm_group_mgt_node *Routing_info::get_src_comp_node() 
{ 
    return comp_comm_group_mgt_mgr->search_global_node(src_comp_full_name); 
}


Comp_comm_group_mgt_node *Routing_info::get_dst_comp_node() 
{ 
    return comp_comm_group_mgt_mgr->search_global_node(true_dst_comp_full_name); 
}


