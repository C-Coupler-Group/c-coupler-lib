/***************************************************************
 *  Copyright (c) 2017, Tsinghua University.
 *  This is a source file of C-Coupler.
 *  This file was initially finished by Dr. Cheng Zhang and then
 *  modified by Dr. Cheng Zhang and Dr. Li Liu. 
 *  If you have any problem, 
 *  please contact Dr. Cheng Zhang via zhangc-cess@tsinghua.edu.cn
 *  or Dr. Li Liu via liuli-cess@tsinghua.edu.cn
 ***************************************************************/



#include "runtime_trans_algorithm.h"
#include "global_data.h"
#include <string.h>
#include <unistd.h>


template <class T> void Runtime_trans_algorithm::pack_segment_data(T *mpi_buf, T *field_data_buf, int segment_start, int segment_size, int field_2D_size, int num_lev, bool is_V1D_sub_grid_after_H2D_sub_grid)
{
    int i, j, offset;

    if (is_V1D_sub_grid_after_H2D_sub_grid) {
        for (i = segment_start, offset = 0; i < segment_size+segment_start; i ++)
            for (j = 0; j < num_lev; j ++)
                mpi_buf[offset++] = field_data_buf[i+j*field_2D_size];
    }
    else {
        for (i = segment_start, offset = 0; i < segment_size+segment_start; i ++)
            for (j = 0; j < num_lev; j ++)
                mpi_buf[offset++] = field_data_buf[i*num_lev+j];        
    }
}


template <class T> void Runtime_trans_algorithm::unpack_segment_data(T *mpi_buf, T *field_data_buf, int segment_start, int segment_size, int field_2D_size, int num_lev, bool is_V1D_sub_grid_after_H2D_sub_grid)
{
    int i, j, offset;

    if (is_V1D_sub_grid_after_H2D_sub_grid) {
        for (i = segment_start, offset = 0; i < segment_size+segment_start; i ++)
            for (j = 0; j < num_lev; j ++)
                field_data_buf[i+j*field_2D_size] = mpi_buf[offset++];
    }
    else {
        for (i = segment_start, offset = 0; i < segment_size+segment_start; i ++)
            for (j = 0; j < num_lev; j ++)
                field_data_buf[i*num_lev+j] = mpi_buf[offset++];        
    }
}


Runtime_trans_algorithm::Runtime_trans_algorithm(bool send_or_receive, int num_transfered_fields, Field_mem_info ** fields_mem, Routing_info ** routers, MPI_Comm comm, int * ranks, int connection_id)
{
    bool only_have_no_decomp_data = true;


    this->send_or_receive = send_or_receive;
    this->num_transfered_fields = num_transfered_fields;
    this->comm_tag = connection_id;
    EXECUTION_REPORT(REPORT_ERROR,-1, num_transfered_fields > 0, "Software error: Runtime_trans_algorithm does not have transfer fields");

    union_comm = comm;
    MPI_Comm_rank(union_comm, &current_proc_id_union_comm);

    this->fields_mem = new Field_mem_info *[num_transfered_fields];
    fields_data_buffers = new void *[num_transfered_fields];
    fields_routers = new Routing_info *[num_transfered_fields];
    last_history_receive_buffer_index = -1;
    last_field_remote_recv_count = -1;
    current_field_local_recv_count = 1;
    last_receive_sender_time = -1;

    for (int i = 0; i < num_transfered_fields; i ++) {
        this->fields_mem[i] = fields_mem[i];
        fields_data_buffers[i] = fields_mem[i]->get_data_buf();
        fields_routers[i] = routers[i];
    }

    if (send_or_receive) {
        local_comp_node = fields_routers[0]->get_src_comp_node();
        remote_comp_node = fields_routers[0]->get_dst_comp_node();
    }
    else {
        local_comp_node = fields_routers[0]->get_dst_comp_node();
        remote_comp_node = fields_routers[0]->get_src_comp_node();
    }
    strcpy(remote_comp_full_name, remote_comp_node->get_comp_full_name());
    remote_comp_node_updated = false;
    timer_not_bypassed = false;
    comp_id = local_comp_node->get_comp_id();
    comp_node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false, "in Runtime_trans_algorithm::Runtime_trans_algorithm");
    current_proc_local_id = local_comp_node->get_current_proc_local_id();
    current_proc_global_id = comp_comm_group_mgt_mgr->get_current_proc_global_id();
    time_mgr = components_time_mgrs->get_time_mgr(comp_id);
    EXECUTION_REPORT(REPORT_ERROR, -1, time_mgr != NULL, "software error in Runtime_trans_algorithm::Runtime_trans_algorithm: wrong time mgr: %x: %d: %d: %s : %s: %s %s %s", comp_id, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false, "C-Coupler native code get time manager")->get_current_proc_local_id(), current_proc_local_id, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false, "C-Coupler native code get time manager")->get_comp_full_name(), local_comp_node->get_comp_full_name(), remote_comp_node->get_comp_full_name(), fields_routers[0]->get_src_comp_node()->get_comp_name(), fields_routers[0]->get_dst_comp_node()->get_comp_name());
    num_remote_procs = remote_comp_node->get_num_procs();
    num_local_procs = local_comp_node->get_num_procs();
    remote_proc_ranks_in_union_comm = new int [num_remote_procs];
    memcpy(remote_proc_ranks_in_union_comm, ranks, num_remote_procs*sizeof(int));
    sender_time_has_matched = false;

#ifndef USE_ONE_SIDED_MPI
    request = new MPI_Request[num_remote_procs];
    is_first_run = true;
#endif
    transfer_size_with_remote_procs = new int [num_remote_procs];
    send_displs_in_remote_procs = new int [num_remote_procs];
    recv_displs_in_current_proc = new int [num_remote_procs];
    field_grids_num_lev = new long [num_transfered_fields];
    fields_data_type_sizes = new int [num_transfered_fields];
    is_V1D_sub_grid_after_H2D_sub_grid =  new bool [num_transfered_fields];

    memset(transfer_size_with_remote_procs, 0, sizeof(int)*num_remote_procs);
    memset(send_displs_in_remote_procs, 0, sizeof(int)*num_remote_procs);
    memset(recv_displs_in_current_proc, 0, sizeof(int)*num_remote_procs);

    only_have_no_decomp_data = true;
    for (int j = 0; j < num_remote_procs; j ++) {
        for (int i = 0; i < num_transfered_fields; i ++) {
            fields_data_type_sizes[i] = get_data_type_size(fields_mem[i]->get_data_type());
            is_V1D_sub_grid_after_H2D_sub_grid[i] = true;
            if (fields_routers[i]->get_num_dimensions() == 0) 
                field_grids_num_lev[i] = 1;
            else {
                field_grids_num_lev[i] = original_grid_mgr->get_num_grid_levels(fields_mem[i]->get_grid_id());
                only_have_no_decomp_data = false;
                transfer_size_with_remote_procs[j] += fields_routers[i]->get_num_elements_transferred_with_remote_proc(send_or_receive, j) * fields_data_type_sizes[i] * field_grids_num_lev[i];
                is_V1D_sub_grid_after_H2D_sub_grid[i] = original_grid_mgr->is_V1D_sub_grid_after_H2D_sub_grid(fields_mem[i]->get_grid_id());
            }    
        }
        if (transfer_size_with_remote_procs[j] > 0)
            index_remote_procs_with_common_data.push_back(j);
    }

    if (only_have_no_decomp_data) {
        if (send_or_receive) {
            num_remote_procs_related = num_remote_procs / num_local_procs;
            if (current_proc_local_id < (num_remote_procs % num_local_procs))
                num_remote_procs_related += 1;
            remote_proc_idx_begin = current_proc_local_id;
        }
        else {
            num_remote_procs_related = 1;
            remote_proc_idx_begin = current_proc_local_id % num_remote_procs;
        }

        for (int i = 0; i < num_remote_procs_related; i ++) {
            int remote_proc_idx = remote_proc_idx_begin + i * num_local_procs;
            for (int j = 0; j < num_transfered_fields; j ++)
                transfer_size_with_remote_procs[remote_proc_idx] += fields_data_type_sizes[j];
            index_remote_procs_with_common_data.push_back(remote_proc_idx);
        }
    }
    else {
        for (int j = 0; j < num_remote_procs; j ++)
            if (transfer_size_with_remote_procs[j] > 0)
                for (int i = 0; i < num_transfered_fields; i ++)
                    if (fields_routers[i]->get_num_dimensions() == 0)
                        transfer_size_with_remote_procs[j] += fields_data_type_sizes[i];
    }

    int * total_transfer_size_with_remote_procs = new int [num_local_procs * num_remote_procs];
    if (send_or_receive) {
        MPI_Allgather(transfer_size_with_remote_procs, num_remote_procs, MPI_INT, total_transfer_size_with_remote_procs, num_remote_procs, MPI_INT, local_comp_node->get_comm_group());
        for (int i = 0; i < current_proc_local_id; i ++) {
            for (int j = 0; j < num_remote_procs; j ++) {
                send_displs_in_remote_procs[j] += total_transfer_size_with_remote_procs[i*num_remote_procs+j] + 4*sizeof(long);
             }
        }
        for (int j = 0; j < num_remote_procs; j ++)
            send_displs_in_remote_procs[j] += sizeof(long)*4;
    }

    recv_displs_in_current_proc[0] = sizeof(long)*4;
    for (int i = 1; i < num_remote_procs; i ++)
        recv_displs_in_current_proc[i] = recv_displs_in_current_proc[i-1] + transfer_size_with_remote_procs[i-1] + 4*sizeof(long);

    current_receive_field_sender_time = -1;
    last_receive_field_sender_time = -1;
    data_buf_size = 0;
    for (int j = 0; j < num_remote_procs; j ++) 
        data_buf_size += transfer_size_with_remote_procs[j];

    total_buf_size = data_buf_size + (4*num_remote_procs + 4) * sizeof(long);
    total_buf = (char*) (new long[(total_buf_size+sizeof(long)-1)/sizeof(long)]);
    send_tag_buf = (long *) total_buf;
    temp_receive_data_buffer = (char*)(new long [(data_buf_size+sizeof(long)-1)/sizeof(long)]);

    for (int i = 0; i < 4; i ++)
        send_tag_buf[i] = -1;
    for (int i = 0; i < num_remote_procs; i ++) {
        tag_buf = (long *) (total_buf + recv_displs_in_current_proc[i]);
        for (int j = 0; j < 2; j ++)
            tag_buf[j] = -1;
    }

    num_recv_procs_related = 0;
    recv_proc_start = -1;
    if (send_or_receive) {
        int color = (index_remote_procs_with_common_data.size() > 0) ? 1 : 0;
        MPI_Comm_split(local_comp_node->get_comm_group(), color, 0, &sub_comm);
        if (index_remote_procs_with_common_data.size() > 0) {
            for (int j = 0; j < num_remote_procs; j ++)
                for (int i = 0; i < num_local_procs; i ++)
                    if (total_transfer_size_with_remote_procs[i*num_remote_procs+j] > 0){
                        index_recv_procs_with_common_data.push_back(j);
                        break;
                    }

            int sub_rank, sub_size;
            MPI_Comm_rank(sub_comm, &sub_rank);
            MPI_Comm_size(sub_comm, &sub_size);
            num_recv_procs_related = index_recv_procs_with_common_data.size() / sub_size;
            if (sub_rank < index_recv_procs_with_common_data.size()%sub_size)
                num_recv_procs_related ++;
            recv_proc_start = sub_rank * num_recv_procs_related;
            if (sub_rank >= index_recv_procs_with_common_data.size()%sub_size)
                recv_proc_start += index_recv_procs_with_common_data.size() % sub_size;
        }
    }

    delete [] total_transfer_size_with_remote_procs;
}


Runtime_trans_algorithm::~Runtime_trans_algorithm()
{
    delete [] fields_mem;
    delete [] fields_data_buffers;
    delete [] fields_routers;
    delete [] field_grids_num_lev;
    delete [] fields_data_type_sizes;
    delete [] is_V1D_sub_grid_after_H2D_sub_grid;
    delete [] total_buf;
    delete [] transfer_size_with_remote_procs;
    delete [] send_displs_in_remote_procs;
    delete [] recv_displs_in_current_proc;
    delete [] remote_proc_ranks_in_union_comm;
    delete [] temp_receive_data_buffer;
#ifndef USE_ONE_SIDED_MPI
    delete [] request;
#endif
}


void Runtime_trans_algorithm::pass_transfer_parameters(long current_remote_fields_time, int bypass_counter)
{
    this->current_remote_fields_time = current_remote_fields_time;
    this->bypass_counter = bypass_counter;
}


bool Runtime_trans_algorithm::set_local_tags()
{
    MPI_Win_lock(MPI_LOCK_SHARED, current_proc_id_union_comm, 0, data_win);
    send_tag_buf[0] = current_field_local_recv_count;
    send_tag_buf[1] = time_mgr->get_current_full_time();
    send_tag_buf[2] = (long) time_mgr->get_runtype_mark();
    send_tag_buf[3] = time_mgr->get_restart_full_time();
    current_field_local_recv_count ++;
    MPI_Win_unlock(current_proc_id_union_comm, data_win);

    return true;
}


bool Runtime_trans_algorithm::is_remote_data_buf_ready(bool bypass_timer)
{
    long temp_field_remote_recv_count = -100;
    double time1, time2;

    if (index_remote_procs_with_common_data.size() == 0)
        return true;

    wtime(&time1);
    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
        int remote_proc_index = index_remote_procs_with_common_data[i];
        if (transfer_size_with_remote_procs[remote_proc_index] > 0) {
             if (remote_comp_node_updated && last_receive_sender_time < remote_comp_node->get_proc_latest_model_time(remote_proc_index)) {
                EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Can bypass MPI_Get for proc %d", remote_proc_index);
                continue;
            }
            int remote_proc_id = remote_proc_ranks_in_union_comm[remote_proc_index];
            MPI_Win_lock(MPI_LOCK_SHARED, remote_proc_id, 0, data_win);
            MPI_Get(send_tag_buf, sizeof(long)*4, MPI_CHAR, remote_proc_id, 0, sizeof(long)*4, MPI_CHAR, data_win);
            MPI_Win_unlock(remote_proc_id, data_win);
            if (remote_comp_node_updated)
                remote_comp_node->set_proc_latest_model_time(remote_proc_index, send_tag_buf[1]);
            if (send_tag_buf[0] != -1 && send_tag_buf[0] != last_field_remote_recv_count + 1)
                return false;
            if (temp_field_remote_recv_count == -100)
                temp_field_remote_recv_count = send_tag_buf[0];
            if (temp_field_remote_recv_count != send_tag_buf[0])
                return false;
        }
    }

    if (temp_field_remote_recv_count == -1) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, last_field_remote_recv_count == -1 || last_field_remote_recv_count == 0, "Software error in Runtime_trans_algorithm::is_remote_data_buf_ready");
        if (last_field_remote_recv_count != -1) 
            return false;
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Remote buffer component \"%s\" is ready for receiving data: %ld vs %ld vs %ld : %d: %d  %ld", remote_comp_full_name, temp_field_remote_recv_count, last_field_remote_recv_count, last_receive_sender_time, bypass_counter, send_tag_buf[2], send_tag_buf[3]);    

    last_field_remote_recv_count ++;

    wtime(&time2);
    local_comp_node->get_performance_timing_mgr()->performance_timing_add(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_SEND_QUERRY, -1, remote_comp_full_name, time2-time1);
    
    return true;
}


void Runtime_trans_algorithm::receive_data_in_temp_buffer()
{
    bool is_ready = true;
    double time1, time2, time3;


    if (index_remote_procs_with_common_data.size() == 0)
        return;

    if (timer_not_bypassed && last_history_receive_buffer_index != -1) {
        int comp_min_remote_lag_seconds = comp_node->get_min_remote_lag_seconds();
        long current_receiver_full_seconds = ((long)time_mgr->get_current_num_elapsed_day())*86400 + time_mgr->get_current_second();
        long current_sender_full_seconds = time_mgr->get_elapsed_day_from_full_time(current_receive_field_sender_time%((long)10000000000000000))*86400 + (current_receive_field_sender_time%((long)100000));
        if (current_sender_full_seconds + 2*comp_min_remote_lag_seconds > current_receiver_full_seconds)
            return;
    }

#ifndef USE_ONE_SIDED_MPI
    local_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV, -1, remote_comp_full_name);
    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
        int remote_proc_index = index_remote_procs_with_common_data[i];
        if (transfer_size_with_remote_procs[remote_proc_index] == 0) 
            continue;
        data_buf = (void *) (total_buf + recv_displs_in_current_proc[remote_proc_index]);
        int remote_proc_id = remote_proc_ranks_in_union_comm[remote_proc_index];
        MPI_Irecv((char *)data_buf, 4*sizeof(long)+transfer_size_with_remote_procs[remote_proc_index], MPI_CHAR, remote_proc_id, comm_tag, union_comm, &request[i]);
    }    
    local_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV, -1, remote_comp_full_name);
    local_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV_WAIT, -1, remote_comp_full_name);
    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
        int remote_proc_index = index_remote_procs_with_common_data[i];
        if (transfer_size_with_remote_procs[remote_proc_index] == 0) 
            continue;
        MPI_Status state;
        MPI_Wait(&request[i], &state);
    }
    local_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV_WAIT, -1, remote_comp_full_name);
#endif

    wtime(&time1);

#ifdef USE_ONE_SIDED_MPI
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, current_proc_id_union_comm, 0, data_win);
#endif
    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
        int remote_proc_index = index_remote_procs_with_common_data[i];
        tag_buf = (long *) (total_buf + recv_displs_in_current_proc[remote_proc_index]);
        if (i == 0) {
            current_receive_field_sender_time = tag_buf[0];
            current_receive_field_usage_time = tag_buf[1];
        }
        else is_ready = is_ready && (current_receive_field_sender_time == tag_buf[0]);
    }
#ifdef USE_ONE_SIDED_MPI
    MPI_Win_unlock(current_proc_id_union_comm, data_win);
#endif

    if (!is_ready) {
#ifndef USE_ONE_SIDED_MPI
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error1 in MPI_send/recv implementation in Runtime_trans_algorithm::receive_data_in_temp_buffer");
#endif
        return;
    }

    if (last_receive_field_sender_time == current_receive_field_sender_time) {
#ifndef USE_ONE_SIDED_MPI
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error2 in MPI_send/recv implementation in Runtime_trans_algorithm::receive_data_in_temp_buffer");
#endif
        return;
    }

    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
        int remote_proc_index = index_remote_procs_with_common_data[i];
        tag_buf = (long *) (total_buf + recv_displs_in_current_proc[remote_proc_index]);
        if (tag_buf[2] != -1) {
            if (tag_buf[2] == RUNTYPE_MARK_INITIAL || tag_buf[2] == RUNTYPE_MARK_HYBRID) 
                EXECUTION_REPORT(REPORT_ERROR, comp_id, time_mgr->get_runtype_mark() == RUNTYPE_MARK_INITIAL || time_mgr->get_runtype_mark() == RUNTYPE_MARK_HYBRID, "Inconsistency of run type between component models is detected: the component model \"%s\" is in an initial run or hybrid run, while the component model \"%s\" is in a continue run or branch run. Please verify.", remote_comp_full_name, local_comp_node->get_comp_full_name());
            else {        
                EXECUTION_REPORT(REPORT_ERROR, comp_id, time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE || time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH, "Inconsistency of run type between component models is detected: the component model \"%s\" is in an initial run or hybrid run, while the component model \"%s\" is in a continue run or branch run. Please verify.", local_comp_node->get_comp_full_name(), remote_comp_full_name);
                if (time_mgr->get_restart_full_time() != -1 && tag_buf[3] != -1)
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, time_mgr->get_restart_full_time() == tag_buf[3], "The restart time between the two component models \"%s\" and \"%s\" are inconsistent: %ld vs %ld. Please verify.", local_comp_node->get_comp_full_name(), remote_comp_full_name, time_mgr->get_restart_full_time(), tag_buf[3]);
            }
        }
    }

#ifdef USE_ONE_SIDED_MPI
    wtime(&time2);    
    local_comp_node->get_performance_timing_mgr()->performance_timing_add(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV_QUERRY, -1, remote_comp_full_name, time2-time1);
#endif    

    int empty_history_receive_buffer_index = -1;
    if (last_history_receive_buffer_index != -1) {
        for (int i = 0; i < history_receive_fields_mem.size(); i ++) {
            int index_iter = (last_history_receive_buffer_index+i) % history_receive_fields_mem.size();
            if (!history_receive_buffer_status[index_iter]) {
                empty_history_receive_buffer_index = index_iter;
                break;
            }
        }
    }
    if (empty_history_receive_buffer_index == -1) {
        std::vector<bool> temp_history_receive_buffer_status;
        std::vector<long> temp_history_receive_sender_time;
        std::vector<long> temp_history_receive_usage_time;
        std::vector<std::vector<Field_mem_info *> > temp_history_receive_fields_mem;
        for (int i = 0; i < history_receive_fields_mem.size(); i ++) {
            int index_iter = (last_history_receive_buffer_index+i) % history_receive_fields_mem.size();
            temp_history_receive_buffer_status.push_back(history_receive_buffer_status[index_iter]);
            temp_history_receive_sender_time.push_back(history_receive_sender_time[index_iter]);
            temp_history_receive_usage_time.push_back(history_receive_usage_time[index_iter]);
            temp_history_receive_fields_mem.push_back(history_receive_fields_mem[index_iter]);
        }
        history_receive_buffer_status.clear();
        history_receive_sender_time.clear();
        history_receive_usage_time.clear();
        history_receive_fields_mem.clear();
        for (int i = 0; i < temp_history_receive_fields_mem.size(); i ++) {
            history_receive_buffer_status.push_back(temp_history_receive_buffer_status[i]);
            history_receive_sender_time.push_back(temp_history_receive_sender_time[i]);
            history_receive_usage_time.push_back(temp_history_receive_usage_time[i]);
            history_receive_fields_mem.push_back(temp_history_receive_fields_mem[i]);
        }
        last_history_receive_buffer_index = 0;
        empty_history_receive_buffer_index = history_receive_buffer_status.size();
        history_receive_buffer_status.push_back(false);
        history_receive_sender_time.push_back(-1);
        history_receive_usage_time.push_back(-1);
        std::vector<Field_mem_info *> new_receive_fields_mem;
        for (int i = 0; i < num_transfered_fields; i ++) {
            new_receive_fields_mem.push_back(memory_manager->alloc_mem(fields_mem[i], BUF_MARK_DATA_TRANSFER, history_receive_fields_mem.size(), fields_mem[i]->get_data_type(), false));
            for (int j = 0; j < history_receive_fields_mem.size(); j ++)
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, history_receive_fields_mem[j][i] != new_receive_fields_mem[i], "Software error in Runtime_trans_algorithm::receive_data_in_temp_buffer");
        }    
        history_receive_fields_mem.push_back(new_receive_fields_mem);
    }

    history_receive_buffer_status[empty_history_receive_buffer_index] = true;
    history_receive_sender_time[empty_history_receive_buffer_index] = current_receive_field_sender_time;
    history_receive_usage_time[empty_history_receive_buffer_index] = current_receive_field_usage_time;
    last_receive_field_sender_time = current_receive_field_sender_time;

#ifdef USE_ONE_SIDED_MPI
    MPI_Win_lock(MPI_LOCK_SHARED, current_proc_id_union_comm, 0, data_win);
#endif
    int offset = 0;
    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
        int remote_proc_index = index_remote_procs_with_common_data[i];
        if (transfer_size_with_remote_procs[remote_proc_index] == 0) 
            continue;
        data_buf = (void *) (total_buf + recv_displs_in_current_proc[remote_proc_index] + 4*sizeof(long));
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, recv_displs_in_current_proc[remote_proc_index] + 4*sizeof(long) >= 0 && recv_displs_in_current_proc[remote_proc_index] + 4*sizeof(long) + transfer_size_with_remote_procs[remote_proc_index] <= total_buf_size, "Software error in Runtime_trans_algorithm::receive_data_in_temp_buffer: %d + %d vs %d", recv_displs_in_current_proc[remote_proc_index] + 4*sizeof(long), transfer_size_with_remote_procs[remote_proc_index], total_buf_size);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, offset >= 0 && offset + transfer_size_with_remote_procs[remote_proc_index] <= data_buf_size, "Software error in Runtime_trans_algorithm::receive_data_in_temp_buffer: %d + %d vs %d", offset, transfer_size_with_remote_procs[remote_proc_index], data_buf_size);
        memcpy(temp_receive_data_buffer+offset, data_buf, transfer_size_with_remote_procs[remote_proc_index]);
        offset += transfer_size_with_remote_procs[remote_proc_index];
    }    
#ifdef USE_ONE_SIDED_MPI
    MPI_Win_unlock(current_proc_id_union_comm, data_win);
#endif
    
    offset = 0;
    for (int i = 0; i < num_remote_procs; i ++) {
        if (transfer_size_with_remote_procs[i] == 0) 
            continue;
        int old_offset = offset;
        //int offset = recv_displs_in_current_proc[i];
        for (int j = 0; j < num_transfered_fields; j ++) {
            if (fields_routers[j]->get_num_dimensions() == 0) {
                memcpy(history_receive_fields_mem[empty_history_receive_buffer_index][j]->get_data_buf(), temp_receive_data_buffer + offset, fields_data_type_sizes[j]);
                offset += fields_data_type_sizes[j];
            }
            else unpack_MD_data(temp_receive_data_buffer, i, j, history_receive_fields_mem[empty_history_receive_buffer_index][j]->get_data_buf(), &offset);
        }    
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, offset - old_offset == transfer_size_with_remote_procs[i], "C-Coupler software error in recv of runtime_trans_algorithm.");
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Get receiving data from component \"%s\" (at time %ld) into temp buffer", remote_comp_full_name, last_receive_field_sender_time);

#ifdef USE_ONE_SIDED_MPI
    set_local_tags();
    wtime(&time3);
    local_comp_node->get_performance_timing_mgr()->performance_timing_add(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV, -1, remote_comp_full_name, time3-time2);
#endif    
}


bool Runtime_trans_algorithm::run(bool bypass_timer)
{
    if (!bypass_timer)
        timer_not_bypassed = true;
    if (send_or_receive)
        return send(bypass_timer);
    else return recv(bypass_timer);
}


bool Runtime_trans_algorithm::send(bool bypass_timer)
{
    if (!remote_comp_node_updated) {
        remote_comp_node = comp_comm_group_mgt_mgr->load_comp_info_from_XML(local_comp_node->get_comp_id(), remote_comp_full_name, local_comp_node->get_comm_group());
        Comp_comm_group_mgt_node *existing_remote_comp_node = comp_comm_group_mgt_mgr->search_global_node(remote_comp_full_name);
        if (existing_remote_comp_node != NULL) {
            delete remote_comp_node;
            remote_comp_node = existing_remote_comp_node;
        }
        remote_comp_node_updated = true;
        remote_comp_node->allocate_proc_latest_model_time();
    }
#ifndef USE_ONE_SIDED_MPI
    comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_SEND_WAIT, -1, remote_comp_full_name);
    if (!is_first_run) {
        for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
            int remote_proc_index = index_remote_procs_with_common_data[i];
            MPI_Status state;
            MPI_Wait(&request[i], &state);
        }
    }
    comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_SEND_WAIT, -1, remote_comp_full_name);
    is_first_run = false;
#endif

    if (index_remote_procs_with_common_data.size() > 0) {
        preprocess();
#ifdef USE_ONE_SIDED_MPI
        if (!is_remote_data_buf_ready(bypass_timer)) {
            inout_interface_mgr->runtime_receive_algorithms_receive_data();
            return false;
        }
#endif
    }

    for (int j = 0; j < num_transfered_fields; j ++) {
        fields_mem[j]->check_field_sum("before sending data");
        fields_mem[j]->use_field_values("before sending data");
    }  

    if (index_remote_procs_with_common_data.size() == 0)
        return true;

    local_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_SEND, -1, remote_comp_full_name);

    long current_full_time = time_mgr->get_current_full_time();
    int offset = 0;
    //for (int i = 0; i < num_remote_procs; i ++) {
    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++) {
        int remote_proc_index = index_remote_procs_with_common_data[i];
        //if (transfer_size_with_remote_procs[remote_proc_index] == 0) continue;

        offset = 0;
        int old_offset = offset;
        data_buf = (void *) (total_buf + recv_displs_in_current_proc[remote_proc_index] + 4*sizeof(long));
        if (transfer_size_with_remote_procs[remote_proc_index] > 0)
            for (int j = 0; j < num_transfered_fields; j ++) {
                if (fields_routers[j]->get_num_dimensions() == 0) {
                    memcpy((char *)data_buf + offset, fields_data_buffers[j], fields_data_type_sizes[j]);
                    offset += fields_data_type_sizes[j];
                }
                else pack_MD_data(remote_proc_index, j, &offset);
            }

        tag_buf = (long *) (total_buf + recv_displs_in_current_proc[remote_proc_index]);
        if (bypass_timer) {
            tag_buf[0] = current_full_time + (bypass_counter%8)*((long)10000000000000000);
            tag_buf[1] = -999;
        }
        else {
            tag_buf[0] = current_full_time;
            tag_buf[1] = current_remote_fields_time;
        }
        tag_buf[2] = (long) time_mgr->get_runtype_mark();
        tag_buf[3] = time_mgr->get_restart_full_time();

        int remote_proc_id = remote_proc_ranks_in_union_comm[remote_proc_index];

#ifndef USE_ONE_SIDED_MPI
        MPI_Isend(tag_buf, 4*sizeof(long)+transfer_size_with_remote_procs[remote_proc_index], MPI_CHAR, remote_proc_id, comm_tag, union_comm, &request[i]);
#else
        MPI_Win_lock(MPI_LOCK_SHARED, remote_proc_id, 0, data_win);
        MPI_Put(tag_buf, 4*sizeof(long)+transfer_size_with_remote_procs[remote_proc_index], MPI_CHAR, remote_proc_id, send_displs_in_remote_procs[remote_proc_index], 4*sizeof(long)+transfer_size_with_remote_procs[remote_proc_index], MPI_CHAR, data_win);
        MPI_Win_unlock(remote_proc_id, data_win);
#endif

        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Set remote tag to component \"%s\": %ld %ld", remote_comp_full_name, tag_buf[0], tag_buf[1]);
    }

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, offset <= data_buf_size, "Software error in Runtime_trans_algorithm::send: wrong data_buf_size: %d vs %d", offset, data_buf_size);

    if (bypass_timer)
        last_receive_sender_time = (bypass_counter%8)*((long)10000000000000000);
    else last_receive_sender_time = current_remote_fields_time;

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish sending data to component \"%s\": %d", remote_comp_full_name, comm_tag);

    local_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_SEND, -1, remote_comp_full_name);

    return true;
}


bool Runtime_trans_algorithm::recv(bool bypass_timer)
{
    bool received_data_ready = false;
    

#ifdef USE_ONE_SIDED_MPI
    local_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV_WAIT, -1, remote_comp_full_name);
#endif
    if (bypass_timer) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Bypass timer to begin to receive data from component \"%s\": %ld: %d: %d", remote_comp_full_name, current_remote_fields_time, bypass_counter, comm_tag);
    }    
    else EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Use timer to begin to receive data from component \"%s\": %ld %d", remote_comp_full_name, current_remote_fields_time, comm_tag);

    if (index_remote_procs_with_common_data.size() > 0) {
        preprocess();
#ifndef USE_ONE_SIDED_MPI
        receive_data_in_temp_buffer();
#else
        while (!received_data_ready) {
            receive_data_in_temp_buffer();
            received_data_ready = last_history_receive_buffer_index != -1 && history_receive_buffer_status[last_history_receive_buffer_index];
            if (!received_data_ready)
                inout_interface_mgr->runtime_receive_algorithms_receive_data();
        }
#endif
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, last_history_receive_buffer_index >= 0, "Software error with last_history_receive_buffer_index: %d", last_history_receive_buffer_index);
        for (int j = 0; j < num_transfered_fields; j ++)
            memcpy(fields_mem[j]->get_data_buf(), history_receive_fields_mem[last_history_receive_buffer_index][j]->get_data_buf(), fields_mem[j]->get_size_of_field()*get_data_type_size(fields_mem[j]->get_data_type()));
    }

    if (index_remote_procs_with_common_data.size() > 0)
        last_receive_sender_time = history_receive_sender_time[last_history_receive_buffer_index];
    else if (bypass_timer)
        last_receive_sender_time = (bypass_counter%8)*((long)10000000000000000);
    else last_receive_sender_time = current_remote_fields_time;

    for (int j = 0; j < num_transfered_fields; j ++) {
         fields_mem[j]->check_field_sum("after receiving data");
         fields_mem[j]->define_field_values(false);
    }    

    if (index_remote_procs_with_common_data.size() > 0) {
        history_receive_buffer_status[last_history_receive_buffer_index] = false;
        last_history_receive_buffer_index = (last_history_receive_buffer_index+1) % history_receive_buffer_status.size();
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish receiving data from component \"%s\" at the remote model time %ld vs %ld", remote_comp_full_name, last_receive_sender_time, current_remote_fields_time);

#ifdef USE_ONE_SIDED_MPI
    local_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_RECV_WAIT, -1, remote_comp_full_name);
#endif

    return true;
}


long Runtime_trans_algorithm::get_history_receive_sender_time()
{
    return last_receive_sender_time;
}


void Runtime_trans_algorithm::preprocess()
{
    for (int i = 0; i < index_remote_procs_with_common_data.size(); i ++)
        transfer_size_with_remote_procs[index_remote_procs_with_common_data[i]] = 0;

    for (int i = 0; i < num_transfered_fields; i ++) {
        for (int j = 0; j < index_remote_procs_with_common_data.size(); j ++) {
            int remote_proc_index = index_remote_procs_with_common_data[j];
            if (fields_routers[i]->get_num_dimensions() == 0)
                transfer_size_with_remote_procs[remote_proc_index] += fields_data_type_sizes[i];
            else transfer_size_with_remote_procs[remote_proc_index] += fields_routers[i]->get_num_elements_transferred_with_remote_proc(send_or_receive, remote_proc_index) * fields_data_type_sizes[i] * field_grids_num_lev[i];
        }
    }
}


void Runtime_trans_algorithm::pack_MD_data(int remote_proc_index, int field_index, int * offset)
{
    int num_segments;
    int *segment_starts, *num_elements_in_segments;
    int i, j;
    int field_2D_size;


    num_segments = fields_routers[field_index]->get_num_local_indx_segments_with_remote_proc(true, remote_proc_index);
    if (num_segments == 0)
        return;

    segment_starts = fields_routers[field_index]->get_local_indx_segment_starts_with_remote_proc(true, remote_proc_index);
    num_elements_in_segments = fields_routers[field_index]->get_local_indx_segment_lengths_with_remote_proc(true, remote_proc_index);
    field_2D_size = fields_routers[field_index]->get_src_decomp_size();
    for (i = 0; i < num_segments; i ++) {
        switch (fields_data_type_sizes[field_index]) {
            case 1:
                pack_segment_data((char*)((char*)data_buf+(*offset)), (char*)fields_data_buffers[field_index], segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            case 2:
                pack_segment_data((short*)((char*)data_buf+(*offset)), (short*)fields_data_buffers[field_index], segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            case 4:
                pack_segment_data((int*)((char*)data_buf+(*offset)), (int*)fields_data_buffers[field_index], segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            case 8:
                pack_segment_data((double*)((char*)data_buf+(*offset)), (double*)fields_data_buffers[field_index], segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            default:
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR,-1, false, "Software error in Runtime_trans_algorithm::pack_MD_data: unsupported data type in runtime transfer algorithm. Please verify.");
                break;
        }
        (*offset) += num_elements_in_segments[i]*field_grids_num_lev[field_index]*fields_data_type_sizes[field_index];
    }
}


void Runtime_trans_algorithm::unpack_MD_data(void *data_buf, int remote_proc_index, int field_index, void *field_data_buffer, int * offset)
{
    int num_segments;
    int *segment_starts, *num_elements_in_segments;
    int i, j;
    int field_2D_size;


    num_segments = fields_routers[field_index]->get_num_local_indx_segments_with_remote_proc(false, remote_proc_index);
    if (num_segments == 0)
        return;

    segment_starts = fields_routers[field_index]->get_local_indx_segment_starts_with_remote_proc(false, remote_proc_index);
    num_elements_in_segments = fields_routers[field_index]->get_local_indx_segment_lengths_with_remote_proc(false, remote_proc_index);
    field_2D_size = fields_routers[field_index]->get_dst_decomp_size();
    for (i = 0; i < num_segments; i ++) {
        switch (fields_data_type_sizes[field_index]) {
            case 1:
                unpack_segment_data((char*)((char*)data_buf+(*offset)), (char*)field_data_buffer, segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            case 2:
                unpack_segment_data((short*)((char*)data_buf+(*offset)), (short*)field_data_buffer, segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            case 4:
                unpack_segment_data((int*)((char*)data_buf+(*offset)), (int*)field_data_buffer, segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            case 8:
                unpack_segment_data((double*)((char*)data_buf+(*offset)), (double*)field_data_buffer, segment_starts[i], num_elements_in_segments[i], field_2D_size, field_grids_num_lev[field_index], is_V1D_sub_grid_after_H2D_sub_grid[field_index]);
                break;
            default:
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR,-1, false, "Software error in Runtime_trans_algorithm::unpack_MD_data: unsupported data type in runtime transfer algorithm. Please verify.");
                break;
        }
        (*offset) += num_elements_in_segments[i]*field_grids_num_lev[field_index]*fields_data_type_sizes[field_index];
    }
}

