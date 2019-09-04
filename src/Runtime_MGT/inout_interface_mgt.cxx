/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "restart_mgt.h"
#include "runtime_cumulate_average_algorithm.h"
#include "runtime_datatype_transformer.h"
#include "inout_interface_mgt.h"
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>


Connection_field_time_info::Connection_field_time_info(Inout_interface *inout_interface, Coupling_timer *timer, int time_step_in_second, int current_year, int current_month, int current_day, int current_second, int inst_or_aver)
{
    this->inout_interface = inout_interface;
    this->timer = timer;
    this->inst_or_aver = inst_or_aver;
    if (IS_TIME_UNIT_SECOND(timer->get_frequency_unit()))
        lag_seconds = timer->get_remote_lag_count();
    else lag_seconds = timer->get_remote_lag_count() * SECONDS_PER_DAY;

    if (!(components_time_mgrs->get_time_mgr(inout_interface->get_comp_id())->get_runtype_mark() == RUNTYPE_MARK_CONTINUE || components_time_mgrs->get_time_mgr(inout_interface->get_comp_id())->get_runtype_mark() == RUNTYPE_MARK_BRANCH)) {
        this->current_year = current_year;
        this->current_month = current_month;
        this->current_day = current_day;
        this->current_second = current_second;
        current_num_elapsed_days = components_time_mgrs->get_time_mgr(inout_interface->get_comp_id())->get_current_num_elapsed_day();
        this->time_step_in_second = time_step_in_second;
        if (components_time_mgrs->get_time_mgr(inout_interface->get_comp_id())->is_timer_on(timer->get_frequency_unit(), timer->get_frequency_count(), timer->get_local_lag_count())) {
            last_timer_num_elapsed_days = current_num_elapsed_days;
			last_timer_date = current_year*10000 + current_month*100 + current_day;
            last_timer_second = current_second;
        }
        else {
            last_timer_num_elapsed_days = -1;
			last_timer_date = -1;
            last_timer_second = -1;
        }
        next_timer_num_elapsed_days = -1;
		next_timer_date = -1;
        next_timer_second = -1;
        timer->get_time_of_next_timer_on(components_time_mgrs->get_time_mgr(inout_interface->get_comp_id()), current_year, current_month, current_day,
                                         current_second, current_num_elapsed_days, time_step_in_second, next_timer_num_elapsed_days, next_timer_date, next_timer_second, true);
    }
}


void Connection_field_time_info::get_time_of_next_timer_on(bool advance)
{
    timer->get_time_of_next_timer_on(components_time_mgrs->get_time_mgr(inout_interface->get_comp_id()), current_year, current_month, current_day,
                                     current_second, current_num_elapsed_days, time_step_in_second, next_timer_num_elapsed_days, next_timer_date, next_timer_second, advance);
}


void Connection_field_time_info::write_restart_mgt_info(Restart_buffer_container *restart_buffer)
{
    restart_buffer->dump_in_data(&current_year, sizeof(int));
    restart_buffer->dump_in_data(&current_month, sizeof(int));
    restart_buffer->dump_in_data(&current_day, sizeof(int));
    restart_buffer->dump_in_data(&current_second, sizeof(int));
    restart_buffer->dump_in_data(&current_num_elapsed_days, sizeof(int));
    restart_buffer->dump_in_data(&last_timer_num_elapsed_days, sizeof(int));
    restart_buffer->dump_in_data(&last_timer_second, sizeof(int));
    restart_buffer->dump_in_data(&next_timer_num_elapsed_days, sizeof(int));
    restart_buffer->dump_in_data(&next_timer_second, sizeof(int));
    restart_buffer->dump_in_data(&time_step_in_second, sizeof(int));
    restart_buffer->dump_in_data(&inst_or_aver, sizeof(int));
    restart_buffer->dump_in_data(&lag_seconds, sizeof(int));
}


void Connection_field_time_info::import_restart_data(Restart_buffer_container *restart_buffer)
{
    int restart_inst_or_aver;

    
    restart_buffer->load_restart_data(&lag_seconds, sizeof(int));
    restart_buffer->load_restart_data(&restart_inst_or_aver, sizeof(int));
    restart_buffer->load_restart_data(&time_step_in_second, sizeof(int));
    restart_buffer->load_restart_data(&next_timer_second, sizeof(int));
    restart_buffer->load_restart_data(&next_timer_num_elapsed_days, sizeof(int));
    restart_buffer->load_restart_data(&last_timer_second, sizeof(int));
    restart_buffer->load_restart_data(&last_timer_num_elapsed_days, sizeof(int));
    restart_buffer->load_restart_data(&current_num_elapsed_days, sizeof(int));
    restart_buffer->load_restart_data(&current_second, sizeof(int));
    restart_buffer->load_restart_data(&current_day, sizeof(int));
    restart_buffer->load_restart_data(&current_month, sizeof(int));
    restart_buffer->load_restart_data(&current_year, sizeof(int));

    EXECUTION_REPORT(REPORT_ERROR, inout_interface->get_comp_id(), inst_or_aver == restart_inst_or_aver, "Error happens when restarting the simulation in a \"continue\" run or a \"branch\" run: the values of \"inst_or aver\" of the coupling interface \"%s\" are not the same with the original value recorded in the restart data file (the original value is %d while the new value used in the restart run is %d)", inout_interface->get_interface_name(), restart_inst_or_aver, inst_or_aver);
}


Connection_coupling_procedure::Connection_coupling_procedure(Inout_interface *inout_interface, Coupling_connection *coupling_connection)
{
    int field_local_index;

    
    this->inout_interface = inout_interface;
    this->coupling_connection = coupling_connection; 
    coupling_connections_dumped = false;
    remote_bypass_counter = -1;
	last_receive_sender_time = CCPL_NULL_LONG;
    is_coupling_time_out_of_execution = false;
    restart_mgr = comp_comm_group_mgt_mgr->search_global_node(inout_interface->get_comp_id())->get_restart_mgr();

    for (int i = 0; i < coupling_connection->fields_name.size(); i ++)
        for (int j=i+1; j < coupling_connection->fields_name.size(); j ++)
            EXECUTION_REPORT(REPORT_ERROR, -1, !words_are_the_same(coupling_connection->fields_name[i], coupling_connection->fields_name[j]), 
                             "Software error in Connection_coupling_procedure::Connection_coupling_procedure: duplicated field name \"%s\" in a coonection", 
                             coupling_connection->fields_name[i]);

    fields_time_info_src = new Connection_field_time_info(inout_interface, coupling_connection->src_timer, coupling_connection->src_time_step_in_second, coupling_connection->src_current_year, coupling_connection->src_current_month, coupling_connection->src_current_day, coupling_connection->src_current_second, -1);
    fields_time_info_dst = new Connection_field_time_info(inout_interface, coupling_connection->dst_timer, coupling_connection->dst_time_step_in_second, coupling_connection->dst_current_year, coupling_connection->dst_current_month, coupling_connection->dst_current_day, coupling_connection->dst_current_second, coupling_connection->dst_inst_or_aver);
    if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT)
        fields_time_info_src->reset_last_timer_info();
    else fields_time_info_dst->reset_last_timer_info();

    for (int i = 0; i < coupling_connection->src_fields_info.size(); i ++) {
        runtime_inner_averaging_algorithm.push_back(NULL);
        runtime_inter_averaging_algorithm.push_back(NULL);
        runtime_remap_algorithms.push_back(NULL);
        runtime_unit_transform_algorithms.push_back(NULL);
        runtime_datatype_transform_algorithms.push_back(NULL);
        if (i < coupling_connection->fields_name.size()) {
            fields_mem_registered.push_back(inout_interface->search_registered_field_instance(coupling_connection->fields_name[i], field_local_index));
            field_interface_local_index.push_back(field_local_index);
        }
        else {
            fields_mem_registered.push_back(coupling_connection->get_bottom_field(inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT, i-coupling_connection->fields_name.size()));
            field_interface_local_index.push_back(-1);
        }
        current_remote_fields_elapsed_time = -1;
		current_remote_fields_time = -1;
        last_remote_fields_time = -1;
        fields_mem_inner_step_averaged.push_back(NULL);
        fields_mem_inter_step_averaged.push_back(NULL);
        if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT) {
            if (fields_time_info_dst->inst_or_aver == USING_AVERAGE_VALUE)
                fields_mem_inner_step_averaged[i] = memory_manager->alloc_mem(fields_mem_registered[i], BUF_MARK_AVERAGED_INNER, coupling_connection->connection_id, NULL, 
                                                                              inout_interface->get_interface_source() == INTERFACE_SOURCE_REGISTER && i < coupling_connection->fields_name.size());
            else fields_mem_inner_step_averaged[i] = fields_mem_registered[i];
            if (fields_time_info_dst->inst_or_aver == USING_AVERAGE_VALUE && !(fields_time_info_dst->lag_seconds == 0 && fields_time_info_dst->timer->is_the_same_with(fields_time_info_src->timer)))
                fields_mem_inter_step_averaged[i] = memory_manager->alloc_mem(fields_mem_registered[i], BUF_MARK_AVERAGED_INTER, coupling_connection->connection_id, NULL, 
                                                                              inout_interface->get_interface_source() == INTERFACE_SOURCE_REGISTER && i < coupling_connection->fields_name.size());
            else fields_mem_inter_step_averaged[i] = fields_mem_inner_step_averaged[i];
        }
        fields_mem_remapped.push_back(NULL);
        fields_mem_datatype_transformed.push_back(NULL);
        fields_mem_unit_transformed.push_back(NULL);
        fields_mem_transfer.push_back(NULL);

        if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT) {
            if (fields_time_info_dst->inst_or_aver == USING_AVERAGE_VALUE)
                runtime_inner_averaging_algorithm[i] = new Runtime_cumulate_average_algorithm(this, fields_mem_registered[i], fields_mem_inner_step_averaged[i]);
            if (fields_mem_inter_step_averaged[i] != fields_mem_inner_step_averaged[i])
                runtime_inter_averaging_algorithm[i] = new Runtime_cumulate_average_algorithm(this, fields_mem_inner_step_averaged[i], fields_mem_inter_step_averaged[i]);
        }
        const char *transfer_data_type = get_data_type_size(coupling_connection->src_fields_info[i]->data_type) <= get_data_type_size(coupling_connection->dst_fields_info[i]->data_type)? 
                                         coupling_connection->src_fields_info[i]->data_type : coupling_connection->dst_fields_info[i]->data_type;
        if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT) {
            if (!words_are_the_same(transfer_data_type, coupling_connection->src_fields_info[i]->data_type)) {
                EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, 
                                 "For field %s, add data type transformation at src from %s to %s\n", 
                                 fields_mem_registered[i]->get_field_name(), coupling_connection->src_fields_info[i]->data_type, transfer_data_type);
                fields_mem_datatype_transformed[i] = memory_manager->alloc_mem(fields_mem_registered[i], BUF_MARK_DATATYPE_TRANS, coupling_connection->connection_id, transfer_data_type, 
                                                                               inout_interface->get_interface_source() == INTERFACE_SOURCE_REGISTER && i < coupling_connection->fields_name.size());
                runtime_datatype_transform_algorithms[i] = new Runtime_datatype_transformer(fields_mem_inter_step_averaged[i], fields_mem_datatype_transformed[i]);
            }    
        }    
        if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT) {
            if (coupling_connection->dst_fields_info[i]->runtime_remapping_weights == NULL || coupling_connection->dst_fields_info[i]->runtime_remapping_weights->get_parallel_remapping_weights() == NULL)
                fields_mem_transfer[i] = memory_manager->alloc_mem(fields_mem_registered[i], BUF_MARK_DATA_TRANSFER, coupling_connection->connection_id, transfer_data_type, 
                                                                   inout_interface->get_interface_source() == INTERFACE_SOURCE_REGISTER && i < coupling_connection->fields_name.size());
            else {
                fields_mem_transfer[i] = memory_manager->alloc_mem(fields_mem_registered[i]->get_field_name(), coupling_connection->dst_fields_info[i]->runtime_remapping_weights->get_src_decomp_info()->get_decomp_id(), 
                                                                   coupling_connection->dst_fields_info[i]->runtime_remapping_weights->get_src_original_grid()->get_grid_id(), BUF_MARK_DATA_TRANSFER^coupling_connection->connection_id, 
                                                                   transfer_data_type, fields_mem_registered[i]->get_unit(), "internal", inout_interface->get_interface_source() == INTERFACE_SOURCE_REGISTER && i < coupling_connection->fields_name.size());
                fields_mem_remapped[i] = memory_manager->alloc_mem(fields_mem_registered[i]->get_field_name(), fields_mem_registered[i]->get_decomp_id(), fields_mem_registered[i]->get_grid_id(), 
                                                                   BUF_MARK_REMAP_NORMAL^coupling_connection->connection_id, transfer_data_type, fields_mem_registered[i]->get_unit(), "internal", 
                                                                   inout_interface->get_interface_source() == INTERFACE_SOURCE_REGISTER && i < coupling_connection->fields_name.size());
                runtime_remap_algorithms[i] = new Runtime_remap_algorithm(coupling_connection->dst_fields_info[i]->runtime_remapping_weights, fields_mem_transfer[i], fields_mem_remapped[i], coupling_connection->connection_id);
            }
            if (!words_are_the_same(transfer_data_type, coupling_connection->dst_fields_info[i]->data_type)) {
                fields_mem_datatype_transformed[i] = memory_manager->alloc_mem(fields_mem_registered[i], BUF_MARK_DATATYPE_TRANS, coupling_connection->connection_id, coupling_connection->dst_fields_info[i]->data_type, 
                                                                               inout_interface->get_interface_source() == INTERFACE_SOURCE_REGISTER && i < coupling_connection->fields_name.size());
                EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, 
                                 "for field %s, add data type transformation at dst from %s to %s: %x %x %x\n", 
                                 fields_mem_registered[i]->get_field_name(), transfer_data_type, coupling_connection->dst_fields_info[i]->data_type, fields_mem_datatype_transformed[i]->get_grid_id(), fields_mem_registered[i]->get_grid_id(), fields_mem_transfer[i]->get_grid_id());
                if (fields_mem_remapped[i] == NULL)
                    runtime_datatype_transform_algorithms[i] = new Runtime_datatype_transformer(fields_mem_transfer[i], fields_mem_datatype_transformed[i]);
                else runtime_datatype_transform_algorithms[i] = new Runtime_datatype_transformer(fields_mem_remapped[i], fields_mem_datatype_transformed[i]);
            }
        }
        if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT) {
            if (fields_mem_remapped[i] != NULL)
                fields_mem_transfer[i] = fields_mem_remapped[i];
            else if (fields_mem_unit_transformed[i] != NULL)
                fields_mem_transfer[i] = fields_mem_unit_transformed[i];
            else if (fields_mem_datatype_transformed[i] != NULL)
                fields_mem_transfer[i] = fields_mem_datatype_transformed[i];
            else fields_mem_transfer[i] = fields_mem_inter_step_averaged[i];
        }
        else {
            Field_mem_info *last_field_instance = fields_mem_transfer[i];
            if (fields_mem_datatype_transformed[i] != NULL)
                last_field_instance = fields_mem_datatype_transformed[i];
            else if (fields_mem_remapped[i] != NULL)
                last_field_instance = fields_mem_remapped[i];
            else last_field_instance = fields_mem_transfer[i];
            runtime_inter_averaging_algorithm[i] = new Runtime_cumulate_average_algorithm(this, last_field_instance, fields_mem_registered[i]);
        }
    }
    
    if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT)
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"Connection_coupling_procedure::Connection_coupling_procedure")->update_min_max_remote_lag_seconds(fields_time_info_dst->lag_seconds);
}


Connection_coupling_procedure::~Connection_coupling_procedure()
{
    delete fields_time_info_src;
    delete fields_time_info_dst;

    for (int i = 0; i < runtime_inner_averaging_algorithm.size(); i ++) {
        if (runtime_inner_averaging_algorithm[i] != NULL)
            delete runtime_inner_averaging_algorithm[i];
        if (runtime_inter_averaging_algorithm[i] != NULL)
            delete runtime_inter_averaging_algorithm[i];
        if (runtime_remap_algorithms[i] != NULL)
            delete runtime_remap_algorithms[i];
        if (runtime_unit_transform_algorithms[i] != NULL)
            delete runtime_unit_transform_algorithms[i];
        if (runtime_datatype_transform_algorithms[i] != NULL)
            delete runtime_datatype_transform_algorithms[i];
    }

    inout_interface_mgr->erase_runtime_receive_algorithm(runtime_data_transfer_algorithm);
    delete runtime_data_transfer_algorithm;
}


void Connection_coupling_procedure::execute(bool bypass_timer, int *field_update_status, const char *annotation)
{
    Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(inout_interface->get_comp_id());
    int lag_seconds;


    finish_status = false;
    transfer_data = false;

    if (!bypass_timer) {    
        Connection_field_time_info *local_fields_time_info, *remote_fields_time_info;        
        if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT) {
            local_fields_time_info = fields_time_info_dst;
            remote_fields_time_info = fields_time_info_src;
            lag_seconds = local_fields_time_info->lag_seconds;
        }
        else {
            local_fields_time_info = fields_time_info_src;
            remote_fields_time_info = fields_time_info_dst;
            lag_seconds = -remote_fields_time_info->lag_seconds;
        }
        time_mgr->get_current_time(local_fields_time_info->current_year, local_fields_time_info->current_month, local_fields_time_info->current_day, local_fields_time_info->current_second, 0, "CCPL internal");
        local_fields_time_info->current_num_elapsed_days = time_mgr->get_current_num_elapsed_day();  
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, local_fields_time_info->inout_interface->get_comp_id(), !((local_fields_time_info->last_timer_num_elapsed_days != -1)) || ((long)local_fields_time_info->current_num_elapsed_days)*100000+local_fields_time_info->current_second >= ((long)local_fields_time_info->last_timer_num_elapsed_days)*100000+local_fields_time_info->last_timer_second,
                         "Software error in Connection_coupling_procedure::execute: current time is earlier than last timer time");
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, local_fields_time_info->inout_interface->get_comp_id(), ((long)local_fields_time_info->current_num_elapsed_days)*100000+local_fields_time_info->current_second <= ((long)local_fields_time_info->next_timer_num_elapsed_days)*100000+local_fields_time_info->next_timer_second,
                         "Error happens when executing the import/export interface \"%s\": it should but not have already been called at any time when its timer is on. Please check the model code with the annotation \"%s\"", 
                         local_fields_time_info->inout_interface->get_interface_name(), annotation_mgr->get_annotation(local_fields_time_info->inout_interface->get_interface_id(), "registering interface"));
        if (time_mgr->is_timer_on(local_fields_time_info->timer->get_frequency_unit(), local_fields_time_info->timer->get_frequency_count(), local_fields_time_info->timer->get_local_lag_count())) {
            if (((long)local_fields_time_info->current_num_elapsed_days)*100000+local_fields_time_info->current_second == ((long)local_fields_time_info->next_timer_num_elapsed_days)*100000+local_fields_time_info->next_timer_second) {
                local_fields_time_info->last_timer_num_elapsed_days = local_fields_time_info->next_timer_num_elapsed_days;
				local_fields_time_info->last_timer_date = local_fields_time_info->next_timer_date;
                local_fields_time_info->last_timer_second = local_fields_time_info->next_timer_second;
                local_fields_time_info->get_time_of_next_timer_on(true);
            }
            while((((long)remote_fields_time_info->current_num_elapsed_days)*((long)SECONDS_PER_DAY))+remote_fields_time_info->current_second+lag_seconds <= (((long)local_fields_time_info->current_num_elapsed_days)*((long)SECONDS_PER_DAY)) + local_fields_time_info->current_second) {
                if (remote_fields_time_info->timer->is_timer_on(remote_fields_time_info->current_year, remote_fields_time_info->current_month, remote_fields_time_info->current_day, remote_fields_time_info->current_second, remote_fields_time_info->current_num_elapsed_days, 
                                                                time_mgr->get_start_year(), time_mgr->get_start_month(), time_mgr->get_start_day(), time_mgr->get_start_second(), time_mgr->get_start_num_elapsed_day())) {
                    remote_fields_time_info->last_timer_num_elapsed_days = remote_fields_time_info->current_num_elapsed_days;
					remote_fields_time_info->last_timer_date = remote_fields_time_info->current_year*10000 + remote_fields_time_info->current_month*100 + remote_fields_time_info->current_day;
                    remote_fields_time_info->last_timer_second = remote_fields_time_info->current_second;
                }    
                time_mgr->advance_time(remote_fields_time_info->current_year, remote_fields_time_info->current_month, remote_fields_time_info->current_day, remote_fields_time_info->current_second, remote_fields_time_info->current_num_elapsed_days,  remote_fields_time_info->time_step_in_second);
            }            
            remote_fields_time_info->get_time_of_next_timer_on(false);
        }
    }
    
    if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT) { 
#ifdef USE_ONE_SIDED_MPI
        ((Runtime_trans_algorithm*)runtime_data_transfer_algorithm)->receive_data_in_temp_buffer();
#endif
        if (bypass_timer) {
            current_remote_fields_elapsed_time = -1;
			current_remote_fields_time = -1;
            if (inout_interface->get_bypass_counter() == 1) {
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, inout_interface->get_comp_id(), !(!words_are_the_same(time_mgr->get_run_type(), RUNTYPE_CONTINUE) && !words_are_the_same(time_mgr->get_run_type(), RUNTYPE_BRANCH)) || last_remote_fields_time == -1, "Software error in Connection_coupling_procedure::execute: wrong last_remote_fields_time 1");
            }
            else if (!words_are_the_same(time_mgr->get_run_type(), RUNTYPE_CONTINUE) && !words_are_the_same(time_mgr->get_run_type(), RUNTYPE_BRANCH))
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, inout_interface->get_comp_id(), (inout_interface->get_bypass_counter() - 1)%8 == remote_bypass_counter, "Software error in Connection_coupling_procedure::execute: wrong last_remote_fields_time 2");
            transfer_data = true;
        }
        else if (!(fields_time_info_dst->current_num_elapsed_days != fields_time_info_dst->last_timer_num_elapsed_days || fields_time_info_dst->current_second != fields_time_info_dst->last_timer_second)) {
            if (fields_time_info_src->last_timer_num_elapsed_days != -1) {
                current_remote_fields_elapsed_time = ((long)fields_time_info_src->last_timer_num_elapsed_days) * 100000 + fields_time_info_src->last_timer_second; 
				current_remote_fields_time = ((long)fields_time_info_src->last_timer_date)*100000 + fields_time_info_src->last_timer_second;
            }	
            if (current_remote_fields_elapsed_time != -1 && !time_mgr->is_time_out_of_execution(current_remote_fields_elapsed_time) && current_remote_fields_elapsed_time != last_remote_fields_time) {
                EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "The import interface \"%s\" will receive remote data at %ld vs %ld", inout_interface->get_interface_name(), current_remote_fields_elapsed_time, last_remote_fields_time);
                last_remote_fields_time = current_remote_fields_elapsed_time;
                transfer_data = true;
            }
             else {
                EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "Do not redundantly receive remote data at %ld vs %ld", last_remote_fields_time, current_remote_fields_elapsed_time);
                if (current_remote_fields_elapsed_time != -1 && time_mgr->is_time_out_of_execution(current_remote_fields_elapsed_time))
                    is_coupling_time_out_of_execution = true;
             }
        }
        if (transfer_data) {
            for (int i = fields_mem_registered.size() - 1; i >= 0; i --)
                if (field_interface_local_index[i] != -1)
                    field_update_status[field_interface_local_index[i]] = transfer_data? 1 : 0;
            bool read_restart_data = (!bypass_timer && !inout_interface->get_is_child_interface() && restart_mgr->is_in_restart_read_window(current_remote_fields_elapsed_time));
            if (!bypass_timer && !inout_interface->get_is_child_interface() && restart_mgr->is_in_restart_read_window(current_remote_fields_elapsed_time)) {
                EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "The import interface \"%s\" will not receive data from the component model \"%s\" that is at the time %ld (the restart time is %ld)", inout_interface->get_interface_name(), coupling_connection->get_src_comp_full_name(), current_remote_fields_elapsed_time, time_mgr->get_restart_full_time());
                for (int i = 0; i < fields_mem_registered.size(); i ++)
                    restart_mgr->read_restart_field_data(fields_mem_registered[i], inout_interface->get_interface_name(), "imported", true, NULL, false, annotation);
                transfer_data = false;
            }
            else {
                runtime_data_transfer_algorithm->pass_transfer_parameters(current_remote_fields_time, inout_interface->get_bypass_counter());
                runtime_data_transfer_algorithm->run(bypass_timer);
                comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, inout_interface->get_interface_name());
                for (int i = fields_mem_registered.size() - 1; i >= 0; i --) {
                        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "data interpolation");
                        if (runtime_remap_algorithms[i] != NULL)
                            runtime_remap_algorithms[i]->run(true);
                        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "data interpolation");
                        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "data type transformation");
                        if (runtime_datatype_transform_algorithms[i] != NULL)
                            runtime_datatype_transform_algorithms[i]->run(true);                                
                        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "data type transformation");
                        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "data average");
                        if (runtime_inter_averaging_algorithm[i] != NULL)
                            runtime_inter_averaging_algorithm[i]->run(true);
                        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "data average");
                }
                comp_comm_group_mgt_mgr->get_global_node_of_local_comp(inout_interface->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, inout_interface->get_interface_name());
                if (!bypass_timer && !inout_interface->get_is_child_interface() && (restart_mgr->is_in_restart_write_window(current_remote_fields_elapsed_time, true))) {
                    EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "Should write the remote data at the remote time %ld and local %ld into the restart data file", current_remote_fields_elapsed_time, time_mgr->get_current_num_elapsed_day()*((long)100000)+time_mgr->get_current_second());
                    for (int i = 0; i < fields_mem_registered.size(); i ++)
                        restart_mgr->write_restart_field_data(fields_mem_registered[i], inout_interface->get_interface_name(), "imported", true);
                }
            }
        }
        finish_status = true;
		if (transfer_data) {
			last_receive_sender_time = runtime_data_transfer_algorithm->get_history_receive_sender_time() % ((long)10000000000000000);
            remote_bypass_counter = runtime_data_transfer_algorithm->get_history_receive_sender_time() / ((long)10000000000000000);
            EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "interface \"%s\" last_receive_sender_time is %ld", inout_interface->get_interface_name(), last_receive_sender_time);
		}
        for (int i = fields_mem_registered.size() - 1; i >= 0; i --) {
            if (!transfer_data)
                continue;
            EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "Bypass counter: remote is %d while local is %d", remote_bypass_counter, inout_interface->get_bypass_counter());
            if (bypass_timer) {
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, inout_interface->get_comp_id(), remote_bypass_counter == (inout_interface->get_bypass_counter()%8), "Error happens when executing the import interface \"%s\" with its timer bypassed (the corresponding input parameter \"bypass_timer\" has been set to true): the data currently obtained by this import interface should be but is not from a timer bypassed execution of the corresponding export interface \"%s\" of the component model \"%s\". Please verify.", inout_interface->get_interface_name(), coupling_connection->src_comp_interfaces[0].second, coupling_connection->src_comp_interfaces[0].first);
            }
            else {
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, inout_interface->get_comp_id(), remote_bypass_counter == 0, "Error happens when executing the import interface \"%s\" with its timer unbypassed (the corresponding input parameter \"bypass_timer\" has been set to false): the data currently obtained by this import interface should be but is not from a timer unbypassed execution of the corresponding export interface \"%s\" of the component model \"%s\". Please verify.", inout_interface->get_interface_name(), coupling_connection->src_comp_interfaces[0].second, coupling_connection->src_comp_interfaces[0].first);
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, inout_interface->get_comp_id(), time_mgr->get_elapsed_day_from_full_time(last_receive_sender_time)*((long)100000) + last_receive_sender_time%100000 == current_remote_fields_elapsed_time, 
                                 "Software error: Error happens when using the timer to call the import interface \"%s\": this interface call does not receive the data from the corresponding export interface \"%s\" from the component model \"%s\" at the right model time (the receiver wants the imported data at %ld but received the imported data at %ld). Please verify. ", 
                                 inout_interface->get_interface_name(), coupling_connection->src_comp_interfaces[0].second, coupling_connection->src_comp_interfaces[0].first, current_remote_fields_elapsed_time, runtime_data_transfer_algorithm->get_history_receive_sender_time());
            }    
        }
        return;
    }
    else {
        for (int i = fields_mem_registered.size() - 1; i >= 0; i --) {
            if (bypass_timer) {
                current_remote_fields_elapsed_time = -1;
				current_remote_fields_time = -1;
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(i == fields_mem_registered.size() - 1 && !words_are_the_same(time_mgr->get_run_type(), RUNTYPE_CONTINUE) && !words_are_the_same(time_mgr->get_run_type(), RUNTYPE_BRANCH)) || last_remote_fields_time == -1, "Software error in Connection_coupling_procedure::execute: wrong last_remote_fields_time");
                transfer_data = true;
                if (runtime_inner_averaging_algorithm[i] != NULL)
                    runtime_inner_averaging_algorithm[i]->run(true);
                if (runtime_inter_averaging_algorithm[i] != NULL)
                    runtime_inter_averaging_algorithm[i]->run(true);
                if (runtime_datatype_transform_algorithms[i] != NULL)
                    runtime_datatype_transform_algorithms[i]->run(true);
            }
            else {
                Coupling_timer *dst_timer = fields_time_info_dst->timer;
                Coupling_timer *src_timer = fields_time_info_src->timer;            
                lag_seconds = -fields_time_info_dst->lag_seconds;
                if (fields_time_info_src->current_num_elapsed_days != fields_time_info_src->last_timer_num_elapsed_days || fields_time_info_src->current_second != fields_time_info_src->last_timer_second) {
                    if (fields_time_info_dst->inst_or_aver == USING_AVERAGE_VALUE)
                        runtime_inner_averaging_algorithm[i]->run(false);
                    continue;
                }
                if (runtime_inner_averaging_algorithm[i] != NULL)
                    runtime_inner_averaging_algorithm[i]->run(true);
                if (((long)fields_time_info_src->current_num_elapsed_days)*SECONDS_PER_DAY+fields_time_info_src->current_second == ((long)fields_time_info_dst->last_timer_num_elapsed_days)*SECONDS_PER_DAY+fields_time_info_dst->last_timer_second+lag_seconds) {
                    current_remote_fields_elapsed_time = ((long)fields_time_info_dst->last_timer_num_elapsed_days)*100000 + fields_time_info_dst->last_timer_second;
					current_remote_fields_time = ((long)fields_time_info_dst->last_timer_date)*100000 + fields_time_info_dst->last_timer_second;
                    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(i == fields_mem_registered.size() - 1) || last_remote_fields_time != current_remote_fields_elapsed_time, "Software error in Connection_coupling_procedure::execute: wrong last_remote_fields_time");
                    last_remote_fields_time = current_remote_fields_elapsed_time;
                    if (runtime_inter_averaging_algorithm[i] != NULL)
                        runtime_inter_averaging_algorithm[i]->run(true);
                    if (runtime_datatype_transform_algorithms[i] != NULL) 
                        runtime_datatype_transform_algorithms[i]->run(false);
                    if (!time_mgr->is_time_out_of_execution(current_remote_fields_elapsed_time))
                        transfer_data = true;
                    continue;
                }
                if ((((long)fields_time_info_dst->next_timer_num_elapsed_days)*((long)SECONDS_PER_DAY))+fields_time_info_dst->next_timer_second+lag_seconds < (((long)fields_time_info_src->next_timer_num_elapsed_days)*((long)SECONDS_PER_DAY)) + fields_time_info_src->next_timer_second) {
                    current_remote_fields_elapsed_time = ((long)fields_time_info_dst->next_timer_num_elapsed_days)*100000 + fields_time_info_dst->next_timer_second;
					current_remote_fields_time = ((long)fields_time_info_dst->next_timer_date)*100000 + fields_time_info_dst->next_timer_second;
                    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(i == fields_mem_registered.size() - 1) || last_remote_fields_time != current_remote_fields_elapsed_time, "Software error in Connection_coupling_procedure::execute: wrong last_remote_fields_time");
                    last_remote_fields_time = current_remote_fields_elapsed_time;
                    if (runtime_inter_averaging_algorithm[i] != NULL)
                        runtime_inter_averaging_algorithm[i]->run(true);
                    if (runtime_datatype_transform_algorithms[i] != NULL) 
                        runtime_datatype_transform_algorithms[i]->run(false);
                    if (!time_mgr->is_time_out_of_execution(current_remote_fields_elapsed_time)) {
                        transfer_data = true;
                    }
                }
                else {
                    if (runtime_inter_averaging_algorithm[i] != NULL)
                        runtime_inter_averaging_algorithm[i]->run(false);
                }    
            }
        }
        if (!bypass_timer && !inout_interface->get_is_child_interface() && transfer_data && restart_mgr->is_in_restart_read_window(current_remote_fields_elapsed_time)) {
            EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "The export interface \"%s\" will not send data to the component model \"%s\" that is at the time %ld (the restart time is %ld)", inout_interface->get_interface_name(), coupling_connection->get_dst_comp_full_name(), current_remote_fields_elapsed_time, time_mgr->get_restart_full_time());
            transfer_data = false;
        }
        if (!transfer_data)
            finish_status = true;
        if (transfer_data)
            ((Runtime_trans_algorithm*)runtime_data_transfer_algorithm)->pass_transfer_parameters(current_remote_fields_time, inout_interface->get_bypass_counter());
    }
}


void Connection_coupling_procedure::send_fields(bool bypass_timer)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT && !finish_status && transfer_data, "Software error in Connection_coupling_procedure::send_fields");
    finish_status = runtime_data_transfer_algorithm->run(bypass_timer);
}


Field_mem_info *Connection_coupling_procedure::get_data_transfer_field_instance(int i)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, fields_mem_transfer[i] != NULL, "Software error in Connection_coupling_procedure::get_data_transfer_field_instance");
    return fields_mem_transfer[i]; 
}


Runtime_remapping_weights *Connection_coupling_procedure::get_runtime_remapping_weights(int i) 
{
    return coupling_connection->dst_fields_info[i]->runtime_remapping_weights; 
}


void Connection_coupling_procedure::write_restart_mgt_info(Restart_buffer_container *restart_buffer)
{
    int temp_int;

    
    fields_time_info_src->write_restart_mgt_info(restart_buffer);
    fields_time_info_dst->write_restart_mgt_info(restart_buffer);
    restart_buffer->dump_in_data(&last_remote_fields_time, sizeof(long));
    restart_buffer->dump_in_data(&current_remote_fields_elapsed_time, sizeof(long));
    for (int i = fields_mem_registered.size()-1; i >=0; i --)
        restart_buffer->dump_in_string(fields_mem_registered[i]->get_field_name(), -1);
    if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT) {
        for (int i = runtime_inner_averaging_algorithm.size()-1; i >= 0; i --) {
            if (runtime_inner_averaging_algorithm[i] != NULL) {
                runtime_inner_averaging_algorithm[i]->restart_write(restart_buffer, "aver_inner");
                temp_int = 1;
            }
            else  temp_int = 0;
            restart_buffer->dump_in_data(&temp_int, sizeof(int));
        }
        for (int i = runtime_inter_averaging_algorithm.size()-1; i >= 0; i --) {
            if (runtime_inter_averaging_algorithm[i] != NULL) {
                runtime_inter_averaging_algorithm[i]->restart_write(restart_buffer, "aver_inter");
                temp_int = 1;
            }
            else  temp_int = 0;
            restart_buffer->dump_in_data(&temp_int, sizeof(int));
        }
    }
    temp_int = fields_mem_registered.size();
    restart_buffer->dump_in_data(&temp_int, sizeof(int));    
}



void Connection_coupling_procedure::import_restart_data(Restart_buffer_container *restart_buffer)
{
    int num_total_fields, temp_int, i, j;
    long str_size, temp_long;
    char restart_field_name[NAME_STR_SIZE];


    restart_buffer->load_restart_data(&num_total_fields, sizeof(int));
    EXECUTION_REPORT(REPORT_ERROR, inout_interface->get_comp_id(), num_total_fields == fields_mem_registered.size(), "Error happens when loading the restart data file \"%s\": it does not match the configuration of the interface \"%s\": the original number of total fields of this interface recorded in the restart data file is %d while the current number is %d. Please check.", restart_buffer->get_input_restart_mgt_info_file(), inout_interface->get_interface_name(), num_total_fields, fields_mem_registered.size());

    if (inout_interface->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT) {
        for (int i = 0; i < runtime_inter_averaging_algorithm.size(); i ++) {
            restart_buffer->load_restart_data(&temp_int, sizeof(int));
            if (temp_int == 0) {
                EXECUTION_REPORT(REPORT_ERROR, -1, runtime_inter_averaging_algorithm[i] == NULL, "Software error1 in Connection_coupling_procedure::import_restart_data");
                continue;
            }            
            EXECUTION_REPORT(REPORT_ERROR, -1, runtime_inter_averaging_algorithm[i] != NULL, "Software error2 in Connection_coupling_procedure::import_restart_data");
            runtime_inter_averaging_algorithm[i]->restart_read(restart_buffer, "aver_inter");
        }
        for (int i = 0; i < runtime_inner_averaging_algorithm.size(); i ++) {
            restart_buffer->load_restart_data(&temp_int, sizeof(int));
            if (temp_int == 0) {                
                EXECUTION_REPORT(REPORT_ERROR, -1, runtime_inner_averaging_algorithm[i] == NULL, "Software error3 in Connection_coupling_procedure::import_restart_data");
                continue;        
            }
            EXECUTION_REPORT(REPORT_ERROR, -1, runtime_inner_averaging_algorithm[i] != NULL, "Software error4 in Connection_coupling_procedure::import_restart_data");
            runtime_inner_averaging_algorithm[i]->restart_read(restart_buffer, "aver_inner");
        }
    }

    for (i = 0; i < num_total_fields; i ++) {
        restart_buffer->load_restart_string(restart_field_name, str_size, NAME_STR_SIZE);
        for (j = 0; j < fields_mem_registered.size(); j ++)
            if (words_are_the_same(restart_field_name, fields_mem_registered[j]->get_field_name()))
                break;
        EXECUTION_REPORT(REPORT_ERROR, inout_interface->get_comp_id(), i == j, "Error happens when loading the restart data file \"%s\": it does not match the configuration of the interface \"%s\": its original %th field recorded in the restart data file is \"%s\" while the current %th field is \"%s\". Please check.", restart_buffer->get_input_restart_mgt_info_file(), inout_interface->get_interface_name(), i, restart_field_name, fields_mem_registered[i]->get_field_name());
    }
    restart_buffer->load_restart_data(&current_remote_fields_elapsed_time, sizeof(long));
    restart_buffer->load_restart_data(&last_remote_fields_time, sizeof(long));
    fields_time_info_dst->import_restart_data(restart_buffer);
    fields_time_info_src->import_restart_data(restart_buffer);
}


Inout_interface::Inout_interface(const char *temp_array_buffer, long &buffer_content_iter)
{
    int num_interfaces;


    interface_id = 0;
    read_data_from_array_buffer(interface_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(comp_full_name, NAME_STR_SIZE, temp_array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(&interface_type, sizeof(int), temp_array_buffer, buffer_content_iter, true);
    read_data_from_array_buffer(&num_interfaces, sizeof(int), temp_array_buffer, buffer_content_iter, true);
    for (int i = 0; i < num_interfaces; i ++) {
        fields_name.push_back(strdup(temp_array_buffer+buffer_content_iter-NAME_STR_SIZE));
        buffer_content_iter -= NAME_STR_SIZE;
    }
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_full_name);
    if (comp_node != NULL)
        comp_id = comp_node->get_local_node_id();
    else comp_id = -1;

    restart_mgr = NULL;
    inversed_dst_fraction = NULL;
}


Inout_interface::Inout_interface(const char *interface_name, int interface_id, int num_fields, int *field_ids_src, int *field_ids_dst, int timer_id, int inst_or_aver, int array_size_src, int array_size_dst, const char *API_label, const char *annotation)
{
    char child_interface_name[NAME_STR_SIZE];


    sprintf(child_interface_name, "%s_child_export", interface_name);
    children_interfaces.push_back(new Inout_interface(child_interface_name, -1, 1, num_fields, field_ids_src, array_size_src, timer_id, inst_or_aver, "field_instance_IDs_source", annotation, API_ID_INTERFACE_REG_NORMAL_REMAP, INTERFACE_SOURCE_REGISTER, true));
    sprintf(child_interface_name, "%s_child_import", interface_name);
    children_interfaces.push_back(new Inout_interface(child_interface_name, -1, 0, num_fields, field_ids_dst, array_size_dst, timer_id, inst_or_aver, "field_instance_IDs_target", annotation, API_ID_INTERFACE_REG_NORMAL_REMAP, INTERFACE_SOURCE_REGISTER, true));
    initialize_data(interface_name, interface_id, 2, timer_id, inst_or_aver, field_ids_src, INTERFACE_SOURCE_REGISTER, annotation);
    this->timer->reset_remote_lag_count();
    children_interfaces[0]->timer->reset_remote_lag_count();
    children_interfaces[1]->timer->reset_remote_lag_count();

    for (int i = 0; i < num_fields; i ++) {    
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(memory_manager->get_field_instance(field_ids_src[i])->get_field_name(),memory_manager->get_field_instance(field_ids_dst[i])->get_field_name()), "Error happens when calling the API \"%s\" to register an interface named \"%s\": the number %d field instance (\"%s\") specified by the parameter \"field_instance_IDs_source\" are not consistent with the number %d field instance (\"%s\") specified by the parameter \"field_instance_IDs_target\" (a source field instance must have the same field name with the same number of target field instance). Please check the model code with the annotation \"%s\".", API_label, interface_name, i+1, memory_manager->get_field_instance(field_ids_src[i])->get_field_name(), i+1, memory_manager->get_field_instance(field_ids_dst[i])->get_field_name(), annotation);        
        EXECUTION_REPORT(REPORT_ERROR, comp_id, memory_manager->get_field_instance(field_ids_src[i])->get_data_buf() != memory_manager->get_field_instance(field_ids_dst[i])->get_data_buf(), "Error happens when calling the API \"%s\" to register an interface named \"%s\": the number %d field instance (\"%s\") specified by the parameter \"field_instance_IDs_source\" have the same model data buffer with the number %d field instance (\"%s\") specified by the parameter \"field_instance_IDs_target\" (the source and target field instances may be the same), which is not allowed. Please check the model code with the annotation \"%s\".", API_label, interface_name, i+1, memory_manager->get_field_instance(field_ids_src[i])->get_field_name(), i+1, memory_manager->get_field_instance(field_ids_dst[i])->get_field_name(), annotation);
    }
}


Inout_interface::Inout_interface(const char *interface_name, int interface_id, int interface_type, int num_fields, int *field_ids, int array_size, int timer_id, int inst_or_aver, const char *field_ids_parameter_name, const char *annotation, int API_id, int interface_source, bool is_child_interface)
{
    char API_label[NAME_STR_SIZE];


    get_API_hint(-1, API_id, API_label);
    
    common_checking_for_interface_registration(num_fields, field_ids, array_size, timer_id, inst_or_aver, interface_type, interface_name, API_id, interface_source, field_ids_parameter_name, annotation);
    initialize_data(interface_name, interface_id, interface_type, timer_id, inst_or_aver, field_ids, interface_source, annotation);
    this->is_child_interface = is_child_interface;    
    for (int i = 0; i < num_fields; i ++) {
        Field_mem_info *field_instance = memory_manager->get_field_instance(field_ids[i]);
        if (!is_child_interface)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, field_instance->is_CPL_field_inst(), "Error happens when calling the API \"%s\" to register an interface named \"%s\" at the model code with the annotation \"%s\": the field instance of \"%s\" cannot not be referred by an import/export interface because it has not been declared as a coupling field instance. Please check the parameter \"usage_tag\" when registering this field instance (at the model code with the annotation \"%s\")", API_label, interface_name, annotation, field_instance->get_field_name(), annotation_mgr->get_annotation(field_instance->get_field_instance_id(), "allocate field instance"));
        fields_mem_registered.push_back(field_instance);
        fields_connected_status.push_back(false);
		fields_coupling_procedures.push_back(NULL);
        if (interface_type == COUPLING_INTERFACE_MARK_IMPORT && !is_child_interface)
            restart_mgr->add_restarted_field_instance(fields_mem_registered[fields_mem_registered.size()-1], true);
    }
    fields_connected_status.push_back(false);
    num_fields_connected = 0;
}


void Inout_interface::initialize_data(const char *interface_name, int interface_id, int interface_type, int timer_id, int inst_or_aver, int *field_ids, int interface_source, const char *annotation)
{
    this->interface_id = interface_id;
    this->interface_type = interface_type;
    this->execution_checking_status = 0;
    this->last_execution_time = -1;
    this->is_child_interface = false;
    Coupling_timer *existing_timer = timer_mgr->get_timer(timer_id);
    this->timer = new Coupling_timer(existing_timer->get_comp_id(), -1, existing_timer);
    timer_mgr->add_timer(this->timer);
    this->comp_id = this->timer->get_comp_id();
    this->interface_source = interface_source;
    this->inversed_dst_fraction = NULL;
    strcpy(this->interface_name, interface_name);
    strcpy(this->comp_full_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"in Inout_interface::initialize_data")->get_full_name());
    this->inst_or_aver = inst_or_aver;
    annotation_mgr->add_annotation(interface_id, "registering interface", annotation);
    time_mgr = components_time_mgrs->get_time_mgr(comp_id);
    this->bypass_counter = 0;
    this->mgt_info_has_been_restarted = false;
    restart_mgr = comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_restart_mgr();
}


Inout_interface::~Inout_interface()
{
    if (inversed_dst_fraction != NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "inout interface %s %lx release %lx", interface_name, this, inversed_dst_fraction);
        delete [] inversed_dst_fraction;
    }

    for (int i = 0; i < fields_name.size(); i ++)
        delete [] fields_name[i];

    for (int i = 0; i < coupling_procedures.size(); i ++)
        delete coupling_procedures[i];

    for (int i = 0; i < children_interfaces.size(); i ++)
        delete children_interfaces[i];
}


void Inout_interface::common_checking_for_interface_registration(int num_fields, int *field_ids, int array_size, int timer_id, int inst_or_aver, int interface_type, const char *interface_name, int API_id, int interface_source, const char *field_ids_parameter_name, const char *annotation)
{
    int comp_id = -1;
    char str[NAME_STR_SIZE], API_label[NAME_STR_SIZE];


    get_API_hint(-1, API_id, API_label);
    
    EXECUTION_REPORT(REPORT_ERROR, -1, num_fields > 0 && num_fields <= 1000, "Error happens when calling the API \"%s\" to register an interface named \"%s\": the parameter \"num_field_instances\" (currently is %d) cannot be smaller than 1 or larger than the maximum number (1000). Please verify the model code with the annotation \"%s\".", API_label, interface_name, num_fields, annotation);
    EXECUTION_REPORT(REPORT_ERROR, -1, num_fields <= array_size, "Error happens when calling the API \"%s\" to register an interface named \"%s\": the array size (currently is %d) of parameter \"%s\" cannot be smaller than the parameter \"num_field_instances\" (currently is %d). Please verify the model code with the annotation \"%s\".", API_label, interface_name, num_fields, field_ids_parameter_name, array_size, annotation);
    for (int i = 0; i < num_fields; i ++) {
        if (interface_source != INTERFACE_SOURCE_IO_WRITE)
            EXECUTION_REPORT(REPORT_ERROR, -1, memory_manager->check_is_legal_field_instance_id(field_ids[i]) && memory_manager->get_field_instance(field_ids[i])->get_is_registered_model_buf(), "Error happens when calling the API \"%s\" to register an interface named \"%s\": the parameter \"%s\" contains wrong field instance ID (the %dth element of the array is wrong). Please verify the model code related to the annotation \"%s\"", API_label, interface_name, field_ids_parameter_name, i+1, annotation);
        if (i == 0)
            comp_id = memory_manager->get_field_instance(field_ids[i])->get_comp_id();
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, comp_id == memory_manager->get_field_instance(field_ids[i])->get_comp_id(), "Error happens when calling the API \"%s\" to register an interface named \"%s\": the field instances specified via the parameter \"%s\" should but not correspond to the same component model currently: the first field instance corresponds to the component model \"%s\" while the %dth field instance corresponds to the component model \"%s\". Please verify the model code with the annotation \"%s\".", API_label, interface_name, field_ids_parameter_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false, "")->get_comp_full_name(), i+1, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(memory_manager->get_field_instance(field_ids[i])->get_comp_id(),false, "")->get_comp_full_name(), annotation);
    }
    EXECUTION_REPORT(REPORT_ERROR, comp_id, timer_mgr->check_is_legal_timer_id(timer_id), "Error happens when calling the API \"%s\" to register an interface named \"%s\": the parameter \"timer_ID\" (currently is 0x%x) is not the legal ID of a timer. Please verify the model code related to the annotation \"%s\"", API_label, interface_name, timer_id, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, comp_id == timer_mgr->get_timer(timer_id)->get_comp_id(), "Error happens when calling the API \"%s\" to register an interface named \"%s\": the parameter \"timer_ID\" and the parameter \"%s\" do not correspond to the same component model (the parameter \"timer_ID\" corresponds to the component model \"%s\" while \"%s\" corresponds to the component model \"%s\"). Please verify the model code related to the annotation \"%s\"", API_label, interface_name, field_ids_parameter_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(timer_mgr->get_timer(timer_id)->get_comp_id(),false, "")->get_comp_full_name(), field_ids_parameter_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false, "")->get_comp_full_name(), annotation);
    if (interface_source != INTERFACE_SOURCE_IO_WRITE && interface_type == COUPLING_INTERFACE_MARK_EXPORT)
        for (int i = 0; i < num_fields; i ++) 
            for (int j = i+1; j < num_fields; j ++)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, !words_are_the_same(memory_manager->get_field_instance(field_ids[i])->get_field_name(),memory_manager->get_field_instance(field_ids[j])->get_field_name()), "Error happens when calling the API \"%s\" to register an interface named \"%s\": the parameter \"%s\" is not allowed to include more than one instance of the same field (field \"%s\"). Please verify the model code related to the annotation \"%s\"", API_label, interface_name, field_ids_parameter_name, memory_manager->get_field_instance(field_ids[i])->get_field_name(), annotation);            

    sprintf(str, "registering an interface named \"%s\"", interface_name);
    synchronize_comp_processes_for_API(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_interface::Inout_interface"), str, annotation);
    
    if (interface_source == INTERFACE_SOURCE_REGISTER)
        comp_comm_group_mgt_mgr->confirm_coupling_configuration_active(comp_id, API_id, true, annotation);    
    check_and_verify_name_format_of_string_for_API(comp_id, interface_name, API_id, "the interface", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"executing an import/export interface"), str, interface_name, "the corresponding interface name", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"in Inout_interface::Inout_interface"), NULL, num_fields, "num_field_instances", annotation);
    sprintf(str, "\"%s\" (the information of the field instances)", field_ids_parameter_name);
    for (int i = 0; i < num_fields; i ++)
        check_API_parameter_field_instance(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_interface::Inout_interface"), "registering an interface", field_ids[i], str, annotation);
    check_API_parameter_timer(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_interface::Inout_interface"), "registering an interface", timer_id, "timer_ID (the information of the timer)", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_interface::Inout_interface"), "registering an interface", inst_or_aver, "inst_or_aver (the tag for using instantaneous or time averaged field value)", annotation);
}


void Inout_interface::report_common_field_instances(const Inout_interface *another_interface)
{
    if (this->interface_source != INTERFACE_SOURCE_REGISTER || another_interface->interface_source != INTERFACE_SOURCE_REGISTER)
        return;

    if (this->interface_type == COUPLING_INTERFACE_MARK_IMPORT && another_interface->interface_type == COUPLING_INTERFACE_MARK_IMPORT) {
        for (int i = 0; i < this->fields_mem_registered.size(); i ++)
            for (int j = 0; j < another_interface->fields_mem_registered.size(); j ++)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, this->fields_mem_registered[i]->get_data_buf() != another_interface->fields_mem_registered[j]->get_data_buf(), "Two import interfaces (\"%s\" and \"%s\") share the same data buffer of the field (field name is \"%s\") which is not allowed. Please check the model code with the annotation \"%s\" and \"%s\".",
                                 this->interface_name, another_interface->interface_name, fields_mem_registered[i]->get_field_name(), annotation_mgr->get_annotation(this->interface_id, "registering interface"), annotation_mgr->get_annotation(another_interface->interface_id, "registering interface"));        
        return;
    }

    if (this->interface_type != COUPLING_INTERFACE_MARK_EXPORT || another_interface->interface_type != COUPLING_INTERFACE_MARK_EXPORT)
        return;

    for (int i = 0; i < this->fields_mem_registered.size(); i ++)
        for (int j = 0; j < another_interface->fields_mem_registered.size(); j ++)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, !words_are_the_same(this->fields_mem_registered[i]->get_field_name(), another_interface->fields_mem_registered[j]->get_field_name()), "Two export interfaces (\"%s\" and \"%s\") provide the same field (field name is \"%s\") which is not allowed. Please check the model code with the annotation \"%s\" and \"%s\".",
                             this->interface_name, another_interface->interface_name, fields_mem_registered[i]->get_field_name(), annotation_mgr->get_annotation(this->interface_id, "registering interface"), annotation_mgr->get_annotation(another_interface->interface_id, "registering interface"));
}


void Inout_interface::get_fields_name(std::vector<const char*> *fields_name)
{
    if (this->fields_mem_registered.size() > 0) {
        for (int i = 0; i < this->fields_mem_registered.size(); i ++)
            fields_name->push_back(this->fields_mem_registered[i]->get_field_name());
    }
    else {
        for (int i = 0; i < this->fields_name.size(); i ++)
            fields_name->push_back(this->fields_name[i]);
    }
}


const char *Inout_interface::get_field_name(int number)
{
    if (number >= fields_name.size() && number >= fields_mem_registered.size())
        return NULL;

    if (number < fields_name.size())
        return fields_name[number];

    return fields_mem_registered[number]->get_field_name();
}


int Inout_interface::get_num_dst_fields()
{
    if (interface_type == COUPLING_INTERFACE_MARK_IMPORT)
        return fields_mem_registered.size();
    if (interface_type == COUPLING_INTERFACE_MARK_EXPORT)
        return 0;
    if (interface_type == COUPLING_INTERFACE_MARK_NORMAL_REMAP)
        return children_interfaces[1]->fields_mem_registered.size();
    if (interface_type == COUPLING_INTERFACE_MARK_FRAC_REMAP)
        return children_interfaces[1]->fields_mem_registered.size() - 1;
    
    return -1;
}


Field_mem_info *Inout_interface::search_registered_field_instance(const char *field_name, int &field_local_index)
{
    field_local_index = -1;
    for (int i = 0; i < fields_mem_registered.size(); i ++)
        if (words_are_the_same(fields_mem_registered[i]->get_field_name(), field_name)) {
            field_local_index = i;
            return fields_mem_registered[i];
        }

    return NULL;
}


void Inout_interface::transform_interface_into_array(char **temp_array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    int temp_int = 0;


    if (interface_type == COUPLING_INTERFACE_MARK_IMPORT && num_fields_connected == fields_mem_registered.size())
        return;
    
    for (int i = fields_mem_registered.size()-1; i >= 0 ; i --) {
        if (interface_type == COUPLING_INTERFACE_MARK_IMPORT && fields_connected_status[i])
            continue;
        write_data_into_array_buffer(fields_mem_registered[i]->get_field_name(), NAME_STR_SIZE, temp_array_buffer, buffer_max_size, buffer_content_size);
        temp_int ++;
    }
    write_data_into_array_buffer(&temp_int, sizeof(int), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&interface_type, sizeof(int), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(comp_full_name, NAME_STR_SIZE, temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(interface_name, NAME_STR_SIZE, temp_array_buffer, buffer_max_size, buffer_content_size);
}


void Inout_interface::write_restart_mgt_info(Restart_buffer_container *restart_buffer)
{
    if (restart_buffer == NULL)
        restart_buffer = restart_mgr->apply_restart_buffer(comp_full_name, RESTART_BUF_TYPE_INTERFACE, interface_name);

    for (int i = coupling_procedures.size() - 1; i >= 0; i --)
        coupling_procedures[i]->write_restart_mgt_info(restart_buffer);
    int temp_int = coupling_procedures.size();
    restart_buffer->dump_in_data(&temp_int, sizeof(int));
    for (int i = children_interfaces.size()-1; i >= 0; i --)
        children_interfaces[i]->write_restart_mgt_info(restart_buffer);
    temp_int = children_interfaces.size();
    restart_buffer->dump_in_data(&temp_int, sizeof(int));
    restart_buffer->dump_in_data(&last_execution_time, sizeof(long));
    timer->write_timer_into_array(restart_buffer->get_buffer_content_ptr(), *(restart_buffer->get_buffer_max_size_ptr()), *(restart_buffer->get_buffer_content_iter_ptr()));
}


void Inout_interface::import_restart_data(Restart_buffer_container *restart_buffer)
{
    int num_children, num_procedures;
    bool successful;

    if (restart_buffer == NULL)
        restart_buffer = restart_mgr->search_restart_buffer(RESTART_BUF_TYPE_INTERFACE, interface_name); 
    EXECUTION_REPORT(REPORT_ERROR, restart_mgr->get_comp_id(), restart_buffer != NULL, "Error happens when loading the restart data file \"%s\" at the model code with the annotation \"%s\": this file does not include the data for restarting the interface \"%s\"", restart_mgr->get_input_restart_mgt_info_file(), restart_mgr->get_restart_read_annotation(), interface_name);
    Coupling_timer *restart_timer = new Coupling_timer(restart_buffer->get_buffer_content(), *(restart_buffer->get_buffer_content_iter_ptr()), comp_id, false, successful);
    EXECUTION_REPORT(REPORT_ERROR, -1, successful, "Fail to load the restart data file \"%s\": its format is wrong, or the information it includes is not complete. Please try a different restart time with complete restart data files.", restart_mgr->get_input_restart_mgt_info_file());
    EXECUTION_REPORT(REPORT_ERROR, comp_id, restart_timer->is_the_same_with(timer), "Error happens when loading the restart data file \"%s\": the timer of the interface \"%s\" in the restart data file is different from the current timer speicifed by the model code. Please verify.", restart_mgr->get_input_restart_mgt_info_file(), interface_name);
    restart_buffer->load_restart_data(&last_execution_time, sizeof(long));
    restart_buffer->load_restart_data(&num_children, sizeof(int));
    EXECUTION_REPORT(REPORT_ERROR, get_comp_id(), num_children == children_interfaces.size(), "Error happens when loading the restart data file \"%s\": it does not match the configuration of the interface \"%s\". Please check.", restart_mgr->get_input_restart_mgt_info_file(), get_interface_name());
    for (int i = 0; i < num_children; i ++)
        children_interfaces[i]->import_restart_data(restart_buffer);
    restart_buffer->load_restart_data(&num_procedures, sizeof(int));
    EXECUTION_REPORT(REPORT_ERROR, get_comp_id(), num_procedures == coupling_procedures.size(), "Error happens when loading the restart data file \"%s\": it does not match the configuration of the interface \"%s\". Please check.", restart_mgr->get_input_restart_mgt_info_file(), get_interface_name());
    for (int i = 0; i < num_procedures; i ++)
        coupling_procedures[i]->import_restart_data(restart_buffer);
}


void Inout_interface::read_restart_fields(int API_id, const char *annotation)
{
    char API_label[NAME_STR_SIZE];


    if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_INITIAL)
        return;

    get_API_hint(comp_id, API_id, API_label);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, restart_mgr->get_restart_read_data_file_name() != NULL, "Error happens when calling the API \"%s\" to read restart fields: the API \"CCPL_start_restart_read_IO\" has not been called before. Please verify the model code corresponding to the annotation %s", API_label, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, !time_mgr->get_time_has_been_advanced(), "Error happens when calling the API \"%s\" to read restart fields: the model time has already been advanced before. Please verify the model code corresponding to the annotation %s", API_label, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, interface_type == COUPLING_INTERFACE_MARK_IMPORT, "Error happens when calling the API \"CCPL_restart_read_fields_interface\": the corresponding coupling interface \"%s\" is not an import interface (this API only reads restart fields for import interfaces). Please verify the model code with the annotation \"%s\"", interface_name, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, (execution_checking_status & 0x2) == 0, "Error happens when calling the API \"CCPL_restart_read_fields_interface\": the corresponding import interface \"%s\" has been executed without bypassing the timer. Please verify the model code with the annotation \"%s\"", interface_name, annotation);
    synchronize_comp_processes_for_API(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, ""), "read restart fields for the given coupling interface", annotation);
    for (int i = 0; i < fields_mem_registered.size(); i ++)
        restart_mgr->read_restart_field_data(fields_mem_registered[i], NULL, NULL, false, NULL, restart_mgr->get_bypass_import_fields_at_read(), annotation);
}


void Inout_interface::add_coupling_procedure(Connection_coupling_procedure *coupling_procedure)
{
    coupling_procedures.push_back(coupling_procedure);
    if (interface_type == COUPLING_INTERFACE_MARK_IMPORT || interface_type == COUPLING_INTERFACE_MARK_EXPORT) {
        EXECUTION_REPORT(REPORT_ERROR, -1, fields_connected_status.size() > 0, "Software error in Inout_interface::add_coupling_procedure: %s", interface_name);
        for (int i = 0; i < coupling_procedure->fields_mem_registered.size(); i ++)
            for (int j = 0; j < fields_mem_registered.size(); j ++)
                if (coupling_procedure->fields_mem_registered[i] == fields_mem_registered[j]) {
                    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Add coupling procedures to the field \"%s\" of import/export interface \"%s", fields_mem_registered[j]->get_field_name(), interface_name);
                    if (interface_type == COUPLING_INTERFACE_MARK_IMPORT)
                        EXECUTION_REPORT(REPORT_ERROR, -1, !fields_connected_status[j], "Software error in Inout_interface::add_coupling_procedure: %s %s %d", interface_name, fields_mem_registered[j]->get_field_name(), j);
                    if (!fields_connected_status[i])
                        num_fields_connected ++;
                    fields_connected_status[j] = true;
					fields_coupling_procedures[j] = coupling_procedure;
                }
    }
}


void Inout_interface::preprocessing_for_frac_based_remapping()
{
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Pre-process the fraction");

    for (int i = 0; i < fields_mem_registered.size()-1; i ++) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, fields_mem_registered[i] != children_interfaces[0]->fields_mem_registered[i], "Software error1 in Inout_interface::preprocessing_for_frac_based_remapping");        
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, fields_mem_registered[i] != children_interfaces[0]->fields_mem_registered[i], "Software error2 in Inout_interface::preprocessing_for_frac_based_remapping");
    }
    
    for (int i = 0; i < fields_mem_registered.size()-1; i ++) {
        if (words_are_the_same(fields_mem_registered[i]->get_data_type(), DATA_TYPE_FLOAT))
            if (words_are_the_same(fields_mem_registered[fields_mem_registered.size()-1]->get_data_type(), DATA_TYPE_FLOAT))
                arrays_multiplication_template((float*)fields_mem_registered[i]->get_data_buf(), (float*)fields_mem_registered[fields_mem_registered.size()-1]->get_data_buf(), (float*)children_interfaces[0]->fields_mem_registered[i]->get_data_buf(), fields_mem_registered[i]->get_size_of_field());
            else arrays_multiplication_template((float*)fields_mem_registered[i]->get_data_buf(), (double*)fields_mem_registered[fields_mem_registered.size()-1]->get_data_buf(), (float*)children_interfaces[0]->fields_mem_registered[i]->get_data_buf(), fields_mem_registered[i]->get_size_of_field());
        else if (words_are_the_same(fields_mem_registered[fields_mem_registered.size()-1]->get_data_type(), DATA_TYPE_FLOAT))
            arrays_multiplication_template((double*)fields_mem_registered[i]->get_data_buf(), (float*)fields_mem_registered[fields_mem_registered.size()-1]->get_data_buf(), (double*)children_interfaces[0]->fields_mem_registered[i]->get_data_buf(), fields_mem_registered[i]->get_size_of_field());
        else arrays_multiplication_template((double*)fields_mem_registered[i]->get_data_buf(), (double*)fields_mem_registered[fields_mem_registered.size()-1]->get_data_buf(), (double*)children_interfaces[0]->fields_mem_registered[i]->get_data_buf(), fields_mem_registered[i]->get_size_of_field());
        children_interfaces[0]->fields_mem_registered[i]->define_field_values(false);
    }    
}


void Inout_interface::postprocessing_for_frac_based_remapping(bool bypass_timer)
{
    Field_mem_info *dst_value_field, *dst_frac_field;


    if (children_interfaces[1]->coupling_procedures[0]->get_runtime_remap_algorithm(0) == NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Does not post-process the fraction becase the runtime alogrithm is NULL");        
        return;
    }
    
    if (!timer->is_timer_on() && !bypass_timer)
        return;

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Post-process the fraction");     

    dst_frac_field = children_interfaces[1]->fields_mem_registered[fields_mem_registered.size()-1];

    if (words_are_the_same(dst_frac_field->get_data_type(), DATA_TYPE_FLOAT)) {
        float *dst_frac_buf = (float*)dst_frac_field->get_data_buf();
        for (int i = dst_frac_field->get_size_of_field()-1; i >= 0; i --)
            if (dst_frac_buf[i] == (float) 0.0)
                ((float*) inversed_dst_fraction)[i] = dst_frac_buf[i];
            else ((float*) inversed_dst_fraction)[i] = ((float)1.0) / dst_frac_buf[i];    
    }
    else {
        double *dst_frac_buf = (double*)dst_frac_field->get_data_buf();
        for (int i = dst_frac_field->get_size_of_field()-1; i >= 0; i --)
            if (dst_frac_buf[i] == (double) 0.0)
                ((double*) inversed_dst_fraction)[i] = dst_frac_buf[i];
            else ((double*) inversed_dst_fraction)[i] = ((double)1.0) / dst_frac_buf[i];    
    }
    
    for (int i = 0; i < fields_mem_registered.size()-1; i ++) {
        dst_value_field = children_interfaces[1]->fields_mem_registered[i];
        if (words_are_the_same(dst_value_field->get_data_type(), DATA_TYPE_FLOAT))
            if (words_are_the_same(dst_frac_field->get_data_type(), DATA_TYPE_FLOAT))
                arrays_multiplication_template((float*)dst_value_field->get_data_buf(), (float*)inversed_dst_fraction, (float*)dst_value_field->get_data_buf(), dst_frac_field->get_size_of_field());
            else arrays_multiplication_template((float*)dst_value_field->get_data_buf(), (double*)inversed_dst_fraction, (float*)dst_value_field->get_data_buf(), dst_frac_field->get_size_of_field());
        else if (words_are_the_same(dst_frac_field->get_data_type(), DATA_TYPE_FLOAT))
            arrays_multiplication_template((double*)dst_value_field->get_data_buf(), (float*)inversed_dst_fraction, (double*)dst_value_field->get_data_buf(), dst_frac_field->get_size_of_field());
        else arrays_multiplication_template((double*)dst_value_field->get_data_buf(), (double*)inversed_dst_fraction, (double*)dst_value_field->get_data_buf(), dst_frac_field->get_size_of_field());
    }    
}


void Inout_interface::execute(bool bypass_timer, int API_id, int *field_update_status, int size_field_update_status, const char *annotation)
{
    bool at_first_normal_step = false;


    if (interface_type == COUPLING_INTERFACE_MARK_IMPORT) {
        if (fields_mem_registered.size() != num_fields_connected)
            for (int i = 0; i < fields_mem_registered.size(); i ++)
                if (!fields_connected_status[i] && (imported_fields_necessity.size() == 0 || imported_fields_necessity[i] == FIELD_NECESSITY_NECESSARY)) {
                    std::vector<const char *> export_comp_full_names, export_interface_names;
                    inout_interface_mgr->get_all_export_interfaces_of_a_field(comp_id, fields_mem_registered[i]->get_field_name(), export_comp_full_names, export_interface_names);
                    char *error_string = NULL;
                    long string_size;
                    if (comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_current_proc_local_id() == 0) {
                        error_string = new char [NAME_STR_SIZE*(export_comp_full_names.size()+2)];
                        error_string[0] = '\0';
                        for (int i = 0; i < export_comp_full_names.size(); i ++)
                            sprintf(error_string, "%s                   %d) Component model is \"%s\", export interface is \"%s\"\n", error_string, i+1, export_comp_full_names[i], export_interface_names[i]);
                        string_size = strlen(error_string) + 1;
                    }
                    bcast_array_in_one_comp(comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_current_proc_local_id(), &error_string, string_size, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comm_group());
                    EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "ERROR happens when executing the import interface \"%s\" (the corresponding code annotation is \"%s\"): the coupling procedures for the import field \"%s\" have not been fully generated. The export interfaces that have been registered with this field are listed as follows. Please make sure the correct coupling generation.\n%s", interface_name, annotation, fields_mem_registered[i]->get_field_name(), error_string);
                    if (error_string != NULL)
                        delete [] error_string;
                }    
        if (fields_mem_registered.size() > size_field_update_status)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Fail to execute the interface \"%s\" corresponding to the model code with the annotation \"%s\": the array size of \"field_update_status\" (%d) is smaller than the number of fields (%d). Please verify.", interface_name, annotation, size_field_update_status, fields_mem_registered.size());
        for (int i = 0; i < fields_mem_registered.size(); i ++)
            field_update_status[i] = 0;
    }    
    else if (interface_type == COUPLING_INTERFACE_MARK_NORMAL_REMAP) {
        if (((int)children_interfaces[0]->fields_mem_registered.size()) > size_field_update_status)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Fail to execute the interface \"%s\" corresponding to the model code with the annotation \"%s\": the array size of \"field_update_status\" (%d) is smaller than the number of fields (%d). Please verify.", interface_name, annotation, size_field_update_status, children_interfaces[0]->fields_mem_registered.size());
        for (int i = 0; i < ((int)children_interfaces[0]->fields_mem_registered.size()); i ++)
            field_update_status[i] = 0;
    }    
    else if (interface_type == COUPLING_INTERFACE_MARK_FRAC_REMAP) {
        if (((int)children_interfaces[0]->fields_mem_registered.size())-1 > size_field_update_status)
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Fail to execute the interface \"%s\" corresponding to the model code with the annotation \"%s\": the array size of \"field_update_status\" (%d) is smaller than the number of fields (%d). Please verify.", interface_name, annotation, size_field_update_status, children_interfaces[0]->fields_mem_registered.size()-1);
        for (int i = 0; i < ((int)children_interfaces[0]->fields_mem_registered.size())-1; i ++)
            field_update_status[i] = 0;
    }

    if (!is_child_interface && !bypass_timer && !mgt_info_has_been_restarted && (time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE || time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH)) {
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Import restart data for the interface \"%s\"\n", interface_name);
        import_restart_data(NULL);
        mgt_info_has_been_restarted = true;
    }
    
    if (bypass_timer)
        bypass_counter ++;
    if (time_mgr->check_is_model_run_finished()) {
        EXECUTION_REPORT(REPORT_WARNING, comp_id, false, "The import/export interface \"%s\" (corresponding to the model code annotation \"%s\") will not execute at time %08d-%05d because the model run has finished",
                         interface_name, annotation_mgr->get_annotation(interface_id, "registering interface"), time_mgr->get_current_date(), time_mgr->get_current_second());
        return;
    }

    if (bypass_timer && (execution_checking_status & 0x2) != 0)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when executing the interface \"%s\" at the model code with the corresponding annotation \"%s\": the timer of this interface cannot be bypassed again because this interface has been executed with the timer unbypassed before (at the model code with the corresponding annotation \"%s\"). Please check and verify.", interface_name, annotation, annotation_mgr->get_annotation(interface_id, "using timer"));
    if ((execution_checking_status & 0x1) == 0 && bypass_timer || (execution_checking_status & 0x2) == 0 && !bypass_timer) {
        synchronize_comp_processes_for_API(comp_id, API_id, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"software error")->get_comm_group(), "executing an import/export interface", annotation);
        check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"executing an import/export interface"), "executing an import/export interface", interface_name, "the corresponding interface name", annotation);
        int bypass_timer_int;
        if (bypass_timer)
            bypass_timer_int = 0;
        else bypass_timer_int = 1;
        check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"executing an import/export interface"), NULL, bypass_timer_int, "the value for specifying whether bypass timers", annotation);
        if (bypass_timer) {
            execution_checking_status = execution_checking_status | 0x1;
            annotation_mgr->add_annotation(interface_id, "bypassing timer", annotation);
        }
        else {
            at_first_normal_step = (execution_checking_status & 0x2) == 0;
            execution_checking_status = execution_checking_status | 0x2;
            annotation_mgr->add_annotation(interface_id, "using timer", annotation);
        }
    }

    long current_execution_time = ((long)time_mgr->get_current_num_elapsed_day())*100000 + components_time_mgrs->get_time_mgr(comp_id)->get_current_second();
    if (current_execution_time == last_execution_time && !bypass_timer && !at_first_normal_step) {
        int current_year, current_month, current_day, current_second;
        components_time_mgrs->get_time_mgr(comp_id)->get_current_time(current_year, current_month, current_day, current_second, 0, "CCPL internal");
        EXECUTION_REPORT(REPORT_WARNING, comp_id, false, "The import/export interface \"%s\", which is called at the model code with the annotation \"%s\", will not be executed again at the time step %04d-%02d-%02d-%05d, because it has been executed at the same time step before.",
                         interface_name, annotation, current_year, current_month, current_day, current_second);
        return;
    }

    last_execution_time = current_execution_time;

    if (interface_type == COUPLING_INTERFACE_MARK_NORMAL_REMAP || interface_type == COUPLING_INTERFACE_MARK_FRAC_REMAP) {
        for (int i = 0; i < fields_mem_registered.size(); i ++)
            fields_mem_registered[i]->check_field_sum("before executing a remap interface");
        if (interface_type == COUPLING_INTERFACE_MARK_FRAC_REMAP)
            preprocessing_for_frac_based_remapping();
        children_interfaces[0]->execute(bypass_timer, API_id, field_update_status, size_field_update_status+1, annotation);
        children_interfaces[1]->execute(bypass_timer, API_id, field_update_status, size_field_update_status+1, annotation);
        if (interface_type == COUPLING_INTERFACE_MARK_FRAC_REMAP)
            postprocessing_for_frac_based_remapping(bypass_timer);
        return;
    }

    if (interface_type == COUPLING_INTERFACE_MARK_EXPORT) {
        for (int i = 0; i < fields_mem_registered.size(); i ++)
            fields_mem_registered[i]->check_field_sum("before executing an export interface");
        EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "after checking all source fields");
    }
    
    for (int i = 0; i < coupling_procedures.size(); i ++)
        coupling_procedures[i]->execute(bypass_timer, field_update_status, annotation);

    if (interface_type == COUPLING_INTERFACE_MARK_EXPORT) {

#ifdef USE_ONE_SIDED_MPI
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_SEND_WAIT, -1, interface_name);
#endif
        bool all_finish = false;
        while (!all_finish) {
            all_finish = true;
            for (int i = 0; i < coupling_procedures.size(); i ++) {
                if (!coupling_procedures[i]->get_finish_status())
                    coupling_procedures[i]->send_fields(bypass_timer);
                all_finish = all_finish && coupling_procedures[i]->get_finish_status();
            }    
        }
#ifdef USE_ONE_SIDED_MPI
        comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMMUNICATION, TIMING_COMMUNICATION_SEND_WAIT, -1, interface_name);
#endif
    }

    if (interface_type == COUPLING_INTERFACE_MARK_IMPORT) {
        for (int i = 0; i < fields_mem_registered.size(); i ++)
            fields_mem_registered[i]->check_field_sum("after executing an export interface");
    }
}


Inout_interface *Inout_interface::get_child_interface(int i)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, i >= 0 && i < children_interfaces.size(), "Software error in Inout_interface::get_child_interface");

    return children_interfaces[i];
}


void Inout_interface::add_remappling_fraction_processing(void *frac_src, void *frac_dst, int size_frac_src, int size_frac_dst, const char *frac_data_type, const char *API_label, const char *annotation)
{
    for (int j = 0; j < 2; j ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, children_interfaces[j]->coupling_procedures.size() == 1 && children_interfaces[j]->coupling_procedures[0]->get_num_runtime_remap_algorithms() == children_interfaces[j]->fields_mem_registered.size(), "Software error in Inout_interface::add_remappling_fraction_processing");
        for (int i = 0; i < children_interfaces[j]->fields_mem_registered.size(); i ++) {
            Original_grid_info *field_grid = original_grid_mgr->get_original_grid(children_interfaces[j]->fields_mem_registered[i]->get_grid_id());
            EXECUTION_REPORT(REPORT_ERROR, comp_id, field_grid != NULL && field_grid->is_H2D_grid(), "Error happens when calling the API \"%s\" to register an interface named \"%s\": field \"%s\" is not on an horizontal grid. Please verify    the model code model with the annotation \"%s\"", API_label, interface_name, children_interfaces[j]->fields_mem_registered[i]->get_field_name(), annotation);
            EXECUTION_REPORT(REPORT_ERROR, comp_id, children_interfaces[j]->coupling_procedures[0]->get_runtime_remap_algorithm(i) == NULL && children_interfaces[j]->coupling_procedures[0]->get_runtime_remap_algorithm(0) == NULL || children_interfaces[j]->coupling_procedures[0]->get_runtime_remap_algorithm(i) != NULL && children_interfaces[j]->coupling_procedures[0]->get_runtime_remap_algorithm(0) != NULL && children_interfaces[j]->coupling_procedures[0]->get_runtime_remap_algorithm(i)->get_runtime_remapping_weights() == children_interfaces[j]->coupling_procedures[0]->get_runtime_remap_algorithm(0)->get_runtime_remapping_weights(), 
                             "Error happens when calling the API \"%s\" to register an interface named \"%s\": The fields (\"%s\" and \"%s\") to be remapped do not share the same remapping weights. Please verify the model code model with the annotation \"%s\".", API_label, interface_name, children_interfaces[j]->fields_mem_registered[0]->get_field_name(), children_interfaces[j]->fields_mem_registered[i]->get_field_name(), annotation);
        }
    }
    Field_mem_info *template_field_src = children_interfaces[0]->fields_mem_registered[0];
    Field_mem_info *template_field_dst = children_interfaces[1]->fields_mem_registered[0];
    EXECUTION_REPORT(REPORT_ERROR, comp_id, template_field_src->get_size_of_field() == size_frac_src, "Error happens when calling the API \"%s\" to register an interface named \"%s\": the array size of the parameter \"frac_src\" is different from the size of each source field instance. Please verify the model code model with the annotation \"%s\".", API_label, interface_name, annotation);
    int has_frac_dst = size_frac_dst == -1? 0 : 1;
    check_API_parameter_int(comp_id, API_ID_INTERFACE_REG_FRAC_REMAP, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"in add_remappling_fraction_processing")->get_comm_group(), "specification (or not)", has_frac_dst, "frac_dst", annotation);
    if (size_frac_dst != -1)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, template_field_dst->get_size_of_field() == size_frac_dst, "Error happens when calling the API \"%s\" to register an interface named \"%s\": the array size of the parameter \"frac_dst\" is different from the size of each target field instance. Please verify the model code model with the annotation \"%s\".", API_label, interface_name, annotation);

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish checking for adding remappling fraction processing for the remapping interface \"%s\"", interface_name);

    Field_mem_info *frac_field_src = memory_manager->alloc_mem("remap_frac", template_field_src->get_decomp_id(), template_field_src->get_comp_or_grid_id(), BUF_MARK_REMAP_FRAC ^ coupling_generator->get_latest_connection_id(), frac_data_type, "unitless", "source fraction for remapping", false);
    frac_field_src->reset_mem_buf(frac_src, true, -1);
    Field_mem_info *frac_field_dst = memory_manager->alloc_mem("remap_frac", template_field_dst->get_decomp_id(), template_field_dst->get_comp_or_grid_id(), BUF_MARK_REMAP_FRAC ^ coupling_generator->get_latest_connection_id(), frac_data_type, "unitless", "target fraction for remapping", false);
    if (size_frac_dst != -1) 
        frac_field_dst->reset_mem_buf(frac_dst, true, -1);
    memset(frac_field_dst->get_data_buf(), 0, frac_field_dst->get_size_of_field()*get_data_type_size(frac_field_dst->get_data_type()));
    interface_type = COUPLING_INTERFACE_MARK_FRAC_REMAP;
    EXECUTION_REPORT(REPORT_ERROR, -1, fields_mem_registered.size() == 0, "Software error in Inout_interface::add_remappling_fraction_processing");
    for (int i = 0; i < children_interfaces[0]->fields_mem_registered.size(); i ++) {
        fields_mem_registered.push_back(children_interfaces[0]->fields_mem_registered[i]);
        children_interfaces[0]->fields_mem_registered[i] = memory_manager->alloc_mem(fields_mem_registered[i], BUF_MARK_REMAP_FRAC, coupling_generator->get_latest_connection_id(), fields_mem_registered[i]->get_data_type(), true);
    }
    fields_mem_registered.push_back(frac_field_src);

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "After field instance allocation for adding remappling fraction processing for the remapping interface \"%s\"", interface_name);
    
    children_interfaces[0]->fields_mem_registered.push_back(frac_field_src);
    children_interfaces[1]->fields_mem_registered.push_back(frac_field_dst);
    delete children_interfaces[0]->coupling_procedures[0];
    delete children_interfaces[1]->coupling_procedures[0];
    for (int i = 0; i < children_interfaces[1]->fields_connected_status.size(); i ++)
        children_interfaces[1]->fields_connected_status[i] = false;
    children_interfaces[1]->num_fields_connected = 0;
    children_interfaces[0]->coupling_procedures.clear();
    children_interfaces[1]->coupling_procedures.clear();
    int num_fields = children_interfaces[0]->fields_mem_registered.size();
    int *field_ids_src = new int [num_fields];
    for (int i = 0; i < num_fields; i ++)
        field_ids_src[i] = children_interfaces[0]->fields_mem_registered[i]->get_field_instance_id();
    
    inout_interface_mgr->generate_remapping_interface_connection(this, num_fields, field_ids_src, true);
    delete [] field_ids_src;

    if (frac_field_dst->get_size_of_field() > 0)
        inversed_dst_fraction = new char [frac_field_dst->get_size_of_field()*get_data_type_size(frac_field_dst->get_data_type())];
}


int Inout_interface::get_h2d_grid_area_in_remapping_weights(const char *interface_name, int field_index, void *output_area_data, int area_array_size, const char *data_type, const char *annotation)
{
    int i, j;
    double *selected_area_array_in_wgts = NULL;
    bool field_has_connection = false;

    
    if (children_interfaces.size() > 0)
        return children_interfaces[0]->get_h2d_grid_area_in_remapping_weights(interface_name, field_index, output_area_data, area_array_size, data_type, annotation);

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, field_index >= 0 && field_index < fields_mem_registered.size(), "ERROR happens when calling the API \"CCPL_get_H2D_grid_area_in_remapping_wgts\" based on the coupling interface \"%s\": the parameter of field index (%d) is out of bounds ([1,%d]). Please verify the model code with the annotation \"%s\"", interface_name, field_index+1, fields_mem_registered.size(), annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, fields_mem_registered[field_index]->get_grid_id() != -1, "ERROR happens when calling the API \"CCPL_get_H2D_grid_area_in_remapping_wgts\" based on the coupling interface \"%s\": the field \"%s\" corresponding to the field index (%d) is not on a grid. Please verify the model code with the annotation \"%s\"", interface_name, fields_mem_registered[field_index]->get_field_name(), field_index+1, fields_mem_registered.size()-1, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->search_grid_info(fields_mem_registered[field_index]->get_grid_id())->is_H2D_grid(), "ERROR happens when calling the API \"CCPL_get_H2D_grid_area_in_remapping_wgts\" based on the coupling interface \"%s\": the field \"%s\" corresponding to the field index (%d) is not on a horizontal grid. Please verify the model code with the annotation \"%s\"", interface_name, fields_mem_registered[field_index]->get_field_name(), field_index+1, annotation);

    for (i = 0; i < coupling_procedures.size(); i ++) {
        for (j = 0; j < coupling_procedures[i]->coupling_connection->fields_name.size(); j ++)
            if (words_are_the_same(fields_mem_registered[field_index]->get_field_name(),coupling_procedures[i]->coupling_connection->fields_name[j])) {
                if (interface_type == COUPLING_INTERFACE_MARK_IMPORT) {
                    if (coupling_procedures[i]->coupling_connection->dst_fields_info[j]->runtime_remapping_weights != NULL && coupling_procedures[i]->coupling_connection->dst_fields_info[j]->runtime_remapping_weights->get_src_H2D_grid_area() != NULL) {
                        selected_area_array_in_wgts = coupling_procedures[i]->coupling_connection->dst_fields_info[j]->runtime_remapping_weights->get_dst_H2D_grid_area();
                        break;
                    }
                }
                else if (interface_type == COUPLING_INTERFACE_MARK_EXPORT) {
                    if (coupling_procedures[i]->coupling_connection->src_fields_info[j]->runtime_remapping_weights != NULL && coupling_procedures[i]->coupling_connection->src_fields_info[j]->runtime_remapping_weights->get_src_H2D_grid_area() != NULL) {
                        selected_area_array_in_wgts = coupling_procedures[i]->coupling_connection->src_fields_info[j]->runtime_remapping_weights->get_src_H2D_grid_area();
                        break;
                    }                    
                }
                else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Inout_interface::get_h2d_grid_area_in_remapping_weights");
                field_has_connection = true;
            }    
        if (selected_area_array_in_wgts != NULL)
            break;
    }

    EXECUTION_REPORT(REPORT_WARNING, comp_id, field_has_connection, "WARNING happens when calling the API \"CCPL_get_H2D_grid_area_in_remapping_wgts\" based on the coupling interface \"%s\": the field \"%s\" corresponding to the field index (%d) has not been used in model coupling. Please verify the model code with the annotation \"%s\"", interface_name, fields_mem_registered[field_index]->get_field_name(), field_index+1, annotation);

    if (selected_area_array_in_wgts == NULL)
        return 0;

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, area_array_size >= fields_mem_registered[field_index]->get_size_of_field(), "ERROR happens when calling the API \"CCPL_get_H2D_grid_area_in_remapping_wgts\" based on the field \"%s\" of the coupling interface \"%s\": the array size (%d) of the parameter \"area_array\" is smaller than required (%d). ", fields_mem_registered[field_index]->get_field_name(), interface_name, area_array_size, fields_mem_registered[field_index]->get_size_of_field());
    Decomp_info *decomp_info = decomps_info_mgr->get_decomp_info(fields_mem_registered[field_index]->get_decomp_id());    
    const int *local_cells_global_index = decomp_info->get_local_cell_global_indx();
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
        float *float_output_area = (float*) output_area_data;
        for (i = 0; i < decomp_info->get_num_local_cells(); i ++) {
            if (local_cells_global_index[i] != CCPL_NULL_INT)
                float_output_area[i] = (float) (selected_area_array_in_wgts[local_cells_global_index[i]]);
            else float_output_area[i] = NULL_COORD_VALUE;
        }
    }
    else if (words_are_the_same(data_type, DATA_TYPE_DOUBLE)) {
        double *double_output_area = (double*) output_area_data;
        for (i = 0; i < decomp_info->get_num_local_cells(); i ++)
            if (local_cells_global_index[i] != CCPL_NULL_INT)
                double_output_area[i] = selected_area_array_in_wgts[local_cells_global_index[i]];
            else double_output_area[i] = NULL_COORD_VALUE;
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Inout_interface::get_h2d_grid_area_in_remapping_weights: wrong data type");

    return 1;
}


void Inout_interface::write_export_info_into_XML_file(TiXmlElement *parent_element)
{
    TiXmlElement *current_element = new TiXmlElement("export_interface");
    parent_element->LinkEndChild(current_element);
    current_element->SetAttribute("interface_name", interface_name);
    for (int i = 0; i < fields_mem_registered.size(); i ++) {
        TiXmlElement *field_element = new TiXmlElement("export_field");
        current_element->LinkEndChild(field_element);
        field_element->SetAttribute("field_name", fields_mem_registered[i]->get_field_name());
    }
}


void Inout_interface::set_fields_necessity(int *necessity, int size_necessity, const char *annotation)
{
    EXECUTION_REPORT(REPORT_ERROR, comp_id, size_necessity >= fields_mem_registered.size(), "ERROR happens when calling the API \"CCPL_register_import_interface\" to register the export interface \"%s\": the array size (currently is %d) of the parameter \"necessity\" is smaller than the number of field instances (currently is %d) of this interface. Please verify the model code with the annotation \"%s\".", interface_name, size_necessity, fields_mem_registered.size(), annotation);
    check_API_parameter_data_array(comp_id, API_ID_INTERFACE_REG_IMPORT, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comm_group(), "registering an import interface", fields_mem_registered.size(), sizeof(int), (const char*)necessity, "necessity", annotation);
    for (int i = 0; i < fields_mem_registered.size(); i ++) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, necessity[i] == FIELD_NECESSITY_NECESSARY || necessity[i] == FIELD_NECESSITY_OPTIONAL, "ERROR happens when calling the API \"CCPL_register_import_interface\" to register the export interface \"%s\": the number %d value (%d) in the parameter \"necessity\" is not either %d (means necessary) nor %d (means optional). Please verify the model code with the annotation \"%s\".", interface_name, i, necessity[i], FIELD_NECESSITY_NECESSARY, FIELD_NECESSITY_OPTIONAL, annotation);
        imported_fields_necessity.push_back(necessity[i]);
    }
}


int Inout_interface::check_is_import_field_connected(int field_instance_id, const char *annotation)
{
    EXECUTION_REPORT(REPORT_ERROR, comp_id, interface_type == COUPLING_INTERFACE_MARK_IMPORT, "ERROR happens when calling the API \"CCPL_check_is_import_field_connected\": the corresponding coupling interface \"%s\" is not an import interface. Please verify the model code with the annotation \"%s\".", interface_name, annotation);
    int i;
    for (i = 0; i < fields_mem_registered.size(); i ++)
        if (fields_mem_registered[i]->get_field_instance_id() == field_instance_id)
            break;
    EXECUTION_REPORT(REPORT_ERROR, comp_id, i < fields_mem_registered.size(), "ERROR happens when calling the API \"CCPL_check_is_import_field_connected\": the parameter \"field_instance_id\" (currently is 0x%x) fails to specify a field instance in the corresponding coupling interface \"%s\". Please verify the model code with the annotation \"%s\".", field_instance_id, interface_name, annotation);

    return (fields_connected_status[i]? 1 : 0);
}


void Inout_interface::get_sender_time(int size_sender_date, int size_sender_elapsed_days, int size_sender_second, int *sender_date, int *sender_elapsed_days, int *sender_second, const char *annotation)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, interface_type == COUPLING_INTERFACE_MARK_IMPORT, "ERROR happens when calling the API \"CCPL_get_import_fields_sender_time\": the corresponding coupling interface \"%s\" is not an import interface. Please verify the model code with the annotation \"%s\".", interface_name, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, fields_mem_registered.size() <= size_sender_elapsed_days, "ERROR happens when calling the API \"CCPL_get_import_fields_sender_time\" to get the current sender time of the field instances imported by the interface \"%s\": the array size of the input parameter \"sender_date\" (%d) is smaller than the number of fields (%d). Please verify the model code with the annotation \"%s\".", interface_name, size_sender_second, fields_mem_registered.size(), annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, fields_mem_registered.size() <= size_sender_second, "ERROR happens when calling the API \"CCPL_get_import_fields_sender_time\" to get the current sender time of the field instances imported by the interface \"%s\": the array size of the input parameter \"sender_second\" (%d) is smaller than the number of fields (%d). Please verify the model code with the annotation \"%s\".", interface_name, size_sender_second, fields_mem_registered.size(), annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, fields_mem_registered.size() <= size_sender_elapsed_days, "ERROR happens when calling the API \"CCPL_get_import_fields_sender_time\" to get the current sender time of the field instances imported by the interface \"%s\": the array size of the input parameter \"sender_elapsed_days\" (%d) is smaller than the number of fields (%d). Please verify the model code with the annotation \"%s\".", interface_name, size_sender_elapsed_days, fields_mem_registered.size(), annotation);	
	
	for (int i = 0; i < fields_mem_registered.size(); i ++) {
		sender_elapsed_days[i] = CCPL_NULL_INT;
		sender_second[i] = CCPL_NULL_INT;
		sender_date[i] = CCPL_NULL_INT;
		if (fields_coupling_procedures[i] != NULL && fields_coupling_procedures[i]->get_last_receive_sender_time() != CCPL_NULL_LONG) {
			sender_date[i] = fields_coupling_procedures[i]->get_last_receive_sender_time() / 100000;
			sender_second[i] = fields_coupling_procedures[i]->get_last_receive_sender_time() % 100000;
			sender_elapsed_days[i] = time_mgr->calculate_elapsed_day(sender_date[i]/10000, (sender_date[i]%10000)/100, sender_date[i]%100);
		}
	}
}


void Inout_interface::dump_active_coupling_connections_into_XML(TiXmlElement *root_element)
{
    TiXmlElement *parent_element = NULL;
    TiXmlElement *interface = NULL;
    
    for (TiXmlNode *child = root_element->FirstChild(); child != NULL; child = child->NextSibling())
        if(interface_type == COUPLING_INTERFACE_MARK_IMPORT && words_are_the_same(child->Value(), "import_interfaces") || interface_type == COUPLING_INTERFACE_MARK_EXPORT && words_are_the_same(child->Value(), "export_interfaces")) {
            parent_element = child->ToElement();
            break;
        }
    EXECUTION_REPORT(REPORT_ERROR, -1, parent_element != NULL, "Software error Inout_interface::dump_active_coupling_connections_into_XML");
    
    for (TiXmlNode *child = parent_element->FirstChild(); child != NULL; child = child->NextSibling()) {
        TiXmlElement *temp = child->ToElement();
        if(words_are_the_same(temp->Attribute("interface_name"), interface_name)) {
            interface = child->ToElement();
            break;
        }
    }
    
    if (interface == NULL) {
        interface = new TiXmlElement("interface");
        parent_element->LinkEndChild(interface);
        interface->SetAttribute("interface_name", interface_name);
    }
    for (int i = 0; i < coupling_procedures.size(); i++) {
        if (coupling_procedures[i]->get_coupling_connections_dumped())
            continue;
        Coupling_connection *coupling_connection = coupling_procedures[i]->get_coupling_connection();
        TiXmlElement *fields = new TiXmlElement("fields");
        interface->LinkEndChild(fields);
        if (interface_type == COUPLING_INTERFACE_MARK_IMPORT){
            fields->SetAttribute("comp_full_name",coupling_connection->src_comp_interfaces[0].first);
            fields->SetAttribute("interface_name",coupling_connection->src_comp_interfaces[0].second);
        }
        if (interface_type == COUPLING_INTERFACE_MARK_EXPORT){
            fields->SetAttribute("comp_full_name",coupling_connection->dst_comp_full_name);
            fields->SetAttribute("interface_name",coupling_connection->dst_interface_name);
        }
        for (int j = 0; j < coupling_connection->fields_name.size(); j++) {
            TiXmlElement *field = new TiXmlElement("field");
            fields->LinkEndChild(field);
            field->SetAttribute("name",coupling_connection->fields_name[j]);
        }
        coupling_procedures[i]->set_coupling_connections_dumped();
    }
}


bool Inout_interface::is_in_restart_write_window()
{
    if (interface_type != COUPLING_INTERFACE_MARK_IMPORT || is_child_interface)
        return false;
        
    for (int i = 0; i < coupling_procedures.size(); i ++)
        if (coupling_procedures[i]->is_in_restart_write_window() || coupling_procedures[i]->get_is_coupling_time_out_of_execution())
            return true;
    return false;
}


void Inout_interface::dump_active_coupling_connections()
{
    char XML_file_name[NAME_STR_SIZE];
    TiXmlElement *root_element;
    TiXmlDocument *XML_file;
    int i;


    if (interface_type != COUPLING_INTERFACE_MARK_IMPORT && interface_type != COUPLING_INTERFACE_MARK_EXPORT || comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_current_proc_local_id() != 0)
        return;

    for(i = 0; i < coupling_procedures.size(); i ++)
        if(!coupling_procedures[i]->get_coupling_connections_dumped())
            break;
    if (i == coupling_procedures.size())
        return;

    sprintf(XML_file_name, "%s/%s.active_coupling_connections.xml", comp_comm_group_mgt_mgr->get_active_coupling_connections_dir(), comp_full_name);
    XML_file = open_XML_file_to_read(comp_id, XML_file_name, MPI_COMM_NULL, false);
    if (XML_file == NULL) {
        XML_file = new TiXmlDocument;
        TiXmlDeclaration *XML_declaration = new TiXmlDeclaration(("1.0"),(""),(""));
        XML_file->LinkEndChild(XML_declaration);
        root_element = new TiXmlElement("Component");
        XML_file->LinkEndChild(root_element);
        root_element->SetAttribute("name", comp_full_name);
        TiXmlElement *import_interfaces = new TiXmlElement("import_interfaces");
        TiXmlElement *export_interfaces = new TiXmlElement("export_interfaces");
        root_element->LinkEndChild(import_interfaces);
        root_element->LinkEndChild(export_interfaces);
    }
    else root_element = XML_file->RootElement();
    dump_active_coupling_connections_into_XML(root_element);
    EXECUTION_REPORT(REPORT_ERROR, -1, XML_file->SaveFile(XML_file_name), "software error in Inout_interface_mgt::dump_active_coupling_connections: fail to write the XML file %s", XML_file_name);
    delete XML_file;
}


Inout_interface_mgt::Inout_interface_mgt(const char *temp_array_buffer, long buffer_content_iter)
{
    while (buffer_content_iter > 0)
        interfaces.push_back(new Inout_interface(temp_array_buffer, buffer_content_iter));
}


Inout_interface_mgt::~Inout_interface_mgt()
{
    for (int i = 0; i < interfaces.size(); i ++)
        delete interfaces[i];
}


void Inout_interface_mgt::generate_remapping_interface_connection(Inout_interface *new_interface, int num_fields, int *field_ids_src, bool has_frac_remapping) 
{
    EXECUTION_REPORT_LOG(REPORT_LOG, new_interface->get_comp_id(), true, "start to generate the coupling connection of the remapping interface \"%s\"", new_interface->get_interface_name());

    coupling_generator->synchronize_latest_connection_id(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(new_interface->get_comp_id(), ""));
    Coupling_connection *coupling_connection = new Coupling_connection(coupling_generator->apply_connection_id());
    Inout_interface *child_interface_export = new_interface->get_child_interface(0);
    Inout_interface *child_interface_import = new_interface->get_child_interface(1);
    std::pair<const char*, const char*> src_comp_interface;

    if (!has_frac_remapping)
        interfaces.push_back(new_interface);
    strcpy(coupling_connection->dst_comp_full_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(new_interface->get_comp_id(),false,"in Inout_interface_mgt::register_normal_remap_interface")->get_full_name());
    strcpy(coupling_connection->dst_interface_name, child_interface_import->get_interface_name());
    for (int i = 0; i < num_fields; i ++)
        coupling_connection->fields_name.push_back(strdup(memory_manager->get_field_instance(field_ids_src[i])->get_field_name()));
    src_comp_interface.first = strdup(coupling_connection->dst_comp_full_name);
    src_comp_interface.second = strdup(child_interface_export->get_interface_name());
    coupling_connection->src_comp_interfaces.push_back(src_comp_interface);

    interfaces.push_back(child_interface_export);
    interfaces.push_back(child_interface_import);

    coupling_connection->generate_a_coupling_procedure(has_frac_remapping);
//    delete coupling_connection;

    interfaces.erase(interfaces.begin()+interfaces.size()-1);
    interfaces.erase(interfaces.begin()+interfaces.size()-1);

    EXECUTION_REPORT_LOG(REPORT_LOG, new_interface->get_comp_id(), true, "Finish generating the coupling connection of the remapping interface \"%s\"", new_interface->get_interface_name());
}


int Inout_interface_mgt::register_normal_remap_interface(const char *interface_name, int num_fields, int *field_ids_src, int *field_ids_dst, int timer_id, int inst_or_aver, int array_size_src, int array_size_dst, const char *API_label, const char *annotation)
{
    Inout_interface *new_interface = new Inout_interface(interface_name, get_next_interface_id(), num_fields, field_ids_src, field_ids_dst, timer_id, inst_or_aver, array_size_src, array_size_dst, API_label, annotation);
    Inout_interface *existing_interface = get_interface(new_interface->get_comp_id(), interface_name);
    if (existing_interface != NULL)
        EXECUTION_REPORT(REPORT_ERROR, new_interface->get_comp_id(), existing_interface == NULL, "Error happens when calling the API \"%s\" to register an interface named \"%s\" at the model code model with the annotation \"%s\": an interface with the same name has already been registered at the model code with the annotation \"%s\". Please verify.", API_label, interface_name, annotation, annotation_mgr->get_annotation(existing_interface->get_interface_id(), "registering interface"));
    generate_remapping_interface_connection(new_interface, num_fields, field_ids_src, false);

    EXECUTION_REPORT_LOG(REPORT_LOG, new_interface->get_comp_id(), true, "Finish generating a normal remapping interface \"%s\"", new_interface->get_interface_name());
    
    return new_interface->get_interface_id();
}


int Inout_interface_mgt::register_frac_based_remap_interface(const char *interface_name, int num_fields, int *field_ids_src, int *field_ids_dst, int timer_id, int inst_or_aver, int array_size_src, int array_size_dst, void *frac_src, void *frac_dst, int size_frac_src, int size_frac_dst, const char *frac_data_type, const char *API_label, const char *annotation)
{
    int new_remap_interface_id = register_normal_remap_interface(interface_name, num_fields, field_ids_src, field_ids_dst, timer_id, inst_or_aver, array_size_src, array_size_dst, API_label, annotation);
    Inout_interface *new_remap_interface = get_interface(new_remap_interface_id);
    EXECUTION_REPORT_LOG(REPORT_LOG, new_remap_interface->get_comp_id(), true, "Finish generating the normal part for a fraction based remapping interface \"%s\"", new_remap_interface->get_interface_name());
    new_remap_interface->add_remappling_fraction_processing(frac_src, frac_dst, size_frac_src, size_frac_dst, frac_data_type, API_label, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, new_remap_interface->get_comp_id(), true, "Adding fraction process for the remapping interface \"%s\"", new_remap_interface->get_interface_name());
    return new_remap_interface->get_interface_id();
}


int Inout_interface_mgt::register_inout_interface(const char *interface_name, int interface_type, int num_fields, int *field_ids, int array_size, int timer_id, int inst_or_aver, const char *annotation, int interface_source)
{
    int API_id = interface_type == COUPLING_INTERFACE_MARK_IMPORT? API_ID_INTERFACE_REG_IMPORT : API_ID_INTERFACE_REG_EXPORT;
    Inout_interface *new_interface = new Inout_interface(interface_name, get_next_interface_id(), interface_type, num_fields, field_ids, array_size, timer_id, inst_or_aver, "field_instance_IDs", annotation, API_id, interface_source, false);
    for (int i = 0; i < interfaces.size(); i ++) {
        if (new_interface->get_comp_id() != interfaces[i]->get_comp_id())
            continue;
        EXECUTION_REPORT(REPORT_ERROR, new_interface->get_comp_id(), !words_are_the_same(interface_name, interfaces[i]->get_interface_name()), 
                         "Fail to register an import/export interface named \"%s\" at the model code with the annotation \"%s\" because an interface with the same name has been registered at the model code with the annotation \"%s\"",
                         interface_name, annotation, annotation_mgr->get_annotation(interfaces[i]->get_interface_id(), "registering interface"));
        //new_interface->report_common_field_instances(interfaces[i]);
        
    }
    interfaces.push_back(new_interface);

    if (interface_type == COUPLING_INTERFACE_MARK_EXPORT)
        write_comp_export_info_into_XML_file(new_interface->get_comp_id());
    
    return new_interface->get_interface_id();
}


int Inout_interface_mgt::get_next_interface_id()
{
    return TYPE_INOUT_INTERFACE_ID_PREFIX|interfaces.size();
}


bool Inout_interface_mgt::is_interface_id_legal(int interface_id)
{
    if ((interface_id & TYPE_ID_PREFIX_MASK) != TYPE_INOUT_INTERFACE_ID_PREFIX)
        return false;

    return (interface_id&TYPE_ID_SUFFIX_MASK) < interfaces.size();
}


Inout_interface *Inout_interface_mgt::get_interface(int interface_id)
{
    if (!is_interface_id_legal(interface_id))
        return NULL;

    return interfaces[interface_id&TYPE_ID_SUFFIX_MASK];
}


Inout_interface *Inout_interface_mgt::get_interface(const char *comp_full_name, const char *interface_name)
{
    for (int i = 0; i < interfaces.size(); i ++)
        if (words_are_the_same(interfaces[i]->get_comp_full_name(),comp_full_name) && words_are_the_same(interfaces[i]->get_interface_name(), interface_name))
            return interfaces[i];

    return NULL;
}


Inout_interface *Inout_interface_mgt::get_interface(int comp_id, const char *interface_name)
{
    for (int i = 0; i < interfaces.size(); i ++)
        if (interfaces[i]->get_comp_id() == comp_id && words_are_the_same(interfaces[i]->get_interface_name(), interface_name))
            return interfaces[i];

    return NULL;
}


void Inout_interface_mgt::get_all_unconnected_inout_interface_fields_info(std::vector<const char*> &all_descendant_real_comp_fullnames, char **temp_array_buffer, long &buffer_content_size, MPI_Comm comm)
{
    char *local_temp_array_buffer = NULL;
    long local_buffer_max_size = 0, local_buffer_content_size = 0;
    int num_total_local_proc, current_proc_local_id, *all_array_size;


    for (int i = 0; i < all_descendant_real_comp_fullnames.size(); i ++) {
        Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(all_descendant_real_comp_fullnames[i]);
        if (comp_node == NULL || comp_node->get_current_proc_local_id() != 0)
            continue;
        for (int j = 0; j < interfaces.size(); j ++)
            if (interfaces[j]->get_comp_id() == comp_node->get_comp_id())
                interfaces[j]->transform_interface_into_array(&local_temp_array_buffer, local_buffer_max_size, local_buffer_content_size);
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_size(comm, &num_total_local_proc) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Comm_rank(comm, &current_proc_local_id) == MPI_SUCCESS);
    all_array_size = new int [num_total_local_proc];
    gather_array_in_one_comp(num_total_local_proc, current_proc_local_id, local_temp_array_buffer, local_buffer_content_size, sizeof(char), all_array_size, (void**)temp_array_buffer, buffer_content_size, comm);
    delete [] all_array_size;
}


void Inout_interface_mgt::get_all_import_interfaces_of_a_component(std::vector<Inout_interface*> &import_interfaces, int comp_id)
{
    import_interfaces.clear();

    for (int i = 0; i < interfaces.size(); i ++)
        if (interfaces[i]->get_comp_id() == comp_id && interfaces[i]->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT)
            import_interfaces.push_back(interfaces[i]);
}


void Inout_interface_mgt::get_all_import_interfaces_of_a_component(std::vector<Inout_interface*> &import_interfaces, const char *comp_full_name)
{
    import_interfaces.clear();

    for (int i = 0; i < interfaces.size(); i ++)
        if (words_are_the_same(interfaces[i]->get_comp_full_name(), comp_full_name) && interfaces[i]->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT)
            import_interfaces.push_back(interfaces[i]);
}


void Inout_interface_mgt::execute_interface(int interface_id, int API_id, bool bypass_timer, int *field_update_status, int size_field_update_status, int *num_dst_fields, const char *annotation)
{
    Inout_interface *inout_interface;
    char API_label[NAME_STR_SIZE];

    
    get_API_hint(-1, API_id, API_label);

    if (!is_interface_id_legal(interface_id))
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when executing an interface through the API \"%s\": the given interface ID 0x%x is illegal. Please check the model code with the annotation \"%s\"", API_label, interface_id, annotation);
    inout_interface = get_interface(interface_id);
    EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "Begin to execute interface \"%s\" (model code annotation is \"%s\")", inout_interface->get_interface_name(), annotation);    
    inout_interface->execute(bypass_timer, API_id, field_update_status, size_field_update_status, annotation);
    *num_dst_fields = inout_interface->get_num_dst_fields();
    EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "Finish executing interface \"%s\" (model code annotation is \"%s\")", inout_interface->get_interface_name(), annotation);

    inout_interface->dump_active_coupling_connections();
}


void Inout_interface_mgt::execute_interface(int comp_id, int API_id, const char *interface_name, bool bypass_timer, int *field_update_status, int size_field_update_status, int *num_dst_fields, const char *annotation)
{
    Inout_interface *inout_interface;
    char API_label[NAME_STR_SIZE];

    
    get_API_hint(-1, API_id, API_label);
    if (!comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id,true))
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when executing an interface through the API \"%s\": the given component model ID 0x%x is illegal. Please check the model code with the annotation \"%s\"", API_label, comp_id, annotation);
    inout_interface = get_interface(comp_id, interface_name);
    if (inout_interface == NULL)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when executing an interface through the API \"%s\": the corresponding component model \"%s\" does not have an interface with the given name \"%s\". Please check the model code with the annotation \"%s\"", API_label, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,false,"")->get_comp_full_name(), interface_name, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "Begin to execute interface \"%s\" (model code annotation is \"%s\")", inout_interface->get_interface_name(), annotation);    
    inout_interface->execute(bypass_timer, API_id, field_update_status, size_field_update_status, annotation);
    *num_dst_fields = inout_interface->get_num_dst_fields();
    EXECUTION_REPORT_LOG(REPORT_LOG, inout_interface->get_comp_id(), true, "Finish executing interface \"%s\" (model code annotation is \"%s\")", inout_interface->get_interface_name(), annotation);

    inout_interface->dump_active_coupling_connections();
}


void Inout_interface_mgt::runtime_receive_algorithms_receive_data()
{
#ifdef USE_ONE_SIDED_MPI
    for (int i = 0; i < all_runtime_receive_algorithms.size(); i ++)
        all_runtime_receive_algorithms[i]->receive_data_in_temp_buffer();
#endif
}


void Inout_interface_mgt::erase_runtime_receive_algorithm(Runtime_trans_algorithm *new_algorithm)
{
    for (int i = 0; i < all_runtime_receive_algorithms.size(); i ++)
        if (all_runtime_receive_algorithms[i] == new_algorithm) {
            all_runtime_receive_algorithms.erase(all_runtime_receive_algorithms.begin()+i);
            break;
        }
}


void Inout_interface_mgt::free_all_MPI_wins()
{
    for (int i = 0; i < all_MPI_wins.size(); i ++) {
        MPI_Win mpi_win = all_MPI_wins[i];
        MPI_Win_free(&mpi_win);
    }
}


void Inout_interface_mgt::write_into_restart_buffers(int comp_id)
{
    char *array_buffer;
    long buffer_max_size, buffer_content_size;

    
    for (int i = 0; i < interfaces.size(); i ++)
        if (interfaces[i]->get_comp_id() == comp_id) {
            if (!interfaces[i]->has_been_executed_with_timer())
                continue;
            interfaces[i]->write_restart_mgt_info(NULL);
        }            
}


void Inout_interface_mgt::write_comp_export_info_into_XML_file(int comp_id)
{
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
    char XML_file_name[NAME_STR_SIZE];


    if (comp_node->get_current_proc_local_id() == 0) {
        TiXmlDocument *XML_file = new TiXmlDocument;
        TiXmlDeclaration *XML_declaration = new TiXmlDeclaration(("1.0"),(""),(""));
        EXECUTION_REPORT(REPORT_ERROR, -1, XML_file != NULL, "Software error: cannot create an xml file");
        XML_file->LinkEndChild(XML_declaration);
        TiXmlElement *root_element = new TiXmlElement("export_interfaces");
        root_element->SetAttribute("comp_full_name", comp_node->get_comp_full_name());
        XML_file->LinkEndChild(root_element);
        for (int i = 0; i < interfaces.size(); i ++)
            if (interfaces[i]->get_interface_type() == COUPLING_INTERFACE_MARK_EXPORT && interfaces[i]->get_comp_id() == comp_id)
                interfaces[i]->write_export_info_into_XML_file(root_element);
            sprintf(XML_file_name, "%s/%s.exports_info.xml", comp_comm_group_mgt_mgr->get_components_exports_dir(), comp_node->get_full_name());
            EXECUTION_REPORT(REPORT_ERROR, -1, XML_file->SaveFile(XML_file_name), "Software error in Inout_interface_mgt::write_comp_export_info_into_XML_file: fail to write the XML file %s", XML_file_name);
            delete XML_file;    
    }
    MPI_Barrier(comp_node->get_comm_group());
}


void Inout_interface_mgt::get_all_export_interfaces_of_a_field(int comp_id, const char *field_name, std::vector<const char*> &export_comp_full_names, std::vector<const char*> &export_interface_names)
{
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
    int line_number;
    char XML_file_name[NAME_STR_SIZE];


    for (int i = 0; i < export_comp_full_names.size(); i ++) {
        delete [] export_comp_full_names[i];
        delete [] export_interface_names[i];
    }
    export_comp_full_names.clear();
    export_interface_names.clear();
     
    if (comp_node->get_current_proc_local_id() == 0) {
        DIR *cur_dir = opendir(comp_comm_group_mgt_mgr->get_components_exports_dir());
        struct dirent *ent = NULL;
        struct stat st;
        EXECUTION_REPORT(REPORT_ERROR, -1, cur_dir != NULL, "Software error in Inout_interface_mgt::get_all_export_interfaces_of_a_field");
        while ((ent = readdir(cur_dir)) != NULL) {
            stat(ent->d_name, &st);
            if (!(strlen(ent->d_name) > 4 && words_are_the_same(ent->d_name+strlen(ent->d_name)-4, ".xml")))
                continue;
            sprintf(XML_file_name, "%s/%s", comp_comm_group_mgt_mgr->get_components_exports_dir(), ent->d_name);
            TiXmlDocument *XML_file = open_XML_file_to_read(comp_id, XML_file_name, MPI_COMM_NULL, false);
            EXECUTION_REPORT(REPORT_ERROR, -1, XML_file != NULL, "Software error in Inout_interface_mgt::get_all_export_interfaces_of_a_field: no XML file");
            TiXmlElement *root_element = XML_file->FirstChildElement();
            const char *comp_full_name = get_XML_attribute(comp_id, -1, root_element, "comp_full_name", XML_file_name, line_number, "the full name of the corresponding component model", "internal XML files generated by C-Coupler", true);
            for (TiXmlNode *export_interface_node = root_element->FirstChildElement(); export_interface_node != NULL; export_interface_node = export_interface_node->NextSibling()) {
                if (export_interface_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                    continue;
                TiXmlElement *export_interface_element = export_interface_node->ToElement();
                const char *interface_name = get_XML_attribute(comp_id, -1, export_interface_element, "interface_name", XML_file_name, line_number, "the name of the export interface", "internal XML files generated by C-Coupler", true);
                for (TiXmlNode *export_field_node = export_interface_element->FirstChildElement(); export_field_node != NULL; export_field_node = export_field_node->NextSibling()) {
                    if (export_field_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                        continue;
                    TiXmlElement *export_field_element = export_field_node->ToElement();
                    const char *XML_field_name = get_XML_attribute(comp_id, -1, export_field_element, "field_name", XML_file_name, line_number, "the name of the export field", "internal XML files generated by C-Coupler", true);
                    if (words_are_the_same(XML_field_name, field_name)) {
                        export_comp_full_names.push_back(strdup(comp_full_name));
                        export_interface_names.push_back(strdup(interface_name));
                    }
                }
            }
            
            delete XML_file;
        }
    }
    MPI_Barrier(comp_node->get_comm_group());
}


Inout_interface *Inout_interface_mgt::search_an_inout_interface_executed_with_timer(int comp_id)
{
    for (int i = 0; i < interfaces.size(); i ++)
        if (interfaces[i]->get_comp_id() == comp_id && interfaces[i]->get_interface_type() < 2 && interfaces[i]->has_been_executed_with_timer())
            return interfaces[i];

    return NULL;
}


int Inout_interface_mgt::get_h2d_grid_area_in_remapping_weights(int interface_id, int field_index, void *output_area_data, int area_array_size, const char *data_type, const char *annotation)
{
    Inout_interface *inout_interface = get_interface(interface_id);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, inout_interface != NULL, "ERROR happens when calling the API \"CCPL_get_H2D_grid_area_in_remapping_wgts\": the parameter of interface ID is wrong. Please verify the model code with the annotation \"%s\"", annotation);
    return inout_interface->get_h2d_grid_area_in_remapping_weights(inout_interface->get_interface_name(), field_index, output_area_data, area_array_size, data_type, annotation);
}


bool Inout_interface_mgt::is_comp_in_restart_write_window(int comp_id)
{
    for (int i = 0; i < interfaces.size(); i ++)
        if (interfaces[i]->get_comp_id() == comp_id && interfaces[i]->is_in_restart_write_window())
            return true;
        
    return false;    
}


