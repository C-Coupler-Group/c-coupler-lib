/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "IO_field_mgt.h"


IO_field::IO_field(int IO_field_id, int field_instance_id, const char *field_IO_name, const char *annotation)
{
    this->IO_field_id = IO_field_id;

    Field_mem_info *field_inst = memory_manager->get_field_instance(field_instance_id);
    this->comp_id = field_inst->get_comp_id();
    this->field_instance_id = field_instance_id;
    if (strlen(field_IO_name) > 0) {
        strcpy(this->field_IO_name, field_IO_name);
        check_and_verify_name_format_of_string_for_API(this->comp_id, field_IO_name, API_ID_FIELD_MGT_REG_IO_FIELD_from_INST, "name of the I/O field in the data file", annotation);
    }
    else strcpy(this->field_IO_name, field_inst->get_field_name());
    strcpy(this->field_long_name, fields_info->search_field_info(field_inst->get_field_name())->field_long_name);
    strcpy(this->field_unit, field_inst->get_unit());
}


IO_field::IO_field(int IO_field_id, int comp_or_grid_id, int decomp_id, int field_size, void *data_buffer, const char * field_IO_name, const char *long_name, const char *unit, const char *data_type, const char * annotation)
{
    Field_mem_info *field_mem;


    this->IO_field_id = IO_field_id;
    check_and_verify_name_format_of_string_for_API(this->comp_id, field_IO_name, API_ID_FIELD_MGT_REG_IO_FIELD_from_BUFFER, "name of the I/O field in the data file", annotation);
    strcpy(this->field_IO_name, field_IO_name);
    strcpy(this->field_unit, unit);
    strcpy(this->field_long_name, long_name);

    if (decomp_id == -1) {
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,true), "The parameter of component id when calling the CCPL interface \"CCPL_register_IO_field\" for registering IO field \"%s\" is wrong. Please verify the model code with the annotation \"%s\"", field_IO_name, annotation);
        this->comp_id = comp_or_grid_id;
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(decomp_id), "The parameter of decomposition id when calling the CCPL interface \"CCPL_register_IO_field\" for registering IO field \"%s\" is wrong. Please verify the model code with the annotation \"%s\"", field_IO_name, annotation);
        this->comp_id = decomps_info_mgr->get_decomp_info(decomp_id)->get_comp_id();
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "The parameter of grid id when calling the CCPL interface \"CCPL_register_IO_field\" for registering IO field \"%s\" is wrong. Please verify the model code with the annotation \"%s\"", field_IO_name, annotation);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_comp_id_of_grid(comp_or_grid_id), "The parameters of grid id and decomposition id when calling the CCPL interface \"CCPL_register_IO_field\" for registering IO field \"%s\" do not belong to the same component. Please verify the model code with the annotation \"%s\"", field_IO_name, annotation);
    }

    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(field_IO_name) > 0, "The parameter of field I/O name when calling the CCPL interface \"CCPL_register_IO_field\" is empty. Please verify the model code with the annotation \"%s\"", annotation);

    field_instance_id = memory_manager->register_external_field_instance(field_IO_name, data_buffer, field_size, decomp_id, comp_or_grid_id, BUF_MARK_IO_FIELD_REG, REG_FIELD_TAG_REST, unit, data_type, annotation);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(field_IO_name) > 0, "The parameter of field I/O name when calling the CCPL interface \"CCPL_register_IO_field\" cannot be an empty string. Please verify the model code with the annotation \"%s\"", annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(long_name) > 0, "The parameter of long name of the I/O field \"%s\" when calling the CCPL interface \"CCPL_register_IO_field\" cannot be an empty string. Please verify the model code with the annotation \"%s\"", field_IO_name, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, strlen(unit) > 0, "The parameter of unit of the I/O field \"%s\" when calling the CCPL interface \"CCPL_register_IO_field\" cannot be an empty string. Please verify the model code with the annotation \"%s\"", field_IO_name, annotation);
}


IO_field *IO_field_mgt::search_IO_field(int comp_id, const char *field_IO_name)
{
    for (int i = 0; i < IO_fields.size(); i ++)
        if (IO_fields[i]->get_comp_id() == comp_id && words_are_the_same(IO_fields[i]->get_field_IO_name(), field_IO_name))
            return IO_fields[i];

    return NULL;
}


void IO_field_mgt::check_for_registering_IO_field(IO_field *new_IO_field, const char *annotation, int API_id)
{
    IO_field *existing_field = search_IO_field(new_IO_field->get_comp_id(), new_IO_field->get_field_IO_name());
    if (existing_field != NULL)
        EXECUTION_REPORT(REPORT_ERROR, new_IO_field->get_comp_id(), false, "IO field \"%s\" has been registered before (the corresponding model code annotation is \"%s\"). It cannot be registered again at the model code with the annotation \"%s\"",
                         new_IO_field->get_field_IO_name(), annotation_mgr->get_annotation(existing_field->get_IO_field_id(), "registering I/O field"), annotation);
    annotation_mgr->add_annotation(new_IO_field->get_IO_field_id(), "registering I/O field", annotation);
    synchronize_comp_processes_for_API(new_IO_field->get_comp_id(), API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(new_IO_field->get_comp_id(),""), "registering an I/O field", annotation);
    IO_fields.push_back(new_IO_field);
}


IO_field_mgt::~IO_field_mgt()
{
    for (int i = 0; i < IO_fields.size(); i ++)
        delete IO_fields[i];
}


int IO_field_mgt::register_IO_field(int field_instance_id, const char *field_IO_name, const char *annotation)
{
    int IO_field_id = TYPE_IO_FIELD_PREFIX | IO_fields.size();
    IO_field *new_IO_field = new IO_field(IO_field_id, field_instance_id, field_IO_name, annotation);
    check_for_registering_IO_field(new_IO_field, annotation, API_ID_FIELD_MGT_REG_IO_FIELD_from_INST);
    check_API_parameter_field_instance(new_IO_field->get_comp_id(), API_ID_FIELD_MGT_REG_IO_FIELD_from_INST, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(new_IO_field->get_comp_id(),""), "registering an I/O field", new_IO_field->get_field_instance_id(), "field_inst_id", annotation);
    return IO_field_id;
}


int IO_field_mgt::register_IO_fields(int num_field_inst, int size_field_inst_ids, int *field_inst_ids, const char *annotation)
{
    int comp_id = memory_manager->get_field_instance(field_inst_ids[0])->get_comp_id();
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"");


    comp_comm_group_mgt_mgr->confirm_coupling_configuration_active(comp_id, API_ID_FIELD_MGT_REG_IO_FIELDs_from_INSTs, true, annotation);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, num_field_inst <= size_field_inst_ids, "Error happers when calling the API \"CCPL_register_IO_fields_from_field_instances\": the array size of the parameter \"field_inst_ids\" cannot be smaller than the parameter \"num_field_inst\". Please check the model code with the annotation \"%s\".", annotation);
    for (int i = 1; i < num_field_inst; i ++)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, comp_id == memory_manager->get_field_instance(field_inst_ids[i])->get_comp_id(), "Error happers when calling the API \"CCPL_register_IO_fields_from_field_instances\": the field instances specified by the parameter \"field_inst_ids\" do not correspond to the same component model. Please check the model code with the annotation \"%s\".", annotation);
    synchronize_comp_processes_for_API(comp_id, API_ID_FIELD_MGT_REG_IO_FIELDs_from_INSTs, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,""), "registering I/O fields", annotation);
    check_API_parameter_int(comp_id, API_ID_FIELD_MGT_REG_IO_FIELDs_from_INSTs, comm, NULL, num_field_inst, "\"num_field_inst\"", annotation);
    for (int i = 0; i < num_field_inst; i ++)
        check_API_parameter_field_instance(comp_id, API_ID_FIELD_MGT_REG_IO_FIELDs_from_INSTs, comm, "registering I/O fields", field_inst_ids[i], "field_inst_ids (detailed field instances)", annotation);    
    for (int i = 0; i < num_field_inst; i ++)
        IO_fields.push_back(new IO_field(TYPE_IO_FIELD_PREFIX|IO_fields.size(), field_inst_ids[i], memory_manager->get_field_instance(field_inst_ids[i])->get_field_name(), annotation));
    
    return 0;
}


int IO_field_mgt::register_IO_field(int comp_or_grid_id, int decomp_id, int field_size, void *data_buffer, const char * field_IO_name, const char *long_name, const char *unit, const char *data_type, const char * annotation)
{
    int IO_field_id = TYPE_IO_FIELD_PREFIX | IO_fields.size();
    IO_field *new_IO_field = new IO_field(IO_field_id, comp_or_grid_id, decomp_id, field_size, data_buffer, field_IO_name, long_name, unit, data_type, annotation);
    check_for_registering_IO_field(new_IO_field, annotation, API_ID_FIELD_MGT_REG_IO_FIELD_from_BUFFER);
    return IO_field_id;
}


IO_output_procedure::~IO_output_procedure()
{
    if (field_update_status != NULL)
        delete [] field_update_status;
}


IO_output_procedure::IO_output_procedure(int comp_id, int procedure_id, Coupling_timer *default_field_timer, Coupling_timer *default_file_timer, bool synchronized_IO)
{
    this->comp_id = comp_id;
    this->procedure_id = procedure_id;
    inst_or_aver = USING_AVERAGE_VALUE;
    import_interface = NULL;
    export_interface = NULL;
    time_mgr = components_time_mgrs->get_time_mgr(comp_id);
    netcdf_file_object = NULL;
    write_grid_name = false;

    include_all_component_io_fields();

    field_update_status = NULL;
    
    if (IO_fields.size() == 0)
        return;

    field_update_status = new int [IO_fields.size()];
    
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id,true), "Software error in IO_output_procedure::IO_output_procedure: wrong comp id");
    
    if (default_field_timer == NULL)
        field_timer = timer_mgr->get_timer(timer_mgr->define_timer(comp_id, FREQUENCY_UNIT_DAYS, 1, 0, 0, "default timer for I/O fields"));
    else {
        field_timer = timer_mgr->get_timer(timer_mgr->define_timer(comp_id, default_field_timer));
        field_timer->reset_remote_lag_count();
    }
    
    if (default_file_timer == NULL)
        file_timer = timer_mgr->get_timer(timer_mgr->define_timer(comp_id, FREQUENCY_UNIT_DAYS, 10, 0, 0, "default timer for I/O fields"));
    else {
        file_timer = timer_mgr->get_timer(timer_mgr->define_timer(comp_id, default_file_timer));
        file_timer->reset_remote_lag_count();
    }

    int *fields_id = new int [IO_fields.size()];
    int field_timer_id = field_timer->get_timer_id();
    for (int i = 0; i < IO_fields.size(); i ++)
        fields_id[i] = IO_fields[i]->get_field_instance_id();

    export_interface = inout_interface_mgr->get_interface(inout_interface_mgr->register_inout_interface("Default_IO_output", 1, IO_fields.size(), fields_id, IO_fields.size(), field_timer_id, inst_or_aver, "register default IO field to output for a component", INTERFACE_SOURCE_IO_OUTPUT));

    if (synchronized_IO) {
        for (int i = 0; i < IO_fields.size(); i ++) {
            Field_mem_info *IO_field_instance, *mirror_field_instance;
            IO_field_instance = memory_manager->get_field_instance(IO_fields[i]->get_field_instance_id());
            const char *data_type = IO_field_instance->get_data_type();
            if (words_are_the_same(data_type, DATA_TYPE_DOUBLE) || words_are_the_same(data_type, DATA_TYPE_FLOAT))
                mirror_field_instance = memory_manager->alloc_mem(IO_field_instance, BUF_MARK_IO_FIELD_MIRROR, inout_interface_mgr->get_next_interface_id(), DATA_TYPE_FLOAT, false);
            else mirror_field_instance = memory_manager->alloc_mem(IO_field_instance, BUF_MARK_IO_FIELD_MIRROR, inout_interface_mgr->get_next_interface_id(), DATA_TYPE_INT, false);
            fields_id[i] = mirror_field_instance->get_field_instance_id();
            data_write_field_insts.push_back(mirror_field_instance);
        }
        import_interface = inout_interface_mgr->get_interface(inout_interface_mgr->register_inout_interface("Default_IO_write", 0, IO_fields.size(), fields_id, IO_fields.size(), field_timer_id, inst_or_aver, "register default IO field to write for a component", INTERFACE_SOURCE_IO_WRITE));
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "asynchronized IO is not supported yet");
    }

    delete [] fields_id;
}


void IO_output_procedure::execute()
{
    if (export_interface != NULL)
        export_interface->execute(false, API_ID_INTERFACE_EXECUTE_WITH_ID, field_update_status, IO_fields.size(), "IO output procedure export interface execute");

    if (import_interface != NULL) {
        import_interface->execute(false, API_ID_INTERFACE_EXECUTE_WITH_ID, field_update_status, IO_fields.size(), "IO output procedure import interface execute");
        if (field_timer->is_timer_on()) {
            if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in IO_output_procedure::execute") == 0) {
                if (netcdf_file_object == NULL || file_timer->is_timer_on()) {
                    char time_string[NAME_STR_SIZE];
                    char full_file_name[NAME_STR_SIZE];
                    char file_header[NAME_STR_SIZE];
                    if (IS_TIME_UNIT_SECOND(field_timer->get_frequency_unit()))
                        sprintf(time_string, "%04d%02d%02d-%05d", time_mgr->get_current_year(), time_mgr->get_current_month(), time_mgr->get_current_day(), time_mgr->get_current_second());
                    else if (IS_TIME_UNIT_DAY(field_timer->get_frequency_unit()))
                        sprintf(time_string, "%04d%02d%02d", time_mgr->get_current_year(), time_mgr->get_current_month(), time_mgr->get_current_day());    
                    else if (IS_TIME_UNIT_MONTH(field_timer->get_frequency_unit()))
                        sprintf(time_string, "%04d%02d", time_mgr->get_current_year(), time_mgr->get_current_month());    
                    else if (IS_TIME_UNIT_YEAR(field_timer->get_frequency_unit()))
                        sprintf(time_string, "%04d", time_mgr->get_current_year());
                    comp_comm_group_mgt_mgr->get_output_data_file_header(comp_id, file_header);
                    sprintf(full_file_name, "%s.%s.h%d.nc",file_header, time_string, procedure_id);
                    if (netcdf_file_object != NULL)
                        delete netcdf_file_object;
                    netcdf_file_object = new IO_netcdf(full_file_name, full_file_name, "w", true);
                    // compset_communicators_info_mgr->write_case_info(netcdf_file_object);   // to be modify shortly
                }
            }
            for (int i = 0; i < data_write_field_insts.size(); i ++) {
                data_write_field_insts[i]->check_field_sum("before writing data into a file");
                fields_gather_scatter_mgr->gather_write_field(netcdf_file_object, data_write_field_insts[i], write_grid_name, time_mgr->get_current_date(), time_mgr->get_current_second(), false);
            }
        }
    }
}


Coupling_connection *IO_output_procedure::generate_coupling_connection(int connection_id)
{
    Coupling_connection *coupling_connection = NULL;
    
    if (import_interface != NULL && export_interface != NULL) {        
        coupling_generator->synchronize_latest_connection_id(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, ""));
        coupling_connection = new Coupling_connection(coupling_generator->apply_connection_id());
        strcpy(coupling_connection->dst_comp_full_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,true,"in IO_output_procedure::generate_coupling_connection")->get_full_name());
        strcpy(coupling_connection->dst_interface_name, import_interface->get_interface_name());
        std::pair<const char*,const char*> src_comp_interface;
        src_comp_interface.first = strdup(coupling_connection->dst_comp_full_name);
        src_comp_interface.second = strdup(export_interface->get_interface_name());
        coupling_connection->src_comp_interfaces.push_back(src_comp_interface);
        std::vector<const char*> import_fields_name;
        import_interface->get_fields_name(&import_fields_name);
        for (int k = 0; k < import_fields_name.size(); k ++)
            coupling_connection->fields_name.push_back(strdup(import_fields_name[k]));
    }

    return coupling_connection;
}


void IO_output_procedure::include_all_component_io_fields()
{
    for (int i = 0; i < IO_fields_mgr->IO_fields.size(); i ++)
        if (IO_fields_mgr->IO_fields[i]->get_comp_id() == comp_id)
            IO_fields.push_back(IO_fields_mgr->IO_fields[i]);
}


Component_IO_output_procedures::Component_IO_output_procedures(int comp_id, const char *xml_file_name, bool synchronized_IO)
{
    this->comp_id = comp_id;
    
    if (xml_file_name == NULL)
        IO_output_procedures.push_back(new IO_output_procedure(comp_id, IO_output_procedures.size(), NULL, NULL, synchronized_IO));
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "asynchronized IO is not supported yet");
    }
}


Component_IO_output_procedures::~Component_IO_output_procedures()
{
    for (int i = 0; i < IO_output_procedures.size(); i ++)
        delete IO_output_procedures[i];
}


void Component_IO_output_procedures::generate_coupling_connection(std::vector<Coupling_connection*> &all_IO_connections, int basic_connection_id)
{
    for (int i = 0; i < IO_output_procedures.size(); i ++) {
        Coupling_connection *coupling_connection = IO_output_procedures[i]->generate_coupling_connection(basic_connection_id+IO_output_procedures.size());
        if (coupling_connection != NULL)
            all_IO_connections.push_back(coupling_connection);
    }
}


void Component_IO_output_procedures::execute()
{
    for (int i = 0; i < IO_output_procedures.size(); i ++)
        IO_output_procedures[i]->execute();
}


Components_IO_output_procedures_mgt::~Components_IO_output_procedures_mgt()
{
    for (int i = 0; i < components_IO_output_procedures.size(); i ++)
        delete components_IO_output_procedures[i];
}


void Components_IO_output_procedures_mgt::add_component_IO_output_procedures(int comp_id, const char *xml_file_name, bool synchronized_IO)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id,true), "Software error in Component_IO_output_procedures::Component_IO_output_procedures: wrong comp id");
    // ...
    int true_comp_id = (comp_id&TYPE_ID_SUFFIX_MASK);
    for (int i = components_IO_output_procedures.size(); i <= true_comp_id; i ++)
        components_IO_output_procedures.push_back(NULL);
    EXECUTION_REPORT(REPORT_ERROR, -1, components_IO_output_procedures[true_comp_id] == NULL, "Software error in Component_IO_output_procedures::Component_IO_output_procedures: wrong");    
    components_IO_output_procedures[true_comp_id] = new Component_IO_output_procedures(comp_id, xml_file_name, synchronized_IO);
}


void Components_IO_output_procedures_mgt::add_all_components_IO_output_procedures()
{
    const int *all_components_ids = comp_comm_group_mgt_mgr->get_all_components_ids();
    for (int i = 1; i < all_components_ids[0]; i ++) {
        if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(all_components_ids[i], "in add_all_components_IO_output_procedures") != -1)
            add_component_IO_output_procedures(all_components_ids[i], NULL, true);
    }
}


Component_IO_output_procedures *Components_IO_output_procedures_mgt::get_component_IO_output_procedures(int comp_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id,true), "Software error in Components_IO_output_procedures_mgt::get_component_IO_output_procedures");
    int true_comp_id = comp_id & TYPE_ID_SUFFIX_MASK;
    EXECUTION_REPORT(REPORT_ERROR, -1, components_IO_output_procedures.size() > true_comp_id && components_IO_output_procedures[true_comp_id] != NULL, "Software error in Components_IO_output_procedures_mgt::get_component_IO_output_procedures");
    return components_IO_output_procedures[true_comp_id];
}


