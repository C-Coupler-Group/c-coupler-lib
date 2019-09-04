/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "compset_communicators_info_mgt.h"
#include <stdio.h>
#include <string.h>
#include "global_data.h"
#include "cor_global_data.h"
#include "quick_sort.h"
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>


#define LOG_BUFFER_MAX_SIZE              ((int)1024*1024*5)
#define LOG_BUFFER_MAX_CONTENT_SIZE      (LOG_BUFFER_MAX_SIZE/5*4)


void output_CCPL_log(const char *log_string, const char *log_file_name, char **log_buffer, int &log_buffer_content_size, bool flush_log_file)
{
    if (*log_buffer == NULL) {
        *log_buffer = new char [LOG_BUFFER_MAX_SIZE];
        *log_buffer[0] = '\0';
        log_buffer_content_size = 0;
    }

    strcat(*log_buffer, log_string);
    log_buffer_content_size += strlen(log_string);

    if (log_buffer_content_size >= LOG_BUFFER_MAX_CONTENT_SIZE)
        flush_log_file = true;

    if (flush_log_file) {
        FILE *log_file = stdout;
        if (log_file_name != NULL)
            log_file = fopen(log_file_name, "a+");
        fprintf(log_file, *log_buffer);
        fflush(log_file);
        if (log_file_name != NULL)
            fclose(log_file);
        *log_buffer[0] = '\0';
        log_buffer_content_size = 0;
    }
}


void recursively_remove_directory()
{
    DIR *cur_dir = opendir(".");
    struct dirent *ent = NULL;
    struct stat st;
 
    if (cur_dir == NULL)
        return;
 
    while ((ent = readdir(cur_dir)) != NULL) {
        if (stat(ent->d_name, &st) != 0)
			continue;
     
        if (words_are_the_same(ent->d_name, ".") || words_are_the_same(ent->d_name, ".."))
            continue;
 
        if (S_ISDIR(st.st_mode)) {
		    char old_path[NAME_STR_SIZE];
		    getcwd(old_path, NAME_STR_SIZE);
            chdir(ent->d_name);
            recursively_remove_directory();
            chdir(old_path);
        }
        remove(ent->d_name);
    }
     
    closedir(cur_dir);
}


void remove_directory(const char *path)
{
    char old_path[NAME_STR_SIZE];
 
    getcwd(old_path, NAME_STR_SIZE);
     
    if (chdir(path) == -1)
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "Fail to open \"%s\" which should be a directory");
    
    recursively_remove_directory();
    chdir(old_path);
}


void create_directory(const char *path, MPI_Comm comm, bool is_root_proc, bool new_dir)
{
    char buffer[NAME_STR_SIZE];
    

    if (is_root_proc) {
        DIR *dir=opendir(path);
        if (dir != NULL && new_dir)
            remove_directory(path);
        if (dir == NULL) {
            umask(0);
            mkdir(path, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
            int retcode = errno;
            dir=opendir(path);
            EXECUTION_REPORT(REPORT_ERROR, -1, dir != NULL, "Directory \"%s\" cannot be created: %s. Please verify.", path, strerror(retcode));
        }
    }    
}


Comp_comm_group_mgt_node::~Comp_comm_group_mgt_node()
{
    if (performance_timing_mgr != NULL)
        delete performance_timing_mgr;
    
    if (log_buffer != NULL) {
        output_log("", true);
        delete [] log_buffer;
    }

    if (temp_array_buffer != NULL)
        delete [] temp_array_buffer;

    if (restart_mgr != NULL)
        delete restart_mgr;

    if (proc_latest_model_time != NULL)
        delete [] proc_latest_model_time;
}


Comp_comm_group_mgt_node::Comp_comm_group_mgt_node(const char *comp_name, const char *comp_type, int comp_id, Comp_comm_group_mgt_node *parent, MPI_Comm &comm, bool enabled_in_parent_coupling_gen, const char *annotation)
{
    std::vector<char*> unique_comp_name;
    int i, j, num_procs, current_proc_local_id_in_parent, *process_comp_id, *processes_global_id;
    MPI_Comm parent_comm;
    char dir[NAME_STR_SIZE];
    Comp_comm_group_mgt_node *ancestor = parent;

    
    strcpy(this->comp_name, comp_name);
    strcpy(this->comp_type, comp_type);
    if (parent != NULL && words_are_the_same(parent->get_comp_type(), COMP_TYPE_PSEUDO_COUPLED) && parent->children.size() > 0)
        EXECUTION_REPORT(REPORT_ERROR, parent->get_comp_id(), false, "Error happens when registering the component model \"%s\": its parent \"%s\" is an inactive component model and already has one child \"%s\". Please note that an inactive component model can have at most one child. Please check the model code with the annotation \"%s\"", comp_name, parent->comp_name, parent->children[0]->comp_name, annotation);

    while(ancestor != NULL && words_are_the_same(ancestor->get_comp_type(), COMP_TYPE_PSEUDO_COUPLED)) 
        ancestor = ancestor->get_parent();
    if (ancestor == NULL || words_are_the_same(ancestor->get_comp_name(), COMP_TYPE_ROOT))
        strcpy(this->full_name, this->comp_name);
    else sprintf(this->full_name, "%s@%s", ancestor->get_full_name(), this->comp_name);
    strcpy(this->annotation_start, annotation);
    this->annotation_end[0] = '\0';
    this->comp_id = comp_id;
    this->parent = parent;
    this->buffer_content_size = 0;
    this->buffer_max_size = 1024;
    this->temp_array_buffer = new char [buffer_max_size];
    this->definition_finalized = false;    
    this->proc_latest_model_time = NULL;
    this->enabled_in_parent_coupling_generation = enabled_in_parent_coupling_gen;
    this->log_buffer = NULL;
    this->performance_timing_mgr = new Performance_timing_mgt(comp_id);
    restart_mgr = new Restart_mgt(this);
    comp_ccpl_log_file_name[0] = '\0';
    comp_model_log_file_name[0] = '\0';
    comp_model_log_file_device = -1;
    min_remote_lag_seconds = 0;
    max_remote_lag_seconds = 0;

    if (comm != MPI_COMM_NULL) {
        comm_group = comm;
        if (parent == NULL)
            synchronize_comp_processes_for_API(-1, API_ID_COMP_MGT_REG_COMP, comm, "checking the given communicator for registering root component", annotation);
        else {
            char tmp_string[NAME_STR_SIZE];
            sprintf(tmp_string, "for checking the given communicator for registering a child component \"%s\"", comp_name);
            synchronize_comp_processes_for_API(parent->get_comp_id(), API_ID_COMP_MGT_REG_COMP, comm, tmp_string, annotation);            
            check_API_parameter_string(parent->get_comp_id(), API_ID_COMP_MGT_REG_COMP, comm, "registering a component model", parent->get_comp_name(), "\"parent_id\" (the parent component model)", annotation);
        }    
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR,-1, parent != NULL, "Software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node for checking parent");
        synchronize_comp_processes_for_API(parent->get_comp_id(), API_ID_COMP_MGT_REG_COMP, parent->get_comm_group(), "checking the communicator of the current component for registering its children component", annotation);
        check_API_parameter_string(parent->get_comp_id(), API_ID_COMP_MGT_REG_COMP, parent->get_comm_group(), "registering a component model", parent->get_comp_name(), "\"parent_id\" (the parent component model)", annotation);
        parent_comm = parent->get_comm_group();
        if ((parent->comp_id&TYPE_ID_SUFFIX_MASK) != 0)
            EXECUTION_REPORT_LOG(REPORT_LOG, parent->comp_id, true, 
                             "Before the MPI_barrier for synchronizing all processes of the parent component \"%s\" for registering its children components including \"%s\" (the corresponding model code annotation is \"%s\")", 
                             parent->get_comp_name(), comp_name, annotation);
        else if (parent->get_current_proc_local_id() == 0) 
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, 
                             "Before the MPI_barrier for synchronizing all processes of the whole coupled model for registering root components including \"%s\" (the corresponding model code annotation is \"%s\")", 
                             comp_name, annotation);    
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Barrier(parent_comm) == MPI_SUCCESS);
        if ((parent->comp_id&TYPE_ID_SUFFIX_MASK) != 0)
            EXECUTION_REPORT_LOG(REPORT_LOG, parent->comp_id, true, 
                             "After the MPI_barrier for synchronizing all processes of the parent component \"%s\" for registering its children components including \"%s\" (the corresponding model code annotation is \"%s\")", 
                             parent->get_comp_name(), comp_name, annotation);
        else if (parent->get_current_proc_local_id() == 0) 
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, 
                             "After the MPI_barrier for synchronizing all processes of the whole coupled model for registering root components including \"%s\" (the corresponding model code annotation is \"%s\")", 
                             comp_name, annotation);    
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_size(parent_comm, &num_procs) == MPI_SUCCESS);
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_rank(parent_comm, &current_proc_local_id_in_parent) == MPI_SUCCESS);
        char *all_comp_name;
        if (current_proc_local_id_in_parent == 0) 
            all_comp_name = new char [NAME_STR_SIZE*num_procs];
        process_comp_id = new int [num_procs];
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Gather((char*)comp_name, NAME_STR_SIZE, MPI_CHAR, all_comp_name, NAME_STR_SIZE, MPI_CHAR, 0, parent_comm) == MPI_SUCCESS);
        if (current_proc_local_id_in_parent == 0) {
            unique_comp_name.push_back(all_comp_name);
            process_comp_id[0] = 0;
            for (i = 1; i < num_procs; i ++) {
                for (j = 0; j < unique_comp_name.size(); j ++)
                    if (words_are_the_same(unique_comp_name[j], all_comp_name+i*NAME_STR_SIZE)) {
                        break;
                    }
                if (j == unique_comp_name.size())
                    unique_comp_name.push_back(all_comp_name+i*NAME_STR_SIZE);
                process_comp_id[i] = j;
            }
        }
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Bcast(process_comp_id, num_procs, MPI_INT, 0, parent_comm)  == MPI_SUCCESS);
        EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_split(parent_comm, process_comp_id[current_proc_local_id_in_parent], 0, &comm_group) == MPI_SUCCESS);
        if (current_proc_local_id_in_parent == 0)
            delete [] all_comp_name;
        delete [] process_comp_id;
        comm = comm_group;
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(comp_type, COMP_TYPE_CPL) || words_are_the_same(comp_type, COMP_TYPE_ATM) || words_are_the_same(comp_type, COMP_TYPE_ATM_CHEM) || words_are_the_same(comp_type, COMP_TYPE_OCN) ||
                     words_are_the_same(comp_type, COMP_TYPE_LND) || words_are_the_same(comp_type, COMP_TYPE_SEA_ICE) || words_are_the_same(comp_type, COMP_TYPE_WAVE) || words_are_the_same(comp_type, COMP_TYPE_ROOT) || 
                     words_are_the_same(comp_type, COMP_TYPE_PSEUDO_COUPLED) || words_are_the_same(comp_type, COMP_TYPE_ACTIVE_COUPLED) || words_are_the_same(comp_type, COMP_TYPE_GLC) || words_are_the_same(comp_type, COMP_TYPE_RUNOFF), 
                     "Error happens when registering the component model \"%s\" at the model code with the annotation is \"%s\": the model type \"%s\" is wrong. Please verify.", comp_name, annotation, comp_type);    
    if (parent != NULL && words_are_the_same(comp_type, COMP_TYPE_PSEUDO_COUPLED))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(parent->comp_type, COMP_TYPE_PSEUDO_COUPLED) || words_are_the_same(parent->comp_type, COMP_TYPE_ROOT), 
                         "Error happens when calling the API \"CCPL_register_component\" to register a component model \"%s\" of type \"pesudo_coupled_system\": the type of its parent component model \"%s\" is \"%s\" but not \"pesudo_coupled_system\". Please check the model code related to the annotation \"%s\". Please verify", 
                         comp_name, parent->comp_name, parent->comp_type, annotation);
/*
    if (parent != NULL && words_are_the_same(parent->comp_type, COMP_TYPE_PSEUDO_COUPLED) && parent->children.size() > 0)
        EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                         "Error happens when calling the API \"CCPL_register_component\" to register a component model \"%s\" at the model code with the annotation \"%s\": its parent \"%s\" is a coupled system (type is \"pesudo_coupled_system\") that can only have one child in a process while it already has a child \"%s\" (registerd at the model code with the annotation \"%s\"). Please verify.", 
                         comp_name, annotation, parent->comp_name, parent->children[0]->comp_name, parent->children[0]->get_annotation_start());
*/
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_rank(comm_group, &current_proc_local_id) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_rank(MPI_COMM_WORLD, &current_proc_global_id) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_size(comm_group, &num_procs) == MPI_SUCCESS);
    processes_global_id = new int [num_procs];
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Allgather(&current_proc_global_id, 1, MPI_INT, processes_global_id, 1, MPI_INT, comm_group) == MPI_SUCCESS);
    for (i = 0; i < num_procs; i ++)
        local_processes_global_ids.push_back(processes_global_id[i]);
    delete [] processes_global_id;

    if (parent != NULL) {
        for (i = 0; i < local_processes_global_ids.size(); i ++) {
            for (j = 0; j < parent->local_processes_global_ids.size(); j ++)
                if (local_processes_global_ids[i] == parent->local_processes_global_ids[j])
                    break;
            if (current_proc_local_id == 0)
                EXECUTION_REPORT(REPORT_ERROR,-1, j < parent->local_processes_global_ids.size(), 
                                 "The processes of component \"%s\" must be a subset of the processes of its parent component \"%s\". Please check the model code related to the annotations \"%s\" and \"%s\"", 
                                 comp_name, parent->get_comp_name(), annotation_start, parent->annotation_start);
        }
        parent->children.push_back(this);
    }
    
    if (ancestor != NULL) {
        if (is_real_component_model()) {
            sprintf(working_dir, "%s/CCPL_dir/run/data/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type);
            create_directory(working_dir, comm_group, get_current_proc_local_id() == 0, false);
            sprintf(working_dir, "%s/CCPL_dir/run/data/%s/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type, full_name);
            create_directory(working_dir, comm_group, get_current_proc_local_id() == 0, false);
        }
        sprintf(dir, "%s/CCPL_dir/run/CCPL_logs/by_components/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type);
        create_directory(dir, comm_group, get_current_proc_local_id() == 0, false);
        sprintf(dir, "%s/CCPL_dir/run/model_logs/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type);
        create_directory(dir, comm_group, get_current_proc_local_id() == 0, false);
        sprintf(dir, "%s/CCPL_dir/run/CCPL_logs/by_components/%s/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type,full_name);
        create_directory(dir, comm_group, get_current_proc_local_id() == 0, false);        
        sprintf(dir, "%s/CCPL_dir/run/model_logs/%s/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type,full_name);
        create_directory(dir, comm_group, get_current_proc_local_id() == 0, false);
        MPI_Barrier(get_comm_group());
        sprintf(comp_ccpl_log_file_name, "%s/CCPL_dir/run/CCPL_logs/by_components/%s/%s/%s.CCPL.log.%d", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type, full_name, get_comp_name(), get_current_proc_local_id());
        sprintf(comp_model_log_file_name, "%s/CCPL_dir/run/model_logs/%s/%s/%s.log.%d", comp_comm_group_mgt_mgr->get_root_working_dir(), comp_type, full_name, get_comp_name(), get_current_proc_local_id());
    }

    if (parent != NULL)
        check_API_parameter_string(parent->get_comp_id(), API_ID_COMP_MGT_REG_COMP, comm_group, "registering a component model", comp_type, "\"comp_type\"", annotation);
    else check_API_parameter_string(-1, API_ID_COMP_MGT_REG_COMP, comm_group, "registering a component model", comp_type, "\"comp_type\"", annotation);

    if (current_proc_local_id == 0) {
        char XML_file_name[NAME_STR_SIZE];
        sprintf(XML_file_name, "%s/%s.basic_info.xml", comp_comm_group_mgt_mgr->get_components_processes_dir(), full_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, !does_file_exist(XML_file_name), 
                         "Error happens when registering a component model \"%s\": another componet model with the same name has already been registered. Please check the model code related to the annotations \"%s\"", 
                         full_name, annotation);
        TiXmlDocument *XML_file = new TiXmlDocument;
        TiXmlDeclaration *XML_declaration = new TiXmlDeclaration(("1.0"),(""),(""));
        EXECUTION_REPORT(REPORT_ERROR, -1, XML_file != NULL, "Software error: cannot create an xml file");
        XML_file->LinkEndChild(XML_declaration);
        TiXmlElement *root_element = new TiXmlElement("Component");
        XML_file->LinkEndChild(root_element);
        write_node_into_XML(root_element);
        EXECUTION_REPORT(REPORT_ERROR, -1, XML_file->SaveFile(XML_file_name), "Software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node: fail to write the XML file %s", XML_file_name);
        delete XML_file;
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering the component model \%s\"", full_name);
}


Comp_comm_group_mgt_node::Comp_comm_group_mgt_node(TiXmlElement *XML_element, const char *specified_full_name, const char *XML_file_name)
{
    int line_number;


    comp_id = -1;
    temp_array_buffer = NULL;
    proc_latest_model_time = NULL;
    comp_model_log_file_device = -1;
    performance_timing_mgr = NULL;
    log_buffer = NULL;
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(XML_element->Value(), "Online_Model"), "Software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node: wrong element name");
    const char *XML_comp_name = get_XML_attribute(comp_id, CCPL_NAME_STR_LEN, XML_element, "comp_name", XML_file_name, line_number, "the name of the component model", "internal configuration file of component information", true);
    const char *XML_full_name = get_XML_attribute(comp_id, 512, XML_element, "full_name", XML_file_name, line_number, "the full name of the component model", "internal configuration file of component information", true);
    const char *XML_comp_type = get_XML_attribute(comp_id, 512, XML_element, "comp_type", XML_file_name, line_number, "the type of the component model", "internal configuration file of component information", true);
    const char *XML_enabled_in_parent_coupling_generation = get_XML_attribute(comp_id, CCPL_NAME_STR_LEN, XML_element, "enabled_in_parent_coupling_generation", XML_file_name, line_number, "enabled_in_parent_coupling_generation", "internal configuration file of component information", true);
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(XML_enabled_in_parent_coupling_generation,"true") || words_are_the_same(XML_enabled_in_parent_coupling_generation,"false"), "software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node: XML");
    if (words_are_the_same(XML_enabled_in_parent_coupling_generation,"true"))
        enabled_in_parent_coupling_generation = true;
    else enabled_in_parent_coupling_generation = false;
    const char *XML_processes = get_XML_attribute(comp_id, -1, XML_element, "processes", XML_file_name, line_number, "global IDs of the processes of the component model", "internal configuration file of component information", true);
    strcpy(this->comp_name, XML_comp_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(specified_full_name, XML_full_name), "Software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node: the full name specified is different from the full name in XML file %s: %s vs %s", XML_file_name, specified_full_name, XML_full_name);
    strcpy(this->full_name, XML_full_name);
    strcpy(this->comp_type, XML_comp_type);
    int segment_start, segment_end;
    for (int i = 1; i < strlen(XML_processes)+1; i ++) {
        if (XML_processes[i-1] == ' ') {
            segment_start = XML_processes[i]-'0';
            segment_end = -1;
            EXECUTION_REPORT(REPORT_ERROR, -1, segment_start >= 0 && segment_start <= 9, "Software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node: wrong format");
        }
        else if (XML_processes[i-1] == '~') {
            segment_end = XML_processes[i]-'0';
            EXECUTION_REPORT(REPORT_ERROR, -1, segment_end >= 0 && segment_end <= 9, "Software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node: wrong format");
        }
        else if (XML_processes[i] == ' ' || XML_processes[i] == '\0') {
            if (segment_end == -1)
                local_processes_global_ids.push_back(segment_start);
            else {
                for (int j = segment_start; j <= segment_end; j ++)
                    local_processes_global_ids.push_back(j);
            }
        }
        else if (XML_processes[i] != '~') {
            int digit = XML_processes[i] - '0';            
            EXECUTION_REPORT(REPORT_ERROR, -1, digit >= 0 && digit <= 9, "Software error in Comp_comm_group_mgt_node::Comp_comm_group_mgt_node: wrong format");
            if (segment_end == -1)
                segment_start = segment_start * 10 + digit;
            else segment_end = segment_end * 10 + digit;
        }
    }

    current_proc_local_id = -1;
    parent = NULL;
    restart_mgr = NULL;
}


void Comp_comm_group_mgt_node::transform_node_into_array()
{
    int num_procs, proc_id, num_children;


    num_procs = local_processes_global_ids.size();
    for (int i = 0; i < num_procs; i ++) {
        proc_id = local_processes_global_ids[i];
        write_data_into_array_buffer(&proc_id, sizeof(int), &temp_array_buffer, buffer_max_size, buffer_content_size);
    }
    write_data_into_array_buffer(&num_procs, sizeof(int), &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&enabled_in_parent_coupling_generation, sizeof(bool), &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&comm_group, sizeof(MPI_Comm), &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(comp_type, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(comp_name, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(full_name, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(comp_model_log_file_name, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(comp_ccpl_log_file_name, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(working_dir, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(annotation_start, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(annotation_end, NAME_STR_SIZE, &temp_array_buffer, buffer_max_size, buffer_content_size);
}


void Comp_comm_group_mgt_node::write_node_into_XML(TiXmlElement *parent_element)
{
    int i, num_segments;
    int *segments_start, *segments_end;
    TiXmlElement * current_element;
    char *string;

    
    current_element = new TiXmlElement("Online_Model");
    parent_element->LinkEndChild(current_element);
    current_element->SetAttribute("comp_name", comp_name);
    current_element->SetAttribute("full_name", full_name);
    current_element->SetAttribute("comp_type", comp_type);
    if (parent != NULL)
        current_element->SetAttribute("parent_full_name", parent->get_comp_full_name()); 
    else current_element->SetAttribute("parent_full_name", "NULL"); 
    
    if (enabled_in_parent_coupling_generation)
        current_element->SetAttribute("enabled_in_parent_coupling_generation", "true");
    else current_element->SetAttribute("enabled_in_parent_coupling_generation", "false");

    segments_start = new int [local_processes_global_ids.size()];
    segments_end = new int [local_processes_global_ids.size()];
    segments_start[0] = local_processes_global_ids[0];
    for (i = 1, num_segments = 1; i < local_processes_global_ids.size(); i ++) {
        if (local_processes_global_ids[i] != local_processes_global_ids[i-1]+1) {
            segments_end[num_segments-1] = local_processes_global_ids[i-1];
            segments_start[num_segments] = local_processes_global_ids[i];
            num_segments ++;
        }
    }
    segments_end[num_segments-1] = local_processes_global_ids[local_processes_global_ids.size()-1];
    string = new char [num_segments*(8*2+1)];
    string[0] = '\0';
    for (i = 0; i < num_segments; i ++) {
        if (segments_start[i] != segments_end[i])
            sprintf(string+strlen(string), " %d~%d", segments_start[i], segments_end[i]);
        else sprintf(string+strlen(string), " %d", segments_start[i]);
    }

    current_element->SetAttribute("processes", string);
    
    delete [] segments_start;
    delete [] segments_end;
    delete [] string;
}


void Comp_comm_group_mgt_node::update_child(const Comp_comm_group_mgt_node *child_old, Comp_comm_group_mgt_node *child_new)
{
    int i;


    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(child_old->full_name, child_new->full_name), "software error in Comp_comm_group_mgt_node::update_child: children names are not the same");
    for (i = 0; i < children.size(); i ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, children[i]->parent == this, "software error in Comp_comm_group_mgt_node::update_child: wrong parent1");
        EXECUTION_REPORT(REPORT_ERROR, -1, children[i]->parent != children[i], "software error in Comp_comm_group_mgt_node::update_child: wrong parent2");
        if (children[i] == child_old) {
            EXECUTION_REPORT_LOG(REPORT_LOG, child_old->comp_id, true, "Link the parent of component \"%s\" to \"%s\"", child_old->full_name, child_old->parent->get_full_name());
            child_new->parent = this;
            children[i] = child_new;
            break;
        }
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, i < children.size(), "software error in Comp_comm_group_mgt_node::update_child");
}


void Comp_comm_group_mgt_node::transfer_data_buffer(Comp_comm_group_mgt_node *new_node)
{
    if (new_node->temp_array_buffer != NULL)
        delete [] new_node->temp_array_buffer;

    new_node->temp_array_buffer = this->temp_array_buffer;
    new_node->buffer_content_iter = this->buffer_content_iter;
    new_node->buffer_content_size = this->buffer_content_size;
    new_node->buffer_max_size = this->buffer_max_size;
    this->temp_array_buffer = NULL;
}


void Comp_comm_group_mgt_node::confirm_coupling_configuration_active(int API_id, bool require_real_model, const char *annotation)
{
    char API_label[NAME_STR_SIZE]; 

    get_API_hint(comp_id, API_id, API_label);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, !definition_finalized, 
                     "ERROR happens when calling the API \"%s\" at the model code with the annotation \"%s\": the coupling configuration stage of the corresponding component model \"%s\" has been ended at the model code with the annotation \"%s\"", 
                     API_label, annotation, comp_name, get_annotation_end());
    if (require_real_model)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, is_real_component_model(), 
                         "ERROR happens when calling the API \"%s\" at the model code with the annotation \"%s\": the corresponding component model \"%s\" cannot handle coupling configuration because it is a pesudo coupled system (its type is \"pesudo_coupled_system\"). Please verify.", 
                         API_label, annotation, comp_name);
}


int Comp_comm_group_mgt_node::get_local_proc_global_id(int local_indx)
{
    if (local_indx < local_processes_global_ids.size())
        return local_processes_global_ids[local_indx];

    EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Comp_comm_group_mgt_node::get_local_proc_global_id");
    
    return -1;
}


bool Comp_comm_group_mgt_node::is_real_component_model()
{ 
    return !words_are_the_same(comp_type, COMP_TYPE_PSEUDO_COUPLED) && !words_are_the_same(comp_type, COMP_TYPE_ROOT); 
}


bool Comp_comm_group_mgt_node::have_local_process(int local_proc_global_id)
{
    for (int i = 0; i < local_processes_global_ids.size(); i ++)
        if (local_proc_global_id == local_processes_global_ids[i])
            return true;
        
    return false;
}


void Comp_comm_group_mgt_node::allocate_proc_latest_model_time()
{
    if (proc_latest_model_time != NULL)
        return;
    
    EXECUTION_REPORT(REPORT_ERROR, -1, get_num_procs() > 0, "Software error in Comp_comm_group_mgt_node::allocate_proc_latest_model_time");
    proc_latest_model_time = new long [get_num_procs()];
    for (int i = 0; i < get_num_procs(); i ++)
        proc_latest_model_time[i] = -1;
}


void Comp_comm_group_mgt_node::set_current_proc_current_time(int days, int second)
{
    allocate_proc_latest_model_time();
    proc_latest_model_time[current_proc_local_id] = ((long)days)*((long)100000) + (long)second;
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "set_proc_latest_model_time %ld: %d %d", proc_latest_model_time[current_proc_local_id], days, second);
}


void Comp_comm_group_mgt_node::set_proc_latest_model_time(int proc_id, long model_time)
{
    allocate_proc_latest_model_time();
    EXECUTION_REPORT(REPORT_ERROR, -1, proc_id >= 0 && proc_id < get_num_procs(), "Software error in set_proc_latest_model_time: wrong proc id: %d vs %d", proc_id, get_num_procs());

    if (model_time > proc_latest_model_time[proc_id])
        proc_latest_model_time[proc_id] = model_time;
}


long Comp_comm_group_mgt_node::get_proc_latest_model_time(int proc_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, proc_id >= 0 && proc_id < get_num_procs(), "Software error in get_proc_latest_model_time: wrong proc id");
    return proc_latest_model_time[proc_id];
}


void Comp_comm_group_mgt_node::get_all_descendant_real_comp_fullnames(int top_comp_id, std::vector<const char*> &all_descendant_real_comp_fullnames, char **temp_array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    char *local_temp_array_buffer = NULL, *gather_temp_array_buffer = NULL;
    long local_buffer_max_size, local_buffer_content_size = 0, gather_buffer_content_size = 0;
    

    for (int i = 0; i < children.size(); i ++)
        if (children[i]->enabled_in_parent_coupling_generation)
            children[i]->get_all_descendant_real_comp_fullnames(top_comp_id, all_descendant_real_comp_fullnames, &local_temp_array_buffer, local_buffer_max_size, local_buffer_content_size);

    if (is_real_component_model() && current_proc_local_id == 0)
        dump_string(full_name, -1, &local_temp_array_buffer, local_buffer_max_size, local_buffer_content_size);

    if (current_proc_local_id != -1) {
        int *all_array_size = new int [get_num_procs()];
        gather_array_in_one_comp(get_num_procs(), current_proc_local_id, local_temp_array_buffer, local_buffer_content_size, sizeof(char), all_array_size, (void**)(&gather_temp_array_buffer), gather_buffer_content_size, comm_group);
        if (current_proc_local_id == 0)
            write_data_into_array_buffer(gather_temp_array_buffer, gather_buffer_content_size, temp_array_buffer, buffer_max_size, buffer_content_size);
        delete [] all_array_size;
    }

    if (local_temp_array_buffer != NULL)
        delete [] local_temp_array_buffer;
    if (gather_temp_array_buffer != NULL)
        delete [] gather_temp_array_buffer;

    if (comp_id == top_comp_id) {
        bcast_array_in_one_comp(current_proc_local_id, temp_array_buffer, buffer_content_size, comm_group);
        char temp_full_name[NAME_STR_SIZE];
        long str_size;
        while (buffer_content_size > 0) {
            load_string(temp_full_name, str_size, NAME_STR_SIZE, *temp_array_buffer, buffer_content_size, "C-Coupler internal");
            all_descendant_real_comp_fullnames.push_back(strdup(temp_full_name));
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, buffer_content_size == 0, "Software error in Comp_comm_group_mgt_node::get_all_descendant_real_comp_fullnames");
        if (*temp_array_buffer != NULL) {
            delete [] *temp_array_buffer;
            *temp_array_buffer = NULL;
        }
        if (current_proc_local_id == 0)
            for (int i = 0; i < all_descendant_real_comp_fullnames.size(); i ++)
                EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "%x comp %d for internal generation %s\n", top_comp_id, i, all_descendant_real_comp_fullnames[i]);
    }
}


int Comp_comm_group_mgt_node::open_comp_model_log_file(int *log_file_device_id)
{
    int log_file_opened;

    
    if (comp_model_log_file_device != -1)
        log_file_opened = 1;
    else log_file_opened = 0;

    comp_model_log_file_device = 100+(comp_id&TYPE_ID_SUFFIX_MASK);
    *log_file_device_id = comp_model_log_file_device;
    
    return log_file_opened;
}


int Comp_comm_group_mgt_node::get_min_remote_lag_seconds()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_proc_local_id != -1, "Software error in Comp_comm_group_mgt_node::get_min_remote_lag_seconds");    
    return min_remote_lag_seconds;
}


int Comp_comm_group_mgt_node::get_max_remote_lag_seconds()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_proc_local_id != -1, "Software error in Comp_comm_group_mgt_node::get_max_remote_lag_seconds");    
    return max_remote_lag_seconds;
}


void Comp_comm_group_mgt_node::update_min_max_remote_lag_seconds(int remote_lag_seconds)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_proc_local_id != -1, "Software error in Comp_comm_group_mgt_node::update_min_max_remote_lag_seconds");
    if (remote_lag_seconds < this->min_remote_lag_seconds)
        this->min_remote_lag_seconds = remote_lag_seconds;
    if (remote_lag_seconds < this->max_remote_lag_seconds)
        this->max_remote_lag_seconds = remote_lag_seconds;
}


void Comp_comm_group_mgt_node::output_log(const char *log_string, bool flush_log_file)
{
    output_CCPL_log(log_string, comp_ccpl_log_file_name, &log_buffer, log_buffer_content_size, flush_log_file);
}


void Comp_comm_group_mgt_node::reset_local_node_id(int new_id) 
{ 
    comp_id = new_id; 
    EXECUTION_REPORT(REPORT_ERROR, -1, restart_mgr == NULL, "Software error in Comp_comm_group_mgt_node::reset_local_node_id");
}



Comp_comm_group_mgt_mgr::Comp_comm_group_mgt_mgr(const char *executable_name)
{
    int i, j, num_procs, proc_id;
    char temp_string[NAME_STR_SIZE];
    std::vector<char*> unique_executable_name;


    unique_comp_id_indx = 0;
    global_node_array.clear();
    global_node_root = NULL;
    definition_finalized = false;
    CCPL_platform_log_dir[0] = '\0';
    log_buffer = NULL; 
    EXECUTION_REPORT(REPORT_ERROR, -1, getcwd(root_working_dir,NAME_STR_SIZE) != NULL, 
                     "Cannot get the current working directory for running the model");

    for (i = strlen(executable_name)-1; i >= 0; i --)
        if (executable_name[i] == '/')
            break;
    i ++;
    EXECUTION_REPORT(REPORT_ERROR,-1, i < strlen(executable_name), "Software error1 in Comp_comm_group_mgt_mgr::Comp_comm_group_mgt_mgr");
    strcpy(this->executable_name, executable_name+i);    
    EXECUTION_REPORT(REPORT_ERROR, -1, !words_are_the_same(this->executable_name, "all"), "Error happens when using the executable \"%s\" to run the coupled system: the name of any executable cannot be all");

    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_rank(MPI_COMM_WORLD, &current_proc_global_id) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_size(MPI_COMM_WORLD, &num_total_global_procs) == MPI_SUCCESS);    
	if (current_proc_global_id == 0) {
		sprintf(temp_string, "%s/CCPL_dir/config", root_working_dir);
		EXECUTION_REPORT(REPORT_ERROR, -1, does_file_exist(temp_string), "Fail to initialize C-Coupler: the directory (\"%s\") for configuration files of C-Coupler does not exit.", temp_string);
	}
    sprintf(temp_string, "%s/CCPL_dir/run", root_working_dir);
    create_directory(temp_string, MPI_COMM_WORLD, current_proc_global_id == 0, false);
    sprintf(temp_string, "%s/CCPL_dir/run/CCPL_logs", root_working_dir);
    create_directory(temp_string, MPI_COMM_WORLD, current_proc_global_id == 0, true);
    sprintf(temp_string, "%s/CCPL_dir/run/model_logs", root_working_dir);
    create_directory(temp_string, MPI_COMM_WORLD, current_proc_global_id == 0, true);
    sprintf(temp_string, "%s/CCPL_dir/run/CCPL_logs/by_components", root_working_dir);
    create_directory(temp_string, MPI_COMM_WORLD, current_proc_global_id == 0, false);
    sprintf(temp_string, "%s/CCPL_dir/run/data", root_working_dir);
    create_directory(temp_string, MPI_COMM_WORLD, current_proc_global_id == 0, false);
    sprintf(temp_string, "%s/CCPL_dir/run/data/all", root_working_dir);
    create_directory(temp_string, MPI_COMM_WORLD, current_proc_global_id == 0, false);
    sprintf(internal_H2D_grids_dir, "%s/CCPL_dir/run/data/all/internal_H2D_grids", root_working_dir);
    create_directory(internal_H2D_grids_dir, MPI_COMM_WORLD, current_proc_global_id == 0, true);
    sprintf(internal_remapping_weights_dir, "%s/CCPL_dir/run/data/all/internal_remapping_weights", root_working_dir);
    create_directory(internal_remapping_weights_dir, MPI_COMM_WORLD, current_proc_global_id == 0, false);
    sprintf(components_processes_dir, "%s/CCPL_dir/run/data/all/components_processes", root_working_dir);
    create_directory(components_processes_dir, MPI_COMM_WORLD, current_proc_global_id == 0, true);
    sprintf(components_exports_dir, "%s/CCPL_dir/run/data/all/components_exports", root_working_dir);
    create_directory(components_exports_dir, MPI_COMM_WORLD, current_proc_global_id == 0, true);
    sprintf(active_coupling_connections_dir, "%s/CCPL_dir/run/data/all/active_coupling_connections", root_working_dir);
    create_directory(active_coupling_connections_dir, MPI_COMM_WORLD, current_proc_global_id == 0, true);    
    sprintf(comps_ending_config_status_dir, "%s/CCPL_dir/run/data/all/comps_ending_config_status", root_working_dir);
    create_directory(comps_ending_config_status_dir, MPI_COMM_WORLD, current_proc_global_id == 0, true);
    sprintf(restart_common_dir, "%s/CCPL_dir/run/data/all/restart", root_working_dir);
    create_directory(restart_common_dir, MPI_COMM_WORLD, current_proc_global_id == 0, false);
    sprintf(runtime_config_root_dir, "%s/CCPL_dir/config", root_working_dir);
    root_comp_config_dir[0] = '\0';
    sprintf(temp_string, "%s/CCPL_dir/run/CCPL_logs/by_executables", root_working_dir);
    create_directory(temp_string, MPI_COMM_WORLD, current_proc_global_id == 0, false);    
    sprintf(exe_log_file_name, "%s/CCPL_dir/run/CCPL_logs/by_executables/%s/%s.CCPL.log.%d", root_working_dir, this->executable_name, this->executable_name, current_proc_global_id);
    MPI_Barrier(MPI_COMM_WORLD);
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_size(MPI_COMM_WORLD, &num_procs) == MPI_SUCCESS);
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_rank(MPI_COMM_WORLD, &proc_id) == MPI_SUCCESS);
    char *all_executable_name, dir[NAME_STR_SIZE];
    if (proc_id == 0)
        all_executable_name = new char [NAME_STR_SIZE*num_procs];
    EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Gather((char*)this->executable_name, NAME_STR_SIZE, MPI_CHAR, all_executable_name, NAME_STR_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD) == MPI_SUCCESS);
    if (proc_id == 0) {
        unique_executable_name.push_back(all_executable_name);
        for (i = 1; i < num_procs; i ++) {
            for (j = 0; j < unique_executable_name.size(); j ++)
                if (words_are_the_same(unique_executable_name[j], all_executable_name+i*NAME_STR_SIZE)) {
                    break;
                }
            if (j == unique_executable_name.size())
                unique_executable_name.push_back(all_executable_name+i*NAME_STR_SIZE);
        }
        for (int j = 0; j < unique_executable_name.size(); j ++) {
            sprintf(dir, "%s/CCPL_dir/run/CCPL_logs/by_executables/%s", root_working_dir, unique_executable_name[j]);
            create_directory(dir, MPI_COMM_WORLD, true, false);
        }
        delete [] all_executable_name;    
    }

    MPI_Barrier(MPI_COMM_WORLD);
}


Comp_comm_group_mgt_mgr::~Comp_comm_group_mgt_mgr()
{
    if (log_buffer != NULL) {
        output_log("", true);
        delete [] log_buffer;
    }

    for (int i = 0; i < global_node_array.size(); i ++)
        delete global_node_array[i];

    for (int i = 0; i < root_comps_full_names.size(); i ++)
        delete root_comps_full_names[i];
}


void Comp_comm_group_mgt_mgr::transform_global_node_tree_into_array(Comp_comm_group_mgt_node *current_global_node, Comp_comm_group_mgt_node **all_global_nodes, int &global_node_id)
{
    all_global_nodes[global_node_id++] = current_global_node;
    for (int i = 0; i < current_global_node->get_num_children(); i ++)
        transform_global_node_tree_into_array(current_global_node->get_child(i), all_global_nodes, global_node_id);
}


bool Comp_comm_group_mgt_mgr::is_legal_local_comp_id(int local_comp_id, bool is_external_call)
{
    int i;

    
    if ((local_comp_id&TYPE_ID_PREFIX_MASK) != TYPE_COMP_LOCAL_ID_PREFIX)
        return false;
    
    for (i = 0; i < global_node_array.size(); i ++)
        if (local_comp_id == global_node_array[i]->get_comp_id())
            break;

    if (i == global_node_array.size())
        return false;

    if (is_external_call)
         return !does_comp_name_include_reserved_prefix(global_node_array[i]->get_comp_name());
    else return true;
}


int Comp_comm_group_mgt_mgr::register_component(const char *comp_name, const char *comp_type, MPI_Comm &comm, int parent_local_id, bool enabled_in_parent_coupling_gen, int change_dir, const char *annotation)
{
    int i;
    Comp_comm_group_mgt_node *root_local_node, *new_comp;
    MPI_Comm global_comm = MPI_COMM_WORLD;
    char hint[NAME_STR_SIZE];
    
    
    if (definition_finalized)
        EXECUTION_REPORT(REPORT_ERROR, -1, !definition_finalized, 
                         "Cannot register component \"%s\" at the model code with the annotation \"%s\" because the stage of registering coupling configurations of the whole coupled model has been ended at the model code with the annotation \"%s\"", 
                         comp_name, annotation, global_node_array[0]->get_annotation_end());
    
    for (i = 0; i < global_node_array.size(); i ++)
        if (words_are_the_same(global_node_array[i]->get_comp_name(), comp_name))
            break;
        
    if (i < global_node_array.size())
        EXECUTION_REPORT(REPORT_ERROR, -1, i == global_node_array.size(), "Error happens when registering a component model named \"%s\" at the model code with the annotation \"%s\": a component model with the same name has already been registered before, at the model code with the annotation \"%s\". Please note that, any two component models on the same MPI process cannot have the same name", comp_name, annotation, global_node_array[i]->get_annotation_start());

    if (parent_local_id == -1) {
        root_local_node = new Comp_comm_group_mgt_node(COMP_TYPE_ROOT, COMP_TYPE_ROOT, (unique_comp_id_indx++)|TYPE_COMP_LOCAL_ID_PREFIX, NULL, global_comm, enabled_in_parent_coupling_gen, annotation);
        global_node_array.push_back(root_local_node);
        global_node_root = root_local_node;
        new_comp = new Comp_comm_group_mgt_node(comp_name, comp_type, (unique_comp_id_indx++)|TYPE_COMP_LOCAL_ID_PREFIX, root_local_node, comm, enabled_in_parent_coupling_gen, annotation);
        global_node_array.push_back(new_comp);
    }
    else {
        Comp_comm_group_mgt_node *parent_comp_node = search_global_node(parent_local_id);
        EXECUTION_REPORT(REPORT_ERROR, -1, !parent_comp_node->is_definition_finalized(), 
                         "Cannot register component \"%s\" at the model code with the annotation \"%s\" because the registration corresponding to the parent \"%s\" has been ended at the model code with the annotation \"%s\"", 
                         comp_name, annotation, parent_comp_node->get_comp_name(), parent_comp_node->get_annotation_end()); // add debug information
        new_comp = new Comp_comm_group_mgt_node(comp_name, comp_type, (unique_comp_id_indx++)|TYPE_COMP_LOCAL_ID_PREFIX, parent_comp_node, comm, enabled_in_parent_coupling_gen, annotation);
        global_node_array.push_back(new_comp);
    }
    sprintf(hint, "regietering a component model \"%s\"", comp_name);
    check_API_parameter_bool(new_comp->get_comp_id(), API_ID_COMP_MGT_REG_COMP, new_comp->get_comm_group(), hint, enabled_in_parent_coupling_gen, "enabled_in_parent_coupling_gen", annotation);

    if (parent_local_id == -1) {
        char *temp_array_buffer = NULL, *gather_array_buffer = NULL;
        long local_buffer_max_size, local_buffer_content_size = 0, gather_buffer_content_size = 0, str_size;
        bool temp_enabled_in_parent_coupling_gen;
        char root_comp_full_name[NAME_STR_SIZE];
        if (new_comp->get_current_proc_local_id() == 0) {
            dump_string(new_comp->get_comp_full_name(), -1, &temp_array_buffer, local_buffer_max_size, local_buffer_content_size);
            write_data_into_array_buffer(&enabled_in_parent_coupling_gen, sizeof(bool), &temp_array_buffer, local_buffer_max_size, local_buffer_content_size);
        }
        gather_array_in_one_comp(num_total_global_procs, current_proc_global_id, temp_array_buffer, local_buffer_content_size, sizeof(char), NULL, (void**)(&gather_array_buffer), gather_buffer_content_size, MPI_COMM_WORLD);
        bcast_array_in_one_comp(current_proc_global_id, &gather_array_buffer, gather_buffer_content_size, MPI_COMM_WORLD);
        while(gather_buffer_content_size > 0) {
            read_data_from_array_buffer(&temp_enabled_in_parent_coupling_gen, sizeof(bool), gather_array_buffer, gather_buffer_content_size, true);
            load_string(root_comp_full_name, str_size, NAME_STR_SIZE, gather_array_buffer, gather_buffer_content_size, NULL);
            root_comps_enabled_in_parent_coupling_generation.push_back(temp_enabled_in_parent_coupling_gen);
            root_comps_full_names.push_back(strdup(root_comp_full_name));
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, gather_buffer_content_size == 0, "Software error in Comp_comm_group_mgt_mgr::register_component");
        if (temp_array_buffer != NULL)
            delete [] temp_array_buffer;
        if (gather_array_buffer != NULL)
            delete [] gather_array_buffer;

        sprintf(root_comp_config_dir, "%s/%s", runtime_config_root_dir, new_comp->get_comp_name());
        CCPL_platform_log_dir[0] = '\0';
        check_API_parameter_int(new_comp->get_comp_id(), API_ID_COMP_MGT_REG_COMP, new_comp->get_comm_group(), "registering a component model", change_dir, "change_dir", annotation);
        if (change_dir == 1) {
            char new_dir[NAME_STR_SIZE];
            sprintf(new_dir, "%s/run/%s/%s/data", root_working_dir, new_comp->get_comp_type(), new_comp->get_comp_name());
            sprintf(CCPL_platform_log_dir, "%s/run/%s/%s/run_logs", root_working_dir, new_comp->get_comp_type(), new_comp->get_comp_name());
            DIR *dir=opendir(new_dir);
            EXECUTION_REPORT(REPORT_ERROR, new_comp->get_comp_id(), dir != NULL, "Fail to change working directory for the first active component model \"%s\": the directory \"%s\" does not exist.", new_comp->get_comp_name(), new_dir);
            chdir(new_dir);
            EXECUTION_REPORT_LOG(REPORT_LOG, new_comp->get_comp_id(), true, "change working directory to \"%s\"", new_dir);
        }
        original_grid_mgr->initialize_CoR_grids();
    }

    EXECUTION_REPORT(REPORT_PROGRESS, new_comp->get_comp_id(), true, "The component model \"%s\" is successfully registered at the model code with the annotation \"%s\".", new_comp->get_full_name(), annotation);

    EXECUTION_REPORT(REPORT_ERROR, new_comp->get_comp_id(), !does_comp_name_include_reserved_prefix(comp_name), "Error happens when registering a component model \"%s\": its name should not include the prefix \"%s\", \"%s\", \"%s\" and \"%s\". Please verify the model code with the annotation \"%s\"", comp_name, COMP_TYPE_ROOT, DATAINST_NAME_PREFIX, DATAMODEL_NAME_PREFIX, ALGMODEL_NAME_PREFIX, annotation);

    return new_comp->get_local_node_id();
}


Comp_comm_group_mgt_node *Comp_comm_group_mgt_mgr::get_global_node_of_local_comp(int local_comp_id, bool is_external_call, const char *annotation)
{    
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_legal_local_comp_id(local_comp_id, is_external_call), "The id of component is wrong when getting the management node of a component. Please check the model code with the annotation \"%s\"", annotation); 

    return search_global_node(local_comp_id);
}


MPI_Comm Comp_comm_group_mgt_mgr::get_comm_group_of_local_comp(int local_comp_id, const char *annotation)
{
    return get_global_node_of_local_comp(local_comp_id,false,annotation)->get_comm_group();
}


void Comp_comm_group_mgt_mgr::get_output_data_file_header(int comp_id, char *data_file_header)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_legal_local_comp_id(comp_id, true), "software error in Comp_comm_group_mgt_mgr::get_data_file_header");
    sprintf(data_file_header, "%s/%s", search_global_node(comp_id)->get_working_dir(), search_global_node(comp_id)->get_comp_name());
}


const char *Comp_comm_group_mgt_mgr::get_comp_ccpl_log_file_name(int comp_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_legal_local_comp_id(comp_id, true), "software error in Comp_comm_group_mgt_mgr::get_comp_ccpl_log_file_name");
    return search_global_node(comp_id)->get_comp_ccpl_log_file_name();
}


const char *Comp_comm_group_mgt_mgr::get_comp_model_log_file(int comp_id, int &device_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_legal_local_comp_id(comp_id, true), "software error in Comp_comm_group_mgt_mgr::get_comp_model_log_file");
    return search_global_node(comp_id)->get_comp_model_log_file(device_id);
}


Comp_comm_group_mgt_node *Comp_comm_group_mgt_mgr::search_comp_with_comp_name(const char *comp_name)
{
    for (int i = 0; i < global_node_array.size(); i ++)
        if (words_are_the_same(global_node_array[i]->get_comp_name(), comp_name) && global_node_array[i]->get_current_proc_local_id() != -1)
            return global_node_array[i];
        
    return NULL;    
}


void Comp_comm_group_mgt_mgr::check_validation()
{
    for (int i = 0; i < global_node_array.size(); i ++) {        
        EXECUTION_REPORT(REPORT_ERROR, -1, (global_node_array[i]->get_comp_id()&TYPE_ID_SUFFIX_MASK) == i, "Software error in Comp_comm_group_mgt_mgr::check_validation: wrong comp id does not match array index");
        EXECUTION_REPORT(REPORT_ERROR, -1, search_global_node(global_node_array[i]->get_comp_id()) == global_node_array[i], "Software error in Comp_comm_group_mgt_mgr::check_validation: wrong comp id");
        if (global_node_array[i]->get_parent() != NULL) {
            EXECUTION_REPORT(REPORT_ERROR, -1, (global_node_array[i]->get_comp_id()&TYPE_ID_SUFFIX_MASK) > (global_node_array[i]->get_parent()->get_comp_id()&TYPE_ID_SUFFIX_MASK), "Software error in Comp_comm_group_mgt_mgr::check_validation: wrong parent1: %s vs %s", global_node_array[i]->get_full_name(), global_node_array[i]->get_parent()->get_full_name());
            char full_name[NAME_STR_SIZE];
            if ( global_node_array[i]->get_parent()->is_real_component_model()) {
                sprintf(full_name, "%s@%s", global_node_array[i]->get_parent()->get_full_name(), global_node_array[i]->get_comp_name());
                EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(full_name, global_node_array[i]->get_full_name()), "Software error in Comp_comm_group_mgt_mgr::check_validation: wrong parent2: %s vs %s", global_node_array[i]->get_full_name(), global_node_array[i]->get_parent()->get_full_name());
            }
            EXECUTION_REPORT(REPORT_ERROR, -1, search_global_node(global_node_array[i]->get_parent()->get_full_name()) == global_node_array[i]->get_parent(), "Software error in Comp_comm_group_mgt_mgr::check_validation: wrong parent3: %s vs %s", global_node_array[i]->get_full_name(), global_node_array[i]->get_parent()->get_full_name());
        }    
    }    
}


Comp_comm_group_mgt_node *Comp_comm_group_mgt_mgr::search_global_node(const char *full_name)
{    
    for (int i = 0; i < global_node_array.size(); i ++)
        if (words_are_the_same(global_node_array[i]->get_full_name(), full_name)) {
            return global_node_array[i];
        }
        
    return NULL;    
}


Comp_comm_group_mgt_node *Comp_comm_group_mgt_mgr::search_global_node(int global_node_id)
{
    for (int i = 0; i < global_node_array.size(); i ++)
        if (global_node_array[i]->get_local_node_id() == global_node_id)
            return global_node_array[i];
        
    return NULL;
}


int Comp_comm_group_mgt_mgr::get_current_proc_id_in_comp(int comp_id, const char *annotation)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_legal_local_comp_id(comp_id, true), 
                     "The component id specified for getting the id of the current process is wrong. Please check the model code with the annotation %s.", 
                     annotation); 
    return search_global_node(comp_id)->get_current_proc_local_id();
}


int Comp_comm_group_mgt_mgr::get_num_proc_in_comp(int comp_id, const char *annotation)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_legal_local_comp_id(comp_id, true), 
                     "The component id specified for getting the the number of processes is wrong. Please check the model code with the annotation %s.", 
                     annotation); 
    EXECUTION_REPORT(REPORT_ERROR, -1, search_global_node(comp_id)->get_current_proc_local_id() != -1, 
                     "The component id specified for getting the the number of processes is wrong. Please check the model code with the annotation %s.", 
                     annotation);
    return search_global_node(comp_id)->get_num_procs();
}


void Comp_comm_group_mgt_mgr::confirm_coupling_configuration_active(int comp_id, int API_id, bool require_real_model, const char *annotation)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id, true), "software error in Comp_comm_group_mgt_mgr::confirm_coupling_configuration_active");
    get_global_node_of_local_comp(comp_id,true,"in Comp_comm_group_mgt_mgr::confirm_coupling_configuration_active")->confirm_coupling_configuration_active(API_id, require_real_model, annotation);
}


const int *Comp_comm_group_mgt_mgr::get_all_components_ids()
{
    int *all_components_ids = new int [global_node_array.size()];


    for (int i = 1; i < global_node_array.size(); i ++) {
        all_components_ids[i] = global_node_array[i]->get_comp_id();
    }
    all_components_ids[0] = global_node_array.size();

    return all_components_ids;
}


bool Comp_comm_group_mgt_mgr::has_comp_ended_configuration(const char *comp_full_name)
{
    char status_file_name[NAME_STR_SIZE];
    sprintf(status_file_name, "%s/%s.end", comps_ending_config_status_dir, comp_full_name);
    FILE *status_file = fopen(status_file_name, "r");
    if (status_file == NULL)
        return false;

    fclose(status_file);
    return true;
}


void Comp_comm_group_mgt_mgr::push_comp_node(Comp_comm_group_mgt_node *comp_node) 
{
    comp_node->reset_local_node_id((unique_comp_id_indx++)|TYPE_COMP_LOCAL_ID_PREFIX);
    global_node_array.push_back(comp_node); 
}


Comp_comm_group_mgt_node *Comp_comm_group_mgt_mgr::pop_comp_node()
{
    Comp_comm_group_mgt_node *top_comp_node = global_node_array[global_node_array.size()-1];
    global_node_array.erase(global_node_array.begin()+global_node_array.size()-1);
    return top_comp_node;
}


void Comp_comm_group_mgt_mgr::set_current_proc_current_time(int comp_id, int days, int second)
{
    get_global_node_of_local_comp(comp_id,true,"Comp_comm_group_mgt_mgr::set_current_proc_time")->set_current_proc_current_time(days, second);
}


Comp_comm_group_mgt_node *Comp_comm_group_mgt_mgr::load_comp_info_from_XML(int host_comp_id, const char *comp_full_name, MPI_Comm comm)
{
    char XML_file_name[NAME_STR_SIZE];
    TiXmlDocument *XML_file;
    int i;


    sprintf(XML_file_name, "%s/%s.basic_info.xml", comp_comm_group_mgt_mgr->get_components_processes_dir(), comp_full_name);
    XML_file = open_XML_file_to_read(host_comp_id, XML_file_name, comm, true);
    TiXmlElement *XML_element = XML_file->FirstChildElement();
    TiXmlElement *Online_Model = XML_element->FirstChildElement();
    Comp_comm_group_mgt_node *pesudo_comp_node = new Comp_comm_group_mgt_node(Online_Model, comp_full_name, XML_file_name);
    delete XML_file;
        
    return pesudo_comp_node;
}


void Comp_comm_group_mgt_mgr::get_root_comps_for_overall_coupling_generation(std::vector<const char *> &all_comp_fullnames_for_coupling_generation)
{
    for (int i = 0; i < root_comps_full_names.size(); i ++) 
        if (root_comps_enabled_in_parent_coupling_generation[i])
            all_comp_fullnames_for_coupling_generation.push_back(strdup(root_comps_full_names[i]));
}


bool Comp_comm_group_mgt_mgr::is_comp_type_coupled(int host_comp_id, const char *comp_type, const char *annotation)
{
    char comp_full_name[NAME_STR_SIZE];


    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(comp_type,COMP_TYPE_CPL) || words_are_the_same(comp_type,COMP_TYPE_ATM) || words_are_the_same(comp_type,COMP_TYPE_GLC) || words_are_the_same(comp_type,COMP_TYPE_ATM_CHEM) || words_are_the_same(comp_type,COMP_TYPE_OCN) || words_are_the_same(comp_type,COMP_TYPE_LND) || words_are_the_same(comp_type,COMP_TYPE_SEA_ICE) || words_are_the_same(comp_type,COMP_TYPE_WAVE) || words_are_the_same(comp_type,COMP_TYPE_RUNOFF), 
                     "ERROR happens when calling the API \"CCPL_is_comp_type_coupled\": the component type \"%s\" is unknown. Please verify the model code with the annotation \"%s\"", comp_type, annotation);

    DIR *cur_dir = opendir(comp_comm_group_mgt_mgr->get_components_processes_dir());
    struct dirent *ent = NULL;
    struct stat st;
    EXECUTION_REPORT(REPORT_ERROR, -1, cur_dir != NULL, "Comp_comm_group_mgt_mgr::is_comp_type_coupled");
    while ((ent = readdir(cur_dir)) != NULL) {
        stat(ent->d_name, &st);
        if (!(strlen(ent->d_name) > strlen(".basic_info.xml") && words_are_the_same(ent->d_name+strlen(ent->d_name)-strlen(".basic_info.xml"), ".basic_info.xml")))
            continue;
        strncpy(comp_full_name, ent->d_name, strlen(ent->d_name)-strlen(".basic_info.xml"));
        comp_full_name[strlen(ent->d_name)-strlen(".basic_info.xml")] = '\0';
        Comp_comm_group_mgt_node *temp_comp_node = load_comp_info_from_XML(host_comp_id, comp_full_name, get_comm_group_of_local_comp(host_comp_id, "is_comp_type_coupled"));
        if (words_are_the_same(temp_comp_node->get_comp_type(), comp_type)) {
            delete temp_comp_node;
            return true;
        }
        delete temp_comp_node;
    }
    
    return false;
}


void Comp_comm_group_mgt_mgr::output_log(const char *log_string, bool flush_log_file)
{
    output_CCPL_log(log_string, exe_log_file_name, &log_buffer, log_buffer_content_size, flush_log_file);
}


void Comp_comm_group_mgt_mgr::output_performance_timing()
{
    for (int i = 0; i < global_node_array.size(); i ++)
        if (global_node_array[i]->is_real_component_model() && global_node_array[i]->get_performance_timing_mgr() != NULL)
            global_node_array[i]->get_performance_timing_mgr()->performance_timing_output();
}


bool Comp_comm_group_mgt_mgr::does_comp_name_include_reserved_prefix(const char *comp_name)
{
    return strncmp(comp_name, COMP_TYPE_ROOT, strlen(COMP_TYPE_ROOT)) == 0 || 
           strncmp(comp_name, DATAINST_NAME_PREFIX, strlen(DATAINST_NAME_PREFIX)) == 0 || 
           strncmp(comp_name, DATAMODEL_NAME_PREFIX, strlen(DATAMODEL_NAME_PREFIX)) == 0 || 
           strncmp(comp_name, ALGMODEL_NAME_PREFIX, strlen(ALGMODEL_NAME_PREFIX)) == 0;
}
