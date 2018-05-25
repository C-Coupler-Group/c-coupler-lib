/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "execution_report.h"
#include "cor_global_data.h"
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>


bool report_external_log_enabled;
bool report_error_enabled;
bool report_progress_enabled;
bool report_internal_log_enabled;
bool flush_log_file;


void import_report_setting()
{
    char XML_file_name[NAME_STR_SIZE];
    int line_number;
    char keywords[5][NAME_STR_SIZE];
    bool report_setting[5];


    report_external_log_enabled = false;
    report_error_enabled = false;
    report_internal_log_enabled = false;
    report_progress_enabled = false;
    flush_log_file = false;

    sprintf(XML_file_name, "%s/all/CCPL_report.xml", comp_comm_group_mgt_mgr->get_config_root_dir());
    TiXmlDocument *XML_file = open_XML_file_to_read(-1, XML_file_name, MPI_COMM_WORLD, false);
    if (XML_file == NULL)
        return;

    sprintf(keywords[0], "report_internal_log");
    sprintf(keywords[1], "report_external_log");
    sprintf(keywords[2], "report_progress");
    sprintf(keywords[3], "report_error");
    sprintf(keywords[4], "flush_log_file");
    
    TiXmlElement *XML_element = XML_file->FirstChildElement();
    for (int i = 0; i < 5; i ++) {
        report_setting[i] = false;
        const char *setting = XML_element->Attribute(keywords[i], &line_number);
        if (setting == NULL)
            continue;
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(setting, "on") || words_are_the_same(setting, "off"), "Error happens when loading the XML file \"%s\": the value of \"%s\" must be \"on\" or \"off\". Please verify the XML file around line number %d", XML_file_name, keywords[i]);
        report_setting[i] = words_are_the_same(setting,"on");
    }

    delete XML_file;

    report_internal_log_enabled = report_setting[0];
    report_external_log_enabled = report_setting[1];
    report_progress_enabled = report_setting[2];
    report_error_enabled = report_setting[3];
    flush_log_file = report_setting[4];
}


void wtime(double *t)
{
    static int sec = -1;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    if (sec < 0) 
        sec = tv.tv_sec;
    *t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}


void report_header(int report_type, int comp_id, bool &condition, char *output_format)
{
    if (comp_id != -1 && (comp_id == comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id() || !comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id, true) || components_time_mgrs->get_time_mgr(comp_id) == NULL))
        comp_id = -1;
    
    output_format[0] = '\0';
    
    switch (report_type) {
        case REPORT_ERROR:
            condition = !condition;
            break;
        case REPORT_LOG:
            condition = report_internal_log_enabled;
            break;
        case REPORT_EXTERNAL_LOG:
            condition = report_external_log_enabled;
            break;
        case REPORT_WARNING:
            condition = !condition;
            break;
        case REPORT_PROGRESS:
            if (comp_id == -1)
                condition = comp_comm_group_mgt_mgr->get_current_proc_global_id() == 0 && report_progress_enabled;
            else condition = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,true,"")->get_current_proc_local_id() == 0 && report_progress_enabled;
            break;
        case REPORT_CONSTANTLY:
            condition = true;
            break;
        default:
            printf("report type %d is not support\n", report_type);
            assert(false);
    }

    if (!condition)
        return;

    switch (report_type) {
        case REPORT_ERROR:
            sprintf(output_format+strlen(output_format), "C-Coupler REPORT ERROR: ");
            break;
        case REPORT_LOG:
            sprintf(output_format+strlen(output_format), "C-Coupler REPORT LOG: ");    
            break;
        case REPORT_EXTERNAL_LOG:
            sprintf(output_format+strlen(output_format), "C-Coupler REPORT LOG: ");
            break;
        case REPORT_WARNING:
            sprintf(output_format+strlen(output_format), "C-Coupler REPORT WARNING: ");
            break;
        case REPORT_PROGRESS:
             sprintf(output_format+strlen(output_format), "C-Coupler REPORT PROGRESS: ");
            break;
        case REPORT_CONSTANTLY:
            sprintf(output_format+strlen(output_format), "C-Coupler REPORT PROGRESS: ");
            break;
        default:
            printf("Software error: report type %d is not supported\n", report_type);
            assert(false);
            break;
    }
    if (comp_comm_group_mgt_mgr != NULL) {
        if (comp_id == -1)
            sprintf(output_format+strlen(output_format)-2, " in the root component model corresponding to the executable named \"%s\": ", comp_comm_group_mgt_mgr->get_executable_name());
        else {
            int current_date = components_time_mgrs->get_time_mgr(comp_id)->get_current_date();
            int current_second = components_time_mgrs->get_time_mgr(comp_id)->get_current_second();
            int current_step_id = components_time_mgrs->get_time_mgr(comp_id)->get_current_step_id();
            sprintf(output_format+strlen(output_format)-2, " in the component model \"%s\" corresponding to the executable named \"%s\", at the current simulation time of %08d-%05d (the current step number is %d): ", 
                    comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,true,"execution report")->get_full_name(), comp_comm_group_mgt_mgr->get_executable_name(), current_date, current_second, current_step_id);
        }
    }
}


void execution_report(int report_type, int comp_id, bool condition, const char *format, ...)
{
    char output_format[NAME_STR_SIZE*16];
    char output_string[NAME_STR_SIZE*16*16];


    report_header(report_type, comp_id, condition, output_format);
    
    if (!condition)
        return;
    
    strcat(output_format, format);
    strcat(output_format, "\n\n");

    if (comp_id != -1 && (comp_id == comp_comm_group_mgt_mgr->get_global_node_root()->get_comp_id() || !comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id, true) || components_time_mgrs->get_time_mgr(comp_id) == NULL))
        comp_id = -1;

    va_list pArgList;
    va_start(pArgList, format);
    vsprintf(output_string, output_format, pArgList);
    va_end(pArgList);    
    if (comp_comm_group_mgt_mgr == NULL) {
        fprintf(stdout, output_string);
        if (flush_log_file || report_type == REPORT_ERROR)
            fflush(stdout);
    }
    else {
        const char *log_file_name1 = comp_comm_group_mgt_mgr->get_exe_log_file_name();
        const char *log_file_name2 = NULL;
        comp_comm_group_mgt_mgr->output_log(output_string, flush_log_file || report_type == REPORT_ERROR);
        if (comp_id != -1) {
            log_file_name2 = comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comp_ccpl_log_file_name();
            comp_comm_group_mgt_mgr->search_global_node(comp_id)->output_log(output_string, flush_log_file || report_type == REPORT_ERROR);
        }
        if (report_type == REPORT_ERROR) {
            if (log_file_name2 == NULL)
                printf("ERROR happens. Please check the log file \"%s\"\n\n", log_file_name1);
            else printf("ERROR happens. Please check the log file \"%s\" or \"%s\"\n\n", log_file_name1, log_file_name2);
        }
    }

    if (report_type == REPORT_ERROR)
        assert(false);
}



void execution_report(int report_type, int comp_id, bool condition) 
{
    execution_report(report_type, comp_id, condition, "report without explicit hint");
}


