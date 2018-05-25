/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "field_info_mgt.h"
#include "common_utils.h"
#include "cor_global_data.h"
#include <string.h>


int Field_info_mgt::get_field_num_dims(const char *field_dim, const char *cfg_name)
{
    if (words_are_the_same(field_dim, FIELD_0_DIM))
        return 0;
    if (words_are_the_same(field_dim, FIELD_2_DIM))
        return 2;
    if (words_are_the_same(field_dim, FIELD_V1_DIM))
        return 1;
    if (words_are_the_same(field_dim, FIELD_3_DIM))
        return 3;
    if (words_are_the_same(field_dim, FIELD_4_DIM))
        return 4;

    EXECUTION_REPORT(REPORT_ERROR,-1, false, "\"%s\" is an undefined description of the number of dimensions of field. Please verify the configuration file %s", field_dim, cfg_name);
    return -1;
}


const field_attr *Field_info_mgt::search_field_info(const char *field_name)
{
    for (int i = 0; i < fields_attr.size(); i ++)
        if (words_are_the_same(field_name, fields_attr[i].field_name))
            return &fields_attr[i];

    return NULL;
}


void Field_info_mgt::add_field_info(const char *field_name, const char *field_long_name, const char *field_unit, const char *field_dim, const char *field_type, int line_number)
{
    field_attr local_attr;


    strcpy(local_attr.field_name, field_name);
    strcpy(local_attr.field_long_name, field_long_name);
    strcpy(local_attr.field_unit, field_unit);
    strcpy(local_attr.field_dim, field_dim);
    strcpy(local_attr.field_type, field_type);
    local_attr.line_number = line_number;
    fields_attr.push_back(local_attr);
    EXECUTION_REPORT(REPORT_ERROR, -1, search_field_info(local_attr.field_name) == &(fields_attr[fields_attr.size()-1]), "Software error in Field_info_mgt::add_field_info");
}


Field_info_mgt::Field_info_mgt()
{
    char XML_file_name[NAME_STR_SIZE];
    int line_number;
    field_attr local_attr;
    

    add_field_info("remap_frac", "fraction used for H2D remapping", "unitless", "H2D", "flux", -1);

    sprintf(XML_file_name, "%s/all/public_field_attribute.xml", comp_comm_group_mgt_mgr->get_config_root_dir());
    TiXmlDocument *XML_file = open_XML_file_to_read(-1, XML_file_name, MPI_COMM_WORLD, false);
    if (XML_file == NULL) {
        if (comp_comm_group_mgt_mgr->get_current_proc_global_id() == 0)
            EXECUTION_REPORT(REPORT_WARNING, -1, false, "There is no configuration file public_field_attribute.xml under the directory \"%s/all\", which indicates that no fields will be coupled among component models.", comp_comm_group_mgt_mgr->get_config_root_dir());
        return;
    }
    
    for (TiXmlNode *field_XML_node = XML_file->FirstChildElement(); field_XML_node != NULL; field_XML_node = field_XML_node->NextSibling()) {
        if (field_XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *field_XML_element = field_XML_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(field_XML_element->Value(),"field"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"field\". Please verify the XML file arround the line number %d.", XML_file_name, field_XML_element->Row());
        const char *field_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, field_XML_element, "name", XML_file_name, line_number, "name of a field", "configuration of the attributes of shared fields for coupling", true);
        const field_attr *existing_field = search_field_info(field_name);
        if (existing_field != NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "Cannot spefify the attributes of field \"%s\" in the XML file \"%s\" around the line number %d again because it has already been specified around the line number %d", field_name, XML_file_name, line_number, existing_field->line_number);
        const char *field_long_name = get_XML_attribute(-1, 1000, field_XML_element, "long_name", XML_file_name, line_number, "long name of a field", "configuration of the attributes of shared fields for coupling", true);
        const char *field_dimensions = get_XML_attribute(-1, -1, field_XML_element, "dimensions", XML_file_name, line_number, "information of dimensions (0D, H2D, V1D or V3D) of a field", "configuration of the attributes of shared fields for coupling", true);
        if (comp_comm_group_mgt_mgr->get_current_proc_global_id() == 0)
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(field_dimensions, FIELD_0_DIM) || words_are_the_same(field_dimensions, FIELD_2_DIM) || words_are_the_same(field_dimensions, FIELD_3_DIM) || words_are_the_same(field_dimensions, FIELD_V1_DIM), "The dimensions of field \"%s\" is wrong (must be \"0D\", \"H2D\", \"V1D\" or \"V3D\"). Please verify the XML file \"%s\" arround the line number %d.", field_name, XML_file_name, field_XML_element->Row());
        const char *default_unit = get_XML_attribute(-1, CCPL_NAME_STR_LEN, field_XML_element, "default_unit", XML_file_name, line_number, "default unit of a field", "configuration of the attributes of shared fields for coupling", true);
        const char *field_type = get_XML_attribute(-1, -1, field_XML_element, "type", XML_file_name, line_number, "default unit of a field", "configuration of the attributes of shared fields for coupling", true);
        add_field_info(field_name, field_long_name, default_unit, field_dimensions, field_type, line_number);
    }

    delete XML_file;
}


const char *Field_info_mgt::get_field_long_name(const char *field_name)
{
    if (search_field_info(field_name) == NULL)
        return NULL;
    
    return search_field_info(field_name)->field_long_name;
}


const char *Field_info_mgt::get_field_unit(const char *field_name)
{
    if (search_field_info(field_name) == NULL)
        return NULL;

    return search_field_info(field_name)->field_unit;
}

