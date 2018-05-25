/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "annotation_mgt.h"
#include "execution_report.h"
#include "remap_common_utils.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
  

Annotation_mgt::Annotation_mgt()
{
    annoation_lookup_table = new Dictionary<const char*>(1024);
}


Annotation_mgt::~Annotation_mgt()
{
    delete annoation_lookup_table;
    for (int i = 0; i < annotations.size(); i ++)
        delete annotations[i];
}


void Annotation_mgt::add_annotation(int object_id, const char *tag, const char *annotation)
{
    char key[NAME_STR_SIZE];
    char *local_annotation = strdup(annotation);


    sprintf(key, "%x @ %s", object_id, tag);
    annoation_lookup_table->insert(key, local_annotation);
    annotations.push_back(local_annotation);
}


const char *Annotation_mgt::get_annotation(int object_id, const char *tag)
{
    char key[NAME_STR_SIZE];
    const char *annotation;


    sprintf(key, "%x @ %s", object_id, tag);
    annotation = annoation_lookup_table->search(key, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, annotation != NULL, "Software error in Annotation_mgt::get_annotation");
    return annotation;
}

