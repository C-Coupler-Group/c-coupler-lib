/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef FIELD_MGT
#define FIELD_MGT

#include "common_utils.h"
#include <vector>


#define FIELD_0_DIM       "0D"
#define FIELD_2_DIM       "H2D"
#define FIELD_V1_DIM      "V1D"
#define FIELD_3_DIM       "V3D"
#define FIELD_4_DIM       "4D"


struct field_attr
{
    char field_name[NAME_STR_SIZE];
    char field_long_name[NAME_STR_SIZE];
    char field_unit[NAME_STR_SIZE];
    char field_dim[NAME_STR_SIZE];           // dimension info: scalar, 1D, 2D, 3D, etc
    char field_type[NAME_STR_SIZE];          // state or flux
    int line_number;
};


class Field_info_mgt
{
private:
    std::vector<field_attr> fields_attr;
    
public:
    Field_info_mgt();
    ~Field_info_mgt() {}
    const field_attr* search_field_info(const char*);
    int get_field_num_dims(const char*, const char*);
    const char *get_field_long_name(const char*);
    const char *get_field_unit(const char*);
    void add_field_info(const char*, const char*, const char*, const char*, const char *, int);
};

#endif
