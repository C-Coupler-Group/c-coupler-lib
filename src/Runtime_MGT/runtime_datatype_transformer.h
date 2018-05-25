/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#ifndef RUNTIME_DATATYPE_TRANSFORMER
#define RUNTIME_DATATYPE_TRANSFORMER


#include <vector>
#include "memory_mgt.h"


class Runtime_datatype_transformer
{
    private:
        std::vector<Field_mem_info*> src_fields;
        std::vector<Field_mem_info*> dst_fields;

    public:
        Runtime_datatype_transformer() {}
        Runtime_datatype_transformer(Field_mem_info*, Field_mem_info*);
        ~Runtime_datatype_transformer() {}
        void transform_fields_datatype();
        bool run(bool);
};

#endif

