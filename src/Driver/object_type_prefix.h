/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef OBJECT_TYPE_PREFIX_H
#define OBJECT_TYPE_PREFIX_H


#define TYPE_ID_PREFIX_MASK                    ((int)0xFF000000)
#define TYPE_ID_SUFFIX_MASK                    ((int)0x00FFFFFF)
#define TYPE_COMP_LOCAL_ID_PREFIX              ((int)0x01000000)
#define TYPE_GRID_LOCAL_ID_PREFIX              ((int)0x03000000)
#define TYPE_GRID_GLOBAL_ID_PREFIX             ((int)0x04000000)
#define TYPE_DECOMP_ID_PREFIX                  ((int)0x06000000)
#define TYPE_FIELD_INST_ID_PREFIX              ((int)0x07000000)
#define TYPE_TIMER_ID_PREFIX                   ((int)0x08000000)
#define TYPE_INOUT_INTERFACE_ID_PREFIX         ((int)0x09000000)
#define TYPE_IO_FIELD_PREFIX                   ((int)0x0a000000)


#endif

