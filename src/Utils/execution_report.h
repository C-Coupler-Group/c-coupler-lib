/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef EXECUTION_REPORT_H
#define EXECUTION_REPORT_H

#define EXECUTION_REPORT                        execution_report
#define EXECUTION_REPORT_LOG                    if (report_internal_log_enabled) execution_report
#define EXECUTION_REPORT_ERROR_OPTIONALLY       if (report_error_enabled) execution_report



#include <stdarg.h>


#define REPORT_ERROR               ((int)1)
#define REPORT_LOG                 ((int)2)
#define REPORT_WARNING             ((int)3)
#define REPORT_PROGRESS            ((int)4)
#define REPORT_CONSTANTLY          ((int)5)
#define REPORT_EXTERNAL_LOG        ((int)6)


#define CCPL_NULL_INT                                   ((int)(0x7FFFFFFF))
#define CCPL_NULL_LONG                                  ((long)(0x7FFFFFFFFFFFFFFF))


extern bool report_external_log_enabled;
extern bool report_error_enabled;
extern bool report_progress_enabled;
extern bool report_internal_log_enabled;
extern bool flush_log_file;


extern void import_report_setting();
extern void wtime(double *t);


extern void execution_report(int, int, bool);
extern void execution_report(int, int, bool, const char *, ...);


#endif
