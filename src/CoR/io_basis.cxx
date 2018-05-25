/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "io_basis.h"
#include <math.h>
#include <string.h>


template <class T1, class T2> void copy_data_values_for_IO(T1 *data_values_in_application, T1 fill_value_in_application, T2 *data_values_in_io, T2 fill_value_in_io, bool *mask, long size, bool is_restart_field)
{
    for (long i = 0; i < size; i ++) {
        data_values_in_io[i] = (T2) data_values_in_application[i];
        if (!is_restart_field && ((mask != NULL &&!mask[i]) || 
            fabs((double)data_values_in_application[i]) >= 1.0e10 ||
            (fabs((double)data_values_in_application[i]) <= 1.5*fabs((double)fill_value_in_application) && fabs((double)data_values_in_application[i]) >= 0.5*fabs((double)fill_value_in_application))))
            data_values_in_io[i] = fill_value_in_io;
    }
}


template <class T> void transfer_data_from_float_to_short(T *data_values_in_application, T fill_value_in_application, short *data_values_in_io, short fill_value_in_io, bool *mask, long size, double add_offset, double scale_factor)
{
    for (long i = 0; i < size; i ++) {
        if ((mask != NULL &&!mask[i]) || 
            fabs((double)data_values_in_application[i]) >= 1.0e10 ||
            (fabs((double)data_values_in_application[i]) <= 1.5*fabs((double)fill_value_in_application) && fabs((double)data_values_in_application[i]) >= 0.5*fabs((double)fill_value_in_application)))
            data_values_in_io[i] = fill_value_in_io;
        else data_values_in_io[i] = (data_values_in_application[i]-add_offset)/(scale_factor);
    }
}


bool IO_basis::match_IO_object(const char *object_name)
{
    return words_are_the_same(object_name, this->object_name);
}


Remap_grid_data_class *IO_basis::generate_field_data_for_IO(Remap_grid_data_class *field_data_in_application, bool is_restart_field)
{
    Remap_grid_data_class *field_data_for_io;


    if (field_data_in_application->get_coord_value_grid() == NULL) {
        strcpy(field_data_in_application->get_grid_data_field()->data_type_in_IO_file, field_data_in_application->get_grid_data_field()->data_type_in_application);
        return field_data_in_application;
    }

    field_data_in_application->interchange_grid_data(field_data_in_application->get_coord_value_grid());
    field_data_for_io = field_data_in_application->duplicate_grid_data_field(field_data_in_application->get_coord_value_grid(), field_data_in_application->get_grid_data_field()->required_data_size/field_data_in_application->get_coord_value_grid()->get_grid_size(), false, false);
    strcpy(field_data_for_io->get_grid_data_field()->data_type_in_IO_file, field_data_in_application->get_grid_data_field()->data_type_in_IO_file);
    if (is_restart_field)
        strcpy(field_data_for_io->get_grid_data_field()->data_type_in_IO_file, field_data_in_application->get_grid_data_field()->data_type_in_application);
        
    field_data_for_io->change_datatype_in_application(field_data_for_io->get_grid_data_field()->data_type_in_IO_file);
    copy_field_data_for_IO(field_data_in_application, field_data_for_io, is_restart_field);

    return field_data_for_io;
}


void IO_basis::copy_field_data_for_IO(Remap_grid_data_class *field_data_in_application, Remap_grid_data_class *field_data_in_io, bool is_restart_field)
{
    Remap_grid_class *mask_sub_grids[256];
    int num_mask_sub_grids;
    Remap_grid_data_class *runtime_mask_field = NULL;
    bool *runtime_mask_values = NULL;
    double fill_value_application = 0, fill_value_io = 0;
    double scale_factor, add_offset;


    field_data_in_application->get_grid_data_field()->read_fill_value();
    if (field_data_in_application->get_grid_data_field()->have_fill_value)
        fill_value_application = field_data_in_application->get_grid_data_field()->fill_value;
    field_data_in_io->get_grid_data_field()->read_fill_value();
    if (field_data_in_io->get_grid_data_field()->have_fill_value)
        fill_value_io = field_data_in_io->get_grid_data_field()->fill_value;    

    if (field_data_in_application->get_coord_value_grid() != NULL) {
        field_data_in_application->get_coord_value_grid()->compute_remap_field_data_runtime_mask(field_data_in_application->get_coord_value_grid(), mask_sub_grids, &num_mask_sub_grids, &runtime_mask_field);
        if (runtime_mask_field != NULL)
            runtime_mask_values = (bool*) runtime_mask_field->get_grid_data_field()->data_buf;
    }

    if (words_are_the_same(field_data_in_application->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE)) 
        if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE)) 
            copy_data_values_for_IO((double*)field_data_in_application->get_grid_data_field()->data_buf, (double)fill_value_application, 
                                   (double*)field_data_in_io->get_grid_data_field()->data_buf, (double)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT))
            copy_data_values_for_IO((double*)field_data_in_application->get_grid_data_field()->data_buf, (double)fill_value_application, 
                                   (float*)field_data_in_io->get_grid_data_field()->data_buf, (float)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_SHORT)) {
            field_data_in_io->get_grid_data_field()->read_scale_factor_and_add_offset(&scale_factor, &add_offset);
            transfer_data_from_float_to_short((double*)field_data_in_application->get_grid_data_field()->data_buf, (double)fill_value_application, 
                                   (short*)field_data_in_io->get_grid_data_field()->data_buf, (short)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, add_offset, scale_factor);
        }
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error1 in copy_field_data_for_IO\n");
    else if (words_are_the_same(field_data_in_application->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT)) 
        if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE))
            copy_data_values_for_IO((float*)field_data_in_application->get_grid_data_field()->data_buf, (float)fill_value_application, 
                                   (double*)field_data_in_io->get_grid_data_field()->data_buf, (double)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT))
            copy_data_values_for_IO((float*)field_data_in_application->get_grid_data_field()->data_buf, (float)fill_value_application, 
                                   (float*)field_data_in_io->get_grid_data_field()->data_buf, (float)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_SHORT)) {
            field_data_in_io->get_grid_data_field()->read_scale_factor_and_add_offset(&scale_factor, &add_offset);
            transfer_data_from_float_to_short((float*)field_data_in_application->get_grid_data_field()->data_buf, (float)fill_value_application, 
                                   (short*)field_data_in_io->get_grid_data_field()->data_buf, (short)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, add_offset, scale_factor);
        }
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error2 in copy_field_data_for_IO\n");
    else if (words_are_the_same(field_data_in_application->get_grid_data_field()->data_type_in_application, DATA_TYPE_LONG)) 
        if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_INT)) 
            copy_data_values_for_IO((long*)field_data_in_application->get_grid_data_field()->data_buf, (long)fill_value_application, 
                                   (int*)field_data_in_io->get_grid_data_field()->data_buf, (int)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error3 in copy_field_data_for_IO\n");    
    else if (words_are_the_same(field_data_in_application->get_grid_data_field()->data_type_in_application, DATA_TYPE_INT)) 
        if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_INT)) 
            copy_data_values_for_IO((int*)field_data_in_application->get_grid_data_field()->data_buf, (int)fill_value_application, 
                                   (int*)field_data_in_io->get_grid_data_field()->data_buf, (int)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error4 in copy_field_data_for_IO\n"); 
    else if (words_are_the_same(field_data_in_application->get_grid_data_field()->data_type_in_application, DATA_TYPE_BOOL)) 
        if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_INT)) 
            copy_data_values_for_IO((bool*)field_data_in_application->get_grid_data_field()->data_buf, (bool)fill_value_application, 
                                   (int*)field_data_in_io->get_grid_data_field()->data_buf, (int)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_BOOL)) 
                    copy_data_values_for_IO((bool*)field_data_in_application->get_grid_data_field()->data_buf, (bool)fill_value_application, 
                                           (bool*)field_data_in_io->get_grid_data_field()->data_buf, (bool)fill_value_io, runtime_mask_values, 
                                           field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error5 in copy_field_data_for_IO\n");
    else if (words_are_the_same(field_data_in_application->get_grid_data_field()->data_type_in_application, DATA_TYPE_SHORT)) 
        if (words_are_the_same(field_data_in_io->get_grid_data_field()->data_type_in_application, DATA_TYPE_SHORT)) 
            copy_data_values_for_IO((short*)field_data_in_application->get_grid_data_field()->data_buf, (short)fill_value_application, 
                                   (short*)field_data_in_io->get_grid_data_field()->data_buf, (short)fill_value_io, runtime_mask_values, 
                                   field_data_in_application->get_grid_data_field()->required_data_size, is_restart_field);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error6 in copy_field_data_for_IO\n");
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "C-Coupler error in copy_field_data_for_IO\n");

    if (runtime_mask_field != NULL)
        delete runtime_mask_field;
}


int IO_basis::get_recorded_grid_num(Remap_grid_class *grid)
{
    for (int i = 0; i < recorded_grids.size(); i ++)
        if (recorded_grids[i] == grid || words_are_the_same(recorded_grids[i]->get_grid_name(), grid->get_grid_name()))
            return i;

    EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in IO_basis::get_recorded_grid_num %s %lx", grid->get_grid_name(), grid);
    return -1;
}

