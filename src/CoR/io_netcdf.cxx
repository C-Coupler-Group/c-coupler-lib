/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "io_netcdf.h"
#include "cor_global_data.h"
#include "execution_report.h"
#include "cor_cpl_interface.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


IO_netcdf::IO_netcdf(int ncfile_id)
{
    this->io_with_time_info = false;
    strcpy(this->object_name, "NULL");
    strcpy(this->file_type, FILE_TYPE_NETCDF);
    strcpy(this->file_name, "NULL");
    strcpy(this->open_format, "NULL");
    this->is_external_file = true;
    this->ncfile_id = ncfile_id;
}


IO_netcdf::IO_netcdf(const char *object_name, const char *file_name, const char *format, bool io_with_time_info)
{
    this->io_with_time_info = io_with_time_info;
    strcpy(this->object_name, object_name);
    strcpy(this->file_type, FILE_TYPE_NETCDF);
    strcpy(this->file_name, file_name);
    strcpy(this->open_format, format);
    this->is_external_file = false;
    if (words_are_the_same(format, "r"))
        rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
    else if (words_are_the_same(format, "w")) {
        rcode = nc_create(file_name, NC_CLOBBER, &ncfile_id);
        report_nc_error();
        nc_enddef(ncfile_id);        
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "the format of openning netcdf file must be read or write (\"r\" or \"w\")\n");
    report_nc_error();

    rcode = nc_close(ncfile_id);
    report_nc_error();
}


IO_netcdf::~IO_netcdf()
{
}


void IO_netcdf::report_nc_error()
{
    EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), file_name);
}


void IO_netcdf::datatype_from_netcdf_to_application(nc_type nc_data_type, char *application_data_type, const char *field_name)
{
    switch (nc_data_type) {
        case NC_BYTE:
            strcpy(application_data_type, DATA_TYPE_BOOL);
            break;
        case NC_CHAR:
            strcpy(application_data_type, DATA_TYPE_CHAR);
            break;
        case NC_INT:
            strcpy(application_data_type, DATA_TYPE_INT);
            break;
        case NC_SHORT:
            strcpy(application_data_type, DATA_TYPE_SHORT);
            break;
        case NC_FLOAT:
            strcpy(application_data_type, DATA_TYPE_FLOAT);
            break;
        case NC_DOUBLE:
            strcpy(application_data_type, DATA_TYPE_DOUBLE);
            break;
        default:
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "nc data type is not supported when reading field \"%s\" from netcdf file\n", field_name);
            break;
    }
}


void IO_netcdf::datatype_from_application_to_netcdf(const char *application_data_type, nc_type *nc_data_type)
{
    if (words_are_the_same(application_data_type, DATA_TYPE_BOOL))
        *nc_data_type = NC_INT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_CHAR))
        *nc_data_type = NC_CHAR;
    else if (words_are_the_same(application_data_type, DATA_TYPE_INT))
        *nc_data_type = NC_INT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_SHORT))
        *nc_data_type = NC_SHORT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_FLOAT))
        *nc_data_type = NC_FLOAT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_DOUBLE))
        *nc_data_type = NC_DOUBLE;
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in datatype_from_application_to_netcdf \"%s\"\n", application_data_type);
}


bool IO_netcdf::read_data(Remap_data_field *read_data_field, int time_pos, bool check_existence)
{
    int i, num_attributes, num_dimensions, variable_id, dimension_ids[256];
    char variable_name[256];
    nc_type nc_data_type;
    unsigned long data_size, dimension_size;
    Remap_field_attribute field_attribute;
    size_t starts[256], counts[256];


    rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
    report_nc_error();

    rcode = nc_inq_varid(ncfile_id, read_data_field->field_name_in_IO_file, &variable_id);
    if (!check_existence && rcode == NC_ENOTVAR) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Does not find the field \"%s\" in the data file \"%s\"", read_data_field->field_name_in_IO_file, file_name);
        rcode = nc_close(ncfile_id);
        report_nc_error();
        return false;
    }
    report_nc_error();
    rcode = nc_inq_var(ncfile_id, variable_id, variable_name, &nc_data_type, &num_dimensions, dimension_ids, &num_attributes);
    report_nc_error();

    for (i = 0, data_size = 1; i < num_dimensions; i ++) {
        rcode = nc_inq_dimlen(ncfile_id, dimension_ids[i], &dimension_size);
        report_nc_error();
        if (time_pos != -1 && i == 0) {
            EXECUTION_REPORT(REPORT_ERROR, -1, time_pos>=0 && time_pos < dimension_size, "C-Coupler error in IO_netcdf::read_data");
            starts[0] = time_pos;
            counts[0] = 1;
        }
        else {
            starts[i] = 0;
            counts[i] = dimension_size;
            data_size *= dimension_size;
        }
    }
    read_data_field->read_data_size = data_size;
    if (read_data_field->data_buf != NULL)
        EXECUTION_REPORT(REPORT_ERROR, -1, read_data_field->required_data_size == read_data_field->read_data_size, 
                     "the data size of field \"%s\" in netcdf file is different from of the data size of field \"%s\" determined by grid\n",
                     read_data_field->field_name_in_IO_file, 
                     read_data_field->field_name_in_application);
    else {
        read_data_field->required_data_size = data_size;
        read_data_field->data_buf = new char [data_size*get_data_type_size(read_data_field->data_type_in_application)];
    }

    /* Record the data type in netcdf file */
    datatype_from_netcdf_to_application(nc_data_type, read_data_field->data_type_in_IO_file, read_data_field->field_name_in_IO_file);

    /* Read and record the attributes of data */
    for (i = 0; i < num_attributes; i ++) {
        rcode = nc_inq_attname(ncfile_id, variable_id, i, field_attribute.attribute_name);
        report_nc_error();
        rcode = nc_inq_att(ncfile_id, variable_id, field_attribute.attribute_name, &nc_data_type, &field_attribute.attribute_size);
        report_nc_error();
        datatype_from_netcdf_to_application(nc_data_type, field_attribute.attribute_type, field_attribute.attribute_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, get_data_type_size(field_attribute.attribute_type)*field_attribute.attribute_size <= 8192, 
                     "value of attribute \"%s\" is out-of-bound (8192 bytes) when reading info of field \"%s\" from netcdf file\n", 
                     field_attribute.attribute_name, read_data_field->field_name_in_IO_file);
        switch (nc_data_type) {
            case NC_BYTE:
            case NC_CHAR:
                rcode = nc_get_att_text(ncfile_id, variable_id, field_attribute.attribute_name, field_attribute.attribute_value);
                field_attribute.attribute_value[field_attribute.attribute_size] = '\0';
                break;
            case NC_SHORT:
                rcode = nc_get_att_short(ncfile_id, variable_id, field_attribute.attribute_name, (short*)field_attribute.attribute_value);
                break;
            case NC_INT:
                rcode = nc_get_att_int(ncfile_id, variable_id, field_attribute.attribute_name, (int*)field_attribute.attribute_value);
                break;
            case NC_FLOAT:
                rcode = nc_get_att_float(ncfile_id, variable_id, field_attribute.attribute_name, (float*)field_attribute.attribute_value);
                break;
            case NC_DOUBLE:
                rcode = nc_get_att_double(ncfile_id, variable_id, field_attribute.attribute_name, (double*)field_attribute.attribute_value);
                break;
            default:
                EXECUTION_REPORT(REPORT_ERROR, -1, false, "Netcdf file may be corrupted: netcdf data type is unknown\n");
        }
        report_nc_error();
        read_data_field->field_attributes.push_back(field_attribute);
    }

    /* Read the data from netcdf file and transform data type when necessary */
    if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_CHAR)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_CHAR), 
                     "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) must be the same (char)\n",
                     read_data_field->field_name_in_application, read_data_field->field_name_in_IO_file);
        rcode = nc_get_vara_uchar(ncfile_id, variable_id, starts, counts, (unsigned char *) read_data_field->data_buf);
    }
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_BOOL)) {
        char *temp_buffer = new char[read_data_field->required_data_size*8];
        if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE)) {
            rcode = nc_get_vara_double(ncfile_id, variable_id, starts, counts, (double*) temp_buffer);            
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((double *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_FLOAT)) {
            rcode = nc_get_vara_float(ncfile_id, variable_id, starts, counts, (float*) temp_buffer);
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((float *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_INT)) {
            rcode = nc_get_vara_int(ncfile_id, variable_id, starts, counts, (int*) temp_buffer);
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((int *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_SHORT)) {
            rcode = nc_get_vara_short(ncfile_id, variable_id, starts, counts, (short*) temp_buffer);
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((short *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_CHAR))
            rcode = nc_get_vara_uchar(ncfile_id, variable_id, starts, counts, (unsigned char *) read_data_field->data_buf);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                          "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) in netcdf is \"%s\", which is not support\n",
                          read_data_field->field_name_in_application, 
                          read_data_field->field_name_in_IO_file,
                          read_data_field->data_type_in_IO_file);

        delete [] temp_buffer;
    }
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_FLOAT)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_FLOAT) || 
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_LONG) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_INT) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_SHORT), 
                     "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) in netcdf file must be float, double, long, int or short\n",
                     read_data_field->field_name_in_application, 
                     read_data_field->field_name_in_IO_file);
        rcode = nc_get_vara_float(ncfile_id, variable_id, starts, counts, (float*) read_data_field->data_buf);
    }
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_INT)) 
        rcode = nc_get_vara_int(ncfile_id, variable_id, starts, counts, (int *) read_data_field->data_buf);
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_SHORT)) 
        rcode = nc_get_vara_short(ncfile_id, variable_id, starts, counts, (short *) read_data_field->data_buf);
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_DOUBLE)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_FLOAT) || 
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_LONG) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_INT) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_SHORT), 
                     "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) in netcdf file must be float, double, long, int or short\n",
                     read_data_field->field_name_in_application, 
                     read_data_field->field_name_in_IO_file);
        rcode = nc_get_vara_double(ncfile_id, variable_id, starts, counts, (double*) read_data_field->data_buf);
    }

    report_nc_error();
    rcode = nc_close(ncfile_id);
    report_nc_error();

    return true;
}


void IO_netcdf::write_grid(Remap_grid_class *associated_grid, bool write_grid_name, bool use_script_format)
{
    int num_sized_sub_grids, num_leaf_grids, num_masked_sub_grids, num_sphere_leaf_grids, i, dim_ncid;
    Remap_grid_class *sized_sub_grids[256], *leaf_grids[256], *masked_sub_grids[256];
    Remap_grid_data_class *grid_data_field;
    char tmp_string[256];


    if (associated_grid == NULL)
        return;

    rcode = nc_open(file_name, NC_WRITE, &ncfile_id);
    report_nc_error();

    rcode = nc_redef(ncfile_id);
    report_nc_error();
    associated_grid->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
    for (i = 0; i < num_sized_sub_grids; i ++)
        if (sized_grids_map.find(sized_sub_grids[i]) == sized_grids_map.end()) {
            if (write_grid_name) {
                if (sized_sub_grids[i]->get_num_dimensions() == 1)
                    sprintf(tmp_string, "dim_%s_%d", sized_sub_grids[i]->get_coord_label(), recorded_grids.size());
                else if (sized_sub_grids[i]->get_is_sphere_grid())
                    sprintf(tmp_string, "dim_H2D_%d", recorded_grids.size());
                else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in IO_netcdf::write_grid");
            }
            else if (sized_sub_grids[i]->get_num_dimensions() == 1)
                sprintf(tmp_string, "%s", sized_sub_grids[i]->get_coord_label()); 
            else sprintf(tmp_string, "grid_size", sized_sub_grids[i]->get_grid_name()); 
            rcode = nc_def_dim(ncfile_id, tmp_string, sized_sub_grids[i]->get_grid_size(), &dim_ncid);
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "define dim %s for grid \"%s\" (%lx) in ncfile %s", tmp_string, sized_sub_grids[i]->get_grid_name(), sized_sub_grids[i], file_name);
            report_nc_error();
            sized_grids_map[sized_sub_grids[i]] = dim_ncid;
            recorded_grids.push_back(sized_sub_grids[i]);
        }
    if (use_script_format) {
        if (sized_grids_map.find(associated_grid) == sized_grids_map.end()) {
            rcode = nc_def_dim(ncfile_id, "grid_size", associated_grid->get_grid_size(), &dim_ncid);
            report_nc_error();
            sized_grids_map[associated_grid] = dim_ncid;
            recorded_grids.push_back(associated_grid);
        }
        rcode = nc_def_dim(ncfile_id, "grid_rank", num_sized_sub_grids, &dim_ncid);
        report_nc_error();
        int grid_dim_size[256], dims_ncid;
        for (int i = 0; i < num_sized_sub_grids; i ++)
            grid_dim_size[i] = sized_sub_grids[i]->get_grid_size();
        rcode = nc_def_var(ncfile_id, "grid_dims", NC_INT, 1, &dim_ncid, &dims_ncid);
        report_nc_error();
        rcode = nc_enddef(ncfile_id);
        report_nc_error();
        nc_put_var_int(ncfile_id, dims_ncid, grid_dim_size);
        report_nc_error();
        rcode = nc_redef(ncfile_id);
        report_nc_error();
    }
    nc_enddef(ncfile_id);
    report_nc_error();

    associated_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, associated_grid);
    for (i = 0; i < num_leaf_grids; i ++) {
        grid_data_field = leaf_grids[i]->get_grid_center_field();
        if (grid_data_field != NULL) 
            write_field_data(grid_data_field, associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format);
        if (leaf_grids[i]->get_sigma_grid_sigma_value_field() != NULL) {
            EXECUTION_REPORT(REPORT_ERROR, -1, grid_data_field == NULL, "Software error in IO_netcdf::write_grid");
            write_field_data(leaf_grids[i]->get_sigma_grid_sigma_value_field(), associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format);
            if (leaf_grids[i]->get_hybrid_grid_coefficient_field() != NULL)
                write_field_data(leaf_grids[i]->get_hybrid_grid_coefficient_field(), associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format);
            sprintf(tmp_string, "P0");
            if (write_grid_name)
                sprintf(tmp_string, "grid_%d_P0", get_recorded_grid_num(leaf_grids[i]));
            double P0 = leaf_grids[i]->get_sigma_grid_top_value();            
            rcode = nc_redef(ncfile_id);
            report_nc_error();
            rcode = nc_put_att_double(ncfile_id, NC_GLOBAL, tmp_string, NC_DOUBLE, 1, &P0);
            report_nc_error();
            nc_enddef(ncfile_id);
            report_nc_error();
        }
        grid_data_field = leaf_grids[i]->get_grid_vertex_field();
        if (grid_data_field != NULL && !grid_data_field->get_coord_value_grid()->get_are_vertex_values_set_in_default()) {
            if (grid_data_field->get_coord_value_grid()->get_num_dimensions() == 1)            
               sprintf(tmp_string, "num_vertexes_%s", grid_data_field->get_coord_value_grid()->get_coord_label());
            else if (use_script_format)
                sprintf(tmp_string, "grid_corners");
            else sprintf(tmp_string, "num_vertexes_H2D");
            rcode = nc_inq_dimid(ncfile_id, tmp_string, &dim_ncid);            
            if (rcode == NC_EBADDIM) {
                rcode = nc_redef(ncfile_id);
                report_nc_error();
                rcode = nc_def_dim(ncfile_id, tmp_string, grid_data_field->get_coord_value_grid()->get_num_vertexes(), &dim_ncid);
                nc_enddef(ncfile_id);
                report_nc_error();
            }
            write_field_data(grid_data_field, associated_grid, true, GRID_VERTEX_LABEL, dim_ncid, write_grid_name, use_script_format);            
        }
    }

    associated_grid->get_masked_sub_grids(&num_masked_sub_grids, masked_sub_grids);
    for (i = 0; i < num_masked_sub_grids; i ++)
        write_field_data(masked_sub_grids[i]->get_grid_mask_field(), associated_grid, true, GRID_MASK_LABEL, -1, write_grid_name, use_script_format);
    if (associated_grid->get_grid_imported_area() != NULL)
        write_field_data(associated_grid->get_grid_imported_area(), associated_grid, true, "area", -1, write_grid_name, use_script_format);

    rcode = nc_close(ncfile_id);
    report_nc_error();
}


void IO_netcdf::write_field_data(Remap_grid_data_class *field_data, 
                                Remap_grid_class *interchange_grid,
                                bool is_grid_data, 
                                const char *grid_field_type, 
                                int dim_ncid_num_vertex,
                                bool write_grid_name,
                                bool use_script_format)
{
    int num_sized_sub_grids, num_dims, i;
    unsigned long io_data_size, dimension_size;
    Remap_grid_class *sized_sub_grids[256];
    char tmp_string[256];
    int var_ncid, dim_ncids[256];
    size_t starts[256], counts[256];
    nc_type nc_data_type;


    if (!is_grid_data)
        field_data->set_masked_cell_to_missing_value();

    tmp_string[0] = '\0';
    if (is_grid_data && write_grid_name) {
        if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application, "mask"))
            sprintf(tmp_string, "grid_%d_%s", get_recorded_grid_num(field_data->get_coord_value_grid()), field_data->get_grid_data_field()->field_name_in_application);
        else {
            if (words_are_the_same(grid_field_type, GRID_VERTEX_LABEL))
                sprintf(tmp_string, "grid_%d_%s_%s", get_recorded_grid_num(field_data->get_coord_value_grid()), grid_field_type, field_data->get_grid_data_field()->field_name_in_application);
            else sprintf(tmp_string, "grid_%d_%s", get_recorded_grid_num(field_data->get_coord_value_grid()), field_data->get_grid_data_field()->field_name_in_application);
        }
    }
    else {
        if (!words_are_the_same(field_data->get_grid_data_field()->field_name_in_IO_file, "\0") && !is_grid_data)
            sprintf(tmp_string, "%s", field_data->get_grid_data_field()->field_name_in_IO_file);
        else {
            if (words_are_the_same(grid_field_type, GRID_VERTEX_LABEL))
                if (!use_script_format)
                    sprintf(tmp_string , "%s_%s", grid_field_type, field_data->get_grid_data_field()->field_name_in_application);
                else if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application,COORD_LABEL_LON))
                    sprintf(tmp_string , SCRIP_VERTEX_LON_LABEL);
                else sprintf(tmp_string , SCRIP_VERTEX_LAT_LABEL);
            else {
                sprintf(tmp_string , "%s", field_data->get_grid_data_field()->field_name_in_application);
                if (use_script_format) {
                    if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application, COORD_LABEL_LON))
                        sprintf(tmp_string , SCRIP_CENTER_LON_LABEL);
                    else if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application, COORD_LABEL_LAT))
                        sprintf(tmp_string , SCRIP_CENTER_LAT_LABEL);
                    else if (words_are_the_same(grid_field_type, GRID_MASK_LABEL))
                        sprintf(tmp_string , SCRIP_MASK_LABEL);
                }
            }
        }
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "IO field name is %s", tmp_string);
    }

    rcode = nc_inq_varid(ncfile_id, tmp_string, &var_ncid);
    if (rcode != NC_ENOTVAR) {
        if (is_grid_data)
            return;
        else EXECUTION_REPORT(REPORT_WARNING, -1, io_with_time_info,
                            "field data \"%s\" has been written to netcdf file \"%s\" before. The old data will be overwritten\n",
                            field_data->get_grid_data_field()->field_name_in_application, file_name);
    }

    if (interchange_grid != NULL) {
        field_data->interchange_grid_data(interchange_grid);        
        if (sized_grids_map.find(field_data->get_coord_value_grid()) != sized_grids_map.end()) {
            num_sized_sub_grids = 1;
            sized_sub_grids[0] = field_data->get_coord_value_grid();
        }
        else field_data->get_coord_value_grid()->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
        for (i = 0; i < num_sized_sub_grids; i ++) {
            EXECUTION_REPORT(REPORT_ERROR, -1, sized_grids_map.find(sized_sub_grids[i]) != sized_grids_map.end(), "remap software error1 in write_field_data\n");
            dim_ncids[num_sized_sub_grids-1-i] = sized_grids_map[sized_sub_grids[i]];
            starts[num_sized_sub_grids-1-i] = 0;
            counts[num_sized_sub_grids-1-i] = sized_sub_grids[i]->get_grid_size();
        }
        num_dims = num_sized_sub_grids;
    }
    else num_dims = 0;
    if (is_grid_data) {
        if (dim_ncid_num_vertex != -1) {
            starts[num_dims] = 0;
            counts[num_dims] = field_data->get_coord_value_grid()->get_num_vertexes();
            dim_ncids[num_dims++] = dim_ncid_num_vertex;
        }
    }
    for (i = 0, io_data_size = 1; i < num_dims; i ++) {
        rcode = nc_inq_dimlen(ncfile_id, dim_ncids[i], &dimension_size);
        report_nc_error();
        io_data_size *= dimension_size;
    }
    EXECUTION_REPORT(REPORT_ERROR, -1, field_data->get_grid_data_field()->required_data_size == io_data_size, "C-Coupler error: the data size in field for writing and IO file must be the same: %ld : %ld", field_data->get_grid_data_field()->required_data_size, io_data_size);
    if (!is_grid_data && io_with_time_info) {
        for (i = num_dims; i > 0; i --) {
            dim_ncids[i] = dim_ncids[i-1];
            starts[i] = starts[i-1];
            counts[i] = counts[i-1];
        }
        num_dims ++;
        dim_ncids[0] = time_dim_id;
        starts[0] = time_count - 1;
        counts[0] = 1;
    }

    rcode = nc_inq_varid(ncfile_id, tmp_string, &var_ncid);
    if (rcode == NC_ENOTVAR) {
        rcode = nc_redef(ncfile_id);
        report_nc_error();
        datatype_from_application_to_netcdf(field_data->get_grid_data_field()->data_type_in_IO_file, &nc_data_type);
        rcode = nc_def_var(ncfile_id, tmp_string, nc_data_type, num_dims, dim_ncids, &var_ncid);
        report_nc_error();
        for (i = 0; i < field_data->get_grid_data_field()->field_attributes.size(); i ++) {
            datatype_from_application_to_netcdf(field_data->get_grid_data_field()->field_attributes[i].attribute_type, &nc_data_type);
            switch (nc_data_type) {
                case NC_BYTE:
                case NC_CHAR:
                    rcode = nc_put_att_text(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name,
                                            field_data->get_grid_data_field()->field_attributes[i].attribute_size, field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_SHORT:
                    rcode = nc_put_att_short(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (short*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_INT:
                    rcode = nc_put_att_int(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (int*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_FLOAT:
                    rcode = nc_put_att_float(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (float*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_DOUBLE:
                    rcode = nc_put_att_double(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (double*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                default:
                    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in write_field_data\n");
            }
            report_nc_error();
        }
        nc_enddef(ncfile_id);
        report_nc_error();
    }

    if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_BOOL)) {
        int *temp_buffer = new int [field_data->get_grid_data_field()->required_data_size];
        for (long i = 0; i < field_data->get_grid_data_field()->required_data_size; i ++)
            if (((bool*)field_data->get_grid_data_field()->data_buf)[i])
                temp_buffer[i] = 1;
            else temp_buffer[i] = 0;
        rcode = nc_put_vara_int(ncfile_id, var_ncid, starts, counts, temp_buffer);
        delete [] temp_buffer;
    }
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_CHAR))
        rcode = nc_put_vara_schar(ncfile_id, var_ncid, starts, counts, (signed char *) field_data->get_grid_data_field()->data_buf);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT)) 
        rcode = nc_put_vara_float(ncfile_id, var_ncid, starts, counts, (float*) field_data->get_grid_data_field()->data_buf);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_INT))
        rcode = nc_put_vara_int(ncfile_id, var_ncid, starts, counts, (int *) field_data->get_grid_data_field()->data_buf);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_SHORT))
        rcode = nc_put_vara_short(ncfile_id, var_ncid, starts, counts, (short *) field_data->get_grid_data_field()->data_buf);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE))
        rcode = nc_put_vara_double(ncfile_id, var_ncid, starts, counts, (double*) field_data->get_grid_data_field()->data_buf);    
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error3 in write_field_data\n");
    report_nc_error(); 
}


void IO_netcdf::write_grided_data(Remap_grid_data_class *grided_data, bool write_grid_name, int date, int datesec, bool is_restart_field)
{
    unsigned long starts, counts, dim_len;
    int current_date, current_datesec;
    int time_var_id, date_var_id, datesec_var_id, tmp_var_id;;
    Remap_grid_data_class *tmp_field_data_for_io;


    if (execution_phase_number == 0)
        return;
    
    rcode = nc_open(file_name, NC_WRITE, &ncfile_id);
    report_nc_error();

    if (!io_with_time_info)
        EXECUTION_REPORT(REPORT_ERROR, -1, date == -1 && datesec == -1, "remap software error in write_grided_data \n");
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, date > 0 && datesec >= 0, "remap software error in write_grided_data \n");
        rcode = nc_inq_dimid(ncfile_id, "time", &time_dim_id); 
        if (rcode == NC_EBADDIM) {
            rcode = nc_redef(ncfile_id);
            report_nc_error();
            rcode = nc_def_dim(ncfile_id, "time", NC_UNLIMITED, &time_dim_id);
            report_nc_error();
            rcode = nc_def_var(ncfile_id, "time", NC_INT, 1, &time_dim_id, &time_var_id);
            report_nc_error();
            rcode = nc_def_var(ncfile_id, "date", NC_INT, 1, &time_dim_id, &date_var_id);
            report_nc_error();
            rcode = nc_def_var(ncfile_id, "datesec", NC_INT, 1, &time_dim_id, &datesec_var_id);
            report_nc_error();
            nc_enddef(ncfile_id);
            report_nc_error();
            time_count = 0;
            current_date = -1;
            current_datesec = -1;
        }
        else {
            rcode = nc_inq_dimlen(ncfile_id, time_dim_id, &dim_len);
            time_count = dim_len;
            report_nc_error();
            starts = time_count - 1;
            counts = 1;
            rcode = nc_inq_varid(ncfile_id, "date", &date_var_id);
            report_nc_error();
            rcode = nc_get_vara_int(ncfile_id, date_var_id, &starts, &counts, &current_date);
            report_nc_error();
            rcode = nc_inq_varid(ncfile_id, "datesec", &datesec_var_id);
            report_nc_error();
            rcode = nc_get_vara_int(ncfile_id, datesec_var_id, &starts, &counts, &current_datesec);
            report_nc_error();
        }
        if (current_date != date || current_datesec != datesec) {
            time_count ++;
            starts = time_count - 1;
            counts = 1;
            rcode = nc_inq_varid(ncfile_id, "time", &time_var_id);  
            report_nc_error();
            rcode = nc_put_vara_int(ncfile_id, time_var_id, &starts, &counts, &time_count);   
            report_nc_error();
            rcode = nc_inq_varid(ncfile_id, "date", &date_var_id);  
            report_nc_error();
            rcode = nc_put_vara_int(ncfile_id, date_var_id, &starts, &counts, &date);   
            report_nc_error();
            rcode = nc_inq_varid(ncfile_id, "datesec", &datesec_var_id);  
            report_nc_error();
            rcode = nc_put_vara_int(ncfile_id, datesec_var_id, &starts, &counts, &datesec);   
            report_nc_error();
        }
    }

    rcode = nc_close(ncfile_id);
    report_nc_error();

    write_grid(grided_data->get_coord_value_grid(), write_grid_name, false);

    rcode = nc_open(file_name, NC_WRITE, &ncfile_id);
    report_nc_error();

    if (strlen(grided_data->get_grid_data_field()->data_type_in_IO_file) == 0)
        strcpy(grided_data->get_grid_data_field()->data_type_in_IO_file, grided_data->get_grid_data_field()->data_type_in_application);
    tmp_field_data_for_io = generate_field_data_for_IO(grided_data, is_restart_field);
    write_field_data(tmp_field_data_for_io, grided_data->get_coord_value_grid(), false, "", -1, write_grid_name, false);
    if (tmp_field_data_for_io != grided_data)
        delete tmp_field_data_for_io;

    rcode = nc_close(ncfile_id);
    report_nc_error();
}


long IO_netcdf::get_dimension_size(const char *dim_name, MPI_Comm comm, bool is_root_proc)
{
    int dimension_id;
    long dimension_size = -1;


    if (is_root_proc) {
        rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
        report_nc_error();
        rcode = nc_inq_dimid(ncfile_id, dim_name, &dimension_id);
        if (rcode == NC_NOERR) {
            rcode = nc_inq_dimlen(ncfile_id, dimension_id, (unsigned long *)(&dimension_size));
            report_nc_error();   
        }
        rcode = nc_close(ncfile_id);
        report_nc_error();
    }

    if (comm != MPI_COMM_NULL)
        MPI_Bcast(&dimension_size, 1, MPI_LONG, 0, comm);

    return dimension_size;
}


void IO_netcdf::write_remap_weights(Remap_weight_of_strategy_class *remap_weights)
{
    Remap_grid_class *remap_grid_src, *remap_grid_dst, *leaf_grids[256];
    int dim_ncid_n_a, dim_ncid_n_b, dim_ncid_n_s, dim_ncid_nv_a, dim_ncid_nv_b;
    int col_id, row_id, S_id, area_a_id, area_b_id, yc_a_id, xc_a_id, yv_a_id, xv_a_id;
    int yc_b_id, xc_b_id, yv_b_id, xv_b_id, mask_a_id, mask_b_id;
    int num_leaf_grids, i;
    Remap_operator_basis *remap_operator;
    Remap_weight_sparse_matrix *weight_sparse_matrix;
    Remap_operator_grid *remap_operator_grid_src, *remap_operator_grid_dst;
    double *area_or_volumn_a, *area_or_volumn_b;
    int *temp_int_values, dim_ncids[2];
    long j;

    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(open_format, "w"), "can not write to netcdf file %s: %s, whose open format is not write\n", object_name, file_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights != NULL, "remap software error1 in write_remap_weights of netcdf file\n");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_grid_manager->search_remap_grid_with_grid_name(remap_weights->get_data_grid_src()->get_grid_name()) != NULL && remap_grid_manager->search_remap_grid_with_grid_name(remap_weights->get_data_grid_dst()->get_grid_name()) != NULL, "Software error in Remap_weight_of_operator_class::write_overall_remapping_weights");

    remap_grid_src = remap_grid_manager->search_remap_grid_with_grid_name(remap_weights->get_data_grid_src()->get_grid_name());
    remap_grid_dst = remap_grid_manager->search_remap_grid_with_grid_name(remap_weights->get_data_grid_dst()->get_grid_name());

	if (remap_weights->get_remap_strategy() != NULL) {
	    EXECUTION_REPORT(REPORT_ERROR, -1, remap_grid_src->are_all_vertex_fields_specified_by_user(),
	                 "all vertex values of coordinates in src grid \"%s\" must have been set by users\n",
	                 remap_grid_src->get_grid_name());
	    EXECUTION_REPORT(REPORT_ERROR, -1, remap_grid_dst->are_all_vertex_fields_specified_by_user(),
	                 "all vertex values of coordinates in dst grid \"%s\" must have been set by users\n",
	                 remap_grid_dst->get_grid_name());   
	}

    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights->get_data_grid_src()->get_is_sphere_grid() &&
                 (remap_weights->get_remap_strategy() == NULL || remap_weights->get_remap_strategy()->get_num_remap_operator() == 1 &&
                 remap_weights->get_remap_strategy()->get_remap_operator(0)->get_num_dimensions() == 2),
                 "for SCRIP format of remap weights, we only support horizontal 2D remap of only one remap algorithm\n");
    if (execution_phase_number == 1) {
        remap_operator = remap_weights->get_unique_remap_operator_of_weights();
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator != NULL,
                     "for SCRIP format of remap weights, we only support horizontal 2D remap of only one 2D remap algorithm\n");
        remap_operator_grid_src = new Remap_operator_grid(remap_weights->get_data_grid_src(), remap_operator, false, false);
        remap_operator_grid_dst = new Remap_operator_grid(remap_weights->get_data_grid_dst(), remap_operator, false, false);
        remap_operator_grid_src->update_operator_grid_data();
        remap_operator_grid_dst->update_operator_grid_data();
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator->get_num_remap_weights_groups() == 1,
                     "for SCRIP format of remap weights, we only support horizontal 2D remap of only one remap algorithm\n");
        weight_sparse_matrix = remap_operator->get_remap_weights_group(0);
        area_or_volumn_a = remap_weights->get_data_grid_src()->get_area_or_volumn();
        area_or_volumn_b = remap_weights->get_data_grid_dst()->get_area_or_volumn();
		if (remap_weights->get_remap_strategy() != NULL) 
	        EXECUTION_REPORT(REPORT_ERROR, -1, area_or_volumn_a != NULL && area_or_volumn_b != NULL, "remap software error2 in write_remap_weights of netcdf file\n");
        rcode = nc_open(file_name, NC_WRITE, &ncfile_id);
        rcode = nc_redef(ncfile_id);
        report_nc_error();
        rcode = nc_def_dim(ncfile_id, "n_a", remap_weights->get_data_grid_src()->get_grid_size(), &dim_ncid_n_a);
        report_nc_error();
        rcode = nc_def_dim(ncfile_id, "n_b", remap_weights->get_data_grid_dst()->get_grid_size(), &dim_ncid_n_b);
        report_nc_error();
        rcode = nc_def_dim(ncfile_id, "n_s", weight_sparse_matrix->get_num_weights(), &dim_ncid_n_s);
        report_nc_error();
        rcode = nc_def_var(ncfile_id, "col", NC_INT, 1, &dim_ncid_n_s, &col_id);
        report_nc_error();
        rcode = nc_def_var(ncfile_id, "row", NC_INT, 1, &dim_ncid_n_s, &row_id);
        report_nc_error();
        rcode = nc_def_var(ncfile_id, "S", NC_DOUBLE, 1, &dim_ncid_n_s, &S_id);
        report_nc_error();
		if (area_or_volumn_a != NULL) {
        	rcode = nc_def_var(ncfile_id, "area_a", NC_DOUBLE, 1, &dim_ncid_n_a, &area_a_id);
	        report_nc_error();            
		}
		if (area_or_volumn_b != NULL) {
        	rcode = nc_def_var(ncfile_id, "area_b", NC_DOUBLE, 1, &dim_ncid_n_b, &area_b_id);
	        report_nc_error();
		}
        rcode = nc_def_var(ncfile_id, "yc_a", NC_DOUBLE, 1, &dim_ncid_n_a, &yc_a_id);
        report_nc_error();
        rcode = nc_def_var(ncfile_id, "xc_a", NC_DOUBLE, 1, &dim_ncid_n_a, &xc_a_id);
        report_nc_error();
        rcode = nc_def_var(ncfile_id, "yc_b", NC_DOUBLE, 1, &dim_ncid_n_b, &yc_b_id);
        report_nc_error();
        rcode = nc_def_var(ncfile_id, "xc_b", NC_DOUBLE, 1, &dim_ncid_n_b, &xc_b_id);
        report_nc_error();

		if (remap_grid_src->are_all_vertex_fields_specified_by_user()) {
	        rcode = nc_def_dim(ncfile_id, "nv_a", remap_operator_grid_src->get_num_vertexes(), &dim_ncid_nv_a);
	        report_nc_error();            
	        dim_ncids[0] = dim_ncid_n_a;
	        dim_ncids[1] = dim_ncid_nv_a;
	        rcode = nc_def_var(ncfile_id, "yv_a", NC_DOUBLE, 2, dim_ncids, &yv_a_id);
	        report_nc_error();
	        rcode = nc_def_var(ncfile_id, "xv_a", NC_DOUBLE, 2, dim_ncids, &xv_a_id);
	        report_nc_error();
		}
		if (remap_grid_dst->are_all_vertex_fields_specified_by_user()) {
	        rcode = nc_def_dim(ncfile_id, "nv_b", remap_operator_grid_dst->get_num_vertexes(), &dim_ncid_nv_b);
	        report_nc_error();
	        dim_ncids[0] = dim_ncid_n_b;
	        dim_ncids[1] = dim_ncid_nv_b;
	        rcode = nc_def_var(ncfile_id, "yv_b", NC_DOUBLE, 2, dim_ncids, &yv_b_id);
	        report_nc_error();
	        rcode = nc_def_var(ncfile_id, "xv_b", NC_DOUBLE, 2, dim_ncids, &xv_b_id);
	        report_nc_error();
		}
        rcode = nc_def_var(ncfile_id, "mask_a", NC_INT, 1, &dim_ncid_n_a, &mask_a_id);
        report_nc_error();
        rcode = nc_def_var(ncfile_id, "mask_b", NC_INT, 1, &dim_ncid_n_b, &mask_b_id);
        report_nc_error();
        rcode = nc_enddef(ncfile_id);
        report_nc_error();
        temp_int_values = new int [weight_sparse_matrix->get_num_weights()];
        for (j = 0; j < weight_sparse_matrix->get_num_weights(); j ++)
            temp_int_values[j] = weight_sparse_matrix->get_indexes_src_grid()[j] + 1;
        rcode = nc_put_var_int(ncfile_id, col_id, temp_int_values);
        report_nc_error();
        for (j = 0; j < weight_sparse_matrix->get_num_weights(); j ++)
            temp_int_values[j] = weight_sparse_matrix->get_indexes_dst_grid()[j] + 1;
        rcode = nc_put_var_int(ncfile_id, row_id, temp_int_values);
        delete [] temp_int_values;
        report_nc_error();
        rcode = nc_put_var_double(ncfile_id, S_id, weight_sparse_matrix->get_weight_values());
        report_nc_error();
		if (area_or_volumn_a != NULL) {
        	rcode = nc_put_var_double(ncfile_id, area_a_id, area_or_volumn_a);
	        report_nc_error();
		}
		if (area_or_volumn_b != NULL) {
        	rcode = nc_put_var_double(ncfile_id, area_b_id, area_or_volumn_b);
	        report_nc_error();
		}
        
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_grid_src->get_center_coord_values()[0] != NULL &&
                         remap_operator_grid_src->get_center_coord_values()[1] != NULL, 
                         "remap software error3 in write_remap_weights of netcdf file\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_grid_dst->get_center_coord_values()[0] != NULL &&
                         remap_operator_grid_dst->get_center_coord_values()[1] != NULL, 
                         "remap software error4 in write_remap_weights of netcdf file\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_grid_src->get_vertex_coord_values()[0] != NULL &&
                         remap_operator_grid_src->get_vertex_coord_values()[1] != NULL, 
                         "remap software error5 in write_remap_weights of netcdf file\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_grid_dst->get_vertex_coord_values()[0] != NULL &&
                         remap_operator_grid_dst->get_vertex_coord_values()[1] != NULL, 
                         "remap software error6 in write_remap_weights of netcdf file\n");
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator_grid_src->get_mask_values() != NULL &&
                         remap_operator_grid_dst->get_mask_values() != NULL, 
                         "remap software error7 in write_remap_weights of netcdf file\n");
        rcode = nc_put_var_double(ncfile_id, yc_a_id, remap_operator_grid_src->get_center_coord_values()[1]);
        report_nc_error();
        rcode = nc_put_var_double(ncfile_id, xc_a_id, remap_operator_grid_src->get_center_coord_values()[0]);
        report_nc_error();
        rcode = nc_put_var_double(ncfile_id, yc_b_id, remap_operator_grid_dst->get_center_coord_values()[1]);
        report_nc_error();
        rcode = nc_put_var_double(ncfile_id, xc_b_id, remap_operator_grid_dst->get_center_coord_values()[0]);
        report_nc_error();
		if (remap_grid_src->are_all_vertex_fields_specified_by_user()) {
	        rcode = nc_put_var_double(ncfile_id, yv_a_id, remap_operator_grid_src->get_vertex_coord_values()[1]);
	        report_nc_error();
	        rcode = nc_put_var_double(ncfile_id, xv_a_id, remap_operator_grid_src->get_vertex_coord_values()[0]);
	        report_nc_error();
		}
		if (remap_grid_dst->are_all_vertex_fields_specified_by_user()) {
	        rcode = nc_put_var_double(ncfile_id, yv_b_id, remap_operator_grid_dst->get_vertex_coord_values()[1]);
    	    report_nc_error();
        	rcode = nc_put_var_double(ncfile_id, xv_b_id, remap_operator_grid_dst->get_vertex_coord_values()[0]);
	        report_nc_error();
    	}
        temp_int_values = new int [remap_operator->get_src_grid()->get_grid_size()];
        for (j = 0; j < remap_operator->get_src_grid()->get_grid_size(); j ++)
            if (remap_operator_grid_src->get_mask_values()[j])
                temp_int_values[j] = 1;
            else temp_int_values[j] = 0;
        rcode = nc_put_var_int(ncfile_id, mask_a_id, temp_int_values);
        report_nc_error();
        delete [] temp_int_values;
        temp_int_values = new int [remap_operator->get_dst_grid()->get_grid_size()];
        for (j = 0; j < remap_operator->get_dst_grid()->get_grid_size(); j ++)
            if (remap_operator_grid_dst->get_mask_values()[j])
                temp_int_values[j] = 1;
            else temp_int_values[j] = 0;
        rcode = nc_put_var_int(ncfile_id, mask_b_id, temp_int_values);
        report_nc_error();
        delete [] temp_int_values;
        rcode = nc_close(ncfile_id);
        report_nc_error();
        delete remap_operator_grid_src;
        delete remap_operator_grid_dst;
    }
}


void IO_netcdf::put_global_attr(const char *text_title, const void *attr_value, const char *local_data_type, const char *nc_data_type, int size)
{
    int nc_datatype;

    
    if (!is_external_file) {
        rcode = nc_open(file_name, NC_WRITE, &ncfile_id);
        report_nc_error();
    }
    rcode = nc_redef(ncfile_id);
    report_nc_error();

    if (words_are_the_same(nc_data_type, DATA_TYPE_STRING))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(local_data_type, DATA_TYPE_STRING), "software error in IO_netcdf::put_global_attr: miss match of data type");
    else if (words_are_the_same(nc_data_type, DATA_TYPE_FLOAT))
        nc_datatype = NC_FLOAT;
    else if (words_are_the_same(nc_data_type, DATA_TYPE_DOUBLE))
        nc_datatype = NC_DOUBLE;
    else if (words_are_the_same(nc_data_type, DATA_TYPE_INT))
        nc_datatype = NC_INT;
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "software error in IO_netcdf::put_global_attr: wrong nc data type %s", nc_datatype);
    
    if (words_are_the_same(local_data_type, DATA_TYPE_STRING))
        rcode = nc_put_att_text(ncfile_id, NC_GLOBAL, text_title, strlen((const char*)attr_value), (const char*)attr_value);
    else if (words_are_the_same(local_data_type, DATA_TYPE_FLOAT))
        rcode = nc_put_att_float(ncfile_id, NC_GLOBAL, text_title, NC_FLOAT, size, (const float*)attr_value);
    else if (words_are_the_same(local_data_type, DATA_TYPE_DOUBLE))
        rcode = nc_put_att_double(ncfile_id, NC_GLOBAL, text_title, NC_DOUBLE, size, (const double*)attr_value);
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "software error in IO_netcdf::put_global_attr: wrong local data type %s", local_data_type);
    report_nc_error();
    nc_enddef(ncfile_id);
    report_nc_error();
    if (!is_external_file) {
        rcode = nc_close(ncfile_id);
        report_nc_error();
    }
}


bool IO_netcdf::get_file_field_attribute(const char *field_name, const char *attribute_name, char *attribute_value, char *data_type)
{
    int variable_id = NC_GLOBAL;
    nc_type nc_data_type;
    unsigned long attribute_size;

    
    rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
    report_nc_error();

    if (field_name != NULL) {
        rcode = nc_inq_varid(ncfile_id, field_name, &variable_id);
        if (rcode != NC_NOERR) {
            rcode = nc_close(ncfile_id);
            report_nc_error();
            return false;
        }
    }

    rcode = nc_inq_att(ncfile_id, variable_id, attribute_name, &nc_data_type, &attribute_size);
    if (rcode != NC_NOERR) {
        rcode = nc_close(ncfile_id);
        report_nc_error();
        return false;    
    }

    switch (nc_data_type) {
        case NC_BYTE:
        case NC_CHAR:
            strcpy(data_type, DATA_TYPE_STRING);
            rcode = nc_get_att_text(ncfile_id, variable_id, attribute_name, attribute_value);
            attribute_value[attribute_size] = '\0';
            break;
        case NC_SHORT:
            rcode = nc_get_att_short(ncfile_id, variable_id, attribute_name, (short*)attribute_value);
            strcpy(data_type, DATA_TYPE_SHORT);
            break;
        case NC_INT:
            strcpy(data_type, DATA_TYPE_INT);
            rcode = nc_get_att_int(ncfile_id, variable_id, attribute_name, (int*)attribute_value);
            break;
        case NC_FLOAT:
            strcpy(data_type, DATA_TYPE_FLOAT);
            rcode = nc_get_att_float(ncfile_id, variable_id, attribute_name, (float*)attribute_value);
            break;
        case NC_DOUBLE:
            strcpy(data_type, DATA_TYPE_DOUBLE);
            rcode = nc_get_att_double(ncfile_id, variable_id, attribute_name, (double*)attribute_value);
            break;
        default:
            rcode = nc_close(ncfile_id);
            report_nc_error();
            return false;
    }    

    if (rcode != NC_NOERR) {
        rcode = nc_close(ncfile_id);
        report_nc_error();
        return false;
    }

    rcode = nc_close(ncfile_id);
    report_nc_error();
    
    return true;
}


bool IO_netcdf::get_file_field_string_attribute(const char *field_name, const char *attribute_name, char *attribute_value, char *data_type, MPI_Comm comm, bool is_root_proc)
{
    int success;

    attribute_value[0] = '\0';
    data_type[0] = '\0';

    if (is_root_proc)
        success = get_file_field_attribute(field_name, attribute_name, attribute_value, data_type)? 1 : 0;
    if (comm != MPI_COMM_NULL)
        MPI_Bcast(&success, 1, MPI_INT, 0, comm);
    if (success == 0)
        return false;
    if (comm != MPI_COMM_NULL) {
        MPI_Bcast(attribute_value, NAME_STR_SIZE, MPI_CHAR, 0, comm);
        MPI_Bcast(data_type, NAME_STR_SIZE, MPI_CHAR, 0, comm);
    }

    return true;
}


void IO_netcdf::read_file_field(const char *field_name, void **data_array_ptr, int *field_size, char *data_type, MPI_Comm comm, bool is_root_proc)
{
    int i, variable_id, *dim_ids, *dim_size, total_size, have_field = 1;
    size_t dim_len;
    nc_type nc_var_type;
    char *data_array = NULL;
    int num_dims = 0;

    *data_array_ptr = NULL;
    *field_size = -1;
    
    if (is_root_proc) {
        rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
        report_nc_error();
        rcode = nc_inq_varid(ncfile_id, field_name, &variable_id);
        if (rcode != NC_NOERR) {
            rcode = nc_close(ncfile_id);
            report_nc_error();
            have_field = 0;
        }
        else {
            rcode = nc_inq_varndims(ncfile_id, variable_id, &num_dims);
            report_nc_error();
            dim_ids = new int [num_dims];
            dim_size = new int [num_dims];
            rcode = nc_inq_vardimid(ncfile_id, variable_id, dim_ids);
            report_nc_error();
            for (i = 0; i < num_dims; i ++) {
                rcode = nc_inq_dimlen(ncfile_id, dim_ids[i], &dim_len);
                report_nc_error();
                dim_size[i] = dim_len;
            }
            rcode = nc_inq_vartype(ncfile_id, variable_id, &nc_var_type);
            report_nc_error();
            datatype_from_netcdf_to_application(nc_var_type, data_type, field_name);

            total_size = 1;
            for (i = 0; i < num_dims; i ++)
                total_size *= dim_size[i];
            data_array = new char [total_size*get_data_type_size(data_type)];
            switch(nc_var_type) {
                case NC_SHORT:
                    rcode = nc_get_var_short(ncfile_id, variable_id, (short*)data_array);
                    break;
                case NC_INT:
                    rcode = nc_get_var_int(ncfile_id, variable_id, (int*)data_array);
                    break;
                case NC_FLOAT:
                    rcode = nc_get_var_float(ncfile_id, variable_id, (float*)data_array);
                    break;
                case NC_DOUBLE:
                    rcode = nc_get_var_double(ncfile_id, variable_id, (double*)data_array);
                    break;
                default:
                    EXECUTION_REPORT(REPORT_ERROR, -1, false, "software error in IO_netcdf::read_file_field: data type %s is not supported", data_type);
                    break;
            }
            report_nc_error();
            *field_size = total_size;
            delete [] dim_ids;
            delete [] dim_size;
            rcode = nc_close(ncfile_id);
            report_nc_error();
        }
    }

    if (comm != MPI_COMM_NULL) {
        MPI_Bcast(&have_field, 1, MPI_INT, 0, comm);
        if (have_field == 0)
            return;
        MPI_Bcast(field_size, 1, MPI_INT, 0, comm);
        MPI_Bcast(data_type, NAME_STR_SIZE, MPI_CHAR, 0, comm);
		if (field_size > 0) {
	        if (data_array == NULL)
	            data_array = new char [(*field_size)*get_data_type_size(data_type)];
	        MPI_Bcast(data_array, (*field_size)*get_data_type_size(data_type), MPI_CHAR, 0, comm);
		}
    }

    *data_array_ptr = data_array;
}


void IO_netcdf::read_remap_weights(Remap_weight_of_strategy_class *remap_weights, Remap_strategy_class *remap_strategy, bool read_weight_values)
{
    double *area, *weight_values;
    int var_id;
    unsigned long grid_size, num_weights, i;
    long *indexes_src, *indexes_dst;
    int *tmp_indexes_src, *tmp_indexes_dst;
    Remap_operator_basis *duplicated_remap_operator;
    Remap_weight_of_operator_instance_class *remap_operator_instance;
    Remap_weight_sparse_matrix *weight_sparse_matrix;


    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(open_format, "r"), "can not read netcdf file %s: %s, whose open format is not read\n", object_name, file_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights->get_data_grid_src()->are_all_vertex_fields_specified_by_user(),
                 "all vertex values of coordinates in src grid \"%s\" must have been set by users\n",
                 remap_weights->get_data_grid_src()->get_grid_name());
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights->get_data_grid_dst()->are_all_vertex_fields_specified_by_user(),
                 "all vertex values of coordinates in dst grid \"%s\" must have been set by users\n",
                 remap_weights->get_data_grid_dst()->get_grid_name());
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights->get_data_grid_src()->get_is_sphere_grid() &&
                 remap_weights->get_remap_strategy()->get_num_remap_operator() == 1,
                 "for SCRIP format of remap weights, we only support horizontal 2D remap of only one remap algorithm\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights->get_data_grid_src()->is_similar_grid_with(remap_weights->get_remap_strategy()->get_remap_operator(0)->get_src_grid()),
                 "the src grid of remap operator does not match the src grid of remap weights\n");
    EXECUTION_REPORT(REPORT_ERROR, -1, remap_weights->get_data_grid_dst()->is_similar_grid_with(remap_weights->get_remap_strategy()->get_remap_operator(0)->get_dst_grid()),
                 "the src grid of remap operator does not match the src grid of remap weights\n");

    if (execution_phase_number == 1) {
        rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
        report_nc_error();
        rcode = nc_inq_dimid(ncfile_id, "n_a", &var_id);
        report_nc_error();
        rcode = nc_inq_dimlen(ncfile_id, var_id, &grid_size);
        report_nc_error();
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_size == remap_weights->get_data_grid_src()->get_grid_size() && remap_weights->get_data_grid_src()->get_area_or_volumn() != NULL, 
                     "the src grid in netcdf file does not match the src grid of remap operator\n");
        area = new double [grid_size];
        rcode = nc_inq_varid(ncfile_id, "area_a", &var_id);
        report_nc_error();
        rcode = nc_get_var_double(ncfile_id, var_id, remap_weights->get_data_grid_src()->get_area_or_volumn());
        report_nc_error();
        rcode = nc_inq_dimid(ncfile_id, "n_b", &var_id);
        report_nc_error();
        rcode = nc_inq_dimlen(ncfile_id, var_id, &grid_size);
        report_nc_error();
        EXECUTION_REPORT(REPORT_ERROR, -1, grid_size == remap_weights->get_data_grid_dst()->get_grid_size() && remap_weights->get_data_grid_dst()->get_area_or_volumn() != NULL, 
                     "the dst grid in netcdf file does not match the dst grid of remap operator\n");
        area = new double [grid_size];
        rcode = nc_inq_varid(ncfile_id, "area_b", &var_id);
        if (rcode == NC_NOERR) {
            rcode = nc_get_var_double(ncfile_id, var_id, remap_weights->get_data_grid_dst()->get_area_or_volumn());
            report_nc_error();
        }
        rcode = nc_inq_dimid(ncfile_id, "n_s", &var_id);
        report_nc_error();
        rcode = nc_inq_dimlen(ncfile_id, var_id, &num_weights);
        report_nc_error();
        if (read_weight_values) {
            indexes_src = new long [num_weights];
            indexes_dst = new long [num_weights];
            tmp_indexes_src = new int [num_weights];
            tmp_indexes_dst = new int [num_weights];
            weight_values = new double [num_weights];
            rcode = nc_inq_varid(ncfile_id, "col", &var_id);
            report_nc_error();
            rcode = nc_get_var_int(ncfile_id, var_id, tmp_indexes_src);
            report_nc_error();
            rcode = nc_inq_varid(ncfile_id, "row", &var_id);
            report_nc_error();
            rcode = nc_get_var_int(ncfile_id, var_id, tmp_indexes_dst);
            report_nc_error();
            rcode = nc_inq_varid(ncfile_id, "S", &var_id);
            report_nc_error();
            rcode = nc_get_var_double(ncfile_id, var_id, weight_values);
            report_nc_error();
            for (i = 0; i < num_weights; i ++) {
                indexes_src[i] = tmp_indexes_src[i] - 1;
                indexes_dst[i] = tmp_indexes_dst[i] - 1;
            }
            duplicated_remap_operator = remap_strategy->get_remap_operator(0)->duplicate_remap_operator(false);
            weight_sparse_matrix = new Remap_weight_sparse_matrix(remap_strategy->get_remap_operator(0),
                                                                  num_weights, indexes_src, indexes_dst, weight_values, 
                                                                  0, NULL);
            duplicated_remap_operator->add_weight_sparse_matrix(weight_sparse_matrix);
        }
        else duplicated_remap_operator = NULL;
        remap_operator_instance = new Remap_weight_of_operator_instance_class(remap_weights->get_data_grid_src(),
                                                                     remap_weights->get_data_grid_dst(),
                                                                     0, remap_strategy->get_remap_operator(0),
                                                                     duplicated_remap_operator);
        remap_weights->add_remap_weight_of_operator_instance(remap_operator_instance, remap_weights->get_data_grid_src(), remap_weights->get_data_grid_dst(), remap_strategy->get_remap_operator(0), remap_strategy->get_remap_operator(0)->get_src_grid(), remap_strategy->get_remap_operator(0)->get_dst_grid());

        if (read_weight_values) {
            delete [] tmp_indexes_src;
            delete [] tmp_indexes_dst;
        }
        delete [] area;
    }
}

