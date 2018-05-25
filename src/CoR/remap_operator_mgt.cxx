/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_operator_mgt.h"
#include "parse_special_words.h"
#include "remap_operator_bilinear.h"
#include "remap_operator_linear.h"
#include "remap_operator_distwgt.h"
#include "remap_operator_conserv_2D.h"
#include "remap_operator_smooth.h"
#include "remap_operator_spline_1D.h"
#include <string.h>
#include <stdio.h>


void Remap_operator_mgt::execute(const char*function, Remap_statement_operand **statement_operands, int num_operands)
{
    int i;
    int num_remap_grids;
    Remap_grid_class *remap_grids[256];
    char temp_string[256];
    int num_leaf_grids_src, num_leaf_grids_dst;
    Remap_grid_class *leaf_grids_src[256], *leaf_grids_dst[256];

    if (words_are_the_same(function, FUNCTION_WORD_NEW_OPERATOR)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3 || num_operands == 4, "function \"%s\" must have one result parameter and two or three input parameters\n", function);
        check_is_parameter_object_type_remap_operator(function, 0, statement_operands[0], "the new remap operator generated");
        check_is_parameter_string_type(function, 1, statement_operands[1], "the name of remap operator");
        for (i = 2, num_remap_grids = 0; i < num_operands; i ++) {
            check_is_parameter_object_type_grid(function, i, statement_operands[i], "the src or dst grid of remap operator");
            remap_grids[num_remap_grids++] = remap_grid_manager->search_remap_grid_with_grid_name(statement_operands[i]->object->object_name);
        }
        if (num_remap_grids == 2) {
            remap_grids[0]->get_leaf_grids(&num_leaf_grids_src, leaf_grids_src, remap_grids[0]);
            remap_grids[1]->get_leaf_grids(&num_leaf_grids_dst, leaf_grids_dst, remap_grids[1]);
            EXECUTION_REPORT(REPORT_ERROR, -1, num_leaf_grids_src == num_leaf_grids_dst,
                         "the src and dst grids of the remap operator must have the same number of dimensions\n");
            for (i = 0; i < num_leaf_grids_src; i ++)
                EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(leaf_grids_src[i]->get_coord_label(), leaf_grids_dst[i]->get_coord_label()),
                             "the coordinate labels of source grid \"%s\" and destination grid \"%s\" must be the same\n",
                             remap_grids[0]->get_grid_name(),
                             remap_grids[1]->get_grid_name());
        }
        if (words_are_the_same(statement_operands[1]->extension_names[0], REMAP_OPERATOR_NAME_BILINEAR)) {
            remap_operators.push_back(new Remap_operator_bilinear(statement_operands[0]->object->object_name, num_remap_grids, remap_grids));
        }
        else if (words_are_the_same(statement_operands[1]->extension_names[0], REMAP_OPERATOR_NAME_CONSERV_2D)) {
            remap_operators.push_back(new Remap_operator_conserv_2D(statement_operands[0]->object->object_name, num_remap_grids, remap_grids));
        }
        else if (words_are_the_same(statement_operands[1]->extension_names[0], REMAP_OPERATOR_NAME_DISTWGT)) {
            remap_operators.push_back(new Remap_operator_distwgt(statement_operands[0]->object->object_name, num_remap_grids, remap_grids));
        }        
        else if (words_are_the_same(statement_operands[1]->extension_names[0], REMAP_OPERATOR_NAME_LINEAR)) {
            remap_operators.push_back(new Remap_operator_linear(statement_operands[0]->object->object_name, num_remap_grids, remap_grids));
        }
        else if (words_are_the_same(statement_operands[1]->extension_names[0], REMAP_OPERATOR_NAME_SPLINE_1D)) {
            remap_operators.push_back(new Remap_operator_spline_1D(statement_operands[0]->object->object_name, num_remap_grids, remap_grids));
        }        
        else if (words_are_the_same(statement_operands[1]->extension_names[0], REMAP_OPERATOR_NAME_SMOOTH)) {
            remap_operators.push_back(new Remap_operator_smooth(statement_operands[0]->object->object_name, num_remap_grids, remap_grids));
        }
        else EXECUTION_REPORT(REPORT_ERROR, -1, false,
                          "the first input parameter \"%s\" of function \"%s\" is not a correct name of existing remap algorithms\n",
                          statement_operands[1]->extension_names[0],
                          function);        
        for (i = 0; i < num_remap_grids; i ++)
            remap_grids[i]->end_grid_definition_stage(remap_operators[remap_operators.size()-1]);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_SET_OPERATOR_PARA)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" must have three input parameters\n", function);
        check_is_parameter_object_type_remap_operator(function, 1, statement_operands[0], "the remap operator to be set parameter");
        check_is_parameter_string_type(function, 2, statement_operands[1], "the name of parameter");
        check_is_parameter_string_type(function, 3, statement_operands[2], "the value of parameter");   
        search_remap_operator(statement_operands[0]->object->object_name)->set_parameter(statement_operands[1]->extension_names[0], 
                                                                                         statement_operands[2]->extension_names[0]);
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "function \"%s\" is not support in the management of remap operators\n", function);
}


Remap_operator_basis *Remap_operator_mgt::search_remap_operator(const char *object_name)
{
    for (int i = 0; i < remap_operators.size(); i ++)
        if (remap_operators[i]->match_remap_operator(object_name))
            return remap_operators[i];

    return NULL;
}


Remap_operator_basis *Remap_operator_mgt::search_remap_operator(Remap_grid_class *grid_src, Remap_grid_class *grid_dst, const char *operator_name)
{
    for (int i = 0; i < remap_operators.size(); i ++)
        if (remap_operators[i]->match_remap_operator(grid_src, grid_dst, operator_name))
            return remap_operators[i];

    return NULL;
}


void Remap_operator_mgt::add_remap_operator(Remap_operator_basis *new_operator)
{
    remap_operators.push_back(new_operator);
}


Remap_operator_mgt::~Remap_operator_mgt()
{
    for (int i = 0; i < remap_operators.size(); i ++)
        delete remap_operators[i];
}


int Remap_operator_mgt::get_remap_operator_num_dim(const char *operator_name)
{
    if (words_are_the_same(operator_name, REMAP_OPERATOR_NAME_BILINEAR) || words_are_the_same(operator_name, REMAP_OPERATOR_NAME_DISTWGT) ||
        words_are_the_same(operator_name, REMAP_OPERATOR_NAME_CONSERV_2D))
        return 2;

    if (words_are_the_same(operator_name, REMAP_OPERATOR_NAME_LINEAR) || words_are_the_same(operator_name, REMAP_OPERATOR_NAME_SPLINE_1D))
        return 1;

    return -1;
}


int Remap_operator_mgt::check_operator_parameter(const char *operator_name, const char *parameter_name, const char *parameter_value, char *error_string)
{
    Remap_operator_basis *remap_operator;
    if (words_are_the_same(operator_name, REMAP_OPERATOR_NAME_BILINEAR))
        remap_operator = new Remap_operator_bilinear();
    else if (words_are_the_same(operator_name, REMAP_OPERATOR_NAME_CONSERV_2D))
        remap_operator = new Remap_operator_conserv_2D();
    else if (words_are_the_same(operator_name, REMAP_OPERATOR_NAME_DISTWGT))
        remap_operator = new Remap_operator_distwgt();
    else if (words_are_the_same(operator_name, REMAP_OPERATOR_NAME_LINEAR)) 
        remap_operator = new Remap_operator_linear();
    else if (words_are_the_same(operator_name, REMAP_OPERATOR_NAME_SPLINE_1D))
        remap_operator = new Remap_operator_spline_1D();
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Remap_operator_mgt::check_operator_parameter");

    int check_result = remap_operator->check_parameter(parameter_name, parameter_value, error_string);
    delete remap_operator;

    return check_result;
}


