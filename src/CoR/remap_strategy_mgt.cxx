/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_strategy_mgt.h"
#include "parse_special_words.h"
#include <string.h>
#include <stdio.h>


void Remap_strategy_mgt::execute(const char*function, Remap_statement_operand **statement_operands, int num_operands)
{
    int i;
    int num_remap_operators = 0;
    Remap_operator_basis *remap_operators[1024];


    if (words_are_the_same(function, FUNCTION_WORD_COMBINE_OPERATORS)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands >= 2, "function \"%s\" must have one result parameter and at least two input parameters\n", function);
        check_is_parameter_object_type_remap_scheme(function, 0, statement_operands[0], "the remap scheme to be generated");
        for (i = 1; i < num_operands; i ++) {
            check_is_parameter_object_type_remap_operator(function, i, statement_operands[i], "one remap operator to generate the remap scheme");
            remap_operators[num_remap_operators++] = remap_operator_manager->search_remap_operator(statement_operands[i]->object->object_name);
        }
        remap_strategies.push_back(new Remap_strategy_class(statement_operands[0]->object->object_name, num_remap_operators, remap_operators));
    }
    else if (words_are_the_same(function, FUNCTION_WORD_REMAP)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" must have non result parameter and three input parameters\n", function);
        check_is_parameter_object_type_remap_scheme(function, 1, statement_operands[0], "the remap scheme used to do remapping");
        check_is_parameter_object_type_field_data(function, 2, statement_operands[1], "the src field of remapping");
        check_is_parameter_object_type_field_data(function, 3, statement_operands[2], "the dst field of remapping");        
        remap_strategy_manager->search_remap_strategy(statement_operands[0]->object->object_name)->remap_fields(statement_operands[1]->object->object_name, 
                                                                                                                statement_operands[2]->object->object_name);
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "\"%s\" is an undefined function\n", function);
}


Remap_strategy_class *Remap_strategy_mgt::search_remap_strategy(const char *strategy_name)
{
    for (int i = 0; i < remap_strategies.size(); i ++)
        if (remap_strategies[i]->match_remap_strategy(strategy_name))
            return remap_strategies[i];

    return NULL;
}


Remap_strategy_mgt::~Remap_strategy_mgt()
{
    for (int i = 0; i < remap_strategies.size(); i ++)
        delete remap_strategies[i];
}


void Remap_strategy_mgt::add_remap_strategy(Remap_strategy_class *remap_strategy)
{
    remap_strategies.push_back(remap_strategy);
}

