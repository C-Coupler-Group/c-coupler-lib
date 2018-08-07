/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_weight_of_strategy_mgt.h"
#include "remap_weight_sparse_matrix.h"
#include "remap_common_utils.h"
#include "parse_special_words.h"


Remap_weight_of_operator_class *Remap_weight_of_operator_mgt::search_remap_weights_of_operator(Remap_grid_class *field_data_grid_src, Remap_grid_class *field_data_grid_dst, Remap_operator_basis *original_remap_operator)
{
    for (int i = 0; i < remap_weights_of_operators.size(); i ++) {
        if (field_data_grid_src->is_similar_grid_with(remap_weights_of_operators[i]->get_field_data_grid_src())&&
            field_data_grid_dst->is_similar_grid_with(remap_weights_of_operators[i]->get_field_data_grid_dst()) &&
            original_remap_operator == remap_weights_of_operators[i]->get_original_remap_operator())
            return remap_weights_of_operators[i];
    }
        
    return NULL;
}


void Remap_weight_of_operator_mgt::add_remap_weights_of_operator(Remap_weight_of_operator_class *remap_weight_of_operator)
{
/*
    EXECUTION_REPORT(REPORT_ERROR, -1, search_remap_weights_of_operator(remap_weight_of_operator->get_field_data_grid_src(), remap_weight_of_operator->get_field_data_grid_dst(), remap_weight_of_operator->get_original_remap_operator()) == NULL,
                     "C-Coupler error1 in add_remap_weights_of_operator of Remap_weight_of_operator_mgt");
*/
    remap_weights_of_operators.push_back(remap_weight_of_operator);
}


Remap_weight_of_operator_mgt::~Remap_weight_of_operator_mgt()
{
    for (int i = 0; i < remap_weights_of_operators.size(); i ++)
        delete remap_weights_of_operators[i];
}


void Remap_weight_of_strategy_mgt::execute(const char*function, Remap_statement_operand **statement_operands, int num_operands)
{
    int i;
    Remap_weight_of_strategy_class *remap_weights;
    Remap_grid_data_class *field_data_src, *field_data_dst;
    

    if (words_are_the_same(function, FUNCTION_WORD_COMPUTE_REMAP_WEIGHTS)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 4, "function \"%s\" must have one result parameter and three input parameters\n", function);
        check_is_parameter_object_type_remap_weights(function, 0, statement_operands[0], "the remap weights computed");
        check_is_parameter_object_type_remap_scheme(function, 1, statement_operands[1], "the remap scheme corresponding to the remap weights");
        check_is_parameter_object_type_grid(function, 2, statement_operands[2], "the src grid of remap weights (the grid of src field)\n");
        check_is_parameter_object_type_grid(function, 3, statement_operands[3], "the dst grid of remap weights (the grid of dst field)\n");
        generate_new_remap_weights(statement_operands[0]->object->object_name, statement_operands[1]->object->object_name,
                                   statement_operands[2]->object->object_name, statement_operands[3]->object->object_name,
                                   NULL, NULL, false);
        remap_weights_of_strategies[remap_weights_of_strategies.size()-1]->add_remap_weight_of_operators_to_manager(false);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_REMAP)){
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 3, "function \"%s\" must have no result parameter and three input parameters\n", function);
        check_is_parameter_object_type_remap_weights(function, 1, statement_operands[0], "the remap weights used for remapping");
        check_is_parameter_object_type_field_data(function, 2, statement_operands[1], "the src field of remapping");
        check_is_parameter_object_type_field_data(function, 3, statement_operands[2], "the dst field of remapping");     
        remap_weights = remap_weights_of_strategy_manager->search_remap_weight_of_strategy(statement_operands[0]->object->object_name);
        field_data_src = remap_field_data_manager->search_remap_field_data(statement_operands[1]->object->object_name);
        field_data_dst = remap_field_data_manager->search_remap_field_data(statement_operands[2]->object->object_name);
        remap_weights->do_remap(-1, field_data_src, field_data_dst);
    }
    else if (words_are_the_same(function, FUNCTION_WORD_READ_REMAP_WEIGHTS)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, num_operands == 6, "function \"%s\" must have one result parameter and five input parameters\n", function);
        check_is_parameter_object_type_remap_weights(function, 0, statement_operands[0], "the remap weights to be read from IO file");
        check_is_parameter_object_type_remap_scheme(function, 1, statement_operands[1], "the remap scheme corresponding to the remap weights");
           check_is_parameter_object_type_grid(function, 2, statement_operands[2], "the src grid of remap weights (the grid of src field)\n");
        check_is_parameter_object_type_grid(function, 3, statement_operands[3], "the dst grid of remap weights (the grid of dst field)\n");
        check_is_parameter_object_type_IO(function, 4, statement_operands[4], "the IO file to read remap weights");
        check_is_parameter_string_type(function, 5, statement_operands[5], "the format of remap weights in IO file (SCRIP or C-Coupler)");
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(statement_operands[5]->extension_names[0], "SCRIP") || words_are_the_same(statement_operands[5]->extension_names[0], "C-Coupler") ,
                     "the fifth input parameter of function \"%s\" must be a string of \"SCRIP\" or \"C-Coupler\"\n", function);
        remap_weights = generate_new_remap_weights(statement_operands[0]->object->object_name, statement_operands[1]->object->object_name,
                                                   statement_operands[2]->object->object_name, statement_operands[3]->object->object_name,
                                                   statement_operands[4]->object->object_name, statement_operands[5]->extension_names[0], true);
        remap_weights_of_strategies[remap_weights_of_strategies.size()-1]->add_remap_weight_of_operators_to_manager(false);
    }
}


Remap_weight_of_strategy_class *Remap_weight_of_strategy_mgt::search_or_add_remap_weight_of_strategy(Remap_grid_class *field_grid_src, Remap_grid_class *field_grid_dst, Remap_strategy_class *remap_strategy, const char *weight_object_name,
                                                                                                             const char *input_IO_file_name, const char *weight_IO_format, bool read_from_io)
{
    char temp_object_name[512];

    
    for (int i = 0; i < remap_weights_of_strategies.size(); i ++) 
        if (field_grid_src == remap_weights_of_strategies[i]->get_data_grid_src() && field_grid_dst == remap_weights_of_strategies[i]->get_data_grid_dst() && remap_strategy == remap_weights_of_strategies[i]->get_remap_strategy()) {
            if (weight_object_name != NULL)
                remap_weights_of_strategies[i]->renew_object_name(weight_object_name);
            return remap_weights_of_strategies[i];
        }

    if (weight_object_name != NULL)
        strcpy(temp_object_name, weight_object_name);
    else sprintf(temp_object_name, "TEMP_WEIGHT_%s_%s_%s", remap_strategy->get_strategy_name(), field_grid_src->get_grid_name(), field_grid_dst->get_grid_name());
    remap_weights_of_strategies.push_back(new Remap_weight_of_strategy_class(temp_object_name,remap_strategy->get_strategy_name(), field_grid_src->get_grid_name(), field_grid_dst->get_grid_name(),
                                                                             input_IO_file_name, weight_IO_format, read_from_io));

    return remap_weights_of_strategies[remap_weights_of_strategies.size()-1];
}


Remap_weight_of_strategy_class *Remap_weight_of_strategy_mgt::generate_new_remap_weights(const char *object_name, const char *remap_strategy_name, const char *data_grid_name_src, const char *data_grid_name_dst,
                                                                                         const char *input_IO_file_name, const char *weight_IO_format, bool read_from_io)
{
    Remap_strategy_class *remap_strategy;
    Remap_grid_class *data_grid_src, *data_grid_dst;
    
        
    remap_strategy = remap_strategy_manager->search_remap_strategy(remap_strategy_name);
    data_grid_src = remap_grid_manager->search_remap_grid_with_grid_name(data_grid_name_src);
    data_grid_dst = remap_grid_manager->search_remap_grid_with_grid_name(data_grid_name_dst);

    EXECUTION_REPORT(REPORT_ERROR, -1, remap_strategy != NULL && data_grid_src != NULL && data_grid_dst != NULL, "C-Coupler error in Remap_weight_of_strategy_class::Remap_weight_of_strategy_class");

    return search_or_add_remap_weight_of_strategy(data_grid_src, data_grid_dst, remap_strategy, object_name, input_IO_file_name, weight_IO_format, read_from_io);
}


Remap_weight_of_strategy_class *Remap_weight_of_strategy_mgt::search_remap_weight_of_strategy(const char *object_name)
{
    for (int i = 0; i < remap_weights_of_strategies.size(); i ++) 
        if (remap_weights_of_strategies[i]->match_object_name(object_name))
            return remap_weights_of_strategies[i];

    return NULL;
}


Remap_weight_of_strategy_mgt::~Remap_weight_of_strategy_mgt()
{
    for (int i = 0; i < remap_weights_of_strategies.size(); i ++)
        delete remap_weights_of_strategies[i];
}


void Remap_weight_of_strategy_mgt::add_remap_weight_of_strategy(Remap_weight_of_strategy_class *remap_weight_of_strategy)
{
    remap_weights_of_strategies.push_back(remap_weight_of_strategy);
}

