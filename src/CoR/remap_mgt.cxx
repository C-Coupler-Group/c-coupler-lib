/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_mgt.h"
#include "remap_parser.h"
#include "parse_special_words.h"
#include "execution_report.h"
#include <string.h>


void Remap_mgt::push_back_words(const char *type, const char *word)
{
    char *localized_word;

    
    localized_word = new char [256];
    strcpy(localized_word, word);
    if (words_are_the_same(type, "function"))
        function_words.push_back(localized_word);
    else if (words_are_the_same(type, "reserved"))
        reserved_words.push_back(localized_word);
}


void Remap_mgt::push_back_all_words()
{
    push_back_words("function", FUNCTION_WORD_ADD_NC_FILE); 
    push_back_words("function", FUNCTION_WORD_ADD_BIN_FILE);
    push_back_words("function", FUNCTION_WORD_NEW_1D_GRID);
    push_back_words("function", FUNCTION_WORD_NEW_PARTIAL_GRID);
    push_back_words("function", FUNCTION_WORD_ADD_GRID_AREA);
    push_back_words("function", FUNCTION_WORD_ADD_AREA_BOUND);
    push_back_words("function", FUNCTION_WORD_COMBINE_GRIDS);
    push_back_words("function", FUNCTION_WORD_NEW_OPERATOR);
    push_back_words("function", FUNCTION_WORD_COMBINE_OPERATORS);
    push_back_words("function", FUNCTION_WORD_WRITE_REMAP_WEIGHTS);
    push_back_words("function", FUNCTION_WORD_READ_REMAP_WEIGHTS);
    push_back_words("function", FUNCTION_WORD_REMAP);
    push_back_words("function", FUNCTION_WORD_WRITE_FIELD);
    push_back_words("function", FUNCTION_WORD_READ_FIELD);
    push_back_words("function", FUNCTION_WORD_SET_BOUNDARY);
    push_back_words("function", FUNCIION_WORD_LEV_COORD_FROM_SIGMA);
    push_back_words("function", FUNCIION_WORD_LEV_COORD_FROM_HYBRID);
    push_back_words("function", FUNCTION_WORD_EXTRACT_MASK);
    push_back_words("function", FUNCTION_WORD_COMPUTE_OCN_MASK);
    push_back_words("function", FUNCTION_WORD_ALLOC_FIELD);
    push_back_words("function", FUNCTION_WORD_READ_DATA);
    push_back_words("function", FUNCTION_WORD_ISPAN);
    push_back_words("function", FUNCTION_WORD_FSPAN);
    push_back_words("function", FUNCTION_WORD_GEN_TEST_DATA);
    push_back_words("function", FUNCTION_WORD_EVALUATE_ERROR);
    push_back_words("function", FUNCTION_WORD_COMPUTE_REMAP_WEIGHTS);
    push_back_words("function", FUNCTION_WORD_SET_OPERATOR_PARA);
    push_back_words("function", FUNCIION_WORD_SET_LEV_GRID_SIGMA_INFO);

    push_back_words("reserved", RESERVED_WORD_QUOTE_MARK);
    push_back_words("reserved", RESERVED_WORD_EQUAL);
    push_back_words("reserved", RESERVED_WORD_LEFT_BRACKET);
    push_back_words("reserved", RESERVED_WORD_RIGHT_BRACKET);
    push_back_words("reserved", RESERVED_WORD_COMMA);
    push_back_words("reserved", RESERVED_WORD_ATTRIBUTE);
}


Remap_mgt::Remap_mgt(const char *cfg_file_name)
{
    int num_words_in_statement;
    int i;
    char **words_in_statement;
    Remap_parser *remap_parser = NULL;


    /* Initialize parser for morphology checking and managers for managing data objects, 
        and push back funcion words and reserved words for semantic  analysis */
    io_manager = new IO_mgt();
    remap_strategy_manager = new Remap_strategy_mgt();
    remap_grid_manager = new Remap_grid_mgt();
    remap_operator_manager = new Remap_operator_mgt();
    remap_field_data_manager = new Remap_field_data_mgt();    
    remap_weights_of_strategy_manager = new Remap_weight_of_strategy_mgt();
    sequential_remap_weight_of_operator_manager = new Remap_weight_of_operator_mgt();
    parallel_remap_weight_of_operator_manager = new Remap_weight_of_operator_mgt();
    push_back_all_words();

    /* Initialize the data structure to keep each word in a statement */
    if (cfg_file_name != NULL) {
        remap_parser = new Remap_parser(cfg_file_name);
        words_in_statement = new char *[256];
        for (i = 0; i < 256; i ++)
            words_in_statement[i] = new char [256];

        /* For each statement, check its syntax, analyze its semantic, execute it and then release it */
        line_number = 1;
        while (remap_parser->get_next_parsed_statement(&num_words_in_statement, words_in_statement)) { 
            parse_statement(num_words_in_statement, words_in_statement);
            process_statement();
            release_statement();
            line_number ++;
        }
        delete remap_parser;
        for (i = 0; i < 256; i ++)
            delete [] words_in_statement[i];
        delete [] words_in_statement;
    }
}


/* Function parse_statement checks the syntax of the current statement, and then transforms
    the current statement into intermediate representation Remap_statement */
void Remap_mgt::parse_statement(const int num_words_in_statement, char **words_in_statement)
{
    int i;
    int matched_function_word_id;
    int matched_reserved_word_id;
    bool has_equal = false;
    bool in_function_stage = false;
    bool in_parameter_stage = false;
    bool enable_next_parameter = false;
    bool require_object_attribute = false;
    bool require_next_parameter = false;
    Remap_statement_object *statement_object;
    Remap_statement_object *last_statement_object = NULL;
    Remap_statement_operand *last_statement_operand = NULL;
    int matched_defined_statement_object_id;


    /* Allocate and initialize remap_statement */
    remap_statement = new Remap_statement;
    remap_statement->line_number = line_number;
    remap_statement->result_operand = NULL;
    strcpy(remap_statement->function, "\0");

    for (i = 0; i < num_words_in_statement; i ++) {
        matched_reserved_word_id = match_reserved_words(words_in_statement[i]);
        if (matched_reserved_word_id != -1) {
            EXECUTION_REPORT(REPORT_ERROR, -1, i > 0, 
                         "reserve word \"%s\" can not be the first word in the statement\n", 
                         words_in_statement[i]);
            
            /* Check the case of reserved word after reserved word */
            if (match_reserved_words(words_in_statement[i-1]) != -1) 
                EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(words_in_statement[i-1], RESERVED_WORD_LEFT_BRACKET) && words_are_the_same(words_in_statement[i], RESERVED_WORD_QUOTE_MARK) ||
                             words_are_the_same(words_in_statement[i-1], RESERVED_WORD_QUOTE_MARK) && words_are_the_same(words_in_statement[i], RESERVED_WORD_RIGHT_BRACKET)||
                             words_are_the_same(words_in_statement[i-1], RESERVED_WORD_QUOTE_MARK) && words_are_the_same(words_in_statement[i], RESERVED_WORD_COMMA) || 
                             words_are_the_same(words_in_statement[i-1], RESERVED_WORD_COMMA) && words_are_the_same(words_in_statement[i], RESERVED_WORD_QUOTE_MARK) ||
                             words_are_the_same(words_in_statement[i-1], RESERVED_WORD_EQUAL) && words_are_the_same(words_in_statement[i], RESERVED_WORD_QUOTE_MARK), 
                             "reserved word \"%s\" has a wrong reserved word \"%s\" before it\n", 
                             words_in_statement[i], 
                             words_in_statement[i-1]);
        }

        /* Enumerate the case of different kinds of words */
        if (words_are_the_same(words_in_statement[i], RESERVED_WORD_QUOTE_MARK)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, i+2 < num_words_in_statement && words_are_the_same(words_in_statement[i+2], RESERVED_WORD_QUOTE_MARK),
                         "reserve word \"%s\" must emerge in pair\n",
                         words_in_statement[i]);
            if (!words_are_the_same(words_in_statement[i-1], RESERVED_WORD_EQUAL))
                EXECUTION_REPORT(REPORT_ERROR, -1, in_parameter_stage && enable_next_parameter,
                             "the parameter \"%s\" specfied by reserve word \"%s\" must be a input parameter of function\n",
                             words_in_statement[i+1],
                             words_in_statement[i]);
            last_statement_operand = new Remap_statement_operand;
            last_statement_operand->object = NULL;
            last_statement_operand->num_extension_names = 1;
            strcpy(last_statement_operand->extension_names[0], words_in_statement[i+1]);
            remap_statement->src_operands.push_back(last_statement_operand);
            i=i+2;
            enable_next_parameter = false;
            require_next_parameter = false;
            continue;
        }

        else if (words_are_the_same(words_in_statement[i], RESERVED_WORD_EQUAL)) { 
            EXECUTION_REPORT(REPORT_ERROR, -1, !in_function_stage,
                         "function word can not be used before reserve word equal \"%s\"\n",
                         words_in_statement[i]);
            has_equal = true;
            last_statement_object = NULL;
        }
        
        else if (words_are_the_same(words_in_statement[i], RESERVED_WORD_LEFT_BRACKET)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, match_function_words(words_in_statement[i-1]) != -1,
                         "reserve word \"%s\" must follow function word while \"%s\" is not a function word\n",
                         words_in_statement[i],
                         words_in_statement[i-1]);            
            in_parameter_stage = true;
            enable_next_parameter = true;
        }

        else if (words_are_the_same(words_in_statement[i], RESERVED_WORD_RIGHT_BRACKET)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, !require_next_parameter && in_parameter_stage && !enable_next_parameter, 
                         "reserve word\"%s\" must be used to specify the parameters of function, the word before it must be reserve word \"%s\" or a parameter of function\n", 
                         RESERVED_WORD_RIGHT_BRACKET, 
                         RESERVED_WORD_LEFT_BRACKET);
            in_parameter_stage = false;
        }

        else if (words_are_the_same(words_in_statement[i], RESERVED_WORD_COMMA)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, in_parameter_stage && 
                         !require_next_parameter &&
                         !require_object_attribute &&
                         !enable_next_parameter, 
                         "reserve word \"%s\" must be used to separate the parameters of function, the word before it must be a parameter of function\n", 
                         RESERVED_WORD_COMMA);
            require_next_parameter = true; 
            enable_next_parameter = true;
        }

        else if (words_are_the_same(words_in_statement[i], RESERVED_WORD_ATTRIBUTE)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, last_statement_object != NULL, 
                         "reserve word\"%s\" is used to specify an extension variable of a defined object, the word before it must be an object name or an extension variable\n", 
                         words_in_statement[i]);
            require_object_attribute = true;
        }

        else if ((matched_function_word_id = match_function_words(words_in_statement[i])) != -1) {
            strcpy(remap_statement->function, words_in_statement[i]);
            in_function_stage = true;
            last_statement_object = NULL;
        }

        else {
            matched_defined_statement_object_id = match_defined_objects(words_in_statement[i]); 
            if (require_object_attribute) {
                strcpy(last_statement_operand->extension_names[last_statement_operand->num_extension_names], words_in_statement[i]);
                last_statement_operand->num_extension_names ++;
            }
            else {
                if (matched_defined_statement_object_id == -1) {
                    EXECUTION_REPORT(REPORT_ERROR, -1, !has_equal, 
                                 "\"%s\" is an undefined object, it can not be on the right side of word \"=\"\n", 
                                 words_in_statement[i]);
                    statement_object = new Remap_statement_object; 
                    strcpy(statement_object->object_name, words_in_statement[i]);
                    strcpy(statement_object->object_type, "\0");
                    statement_object->object_pointer = NULL;
                    last_statement_object = statement_object;
                    defined_objects.push_back(statement_object); 
                }
                else {
                    if (!has_equal && !in_function_stage)
                        EXECUTION_REPORT(REPORT_ERROR, -1, num_words_in_statement > i+1 && words_are_the_same(words_in_statement[i+1],RESERVED_WORD_ATTRIBUTE), 
                                     "\"%s\" is a defined object on the left side of word \"=\", it must have reserve word \"%s\" as well as extension variable following it\n", 
                                     words_in_statement[i], 
                                     RESERVED_WORD_ATTRIBUTE);
                    if (in_function_stage)
                        EXECUTION_REPORT(REPORT_ERROR, -1, enable_next_parameter, 
                                     "it may require \"%s\" before word \"%s\"\n", 
                                     RESERVED_WORD_COMMA, 
                                     words_in_statement[i]);
                    last_statement_object = defined_objects[matched_defined_statement_object_id];
                }
                last_statement_operand = new Remap_statement_operand;
                last_statement_operand->object = last_statement_object;
                last_statement_operand->num_extension_names = 0;
                if (!has_equal && !in_function_stage)
                    remap_statement->result_operand = last_statement_operand;
                else remap_statement->src_operands.push_back(last_statement_operand); 
            }
            require_next_parameter = false;
            enable_next_parameter = false;
            require_object_attribute = false;
        }
    }
}


/* Function process_statement checks the semantics of the current statement and then executes the current statement */
void Remap_mgt::process_statement()
{
    int i, j; 
    Remap_statement_operand *statement_operands[256];
    int num_operands = 0;


    if (strlen(remap_statement->function) == 0) {
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_statement->result_operand == NULL && remap_statement->src_operands.size() == 0, 
                     "this line does not have function word, it must be empty line\n");    
        return;
    }

    /* Determine the class type of defined object according to function */
    if (remap_statement->result_operand != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, strlen(remap_statement->result_operand->object->object_type) == 0 && remap_statement->result_operand->num_extension_names == 0 ||
                     strlen(remap_statement->result_operand->object->object_type) > 0 && remap_statement->result_operand->num_extension_names > 0, 
                     "remap software error: the object type of \"%s\" should not be set in this program step\n", 
                     remap_statement->result_operand->object->object_name);    
        if (remap_statement->result_operand->num_extension_names > 0) {
            if (words_are_the_same(remap_statement->function, FUNCTION_WORD_READ_FIELD) || 
                words_are_the_same(remap_statement->function, FUNCTION_WORD_SET_BOUNDARY) ||
                words_are_the_same(remap_statement->function, FUNCTION_WORD_COMPUTE_OCN_MASK) ||
                words_are_the_same(remap_statement->function, FUNCIION_WORD_LEV_COORD_FROM_SIGMA) ||
                words_are_the_same(remap_statement->function, FUNCIION_WORD_LEV_COORD_FROM_HYBRID) ||
                words_are_the_same(remap_statement->function, FUNCTION_WORD_EXTRACT_MASK) ||
                words_are_the_same(remap_statement->function, FUNCTION_WORD_FSPAN) ||
                words_are_the_same(remap_statement->function, FUNCTION_WORD_ISPAN))
                EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(remap_statement->result_operand->object->object_type, OBJECT_TYPE_GRID), 
                             "the object type of result parameter %s is %s (should be grid), which does not match function %s\n",
                             remap_statement->result_operand->object->object_name, remap_statement->result_operand->object->object_type, 
                             remap_statement->function);
            else EXECUTION_REPORT(REPORT_ERROR, -1, false, "function %s does not have result parameter with extension names\n", remap_statement->function);
        }
        else {
            if (words_are_the_same(remap_statement->function, FUNCTION_WORD_ADD_NC_FILE) ||
                words_are_the_same(remap_statement->function, FUNCTION_WORD_ADD_BIN_FILE)) {
                strcpy(remap_statement->result_operand->object->object_type, OBJECT_TYPE_IO);
            }
            else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_NEW_1D_GRID) ||
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_COMBINE_GRIDS) ||
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_NEW_PARTIAL_GRID)) {
                strcpy(remap_statement->result_operand->object->object_type, OBJECT_TYPE_GRID);
            }
            else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_NEW_OPERATOR)) {
                strcpy(remap_statement->result_operand->object->object_type, OBJECT_TYPE_REMAP_OPERATOR);
            }
            else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_COMBINE_OPERATORS)) {
                strcpy(remap_statement->result_operand->object->object_type, OBJECT_TYPE_REMAP_STRATEGY);
            }
            else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_COMPUTE_REMAP_WEIGHTS) ||
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_READ_REMAP_WEIGHTS)) {
                strcpy(remap_statement->result_operand->object->object_type, OBJECT_TYPE_REMAP_WEIGHTS);
            }
            else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_READ_FIELD) ||
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_ISPAN) ||
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_FSPAN) ||
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_ALLOC_FIELD) || 
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_READ_DATA) ||
                     words_are_the_same(remap_statement->function, FUNCTION_WORD_EVALUATE_ERROR)) {
                strcpy(remap_statement->result_operand->object->object_type, OBJECT_TYPE_FIELD_DATA);
            }
            else EXECUTION_REPORT(REPORT_ERROR, -1, false, "function %s does not have result parameter\n", remap_statement->function);
            remap_statement->result_operand->object->object_pointer = new char [256];
            sprintf((char*) remap_statement->result_operand->object->object_pointer, "%s(%s)", remap_statement->result_operand->object->object_name, remap_statement->result_operand->object->object_type);
        }
    }
    else {
        /*check the data type of the first input parameter for non-result function */
        if (words_are_the_same(remap_statement->function, FUNCTION_WORD_SET_OPERATOR_PARA))
            check_is_parameter_object_type_remap_operator(remap_statement->function, 1, remap_statement->src_operands[0], "the remap operator to be set parameter");
        else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_REMAP)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, remap_statement->src_operands[0]->object != NULL &&
                         (words_are_the_same(remap_statement->src_operands[0]->object->object_type, OBJECT_TYPE_REMAP_STRATEGY) ||
                          words_are_the_same(remap_statement->src_operands[0]->object->object_type, OBJECT_TYPE_REMAP_WEIGHTS)),
                          "the first parameter \"%s\" of function \"%s\" must be a defined remap scheme object or a defined remap weight object\n",
                          remap_statement->src_operands[0]->object->object_name,
                          remap_statement->function);
        }
        else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_WRITE_REMAP_WEIGHTS) ||
                 words_are_the_same(remap_statement->function, FUNCTION_WORD_WRITE_FIELD)) 
            check_is_parameter_object_type_IO(remap_statement->function, 1, remap_statement->src_operands[0], "the IO file to record the weight data");
        else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_ADD_GRID_AREA)) 
            check_is_parameter_object_type_grid(remap_statement->function, 1, remap_statement->src_operands[0], "the partial grid to be added the area");   
        else if (words_are_the_same(remap_statement->function, FUNCIION_WORD_SET_LEV_GRID_SIGMA_INFO)) 
            check_is_parameter_object_type_grid(remap_statement->function, 1, remap_statement->src_operands[0], "the level grid (vertical grid) to be set the sigma information");   
        else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_ADD_AREA_BOUND))
            check_is_parameter_object_type_grid(remap_statement->function, 1, remap_statement->src_operands[0], "the partial grid to be added the area bounds");  
        else if (words_are_the_same(remap_statement->function, FUNCTION_WORD_GEN_TEST_DATA) ||
                 words_are_the_same(remap_statement->function, FUNCTION_WORD_EVALUATE_ERROR))
            check_is_parameter_object_type_field_data(remap_statement->function, 1, remap_statement->src_operands[0], "which records the result of evaluation");
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in process_statement");
    }

    /* Pack the operands of statement */
    if (remap_statement->result_operand != NULL) 
        statement_operands[num_operands++] = remap_statement->result_operand;
    for (i = 0; i < remap_statement->src_operands.size(); i ++) 
        statement_operands[num_operands++] = remap_statement->src_operands[i];

    /* Call the execute function of the corresponding manager to check the semantics execute the statement */
    if (words_are_the_same(statement_operands[0]->object->object_type, OBJECT_TYPE_IO)) 
        io_manager->execute(remap_statement->function, statement_operands, num_operands);
    else if (words_are_the_same(statement_operands[0]->object->object_type, OBJECT_TYPE_GRID)) 
        remap_grid_manager->execute(remap_statement->function, statement_operands, num_operands);  
    else if (words_are_the_same(statement_operands[0]->object->object_type, OBJECT_TYPE_REMAP_OPERATOR)) 
        remap_operator_manager->execute(remap_statement->function, statement_operands, num_operands);
    else if (words_are_the_same(statement_operands[0]->object->object_type, OBJECT_TYPE_REMAP_STRATEGY)) 
        remap_strategy_manager->execute(remap_statement->function, statement_operands, num_operands);
    else if (words_are_the_same(statement_operands[0]->object->object_type, OBJECT_TYPE_REMAP_WEIGHTS))
        remap_weights_of_strategy_manager->execute(remap_statement->function, statement_operands, num_operands);
    else if (words_are_the_same(statement_operands[0]->object->object_type, OBJECT_TYPE_FIELD_DATA)) 
        remap_field_data_manager->execute(remap_statement->function, statement_operands, num_operands);
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error: \"%s\" is an unknown object type\n", remap_statement->src_operands[0]->object->object_type);
}


int Remap_mgt::match_function_words(const char *statement_word)
{
    for (int i = 0; i < function_words.size(); i ++) 
        if (words_are_the_same(function_words[i], statement_word)) 
            return i;

    return -1;
}

int Remap_mgt::match_reserved_words(const char *statement_word)
{
    for (int i = 0; i < reserved_words.size(); i ++) 
        if (words_are_the_same(reserved_words[i], statement_word)) 
            return i;

    return -1;
}


int Remap_mgt::match_defined_objects(const char *statement_word)
{
    for (int i = 0; i < defined_objects.size(); i ++) 
        if (words_are_the_same(defined_objects[i]->object_name, statement_word)) 
            return i;

    return -1;
}


void Remap_mgt::release_statement()
{
    if (remap_statement->result_operand != NULL) 
        delete remap_statement->result_operand;
    for (int j = 0; j < remap_statement->src_operands.size(); j ++)
        delete remap_statement->src_operands[j];
    delete remap_statement;
}


Remap_mgt::~Remap_mgt()
{
    for (int i = 0; i < function_words.size(); i ++)
        delete [] function_words[i];
    for (int i = 0; i < reserved_words.size(); i ++)
        delete [] reserved_words[i];
    for (int i = 0; i < defined_objects.size(); i ++) {
        delete [] defined_objects[i]->object_pointer;
        delete defined_objects[i];
    }

    delete io_manager;
    delete remap_weights_of_strategy_manager;
    delete sequential_remap_weight_of_operator_manager;
    delete parallel_remap_weight_of_operator_manager;
    delete remap_strategy_manager;
    delete remap_grid_manager;
    delete remap_operator_manager;
    delete remap_field_data_manager; 
}

