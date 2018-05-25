/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "stdlib.h"
#include "ensemble_mgt.h"


Ensemble_mgt::Ensemble_mgt()
{
    ensemble_member_id = -1;
    have_random_seed_for_perturbation = -1;
    root_random_seed_for_perturbation = -1;
    ensemble_random_seed_for_perturbation = -1;
}


void Ensemble_mgt::Initialize(int ensemble_member_id, int have_random_seed_for_perturbation, int root_random_seed_for_perturbation, const char *perturbation_type)
{
    EXECUTION_REPORT(REPORT_ERROR,-1, ensemble_member_id > 0, "The number of an ensemble member must be a positive integer\n");
    EXECUTION_REPORT_LOG(REPORT_LOG,-1, true, "Ensemble is enabled. The ensemble member id is %d", ensemble_member_id);

    this->ensemble_member_id = ensemble_member_id;
    this->have_random_seed_for_perturbation = (have_random_seed_for_perturbation == 1);
    this->root_random_seed_for_perturbation = root_random_seed_for_perturbation;
    if (have_random_seed_for_perturbation) {
        if (words_are_the_same("set_last_bit_to_1", perturbation_type))
            this->perturbation_type_id = 1;
        else if (words_are_the_same("set_last_bit_to_0", perturbation_type))
            this->perturbation_type_id = 2;
        else if (words_are_the_same("reverse_last_bit", perturbation_type))
            this->perturbation_type_id = 3;
        else if (words_are_the_same("xor_last_bit_with_a_bit", perturbation_type))
            this->perturbation_type_id = 4;
        else EXECUTION_REPORT(REPORT_ERROR,-1, false, "\"%s\" is not a right selection for perturbing the roundoff errors of the registered fields. Existing selections include set_last_bit_to_1, set_last_bit_to_0, reverse_last_bit and xor_last_bit_with_a_bit", perturbation_type);
        srand(root_random_seed_for_perturbation);
        for (int i = 0; i < ensemble_member_id; i ++)
            this->ensemble_random_seed_for_perturbation = rand();
        srand(this->ensemble_random_seed_for_perturbation);
        EXECUTION_REPORT_LOG(REPORT_LOG,-1, true, "In the ensemble experiment of perturbing roundoff errors, for the ensemble member %d, the seed of random number generation is %d", ensemble_member_id, this->ensemble_random_seed_for_perturbation);
    }
}


void Ensemble_mgt::register_a_field_for_perturbation(void *data_buf)
{
    Field_mem_info *registered_field;


    registered_field = memory_manager->search_field_via_data_buf(data_buf, false);
    EXECUTION_REPORT(REPORT_ERROR,-1, registered_field != NULL && registered_field->get_is_registered_model_buf(), 
                     "The field data buffer for perturbing the roundoff errors has not been registered to C-Coupler before. Please check.");
    EXECUTION_REPORT(REPORT_ERROR,-1, words_are_the_same(registered_field->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT) || words_are_the_same(registered_field->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE),
                     "The data type of field %s is not real4 or real8. It cannot be used for perturbing the roundoff errors. Please check.", registered_field->get_field_name());
    registered_fields_for_perturbation.push_back(registered_field);
}


void Ensemble_mgt::perturb_a_field_through_set_last_bit_to_1(void *field_data_buf, const char *data_type, long field_size, int current_random_number)
{
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
        for (long i = 0; i < field_size; i ++) 
            ((int*) field_data_buf)[i] = (((int*) field_data_buf)[i] | ((int)1));
    else if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
        for (long i = 0; i < field_size; i ++) 
            ((long*) field_data_buf)[i] = (((long*) field_data_buf)[i] | ((long)1));
    else EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler software error in perturb_a_field_through_set_last_bit_to_1");
}


void Ensemble_mgt::perturb_a_field_through_set_last_bit_to_0(void *field_data_buf, const char *data_type, long field_size, int current_random_number)
{
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
        for (long i = 0; i < field_size; i ++) 
            ((int*) field_data_buf)[i] = (((int*) field_data_buf)[i] & ((int)0xFFFFFFFE));
    else if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
        for (long i = 0; i < field_size; i ++) 
            ((long*) field_data_buf)[i] = (((long*) field_data_buf)[i] & ((long)0xFFFFFFFFFFFFFFFE));
    else EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler software error in perturb_a_field_through_set_last_bit_to_0");
}


void Ensemble_mgt::perturb_a_field_through_reverse_last_bit(void *field_data_buf, const char *data_type, long field_size, int current_random_number)
{
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
        for (long i = 0; i < field_size; i ++)
            ((int*) field_data_buf)[i] = ((int*) field_data_buf)[i] ^ ((int)1);
    else if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
        for (long i = 0; i < field_size; i ++) 
            ((long*) field_data_buf)[i] = ((long*) field_data_buf)[i] ^ ((long)1);
    else EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler software error in perturb_a_field_through_reverse_last_bit");
}


void Ensemble_mgt::perturb_a_field_through_xor_last_bit_with_a_bit(void *field_data_buf, const char *data_type, long field_size, int current_random_number)
{    
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
        unsigned int number_perturbing_bit = (((unsigned int) current_random_number) >> 1)%32;
        unsigned int perturbing_bitmap = (((unsigned int)1) << number_perturbing_bit);
        for (long i = 0; i < field_size; i ++) {
            unsigned int perturbing_bit_value = ((((unsigned int*)field_data_buf)[i] & perturbing_bitmap) >> number_perturbing_bit);
            if (!(perturbing_bit_value == 0 || perturbing_bit_value == 1))
                EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler software error1 in perturb_a_field_through_xor_last_bit_with_a_bit");
            ((unsigned int*) field_data_buf)[i] = ((unsigned int*) field_data_buf)[i] ^ perturbing_bit_value;
        }
    }
    else if (words_are_the_same(data_type, DATA_TYPE_DOUBLE)) {
        unsigned int number_perturbing_bit = (((unsigned int) current_random_number) >> 1)%64;
        unsigned long perturbing_bitmap = (((unsigned long)1) << number_perturbing_bit);
        for (long i = 0; i < field_size; i ++) {
            unsigned long perturbing_bit_value = ((((unsigned long*)field_data_buf)[i] & perturbing_bitmap) >> number_perturbing_bit);
            if (!(perturbing_bit_value == 0 || perturbing_bit_value == 1))
                EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler software error2 in perturb_a_field_through_xor_last_bit_with_a_bit");
            ((unsigned long*) field_data_buf)[i] = ((unsigned long*) field_data_buf)[i] ^ perturbing_bit_value;
        }
    }
    else EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler software error3 in perturb_a_field_through_xor_last_bit_with_a_bit");
}


void Ensemble_mgt::perturb_an_array(void *field_data_buf, const char *data_type, long field_size, int current_random_number)
{
    if (perturbation_type_id == 1)
        perturb_a_field_through_set_last_bit_to_1(field_data_buf, data_type, field_size, current_random_number);
    else if (perturbation_type_id == 2)
        perturb_a_field_through_set_last_bit_to_0(field_data_buf, data_type, field_size, current_random_number);
    else if (perturbation_type_id == 3)
        perturb_a_field_through_reverse_last_bit(field_data_buf, data_type, field_size, current_random_number);
    else if (perturbation_type_id == 4)
        perturb_a_field_through_xor_last_bit_with_a_bit(field_data_buf, data_type, field_size, current_random_number);
    else EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler software error in Ensemble_mgt::run");
}


void Ensemble_mgt::perturb_a_model_array(void *field_data_buf, const char *data_type, long field_size)
{
    int current_random_number;

    
    if (ensemble_member_id <= 0 || !have_random_seed_for_perturbation)
        return;

    current_random_number = rand();
    if ((current_random_number & (0x000000001)) == 0)
        return;

    EXECUTION_REPORT_LOG(REPORT_LOG,-1, true, "Perturb the values of a model array with random roundoff errors");
    perturb_an_array(field_data_buf, data_type, field_size, current_random_number);
}


void Ensemble_mgt::run()
{
    int current_random_number;

    
    if (ensemble_member_id <= 0 || !have_random_seed_for_perturbation || registered_fields_for_perturbation.size() == 0)
        return;

    current_random_number = rand();
    if ((current_random_number & (0x000000001)) == 0)
        return;

    for (int i = 0; i < registered_fields_for_perturbation.size(); i ++) {
        EXECUTION_REPORT_LOG(REPORT_LOG,-1, true, "Perturb the values of field %s (on grid %s) with random roundoff errors", registered_fields_for_perturbation[i]->get_field_name(), registered_fields_for_perturbation[i]->get_grid_name());
        perturb_an_array(registered_fields_for_perturbation[i]->get_data_buf(), registered_fields_for_perturbation[i]->get_field_data()->get_grid_data_field()->data_type_in_application,
                         registered_fields_for_perturbation[i]->get_size_of_field(), current_random_number);
    }
}

