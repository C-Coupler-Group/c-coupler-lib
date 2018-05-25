/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef ENSEMBLE_MGT_H
#define ENSEMBLE_MGT_H


#include <vector>
#include "memory_mgt.h"


class Ensemble_mgt
{
    private:
        int ensemble_member_id;
        bool have_random_seed_for_perturbation;
        int root_random_seed_for_perturbation;
        int ensemble_random_seed_for_perturbation;
        int perturbation_type_id;
        std::vector<Field_mem_info *> registered_fields_for_perturbation;
        
        void perturb_a_field_through_set_last_bit_to_1(void*, const char*, long, int);
        void perturb_a_field_through_set_last_bit_to_0(void*, const char*, long, int);
        void perturb_a_field_through_reverse_last_bit(void*, const char*, long, int);
        void perturb_a_field_through_xor_last_bit_with_a_bit(void*, const char*, long, int);
        void perturb_an_array(void*, const char*, long, int);

    public:
        Ensemble_mgt();
        void Initialize(int, int, int, const char*);
        ~Ensemble_mgt() {}
        void register_a_field_for_perturbation(void *);
        void perturb_fields_with_roundoff_errors();
        void perturb_a_model_array(void*, const char*, long);
        void run();
};


#endif
