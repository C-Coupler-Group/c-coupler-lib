/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_MGT
#define REMAP_MGT


#include "remap_statement_operand.h"
#include <vector>


class Remap_mgt
{
    private:

        struct Remap_statement 
        {
            Remap_statement_operand *result_operand;
            char function[256];
            int line_number;
            std::vector<Remap_statement_operand*> src_operands;
        };

        std::vector<char*> function_words;
        std::vector<char*> reserved_words;
        std::vector<Remap_statement_object*> defined_objects;
        Remap_statement *remap_statement;

        void parse_statement(const int, char **);
        void push_back_words(const char*, const char*);
        void push_back_all_words();
        int match_function_words(const char*);
        int match_reserved_words(const char*);
        int match_defined_objects(const char*);
        void release_statement();
        void process_statement();

    public:
        Remap_mgt(const char *);
        ~Remap_mgt();
};

#endif
