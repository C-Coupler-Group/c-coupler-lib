/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef REMAP_PARSER
#define REMAP_PARSER


#include <vector>
#include <stdio.h>


class Remap_parser
{
    private:
        FILE *script_fp;
        std::vector<char*> seperator_words;
        std::vector<char*> reserved_words;

        void push_back_words(const char*, const char*);
        bool get_next_line_of_script_file(char*, FILE*);
        bool get_next_word_in_line(char*, char**);
        int match_seperator_words(const char*);
        int match_reserved_words(const char*);

    public: 
        Remap_parser(const char *);
        ~Remap_parser();
        bool get_next_parsed_statement(int*, char**);
};

#endif
