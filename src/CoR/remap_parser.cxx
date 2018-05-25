/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "remap_parser.h"
#include "parse_special_words.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>


void Remap_parser::push_back_words(const char *type, const char *word)
{
    char *localized_word;

    
    localized_word = new char [256];
    strcpy(localized_word, word);
    if (words_are_the_same(type, "seperator"))
        seperator_words.push_back(localized_word);
    else if (words_are_the_same(type, "reserved")) 
        reserved_words.push_back(localized_word);
}


Remap_parser::Remap_parser(const char* script_file_name)
{
    if((script_fp = fopen(script_file_name, "rb")) == NULL){
        printf("Config file %s can not be opened\n", script_file_name);
        assert(false);
    }     

    /* Push back all seperator words */
    push_back_words("seperator", SEPERATE_WORD_SPACE);
    push_back_words("seperator", SEPERATE_WORD_TAB);
    push_back_words("seperator", RESERVED_WORD_QUOTE_MARK);
    push_back_words("seperator", RESERVED_WORD_EQUAL);
    push_back_words("seperator", RESERVED_WORD_LEFT_BRACKET);
    push_back_words("seperator", RESERVED_WORD_RIGHT_BRACKET);
    push_back_words("seperator", RESERVED_WORD_COMMA);
    push_back_words("seperator", RESERVED_WORD_ATTRIBUTE);

    /* Push back all reserved words */
    push_back_words("reserved", RESERVED_WORD_QUOTE_MARK);
    push_back_words("reserved", RESERVED_WORD_EQUAL);
    push_back_words("reserved", RESERVED_WORD_LEFT_BRACKET);
    push_back_words("reserved", RESERVED_WORD_RIGHT_BRACKET);
    push_back_words("reserved", RESERVED_WORD_COMMA);
    push_back_words("reserved", RESERVED_WORD_ATTRIBUTE);
}


Remap_parser::~Remap_parser()
{
    for (int i = 0; i < reserved_words.size(); i ++)
        delete [] reserved_words[i];
    for (int i = 0; i < seperator_words.size(); i ++)
        delete [] seperator_words[i];
}


/* Get the current line in the script file. Return false when touching the end of file */
bool Remap_parser::get_next_line_of_script_file(char *line, FILE *fp)
{
    char c;
    int iter = 0;
    

    while (!feof(fp) && (c = getc(fp)) != -1) {
        if (c == '\n') 
            break;
        line[iter ++] = c;
    }
    line[iter ++] = '\0';
    if (feof(fp))
        return false;

    return true;
}


/* Get the current word in the current line. Words are seperated by seperator words, reserved words, 
    space word and tab word. Reserved words will be recorded as a enssential word while other kinds of 
    words for seperating will be neglected. Annotation words are used to mark the notes in a line, they 
    and the characters after them will be neglected. */
bool Remap_parser::get_next_word_in_line(char *essential_word, char **line)
{
    char *word = essential_word;
    int matched_seperator_word_id;
    int matched_reserved_word_id;


    if ((*line)[0] == '\0' || (*line)[0] == RESERVED_WORD_ANNOTATION) {
        (*line) = NULL;
        return false;
    }

    while (1) {
        matched_seperator_word_id = match_seperator_words(*line);
        matched_reserved_word_id = match_reserved_words(*line); 
        if (matched_reserved_word_id != -1) {
            strcpy(essential_word, reserved_words[matched_reserved_word_id]);
            *line = *line + strlen(reserved_words[matched_reserved_word_id]); 
            return true;
        }
        if (matched_seperator_word_id == -1)
            break;
        *line = *line + strlen(seperator_words[matched_seperator_word_id]);
    }    

    while ((*line)[0] != '\0' && match_seperator_words(*line) == -1) {
        *essential_word = (*line)[0];
        essential_word ++;
        (*line) ++;
    }
    
    *essential_word = '\0';

    return strlen(word) != 0;
}


/* Search a seperator word which is the same with a prefix of given statement string. When there are
    several seperator words match some prefixes of statment string, the longest seperator will be returned */
int Remap_parser::match_seperator_words(const char *statement_string)
{
    int i;
    int matched_id;


    for (i = 0, matched_id = -1; i < seperator_words.size(); i ++) {
        if (strncmp(seperator_words[i], statement_string, strlen(seperator_words[i])) == 0) {
            if (matched_id == -1)
                matched_id = i;
            else if (strlen(seperator_words[matched_id]) < strlen(seperator_words[i]))
                matched_id = i;
        }
    }

    return matched_id;
}


/* Search a reserved word which is the same with a prefix of given statement string. When there are
    several reserved words match some prefixes of statment string, the longest seperator will be returned */
int Remap_parser::match_reserved_words(const char *statement_string)
{
    int i;
    int matched_id;


    for (i = 0, matched_id = -1; i < reserved_words.size(); i ++) {
        if (strncmp(reserved_words[i], statement_string, strlen(reserved_words[i])) == 0) {
            if (matched_id == -1)
                matched_id = i;
            else if (strlen(reserved_words[matched_id]) < strlen(reserved_words[i]))
                matched_id = i;
        }
    }

    return matched_id;
}


/* Get next statement and check syntax at the same time */
bool Remap_parser::get_next_parsed_statement(int *num_words, char **words)
{
    char statement_string[NAME_STR_SIZE];
    char *line_iter;


    if (!get_next_line_of_script_file(statement_string, script_fp)) 
        return false;

    line_iter = statement_string;
    *num_words = 0;
    while (get_next_word_in_line(words[*num_words], &line_iter)) 
        (*num_words) ++;
    return true;
}
