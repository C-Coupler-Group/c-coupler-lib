/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef DICTIONARY_H
#define DICTIONARY_H


#include "stdlib.h"
#include <string.h>
#include <stdio.h>
#include "execution_report.h"
#include "remap_common_utils.h"



#define INITIAL_SIZE (1024)
#define GROWTH_FACTOR (2)
#define MAX_LOAD_FACTOR (2)
#define MULTIPLIER (97)


template <class T>
struct Dictionary_node {
    struct Dictionary_node *next;
    char *key;
    T value;
};


template <class T>
class Dictionary
{
    private:
        int num_hashing_entries;
        int num_elements;
        struct Dictionary_node<T> **hashing_table;
        
        unsigned long hash_function(const char*);    
        void initialize_hashing_table(int);
        void delete_hashing_table();
        void increase_hashing_entries();

    public:
        Dictionary(int);
        ~Dictionary();
        T search(const char*, bool);
        void insert(const char*, T);
        void remove(const char*);
};


template <class T>  
Dictionary<T>::Dictionary(int size)
{
    int true_size = size / 4;
  
  
    if (true_size < INITIAL_SIZE)
        true_size = INITIAL_SIZE;
    if (true_size > size)
        true_size = size;
  
    num_elements = 0;
    initialize_hashing_table(true_size);
}
  

template <class T>  
void Dictionary<T>::initialize_hashing_table(int size)
{
    num_hashing_entries = size;
    hashing_table = (Dictionary_node<T>**) new Dictionary_node<T> *[size];
    for(int i = 0; i < size; i++) 
        hashing_table[i] = NULL;
}


template <class T>
void Dictionary<T>::delete_hashing_table()
{
    int i;
    Dictionary_node<T> *e, *next;


    for(i = 0; i < num_hashing_entries; i++) {
        for(e = hashing_table[i]; e != NULL; e = next) {
            next = e->next;
            delete [] e->key;
            delete e;
        }
    }
  
    delete [] hashing_table;
}
  

template <class T>  
Dictionary<T>::~Dictionary()
{
    delete_hashing_table();
}
  

template <class T>
unsigned long Dictionary<T>::hash_function(const char *keyword)
{
    unsigned const char *us;
    unsigned long h;
  
    h = 0;
  
    for(us = (unsigned const char *) keyword; *us; us++) {
        h = (h << 6) + (h << 16) - h + *us;
    }

    return h % num_hashing_entries;
}
  

template <class T>
void Dictionary<T>::increase_hashing_entries()
{   
    int old_num_hashing_entries = num_hashing_entries, new_num_hashing_entries;
    Dictionary_node<T> **old_hashing_table = hashing_table, **new_hashing_table;
    Dictionary_node<T> *e, *next;
      
  
    initialize_hashing_table(num_hashing_entries * GROWTH_FACTOR);
    new_hashing_table = hashing_table;
    new_num_hashing_entries = num_hashing_entries;
    num_elements = 0;
  
    for(int i = 0; i < old_num_hashing_entries; i++) {
        for(e = old_hashing_table[i]; e != NULL; e = e->next)
            insert(e->key, e->value);
    }
  
    hashing_table = old_hashing_table;
    num_hashing_entries = old_num_hashing_entries;
    delete_hashing_table();
    hashing_table = new_hashing_table;
    num_hashing_entries = new_num_hashing_entries;
}
  

template <class T>
void Dictionary<T>::insert(const char *key, T value)
{
    Dictionary_node<T> *e;
    unsigned long h;
      
  
    EXECUTION_REPORT(REPORT_ERROR, !words_are_the_same(key,""), "The key to be inserted into the dictionary cannot be empty!");
    EXECUTION_REPORT(REPORT_ERROR, search(key, false) == NULL, "The key \"%s\" has been inserted into the dictionary before. It cannot be inserted again. Please check.");
    e = new Dictionary_node<T>;
    e->key = strdup(key);
    e->value = value;
    h = hash_function(key);
    e->next = hashing_table[h];
    hashing_table[h] = e;
    num_elements ++;
  
    /* increase_hashing_entries hashing_table if there is not enough room */
    if(num_elements >= num_hashing_entries * MAX_LOAD_FACTOR)
    increase_hashing_entries();
}

  
/* return the most recently inserted value associated with a key */
/* or 0 if no matching key is present */
template <class T>
T Dictionary<T>::search(const char *key, bool check)
{
    Dictionary_node<T> *e;
  
  
    for(e = hashing_table[hash_function(key)]; e != NULL; e = e->next) {
        if(strcmp(e->key, key) == 0) {
            /* got it */
            return e->value;
        }
    }
  
    if (check)
        EXECUTION_REPORT(REPORT_ERROR, false, "Cannot find the entry for the keyword \"%s\" in the dictionary", key);
  
    return (T) 0;
}
  

template <class T>
void Dictionary<T>::remove(const char *key)
{
    Dictionary_node<T> **prev, *e;  
  
  
    for(prev = &(hashing_table[hash_function(key)]); *prev != 0; prev = &((*prev)->next)) {
        if(strcmp((*prev)->key, key) == 0) {
            /* got it */
            e = *prev;
            *prev = e->next;
            delete [] e->key;
            delete e;
            num_elements --;
            break;
        }
    }
}
  


#endif

