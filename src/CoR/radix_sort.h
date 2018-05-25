/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef RADIX_SORT_H
#define RADIX_SORT_H


#include "quick_sort.h"
#include <math.h>


template <class T1, class T2> class Radix_sort
{
    public:
        int num_radixes;
        T1 *radix_values[256];
        T2 *content;
        long num_elements;
        T1 tolerable_error;

        Radix_sort(T1**, int, T2*, long, T1);
        void recursively_do_radix_sort(long, long, int);
        void do_radix_sort();
        void do_quick_sort(long, long, int);
        void do_ins_sort(long, long, int);
        long partition(long, long, int, T1);
        void radix_sort_swap(long, long);
        ~Radix_sort();
};


template <class T1, class T2> Radix_sort<T1,T2>::Radix_sort(T1 **radix_values, 
                                                           int num_radixes, 
                                                           T2 *content, 
                                                           long num_elements,
                                                           T1 tolerable_error)
{
    long i, j;


    this->num_radixes = num_radixes;
    this->num_elements = num_elements;
    this->tolerable_error = fabs(tolerable_error);
    this->content = new T2 [num_elements];
    for (i = 0; i < num_elements; i ++)
        this->content[i] = content[i];
    for (i = 0; i < num_radixes; i ++) {
        this->radix_values[i] = new T1 [num_elements];
        for (j = 0; j < num_elements; j ++)
            this->radix_values[i][j] = radix_values[i][j];
    }
}


template <class T1, class T2> Radix_sort<T1,T2>::~Radix_sort()
{
    delete [] content;
    for(int i = 0; i < num_radixes; i++)
        delete [] radix_values[i];
}


template <class T1, class T2> void Radix_sort<T1,T2>::do_radix_sort()
{
    recursively_do_radix_sort(0, num_elements-1, 0);
}


template <class T1, class T2> void Radix_sort<T1,T2>::recursively_do_radix_sort(long segment_start, long segment_end, int radix_id)
{
    long i;
    long new_segment_start, new_segment_end;
    T1 last_radix_value;


    if (radix_id >= num_radixes)
        return;

    do_quick_sort(segment_start, segment_end, radix_id);

    last_radix_value = radix_values[radix_id][segment_start];
    new_segment_start = segment_start;
    for (i = segment_start+1; i < segment_end; i ++) {
        if (fabs(radix_values[radix_id][i] - last_radix_value) > tolerable_error) {
            new_segment_end = i - 1;
            recursively_do_radix_sort(new_segment_start, new_segment_end, radix_id+1);
            new_segment_start = i;
            last_radix_value = radix_values[radix_id][i];
        }
    }
    new_segment_end = segment_end;
    recursively_do_radix_sort(new_segment_start, new_segment_end, radix_id+1);
}


template <class T1, class T2> void Radix_sort<T1,T2>::radix_sort_swap(long pos1, long pos2)
{
    for (int k = 0; k < num_radixes; k ++)
        swap(radix_values[k]+pos1, radix_values[k]+pos2);
    swap(content+pos1, content+pos2);
}


template <class T1, class T2> void Radix_sort<T1,T2>::do_quick_sort(long segment_start, long segment_end, int radix_id)
{
    long pivotindex, partition_pos;
    
    
    if (segment_end-segment_start < QUICK_SORTHRESH) 
        do_ins_sort(segment_start, segment_end, radix_id);
    else {
        pivotindex = (segment_start+segment_end) / 2;
        radix_sort_swap(pivotindex, segment_end);
        partition_pos = partition(segment_start-1, segment_end, radix_id, radix_values[radix_id][segment_end]);
        radix_sort_swap(partition_pos, segment_end);
        if(partition_pos - segment_start > 1) 
            do_quick_sort(segment_start, partition_pos-1, radix_id);
        if(segment_end-partition_pos > 1) 
            do_quick_sort(partition_pos+1, segment_end, radix_id);
    }
}


template <class T1, class T2> void Radix_sort<T1,T2>::do_ins_sort(long segment_start, long segment_end, int radix_id)
{
    for (long i = segment_start; i <= segment_end; i ++)
        for (long j = i; j > segment_start && (radix_values[radix_id][j] > radix_values[radix_id][j-1]); j--)
            radix_sort_swap(j, j-1);
}


template <class T1, class T2> long Radix_sort<T1,T2>::partition(long segment_start, long segment_end, int radix_id, T1 pivot)
{
    do {
        while (radix_values[radix_id][++segment_start] > pivot);
        while (segment_end > 0 && radix_values[radix_id][--segment_end] < pivot);
        radix_sort_swap(segment_start, segment_end);
    } while (segment_start < segment_end);
    radix_sort_swap(segment_start, segment_end);
    return segment_start;
}


#endif

