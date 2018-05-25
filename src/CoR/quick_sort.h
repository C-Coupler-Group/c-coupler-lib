/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef QUICK_SORT_H
#define QUICK_SORT_H


#define QUICK_SORTHRESH 4


template <class T> void swap(T* k, T* j)
{
    T temp;
    temp = *j;
    *j = *k;
    *k = temp;
}


template <class T1, class T2> void do_ins_sort(T1 *sorted_values,
                                               T2 *content_values,
                                               long segment_start, 
                                               long segment_end)
{
    long i, j;
    for (i = segment_start+1; i <= segment_end; i ++) {
        for (j = i; j > segment_start && (sorted_values[j] < sorted_values[j-1]); j--) {
            swap(sorted_values+j, sorted_values+j-1);
            if (content_values != NULL)
                swap(content_values+j, content_values+j-1);
        }
        if (content_values != NULL)
            for (; j > segment_start && (sorted_values[j] == sorted_values[j-1]) && (content_values[j] < content_values[j-1]); j--)
                swap(content_values+j, content_values+j-1);
    }
}


template <class T1, class T2> long partition(T1 *sorted_values,
                                             T2 *content_values,
                                             long segment_start, 
                                             long segment_end, 
                                             T1 pivot)
{
    do {
        while (sorted_values[++segment_start] < pivot);
        while (segment_end > 0 && sorted_values[--segment_end] > pivot);
        swap(sorted_values+segment_start, sorted_values+segment_end);
        if (content_values != NULL)
            swap(content_values+segment_start, content_values+segment_end);
    } while (segment_start < segment_end);
    swap(sorted_values+segment_start, sorted_values+segment_end);
    if (content_values != NULL)
        swap(content_values+segment_start, content_values+segment_end);

    return segment_start;
}


template <class T1, class T2> void do_quick_sort(T1 *sorted_values,
                                                 T2 *content_values,
                                                 long segment_start, 
                                                 long segment_end)
{
    long pivotindex, partition_pos, i, j;
    
    
    if (segment_end-segment_start < QUICK_SORTHRESH) 
        do_ins_sort(sorted_values, content_values, segment_start, segment_end);
    else {
        pivotindex = (segment_start+segment_end) / 2;
        swap(sorted_values+pivotindex, sorted_values+segment_end);
        if (content_values != NULL)
            swap(content_values+pivotindex, content_values+segment_end);
        partition_pos = partition(sorted_values, content_values, segment_start-1, segment_end, sorted_values[segment_end]);
        swap(sorted_values+partition_pos, sorted_values+segment_end);
        if (content_values != NULL)
            swap(content_values+partition_pos, content_values+segment_end);
        if(partition_pos - segment_start > 1) 
            do_quick_sort(sorted_values, content_values, segment_start, partition_pos-1);
        if(segment_end-partition_pos > 1) 
            do_quick_sort(sorted_values, content_values, partition_pos+1, segment_end);
        if (content_values != NULL) {
            for (i = partition_pos-1; i >= segment_start && sorted_values[i] == sorted_values[partition_pos]; i --);
            for (j = partition_pos+1; j <= segment_end && sorted_values[j] == sorted_values[partition_pos]; j ++);
            if ((++i) < (--j))
                do_quick_sort(content_values, (long*)NULL, i, j);
        }
    }
}


#endif

