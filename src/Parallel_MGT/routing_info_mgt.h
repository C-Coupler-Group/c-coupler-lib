/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu and then
  *  modified by Dr. Cheng Zhang and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn or
  *  Dr. Cheng Zhang via zhangc-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef ROUTER_MGT_H
#define ROUTER_MGT_H


#include "common_utils.h" 
#include "decomp_info_mgt.h"
#include "compset_communicators_info_mgt.h"
#include <vector>


struct Routing_info_with_one_process
{
    int remote_proc_global_id;
    int num_elements_transferred;
    int num_local_indx_segments;
    int *local_indx_segment_starts;
    int *local_indx_segment_lengths;
    bool send_or_recv;                           // true is send and false is recv
};


class Routing_info
{
    private:
        int src_comp_id;
        int dst_comp_id;
        char src_comp_full_name[NAME_STR_SIZE];
        char index_dst_comp_full_name[NAME_STR_SIZE];
        char true_dst_comp_full_name[NAME_STR_SIZE];
        char src_decomp_name[NAME_STR_SIZE];
        char dst_decomp_name[NAME_STR_SIZE];
        int src_decomp_size;
        int dst_decomp_size;
        Decomp_info *src_decomp_info;
        Decomp_info *dst_decomp_info;
        Comp_comm_group_mgt_node * src_comp_node;
        Comp_comm_group_mgt_node * dst_comp_node;
        int src_comp_node_id;
        int dst_comp_node_id;
        int current_proc_id_src_comp;
        int current_proc_id_dst_comp;
        char remote_comp_name[NAME_STR_SIZE];
        char local_decomp_name[NAME_STR_SIZE];
        char remote_decomp_name[NAME_STR_SIZE];
        int num_dimensions;
        int total_num_transferred_cells;
        long local_decomp_size;
        long remap_decomp_size;
        std::vector<Routing_info_with_one_process *> recv_from_remote_procs_routing_info;
        std::vector<Routing_info_with_one_process *> send_to_remote_procs_routing_info;

    public:
        Routing_info(const int, const int, const char*, const char*);
        ~Routing_info();
        Routing_info_with_one_process *get_routing_info(bool, int);
        int get_num_elements_transferred_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->num_elements_transferred; }
        int *get_local_indx_segment_starts_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->local_indx_segment_starts; }
        int *get_local_indx_segment_lengths_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->local_indx_segment_lengths; }
        int get_num_local_indx_segments_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->num_local_indx_segments; }
        bool match_router(const int, const int, const char*, const char*);
        int get_num_dimensions() { return num_dimensions; }
        long get_local_decomp_size() { return local_decomp_size; }
        long get_remap_decomp_size() { return remap_decomp_size; }
        Comp_comm_group_mgt_node *get_src_comp_node();
        Comp_comm_group_mgt_node *get_dst_comp_node();
        long get_src_decomp_size() { return src_decomp_size; }
        long get_dst_decomp_size() { return dst_decomp_size; }
        
    private:
        void build_2D_router();
        Routing_info_with_one_process *compute_routing_info_between_decomps(int, const int*, int, const int*, int, int, int, bool);
};


class Routing_info_mgt
{
    private:
        std::vector<Routing_info *> routers;
    
    public:
        Routing_info_mgt() {}
        ~Routing_info_mgt();
        Routing_info *search_router(const int, const int, const char*, const char*);
        Routing_info *search_or_add_router(const int, const int, const char*, const char*);
};

#endif

