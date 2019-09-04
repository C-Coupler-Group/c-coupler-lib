/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef DECOMP_GRID_MGT
#define DECOMP_GRID_MGT


#include "common_utils.h"
#include "decomp_info_mgt.h"
#include "cor_global_data.h"
#include <vector>


class Decomp_grid_info
{
    private:
        char decomp_name[NAME_STR_SIZE*2];
        int decomp_id;
        Remap_grid_class *decomp_grid;
        Remap_grid_class *original_grid;

    public:
        Decomp_grid_info(int, Remap_grid_class*);
        ~Decomp_grid_info();
        Remap_grid_class *get_decomp_grid() { return decomp_grid; }
        Remap_grid_class *get_original_grid() { return original_grid; }
        const char *get_decomp_name() { return decomp_name; }
        int get_decomp_id() { return decomp_id; }
};


class Decomp_grid_mgt
{
    private:
        std::vector<Decomp_grid_info*> decomp_grids_info;

    public:
        Decomp_grid_mgt() {}
        Decomp_grid_info *search_decomp_grid_info(int, Remap_grid_class*, bool);
        ~Decomp_grid_mgt();
		void set_decomp_grids_using_3D_level_coord(Remap_grid_class *);
		Remap_grid_class *search_decomp_grid_original_grid(int, Remap_grid_class*);
};

#endif
