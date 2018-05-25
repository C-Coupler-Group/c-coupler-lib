/***************************************************************
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef GRID_CELL_SEARCH_ENGINE
#define GRID_CELL_SEARCH_ENGINE


#include "remap_grid_class.h"


#define TILE_DIVIDE_FACTOR         2
#define MAX_NUM_CELLS_IN_TILE      8

#define EDGE_TYPE_LATLON           1
#define EDGE_TYPE_GREAT_ARC        2
#define EDGE_TYPE_XY               3       // not supported currently


class H2D_grid_cell_search_cell;

struct H2D_grid_cell_cartesian_coord
{
    protected:
        friend class H2D_grid_cell_search_cell;
        double center_x;
        double center_y;
        double center_z;
        double *vertex_x;
        double *vertex_y;
        double *vertex_z;

    public:
        H2D_grid_cell_cartesian_coord();
        ~H2D_grid_cell_cartesian_coord();
};


class H2D_grid_cell_search_cell
{
    private:
        int cell_index;
        double center_lon;
        double center_lat;
        int num_vertex;
        double *vertex_lons;
        double *vertex_lats;
        double bounding_circle_center_lon;
        double bounding_circle_center_lat;
        double bounding_circle_radius;
        bool mask;
        int edge_type;   // 1: latlon edge; 2: great arc edge
        H2D_grid_cell_cartesian_coord *cartesian_coord;

        bool check_overlapping_for_latlon_grid(const H2D_grid_cell_search_cell*) const;
        bool check_overlapping_for_cartesian_grid(const H2D_grid_cell_search_cell*) const;
        bool check_possible_overlapping(const H2D_grid_cell_search_cell*) const;
        bool is_point_in_latlon_coord_cell(double, double) const;
        bool is_point_in_cartesian_coord_cell(double, double) const;

    public:
        H2D_grid_cell_search_cell(int, double, double, bool, int, const double*, const double*, int);
        H2D_grid_cell_search_cell() {}
        ~H2D_grid_cell_search_cell();
        double get_center_lon() const { return center_lon; }
        double get_center_lat() const { return center_lat; }
        int get_num_vertex() const { return num_vertex; }
        int get_cell_index() const { return cell_index; }
        double get_vertex_lon(int indx) const ;
        double get_vertex_lat(int indx) const ;
        bool get_mask() const { return mask; }
        void set_mask(bool new_mask) { mask = new_mask; } 
        int get_edge_type() const { return edge_type; }
        bool check_overlapping(const H2D_grid_cell_search_cell*, bool) const;
        double get_bounding_circle_radius() const { return bounding_circle_radius; }
        double get_bounding_circle_center_lon() const { return bounding_circle_center_lon; }
        double get_bounding_circle_center_lat() const { return bounding_circle_center_lat; }
};


class H2D_grid_cell_search_tile
{
    private:
        int num_cells;
        H2D_grid_cell_search_cell **cells;
        H2D_grid_cell_search_cell **cells_buffer;
        H2D_grid_cell_search_tile *children[TILE_DIVIDE_FACTOR*TILE_DIVIDE_FACTOR];
        H2D_grid_cell_search_tile *parent;
        long *index_buffer;
        double center_lon;
        double center_lat;
        double circle_radius;
        double dlon;
        double dlat;

    public:
        bool has_cell_index(int);
        H2D_grid_cell_search_tile(int, H2D_grid_cell_search_cell**, H2D_grid_cell_search_cell**, long*, H2D_grid_cell_search_tile*, double, double, double, double);
        ~H2D_grid_cell_search_tile();
        void divide_tile();
        void compute_bounding_circle();
        bool search_points_within_distance(double, double, double, int&, long*, double*, bool);
        void search_overlapping_cells(int&, long*, const H2D_grid_cell_search_cell*, bool, bool);
};


class H2D_grid_cell_search_engine
{
    private:
        const Remap_grid_class *remap_grid;
        H2D_grid_cell_search_cell **cells;
        H2D_grid_cell_search_cell **cells_ptr;
        H2D_grid_cell_search_cell **cells_buffer;
        long *index_buffer;
        double *dist_buffer;
        H2D_grid_cell_search_tile *root_tile;
        double dist_threshold;
        int num_cells;
        
    public:
        H2D_grid_cell_search_engine(const Remap_grid_class*, const double*, const double*, const bool*, const bool*, int, const double*, const double*, int, bool);
        ~H2D_grid_cell_search_engine();
        void recursively_search_initial_boundary(double, double, double, double, double&, double&, double&, double&);
        void search_nearest_points_var_number(int, double, double, int&, long*, double *, bool);
        void search_nearest_points_var_distance(double, double, double, int&, long*, double*, bool);
        void search_overlapping_cells(int&, long*, const H2D_grid_cell_search_cell*, bool, bool) const;
        int search_cell_of_locating_point(double, double, bool) const;
        const H2D_grid_cell_search_cell* get_cell(int) const;
        void update(const bool*);
};

#endif

