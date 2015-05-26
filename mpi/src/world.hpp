#ifndef WORLD_HPP
#define WORLD_HPP

struct world_info
{
    int size;
    int dims_x;
    int dims_y;
    
    int left_proc;
    int right_proc;
    int top_proc;
    int bottom_proc;
    
    int rank;
    int cart_rank;
    int coord_x;
    int coord_y;
};

#endif // WORLD_HPP
