#include <iostream>
#include <vector>
#include <mpi.h>
#include <assert.h>
#include <iomanip>


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



int printvec2d(std::vector<int> vec, int rows, int cols) {
    for (int i=0; i<rows+2; ++i) {
        for (int j=0; j<cols; ++j) {
            std::cout << std::setw(3) << vec[i*cols + j];// << " ";
        }
        std::cout << "\n";
    }
}


//int printarr2d(int* arr, int rows, int cols) {
//    for (int i=0; i<rows; ++i) {
//        for (int j=0; j<cols; ++j) {
//            std::cout << arr[i*cols+j] << " ";
//        }
//        std::cout << "\n";
//    }
//}


int main(int argc, char* argv[])
{

    int N = 6;
    
    // initialize MPI domain
    MPI_Init(&argc, &argv);
    
    world_info world;
    
    MPI_Comm_size(MPI_COMM_WORLD, &world.size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world.rank);
    
    int dims[2] = {0,1};
    MPI_Dims_create(world.size, 2, dims);
    world.dims_x = dims[0];
    world.dims_y = dims[1];
    
    if (world.rank == 0)
        std::cout
        << "processes: " << world.size << "\n"
        << "dims_x: " << world.dims_x << "\n"
        << "dims_y: " << world.dims_y << "\n"
        << std::endl;
    
//    if (world.rank == 0) {
//        printf("processes: %d\ndims_x: %d\ndims_y: %d\n", world.size, world.dims_x, world.dims_y);
//    }
    
    int tmp = N;
    N = tmp % world.dims_x == 0 ? tmp : tmp + (world.dims_x - tmp % world.dims_x);
    assert(N % world.dims_x == 0);
    
    
    
    
    MPI_Comm cart_comm;
    
    // build periodic process geometry with cartesian communicator
    int periods[2] = {false, false};
    int dims_grid[2] = {world.dims_x, world.dims_y};
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims_grid, periods, true, &cart_comm);
    
    MPI_Comm_rank(cart_comm, &world.cart_rank);
    
    MPI_Cart_shift(cart_comm, 0, 1, &world.top_proc, &world.bottom_proc);
    
    int coords[2];
    MPI_Cart_coords(cart_comm, world.cart_rank, 2, coords);
    
    world.coord_x = coords[0];
    world.coord_y = coords[1];
    
    // global grid
    int Nx_glo = N;
    int Ny_glo = N;
    int NN_glo = Nx_glo * Ny_glo;
    
    // local grid
    int Nx_loc = Nx_glo / world.dims_x;
    int Ny_loc = Ny_glo / world.dims_y;
    int NN_loc = Nx_loc * Ny_loc;
    
    
    // datatypes
    MPI_Datatype bottom_boundary, top_boundary, block, block_resized;
    
    // build contiguous vectors for boundaries (rows)
    // each process has multiple rows in the grid
    MPI_Type_contiguous(Ny_loc, MPI_INT, &bottom_boundary);
    MPI_Type_commit(&bottom_boundary);
    
    MPI_Type_contiguous(Ny_loc, MPI_INT, &top_boundary);
    MPI_Type_commit(&top_boundary);
    
    // build blocks for transpose
    int sizes[2]    = {Nx_loc, Ny_loc}; // size of global vector
    int subsizes[2] = {Nx_loc, Nx_loc}; // size of sub-region (square)
    int starts[2]   = {0,0};
    
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &block);
    // or use vector -> test
//    MPI_Type_vector(Ny_loc, Nx_loc, Ny_loc, MPI_INT, &block);
    
    MPI_Type_commit(&block);

    // resize data structure, so that it is contagious (for alltoall)
    MPI_Type_create_resized(block, 0, Nx_loc*sizeof(int), &block_resized);
    MPI_Type_free(&block);
    MPI_Type_commit(&block_resized);
    
    

    // initialize data structure
    std::vector<int> send((Nx_loc+2)*Ny_loc,0);
    std::vector<int> recv((Nx_loc+2)*Ny_loc,0);
    
    int ind = 0;
    for (int i=0; i<Nx_loc+2; ++i) {
        for (int j=0; j<Ny_loc; ++j) {
            send[i*Ny_loc+j] = (world.coord_x*Nx_loc+i)*Ny_glo + j;
        }
    }

    // ensure that all procs finished with the initialization    
    MPI_Barrier(MPI_COMM_WORLD);

    
    




    // transpose globally block-wise
    MPI_Alltoall(&send[Ny_loc], 1, block_resized, &recv[Ny_loc], 1, block_resized, MPI_COMM_WORLD);
    
    // transpose globally block-wise
//    MPI_Alltoall(&send[0], 1, block_resized, &recv[0], 1, block_resized, MPI_COMM_WORLD);
    
    // locally transpose blocks, WORKS as expected
    // loop over blocks TODO parallelize block loop with openmp
    int Nb_loc = Ny_loc/Nx_loc;
    int ind1, ind2;
    int tmp2;
    for (int b=0; b<Nb_loc; ++b) {
        for (int i=0; i<Nx_loc; ++i) {
            for (int j=0; j<i; ++j) {
                ind1 = (i+1)*Ny_loc + j + b*Nx_loc; // regular index + offset of block
                ind2 = (j+1)*Ny_loc + i + b*Nx_loc; // switch i and j
                
                std::swap(recv[ind1], recv[ind2]);
            }
        }
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    // exchange boundaries
    MPI_Request request[4];
    MPI_Status status[4];
    
    if (world.coord_x % 2 == 0) { // first send top, then botton
        MPI_Isend(&recv[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, 0, cart_comm, &request[0]);
        MPI_Irecv(&recv[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, 0, cart_comm, &request[1]);
        MPI_Isend(&recv[(Ny_loc)],            1, top_boundary,    world.top_proc,    0, cart_comm, &request[2]);
        MPI_Irecv(&recv[0],                   1, top_boundary,    world.top_proc,    0, cart_comm, &request[3]);
    }
    else { // first send botton, then top
        
        MPI_Irecv(&recv[0],                   1, top_boundary,    world.top_proc,    0, cart_comm, &request[0]);
        MPI_Isend(&recv[(Ny_loc)],            1, top_boundary,    world.top_proc,    0, cart_comm, &request[1]);
        MPI_Irecv(&recv[(Nx_loc+1)*(Ny_loc)], 1, bottom_boundary, world.bottom_proc, 0, cart_comm, &request[2]);
        MPI_Isend(&recv[(Nx_loc)*(Ny_loc)],   1, bottom_boundary, world.bottom_proc, 0, cart_comm, &request[3]);
    }
    
    MPI_Waitall(4,request,status);
    
    
    // print result
    for (int i=0; i<world.dims_x; ++i) {
        if (world.rank != i) continue;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "\nSend - Proc " << i << "\n";
        printvec2d(send,Nx_loc,Ny_loc);
    }
    for (int i=0; i<world.dims_x; ++i) {
        if (world.rank != i) continue;
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "\nRecv - Proc " << i << "\n";
        printvec2d(recv,Nx_loc,Ny_loc);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Finalize();
    
    return 0;
}
