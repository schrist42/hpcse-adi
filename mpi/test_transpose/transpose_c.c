/* file: transpose.c 
 * 
 * Description: Transpose operation using MPI derived datatypes with 
 * MPI collective communications. 
 * 
 * Author: Spenser Gilliland <spenser_at_[hidden]> 
 */ 
#include <mpi.h> 
#include <stdio.h> 
#include <unistd.h> 
#define N (8) 
float matrix[N][N]; 
void print_matrix(int wrank, int wrows) { 
    int i , j; 
    MPI_Barrier(MPI_COMM_WORLD); 
    if(wrank == 0) printf("Matrix = \n"); 
    MPI_Barrier(MPI_COMM_WORLD); 
    for(i = 0; i < N; i++) { 
        if(i >= wrank*wrows && i < (wrank+1)*wrows) { 
            printf("%2d:", wrank); 
            for(j = 0; j < N; j++) { 
                printf("%6.2g", matrix[i][j]); 
            } 
            printf("\n"); 
        } 
        usleep(1); 
        MPI_Barrier(MPI_COMM_WORLD); 
    } 
} 
int main(int argc, char *argv[]) { 
    int wsize, wrank, wrows; 
    int i, j, k; 
    int row, col; 
    float temp; 
    MPI_Datatype mpi_all_unaligned_t, mpi_all_t; 
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &wsize); 
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank); 
    wrows = N / wsize; 
    MPI_Type_vector(N, wrows, N, MPI_FLOAT, &mpi_all_unaligned_t); 
    MPI_Type_create_resized(mpi_all_unaligned_t, 0, wrows*sizeof(float), &mpi_all_t); 
    MPI_Type_free(&mpi_all_unaligned_t); 
    MPI_Type_commit(&mpi_all_t); 
    for(i = 0; i < N; i++) { 
        /* Initialize data on the rows of the matrix owned by this rank */ 
        if (i >= wrank*wrows && i < (wrank+1)*wrows) { 
            for(j = 0; j < N; j++) { 
                matrix[i][j] = i*N + j; 
            } 
        } 
    } 
    print_matrix(wrank, wrows); 
    /* Local Transpose */ 
    row = wrank*wrows; 
    for(k = 0; k < wsize; k++) { 
        col = k*wrows; 
        for( i = 0; i < wrows; i++) { 
            for(j =i+1; j < wrows; j++) { 
                temp = matrix[row+i][col + j]; 
                matrix[row + i][col + j] = matrix[row + j][col + i]; 
                matrix[row + j][col + i] = temp; 
            } 
        } 
    } 
    print_matrix(wrank, wrows); 
    /* Global Transpose */ 
    MPI_Alltoall(matrix[wrank*wrows], 1, mpi_all_t, 
                 matrix[wrank*wrows], 1, mpi_all_t, 
                 MPI_COMM_WORLD); 
    print_matrix(wrank, wrows); 
    MPI_Finalize(); 
    return 0; 
} 
