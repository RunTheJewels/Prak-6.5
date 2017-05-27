#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>
#include "mpi.h"
#include <omp.h>
#include <stdio.h>

using namespace std;
typedef std::complex<float> complexd;
uint threads;

complexd operator*(complexd a, complexd b)
{
  return complexd(a.real()*b.real()-a.imag()*b.imag(),a.real()*b.imag()+a.imag()*b.real());
}

const complexd pi = (3.14159265358,0);

int rank = 0, comm_size;

void cubit(vector<complexd>& a, vector<complexd>& b, int loc_size, int N, int K, complexd H[2][2]);
void cubit2(vector<complexd>& a, vector<complexd>& b, int loc_size, int N, int K1, int K2, complexd H[4][4]);

complexd AD[2][2] = {0.707107,0.707107,0.707107,-0.707107};
complexd H[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
complexd swapqbits[4][4] = {{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};

void qft(vector<complexd>& a, vector<complexd>& b, int loc_size, int N)
{

  for(uint i=1; i<=N; ++i) {
      for(uint j=1; j<i; ++j) {
          H[3][3] = exp(complexd(0,1) * pi / pow((float)2, (float)j));
          cubit2(a, b, loc_size, N, i, j, H);
          a.swap(b);
      }
      cubit(a, b,loc_size, N, i, AD);
      a.swap(b);
  }

  for(uint i=1; i<N/2; ++i) {
      cubit2(a, b, loc_size, N, i, N-i+1, swapqbits);
      a.swap(b);
  }
}


void cubit(vector<complexd>& a, vector<complexd>& b, int loc_size, int N, int K, complexd H[2][2])
{
    int P = N - K;
    int stride = 1 << P;
    if (stride < loc_size)
    {
        #pragma omp parallel for num_threads(threads) shared(a,b)
        for (int i = 0; i < loc_size; ++i)
        {
            int j_1 = i & ~stride;
            int j_2 = i | stride;
            int u_i = !((i & stride) == 0);
            b[i] = (H[u_i][0] * a[j_1]) + (H[u_i][1] * a[j_2]);
            printf("b[%d] = H[%d][0] * a[%d] + H[%d][1] * a[%d] = ",i,u_i,j_1,u_i,j_2);
            cout << H[u_i][0] << " * " << a[j_1] << " + " << H[u_i][1] << " * " <<
             a[j_2] << " = " << b[i] << endl;
        }
    } else
    {
        int proc_stride = stride / loc_size;
        vector<complexd> tmp(loc_size);

        MPI_Send(a.data(),loc_size,MPI::COMPLEX,proc_stride^rank,0,MPI_COMM_WORLD);

        MPI_Recv(tmp.data(),loc_size,MPI::COMPLEX,proc_stride^rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < loc_size; ++i)
        {
            if (!(rank & proc_stride))
                b[i] = H[0][0] * a[i] + H[0][1] * tmp[i];
            else
                b[i] = -(H[1][0] * a[i] + H[1][1] * tmp[i]);
        }
    }
}

void cubit2(vector<complexd>& a, vector<complexd>& b, int loc_size, int N, int K1, int K2, complexd H[4][4])
{
    int P1 = N - K1;
    int stride1 = 1 << P1;
    int P2 = N - K2;
    int stride2 = 1 << P2;

    if (stride1 < loc_size and stride2 < loc_size)
    {
        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < loc_size; ++i)
        {
            int i00 = i & ~stride1 & ~stride2;
            int i01 = i & ~stride1 | stride2;
            int i10 = (i | stride1) & ~stride2;
            int i11 = i | stride1 | stride2;
            int ik1 = (i & stride1) >> P1;
            int ik2 = (i & stride2) >> P2;
            int ik = (ik1 << 1) + ik2;
            b[i] = H[ik][0] * a[i00] + H[ik][1] * a[i01] +
            H[ik][2] * a[i10] + H[ik][3] * a[i11];
        }
    } else if (stride1 < loc_size)
    {
        int proc_stride = stride1 / loc_size;
        vector<complexd> tmp(loc_size);;

        MPI_Sendrecv(a.data(),loc_size,MPI::COMPLEX,proc_stride^rank,0,tmp.data(),
        loc_size,MPI::COMPLEX,proc_stride^rank,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < loc_size; ++i)
        {
            int i00 = i & ~stride1 & ~stride2;
            int i01 = i & ~stride1;
            int i10 = (i | stride1) & ~stride2;
            int i11 = i | stride1;
            int ik1 = (i & stride1) >> P1;
            int ik2 = (i & stride2) >> P2;
            int ik = (ik1 << 1) + ik2;

            b[i] = H[ik][0] * a[i00] + H[ik][1] * tmp[i01] +
            H[ik][2] * a[i10] + H[ik][3] * tmp[i11];

        }
    } else if (stride2 < loc_size)
    {
        int proc_stride = stride2 / loc_size;
        vector<complexd> tmp(loc_size);

        MPI_Sendrecv(a.data(),loc_size,MPI::COMPLEX,proc_stride^rank,0,tmp.data(),
        loc_size,MPI::COMPLEX,proc_stride^rank,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < loc_size; ++i)
        {
            int i00 = i & ~stride1 & ~stride2;
            int i01 = i & ~stride1 | stride2;
            int i10 = i & ~stride2;
            int i11 = i | stride2;
            int ik1 = (i & stride1) >> P1;
            int ik2 = (i & stride2) >> P2;
            int ik = (ik1 << 1) + ik2;

            b[i] = H[ik][0] * a[i00] + H[ik][1] * a[i01] +
            H[ik][2] * tmp[i10] + H[ik][3] * tmp[i11];

        }
    } else
    {
        int proc_stride1 = stride1 / loc_size;
        int proc_stride2 = stride2 / loc_size;
        vector<complexd> tmp1(loc_size),tmp2(loc_size),tmp3(loc_size);

        MPI_Sendrecv(a.data(),loc_size,MPI::COMPLEX,proc_stride1^rank,0,tmp1.data(),
        loc_size,MPI::COMPLEX,proc_stride1^rank,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(a.data(),loc_size,MPI::COMPLEX,proc_stride2^rank,0,tmp2.data(),
        loc_size,MPI::COMPLEX,proc_stride2^rank,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(a.data(),loc_size,MPI::COMPLEX,proc_stride1^proc_stride2^rank,0,tmp3.data(),
        loc_size,MPI::COMPLEX,proc_stride1^proc_stride2^rank,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < loc_size; ++i)
        {
            int i00 = i & ~stride1 & ~stride2;
            int i01 = i & ~stride1;
            int i10 = i & ~stride2;
            int i11 = i;
            int ik1 = (i & stride1) >> P1;
            int ik2 = (i & stride2) >> P2;
            int ik = (ik1 << 1) + ik2;
            b[i] = H[ik][0] * a[i00] + H[ik][1] * tmp2[i01] +
            H[ik][2] * tmp1[i10] + H[ik][3] * tmp3[i11];

        }
    }
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

    try
    {
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (argc != 5)
            throw string("Wrong arguments");

        uint N = atoi(argv[1]);
        // uint  = atoi(argv[2]);
        // uint K2 = atoi(argv[3]);
        threads = atoi(argv[4]);

        uint vec_size = (uint) pow((float) 2, (float) N);
        if (vec_size % comm_size != 0 || (uint) comm_size > vec_size)
            throw string("Wrong number of processors");

        uint loc_size = vec_size / comm_size;

        int seed = time(0);
        double start_time, end_time, time_diff_comp = 0;
        double all_times_comp[comm_size];

        vector<complexd> a(loc_size), b(loc_size);

        MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        srand(seed + rank);

        MPI_File fh;

        MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_read_at(fh, rank*loc_size*sizeof(complexd), a.data(), loc_size, MPI_COMPLEX, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);



        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();

        qft(a, b, loc_size, N);


        end_time = MPI_Wtime();
        time_diff_comp = end_time - start_time;


        MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        MPI_File_write_at(fh, rank*loc_size*sizeof(complexd), a.data(), loc_size, MPI_COMPLEX, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);


        MPI_Gather(&time_diff_comp, 1, MPI_DOUBLE, all_times_comp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            cout << comm_size << " " << threads << " " << N << " "  << *std::max_element(&all_times_comp[0], &all_times_comp[comm_size]) << endl;
        }
    }
    catch (const string& e) {
        if (rank == 0)
        cerr << e << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    catch (const MPI::Exception& e) {
        if (rank == 0)
        cerr << "MPI Exception thrown" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    MPI_Finalize();
    return 0;
}
