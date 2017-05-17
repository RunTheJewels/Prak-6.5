#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>
#include "mpi.h"
#include <omp.h>


using namespace std;
typedef std::complex<float> complexd;
uint threads;

int rank = 0, comm_size;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

    try
    {
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (argc < 3) throw "Not nuff args";

        int N = atoi(argv[1]);

        uint vec_size = (uint) pow((float) 2, (float) N);
        if (vec_size % comm_size != 0 || (uint) comm_size > vec_size)
            throw string("Wrong number of processors");

        uint loc_size = vec_size / comm_size;

        int seed = time(0);
        vector<complexd> a(loc_size);

        MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        srand(seed + rank);
        
        float sum = 0.0;
        float all_sum = 0.0;
        for (uint i = 0; i < loc_size; ++i)
        {
            a[i] = complexd((float) std::rand()/RAND_MAX, (float) std::rand()/RAND_MAX);
            sum += std::abs(a[i])*std::abs(a[i]);
        }
        
        
        MPI_Allreduce(&sum, &all_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        for (uint i = 0; i < loc_size; ++i)
        {
            a[i] /= all_sum;
        }

        MPI_File fh;

        MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        MPI_File_write_at(fh, rank*loc_size*sizeof(complexd), a.data(), loc_size, MPI_COMPLEX, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
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
