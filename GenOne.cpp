#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>
#include <stdlib.h>

using namespace std;
typedef std::complex<float> complexd;

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        cout << "Filename needed\n";
        return -1;
    }
    FILE* f = 0;
    f = fopen(argv[1],"wb");
    if (f == 0)
    {
        cout << "Writing error\n";
        return -1;
    }
    complexd zero = 0, one = 1;
    uint vec_size = (uint) pow((float) 2, (float) atoi(argv[2]));
    uint elem_one = atoi(argv[3]);
    for (uint i = 0; i < vec_size; i++)
    {
      fwrite(&zero,sizeof(complexd),1,f);
    }
    fseek(f,elem_one*sizeof(complexd),SEEK_SET);
    fwrite(&one,sizeof(complexd),1,f);
    fclose(f);

    return 0;
}
