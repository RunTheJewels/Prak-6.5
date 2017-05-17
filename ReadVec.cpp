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
    if (argc < 2)
    {
        cout << "Filename needed\n";
        return -1;
    }
    FILE* f = 0;
    f = fopen(argv[1],"rb");
    if (f == 0)
    {
        cout << "No such file\n";
        return -1;
    }
    complexd c;
    while (fread(&c,sizeof(complexd),1,f) != 0)
        cout << c << " ";
    fclose(f);

    return 0;
}
