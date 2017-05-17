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
    if (argv < 2)
    {
        cout << "Filename needed\n";
        return -1;
    }
    FILE* f = 0;
    f = fopen(argv[1]);
    if (f == 0)
    {
        cout << "No such file\n";
        return -1;
    }
    complexd c;
    bool read;
    while (fread(&c,sizeof(comlpexd),1,f) != 0)
        cout << c << " ";
    fclose(f);

    return 0;
}
