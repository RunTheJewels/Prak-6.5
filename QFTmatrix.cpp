#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;
typedef std::complex<float> complexd;

const float pi = 3.14159265358;

int main(int argc, char** argv)
{
    if (argc < 4) 
    {
        cout << "Error";
        return -1;
    }
    
    uint N = atoi(argv[1]);
    uint vec_size = (uint) pow((float) 2, (float) N);
    
    complexd a[vec_size],b[vec_size];
    
    FILE* f = 0;
    f = fopen(argv[2],"rb");
    if (f == 0)
    {
        cout << "No such file\n";
        return -1;
    }

    for (int i = 0; i < vec_size; i++)
    {
        fread((a+i),sizeof(complexd),1,f);
    }
    fclose(f);
        
    complexd w = exp(2 * pi * complexd(0,1) / (float) vec_size);
    
    for (int i = 0; i < vec_size; i++)
    {
        b[i] = 0;
        for (int j = 0; j < vec_size; j++)
        {
            b[i] += a[j] * pow(w, i+j) / sqrt((float) vec_size);
        }
    }
    
    f = fopen(argv[3],"wb");
    for (int i = 0; i < vec_size; i++)
    {
        fwrite((b+i),sizeof(complexd),1,f);
    }
    fclose(f);
    
    return 0;
}