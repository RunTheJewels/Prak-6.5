#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;
typedef std::complex<float> complexd;

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
    
    f = fopen(argv[3],"rb");
    if (f == 0)
    {
        cout << "No such file\n";
        return -1;
    }

    for (int i = 0; i < vec_size; i++)
    {
        fread((b+i),sizeof(complexd),1,f);
    }
    fclose(f);
    
    float diff = 0;
    for (int i = 0; i < vec_size; i++)
    {
        diff += abs(a[i]-b[i]);
    }
    
    cout << diff << endl;
    
    return 0;
}