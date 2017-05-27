#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <complex>
#include <string.h>
#include <ctime>
#include <stdlib.h>

using namespace std;
typedef std::complex<float> complexd;

const float pi = 3.141592653589793238462643383279502;

complexd operator*(complexd a, complexd b)
{
  return complexd(a.real()*b.real()-a.imag()*b.imag(),a.real()*b.imag()+a.imag()*b.real());
}

complexd operator+(complexd a, complexd b)
{
  return complexd(a.real()+b.real(),b.imag()+a.imag());
}

complexd pow(complexd a, int p)
{
  if (p == 0) return complexd(1,0);
  else if (p == 1) return a;
  else return a*pow(a,p-1);
}

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

    fread(a,sizeof(complexd),vec_size,f);
    fclose(f);

    float sqv = sqrt((float) vec_size);
    complexd w = exp(2 * pi * complexd(0,1) / (float) vec_size);
    // cout << w << " " << sqv << endl;

    for (int i = 0; i < vec_size; i++)
    {
        b[i] = 0;
        // cout << i << endl;
        for (int j = 0; j < vec_size; j++)
        {
            b[i] += a[j] * pow(w, (i*j)%vec_size);
            // cout << a[j] << " " << pow(w, (i*j)%vec_size) << " " << a[j] * pow(w, (i*j)%vec_size) << ", ";
        }
        // cout << endl;
        b[i]/=sqv;
    }

    f = fopen(argv[3],"wb");
    fwrite(b,sizeof(complexd),vec_size,f);
    fclose(f);

    return 0;
}
