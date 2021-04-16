#include <cstdio>
#include <cmath>
#include <complex>

typedef std::complex<double> dcomplex;

int main()
{
    const int n = 256;
    dcomplex* e_r = new dcomplex[(n+1)*(n+1)];
    FILE* file = fopen("mri_data_256.dat", "rb");
    fread(e_r, sizeof(dcomplex), (n+1)*(n+1), file);
    
    for (int i=0; i<n+1; ++i) {
        for (int j=0; j<n+1; ++j) {
            printf("%g ", std::abs(e_r[(n+1)*i+j]));
        }
        puts("");
    }

    fclose(file);
    delete [] e_r;
}
