#include <cstdio>
#include <fstream>
#include <complex>
int main(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        std::ifstream f(argv[i], std::ios::binary);
        char m[4]; f.read(m, 4); int r, c; f.read((char*)&r, 4); f.read((char*)&c, 4);
        std::complex<double> v[7]; f.read((char*)v, 7 * sizeof(std::complex<double>));
        printf("%s: range_ref=%.0f range=%.0f doppler_ref=%.0f doppler=%.0f az_ref=%.6f az=%.6f sinr=%.15g\n",
            argv[i], v[0].real(), v[1].real(), v[2].real(), v[3].real(), v[4].real(), v[5].real(), v[6].real());
    }
}