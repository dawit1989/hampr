#include <hampr/core/types.hpp>
#include <hampr/core/config.hpp>
#include <hampr/utils/math_utils.hpp>
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using namespace hampr;

    // Test shift
    {
        array a = {complex(1,0), complex(2,0), complex(3,0), complex(4,0)};
        auto shifted = shift(a, 1);
        assert(shifted[3].real() == 3);
        assert(shifted[0].real() == 0);
    }

    // Test transpose
    {
        mat m = {{complex(1,0), complex(2,0)}, {complex(3,0), complex(4,0)}};
        auto t = transpose(m);
        assert(t[0][1].real() == 3);
        assert(t[1][0].real() == 2);
    }

    // Test conjugate
    {
        array a = {complex(1,2), complex(3,4)};
        auto c = conjugate(a);
        assert(c[0].imag() == -2);
        assert(c[1].imag() == -4);
    }

    // Test dot product
    {
        array a = {complex(1,0), complex(2,0)};
        array b = {complex(3,0), complex(4,0)};
        auto d = dot(a, b);
        assert(d.real() == 11);
    }

    // Test reshape
    {
        array v = {complex(1,0), complex(2,0), complex(3,0), complex(4,0)};
        auto m = reshape(v, 2, 2);
        assert(m[0][0].real() == 1);
        assert(m[1][1].real() == 4);
    }

    // Test convolve
    {
        array a = {complex(1,0), complex(2,0), complex(3,0)};
        array b = {complex(1,0), complex(1,0)};
        auto c = convolve(a, b);
        assert(c.size() == 4);
        assert(c[0].real() == 1);
        assert(c[1].real() == 3);
        assert(c[2].real() == 5);
        assert(c[3].real() == 3);
    }

    // Test config defaults
    {
        Config config;
        assert(config.filter_taps == 128);
        assert(config.num_antennas == 3);
        assert(config.antenna_spacing == 0.528);
    }

    std::cout << "All core tests passed!" << std::endl;
    return 0;
}

