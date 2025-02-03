

// g++ -O3 -Wall -shared -std=c++11 -I"C:\Python310\include" -I"C:\Python310\Lib\site-packages\pybind11\include" resample.cpp -o resampler.pyd -L"C:\Python310\libs" -lpython310 -static-libgcc -static-libstdc++

// g++ -O3 -Wall -shared -std=c++11 `python3-config --cflags --ldflags` resample.cpp -o resampler.pyd -fPIC

// g++ -O3 -Wall -shared -std=c++11 `python3-config --cflags --ldflags` resample.cpp -o resampler.so -fPIC


#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Function definition
std::vector<int> resampstr(const std::vector<double>& p, int m = -1) {
    int n = p.size();
    if (m == -1) m = n;

    // Normalize p to sum to m
    std::vector<double> pn(n);
    double sum_p = std::accumulate(p.begin(), p.end(), 0.0);
    for (int i = 0; i < n; ++i) {
        pn[i] = p[i] / sum_p * m;
    }

    std::vector<int> s(m, 0);
    std::vector<double> r(n);
    
    // Generate random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        r[i] = dis(gen);
    }

    int k = 0;
    double c = 0.0;

    for (int i = 0; i < n; ++i) {
        c += pn[i];
        if (c >= 1.0) {
            int a = std::floor(c);
            c -= a;
            for (int j = 0; j < a; ++j) {
                s[k] = i + 1; // MATLAB uses 1-based indexing
                ++k;
            }
        }
        if (k < m && c >= r[k]) {
            c -= 1.0;
            s[k] = i + 1;
            ++k;
        }
    }

    if (s[m - 1] == 0) {
        auto mode_val = std::max_element(s.begin(), s.end());
        s[m - 1] = *mode_val;
    }

    if (static_cast<int>(s.size()) > m) {//if (s.size() > m) {
        s.resize(m);
    }

    return s;
}

// Pybind11 binding
PYBIND11_MODULE(resampler, m) {
    m.def("resampstr", &resampstr, "Resample function",
          pybind11::arg("p"), pybind11::arg("m") = -1);
}



