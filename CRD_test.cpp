#include <iostream>
#include "headers/ThechycoGlobalVar.hpp"


inline void testCRD() {
    int current_line = 0;
    int next_line = 0;
    int count = 0;
    for (size_t i = 0; i < 2 * mf; ++i) {
        next_line = crd[1][i / 2];
        if (i % 2 == 0 && (next_line - current_line) == 0) {
            std::cout << "current_line = " << current_line << "  " << i / 2 - count << "    " << crd[i % 2][i / 2] - crd[i % 2][i / 2 - 1] << std::endl;
        }
        if ((next_line - current_line) != 0) { std::cout << "/********************/" << std::endl; count = i / 2; }
        current_line = next_line;
    }
}