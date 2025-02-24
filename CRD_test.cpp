#include <iostream>
#include "headers/ThechycoGlobalVar.hpp"


inline void testCRD() {
    int current_line = 0;
    int next_line = 0;
    int first_elem_line = 0;
    int number_rod = 0;

    for (size_t i = 0; i < 2 * mf; ++i) {
        next_line = crd[1][i / 2];
        if (i % 2 == 0 && (next_line - current_line) == 0) {
            std::cout << "current line = " << current_line << "  " 
                      << "num of elem = " << i / 2 - first_elem_line << "    "
                                          << crd[i % 2][i / 2] - crd[i % 2][i / 2 - 1] << std::endl;

            if (crd[i % 2][i / 2] - crd[i % 2][i / 2 - 1] == 4) {
                ++number_rod;
            }
        }
        if ((next_line - current_line) != 0) { 
            std::cout << "\n" << "number of control rod = " << number_rod << std::endl;
            number_rod = 0;
            std::cout << "/********************/" << std::endl; 
            first_elem_line = i / 2; 
        }
        current_line = next_line;
    }
}