#ifndef MATRIX_H
#define MATRIX_H
#include <vector>

void MatVect(int N, const std::vector<double>& A, const std::vector<int>& NC, const std::vector<int>& NE, const std::vector<double>& X, std::vector<double>& B);
void MatTVect(int N, const std::vector<double>& A, const std::vector<int>& NC, const std::vector<int>& NE, const std::vector<double>& X, std::vector<double>& B);
double VectNorm2(int N, const std::vector<double>& A);
void ResCalc(int N, const std::vector<double>& A, const std::vector<int>& NC, const std::vector<int>& NE, const std::vector<double>& B, std::vector<double>& x0, std::vector<double>& x1, double eps);

#endif