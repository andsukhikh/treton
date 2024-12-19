#include <vector>
#include <cmath>
#include <iostream>

//!!!!!!!!!!!!!!здесь возможна ошибка в логике 
//! 
//! 
// 		Расчет произведения матрицы А[NxN] на вектор х[N]. Результат в векторе В[N]. 
// 		Структура матрицы задана с помощью массивов связей NC и NE.
// 		Хранение только отличных от нуля членов. Сначала в матрице А расположены члены с главной диагонали

void MatVect(int N, const std::vector<double>& A, const std::vector<int>& NC, 
                const std::vector<int>& NE, const std::vector<double>& X, std::vector<double>& B) {
    for (int i = 0; i < N; ++i) {
        double s = A[i] * X[i]; 
        for (int j = NC[i]; j < NC[i + 1]; ++j) {
            //std::cout << "NE[" << j << "] = " << NE[j] << std::endl;
            s += A[j + N] * X[NE[j]]; 
        }
        B[i] = s;
    }
}


// 		Расчет произведения матрицы АT[NxN] (транспонированной матрицы по отношению к матрице А) на вектор х[N].  
// 		Результат в векторе в[N]. Структура матрицы задана с помощью массивов связей NC и NE. 
// 		Хранение только отличных от нуля членов. Сначала в матрице А расположены члены с главной диагонали


void MatTVect(int N, const std::vector<double>& A, const std::vector<int>& NC, 
                const std::vector<int>& NE, const std::vector<double>& X, std::vector<double>& B) {
    for (int i = 0; i < N; ++i) {
        double s = A[i] * X[i];
        for (int j = NC[i]; j < NC[i + 1]; ++j) {
            int k = NE[j];
            for (int l = NC[k]; l < NC[k + 1]; ++l) {
                if (NE[l] == i) {
                    s += A[N + l] * X[k];
                    break;
                }
            }
        }
        B[i] = s;
    }
}

//		Расчет квадрата нормы вектора А

double VectNorm2(int N, const std::vector<double>& A) {
    double s = 0.0;
    for (int i = 0; i < N; ++i) {
        s += A[i] * A[i];
    }
    return s;
}


// 		Решение системы уравнений Ах = в методои сопряженных градиентов. Система преобразуется к виду (АТА)х = АТв.
// 		Структура матрицы задана с помощью массивов связей NC и NE. Заданная погрешность расчета - eps, начальное
// 		приближение х0, результат - Х1.

void ResCalc(int N, const std::vector<double>& A, const std::vector<int>& NC, const std::vector<int>& NE, const std::vector<double>& B, 
             std::vector<double>& x0, std::vector<double>& x1, double eps) {
    std::vector<double> R(N), P(N), Ar(N), Ap(N);
    int ni = 10000;
    int k = 0;

    MatVect(N, A, NC, NE, x0, Ar);
    for (int i = 0; i < N; ++i) {
        R[i] = B[i] - Ar[i];
    }

    MatTVect(N, A, NC, NE, R, P);
    /*for (int i = 0; i < N; ++i) {
        Ar[i] = P[i];
    }*/

    Ar = P;

    MatVect(N, A, NC, NE, P, Ap);

    while (k < ni) {
        ++k;
        double ArNorm2 = VectNorm2(N, Ar);
        double alf = ArNorm2 / VectNorm2(N, Ap);
        double amax = 0.0;

        for (int i = 0; i < N; ++i) {
            double dx = alf * P[i];
            x1[i] = x0[i] + dx;
            double adx = std::abs(dx);
            if (adx > amax) amax = adx;
        }

        if (amax < eps) {
            break;
        }

        for (int i = 0; i < N; ++i) {
            R[i] = R[i] - alf * Ap[i];
        }

        MatTVect(N, A, NC, NE, R, Ar);
        double betta = VectNorm2(N, Ar) / ArNorm2;

        for (int i = 0; i < N; ++i) {
            P[i] = Ar[i] + betta * P[i];
            x0[i] = x1[i];
        }

        MatVect(N, A, NC, NE, P, Ap);
    }

    if (k == ni) {
        std::cout << "Превышено максимальное число итераций ni = " << ni << std::endl;
    }
}

