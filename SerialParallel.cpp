#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
using namespace std;
template<typename Type>
void SerialProduct(Type* a, Type* b, Type* c, int n) {
    for (int k1 = 0; k1 < n; k1++) {
        for (int k2 = 0; k2 < n; k2++) {
            c[k1] += a[k1 + k2 * n] * b[k2];
        }
    }
}
template<typename Type1>
void ParallelProductRow(Type1* a, Type1* b, Type1* c, int n) {
#pragma omp  parallel for  
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[i] += a[i + j * n] * b[j];
    }
    }
}
template<typename Type2>
void ParallelProductCol(Type2* a, Type2* b, Type2* c, int n) { 
/*#pragma omp parallel 
    {   double c_prv[10000] = {};
        int i, j;
        for (i = 0; i < n; i++) {
            #pragma omp for
            for (j = 0; j < n; j++) {
                c_prv[i] += a[i + j * n] * b[j];
            }
        }
#pragma omp critical
        {
            for (i = 0; i < n; i++) c[i] += c_prv[i];
        }
    }*/
    int sum = 0;
    for (int i = 0; i < n; i++) {
        sum = 0;
#pragma omp parallel for reduction (+:sum)
        for (int j = 0; j < n; j++) {
            sum += a[i + j * n] * b[j];
        }
        c[i] = sum;
}
}
int main()
{
    ofstream fout("Data.csv");
    if (!fout) {
        cout << "Error!";
    }
    int i = 0;
    int n = 10000;
    int* a = new int[n * n];
    int* b = new int[n];
    int* c = new int[n];
    std::string ss = "";
    for (int ni = 100; ni <= 1000; ni += 100) {
        for (int i = 0; i < ni; i++) {
            b[i] = 1;
            c[i] = 0;
            for (int j = 0; j < ni; j++) {
                a[i + j * ni] = 1;
            }
        }
        double t = omp_get_wtime();
        SerialProduct(a, b, c, ni);
        t = omp_get_wtime() - t;
        ss += std::to_string(ni) + ";" + std::to_string(t) + ";";
        printf_s("Time serial= %f\n", t);
       /* printf_s("Matrix a=\n");
        for (int i = 0; i < ni; i++) {
            printf_s("\n");
            for (int j = 0; j < ni; j++) {
                printf_s("%d ", a[i + j * ni]);
            }
        }
        printf_s("\n");
        printf_s("Vector b=\n");
        for (int i = 0; i < ni; i++) {
            printf_s("%d ", b[i]);

        }*/
        printf_s("\n");
        printf_s("Vector c=\n");
        for (int i = 0; i < ni; i++) {
            printf_s("%d ", c[i]);
            c[i] = 0;
        }
        t = omp_get_wtime();
        ParallelProductRow(a, b, c, ni);
        t = omp_get_wtime() - t;
        ss += std::to_string(t) + ";";
        printf_s("\n");
        printf_s("Time parallelRow= %f\n", t);
        printf_s("Vector c=\n");
        for (int i = 0; i < ni; i++) {
            printf_s("%d ", c[i]);
            c[i] = 0;
        }
        t = omp_get_wtime();
        ParallelProductCol(a, b, c, ni);
        t = omp_get_wtime() - t;
        ss += std::to_string(t) + ";";
        fout << ss;
        fout << "\n";
        ss = "";
        printf_s("\n");
        printf_s("Time parallelCol= %f\n", t);
        printf_s("Vector c=\n");
        for (int i = 0; i < ni; i++) {
            printf_s("%d ", c[i]);
            c[i] = 0;
        }  
    }
    fout.close();
}

