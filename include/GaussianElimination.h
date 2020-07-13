#ifndef BCHCODES_INCLUDE_GAUSSIANELIMINATION_H_
#define BCHCODES_INCLUDE_GAUSSIANELIMINATION_H_

#include <algorithm>
#include <cstring>
int upperEchelonForm(int k, int n, int **newa);
int gaussianElimination(int k, int n, int **a, int ***x);

#endif //BCHCODES_INCLUDE_GAUSSIANELIMINATION_H_
