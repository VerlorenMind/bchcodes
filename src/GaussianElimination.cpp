#include "GaussianElimination.h"

int upperEchelonForm(int k, int n, int **newa) {
  // Transforming matrix to an upper-echelon form
  int r = 0;
  int c = 0;
  while(r < k && c < n)
  {
    int index = -1;
    for(unsigned int i=r; i<k; ++i)
    {
      if(newa[i][c])
      {
        index = i;
        break;
      }
    }
    if(index == -1)
    {
      ++c;
    }
    else
    {
      if(index != r)
      {
        for (unsigned int i = 0; i < n; ++i)
        {
          std::swap(newa[index][i], newa[r][i]);
        }
      }
      for(unsigned int i=r+1; i<k; ++i)
      {
        if(newa[i][c])
        {
          for(unsigned int j=c; j<n; ++j)
          {
            newa[i][j] ^= newa[r][j];
          }
        }
      }
      ++r;
      ++c;
    }
  }
  return r;
}
int gaussianElimination(int k, int n, int **a, int ***x) {
  // Memory init
  int **newa = new int*[k];
  for(unsigned int i=0; i<k; ++i)
  {
    newa[i] = new int [n];
    memcpy(newa[i], a[i], sizeof(int)*n);
  }
#ifdef DEBUG
  std::string tempstring = matrix_to_sstream(k, n, newa).str();
    tempstring = matrix_to_sstream(k, n, a).str();
#endif
  // Main algorithm body
  int r = upperEchelonForm(k, n, newa);
  (*x) = new int*[n-r];
  for(unsigned int i=0; i<n-r; ++i)
  {
    (*x)[i] = new int[n];
    memset((*x)[i], 0, sizeof(int)*n);
  }
  int *order = new int[n];
  for(unsigned int i=0; i<n; ++i)
  {
    order[i] = i;
  }
  for(unsigned int i=0; i<n-r; ++i)
  {
    (*x)[i][r+i] = 1;
  }
#ifdef DEBUG
  tempstring = matrix_to_sstream<int>(n-r, n, *x).str();
#endif
  // Permuting columns
  for(unsigned int i=0; i<r; ++i)
  {
    if(!newa[i][i])
    {
      for(unsigned int j=i+1; j<n; ++j)
      {
        if(newa[i][j])
        {
          std::swap(order[i], order[j]);
          for(unsigned int l=0; l<k; ++l)
          {
            std::swap(newa[l][i], newa[l][j]);
#ifdef DEBUG
            tempstring = matrix_to_sstream(k, n, newa).str();
#endif
          }
          break;
        }
      }
    }
  }

#ifdef DEBUG
  tempstring = matrix_to_sstream(k, n, newa).str();
#endif
  // Computing basis of the solution
  for(int i=r-1; i>=0; --i)
  {
    for(unsigned int j=0; j<n-r; ++j)
    {
      for(unsigned int l=i+1; l<n; ++l)
      {
        (*x)[j][i] ^= newa[i][l] & (*x)[j][l];
      }
    }
  }
#ifdef DEBUG
  tempstring = matrix_to_sstream<int>(n-r, n, *x).str();
#endif
  // Returning solution to original ordering of columns
  for(unsigned int i=0; i<n; ++i)
  {
    while(i != order[i])
    {
      for (unsigned int j = 0; j < n-r; ++j)
      {
        std::swap((*x)[j][i], (*x)[j][order[i]]);
      }
      std::swap(order[i], order[order[i]]);

    }
  }
#ifdef DEBUG
  tempstring = matrix_to_sstream<int>(n-r, n, *x).str();
#endif

  // Clean-up
  for(unsigned int i=0; i<k; ++i)
  {
    delete[] newa[i];
  }
  delete[] newa;
  delete[] order;

  return n-r;
}
