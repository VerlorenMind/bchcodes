#include <iostream>
#include <cstring>
#include <GaussianElimination.h>
#include "GF.h"
#include "../include/BCHCodes.h"

void build_bch_matrices(uint64_t n, uint64_t d, std::ostream& output)
{
    if(__builtin_popcount(n+1) > 1)
    {
        std::cerr << "Length is not a 2^m-1"<<std::endl;
    }
    uint64_t pow = 0;
    uint128 tmp = n;
    tmp += 1;
    while(tmp > 0)
    {
        ++pow;
        tmp >>= 1;
    }
    --pow;
    GF2 field(pow);
    uint128 b = 2;
    uint128 bpow = b;
    uint128 gen = field.min_poly(b);
    for(uint64_t i=2; i<d; ++i)
    {
        bpow = field.mul(bpow, b);
        gen = lcm_poly(gen, field.min_poly(bpow));
    }
    tmp = gen;
    uint64_t k = 0;
    while(tmp > 0)
    {
        ++k;
        tmp >>= 1;
    }
    --k;
    k = n - k;
    uint128 hgen, res;
    div_poly((uint128(1)<<n)+1, gen, hgen, res);
    output<<"Built a ("<<n<<", "<<k<<") coder with gen poly\n"
             <<std::bitset<128>(gen) <<"\nand check poly\n"
             <<std::bitset<128>(hgen)<<std::endl;
    uint128* g = new uint128[k];
    uint128* h = new uint128[n-k];
    uint64_t hdeg = 0;
    tmp = hgen;
    while(tmp > 0)
    {
        ++hdeg;
        tmp >>= 1;
    }
    --hdeg;

    for(uint64_t i=0; i<k; ++i)
    {
        g[i] = gen << i;
    }
    for(uint64_t i=0; i<n-k; ++i)
    {
        h[i] = 0;
        for(uint64_t j=0; j<=hdeg; ++j)
        {
            h[i] ^= ((hgen & (uint128(1) << (hdeg-j))) >> (hdeg-j)) << (j + i);
        }
    }
    output<<"Gen matrix:\n";
    for(uint64_t i=0; i<k; ++i)
    {
        for(uint64_t j=0; j<n; ++j)
            output<<(g[i] & (uint128(1)<<j) ? 1 : 0)<<" ";
        output<<"\n";
    }
    output<<"Check matrix:\n";
    for(uint64_t i=0; i<n-k; ++i)
    {
        for(uint64_t j=0; j<n; ++j)
            output<<(h[i] & (uint128(1)<<j) ? 1 : 0)<<" ";
        output<<"\n";
    }
    delete[] g;
    delete[] h;
}

void build_bch_matrices_bit_order(uint64_t n, uint64_t d, std::ostream& output)
{
    if(__builtin_popcount(n+1) > 1)
    {
        std::cerr << "Length is not a 2^m-1"<<std::endl;
    }
    uint64_t pow = 0;
    uint128 tmp = n;
    tmp += 1;
    while(tmp > 0)
    {
        ++pow;
        tmp >>= 1;
    }
    --pow;
    GF2 field(pow);
    uint128 b = 2;
    uint128 bpow = b;
    uint128 gen = field.min_poly(b);
    int *powers = new int[d];
    powers[0] = 0;
    powers[1] = 1;
    int indPowers = 2;
    for(uint64_t i=2; i<d; ++i)
    {
        bpow = field.mul(bpow, b);
        uint128 tmp = lcm_poly(gen, field.min_poly(bpow));
        if(tmp != gen)
        {
            powers[indPowers++] = i;
        }
        gen = tmp;
    }
    tmp = gen;
    uint64_t k = 0;
    while(tmp > 0)
    {
        ++k;
        tmp >>= 1;
    }
    --k;
    k = n - k;
    uint128 hgen, res;
    div_poly((uint128(1)<<n)+1, gen, hgen, res);
    // output<<"Built a ("<<n<<", "<<k<<") coder with gen poly\n"
    //       <<std::bitset<128>(gen) <<"\nand check poly\n"
    //      <<std::bitset<128>(hgen)<<std::endl;
    uint128* h = new uint128[indPowers*pow];
    memset(h, 0, sizeof(uint128)*(indPowers*pow));
    uint64_t elem = 0;
    uint64_t temp;
    ++n;
    for(unsigned int i=0; i<n; ++i)
    {
        for(unsigned int j=0; j<indPowers; ++j)
        {
            temp = field.pow(elem, powers[j]);
            for(unsigned int l=0; l<pow; ++l)
            {
                h[j*pow+l] ^= (temp & (uint128(1) << l)) ? uint128(1) << i : 0;
            }
        }
        ++elem;
    }
    int **h_buf = new int*[indPowers*pow];
    for(uint64_t i=0; i < indPowers*pow; ++i) {
      h_buf[i] = new int[n];
      for (unsigned int j = 0; j < n; ++j) {
        h_buf[i][j] = (h[i] & (uint128(1) << j) ? 1 : 0);
      }
    }
    int ***g = new int**;
    upperEchelonForm(indPowers*pow, n, h_buf);
    gaussianElimination(n-k, n, h_buf, g);
    output<<"EBCH("<<n<<", "<<k<<")\n\n";
    output<<n<<" "<<k<<"\n\n";
    for(uint64_t i=0; i<k; ++i)
    {
      for(uint64_t j=0; j<n; ++j)
        output<<(*g)[i][j]<<" ";
      output<<"\n";
    }
    output<<"\n";
    for(uint64_t i=0; i<n-k; ++i)
    {
        for(uint64_t j=0; j<n; ++j)
            output<<h_buf[i][j]<<" ";
        output<<"\n";
    }
    delete[] h;
    for(unsigned int i=0; i<indPowers*pow; ++i) {
      delete[] h_buf[i];
    }
    delete[] h_buf;
    for(unsigned int i=0; i<k; ++i) {
      delete[] (*g)[i];
    }
    delete[] (*g);
    delete g;
    delete[] powers;
}
