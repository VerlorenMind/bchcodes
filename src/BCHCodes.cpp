#include <iostream>
#include "GF.h"
#include "../include/BCHCodes.h"

void build_bch_matrices(uint64_t n, uint64_t d, std::ostream& output)
{
    if(__builtin_popcount(n+1) > 1)
    {
        std::cerr << "Length is not a 2^m-1"<<std::endl;
    }
    uint64_t pow = 0;
    uint64_t tmp = n;
    tmp += 1;
    while(tmp > 0)
    {
        ++pow;
        tmp >>= 1;
    }
    --pow;
    GF2 field(pow);
    uint64_t b = 2;
    uint64_t bpow = b;
    uint64_t gen = field.min_poly(b);
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
    uint64_t hgen, res;
    div_poly((uint64_t(1)<<n)+1, gen, hgen, res);
    output<<"Built a ("<<n<<", "<<k<<") coder with gen poly\n"
             <<std::bitset<64>(gen) <<"\nand check poly\n"
             <<std::bitset<64>(hgen)<<std::endl;
    uint64_t* g = new uint64_t[k];
    uint64_t* h = new uint64_t[n-k];
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
            h[i] ^= ((hgen & (uint64_t(1) << hdeg-j)) >> (hdeg-j)) << (j + i);
        }
    }
    output<<"Gen matrix:\n";
    for(uint64_t i=0; i<k; ++i)
    {
        for(uint64_t j=0; j<n; ++j)
            output<<(g[i] & (uint64_t(1)<<j) ? 1 : 0)<<" ";
        output<<"\n";
    }
    output<<"Check matrix:\n";
    for(uint64_t i=0; i<n-k; ++i)
    {
        for(uint64_t j=0; j<n; ++j)
            output<<(h[i] & (uint64_t(1)<<j) ? 1 : 0)<<" ";
        output<<"\n";
    }
    delete[] g;
    delete[] h;
}