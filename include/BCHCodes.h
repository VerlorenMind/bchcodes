#ifndef BCHCODES_BCHCODES_H
#define BCHCODES_BCHCODES_H

#include <cstdint>


void build_bch_matrices(uint64_t n, uint64_t d, std::ostream& output);
void build_bch_matrices_bit_order(uint64_t n, uint64_t d, std::ostream& output);

#endif
