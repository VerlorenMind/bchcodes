#ifndef BCH_GF_H
#define BCH_GF_H

#include <iostream>
#include <bitset>
#include <vector>
#include <string>

const uint64_t PRIMITIVE_POLYS[17] =
        {
                0,
                0,
                0b111,
                0b1011,
                0b10011,
                0b100101,
                0b1000011,
                0b10001001,
                0b100011101,
                0b1000010001,
                0b10000001001,
                0b100000000101,
                0b1000001010011,
                0b10000000011011,
                0b100010001000011,
                0b1000000000000011,
                0b10001000000001011
        };


class GF2 {
/*
 *   Class for elements from GF(2^m)
 */
private:
    uint64_t deg; // Power of a primitive polynomial
    uint64_t primitive; // Primitive polynomial
    uint64_t num; // Number of non-zero elements
    std::pair<uint64_t, uint64_t>* logs; // Table of logs as <element : generator power>
public:
    GF2(uint64_t pow);
    GF2();
    GF2(const GF2& x);
    ~GF2();
    uint64_t add(const uint64_t& a, const uint64_t& b) const;
    uint64_t mul(const uint64_t& a, const uint64_t& b) const;
    uint64_t inv(const uint64_t& a) const;
    uint64_t min_poly(const uint64_t& a) const;
    uint64_t get_elem_deg(uint64_t elem) const;
    uint64_t pow(const uint64_t& a, int pow) const;
    bool is_in(const uint64_t& x) const;
    uint64_t get_deg() const;
    uint64_t get_num() const;
    GF2& operator=(const GF2& a);
};

class GF2X {
    /*
     * Class for polynomials with coefficients from GF(2^m)
     */
private:
    GF2 F; // Field corresponding to coefficients
    uint64_t* c; // Coefficients
    uint64_t deg; // Polynomial degree
    void cut(); // Lowering the degree if highest coeff is zero
public:
    GF2X();
    GF2X(const GF2& field, uint64_t* coeffs, const uint64_t& pow);
    GF2X(const GF2& field, const std::vector<uint64_t>& coeffs);
    GF2X(const GF2& field, uint64_t poly); // Converting element from GF(2^m)
    GF2X(const GF2X& x);
    ~GF2X();
    GF2X operator+(const GF2X& a) const;
    GF2X operator*(const GF2X& b) const;
    GF2X operator*(uint64_t x) const;
    uint64_t& operator[](uint64_t i);
    uint64_t eval(const uint64_t& x) const;
    GF2X& operator=(const GF2X& a);
    uint64_t get_deg();
    std::string print();
    uint64_t* roots(uint64_t& num);

};

uint64_t mul_poly(uint64_t a, uint64_t b);
void div_poly(const uint64_t& a, const uint64_t& b, uint64_t &q, uint64_t &r);
uint64_t gcd_poly(uint64_t a, uint64_t b);
uint64_t lcm_poly(uint64_t a, uint64_t b);
uint64_t mod(uint64_t a, uint64_t b);

#endif //BCH_GF_H
