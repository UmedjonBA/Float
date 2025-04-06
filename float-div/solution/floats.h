#pragma once

// #include <cstdbool.h>
#include <cstdint>
#include <memory>
#include "Conversion.hpp"
#include "DynNum.hpp"
#include "LongInt.hpp"

using std::size_t;

class Float {
  public:
    using exponent_t = DynamicLongInt<true, uint64_t>;
    using mantissa_t = DynamicLongInt<false, uint64_t>;

  private:
    size_t exponent_size;
    size_t mantissa_size;
    bool sign;
    std::unique_ptr<uint64_t, std::default_delete<uint64_t[]>> exponent;
    std::unique_ptr<uint64_t, std::default_delete<uint64_t[]>> mantissa;

    size_t real_exponent_size;
    size_t real_mantissa_size;

    exponent_t get_exponent() const;
    mantissa_t get_mantissa(size_t size) const;

    template <bool is_signed, typename T>
    static void one_shift(DynamicLongInt<is_signed, T>& long_int, size_t shift);

    void get_real_exp_and_man(exponent_t& exponent, mantissa_t& mantissa) const;

    void get_nan();
    void get_inf();
    void get_almost_inf();

    void dump_exponent(const exponent_t& exponent);
    void dump_mantissa(const mantissa_t& mantissa);


    template <bool is_signed, typename T>
    static size_t get_bit_size(const DynamicLongInt<is_signed, T> long_int);

    void normalize(mantissa_t& mantissa, exponent_t& exponent);

    Float sum_helper(const Float& other, bool is_sub) const;

  public:
    Float();
    Float(int exponent_bits, int mantissa_bits);

    Float(const Float& other);
    Float(Float&& other);

    bool is_inf() const;
    bool is_nan() const;
    bool is_zero() const;

    Float& operator=(const Float& other);
    Float& operator=(Float&& other);

    Float& operator+=(const Float& other);
    Float& operator-=(const Float& other);
    Float& operator*=(const Float& other);
    Float& operator/=(const Float& other);
    Float operator+(const Float& other) const;
    Float operator-(const Float& other) const;
    Float operator*(const Float& other) const;
    Float operator/(const Float& other) const;
    Float operator-() const;

    int float_get_exponent_bits() const;
    int float_get_mantissa_bits() const;

    bool float_get_sign() const;
    void float_get_exponent(void* target) const;
    void float_get_mantissa(void* target) const;

    void float_set_sign(bool sign);
    void float_set_exponent(const void* exponent);
    void float_set_mantissa(const void* mantissa);

    Float& parse(std::string string);
    std::string to_string() const;

    Float& prev();
    Float& next();
};

int float_init(Float* self, int exponent_bits, int mantissa_bits);
void float_destroy(Float* self);

int float_get_exponent_bits(Float* self);
int float_get_mantissa_bits(Float* self);

bool float_get_sign(const Float* self);
void float_get_exponent(const Float* self, void* target);
void float_get_mantissa(const Float* self, void* target);

void float_set_sign(Float* self, bool sign);
void float_set_exponent(Float* self, const void* exponent);
void float_set_mantissa(Float* self, const void* mantissa);

DynamicLongInt<true, uint64_t> get_exponent(const Float*);
DynamicLongInt<false, uint64_t> get_mantissa(const Float*, size_t);

void float_add(Float* result, const Float* a, const Float* b);
void float_sub(Float* result, const Float* a, const Float* b);
void float_mul(Float* result, const Float* a, const Float* b);
void float_div(Float* result, const Float* a, const Float* b);

void float_next(Float* self);
void float_prev(Float* self);

void float_parse(Float* self, const char* string);
int float_string(const Float* self, char* string, int n);

