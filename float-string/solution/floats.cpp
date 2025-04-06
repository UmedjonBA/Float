
#include "floats.h"
#include <algorithm>
#include <bitset>

Float::Float() {}
Float::Float(int exponent_bits, int mantissa_bits) {
    exponent_size = exponent_bits;
    mantissa_size = mantissa_bits;
    sign = false;
    real_exponent_size = (exponent_bits + 63) / 64;
    real_mantissa_size = (mantissa_bits + 63) / 64;
    exponent.reset(new uint64_t[real_exponent_size]);
    mantissa.reset(new uint64_t[real_mantissa_size]);
}

Float::Float(const Float& other) : Float(other.exponent_size, other.mantissa_size) {
    sign = other.sign;
    std::copy(other.exponent.get(), other.exponent.get() + real_exponent_size, exponent.get());
    std::copy(other.mantissa.get(), other.mantissa.get() + real_mantissa_size, mantissa.get());
}

Float::Float(Float&& other)
    : exponent_size(other.exponent_size),
      mantissa_size(other.mantissa_size),
      sign(other.sign),
      real_exponent_size(other.real_exponent_size),
      real_mantissa_size(other.real_mantissa_size) {
    exponent = std::move(other.exponent);
    mantissa = std::move(other.mantissa);
}

Float& Float::operator=(const Float& other) {
    if (this == &other) {
        return *this;
    }
    operator=(Float(other));
    return *this;
}

Float& Float::operator=(Float&& other) {
    if (this == &other) {
        return *this;
    }
    exponent_size = other.exponent_size;
    mantissa_size = other.mantissa_size;
    sign = other.sign;
    real_exponent_size = other.real_exponent_size;
    real_mantissa_size = other.real_mantissa_size;
    exponent = std::move(other.exponent);
    mantissa = std::move(other.mantissa);
    return *this;
}

void Float::get_real_exp_and_man(exponent_t& exponent, mantissa_t& mantissa) const {
    auto one = mantissa_t::as(mantissa);
    one_shift(one, mantissa_size);
    auto false_zero = exponent_t::as(exponent);
    one_shift(false_zero, exponent_size - 1);
    --false_zero;
    if (!exponent.is_zero()) {
        mantissa += one;
    } else {
        ++exponent;
    }
    exponent -= false_zero;
}

Float& Float::operator+=(const Float& other) {
    operator=(operator+(other));
    return *this;
}
Float& Float::operator-=(const Float& other) {
    operator=(operator-(other));
    return *this;
}
Float& Float::operator*=(const Float& other) {
    operator=(operator*(other));
    return *this;
}
Float& Float::operator/=(const Float& other) {
    operator=(operator/(other));
    return *this;
}

Float Float::operator-() const {
    Float result = *this;
    result.sign = !sign;
    return result;
}

Float Float::sum_helper(const Float& other, bool is_sub) const {
    if (sign != other.sign) {
        return sum_helper(-other, !is_sub);
    }

    Float result(exponent_size, mantissa_size);

    if (is_nan() || other.is_nan()) {
        result.get_nan();
        return result;
    }
    if (is_sub && is_inf() && other.is_inf()) {
        result.get_nan();
        return result;
    }
    if (is_inf()) {
        result.sign = sign;
        result.get_inf();
        return result;
    }
    if (other.is_inf()) {
        result.sign = is_sub ? !other.sign : other.sign;
        result.get_inf();
        return result;
    }

    result.sign = sign;
    auto exponent_a = get_exponent();
    auto exponent_b = other.get_exponent();

    ssize_t diff = (exponent_a - exponent_b).to_int();
    if (abs(diff) > mantissa_size + 1) {
        result = diff > 0 ? *this : other;
        if (is_sub && diff < 0) {
            result.sign = !result.sign;
        }
        return result;
    }

    size_t buffer_size = mantissa_size;
    buffer_size +=
        std::max(exponent_a, exponent_b).to_int() - std::min(exponent_a, exponent_b).to_int();
    buffer_size = (buffer_size + 63) / 64;
    buffer_size += 2;

    auto mantissa_a = get_mantissa(buffer_size);
    auto mantissa_b = other.get_mantissa(buffer_size);

    get_real_exp_and_man(exponent_a, mantissa_a);
    other.get_real_exp_and_man(exponent_b, mantissa_b);

    if (exponent_a < exponent_b || (exponent_a == exponent_b && mantissa_a < mantissa_b)) {
        if (is_sub) {
            result.sign = !sign;
        }
        std::swap(exponent_a, exponent_b);
        std::swap(mantissa_a, mantissa_b);
    }

    while (exponent_a > exponent_b) {
        mantissa_a <<= 1;
        --exponent_a;
    }
    auto exponent = exponent_t::as(exponent_a);
    exponent = exponent_a;
    auto mantissa = mantissa_t::as(mantissa_a);
    mantissa = mantissa_a;
    if (is_sub) {
        mantissa -= mantissa_b;
    } else {
        mantissa += mantissa_b;
    }
    result.normalize(mantissa, exponent);
    return result;
}

Float Float::operator+(const Float& other) const {
    return sum_helper(other, false);
}

Float Float::operator-(const Float& other) const {
    return sum_helper(other, true);
}

Float Float::operator*(const Float& other) const {
    Float result(exponent_size, mantissa_size);
    if (is_nan() || other.is_nan()) {
        result.get_nan();
        return result;
    }
    if (is_inf() && other.is_zero()) {
        result.get_nan();
        return result;
    }
    if (is_zero() && other.is_inf()) {
        result.get_nan();
        return result;
    }

    auto exponent_a = get_exponent();
    auto exponent_b = other.get_exponent();

    auto mantissa_a = get_mantissa(real_mantissa_size * 2 + 1);
    auto mantissa_b = other.get_mantissa(real_mantissa_size * 2 + 1);

    get_real_exp_and_man(exponent_a, mantissa_a);
    other.get_real_exp_and_man(exponent_b, mantissa_b);

    auto exponent = exponent_t::as(exponent_a);
    exponent = exponent_a;
    exponent += exponent_b;
    exponent -= (exponent_t::as(exponent_a) = mantissa_size);
    auto mantissa = mantissa_t::as(mantissa_a);
    mantissa = mantissa_a;
    mantissa *= mantissa_b;
    result.sign = sign ^ other.sign;
    result.normalize(mantissa, exponent);
    return result;
}

Float Float::operator/(const Float& other) const {
    Float result(exponent_size, mantissa_size);
    if (is_nan() || other.is_nan()) {
        result.get_nan();
        return result;
    }

    if (is_inf() && other.is_zero()) {
        result.get_nan();
        return result;
    }

    if (is_zero() && other.is_zero()) {
        result.get_nan();
        return result;
    }

    if (other.is_zero()) {
        result.get_inf();
        result.sign = sign ^ other.sign;
        return result;
    }

    size_t buffer_size = mantissa_size;
    buffer_size *= 3;
    buffer_size = (buffer_size + 63) / 64;
    buffer_size += 2;

    auto exponent_a = get_exponent();
    auto exponent_b = other.get_exponent();

    auto mantissa_a = get_mantissa(buffer_size);
    auto mantissa_b = other.get_mantissa(buffer_size);

    if (is_zero()) {
        result.sign = sign ^ other.sign;
        result.dump_exponent(exponent_a);
        result.dump_mantissa(mantissa_a);
        return result;
    }

    get_real_exp_and_man(exponent_a, mantissa_a);
    other.get_real_exp_and_man(exponent_b, mantissa_b);

    for (size_t i = 0; i < 2 * mantissa_size; ++i) {
        --exponent_a;
        mantissa_a <<= 1;
    }

    auto exponent = exponent_t::as(exponent_a);
    exponent = exponent_a;
    exponent -= exponent_b;
    exponent += (exponent_t::as(exponent_a) = mantissa_size);
    auto mantissa = mantissa_t::as(mantissa_a);
    mantissa = mantissa_a;
    mantissa /= mantissa_b;
    result.sign = sign ^ other.sign;
    result.normalize(mantissa, exponent);
    return result;
}

Float& Float::prev() {
    if (is_nan()) {
        return *this;
    }
    if (is_zero() && !sign) {
        sign = true;
        return *this;
    }
    if (is_inf() && !sign) {
        get_almost_inf();
        return *this;
    }

    auto mantissa = get_mantissa(real_mantissa_size + 1);
    auto exponent = get_exponent();
    get_real_exp_and_man(exponent, mantissa);

    if (sign) {
        ++mantissa;
    } else {
        --mantissa;
    }
    normalize(mantissa, exponent);
    return *this;
}

Float& Float::next() {
    if (is_nan()) {
        return *this;
    }
    if (is_zero() && sign) {
        sign = false;
        return *this;
    }
    if (is_inf() && sign) {
        get_almost_inf();
        return *this;
    }

    auto mantissa = get_mantissa(real_mantissa_size + 1);
    auto exponent = get_exponent();
    get_real_exp_and_man(exponent, mantissa);

    if (sign) {
        --mantissa;
    } else {
        ++mantissa;
    }
    normalize(mantissa, exponent);
    return *this;
}

Float& Float::parse(std::string string) {
    DynNum string_num(string.data());

    std::string temp_string = string;
    std::transform(temp_string.begin(), temp_string.end(), temp_string.begin(),
                   [](char c) { return std::tolower(c); });

    if (temp_string == "inf") {
        get_inf();
        return *this;
    }

    if (temp_string == "-inf") {
        get_inf();
        sign = true;
        return *this;
    }

    if (temp_string == "nan") {
        get_nan();
        return *this;
    }

    get_almost_inf();

    sign = string_num.sign;
    std::string test_str = to_string();
    DynNum test(test_str.data());

    if (sign ? test > string_num : test < string_num) {
        get_inf();
        return *this;
    }

    size_t mantissa_per = 2;
    size_t exponent_per = 2;
    size_t buff_size = (mantissa_size + mantissa_per + 63) / 64 + 1;
    auto mantissa_max = get_mantissa(buff_size);
    auto mantissa_min = get_mantissa(buff_size);
    float_init(this, exponent_size + exponent_per, mantissa_size + mantissa_per);
    sign = string_num.sign;

    auto exponent_max = get_exponent();
    auto exponent_min = get_exponent();

    mantissa_min = 0;
    one_shift(exponent_max, exponent_size);
    exponent_min = 0;

    dump_mantissa(mantissa_min);
    dump_exponent(exponent_min);
    std::string rand_string = to_string();
    DynNum min_dyn(rand_string.data());

    dump_exponent(exponent_max);
    rand_string = to_string();
    DynNum max_dyn(rand_string.data());

    auto false_zero = exponent_t::as(exponent_max);
    one_shift(false_zero, exponent_size - 1);
    --false_zero;

    DynNum half("0.5");

    while (++exponent_min < exponent_max) {
        --exponent_min;
        auto middle = exponent_min + exponent_max;
        middle >>= 1;

        bool is_zero = middle.is_zero();
        DynNum temp = DynNum(middle.is_zero() ? "0" : "1");
        if (is_zero) {
            ++middle;
        }
        middle -= false_zero;
        temp <<= middle.to_int();
        middle += false_zero;
        if (sign) {
            temp = -temp;
        }
        if (is_zero) {
            --middle;
        }

        if (sign ? temp >= string_num : temp <= string_num) {
            exponent_min = middle;
        } else {
            exponent_max = middle;
        }
    }
    --exponent_min;

    one_shift(mantissa_max, mantissa_size);

    dump_exponent(exponent_min);
    dump_mantissa(mantissa_min);
    rand_string = to_string();
    min_dyn = DynNum(rand_string.data());
    dump_mantissa(mantissa_max);
    rand_string = to_string();
    max_dyn = DynNum(rand_string.data());

    while (++mantissa_min < mantissa_max) {

        --mantissa_min;
        auto middle = mantissa_min + mantissa_max;
        middle >>= 1;
        DynNum temp = min_dyn + max_dyn;
        temp *= half;
        if (sign ? temp >= string_num : temp <= string_num) {
            mantissa_min = middle;
            min_dyn = temp;
        } else {
            mantissa_max = middle;
            max_dyn = temp;
        }
    }
    --mantissa_min;

    dump_exponent(exponent_min);
    dump_mantissa(mantissa_min);

    get_real_exp_and_man(exponent_min, mantissa_min);

    exponent_size -= exponent_per;
    mantissa_size -= mantissa_per;
    exponent_min -= mantissa_per;
    normalize(mantissa_min, exponent_min);
    return *this;
}

std::string Float::to_string() const {
    std::string res = "";
    if (is_inf()) {
        res = "Inf";
        if (sign) {
            res = "-" + res;
        }
    } else if (is_nan()) {
        res = "NaN";
    } else if (is_zero()) {
        res = "0";
        if (sign) {
            res = "-" + res;
        }
    } else {
        auto exponent = get_exponent();
        auto mantissa = get_mantissa(real_mantissa_size);
        DynNum temp;
        temp = mantissa;
        temp >>= mantissa_size;
        temp += DynNum(exponent.is_zero() ? "0" : "1");
        if (exponent.is_zero()) {
            ++exponent;
        }
        if (sign) {
            temp = -temp;
        }

        auto false_zero = exponent_t::as(exponent);
        one_shift(false_zero, exponent_size - 1);
        --false_zero;
        exponent -= false_zero;

        temp <<= exponent.to_int();
        res = Conversion::string(temp);
    }
    return res;
}

int float_init(Float* self, int exponent_bits, int mantissa_bits) {
    *self = Float(exponent_bits, mantissa_bits);
    return 0;
}

void float_destroy(Float* self) {}

int Float::float_get_exponent_bits() const {
    return exponent_size;
}

int Float::float_get_mantissa_bits() const {
    return mantissa_size;
}

bool Float::float_get_sign() const {
    return sign;
}

void Float::float_get_exponent(void* target) const {
    size_t count = 0;
    for (size_t i = 0; i < real_exponent_size; ++i) {
        for (size_t j = 0; j < sizeof(uint64_t); ++j) {
            if (count >= exponent_size) {
                break;
            }
            *(static_cast<uint8_t*>(target) + i * sizeof(uint64_t) + j) =
                (exponent.get())[i] >> (j * 8);
            count += 8;
        }
    }
}

void Float::float_get_mantissa(void* target) const {
    size_t count = 0;
    for (size_t i = 0; i < real_mantissa_size; ++i) {
        for (size_t j = 0; j < sizeof(uint64_t); ++j) {
            if (count >= mantissa_size) {
                break;
            }
            *(static_cast<uint8_t*>(target) + i * sizeof(uint64_t) + j) =
                (mantissa.get())[i] >> (j * 8);
            count += 8;
        }
    }
}

void Float::float_set_sign(bool sign) {
    Float::sign = sign;
}

void Float::float_set_exponent(const void* exponent) {
    std::fill(Float::exponent.get(), Float::exponent.get() + (exponent_size + 63) / 64, 0);
    size_t count = 0;
    size_t buffer_size = (exponent_size + 7) / 8;
    for (size_t i = 0; i < real_exponent_size; ++i) {
        Float::exponent.get()[i] = static_cast<const uint64_t*>(exponent)[i];
    }
}

void Float::float_set_mantissa(const void* mantissa) {
    std::fill(Float::mantissa.get(), Float::mantissa.get() + (mantissa_size + 63) / 64, 0);
    size_t count = 0;
    size_t buffer_size = (mantissa_size + 7) / 8;
    for (size_t i = 0; i < real_mantissa_size; ++i) {
        Float::mantissa.get()[i] = static_cast<const uint64_t*>(mantissa)[i];
    }
}

int float_get_exponent_bits(Float* self) {
    return self->float_get_exponent_bits();
}
int float_get_mantissa_bits(Float* self) {
    return self->float_get_mantissa_bits();
}

bool float_get_sign(const Float* self) {
    return self->float_get_sign();
}
void float_get_exponent(const Float* self, void* target) {
    self->float_get_exponent(target);
}
void float_get_mantissa(const Float* self, void* target) {
    self->float_get_mantissa(target);
}

void float_set_sign(Float* self, bool sign) {
    self->float_set_sign(sign);
}
void float_set_exponent(Float* self, const void* exponent) {
    self->float_set_exponent(exponent);
}
void float_set_mantissa(Float* self, const void* mantissa) {
    self->float_set_mantissa(mantissa);
}

Float::mantissa_t Float::get_mantissa(size_t size) const {
    mantissa_t mantissa = mantissa_t::with_digits(size);
    mantissa = 0;
    mantissa_t tmp = mantissa_t::with_digits(size);
    for (ssize_t i = real_mantissa_size - 1; i >= 0; --i) {
        mantissa <<= 64;
        tmp = (Float::mantissa.get())[i];
        mantissa += tmp;
    }
    return mantissa;
}

Float::exponent_t Float::get_exponent() const {
    exponent_t exponent = exponent_t::with_digits(real_exponent_size);
    exponent = 0;
    exponent_t tmp = exponent_t::with_digits(real_exponent_size);
    for (ssize_t i = real_exponent_size - 1; i >= 0; --i) {
        exponent <<= 64;
        tmp = (Float::exponent.get())[i];
        exponent += tmp;
    }
    return exponent;
}

template <bool is_signed, typename T>
size_t Float::get_bit_size(const DynamicLongInt<is_signed, T> long_int) {
    return long_int.get_digit_count() * sizeof(T) * CHAR_BIT - long_int.countl_zero();
}

void Float::dump_mantissa(const mantissa_t& mantissa) {
    std::copy(mantissa.get_digits(), mantissa.get_digits() + real_mantissa_size,
              Float::mantissa.get());
}

void Float::dump_exponent(const exponent_t& exponent) {
    std::copy(exponent.get_digits(), exponent.get_digits() + real_exponent_size,
              Float::exponent.get());
}


template <bool is_signed, typename T>
void Float::one_shift(DynamicLongInt<is_signed, T>& long_int, size_t shift) {
    long_int = 1;
    long_int <<= shift;
}

void Float::normalize(mantissa_t& mantissa, exponent_t& exponent) {
    if (mantissa.is_zero()) {
        exponent = 0;
        dump_exponent(exponent);
        dump_mantissa(mantissa);
        return;
    }
    size_t mantissa_size = get_bit_size(mantissa);
    exponent_t max_exponent = exponent_t::as(exponent);
    one_shift(max_exponent, exponent_size);
    --max_exponent;
    auto false_zero = exponent_t::as(exponent);
    one_shift(false_zero, exponent_size - 1);
    --false_zero;
    exponent += false_zero;
    bool round_up = false;
    bool more_one = false;
    while (mantissa_size > Float::mantissa_size + 1 && exponent < max_exponent) {
        if (round_up) {
            more_one = true;
        }
        round_up = !mantissa.is_even();
        mantissa >>= 1;
        --mantissa_size;
        ++exponent;
    }

    while (mantissa_size < Float::mantissa_size + 1 && !exponent.is_zero() &&
           !exponent.is_negative()) {
        mantissa <<= 1;
        ++mantissa_size;
        --exponent;
    }

    while (exponent.is_negative() && mantissa_size > 0) {
        if (round_up) {
            more_one = true;
        }
        round_up = !mantissa.is_even();
        mantissa >>= 1;
        --mantissa_size;
        ++exponent;
    }

    if (exponent.is_zero()) {
        if (round_up) {
            more_one = true;
        }
        round_up = !mantissa.is_even();
        mantissa >>= 1;
    }

    if (round_up) {
        if (mantissa.is_even()) {
            if (more_one) {
                ++mantissa;
            }
        } else {
            ++mantissa;
        }
        if (get_bit_size(mantissa) > Float::mantissa_size + 1) {
            mantissa >>= 1;
            ++exponent;
        }
    }

    if (!exponent.is_zero()) {
        mantissa_t one = mantissa_t::as(mantissa);
        one_shift(one, Float::mantissa_size);
        mantissa -= one;
    }

    if (exponent.is_negative()) {
        mantissa = 0;
        exponent = 0;
    }

    if (exponent >= max_exponent) {
        mantissa = 0;
        exponent = max_exponent;
    }

    dump_exponent(exponent);
    dump_mantissa(mantissa);
}

bool Float::is_nan() const {
    exponent_t max_exponent = exponent_t::with_digits(real_exponent_size + 1);
    one_shift(max_exponent, exponent_size);
    --max_exponent;
    exponent_t exponent = get_exponent();
    mantissa_t mantissa = get_mantissa(real_mantissa_size);
    if (exponent == max_exponent && !mantissa.is_zero()) {
        return true;
    }
    return false;
}

bool Float::is_inf() const {
    exponent_t max_exponent = exponent_t::with_digits(real_exponent_size + 1);
    one_shift(max_exponent, exponent_size);
    --max_exponent;
    exponent_t exponent = get_exponent();
    mantissa_t mantissa = get_mantissa(real_mantissa_size);
    if (exponent == max_exponent && mantissa.is_zero()) {
        return true;
    }
    return false;
}

bool Float::is_zero() const {
    exponent_t exponent = get_exponent();
    mantissa_t mantissa = get_mantissa(real_mantissa_size);
    if (exponent.is_zero() && mantissa.is_zero()) {
        return true;
    }
    return false;
}

void Float::get_nan() {
    exponent_t max_exponent = exponent_t::with_digits(real_exponent_size + 1);
    one_shift(max_exponent, exponent_size);
    --max_exponent;
    exponent_t exponent = max_exponent;
    mantissa_t mantissa = mantissa_t::with_digits(real_mantissa_size);
    one_shift(mantissa, mantissa_size - 1);
    dump_exponent(exponent);
    dump_mantissa(mantissa);
}

void Float::get_inf() {
    exponent_t max_exponent = exponent_t::with_digits(real_exponent_size + 1);
    one_shift(max_exponent, exponent_size);
    --max_exponent;
    exponent_t exponent = max_exponent;
    mantissa_t mantissa = mantissa_t::with_digits(real_mantissa_size);
    mantissa = 0;
    dump_exponent(exponent);
    dump_mantissa(mantissa);
}

void Float::get_almost_inf() {
    exponent_t max_exponent = exponent_t::with_digits(real_exponent_size + 1);
    one_shift(max_exponent, exponent_size);
    --max_exponent;
    --max_exponent;
    exponent_t exponent = max_exponent;
    mantissa_t mantissa = mantissa_t::with_digits(real_mantissa_size);
    one_shift(mantissa, mantissa_size);
    --mantissa;
    dump_exponent(exponent);
    dump_mantissa(mantissa);
}

void float_add(Float* result, const Float* a, const Float* b) {
    *result = (*a + *b);
}

void float_sub(Float* result, const Float* a, const Float* b) {
    *result = (*a - *b);
}

void float_mul(Float* result, const Float* a, const Float* b) {
    *result = (*a * *b);
}

void float_div(Float* result, const Float* a, const Float* b) {
    *result = (*a / *b);
}

void float_next(Float* self) {
    self->next();
}

void float_prev(Float* self) {
    self->prev();
}

int float_string(const Float* self, char* string, int n) {
    std::string res = self->to_string();
    if (n != 0) {
        if (res.size() > n) {
            std::copy(res.begin(), res.begin() + n, string);
        } else {
            std::copy(res.begin(), res.end(), string);
        }
    }
    return res.size();
}

void float_parse(Float* self, const char* string) {
    self->parse(string);
}
