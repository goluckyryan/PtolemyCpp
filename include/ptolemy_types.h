#pragma once

#include <cstdint>
#include <complex>
#include <cstring>

// Fortran type aliases for verbatim transliteration
using real4 = float;
using integer2 = int16_t;
using integer4 = int32_t;
using logical = int;       // Fortran LOGICAL is typically 4 bytes
using logical1 = int8_t;   // LOGICAL*1
using complex16 = std::complex<double>;

// Fortran logical constants
constexpr int TRUE_F  = 1;
constexpr int FALSE_F = 0;

// Helper to convert Fortran logical to C++ bool
inline bool ftobool(int val) { return val != 0; }

// Fortran CHARACTER*8 helper
struct char8 {
    char data[8];
    char8() { std::memset(data, ' ', 8); }
    char8(const char* s) {
        std::memset(data, ' ', 8);
        if (s) {
            // Use strnlen to avoid reading past non-null-terminated 8-byte buffers
            size_t len = strnlen(s, 8);
            std::memcpy(data, s, len);
        }
    }
    bool operator==(const char8& other) const { return std::memcmp(data, other.data, 8) == 0; }
    bool operator!=(const char8& other) const { return std::memcmp(data, other.data, 8) != 0; }
};

