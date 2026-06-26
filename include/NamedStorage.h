#pragma once
// NamedStorage — named dynamic-array storage.
//
// Stores named std::vector<double> arrays keyed by an 8-char-style name.
// Used for:
//   - DEFINE input cards (user-supplied named arrays)
//   - linkule SHAPE potential-shape buffers
//   - β-deformation arrays read back by the SHAPE linkule
//
// Names are normalized to exactly 8 chars (truncated or space-padded) to
// match the Fortran CHARACTER*8 convention the legacy pool keyed on.
// Single-threaded, no cross-call lifetime constraint, no fragmentation.

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

class NamedStorage {
public:
    // Lookup. Returns nullptr if not found.
    const std::vector<double>* find(std::string_view name) const {
        auto it = store_.find(normalize(name));
        return (it == store_.end()) ? nullptr : &it->second;
    }
    std::vector<double>* find(std::string_view name) {
        auto it = store_.find(normalize(name));
        return (it == store_.end()) ? nullptr : &it->second;
    }

    // Create or replace a named array of size n (zero-initialized).
    std::vector<double>& alloc(std::string_view name, std::size_t n) {
        auto& v = store_[normalize(name)];
        v.assign(n, 0.0);
        return v;
    }

    // Store/copy user-input values under name (replaces any prior entry).
    std::vector<double>& define(std::string_view name, std::vector<double> values) {
        auto& v = store_[normalize(name)];
        v = std::move(values);
        return v;
    }

    bool has(std::string_view name) const {
        return store_.find(normalize(name)) != store_.end();
    }

    void erase(std::string_view name) {
        store_.erase(normalize(name));
    }

private:
    static std::string normalize(std::string_view name) {
        std::string s(name.substr(0, 8));
        if (s.size() < 8) s.append(8 - s.size(), ' ');
        return s;
    }

    std::unordered_map<std::string, std::vector<double>> store_;
};
