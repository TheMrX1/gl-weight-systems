#pragma once

#include "chord/chord_diagram.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace chord {

class EvaluationBudgetExceeded : public std::runtime_error {
public:
    EvaluationBudgetExceeded() : std::runtime_error("gl recurrence-node budget exceeded") {}
};

inline constexpr std::size_t kResidueLanes = 3;

struct Residues {
    std::array<std::uint64_t, kResidueLanes> value{};

    friend bool operator==(const Residues&, const Residues&) = default;
};

struct EvaluatorStats {
    std::uint64_t calls{};
    std::uint64_t cache_hits{};
    std::uint64_t recurrence_nodes{};
    std::uint64_t decompositions{};
    std::size_t cache_entries{};
    std::size_t max_depth{};
};

class GlEvaluator {
public:
    explicit GlEvaluator(std::uint64_t assignment_seed = 0x9e3779b97f4a7c15ULL,
                         bool traceless_specialization = false);

    void assign_moduli(const std::array<std::uint64_t,kResidueLanes>& moduli) {
        primes_ = moduli;
        clear_cache();
    }
    void assign_casimirs(const std::array<std::array<std::uint64_t,23>,kResidueLanes>& values) {
        casimirs_ = values;
        clear_cache();
    }
    Residues evaluate_chord_word(const Word& word);
    Residues evaluate_cycles(const Cycles& cycles);
    void clear_cache();
    void clear_memo();
    void set_recurrence_budget(std::uint64_t nodes);
    [[nodiscard]] EvaluatorStats stats() const;
    [[nodiscard]] const std::array<std::uint64_t, kResidueLanes>& primes() const;
    [[nodiscard]] const std::array<std::array<std::uint64_t, 23>, kResidueLanes>&
    casimir_values() const;

private:
    struct NormalizedState {
        Cycles cycles;
        bool decomposable{};
        Cycles left;
        Cycles right;
        std::size_t empty_cycles{};
        std::string key;
    };

    Residues evaluate_impl(const Cycles& cycles, std::size_t depth);
    NormalizedState normalize(const Cycles& cycles) const;
    Residues casimir(std::size_t index) const;
    Residues one() const;
    Residues zero() const;
    bool is_zero(const Residues& value) const;
    Residues add(const Residues& lhs, const Residues& rhs) const;
    Residues subtract(const Residues& lhs, const Residues& rhs) const;
    Residues multiply(const Residues& lhs, const Residues& rhs) const;
    Residues power(Residues base, std::size_t exponent) const;

    std::array<std::uint64_t, kResidueLanes> primes_{};
    std::array<std::array<std::uint64_t, 23>, kResidueLanes> casimirs_{};
    std::unordered_map<std::string, Residues> cache_;
    EvaluatorStats stats_{};
    std::uint64_t recurrence_budget_{};
    std::uint64_t recurrence_deadline_{};
};

std::string residues_to_string(const Residues& residues);

}  
