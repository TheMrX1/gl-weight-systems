#include "chord/gl_weight.hpp"

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace chord {
namespace {

std::uint64_t splitmix64(std::uint64_t& state) {
    state += 0x9e3779b97f4a7c15ULL;
    auto value = state;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

void rotate_cycle_to_minimum(Cycle& cycle) {
    if (cycle.empty()) {
        return;
    }
    const auto minimum = std::min_element(cycle.begin(), cycle.end());
    std::rotate(cycle.begin(), minimum, cycle.end());
}

Cycle rotate_left(const Cycle& cycle, std::size_t amount) {
    if (cycle.empty()) {
        return {};
    }
    Cycle result = cycle;
    amount %= result.size();
    std::rotate(result.begin(), result.begin() + static_cast<std::ptrdiff_t>(amount), result.end());
    return result;
}

std::string serialize(const Cycles& cycles) {
    std::string key;
    std::size_t size = 1U;
    for (const auto& cycle : cycles) {
        size += cycle.size() + 1U;
    }
    key.reserve(size);
    key.push_back(static_cast<char>(cycles.size()));
    for (const auto& cycle : cycles) {
        key.push_back(static_cast<char>(cycle.size()));
        for (const auto label : cycle) {
            key.push_back(static_cast<char>(label));
        }
    }
    return key;
}

bool is_increasing_identity_cycle(const Cycle& cycle) {
    for (std::size_t index = 0; index < cycle.size(); ++index) {
        if (cycle[index] != index) {
            return false;
        }
    }
    return true;
}

}  

GlEvaluator::GlEvaluator(std::uint64_t assignment_seed, bool traceless_specialization)
    : primes_{1'000'000'007ULL, 1'000'000'009ULL, 998'244'353ULL} {
    for (std::size_t lane = 0; lane < kResidueLanes; ++lane) {
        auto state = assignment_seed ^ (0xd1b54a32d192ed03ULL * (lane + 1U));
        for (std::size_t index = 0; index < casimirs_[lane].size(); ++index) {
            casimirs_[lane][index] = 1U + splitmix64(state) % (primes_[lane] - 1U);
        }
        if (traceless_specialization) {
            casimirs_[lane][1U] = 0U;
        }
    }
    cache_.reserve(1U << 18U);
}

Residues GlEvaluator::evaluate_chord_word(const Word& word) {
    return evaluate_cycles(word_to_involution_cycles(word));
}

Residues GlEvaluator::evaluate_cycles(const Cycles& cycles) {
    recurrence_deadline_ = recurrence_budget_ == 0U
                               ? std::numeric_limits<std::uint64_t>::max()
                               : stats_.recurrence_nodes + recurrence_budget_;
    return evaluate_impl(cycles, 0U);
}

void GlEvaluator::clear_cache() {
    cache_.clear();
    stats_ = {};
}

void GlEvaluator::clear_memo() {
    cache_.clear();
}

void GlEvaluator::set_recurrence_budget(std::uint64_t nodes) {
    recurrence_budget_ = nodes;
}

EvaluatorStats GlEvaluator::stats() const {
    auto result = stats_;
    result.cache_entries = cache_.size();
    return result;
}

const std::array<std::uint64_t, kResidueLanes>& GlEvaluator::primes() const {
    return primes_;
}

const std::array<std::array<std::uint64_t, 23>, kResidueLanes>&
GlEvaluator::casimir_values() const {
    return casimirs_;
}

Residues GlEvaluator::one() const {
    Residues result;
    result.value.fill(1U);
    return result;
}

Residues GlEvaluator::zero() const {
    return {};
}

bool GlEvaluator::is_zero(const Residues& value) const {
    return std::all_of(value.value.begin(), value.value.end(),
                       [](std::uint64_t lane) { return lane == 0U; });
}

Residues GlEvaluator::casimir(std::size_t index) const {
    if (index >= casimirs_[0].size()) {
        throw std::out_of_range("Casimir index exceeds configured maximum 22");
    }
    Residues result;
    for (std::size_t lane = 0; lane < kResidueLanes; ++lane) {
        result.value[lane] = casimirs_[lane][index];
    }
    return result;
}

Residues GlEvaluator::add(const Residues& lhs, const Residues& rhs) const {
    Residues result;
    for (std::size_t lane = 0; lane < kResidueLanes; ++lane) {
        const auto sum = lhs.value[lane] + rhs.value[lane];
        result.value[lane] = sum >= primes_[lane] ? sum - primes_[lane] : sum;
    }
    return result;
}

Residues GlEvaluator::subtract(const Residues& lhs, const Residues& rhs) const {
    Residues result;
    for (std::size_t lane = 0; lane < kResidueLanes; ++lane) {
        result.value[lane] = lhs.value[lane] >= rhs.value[lane]
                                 ? lhs.value[lane] - rhs.value[lane]
                                 : primes_[lane] - (rhs.value[lane] - lhs.value[lane]);
    }
    return result;
}

Residues GlEvaluator::multiply(const Residues& lhs, const Residues& rhs) const {
    Residues result;
    for (std::size_t lane = 0; lane < kResidueLanes; ++lane) {
        const auto product = static_cast<unsigned __int128>(lhs.value[lane]) * rhs.value[lane];
        result.value[lane] = static_cast<std::uint64_t>(product % primes_[lane]);
    }
    return result;
}

Residues GlEvaluator::power(Residues base, std::size_t exponent) const {
    auto result = one();
    while (exponent != 0U) {
        if ((exponent & 1U) != 0U) {
            result = multiply(result, base);
        }
        exponent >>= 1U;
        if (exponent != 0U) {
            base = multiply(base, base);
        }
    }
    return result;
}

GlEvaluator::NormalizedState GlEvaluator::normalize(const Cycles& input) const {
    NormalizedState result;
    Cycles cycles;
    std::vector<std::uint8_t> labels;
    for (const auto& cycle : input) {
        if (cycle.empty()) {
            ++result.empty_cycles;
            continue;
        }
        cycles.push_back(cycle);
        labels.insert(labels.end(), cycle.begin(), cycle.end());
    }
    if (cycles.empty()) {
        return result;
    }

    std::sort(labels.begin(), labels.end());
    if (std::adjacent_find(labels.begin(), labels.end()) != labels.end()) {
        throw std::invalid_argument("permutation cycles must use every label exactly once");
    }

    std::array<std::uint8_t, 256> compression;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        compression[labels[index]] = static_cast<std::uint8_t>(index);
    }
    for (auto& cycle : cycles) {
        for (auto& label : cycle) {
            label = static_cast<std::uint8_t>(compression[label]);
        }
    }

    const auto vertex_count = labels.size();
    Cycles best;
    bool initialized = false;
    auto candidate = cycles;

    for (std::size_t shift = 0; shift < vertex_count; ++shift) {
        if (shift != 0U) {
            for (auto& cycle : candidate) {
                for (auto& label : cycle) {
                    label = static_cast<std::uint8_t>(
                        label == 0U ? vertex_count - 1U
                                    : static_cast<std::size_t>(label) - 1U);
                }
            }
        }
        for (auto& cycle : candidate) {
            rotate_cycle_to_minimum(cycle);
        }
        std::sort(candidate.begin(), candidate.end());

        std::size_t prefix_vertices = 0U;
        std::size_t prefix_max = 0U;
        for (std::size_t cycle_index = 0; cycle_index + 1U < candidate.size(); ++cycle_index) {
            prefix_vertices += candidate[cycle_index].size();
            for (const auto label : candidate[cycle_index]) {
                prefix_max = std::max(prefix_max, static_cast<std::size_t>(label));
            }
            if (prefix_max + 1U == prefix_vertices) {
                result.decomposable = true;
                result.left.assign(candidate.begin(),
                                   candidate.begin() + static_cast<std::ptrdiff_t>(cycle_index + 1U));
                result.right.assign(candidate.begin() + static_cast<std::ptrdiff_t>(cycle_index + 1U),
                                    candidate.end());
                for (auto& cycle : result.right) {
                    for (auto& label : cycle) {
                        label = static_cast<std::uint8_t>(
                            static_cast<std::size_t>(label) - prefix_vertices);
                    }
                }
                return result;
            }
        }

        if (!initialized || candidate < best) {
            best = candidate;
            initialized = true;
        }
    }

    result.cycles = std::move(best);
    result.key = serialize(result.cycles);
    return result;
}

Residues GlEvaluator::evaluate_impl(const Cycles& cycles, std::size_t depth) {
    ++stats_.calls;
    stats_.max_depth = std::max(stats_.max_depth, depth);
    const auto state = normalize(cycles);
    const auto restore_empty_cycles = [this, count = state.empty_cycles](Residues value) {
        return count == 0U ? value : multiply(power(casimir(0U), count), value);
    };
    if (state.decomposable) {
        ++stats_.decompositions;
        const auto left_value = evaluate_impl(state.left, depth + 1U);
        if (is_zero(left_value)) {
            return zero();
        }
        return restore_empty_cycles(
            multiply(left_value, evaluate_impl(state.right, depth + 1U)));
    }
    if (state.cycles.empty()) {
        return restore_empty_cycles(one());
    }

    if (const auto found = cache_.find(state.key); found != cache_.end()) {
        ++stats_.cache_hits;
        return restore_empty_cycles(found->second);
    }

    if (state.cycles.size() == 1U && is_increasing_identity_cycle(state.cycles.front())) {
        const auto result = casimir(state.cycles.front().size());
        cache_.emplace(state.key, result);
        return restore_empty_cycles(result);
    }

    const auto& first = state.cycles.front();
    std::size_t k = 0U;
    while (k < first.size() && first[k] == k) {
        ++k;
    }
    if (k == first.size() || first[k] <= k) {
        throw std::logic_error("normalized irreducible state has no admissible recurrence position");
    }
    const auto p = first[k];
    const auto predecessor = static_cast<std::uint8_t>(p - 1U);

    std::size_t predecessor_cycle = 0U;
    std::size_t predecessor_position = 0U;
    bool located = false;
    for (std::size_t cycle_index = 0; cycle_index < state.cycles.size() && !located; ++cycle_index) {
        const auto& cycle = state.cycles[cycle_index];
        const auto iterator = std::find(cycle.begin(), cycle.end(), predecessor);
        if (iterator != cycle.end()) {
            predecessor_cycle = cycle_index;
            predecessor_position = static_cast<std::size_t>(iterator - cycle.begin());
            located = true;
        }
    }
    if (!located) {
        throw std::logic_error("predecessor label is missing from permutation");
    }

    auto swapped = state.cycles;
    for (auto& cycle : swapped) {
        for (auto& label : cycle) {
            if (label == p) {
                label = predecessor;
            } else if (label == predecessor) {
                label = p;
            }
        }
    }

    Cycles positive;
    Cycles negative;
    if (predecessor_cycle != 0U) {
        const auto alpha = rotate_left(first, k);
        const auto beta = rotate_left(state.cycles[predecessor_cycle], predecessor_position);
        for (std::size_t index = 1U; index < state.cycles.size(); ++index) {
            if (index != predecessor_cycle) {
                positive.push_back(state.cycles[index]);
                negative.push_back(state.cycles[index]);
            }
        }
        Cycle joined_positive = alpha;
        joined_positive.insert(joined_positive.end(), beta.begin() + 1, beta.end());
        Cycle joined_negative(alpha.begin() + 1, alpha.end());
        joined_negative.insert(joined_negative.end(), beta.begin(), beta.end());
        positive.insert(positive.begin(), std::move(joined_positive));
        negative.insert(negative.begin(), std::move(joined_negative));
    } else {
        if (predecessor_position <= k) {
            throw std::logic_error("recurrence predecessor must occur after the pivot");
        }
        const auto rotated = rotate_left(first, k);
        const auto distance = predecessor_position - k;
        Cycle alpha(rotated.begin(), rotated.begin() + static_cast<std::ptrdiff_t>(distance));
        Cycle beta(rotated.begin() + static_cast<std::ptrdiff_t>(distance), rotated.end());

        positive.push_back(alpha);
        positive.emplace_back(beta.begin() + 1, beta.end());
        negative.emplace_back(alpha.begin() + 1, alpha.end());
        negative.push_back(beta);
        positive.insert(positive.end(), state.cycles.begin() + 1, state.cycles.end());
        negative.insert(negative.end(), state.cycles.begin() + 1, state.cycles.end());
    }

    if (stats_.recurrence_nodes >= recurrence_deadline_) {
        throw EvaluationBudgetExceeded();
    }
    ++stats_.recurrence_nodes;
    auto result = add(evaluate_impl(swapped, depth + 1U),
                      subtract(evaluate_impl(positive, depth + 1U),
                               evaluate_impl(negative, depth + 1U)));
    cache_.emplace(state.key, result);
    return restore_empty_cycles(result);
}

std::string residues_to_string(const Residues& residues) {
    std::ostringstream stream;
    for (std::size_t lane = 0; lane < kResidueLanes; ++lane) {
        if (lane != 0U) {
            stream << ',';
        }
        stream << residues.value[lane];
    }
    return stream.str();
}

}  
