#pragma once
#include "chord/gl_weight.hpp"
#include <memory>
#include <limits>

namespace chord {

using ExactInt = __int128;
using MonomialKey = unsigned __int128;
struct ExactTerm { MonomialKey monomial; ExactInt coefficient; };
using ExactPolynomial = std::vector<ExactTerm>;
ExactPolynomial polynomial_sum(const ExactPolynomial&, const ExactPolynomial&, int sign = 1);
std::string exact_integer_string(ExactInt);
std::string polynomial_string(const ExactPolynomial&);
std::string polynomial_json(const ExactPolynomial&);

class SymbolicGlEvaluator {
public:
    SymbolicGlEvaluator();
    ExactPolynomial evaluate_chord_word(const Word&);
private:
    using NodeId = std::uint32_t;
    enum class Op : std::uint8_t { Zero, One, Casimir, Add, Subtract, Multiply };
    struct Node { Op op; NodeId left{}, right{}; };
    struct NormalizedState {
        Cycles cycles; bool decomposable{}; Cycles left, right;
        std::size_t empty_cycles{}; std::string key;
    };
    std::vector<Node> nodes_;
    std::unordered_map<std::string, NodeId> cache_;
    EvaluatorStats stats_{};
    std::uint64_t recurrence_deadline_ = std::numeric_limits<std::uint64_t>::max();
    std::array<NodeId, 23> leaves_{};
    NodeId evaluate_impl(const Cycles&, std::size_t);
    NormalizedState normalize(const Cycles&) const;
    NodeId one() const { return 1; }
    NodeId zero() const { return 0; }
    bool is_zero(NodeId v) const { return v == 0; }
    NodeId casimir(std::size_t k) const { return leaves_.at(k); }
    NodeId add(NodeId, NodeId);
    NodeId subtract(NodeId, NodeId);
    NodeId multiply(NodeId, NodeId);
    NodeId power(NodeId, std::size_t);
    NodeId node(Op, NodeId, NodeId);
    ExactPolynomial expand(NodeId);
};
}
