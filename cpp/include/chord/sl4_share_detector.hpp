#pragma once

#include "chord/chord_diagram.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace chord {

struct TwoArcShare {
    Word first;
    Word second;
};

struct Block2x2 {
    std::array<std::uint32_t, 4> value{};

    friend bool operator==(const Block2x2&, const Block2x2&) = default;
};

struct ShareEvaluationStats {
    std::size_t path_width{};
    std::size_t peak_states{};
    std::uint64_t transitions{};
};

struct ChannelOddBlock {
    std::size_t channel{};
    std::array<int, 3> highest_weight{};
    std::size_t multiplicity{};
    std::vector<std::uint32_t> value;
};

class Sl4ShareDetector {
public:
    static constexpr std::uint32_t modulus = 65'521U;

    Sl4ShareDetector();

    [[nodiscard]] Block2x2 evaluate(const TwoArcShare& share);
    [[nodiscard]] Block2x2 evaluate_reference(const TwoArcShare& share);
    [[nodiscard]] Block2x2 swap_odd(const TwoArcShare& share);
    [[nodiscard]] Block2x2 mutation_odd(const TwoArcShare& share, int mutation_kind);
    [[nodiscard]] std::vector<ChannelOddBlock> channel_swap_odds(
        const TwoArcShare& share);
    [[nodiscard]] std::vector<ChannelOddBlock> channel_summary() const;
    [[nodiscard]] Block2x2 tensor_square_block(const TwoArcShare& share);
    [[nodiscard]] Block2x2 tensor_square_block_sparse(const TwoArcShare& share);
    [[nodiscard]] Block2x2 tensor_square_mutation_odd(const TwoArcShare& share,
                                                      int mutation_kind);
    [[nodiscard]] Block2x2 tensor_square_mutation_odd_sparse(const TwoArcShare& share,
                                                              int mutation_kind);
    [[nodiscard]] ShareEvaluationStats stats() const;
    void set_sparse_state_cap(std::size_t cap);
    void self_test() const;

private:
    static constexpr std::size_t dimension = 20U;
    static constexpr std::size_t matrix_size = dimension * dimension;

    struct DenseMatrix {
        std::array<std::uint16_t, matrix_size> value{};
    };

    struct MatrixState {
        std::array<std::uint16_t, matrix_size> value{};
    };

    struct SparseEntry {
        std::uint8_t row{};
        std::uint8_t column{};
        std::uint16_t value{};
    };

    struct ColumnEvaluation {
        std::array<std::uint32_t, 2> value{};
        ShareEvaluationStats stats{};
    };

    struct StateEvaluation {
        MatrixState value{};
        ShareEvaluationStats stats{};
    };

    struct MultiplicityChannel {
        std::array<int, 3> highest_weight{};
        std::vector<DenseMatrix> positive;
        std::vector<DenseMatrix> negative;
        std::vector<std::size_t> projection_rows;
        std::vector<std::uint32_t> projection_inverse;
        std::vector<std::uint32_t> gram;
        std::vector<std::uint32_t> gram_inverse;
    };

    using SparseMatrix = std::vector<SparseEntry>;

    std::array<DenseMatrix, 16> generators_{};
    std::array<SparseMatrix, 15> primal_generators_{};
    std::array<SparseMatrix, 15> dual_generators_{};
    DenseMatrix adjoint_linear_{};
    DenseMatrix adjoint_quadratic_{};
    std::size_t projection_first_{};
    std::size_t projection_second_{};
    std::uint32_t projection_inverse_determinant_{};
    Block2x2 adjoint_gram_{};
    Block2x2 adjoint_gram_inverse_{};
    std::vector<MultiplicityChannel> channels_;
    std::array<DenseMatrix, 2> tensor_square_highest_{};
    std::array<std::size_t, 2> tensor_square_projection_rows_{};
    Block2x2 tensor_square_projection_inverse_{};
    ShareEvaluationStats stats_{};
    std::size_t sparse_state_cap_{4'000'000U};

    [[nodiscard]] MatrixState apply_generator(const MatrixState& input,
                                              const SparseMatrix& generator,
                                              bool left,
                                              bool tensor_square = false) const;
    [[nodiscard]] std::array<std::uint32_t, 2> project_matrix(
        const MatrixState& matrix) const;
    [[nodiscard]] ColumnEvaluation evaluate_column(const TwoArcShare& share,
                                                   bool dense_states,
                                                   std::size_t input);
    [[nodiscard]] StateEvaluation contract(const TwoArcShare& share,
                                           bool dense_states,
                                           const MatrixState& initial,
                                           bool tensor_square = false);
    void build_multiplicity_channels();
    void build_tensor_square_channel();
    [[nodiscard]] Block2x2 tensor_square_block_impl(const TwoArcShare& share,
                                                     bool dense_states);
    [[nodiscard]] Block2x2 tensor_square_mutation_odd_impl(const TwoArcShare& share,
                                                            int mutation_kind,
                                                            bool dense_states);
    [[nodiscard]] Block2x2 evaluate_impl(const TwoArcShare& share, bool dense_states);
};

TwoArcShare mutate_two_arc_share(const TwoArcShare& share, int mutation_kind);
TwoArcShare canonical_two_arc_share(const TwoArcShare& share);
std::vector<TwoArcShare> enumerate_two_arc_shares(std::size_t chord_count);
std::size_t two_arc_share_path_width(const TwoArcShare& share,
                                     bool tensor_square = false);
std::uint32_t tensor_square_sewing_score(const Block2x2& difference,
                                         const Block2x2& complement);

}  // namespace chord
