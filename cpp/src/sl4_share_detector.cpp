#include "chord/sl4_share_detector.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace chord {
namespace {

constexpr std::size_t kAmbientDimension = 24U;
constexpr std::size_t kRepresentationDimension = 20U;
constexpr std::uint32_t kModulus = Sl4ShareDetector::modulus;

std::uint32_t add_mod(std::uint32_t left, std::uint32_t right) {
    const auto sum = left + right;
    return sum >= kModulus ? sum - kModulus : sum;
}

std::uint32_t subtract_mod(std::uint32_t left, std::uint32_t right) {
    return left >= right ? left - right : kModulus - (right - left);
}

std::uint32_t multiply_mod(std::uint32_t left, std::uint32_t right) {
    return static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(left) * right) % kModulus);
}

std::uint32_t power_mod(std::uint32_t base, std::uint32_t exponent) {
    std::uint32_t result = 1U;
    while (exponent != 0U) {
        if ((exponent & 1U) != 0U) {
            result = multiply_mod(result, base);
        }
        exponent >>= 1U;
        if (exponent != 0U) {
            base = multiply_mod(base, base);
        }
    }
    return result;
}

std::vector<std::uint32_t> invert_matrix(std::vector<std::uint32_t> matrix,
                                         std::size_t size) {
    std::vector<std::uint32_t> augmented(size * 2U * size, 0U);
    const auto at = [size](std::size_t row, std::size_t column) {
        return row * (2U * size) + column;
    };
    for (std::size_t row = 0U; row < size; ++row) {
        for (std::size_t column = 0U; column < size; ++column) {
            augmented[at(row, column)] = matrix[row * size + column];
        }
        augmented[at(row, size + row)] = 1U;
    }
    for (std::size_t column = 0U; column < size; ++column) {
        std::size_t pivot = column;
        while (pivot < size && augmented[at(pivot, column)] == 0U) {
            ++pivot;
        }
        if (pivot == size) {
            throw std::logic_error("singular modular matrix");
        }
        for (std::size_t item = 0U; item < 2U * size; ++item) {
            std::swap(augmented[at(column, item)], augmented[at(pivot, item)]);
        }
        const auto inverse = power_mod(augmented[at(column, column)], kModulus - 2U);
        for (std::size_t item = 0U; item < 2U * size; ++item) {
            augmented[at(column, item)] = multiply_mod(augmented[at(column, item)], inverse);
        }
        for (std::size_t row = 0U; row < size; ++row) {
            if (row == column || augmented[at(row, column)] == 0U) {
                continue;
            }
            const auto factor = augmented[at(row, column)];
            for (std::size_t item = 0U; item < 2U * size; ++item) {
                augmented[at(row, item)] = subtract_mod(
                    augmented[at(row, item)], multiply_mod(factor, augmented[at(column, item)]));
            }
        }
    }
    std::vector<std::uint32_t> result(size * size);
    for (std::size_t row = 0U; row < size; ++row) {
        for (std::size_t column = 0U; column < size; ++column) {
            result[row * size + column] = augmented[at(row, size + column)];
        }
    }
    return result;
}

std::vector<std::vector<std::uint32_t>> nullspace(
    std::vector<std::vector<std::uint32_t>> matrix, std::size_t columns) {
    std::vector<std::size_t> pivot_columns;
    std::size_t pivot_row = 0U;
    for (std::size_t column = 0U; column < columns && pivot_row < matrix.size(); ++column) {
        std::size_t pivot = pivot_row;
        while (pivot < matrix.size() && matrix[pivot][column] == 0U) {
            ++pivot;
        }
        if (pivot == matrix.size()) {
            continue;
        }
        std::swap(matrix[pivot_row], matrix[pivot]);
        const auto inverse = power_mod(matrix[pivot_row][column], kModulus - 2U);
        for (std::size_t item = column; item < columns; ++item) {
            matrix[pivot_row][item] = multiply_mod(matrix[pivot_row][item], inverse);
        }
        for (std::size_t row = 0U; row < matrix.size(); ++row) {
            if (row == pivot_row || matrix[row][column] == 0U) {
                continue;
            }
            const auto factor = matrix[row][column];
            for (std::size_t item = column; item < columns; ++item) {
                matrix[row][item] = subtract_mod(
                    matrix[row][item], multiply_mod(factor, matrix[pivot_row][item]));
            }
        }
        pivot_columns.push_back(column);
        ++pivot_row;
    }
    std::vector<bool> is_pivot(columns, false);
    for (const auto column : pivot_columns) {
        is_pivot[column] = true;
    }
    std::vector<std::vector<std::uint32_t>> result;
    for (std::size_t free_column = 0U; free_column < columns; ++free_column) {
        if (is_pivot[free_column]) {
            continue;
        }
        std::vector<std::uint32_t> vector(columns, 0U);
        vector[free_column] = 1U;
        for (std::size_t row = 0U; row < pivot_columns.size(); ++row) {
            const auto value = matrix[row][free_column];
            vector[pivot_columns[row]] = value == 0U ? 0U : kModulus - value;
        }
        result.push_back(std::move(vector));
    }
    return result;
}

int signed_residue(std::uint32_t value) {
    return value <= kModulus / 2U ? static_cast<int>(value)
                                 : static_cast<int>(value) - static_cast<int>(kModulus);
}

Block2x2 subtract_blocks(const Block2x2& left, const Block2x2& right) {
    Block2x2 result;
    for (std::size_t index = 0U; index < result.value.size(); ++index) {
        result.value[index] = subtract_mod(left.value[index], right.value[index]);
    }
    return result;
}

std::size_t generator_index(std::size_t row, std::size_t column) {
    return 4U * row + column;
}

struct AmbientBasisVector {
    std::uint8_t tensor{};
    std::uint8_t wedge_first{};
    std::uint8_t wedge_second{};
};

std::array<AmbientBasisVector, kAmbientDimension> ambient_basis() {
    std::array<AmbientBasisVector, kAmbientDimension> result{};
    std::size_t position = 0U;
    for (std::size_t tensor = 0U; tensor < 4U; ++tensor) {
        for (std::size_t first = 0U; first < 4U; ++first) {
            for (std::size_t second = first + 1U; second < 4U; ++second) {
                result[position++] = AmbientBasisVector{
                    static_cast<std::uint8_t>(tensor),
                    static_cast<std::uint8_t>(first),
                    static_cast<std::uint8_t>(second)};
            }
        }
    }
    return result;
}

std::size_t ambient_index(std::size_t tensor, std::size_t first, std::size_t second) {
    if (first >= second) {
        throw std::logic_error("wedge basis indices must be increasing");
    }
    std::size_t wedge_index = 0U;
    for (std::size_t a = 0U; a < 4U; ++a) {
        for (std::size_t b = a + 1U; b < 4U; ++b) {
            if (a == first && b == second) {
                return 6U * tensor + wedge_index;
            }
            ++wedge_index;
        }
    }
    throw std::logic_error("invalid wedge basis index");
}

int permutation_sign(std::array<std::size_t, 3> values) {
    if (values[0] == values[1] || values[0] == values[2] || values[1] == values[2]) {
        return 0;
    }
    int inversions = 0;
    for (std::size_t left = 0U; left < values.size(); ++left) {
        for (std::size_t right = left + 1U; right < values.size(); ++right) {
            inversions += values[left] > values[right] ? 1 : 0;
        }
    }
    return inversions % 2 == 0 ? 1 : -1;
}

std::size_t wedge3_index(std::array<std::size_t, 3> values) {
    std::sort(values.begin(), values.end());
    for (std::size_t missing = 0U; missing < 4U; ++missing) {
        std::array<std::size_t, 3> expected{};
        std::size_t position = 0U;
        for (std::size_t value = 0U; value < 4U; ++value) {
            if (value != missing) {
                expected[position++] = value;
            }
        }
        if (values == expected) {
            return missing;
        }
    }
    throw std::logic_error("invalid exterior-cube basis vector");
}

Word canonical_share_word(const TwoArcShare& share) {
    Word joined = share.first;
    joined.push_back(255U);
    joined.insert(joined.end(), share.second.begin(), share.second.end());
    std::array<int, 256> mapping{};
    mapping.fill(-1);
    std::uint8_t next = 0U;
    for (auto& label : joined) {
        if (label == 255U) {
            continue;
        }
        if (mapping[label] < 0) {
            mapping[label] = next++;
        }
        label = static_cast<std::uint8_t>(mapping[label]);
    }
    return joined;
}

TwoArcShare share_from_key(const Word& key) {
    const auto separator = std::find(key.begin(), key.end(), 255U);
    if (separator == key.end()) {
        throw std::logic_error("two-arc share key has no separator");
    }
    return TwoArcShare{Word(key.begin(), separator), Word(separator + 1, key.end())};
}

void enumerate_matchings(Word& positions, std::uint8_t next_label,
                         const std::function<void(const Word&)>& visitor) {
    const auto empty = std::find(positions.begin(), positions.end(), 255U);
    if (empty == positions.end()) {
        visitor(positions);
        return;
    }
    const auto first = static_cast<std::size_t>(empty - positions.begin());
    positions[first] = next_label;
    for (std::size_t second = first + 1U; second < positions.size(); ++second) {
        if (positions[second] != 255U) {
            continue;
        }
        positions[second] = next_label;
        enumerate_matchings(positions, static_cast<std::uint8_t>(next_label + 1U), visitor);
        positions[second] = 255U;
    }
    positions[first] = 255U;
}

}  

TwoArcShare mutate_two_arc_share(const TwoArcShare& share, int mutation_kind) {
    auto reverse = [](Word word) {
        std::reverse(word.begin(), word.end());
        return word;
    };
    switch (mutation_kind) {
        case 0:
            return share;
        case 1:
            return TwoArcShare{reverse(share.first), reverse(share.second)};
        case 2:
            return TwoArcShare{share.second, share.first};
        case 3:
            return TwoArcShare{reverse(share.second), reverse(share.first)};
        default:
            throw std::invalid_argument("share mutation kind must be in {0,1,2,3}");
    }
}

TwoArcShare canonical_two_arc_share(const TwoArcShare& share) {
    Word best;
    bool initialized = false;
    for (int mutation = 0; mutation < 4; ++mutation) {
        auto candidate = canonical_share_word(mutate_two_arc_share(share, mutation));
        if (!initialized || candidate < best) {
            best = std::move(candidate);
            initialized = true;
        }
    }
    return share_from_key(best);
}

std::vector<TwoArcShare> enumerate_two_arc_shares(std::size_t chord_count) {
    if (chord_count == 0U || chord_count > 15U) {
        throw std::invalid_argument("share enumeration supports 1..15 chords");
    }
    std::map<Word, TwoArcShare> unique;
    Word matching(2U * chord_count, 255U);
    enumerate_matchings(matching, 0U, [&](const Word& word) {
        for (std::size_t cut = 1U; cut < word.size(); ++cut) {
            const TwoArcShare share{Word(word.begin(), word.begin() + static_cast<std::ptrdiff_t>(cut)),
                                    Word(word.begin() + static_cast<std::ptrdiff_t>(cut), word.end())};
            auto canonical = canonical_two_arc_share(share);
            unique.emplace(canonical_share_word(canonical), std::move(canonical));
        }
    });
    std::vector<TwoArcShare> result;
    result.reserve(unique.size());
    for (auto& [key, share] : unique) {
        static_cast<void>(key);
        result.push_back(std::move(share));
    }
    return result;
}

std::size_t two_arc_share_path_width(const TwoArcShare& share, bool tensor_square) {
    Word joined = share.first;
    joined.insert(joined.end(), share.second.begin(), share.second.end());
    const auto normalized = relabel_by_first_occurrence(joined);
    if (normalized.empty() || normalized.size() % 2U != 0U) {
        throw std::invalid_argument("invalid two-arc share");
    }
    const auto chord_count = normalized.size() / 2U;
    if (chord_count > 15U) {
        throw std::invalid_argument("share path-width supports at most fifteen chords");
    }
    std::vector<std::size_t> counts(chord_count, 0U);
    for (const auto label : normalized) {
        if (label >= chord_count) {
            throw std::invalid_argument("invalid share chord label");
        }
        ++counts[label];
    }
    if (std::any_of(counts.begin(), counts.end(), [](std::size_t count) { return count != 2U; })) {
        throw std::invalid_argument("each share chord must occur exactly twice");
    }

    Word left(normalized.begin(),
              normalized.begin() + static_cast<std::ptrdiff_t>(share.first.size()));
    std::reverse(left.begin(), left.end());
    Word right(normalized.begin() + static_cast<std::ptrdiff_t>(share.first.size()),
               normalized.end());
    if (tensor_square) {
        std::reverse(right.begin(), right.end());
    }
    const auto infinity = std::numeric_limits<std::size_t>::max() / 4U;
    std::vector<std::vector<std::size_t>> cost(
        left.size() + 1U, std::vector<std::size_t>(right.size() + 1U, infinity));
    cost[0][0] = 0U;
    const auto open_width = [&](std::size_t left_count, std::size_t right_count) {
        std::uint32_t mask = 0U;
        for (std::size_t index = 0U; index < left_count; ++index) {
            mask ^= std::uint32_t{1} << left[index];
        }
        for (std::size_t index = 0U; index < right_count; ++index) {
            mask ^= std::uint32_t{1} << right[index];
        }
        return static_cast<std::size_t>(std::popcount(mask));
    };
    for (std::size_t left_count = 0U; left_count <= left.size(); ++left_count) {
        for (std::size_t right_count = 0U; right_count <= right.size(); ++right_count) {
            if (cost[left_count][right_count] == infinity) {
                continue;
            }
            if (left_count < left.size()) {
                cost[left_count + 1U][right_count] = std::min(
                    cost[left_count + 1U][right_count],
                    std::max(cost[left_count][right_count],
                             open_width(left_count + 1U, right_count)));
            }
            if (right_count < right.size()) {
                cost[left_count][right_count + 1U] = std::min(
                    cost[left_count][right_count + 1U],
                    std::max(cost[left_count][right_count],
                             open_width(left_count, right_count + 1U)));
            }
        }
    }
    return cost[left.size()][right.size()];
}

std::uint32_t tensor_square_sewing_score(const Block2x2& difference,
                                         const Block2x2& complement) {
    std::uint32_t result = 0U;
    for (std::size_t row = 0U; row < 2U; ++row) {
        for (std::size_t inner = 0U; inner < 2U; ++inner) {
            const auto product = multiply_mod(difference.value[inner * 2U + row],
                                              complement.value[row * 2U + inner]);
            result = row == 0U ? add_mod(result, product) : subtract_mod(result, product);
        }
    }
    return result;
}

Sl4ShareDetector::Sl4ShareDetector() {
    const auto basis = ambient_basis();

    std::array<std::array<std::uint32_t, kAmbientDimension>, 4> wedge_map{};
    for (std::size_t column = 0U; column < basis.size(); ++column) {
        const auto item = basis[column];
        const std::array<std::size_t, 3> triple{
            item.tensor, item.wedge_first, item.wedge_second};
        const auto sign = permutation_sign(triple);
        if (sign != 0) {
            wedge_map[wedge3_index(triple)][column] =
                sign > 0 ? 1U : kModulus - 1U;
        }
    }

    std::array<std::size_t, 4> pivots{};
    std::size_t pivot_count = 0U;
    for (std::size_t column = 0U; column < kAmbientDimension && pivot_count < 4U; ++column) {
        std::size_t pivot = pivot_count;
        while (pivot < 4U && wedge_map[pivot][column] == 0U) {
            ++pivot;
        }
        if (pivot == 4U) {
            continue;
        }
        std::swap(wedge_map[pivot], wedge_map[pivot_count]);
        const auto inverse = power_mod(wedge_map[pivot_count][column], kModulus - 2U);
        for (auto& value : wedge_map[pivot_count]) {
            value = multiply_mod(value, inverse);
        }
        for (std::size_t row = 0U; row < 4U; ++row) {
            if (row == pivot_count || wedge_map[row][column] == 0U) {
                continue;
            }
            const auto factor = wedge_map[row][column];
            for (std::size_t item = 0U; item < kAmbientDimension; ++item) {
                wedge_map[row][item] = subtract_mod(
                    wedge_map[row][item], multiply_mod(factor, wedge_map[pivot_count][item]));
            }
        }
        pivots[pivot_count++] = column;
    }
    if (pivot_count != 4U) {
        throw std::logic_error("exterior multiplication map must have rank four");
    }

    std::array<bool, kAmbientDimension> is_pivot{};
    for (const auto pivot : pivots) {
        is_pivot[pivot] = true;
    }
    std::array<std::size_t, kRepresentationDimension> free_columns{};
    std::size_t free_count = 0U;
    for (std::size_t column = 0U; column < kAmbientDimension; ++column) {
        if (!is_pivot[column]) {
            free_columns[free_count++] = column;
        }
    }
    if (free_count != kRepresentationDimension) {
        throw std::logic_error("[2,1] representation must have dimension twenty");
    }

    std::array<std::array<std::uint32_t, kRepresentationDimension>, kAmbientDimension> kernel{};
    for (std::size_t basis_index = 0U; basis_index < kRepresentationDimension; ++basis_index) {
        const auto free_column = free_columns[basis_index];
        kernel[free_column][basis_index] = 1U;
        for (std::size_t row = 0U; row < 4U; ++row) {
            kernel[pivots[row]][basis_index] =
                wedge_map[row][free_column] == 0U
                    ? 0U
                    : kModulus - wedge_map[row][free_column];
        }
    }

    for (std::size_t matrix_row = 0U; matrix_row < 4U; ++matrix_row) {
        for (std::size_t matrix_column = 0U; matrix_column < 4U; ++matrix_column) {
            std::array<std::array<std::uint32_t, kAmbientDimension>, kAmbientDimension> ambient{};
            for (std::size_t source = 0U; source < basis.size(); ++source) {
                const auto item = basis[source];
                if (item.tensor == matrix_column) {
                    ambient[ambient_index(matrix_row, item.wedge_first, item.wedge_second)][source] =
                        add_mod(ambient[ambient_index(matrix_row, item.wedge_first, item.wedge_second)][source], 1U);
                }
                if (item.wedge_first == matrix_column && matrix_row != item.wedge_second) {
                    const auto first = std::min<std::size_t>(matrix_row, item.wedge_second);
                    const auto second = std::max<std::size_t>(matrix_row, item.wedge_second);
                    const auto sign = matrix_row < item.wedge_second ? 1U : kModulus - 1U;
                    const auto target = ambient_index(item.tensor, first, second);
                    ambient[target][source] = add_mod(ambient[target][source], sign);
                }
                if (item.wedge_second == matrix_column && matrix_row != item.wedge_first) {
                    const auto first = std::min<std::size_t>(item.wedge_first, matrix_row);
                    const auto second = std::max<std::size_t>(item.wedge_first, matrix_row);
                    const auto sign = item.wedge_first < matrix_row ? 1U : kModulus - 1U;
                    const auto target = ambient_index(item.tensor, first, second);
                    ambient[target][source] = add_mod(ambient[target][source], sign);
                }
            }

            auto& generator = generators_[generator_index(matrix_row, matrix_column)];
            for (std::size_t output = 0U; output < kRepresentationDimension; ++output) {
                const auto ambient_output = free_columns[output];
                for (std::size_t input = 0U; input < kRepresentationDimension; ++input) {
                    std::uint64_t sum = 0U;
                    for (std::size_t middle = 0U; middle < kAmbientDimension; ++middle) {
                        sum += static_cast<std::uint64_t>(ambient[ambient_output][middle]) *
                               kernel[middle][input];
                    }
                    generator.value[output * dimension + input] =
                        static_cast<std::uint16_t>(sum % kModulus);
                }
            }
        }
    }

    const auto sparse = [&](const DenseMatrix& matrix) {
        SparseMatrix result;
        for (std::size_t row = 0U; row < dimension; ++row) {
            for (std::size_t column = 0U; column < dimension; ++column) {
                const auto value = matrix.value[row * dimension + column];
                if (value != 0U) {
                    result.push_back(SparseEntry{static_cast<std::uint8_t>(row),
                                                 static_cast<std::uint8_t>(column), value});
                }
            }
        }
        return result;
    };
    const auto combination = [&](const std::vector<std::pair<std::size_t, std::uint32_t>>& terms) {
        DenseMatrix result;
        for (std::size_t index = 0U; index < matrix_size; ++index) {
            std::uint32_t value = 0U;
            for (const auto [generator, coefficient] : terms) {
                value = add_mod(value, multiply_mod(coefficient, generators_[generator].value[index]));
            }
            result.value[index] = static_cast<std::uint16_t>(value);
        }
        return result;
    };

    std::size_t sl_index = 0U;
    for (std::size_t row = 0U; row < 4U; ++row) {
        for (std::size_t column = 0U; column < 4U; ++column) {
            if (row == column) {
                continue;
            }
            primal_generators_[sl_index] = sparse(generators_[generator_index(row, column)]);
            dual_generators_[sl_index] = sparse(generators_[generator_index(column, row)]);
            ++sl_index;
        }
    }
    const auto minus_one = kModulus - 1U;
    std::array<DenseMatrix, 3> cartan{};
    for (std::size_t index = 0U; index < 3U; ++index) {
        cartan[index] = combination({{generator_index(index, index), 1U},
                                     {generator_index(index + 1U, index + 1U), minus_one}});
        primal_generators_[sl_index + index] = sparse(cartan[index]);
    }
    const auto inverse_four = power_mod(4U, kModulus - 2U);
    constexpr std::array<std::array<std::uint32_t, 3>, 3> inverse_cartan_numerator{{
        {{3U, 2U, 1U}}, {{2U, 4U, 2U}}, {{1U, 2U, 3U}}}};
    for (std::size_t dual = 0U; dual < 3U; ++dual) {
        DenseMatrix matrix;
        for (std::size_t index = 0U; index < matrix_size; ++index) {
            std::uint32_t value = 0U;
            for (std::size_t primal = 0U; primal < 3U; ++primal) {
                const auto coefficient = multiply_mod(
                    inverse_cartan_numerator[dual][primal], inverse_four);
                value = add_mod(value, multiply_mod(coefficient, cartan[primal].value[index]));
            }
            matrix.value[index] = static_cast<std::uint16_t>(value);
        }
        dual_generators_[sl_index + dual] = sparse(matrix);
    }

    adjoint_linear_ = generators_[generator_index(0U, 3U)];
    for (std::size_t middle = 0U; middle < 4U; ++middle) {
        const auto& left = generators_[generator_index(0U, middle)];
        const auto& right = generators_[generator_index(middle, 3U)];
        for (std::size_t row = 0U; row < dimension; ++row) {
            for (std::size_t column = 0U; column < dimension; ++column) {
                std::uint64_t sum = adjoint_quadratic_.value[row * dimension + column];
                for (std::size_t inner = 0U; inner < dimension; ++inner) {
                    sum += static_cast<std::uint64_t>(left.value[row * dimension + inner]) *
                           right.value[inner * dimension + column];
                }
                adjoint_quadratic_.value[row * dimension + column] =
                    static_cast<std::uint16_t>(sum % kModulus);
            }
        }
    }

    bool found_projection = false;
    for (std::size_t first = 0U; first < matrix_size && !found_projection; ++first) {
        for (std::size_t second = first + 1U; second < matrix_size; ++second) {
            const auto determinant = subtract_mod(
                multiply_mod(adjoint_linear_.value[first], adjoint_quadratic_.value[second]),
                multiply_mod(adjoint_linear_.value[second], adjoint_quadratic_.value[first]));
            if (determinant != 0U) {
                projection_first_ = first;
                projection_second_ = second;
                projection_inverse_determinant_ = power_mod(determinant, kModulus - 2U);
                found_projection = true;
                break;
            }
        }
    }
    if (!found_projection) {
        throw std::logic_error("the two adjoint embeddings are linearly dependent");
    }

    DenseMatrix opposite_linear = generators_[generator_index(3U, 0U)];
    DenseMatrix opposite_quadratic;
    for (std::size_t middle = 0U; middle < 4U; ++middle) {
        const auto& left = generators_[generator_index(3U, middle)];
        const auto& right = generators_[generator_index(middle, 0U)];
        for (std::size_t row = 0U; row < dimension; ++row) {
            for (std::size_t column = 0U; column < dimension; ++column) {
                std::uint64_t sum = opposite_quadratic.value[row * dimension + column];
                for (std::size_t inner = 0U; inner < dimension; ++inner) {
                    sum += static_cast<std::uint64_t>(left.value[row * dimension + inner]) *
                           right.value[inner * dimension + column];
                }
                opposite_quadratic.value[row * dimension + column] =
                    static_cast<std::uint16_t>(sum % kModulus);
            }
        }
    }
    const std::array<const DenseMatrix*, 2> positive{&adjoint_linear_, &adjoint_quadratic_};
    const std::array<const DenseMatrix*, 2> negative{&opposite_linear, &opposite_quadratic};
    const auto block_at = [](std::size_t row, std::size_t column) {
        return column * 2U + row;
    };
    for (std::size_t row = 0U; row < 2U; ++row) {
        for (std::size_t column = 0U; column < 2U; ++column) {
            std::uint64_t trace = 0U;
            for (std::size_t i = 0U; i < dimension; ++i) {
                for (std::size_t j = 0U; j < dimension; ++j) {
                    trace += static_cast<std::uint64_t>(positive[row]->value[i * dimension + j]) *
                             negative[column]->value[j * dimension + i];
                }
            }
            adjoint_gram_.value[block_at(row, column)] =
                static_cast<std::uint32_t>(trace % kModulus);
        }
    }
    const auto determinant = subtract_mod(
        multiply_mod(adjoint_gram_.value[0], adjoint_gram_.value[3]),
        multiply_mod(adjoint_gram_.value[1], adjoint_gram_.value[2]));
    if (determinant == 0U) {
        throw std::logic_error("adjoint multiplicity Gram matrix is singular");
    }
    const auto inverse = power_mod(determinant, kModulus - 2U);
    adjoint_gram_inverse_.value = {
        multiply_mod(adjoint_gram_.value[3], inverse),
        multiply_mod(subtract_mod(0U, adjoint_gram_.value[1]), inverse),
        multiply_mod(subtract_mod(0U, adjoint_gram_.value[2]), inverse),
        multiply_mod(adjoint_gram_.value[0], inverse)};
    build_multiplicity_channels();
    build_tensor_square_channel();
    self_test();
}

void Sl4ShareDetector::build_multiplicity_channels() {
    using Weight = std::array<std::uint16_t, 3>;
    std::array<Weight, dimension> representation_weights{};
    for (std::size_t basis = 0U; basis < dimension; ++basis) {
        for (std::size_t cartan = 0U; cartan < 3U; ++cartan) {
            const auto& first = generators_[generator_index(cartan, cartan)];
            const auto& second = generators_[generator_index(cartan + 1U, cartan + 1U)];
            representation_weights[basis][cartan] = static_cast<std::uint16_t>(subtract_mod(
                first.value[basis * dimension + basis],
                second.value[basis * dimension + basis]));
            for (std::size_t column = 0U; column < dimension; ++column) {
                if (column != basis &&
                    (first.value[basis * dimension + column] != 0U ||
                     second.value[basis * dimension + column] != 0U)) {
                    throw std::logic_error("the constructed representation basis is not a weight basis");
                }
            }
        }
    }

    std::map<Weight, std::vector<std::size_t>> weight_positions;
    for (std::size_t row = 0U; row < dimension; ++row) {
        for (std::size_t column = 0U; column < dimension; ++column) {
            Weight weight{};
            for (std::size_t cartan = 0U; cartan < 3U; ++cartan) {
                weight[cartan] = static_cast<std::uint16_t>(subtract_mod(
                    representation_weights[row][cartan], representation_weights[column][cartan]));
            }
            weight_positions[weight].push_back(row * dimension + column);
        }
    }

    const auto extremal_vectors = [&](const Weight& weight, bool highest) {
        const auto found = weight_positions.find(weight);
        if (found == weight_positions.end()) {
            return std::vector<DenseMatrix>{};
        }
        const auto& positions = found->second;
        std::vector<std::vector<std::uint32_t>> equations(
            3U * matrix_size, std::vector<std::uint32_t>(positions.size(), 0U));
        for (std::size_t variable = 0U; variable < positions.size(); ++variable) {
            const auto row = positions[variable] / dimension;
            const auto column = positions[variable] % dimension;
            for (std::size_t root = 0U; root < 3U; ++root) {
                const auto generator = highest ? generator_index(root, root + 1U)
                                               : generator_index(root + 1U, root);
                const auto& matrix = generators_[generator];
                for (std::size_t output_row = 0U; output_row < dimension; ++output_row) {
                    const auto coefficient = matrix.value[output_row * dimension + row];
                    if (coefficient != 0U) {
                        auto& target = equations[root * matrix_size +
                                                 output_row * dimension + column][variable];
                        target = add_mod(target, coefficient);
                    }
                }
                for (std::size_t output_column = 0U; output_column < dimension; ++output_column) {
                    const auto coefficient = matrix.value[column * dimension + output_column];
                    if (coefficient != 0U) {
                        auto& target = equations[root * matrix_size +
                                                 row * dimension + output_column][variable];
                        target = subtract_mod(target, coefficient);
                    }
                }
            }
        }
        equations.erase(std::remove_if(equations.begin(), equations.end(),
                                       [](const auto& row) {
                                           return std::all_of(row.begin(), row.end(),
                                                              [](auto value) { return value == 0U; });
                                       }),
                        equations.end());
        const auto kernel = nullspace(std::move(equations), positions.size());
        std::vector<DenseMatrix> result;
        result.reserve(kernel.size());
        for (const auto& vector : kernel) {
            DenseMatrix matrix;
            for (std::size_t index = 0U; index < positions.size(); ++index) {
                matrix.value[positions[index]] = static_cast<std::uint16_t>(vector[index]);
            }
            result.push_back(std::move(matrix));
        }
        return result;
    };

    for (const auto& [weight, positions] : weight_positions) {
        static_cast<void>(positions);
        auto positive = extremal_vectors(weight, true);
        if (positive.size() < 2U) {
            continue;
        }
        Weight opposite_weight{};
        for (std::size_t item = 0U; item < 3U; ++item) {
            opposite_weight[item] = weight[item] == 0U
                                        ? 0U
                                        : static_cast<std::uint16_t>(kModulus - weight[item]);
        }
        auto negative = extremal_vectors(opposite_weight, false);
        if (negative.size() != positive.size()) {
            throw std::logic_error("highest/lowest multiplicities do not match");
        }
        const auto multiplicity = positive.size();

        std::vector<std::vector<std::uint32_t>> basis_rows(
            multiplicity, std::vector<std::uint32_t>(matrix_size));
        for (std::size_t basis = 0U; basis < multiplicity; ++basis) {
            for (std::size_t position = 0U; position < matrix_size; ++position) {
                basis_rows[basis][position] = positive[basis].value[position];
            }
        }
        std::vector<std::size_t> projection_rows;
        std::size_t pivot_row = 0U;
        for (std::size_t column = 0U; column < matrix_size && pivot_row < multiplicity; ++column) {
            std::size_t pivot = pivot_row;
            while (pivot < multiplicity && basis_rows[pivot][column] == 0U) {
                ++pivot;
            }
            if (pivot == multiplicity) {
                continue;
            }
            std::swap(basis_rows[pivot_row], basis_rows[pivot]);
            const auto inverse = power_mod(basis_rows[pivot_row][column], kModulus - 2U);
            for (std::size_t item = column; item < matrix_size; ++item) {
                basis_rows[pivot_row][item] = multiply_mod(basis_rows[pivot_row][item], inverse);
            }
            for (std::size_t row = 0U; row < multiplicity; ++row) {
                if (row == pivot_row || basis_rows[row][column] == 0U) {
                    continue;
                }
                const auto factor = basis_rows[row][column];
                for (std::size_t item = column; item < matrix_size; ++item) {
                    basis_rows[row][item] = subtract_mod(
                        basis_rows[row][item],
                        multiply_mod(factor, basis_rows[pivot_row][item]));
                }
            }
            projection_rows.push_back(column);
            ++pivot_row;
        }
        if (projection_rows.size() != multiplicity) {
            throw std::logic_error("highest-weight vectors are linearly dependent");
        }

        std::vector<std::uint32_t> projection(multiplicity * multiplicity);
        std::vector<std::uint32_t> gram(multiplicity * multiplicity);
        for (std::size_t row = 0U; row < multiplicity; ++row) {
            for (std::size_t column = 0U; column < multiplicity; ++column) {
                projection[row * multiplicity + column] =
                    positive[column].value[projection_rows[row]];
                std::uint64_t trace = 0U;
                for (std::size_t i = 0U; i < dimension; ++i) {
                    for (std::size_t j = 0U; j < dimension; ++j) {
                        trace += static_cast<std::uint64_t>(positive[row].value[i * dimension + j]) *
                                 negative[column].value[j * dimension + i];
                    }
                }
                gram[row * multiplicity + column] =
                    static_cast<std::uint32_t>(trace % kModulus);
            }
        }
        channels_.push_back(MultiplicityChannel{
            {signed_residue(weight[0]), signed_residue(weight[1]), signed_residue(weight[2])},
            std::move(positive), std::move(negative), std::move(projection_rows),
            invert_matrix(std::move(projection), multiplicity), gram,
            invert_matrix(std::move(gram), multiplicity)});
    }
}

void Sl4ShareDetector::build_tensor_square_channel() {
    using Weight = std::array<std::uint16_t, 3>;
    std::array<Weight, dimension> representation_weights{};
    for (std::size_t basis = 0U; basis < dimension; ++basis) {
        for (std::size_t cartan = 0U; cartan < 3U; ++cartan) {
            representation_weights[basis][cartan] = static_cast<std::uint16_t>(subtract_mod(
                generators_[generator_index(cartan, cartan)]
                    .value[basis * dimension + basis],
                generators_[generator_index(cartan + 1U, cartan + 1U)]
                    .value[basis * dimension + basis]));
        }
    }

    constexpr Weight target_weight{1U, 1U, 1U};
    std::vector<std::size_t> positions;
    for (std::size_t row = 0U; row < dimension; ++row) {
        for (std::size_t column = 0U; column < dimension; ++column) {
            Weight weight{};
            for (std::size_t cartan = 0U; cartan < 3U; ++cartan) {
                weight[cartan] = static_cast<std::uint16_t>(add_mod(
                    representation_weights[row][cartan],
                    representation_weights[column][cartan]));
            }
            if (weight == target_weight) {
                positions.push_back(row * dimension + column);
            }
        }
    }

    std::vector<std::vector<std::uint32_t>> equations(
        3U * matrix_size, std::vector<std::uint32_t>(positions.size(), 0U));
    for (std::size_t variable = 0U; variable < positions.size(); ++variable) {
        const auto row = positions[variable] / dimension;
        const auto column = positions[variable] % dimension;
        for (std::size_t root = 0U; root < 3U; ++root) {
            const auto& matrix = generators_[generator_index(root, root + 1U)];
            for (std::size_t output_row = 0U; output_row < dimension; ++output_row) {
                const auto coefficient = matrix.value[output_row * dimension + row];
                if (coefficient != 0U) {
                    auto& target = equations[root * matrix_size +
                                             output_row * dimension + column][variable];
                    target = add_mod(target, coefficient);
                }
            }
            for (std::size_t output_column = 0U; output_column < dimension; ++output_column) {
                const auto coefficient = matrix.value[output_column * dimension + column];
                if (coefficient != 0U) {
                    auto& target = equations[root * matrix_size +
                                             row * dimension + output_column][variable];
                    target = add_mod(target, coefficient);
                }
            }
        }
    }
    equations.erase(std::remove_if(equations.begin(), equations.end(),
                                   [](const auto& row) {
                                       return std::all_of(row.begin(), row.end(),
                                                          [](auto value) { return value == 0U; });
                                   }),
                    equations.end());
    const auto kernel = nullspace(std::move(equations), positions.size());
    if (kernel.size() != 2U) {
        throw std::logic_error("the (1,1,1) highest-weight space in [2,1] tensor square must be two-dimensional");
    }

    std::array<bool, 2> found{};
    for (const auto& vector : kernel) {
        DenseMatrix raw;
        for (std::size_t index = 0U; index < positions.size(); ++index) {
            raw.value[positions[index]] = static_cast<std::uint16_t>(vector[index]);
        }
        for (std::size_t parity = 0U; parity < 2U; ++parity) {
            if (found[parity]) {
                continue;
            }
            DenseMatrix candidate;
            bool nonzero = false;
            for (std::size_t row = 0U; row < dimension; ++row) {
                for (std::size_t column = 0U; column < dimension; ++column) {
                    const auto direct = raw.value[row * dimension + column];
                    const auto flipped = raw.value[column * dimension + row];
                    const auto value = parity == 0U ? add_mod(direct, flipped)
                                                    : subtract_mod(direct, flipped);
                    candidate.value[row * dimension + column] =
                        static_cast<std::uint16_t>(value);
                    nonzero = nonzero || value != 0U;
                }
            }
            if (!nonzero) {
                continue;
            }
            const auto first = std::find_if(candidate.value.begin(), candidate.value.end(),
                                            [](auto value) { return value != 0U; });
            const auto inverse = power_mod(*first, kModulus - 2U);
            for (auto& value : candidate.value) {
                value = static_cast<std::uint16_t>(multiply_mod(value, inverse));
            }
            tensor_square_highest_[parity] = std::move(candidate);
            found[parity] = true;
        }
    }
    if (!found[0] || !found[1]) {
        throw std::logic_error("failed to split the repeated 64 channel into symmetric and antisymmetric copies");
    }

    bool projection_found = false;
    for (std::size_t first = 0U; first < matrix_size && !projection_found; ++first) {
        for (std::size_t second = first + 1U; second < matrix_size; ++second) {
            const auto a = tensor_square_highest_[0].value[first];
            const auto b = tensor_square_highest_[1].value[first];
            const auto c = tensor_square_highest_[0].value[second];
            const auto d = tensor_square_highest_[1].value[second];
            const auto determinant = subtract_mod(multiply_mod(a, d), multiply_mod(b, c));
            if (determinant == 0U) {
                continue;
            }
            const auto inverse = power_mod(determinant, kModulus - 2U);
            tensor_square_projection_rows_ = {first, second};
            tensor_square_projection_inverse_.value = {
                multiply_mod(d, inverse),
                multiply_mod(subtract_mod(0U, b), inverse),
                multiply_mod(subtract_mod(0U, c), inverse),
                multiply_mod(a, inverse)};
            projection_found = true;
            break;
        }
    }
    if (!projection_found) {
        throw std::logic_error("tensor-square highest-weight vectors are linearly dependent");
    }
}

Sl4ShareDetector::MatrixState Sl4ShareDetector::apply_generator(
    const MatrixState& input, const SparseMatrix& generator, bool left,
    bool tensor_square) const {
    MatrixState output;
    const auto scaled = [](std::uint32_t coefficient, std::uint32_t value) {
        if (coefficient == 1U) {
            return value;
        }
        if (coefficient == 2U) {
            return add_mod(value, value);
        }
        if (coefficient == kModulus - 1U) {
            return value == 0U ? 0U : kModulus - value;
        }
        if (coefficient == kModulus - 2U) {
            const auto twice = add_mod(value, value);
            return twice == 0U ? 0U : kModulus - twice;
        }
        return multiply_mod(coefficient, value);
    };
    for (const auto& entry : generator) {
        if (left) {
            for (std::size_t column = 0U; column < dimension; ++column) {
                const auto target = entry.row * dimension + column;
                output.value[target] = static_cast<std::uint16_t>(add_mod(
                    output.value[target],
                    scaled(entry.value, input.value[entry.column * dimension + column])));
            }
        } else if (!tensor_square) {
            for (std::size_t row = 0U; row < dimension; ++row) {
                const auto target = row * dimension + entry.column;
                output.value[target] = static_cast<std::uint16_t>(add_mod(
                    output.value[target],
                    scaled(entry.value, input.value[row * dimension + entry.row])));
            }
        } else {
            for (std::size_t row = 0U; row < dimension; ++row) {
                const auto target = row * dimension + entry.row;
                output.value[target] = static_cast<std::uint16_t>(add_mod(
                    output.value[target],
                    scaled(entry.value, input.value[row * dimension + entry.column])));
            }
        }
    }
    return output;
}

std::array<std::uint32_t, 2> Sl4ShareDetector::project_matrix(
    const MatrixState& matrix) const {
    const auto a1 = adjoint_linear_.value[projection_first_];
    const auto a2 = adjoint_linear_.value[projection_second_];
    const auto b1 = adjoint_quadratic_.value[projection_first_];
    const auto b2 = adjoint_quadratic_.value[projection_second_];
    const auto y1 = matrix.value[projection_first_];
    const auto y2 = matrix.value[projection_second_];
    const auto linear = multiply_mod(
        subtract_mod(multiply_mod(y1, b2), multiply_mod(y2, b1)),
        projection_inverse_determinant_);
    const auto quadratic = multiply_mod(
        subtract_mod(multiply_mod(a1, y2), multiply_mod(a2, y1)),
        projection_inverse_determinant_);
    for (std::size_t index = 0U; index < matrix_size; ++index) {
        const auto expected = add_mod(multiply_mod(linear, adjoint_linear_.value[index]),
                                      multiply_mod(quadratic, adjoint_quadratic_.value[index]));
        if (matrix.value[index] != expected) {
            throw std::logic_error("share operator left the adjoint multiplicity space");
        }
    }
    return {linear, quadratic};
}

Block2x2 Sl4ShareDetector::evaluate(const TwoArcShare& share) {
    return evaluate_impl(share, true);
}

Block2x2 Sl4ShareDetector::evaluate_reference(const TwoArcShare& share) {
    return evaluate_impl(share, false);
}

Block2x2 Sl4ShareDetector::evaluate_impl(const TwoArcShare& share, bool dense_states) {
    const auto first = evaluate_column(share, dense_states, 0U);
    const auto second = evaluate_column(share, dense_states, 1U);
    stats_.path_width = std::max(first.stats.path_width, second.stats.path_width);
    stats_.peak_states = std::max(first.stats.peak_states, second.stats.peak_states);
    stats_.transitions = first.stats.transitions + second.stats.transitions;
    return Block2x2{{first.value[0], first.value[1], second.value[0], second.value[1]}};
}

Sl4ShareDetector::ColumnEvaluation Sl4ShareDetector::evaluate_column(
    const TwoArcShare& share, bool dense_states, std::size_t input) {
    if (input > 1U) {
        throw std::invalid_argument("adjoint input column must be zero or one");
    }
    MatrixState initial;
    const auto& initial_matrix = input == 0U ? adjoint_linear_ : adjoint_quadratic_;
    std::copy(initial_matrix.value.begin(), initial_matrix.value.end(), initial.value.begin());
    auto result = contract(share, dense_states, initial);
    return ColumnEvaluation{project_matrix(result.value), result.stats};
}

Sl4ShareDetector::StateEvaluation Sl4ShareDetector::contract(
    const TwoArcShare& share, bool dense_states, const MatrixState& initial,
    bool tensor_square) {
    Word joined = share.first;
    joined.insert(joined.end(), share.second.begin(), share.second.end());
    const auto normalized = relabel_by_first_occurrence(joined);
    if (normalized.size() != joined.size() || normalized.size() % 2U != 0U) {
        throw std::invalid_argument("invalid two-arc share");
    }
    const auto chord_count = normalized.size() / 2U;
    if (chord_count > 7U) {
        throw std::invalid_argument("the current share detector supports at most seven chords");
    }
    std::vector<std::size_t> counts(chord_count, 0U);
    for (const auto label : normalized) {
        if (label >= chord_count) {
            throw std::invalid_argument("invalid share chord label");
        }
        ++counts[label];
    }
    if (std::any_of(counts.begin(), counts.end(), [](std::size_t count) { return count != 2U; })) {
        throw std::invalid_argument("each share chord must occur exactly twice");
    }

    TwoArcShare canonical_input;
    canonical_input.first.assign(normalized.begin(),
                                 normalized.begin() + static_cast<std::ptrdiff_t>(share.first.size()));
    canonical_input.second.assign(
        normalized.begin() + static_cast<std::ptrdiff_t>(share.first.size()), normalized.end());

    Word left_sequence = canonical_input.first;
    std::reverse(left_sequence.begin(), left_sequence.end());
    Word right_sequence = canonical_input.second;
    if (tensor_square) {

        std::reverse(right_sequence.begin(), right_sequence.end());
    }
    const auto left_size = left_sequence.size();
    const auto right_size = right_sequence.size();
    const auto infinity = std::numeric_limits<std::size_t>::max() / 4U;
    std::vector<std::vector<std::size_t>> cost(left_size + 1U,
                                               std::vector<std::size_t>(right_size + 1U, infinity));
    cost[0][0] = 0U;
    const auto open_width = [&](std::size_t left, std::size_t right) {
        std::uint32_t mask = 0U;
        for (std::size_t index = 0U; index < left; ++index) {
            mask ^= std::uint32_t{1} << left_sequence[index];
        }
        for (std::size_t index = 0U; index < right; ++index) {
            mask ^= std::uint32_t{1} << right_sequence[index];
        }
        return static_cast<std::size_t>(std::popcount(mask));
    };
    for (std::size_t left = 0U; left <= left_size; ++left) {
        for (std::size_t right = 0U; right <= right_size; ++right) {
            if (cost[left][right] == infinity) {
                continue;
            }
            if (left < left_size) {
                const auto candidate = std::max(cost[left][right], open_width(left + 1U, right));
                if (candidate < cost[left + 1U][right]) {
                    cost[left + 1U][right] = candidate;
                }
            }
            if (right < right_size) {
                const auto candidate = std::max(cost[left][right], open_width(left, right + 1U));
                if (candidate < cost[left][right + 1U]) {
                    cost[left][right + 1U] = candidate;
                }
            }
        }
    }

    const auto target_width = cost[left_size][right_size];
    std::array<std::uint64_t, 8> transition_powers{};
    transition_powers[0] = 1U;
    for (std::size_t index = 1U; index < transition_powers.size(); ++index) {
        transition_powers[index] = 15U * transition_powers[index - 1U];
    }
    const auto work_infinity = std::numeric_limits<std::uint64_t>::max();
    std::vector<std::vector<std::uint64_t>> work(
        left_size + 1U, std::vector<std::uint64_t>(right_size + 1U, work_infinity));
    std::vector<std::vector<char>> parent(left_size + 1U,
                                          std::vector<char>(right_size + 1U, 0));
    work[0][0] = 0U;
    for (std::size_t left = 0U; left <= left_size; ++left) {
        for (std::size_t right = 0U; right <= right_size; ++right) {
            if (work[left][right] == work_infinity) {
                continue;
            }
            const auto before = open_width(left, right);
            const auto update = [&](std::size_t next_left, std::size_t next_right, char step) {
                const auto after = open_width(next_left, next_right);
                if (after > target_width) {
                    return;
                }
                const auto candidate = work[left][right] +
                                       transition_powers[std::max(before, after)];
                if (candidate < work[next_left][next_right]) {
                    work[next_left][next_right] = candidate;
                    parent[next_left][next_right] = step;
                }
            };
            if (left < left_size) {
                update(left + 1U, right, 'L');
            }
            if (right < right_size) {
                update(left, right + 1U, 'R');
            }
        }
    }

    std::vector<std::pair<std::uint8_t, bool>> path;
    std::size_t left = left_size;
    std::size_t right = right_size;
    while (left != 0U || right != 0U) {
        if (parent[left][right] == 'L') {
            path.emplace_back(left_sequence[left - 1U], true);
            --left;
        } else if (parent[left][right] == 'R') {
            path.emplace_back(right_sequence[right - 1U], false);
            --right;
        } else {
            throw std::logic_error("failed to reconstruct contraction path");
        }
    }
    std::reverse(path.begin(), path.end());

    ShareEvaluationStats local_stats;
    local_stats.path_width = target_width;
    local_stats.peak_states = 1U;

    if (!dense_states) {
        using SparseStates = std::unordered_map<std::uint32_t, MatrixState>;
        const auto is_nonzero_state = [](const MatrixState& value) {
            return std::any_of(value.value.begin(), value.value.end(),
                               [](std::uint16_t entry) { return entry != 0U; });
        };
        SparseStates states;
        if (is_nonzero_state(initial)) {
            states.emplace(0U, initial);
        }
        const auto run_sparse = [&](auto&& self, std::size_t step,
                                    std::uint32_t current_mask,
                                    SparseStates current) -> MatrixState {
            for (; step < path.size(); ++step) {
                const auto [label, is_left] = path[step];
                const auto bit = std::uint32_t{1} << label;
                if ((current_mask & bit) == 0U &&
                    current.size() > sparse_state_cap_ / 15U) {

                    SparseStates other;
                    other.reserve(current.size() / 2U + 1U);
                    for (auto iterator = current.begin();
                         other.size() < current.size(); ) {
                        auto selected = iterator++;
                        other.insert(current.extract(selected));
                    }
                    auto first = self(self, step, current_mask, std::move(current));
                    const auto second = self(self, step, current_mask, std::move(other));
                    for (std::size_t index = 0U; index < first.value.size(); ++index) {
                        first.value[index] = static_cast<std::uint16_t>(
                            add_mod(first.value[index], second.value[index]));
                    }
                    return first;
                }
                SparseStates next;
                if ((current_mask & bit) == 0U) {
                    next.reserve(current.size() * 15U);
                    for (const auto& [key, matrices] : current) {
                        for (std::size_t generator = 0U; generator < 15U; ++generator) {
                            const auto next_key =
                                key | (static_cast<std::uint32_t>(generator) << (4U * label));
                            auto transformed = apply_generator(
                                matrices, primal_generators_[generator], is_left, tensor_square);
                            if (is_nonzero_state(transformed)) {
                                next.emplace(next_key, std::move(transformed));
                            }
                            ++local_stats.transitions;
                        }
                    }
                    current_mask |= bit;
                } else {
                    next.reserve(current.size() / 8U + 1U);
                    for (const auto& [key, matrices] : current) {
                        const auto generator = (key >> (4U * label)) & 15U;
                        const auto next_key = key & ~(std::uint32_t{15} << (4U * label));
                        auto transformed = apply_generator(
                            matrices, dual_generators_[generator], is_left, tensor_square);
                        if (!is_nonzero_state(transformed)) {
                            ++local_stats.transitions;
                            continue;
                        }
                        auto [iterator, inserted] = next.try_emplace(next_key);
                        if (inserted) {
                            iterator->second = std::move(transformed);
                        } else {
                            for (std::size_t index = 0U; index < transformed.value.size(); ++index) {
                                iterator->second.value[index] =
                                    static_cast<std::uint16_t>(add_mod(
                                        iterator->second.value[index], transformed.value[index]));
                            }
                        }
                        ++local_stats.transitions;
                    }
                    current_mask &= ~bit;
                    for (auto iterator = next.begin(); iterator != next.end();) {
                        if (!is_nonzero_state(iterator->second)) {
                            iterator = next.erase(iterator);
                        } else {
                            ++iterator;
                        }
                    }
                }
                current = std::move(next);
                local_stats.peak_states = std::max(local_stats.peak_states, current.size());
            }
            if (current_mask != 0U) {
                throw std::logic_error("sparse share contraction did not close all chord indices");
            }
            if (current.empty()) {
                return {};
            }
            if (current.size() != 1U || current.begin()->first != 0U) {
                throw std::logic_error("sparse share contraction did not reduce to one state");
            }
            return std::move(current.begin()->second);
        };
        return StateEvaluation{run_sparse(run_sparse, 0U, 0U, std::move(states)), local_stats};
    }

    constexpr std::size_t base = 15U;
    std::vector<MatrixState> states;
    std::vector<std::uint8_t> open_labels;
    std::array<std::size_t, 8> powers{};
    powers[0] = 1U;
    for (std::size_t exponent = 1U; exponent < powers.size(); ++exponent) {
        powers[exponent] = base * powers[exponent - 1U];
    }
    states.reserve(powers[local_stats.path_width]);
    states.push_back(initial);
    for (const auto [label, is_left] : path) {
        const auto found = std::find(open_labels.begin(), open_labels.end(), label);
        if (found == open_labels.end()) {
            const auto old_size = states.size();
            states.resize(old_size * base);
            for (std::size_t state = old_size; state-- > 0U;) {
                const auto source = states[state];
                for (std::size_t generator = 0U; generator < base; ++generator) {
                    states[state * base + generator] =
                        apply_generator(source, primal_generators_[generator],
                                        is_left, tensor_square);
                    ++local_stats.transitions;
                }
            }
            open_labels.push_back(label);
        } else {
            const auto position = static_cast<std::size_t>(found - open_labels.begin());
            const auto inner_count = powers[open_labels.size() - position - 1U];
            const auto outer_count = powers[position];
            for (std::size_t outer = 0U; outer < outer_count; ++outer) {
                for (std::size_t inner = 0U; inner < inner_count; ++inner) {
                    MatrixState target;
                    for (std::size_t generator = 0U; generator < base; ++generator) {
                        const auto source = (outer * base + generator) * inner_count + inner;
                        const auto transformed = apply_generator(
                            states[source], dual_generators_[generator], is_left, tensor_square);
                        for (std::size_t index = 0U; index < transformed.value.size(); ++index) {
                            target.value[index] = static_cast<std::uint16_t>(add_mod(
                                target.value[index], transformed.value[index]));
                        }
                        ++local_stats.transitions;
                    }
                    states[outer * inner_count + inner] = std::move(target);
                }
            }
            states.resize(outer_count * inner_count);
            open_labels.erase(found);
        }
        local_stats.peak_states = std::max(local_stats.peak_states, states.size());
    }
    if (states.size() != 1U || !open_labels.empty()) {
        throw std::logic_error("dense share contraction did not close all chord indices");
    }
    return StateEvaluation{std::move(states.front()), local_stats};
}

Block2x2 Sl4ShareDetector::mutation_odd(const TwoArcShare& share, int mutation_kind) {
    if (mutation_kind == 2) {
        return swap_odd(share);
    }
    const auto original = evaluate(share);
    const auto mutated = evaluate(mutate_two_arc_share(share, mutation_kind));
    Block2x2 result;
    for (std::size_t index = 0U; index < result.value.size(); ++index) {
        result.value[index] = subtract_mod(original.value[index], mutated.value[index]);
    }
    return result;
}

Block2x2 Sl4ShareDetector::tensor_square_block(const TwoArcShare& share) {
    return tensor_square_block_impl(share, true);
}

Block2x2 Sl4ShareDetector::tensor_square_block_sparse(const TwoArcShare& share) {
    return tensor_square_block_impl(share, false);
}

Block2x2 Sl4ShareDetector::tensor_square_block_impl(const TwoArcShare& share,
                                                     bool dense_states) {
    Block2x2 result;
    stats_ = {};
    for (std::size_t input = 0U; input < 2U; ++input) {
        MatrixState initial;
        std::copy(tensor_square_highest_[input].value.begin(),
                  tensor_square_highest_[input].value.end(), initial.value.begin());
        auto evaluated = contract(share, dense_states, initial, true);
        stats_.path_width = std::max(stats_.path_width, evaluated.stats.path_width);
        stats_.peak_states = std::max(stats_.peak_states, evaluated.stats.peak_states);
        stats_.transitions += evaluated.stats.transitions;
        const std::array<std::uint32_t, 2> sampled{
            evaluated.value.value[tensor_square_projection_rows_[0]],
            evaluated.value.value[tensor_square_projection_rows_[1]]};
        for (std::size_t output = 0U; output < 2U; ++output) {
            auto coordinate = std::uint32_t{0};
            for (std::size_t row = 0U; row < 2U; ++row) {
                coordinate = add_mod(
                    coordinate,
                    multiply_mod(tensor_square_projection_inverse_.value[output * 2U + row],
                                 sampled[row]));
            }
            result.value[input * 2U + output] = coordinate;
        }
        for (std::size_t position = 0U; position < matrix_size; ++position) {
            const auto expected = add_mod(
                multiply_mod(result.value[input * 2U],
                             tensor_square_highest_[0].value[position]),
                multiply_mod(result.value[input * 2U + 1U],
                             tensor_square_highest_[1].value[position]));
            if (evaluated.value.value[position] != expected) {
                throw std::logic_error("share operator left the repeated 64 tensor-square channel");
            }
        }
    }
    return result;
}

Block2x2 Sl4ShareDetector::tensor_square_mutation_odd(
    const TwoArcShare& share, int mutation_kind) {
    return tensor_square_mutation_odd_impl(share, mutation_kind, true);
}

Block2x2 Sl4ShareDetector::tensor_square_mutation_odd_sparse(
    const TwoArcShare& share, int mutation_kind) {
    return tensor_square_mutation_odd_impl(share, mutation_kind, false);
}

Block2x2 Sl4ShareDetector::tensor_square_mutation_odd_impl(
    const TwoArcShare& share, int mutation_kind, bool dense_states) {
    if (mutation_kind < 1 || mutation_kind > 3) {
        throw std::invalid_argument("share mutation kind must be in {1,2,3}");
    }
    const auto original = tensor_square_block_impl(share, dense_states);
    const auto first_stats = stats_;
    const auto mutated = tensor_square_block_impl(
        mutate_two_arc_share(share, mutation_kind), dense_states);
    stats_.path_width = std::max(first_stats.path_width, stats_.path_width);
    stats_.peak_states = std::max(first_stats.peak_states, stats_.peak_states);
    stats_.transitions += first_stats.transitions;
    return subtract_blocks(original, mutated);
}

std::vector<ChannelOddBlock> Sl4ShareDetector::channel_summary() const {
    std::vector<ChannelOddBlock> result;
    result.reserve(channels_.size());
    for (std::size_t index = 0U; index < channels_.size(); ++index) {
        result.push_back(ChannelOddBlock{
            index, channels_[index].highest_weight, channels_[index].positive.size(), {}});
    }
    return result;
}

std::vector<ChannelOddBlock> Sl4ShareDetector::channel_swap_odds(
    const TwoArcShare& share) {
    stats_ = {};
    std::vector<ChannelOddBlock> result;
    result.reserve(channels_.size());
    for (std::size_t channel_index = 0U; channel_index < channels_.size(); ++channel_index) {
        const auto& channel = channels_[channel_index];
        const auto multiplicity = channel.positive.size();
        std::vector<std::uint32_t> block(multiplicity * multiplicity, 0U);
        std::vector<std::uint32_t> mutated_block(multiplicity * multiplicity, 0U);
        for (std::size_t input = 0U; input < multiplicity; ++input) {
            MatrixState initial;
            std::copy(channel.positive[input].value.begin(), channel.positive[input].value.end(),
                      initial.value.begin());
            auto evaluated = contract(share, true, initial);
            auto evaluated_mutated = contract(mutate_two_arc_share(share, 2), true, initial);
            stats_.path_width = std::max(
                {stats_.path_width, evaluated.stats.path_width, evaluated_mutated.stats.path_width});
            stats_.peak_states = std::max(
                {stats_.peak_states, evaluated.stats.peak_states, evaluated_mutated.stats.peak_states});
            stats_.transitions += evaluated.stats.transitions + evaluated_mutated.stats.transitions;

            std::vector<std::uint32_t> sampled(multiplicity);
            std::vector<std::uint32_t> sampled_mutated(multiplicity);
            for (std::size_t row = 0U; row < multiplicity; ++row) {
                sampled[row] = evaluated.value.value[channel.projection_rows[row]];
                sampled_mutated[row] =
                    evaluated_mutated.value.value[channel.projection_rows[row]];
            }
            for (std::size_t output = 0U; output < multiplicity; ++output) {
                std::uint32_t coordinate = 0U;
                std::uint32_t mutated_coordinate = 0U;
                for (std::size_t row = 0U; row < multiplicity; ++row) {
                    coordinate = add_mod(
                        coordinate,
                        multiply_mod(channel.projection_inverse[output * multiplicity + row],
                                     sampled[row]));
                    mutated_coordinate = add_mod(
                        mutated_coordinate,
                        multiply_mod(channel.projection_inverse[output * multiplicity + row],
                                     sampled_mutated[row]));
                }
                block[output * multiplicity + input] = coordinate;
                mutated_block[output * multiplicity + input] = mutated_coordinate;
            }
            for (std::size_t position = 0U; position < matrix_size; ++position) {
                std::uint32_t expected = 0U;
                std::uint32_t expected_mutated = 0U;
                for (std::size_t output = 0U; output < multiplicity; ++output) {
                    expected = add_mod(
                        expected,
                        multiply_mod(block[output * multiplicity + input],
                                     channel.positive[output].value[position]));
                    expected_mutated = add_mod(
                        expected_mutated,
                        multiply_mod(mutated_block[output * multiplicity + input],
                                     channel.positive[output].value[position]));
                }
                if (evaluated.value.value[position] != expected ||
                    evaluated_mutated.value.value[position] != expected_mutated) {
                    throw std::logic_error("share operator left a highest-weight multiplicity space");
                }
            }
        }

        std::vector<std::uint32_t> odd(block.size());
        for (std::size_t index = 0U; index < odd.size(); ++index) {
            odd[index] = subtract_mod(block[index], mutated_block[index]);
        }
        result.push_back(ChannelOddBlock{channel_index, channel.highest_weight,
                                         multiplicity, std::move(odd)});
    }
    return result;
}

Block2x2 Sl4ShareDetector::swap_odd(const TwoArcShare& share) {
    const auto block = evaluate(share);
    const auto mutated = evaluate(mutate_two_arc_share(share, 2));
    return subtract_blocks(block, mutated);
}

ShareEvaluationStats Sl4ShareDetector::stats() const {
    return stats_;
}

void Sl4ShareDetector::set_sparse_state_cap(std::size_t cap) {
    if (cap < 15U) {
        throw std::invalid_argument("sparse state cap must be at least 15");
    }
    sparse_state_cap_ = cap;
}

void Sl4ShareDetector::self_test() const {
    for (std::size_t a = 0U; a < 4U; ++a) {
        for (std::size_t b = 0U; b < 4U; ++b) {
            for (std::size_t c = 0U; c < 4U; ++c) {
                for (std::size_t d = 0U; d < 4U; ++d) {
                    const auto& left = generators_[generator_index(a, b)];
                    const auto& right = generators_[generator_index(c, d)];
                    const auto& expected_left = generators_[generator_index(a, d)];
                    const auto& expected_right = generators_[generator_index(c, b)];
                    for (std::size_t row = 0U; row < dimension; ++row) {
                        for (std::size_t column = 0U; column < dimension; ++column) {
                            std::uint64_t lr = 0U;
                            std::uint64_t rl = 0U;
                            for (std::size_t inner = 0U; inner < dimension; ++inner) {
                                lr += static_cast<std::uint64_t>(left.value[row * dimension + inner]) *
                                      right.value[inner * dimension + column];
                                rl += static_cast<std::uint64_t>(right.value[row * dimension + inner]) *
                                      left.value[inner * dimension + column];
                            }
                            auto expected = std::uint32_t{0};
                            if (b == c) {
                                expected = add_mod(expected, expected_left.value[row * dimension + column]);
                            }
                            if (a == d) {
                                expected = subtract_mod(expected, expected_right.value[row * dimension + column]);
                            }
                            if (subtract_mod(static_cast<std::uint32_t>(lr % kModulus),
                                             static_cast<std::uint32_t>(rl % kModulus)) != expected) {
                                throw std::logic_error("[2,1] matrices violate gl4 commutation relations");
                            }
                        }
                    }
                }
            }
        }
    }
}

}  
