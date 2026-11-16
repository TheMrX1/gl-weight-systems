#include "chord/chord_diagram.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace chord {
namespace {

Word rotate_word(const Word& word, std::size_t shift) {
    Word result;
    result.reserve(word.size());
    for (std::size_t index = 0; index < word.size(); ++index) {
        result.push_back(word[(index + shift) % word.size()]);
    }
    return result;
}

Word reverse_word(Word word) {
    std::reverse(word.begin(), word.end());
    return word;
}

Word concatenate(const ShareBlocks& blocks, int mutation_kind) {
    Word result;
    result.reserve(blocks.a.size() + blocks.b.size() + blocks.c.size() + blocks.d.size());

    const auto append = [&result](const Word& block) {
        result.insert(result.end(), block.begin(), block.end());
    };

    switch (mutation_kind) {
        case 0:
            append(blocks.a);
            append(blocks.b);
            append(blocks.c);
            append(blocks.d);
            break;
        case 1:
            append(reverse_word(blocks.a));
            append(blocks.b);
            append(reverse_word(blocks.c));
            append(blocks.d);
            break;
        case 2:
            append(blocks.c);
            append(blocks.b);
            append(blocks.a);
            append(blocks.d);
            break;
        case 3:
            append(reverse_word(blocks.c));
            append(blocks.b);
            append(reverse_word(blocks.a));
            append(blocks.d);
            break;
        default:
            throw std::invalid_argument("mutation kind must be in {0,1,2,3}");
    }
    return result;
}

std::vector<std::uint32_t> labeled_intersection_graph(const Word& word) {
    if (word.size() % 2U != 0U) {
        throw std::invalid_argument("a chord word must have even length");
    }
    const std::size_t chord_count = word.size() / 2U;
    std::vector<std::array<std::size_t, 2>> positions(chord_count);
    std::vector<std::size_t> counts(chord_count, 0U);
    for (std::size_t position = 0; position < word.size(); ++position) {
        const auto label = static_cast<std::size_t>(word[position]);
        if (label >= chord_count || counts[label] >= 2U) {
            throw std::invalid_argument("chord labels must be 0..n-1 and occur exactly twice");
        }
        positions[label][counts[label]++] = position;
    }
    if (std::any_of(counts.begin(), counts.end(), [](std::size_t count) { return count != 2U; })) {
        throw std::invalid_argument("every chord label must occur exactly twice");
    }

    std::vector<std::uint32_t> adjacency(chord_count, 0U);
    for (std::size_t left = 0; left < chord_count; ++left) {
        const auto [a, b] = positions[left];
        for (std::size_t right = left + 1U; right < chord_count; ++right) {
            const auto [c, d] = positions[right];
            const bool crosses = (a < c && c < b && b < d) || (c < a && a < d && d < b);
            if (crosses) {
                adjacency[left] |= std::uint32_t{1} << right;
                adjacency[right] |= std::uint32_t{1} << left;
            }
        }
    }
    return adjacency;
}

Word shuffled_double_word(std::size_t first_label, std::size_t chord_count,
                          std::mt19937_64& rng) {
    Word result;
    result.reserve(2U * chord_count);
    for (std::size_t label = 0; label < chord_count; ++label) {
        const auto value = static_cast<std::uint8_t>(first_label + label);
        result.push_back(value);
        result.push_back(value);
    }
    std::shuffle(result.begin(), result.end(), rng);
    return result;
}

}  

Word relabel_by_first_occurrence(const Word& word) {
    std::array<int, 256> mapping{};
    mapping.fill(-1);
    std::uint16_t next_label = 0;
    Word result;
    result.reserve(word.size());
    for (const auto label : word) {
        if (mapping[label] < 0) {
            if (next_label > std::numeric_limits<std::uint8_t>::max()) {
                throw std::overflow_error("too many chord labels");
            }
            mapping[label] = static_cast<int>(next_label++);
        }
        result.push_back(static_cast<std::uint8_t>(mapping[label]));
    }
    return result;
}

Word canonical_oriented(const Word& word) {
    if (word.empty()) {
        return {};
    }
    Word best;
    bool initialized = false;
    for (std::size_t shift = 0; shift < word.size(); ++shift) {
        auto candidate = relabel_by_first_occurrence(rotate_word(word, shift));
        if (!initialized || candidate < best) {
            best = std::move(candidate);
            initialized = true;
        }
    }
    return best;
}

Word canonical_unoriented(const Word& word) {
    auto best = canonical_oriented(word);
    auto reversed = word;
    std::reverse(reversed.begin(), reversed.end());
    auto reflected = canonical_oriented(reversed);
    return std::min(best, reflected);
}

Cycles word_to_involution_cycles(const Word& word) {
    const auto normalized = relabel_by_first_occurrence(word);
    if (normalized.size() % 2U != 0U) {
        throw std::invalid_argument("a chord word must have even length");
    }
    const std::size_t chord_count = normalized.size() / 2U;
    std::vector<Cycle> positions(chord_count);
    for (std::size_t position = 0; position < normalized.size(); ++position) {
        const auto label = static_cast<std::size_t>(normalized[position]);
        if (label >= chord_count) {
            throw std::invalid_argument("invalid number of chord labels");
        }
        positions[label].push_back(static_cast<std::uint8_t>(position));
    }
    if (std::any_of(positions.begin(), positions.end(),
                    [](const Cycle& cycle) { return cycle.size() != 2U; })) {
        throw std::invalid_argument("every chord label must occur exactly twice");
    }
    return positions;
}

std::vector<std::uint32_t> intersection_graph(const Word& word) {
    return labeled_intersection_graph(relabel_by_first_occurrence(word));
}

bool same_labeled_intersection_graph(const Word& lhs, const Word& rhs) {
    return labeled_intersection_graph(lhs) == labeled_intersection_graph(rhs);
}

std::string word_to_string(const Word& word) {
    std::ostringstream stream;
    for (std::size_t index = 0; index < word.size(); ++index) {
        if (index != 0U) {
            stream << ' ';
        }
        stream << static_cast<unsigned int>(word[index]);
    }
    return stream.str();
}

RandomSharePair random_share_diagram(std::size_t chord_count,
                                     std::size_t share_chords,
                                     std::mt19937_64& rng) {
    if (chord_count < 2U || share_chords == 0U || share_chords >= chord_count) {
        throw std::invalid_argument("share size must be between 1 and n-1");
    }
    if (chord_count > 127U) {
        throw std::invalid_argument("this implementation supports at most 127 chords");
    }

    auto share = shuffled_double_word(0U, share_chords, rng);
    auto complement = shuffled_double_word(share_chords, chord_count - share_chords, rng);
    std::uniform_int_distribution<std::size_t> share_cut(1U, share.size() - 1U);
    std::uniform_int_distribution<std::size_t> complement_cut(1U, complement.size() - 1U);
    const auto a_length = share_cut(rng);
    const auto b_length = complement_cut(rng);

    ShareBlocks blocks;
    blocks.a.assign(share.begin(), share.begin() + static_cast<std::ptrdiff_t>(a_length));
    blocks.c.assign(share.begin() + static_cast<std::ptrdiff_t>(a_length), share.end());
    blocks.b.assign(complement.begin(), complement.begin() + static_cast<std::ptrdiff_t>(b_length));
    blocks.d.assign(complement.begin() + static_cast<std::ptrdiff_t>(b_length), complement.end());

    return RandomSharePair{blocks, concatenate(blocks, 0), share_chords};
}

Word mutate_share(const ShareBlocks& blocks, int mutation_kind) {
    if (mutation_kind < 1 || mutation_kind > 3) {
        throw std::invalid_argument("nontrivial mutation kind must be in {1,2,3}");
    }
    return concatenate(blocks, mutation_kind);
}

}  

