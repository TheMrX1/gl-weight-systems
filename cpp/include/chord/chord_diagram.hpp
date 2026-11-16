#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace chord {

using Word = std::vector<std::uint8_t>;
using Cycle = std::vector<std::uint8_t>;
using Cycles = std::vector<Cycle>;

struct ShareBlocks {
    Word a;
    Word b;
    Word c;
    Word d;
};

struct RandomSharePair {
    ShareBlocks blocks;
    Word original;
    std::size_t share_chords{};
};

Word relabel_by_first_occurrence(const Word& word);
Word canonical_oriented(const Word& word);
Word canonical_unoriented(const Word& word);
Cycles word_to_involution_cycles(const Word& word);
std::vector<std::uint32_t> intersection_graph(const Word& word);
bool same_labeled_intersection_graph(const Word& lhs, const Word& rhs);
std::string word_to_string(const Word& word);

RandomSharePair random_share_diagram(std::size_t chord_count,
                                     std::size_t share_chords,
                                     std::mt19937_64& rng);
Word mutate_share(const ShareBlocks& blocks, int mutation_kind);

}  

