#include "chord/chord_diagram.hpp"
#include "chord/gl_weight.hpp"
#include "chord/sl4_share_detector.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

bool nonzero(const chord::Block2x2& block) {
    return std::any_of(block.value.begin(), block.value.end(),
                       [](auto value) { return value != 0U; });
}

bool connected(const chord::TwoArcShare& share) {
    chord::Word word = share.first;
    word.insert(word.end(), share.second.begin(), share.second.end());
    const auto adjacency = chord::intersection_graph(word);
    if (adjacency.empty()) {
        return false;
    }
    std::uint32_t reached = 1U;
    for (;;) {
        auto expanded = reached;
        for (std::size_t vertex = 0U; vertex < adjacency.size(); ++vertex) {
            if ((reached & (std::uint32_t{1} << vertex)) != 0U) {
                expanded |= adjacency[vertex];
            }
        }
        if (expanded == reached) {
            break;
        }
        reached = expanded;
    }
    return std::popcount(reached) == adjacency.size();
}

chord::Word shifted(const chord::Word& word, std::size_t offset) {
    chord::Word result = word;
    for (auto& label : result) {
        label = static_cast<std::uint8_t>(label + offset);
    }
    return result;
}

chord::ShareBlocks glue_blocks(const chord::TwoArcShare& share,
                               const chord::TwoArcShare& complement,
                               std::size_t share_size) {
    return chord::ShareBlocks{share.first, shifted(complement.first, share_size),
                              share.second, shifted(complement.second, share_size)};
}

chord::Word original_word(const chord::ShareBlocks& blocks) {
    chord::Word result;
    for (const auto* block : {&blocks.a, &blocks.b, &blocks.c, &blocks.d}) {
        result.insert(result.end(), block->begin(), block->end());
    }
    return result;
}

}  

int main() {
    try {
        std::cout << std::unitbuf;
        auto shares_five = chord::enumerate_two_arc_shares(5U);
        auto shares_six = chord::enumerate_two_arc_shares(6U);
        std::stable_sort(shares_six.begin(), shares_six.end(), [](const auto& left, const auto& right) {
            return chord::two_arc_share_path_width(left, true) <
                   chord::two_arc_share_path_width(right, true);
        });

        chord::Sl4ShareDetector detector;
        struct SensitiveShare {
            chord::TwoArcShare share;
            chord::Block2x2 odd;
        };
        std::vector<SensitiveShare> sensitive;
        for (const auto& share : shares_five) {
            if (!connected(share)) {
                continue;
            }
            const auto odd = detector.tensor_square_mutation_odd(share, 2);
            if (nonzero(odd)) {
                sensitive.push_back(SensitiveShare{share, odd});
                if (sensitive.size() == 32U) {
                    break;
                }
            }
        }
        std::cout << "SENSITIVE_FIVE count=" << sensitive.size() << '\n';
        if (sensitive.empty()) {
            return 2;
        }

        const auto shares_four = chord::enumerate_two_arc_shares(4U);
        std::vector<chord::Block2x2> sensitive_four;
        for (const auto& share : shares_four) {
            const auto odd = detector.tensor_square_mutation_odd(share, 2);
            if (nonzero(odd)) {
                sensitive_four.push_back(odd);
            }
        }
        std::vector<chord::Block2x2> complement_five;
        for (const auto& share : shares_five) {
            if (connected(share)) {
                complement_five.push_back(detector.tensor_square_block(share));
                if (complement_five.size() == 256U) {
                    break;
                }
            }
        }
        std::size_t low_degree_checks = 0U;
        for (const auto& difference : sensitive_four) {
            for (const auto& complement : complement_five) {
                ++low_degree_checks;
                if (chord::tensor_square_sewing_score(difference, complement) != 0U) {
                    throw std::logic_error("corrected sewing failed a degree 4+5 regression");
                }
            }
        }
        for (const auto& item : sensitive) {
            for (const auto& complement : complement_five) {
                ++low_degree_checks;
                if (chord::tensor_square_sewing_score(item.odd, complement) != 0U) {
                    throw std::logic_error("corrected sewing failed a degree 5+5 regression");
                }
            }
        }
        std::cout << "LOW_DEGREE_CHECKS count=" << low_degree_checks << " failures=0\n";

        chord::GlEvaluator evaluator(0x74656e736f725f36ULL, true);
        evaluator.set_recurrence_budget(100'000U);
        std::size_t predicted = 0U;
        std::size_t verified = 0U;
        std::size_t budget_skips = 0U;
        std::size_t skipped_width_six = 0U;
        for (std::size_t complement_index = 0U;
             complement_index < shares_six.size();
             ++complement_index) {
            if (complement_index != 0U && complement_index % 1'000U == 0U) {
                std::cout << "PROGRESS complement_index=" << complement_index
                          << " predicted=" << predicted << " verified=" << verified << '\n';
            }
            if (!connected(shares_six[complement_index])) {
                continue;
            }
            if (chord::two_arc_share_path_width(shares_six[complement_index], true) > 5U) {
                ++skipped_width_six;
                continue;
            }
            const auto complement_block = detector.tensor_square_block(shares_six[complement_index]);
            for (std::size_t share_index = 0U; share_index < sensitive.size(); ++share_index) {
                if (chord::tensor_square_sewing_score(
                        sensitive[share_index].odd, complement_block) == 0U) {
                    continue;
                }
                ++predicted;
                if (predicted <= 5U) {
                    std::cout << "PREDICTED share_index=" << share_index
                              << " complement_index=" << complement_index
                              << " share_first=" << chord::word_to_string(sensitive[share_index].share.first)
                              << " share_second=" << chord::word_to_string(sensitive[share_index].share.second)
                              << " complement_first=" << chord::word_to_string(shares_six[complement_index].first)
                              << " complement_second=" << chord::word_to_string(shares_six[complement_index].second)
                              << '\n';
                }
                const auto blocks = glue_blocks(
                    sensitive[share_index].share, shares_six[complement_index], 5U);
                const auto original = original_word(blocks);
                const auto mutant = chord::mutate_share(blocks, 2);
                if (!chord::same_labeled_intersection_graph(original, mutant) ||
                    chord::canonical_oriented(original) == chord::canonical_oriented(mutant)) {
                    continue;
                }
                try {
                    evaluator.clear_memo();
                    const auto original_value = evaluator.evaluate_chord_word(original);
                    evaluator.clear_memo();
                    const auto mutant_value = evaluator.evaluate_chord_word(mutant);
                    ++verified;
                    if (original_value != mutant_value) {
                        std::cout << "FOUND share_index=" << share_index
                                  << " complement_index=" << complement_index
                                  << " original=" << chord::word_to_string(original)
                                  << " mutant=" << chord::word_to_string(mutant)
                                  << " value_original=" << chord::residues_to_string(original_value)
                                  << " value_mutant=" << chord::residues_to_string(mutant_value)
                                  << '\n';
                        return 0;
                    }
                } catch (const chord::EvaluationBudgetExceeded&) {
                    ++budget_skips;
                }
                if (verified >= 200U || budget_skips >= 200U) {
                    std::cout << "NO_WITNESS predicted=" << predicted
                              << " verified=" << verified
                              << " budget_skips=" << budget_skips << '\n';
                    return 3;
                }
            }
        }
        std::cout << "NO_WITNESS predicted=" << predicted << " verified=" << verified
                  << " budget_skips=" << budget_skips
                  << " skipped_width_six=" << skipped_width_six << '\n';
        return 3;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
