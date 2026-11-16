#include "chord/chord_diagram.hpp"
#include "chord/gl_weight.hpp"
#include "chord/sl4_share_detector.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        const chord::ShareBlocks blocks{
            {0U, 1U, 2U, 0U, 3U, 1U, 2U},
            {4U, 5U, 6U, 7U, 8U, 5U, 9U, 4U, 8U, 10U},
            {3U},
            {7U, 9U, 10U, 6U}};
        chord::Word original;
        for (const auto* part : {&blocks.a, &blocks.b, &blocks.c, &blocks.d}) {
            original.insert(original.end(), part->begin(), part->end());
        }
        const auto mutant = chord::mutate_share(blocks, 2);
        if (!chord::same_labeled_intersection_graph(original, mutant)) {
            throw std::logic_error("candidate does not preserve its labeled intersection graph");
        }
        if (chord::canonical_unoriented(original) == chord::canonical_unoriented(mutant)) {
            throw std::logic_error("candidate diagrams are isomorphic");
        }
        chord::Sl4ShareDetector detector;
        const chord::TwoArcShare four{blocks.a, blocks.c};
        const chord::TwoArcShare seven{{0U, 1U, 2U, 3U, 4U, 1U, 5U, 0U, 4U, 6U},
                                        {3U, 5U, 6U, 2U}};
        const auto score = chord::tensor_square_sewing_score(
            detector.tensor_square_mutation_odd(four, 2),
            detector.tensor_square_block(seven));
        if (score == 0U) {
            throw std::logic_error("specialized 4+7 signal vanished");
        }
        chord::GlEvaluator evaluator(0x4c414e444fULL ^ 0x9e3779b97f4a7c15ULL, true);
        const auto first = evaluator.evaluate_chord_word(original);
        evaluator.clear_memo();
        const auto second = evaluator.evaluate_chord_word(mutant);
        if (first == second) {
            throw std::logic_error("universal gl evaluations unexpectedly agree");
        }
        std::cout << "original=" << chord::word_to_string(original) << '\n'
                  << "mutant=" << chord::word_to_string(mutant) << '\n'
                  << "same_labeled_graph=1 nonisomorphic=1 sl4_score=" << score << '\n'
                  << "first=" << chord::residues_to_string(first) << '\n'
                  << "second=" << chord::residues_to_string(second) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
