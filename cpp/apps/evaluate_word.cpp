#include "chord/chord_diagram.hpp"
#include "chord/gl_weight.hpp"

#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

chord::Word parse_word(const std::string& text) {
    std::string normalized = text;
    for (auto& character : normalized) {
        if (character == ',') {
            character = ' ';
        }
    }
    std::istringstream stream(normalized);
    chord::Word word;
    unsigned int label = 0U;
    while (stream >> label) {
        if (label > 255U) {
            throw std::invalid_argument("labels must fit into one byte");
        }
        word.push_back(static_cast<std::uint8_t>(label));
    }
    if (word.empty()) {
        throw std::invalid_argument("empty chord word");
    }
    return word;
}

}  

int main(int argc, char** argv) {
    try {
        bool traceless = false;
        std::string word_text;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--specialization") {
                if (++index >= argc) {
                    throw std::invalid_argument("missing value for --specialization");
                }
                const std::string value = argv[index];
                if (value == "sl") {
                    traceless = true;
                } else if (value == "gl") {
                    traceless = false;
                } else {
                    throw std::invalid_argument("specialization must be gl or sl");
                }
            } else if (word_text.empty()) {
                word_text = argument;
            } else {
                throw std::invalid_argument("only one chord word may be supplied");
            }
        }
        if (word_text.empty()) {
            std::cerr << "Usage: evaluate_word [--specialization gl|sl] '0 1 0 1'\n";
            return 1;
        }
        constexpr std::uint64_t search_seed = 0x4c414e444fULL;
        chord::GlEvaluator evaluator(search_seed ^ 0x9e3779b97f4a7c15ULL, traceless);
        const auto word = parse_word(word_text);
        const auto value = evaluator.evaluate_chord_word(word);
        std::cout << "{\"specialization\":\"" << (traceless ? "sl" : "gl")
                  << "\",\"values\":[";
        for (std::size_t lane = 0; lane < chord::kResidueLanes; ++lane) {
            if (lane != 0U) {
                std::cout << ',';
            }
            std::cout << value.value[lane];
        }
        std::cout << "],\"primes\":[";
        for (std::size_t lane = 0; lane < chord::kResidueLanes; ++lane) {
            if (lane != 0U) {
                std::cout << ',';
            }
            std::cout << evaluator.primes()[lane];
        }
        std::cout << "],\"casimirs\":[";
        for (std::size_t lane = 0; lane < chord::kResidueLanes; ++lane) {
            if (lane != 0U) {
                std::cout << ',';
            }
            std::cout << '[';
            for (std::size_t index = 0; index < 23U; ++index) {
                if (index != 0U) {
                    std::cout << ',';
                }
                std::cout << evaluator.casimir_values()[lane][index];
            }
            std::cout << ']';
        }
        std::cout << "]}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
