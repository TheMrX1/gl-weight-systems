#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Poly = std::array<std::int64_t, 10>;

Poly chromatic(const std::string& text) {
    if (text.empty()) throw std::runtime_error("Empty graph6 input");
    const int n = static_cast<unsigned char>(text[0]) - 63;
    if (n < 0 || n > 9) throw std::runtime_error("Supported order is 0..9");
    if (text.size() != 1 + (n*(n-1)/2 + 5)/6)
        throw std::runtime_error("Wrong graph6 length");
    for (unsigned char c : text)
        if (c < 63 || c > 126) throw std::runtime_error("Invalid graph6 byte");
    std::array<unsigned, 9> adjacency{};
    int bit = 0;
    for (int j=1;j<n;++j) for (int i=0;i<j;++i,++bit) {
        if (((text[1+bit/6]-63) >> (5-bit%6)) & 1) {
            adjacency[i] |= 1u << j;
            adjacency[j] |= 1u << i;
        }
    }
    const unsigned size = 1u << n;
    std::vector<bool> independent(size, false);
    independent[0] = true;
    for (unsigned mask=1;mask<size;++mask) {
        const unsigned v = __builtin_ctz(mask);
        const unsigned rest = mask & (mask-1);
        independent[mask] = independent[rest] && !(adjacency[v] & rest);
    }
    
    std::vector<Poly> dp(size);
    dp[0][0] = 1;
    for (unsigned mask=1;mask<size;++mask) {
        const unsigned anchor = mask & -mask;
        const unsigned rest = mask ^ anchor;
        unsigned sub = rest;
        while (true) {
            const unsigned block = sub | anchor;
            if (independent[block]) {
                const auto& tail = dp[mask ^ block];
                for (int k=1;k<=n;++k) dp[mask][k] += tail[k-1];
            }
            if (!sub) break;
            sub = (sub-1) & rest;
        }
    }
    
    Poly result{}, falling{};
    falling[0] = 1;
    for (int k=0;k<=n;++k) {
        for (int j=0;j<=k;++j) result[j] += dp.back()[k]*falling[j];
        if (k==n) break;
        Poly next{};
        for (int j=0;j<=k;++j) {
            next[j] -= k*falling[j];
            next[j+1] += falling[j];
        }
        falling = next;
    }
    return result;
}

int main() {
    try {
        std::string line;
        while (std::getline(std::cin,line)) {
            if (!line.empty() && line.back()=='\r') line.pop_back();
            auto p = chromatic(line);
            const int n = static_cast<unsigned char>(line[0])-63;
            std::cout << line;
            for (int k=0;k<=n;++k) std::cout << ' ' << p[k];
            std::cout << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
