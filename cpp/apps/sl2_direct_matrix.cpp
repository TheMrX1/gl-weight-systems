#include <boost/json/src.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
using Big = boost::multiprecision::cpp_int;
using Wide = __int128;
namespace json = boost::json;

Wide moment(const std::vector<unsigned>& word, unsigned a) {
    auto dfs = [&](auto&& self, unsigned pos, unsigned state, unsigned seen,
                   unsigned Es, unsigned Fs, Wide product) -> Wide {
        if (pos == word.size()) return state == 0 ? product : 0;
        unsigned bit = 1u << word[pos];
        if (!(seen & bit)) {
            Wide sum = 0;
            if (state) sum += self(self, pos+1, state-1, seen|bit, Es|bit, Fs,
                                   product * (2*(a-state+1)));
            if (state < a) sum += self(self, pos+1, state+1, seen|bit, Es, Fs|bit,
                                       product * (2*(state+1)));
            int h = int(a)-2*int(state);
            if (h) sum += self(self, pos+1, state, seen|bit, Es, Fs, product*h);
            return sum;
        }
        if (Es & bit) return state < a ? self(self,pos+1,state+1,seen,Es,Fs,product*(state+1)) : 0;
        if (Fs & bit) return state ? self(self,pos+1,state-1,seen,Es,Fs,product*(a-state+1)) : 0;
        int h = int(a)-2*int(state);
        return h ? self(self,pos+1,state,seen,Es,Fs,product*h) : 0;
    };
    return dfs(dfs,0,0,0,0,0,1);
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("Expected graph JSON and highest weight");
        unsigned a = std::stoul(argv[2]);
        if (a > 8) throw std::runtime_error("Proved int128 bound restricts highest weight to <=8");
        std::ifstream f(argv[1]); std::ostringstream buffer; buffer << f.rdbuf();
        auto input = json::parse(buffer.str());
        for (const auto& item : input.as_array()) {
            const auto& g = item.as_object();
            unsigned n = g.at("n").to_number<unsigned>();
            if (n > 13) throw std::runtime_error("At most 13 chords");
            std::vector<unsigned> word;
            for (const auto& x : g.at("word").as_array()) word.push_back(x.to_number<unsigned>());
            
            std::reverse(word.begin(),word.end());
            std::vector<Wide> moments(1u<<n); moments[0] = 1;
            std::vector<Big> primitive(1u<<n);
            for (unsigned mask=1; mask<(1u<<n); ++mask) {
                std::vector<unsigned> sub;
                for (auto label : word) if (mask>>label & 1) sub.push_back(label);
                moments[mask] = moment(sub,a);
                primitive[mask] = moments[mask];
                unsigned anchor=mask&-mask,rest=mask^anchor;
                if (rest) {
                    unsigned s=(rest-1)&rest;
                    while (true) {
                        unsigned block=s|anchor;
                        primitive[mask] -= primitive[block]*moments[mask^block];
                        if (!s) break;
                        s=(s-1)&rest;
                    }
                }
            }
            Big expected=0,power=1; unsigned t=a*(a+2);
            for (const auto& c : g.at("primitive_scaled").as_array()) {
                expected += Big(std::string(c.as_string()))*power; power*=t;
            }
            if (primitive.back()!=expected) throw std::runtime_error("Direct matrix check disagrees");
            std::cout << json::serialize(json::object{{"id",g.at("id")},{"highest_weight",a},
                {"t",t},{"primitive_scaled_value",primitive.back().str()},{"passed",true}}) << '\n' << std::flush;
        }
    } catch (const std::exception& e) {std::cerr << e.what() << '\n';return 1;}
}
