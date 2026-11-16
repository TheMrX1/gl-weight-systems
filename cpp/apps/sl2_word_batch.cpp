#include "sl2_exact.hpp"
#include <boost/json/src.hpp>
#include <iostream>
#include <sstream>
namespace json=boost::json;
using namespace primitive_sl2;
int main() {
    Evaluator evaluator;std::string line;
    while(std::getline(std::cin,line)) {
      try {
        std::istringstream in(line);unsigned label;Word w;
        while(in>>label){if(label>255)throw std::runtime_error("Label exceeds 255");w.push_back(label);}
        if(!in.eof() || w.empty())throw std::runtime_error("Expected chord word");
        std::map<unsigned,unsigned> counts;for(auto a:w)++counts[a];
        if(counts.size()>16)throw std::runtime_error("At most 16 chords");
        for(auto [a,k]:counts)if(k!=2)throw std::runtime_error("Each label must occur twice");
        w=chord::relabel_by_first_occurrence(w);
        auto p=evaluator.primitive(w);json::array encoded;
        for(const auto& a:p)encoded.push_back(json::value(a.str()));
        std::cout<<json::serialize(json::object{{"n",counts.size()},{"primitive_scaled",encoded}})<<'\n'<<std::flush;
      }catch(const std::exception& e){std::cout<<json::serialize(json::object{{"error",e.what()}})<<'\n'<<std::flush;}
    }
}
