#include "sl2_exact.hpp"
#include <boost/json/src.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

namespace json=boost::json;
using namespace primitive_sl2;
json::array encode(const Poly& p) {json::array out;for(const auto& a:p)out.push_back(json::value(a.str()));return out;}
int main(int argc,char**argv) {
  try {
    Evaluator evaluator;
    if(argc==3 && std::string(argv[1])=="--word") {
        std::istringstream in(argv[2]);Word w;unsigned a;
        while(in>>a){if(a>255)throw std::runtime_error("Label too large");w.push_back(a);}
        if(!in.eof())throw std::runtime_error("Invalid input");
        std::map<unsigned,unsigned> counts;for(auto a:w)++counts[a];
        if(counts.empty())throw std::runtime_error("Expected a nonempty paired word");
        for(auto [a,k]:counts)if(k!=2)throw std::runtime_error("Each label must occur twice");
        if(counts.size()>16)throw std::runtime_error("At most 16 chords for this executable");
        w=chord::relabel_by_first_occurrence(w);
        auto p=evaluator.primitive(w),l=evaluator.primitive(leaf(w));
        std::cout<<json::serialize(json::object{{"n",w.size()/2},{"moment",encode(evaluator.weight(w))},
            {"primitive",encode(p)},{"leaf_primitive",encode(l)},
            {"normalization","H_D(t)=8^n*w_D(t/8); coefficients are exact integers"}})<<'\n';return 0;
    }
    if(argc!=2)throw std::runtime_error("Expected certificate JSON or --word WORD");
    std::ifstream in(argv[1]);std::ostringstream text;text<<in.rdbuf();
    auto data=json::parse(text.str());const auto& entries=data.as_object().at("entries").as_array();
    Poly sum{0},leaf_sum{0};unsigned count=0;
    auto begin=std::chrono::steady_clock::now();
    for(const auto& item:entries) {
        const auto& entry=item.as_object();Word w;
        for(const auto& a:entry.at("word").as_array())w.push_back(a.to_number<unsigned>());
        auto coeff=std::string(entry.at("coefficient").as_string());auto slash=coeff.find('/');
        Integer num(coeff.substr(0,slash)),den=slash==std::string::npos?Integer(1):Integer(coeff.substr(slash+1));
        if(2%den!=0)throw std::runtime_error("Certificate denominator does not divide 2");
        Integer scale=num*(2/den);
        auto p=evaluator.primitive(w),lp=evaluator.primitive(leaf(w));
        sum.resize(std::max(sum.size(),p.size()));leaf_sum.resize(std::max(leaf_sum.size(),lp.size()));
        for(unsigned k=0;k<p.size();++k)sum[k]+=scale*p[k];
        for(unsigned k=0;k<lp.size();++k)leaf_sum[k]+=scale*lp[k];
        auto expected=p;for(auto& a:expected)a*=-4;
        if(expected!=lp)throw std::runtime_error("Leaf identity failed at entry "+std::to_string(count));
        std::cout<<json::serialize(json::object{{"id",entry.at("id")},
                        {"moment_scaled",encode(evaluator.weight(w))},
                        {"primitive_scaled",encode(p)},{"leaf_primitive_scaled",encode(lp)}})<<'\n';
        if(++count%100==0) {
            std::cout.flush();std::cerr<<"checked="<<count<<" seconds="
            <<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()
            <<" moments="<<evaluator.moment_entries()<<" primitives="<<evaluator.primitive_entries()<<'\n';
        }
    }
    std::cout<<json::serialize(json::object{{"summary",true},{"entries",count},
        {"twice_scaled_primitive_sum_n9",encode(sum)},
        {"twice_scaled_leaf_sum_n10",encode(leaf_sum)},
        {"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()}})<<'\n';
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
