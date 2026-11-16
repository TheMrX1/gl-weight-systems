#include "chord/gl_weight.hpp"
#include <boost/json/src.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

namespace json=boost::json;
constexpr std::uint64_t P=1000000007;
std::uint64_t power(std::uint64_t a, unsigned long long b) {
    std::uint64_t r=1;
    while(b) {if(b&1)r=r*a%P;a=a*a%P;b>>=1;}
    return r;
}
std::uint64_t coefficient(const std::string& s) {
    auto slash=s.find('/');
    auto a=std::stoll(s.substr(0,slash));
    std::uint64_t r=(a%static_cast<long long>(P)+P)%P;
    if(slash!=std::string::npos)r=r*power(std::stoull(s.substr(slash+1)),P-2)%P;
    return r;
}
int main(int argc,char**argv) {
  try {
    if(argc!=2)throw std::runtime_error("Expected certificate JSON");
    chord::GlEvaluator e;
    auto c=e.casimir_values();
    std::array<unsigned,110> ranks{},xs{};
    for(unsigned lane=0;lane<110;++lane) {
        unsigned N=lane<100?lane/20+1:6;
        unsigned x=lane<100?lane%20:lane-99;
        ranks[lane]=N;xs[lane]=x;
        c[lane][0]=N;c[lane][1]=x;
        for(unsigned k=2;k<23;++k)c[lane][k]=c[lane][k-1]*N%P;
    }
    e.assign_casimirs(c);
    auto small=e.evaluate_chord_word({0,1,0,1});
    for(unsigned l=0;l<110;++l) {
        auto N=ranks[l],x=xs[l];
        auto expected=((std::uint64_t(x)*x+P-x)*N*N+std::uint64_t(x)*x)%P;
        if(small.value[l]!=expected)throw std::runtime_error("Crossing control failed");
    }
    std::ifstream in(argv[1]);std::ostringstream text;text<<in.rdbuf();
    auto data=json::parse(text.str());auto& entries=data.as_object().at("entries").as_array();
    std::array<std::uint64_t,110> totals{};
    auto begin=std::chrono::steady_clock::now();
    unsigned count=0;
    for(const auto& entry:entries) {
        auto& obj=entry.as_object();chord::Word w;
        for(const auto& label:obj.at("word").as_array())w.push_back(label.to_number<unsigned>());
        auto r=w;std::reverse(r.begin(),r.end());r=chord::relabel_by_first_occurrence(r);
        w=chord::relabel_by_first_occurrence(w);if(r<w)w=std::move(r);
        if(e.stats().cache_entries>4000000)e.clear_memo();
        auto result=e.evaluate_chord_word(w);
        auto alpha=coefficient(std::string(obj.at("coefficient").as_string()));
        json::array values;
        for(unsigned l=0;l<110;++l) {
            totals[l]=(totals[l]+result.value[l]*alpha)%P;
            values.push_back(result.value[l]);
        }
        std::cout<<json::serialize(json::object{{"id",obj.at("id")},{"values",values}})<<'\n';
        if(++count%100==0) {
            std::cout.flush();
            std::cerr<<"processed="<<count<<" seconds="
                     <<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()
                     <<" nodes="<<e.stats().recurrence_nodes<<'\n';
        }
    }
    json::array lanes;
    for(unsigned l=0;l<110;++l)lanes.push_back(json::object{{"N",ranks[l]},{"x",xs[l]},{"value",totals[l]}});
    std::cout<<json::serialize(json::object{{"summary",true},{"prime",P},{"entries",count},{"lanes",lanes}})<<'\n';
  } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
