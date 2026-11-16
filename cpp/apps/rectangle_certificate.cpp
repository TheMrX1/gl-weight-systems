#include "chord/gl_weight.hpp"
#include <boost/json/src.hpp>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <unordered_map>

namespace json=boost::json;
using chord::Residues;
using chord::Word;

class Probe {
    chord::GlEvaluator evaluator;
    std::unordered_map<std::string,Residues> words;
    std::size_t cap;
public:
    Probe(unsigned n,unsigned a,unsigned b,std::size_t cap,bool wide=false,
          std::array<unsigned,3> lane_n={}):cap(cap) {
        if (!a || !b || b>n || n>128 || a>10000)
            throw std::invalid_argument("Require 1<=height<=N<=128 and 1<=width<=10000");
        if(wide) evaluator.assign_moduli({2305843009213693951ULL,2305843009213693921ULL,2305843009213693907ULL});
        auto c=evaluator.casimir_values();
        for (unsigned lane=0;lane<3;++lane) {
            const unsigned rank=lane_n[lane]?lane_n[lane]:n;
            if(rank<b || rank>128)throw std::invalid_argument("Invalid lane rank");
            const auto p=evaluator.primes()[lane];
            c[lane][0]=rank;
            c[lane][1]=std::uint64_t(a)*b%p;
            for (unsigned k=2;k<23;++k) c[lane][k]=std::uint64_t((unsigned __int128)c[lane][k-1]*(rank+a-b)%p);
        }
        evaluator.assign_casimirs(c);
    }
    const auto& primes() const {return evaluator.primes();}
    auto stats() const {return evaluator.stats();}
    std::size_t word_count() const {return words.size();}
    Residues add(Residues x,const Residues& y) const {
        for (unsigned l=0;l<3;++l) x.value[l]=(x.value[l]+y.value[l])%primes()[l];
        return x;
    }
    Residues subtract(Residues x,const Residues& y) const {
        for (unsigned l=0;l<3;++l) x.value[l]=(x.value[l]+primes()[l]-y.value[l])%primes()[l];
        return x;
    }
    Residues multiply(Residues x,const Residues& y) const {
        for (unsigned l=0;l<3;++l) x.value[l]=std::uint64_t((unsigned __int128)x.value[l]*y.value[l]%primes()[l]);
        return x;
    }
    Residues moment(const Word& input) {
        if (input.empty()) {Residues one;one.value.fill(1);return one;}
        auto w=chord::relabel_by_first_occurrence(input);
        
        Word r(w.rbegin(),w.rend());r=chord::relabel_by_first_occurrence(r);
        if (r<w) w=std::move(r);
        std::string key(w.begin(),w.end());
        if (auto it=words.find(key);it!=words.end()) return it->second;
        auto result=evaluator.evaluate_chord_word(w);
        words.emplace(std::move(key),result);
        return result;
    }
    std::pair<Residues,Residues> evaluate(const Word& input,bool primitive) {
        if(input.size()%2 || input.size()>22) throw std::invalid_argument("Expected at most 11 paired chords");
        auto w=chord::relabel_by_first_occurrence(input);
        std::map<unsigned,unsigned> counts;
        for(auto c:w)++counts[c];
        for(auto [c,count]:counts) if(count!=2)throw std::invalid_argument("Each chord must occur twice");
        if (evaluator.stats().cache_entries>cap) evaluator.clear_memo();
        auto ordinary=moment(w);
        if (!primitive) return {ordinary,{}};
        if(w.empty())return {ordinary,{}};
        unsigned n=counts.size(),size=1u<<n;
        std::vector<Residues> moments(size),cumulants(size);
        moments[0].value.fill(1);
        for (unsigned mask=1;mask<size;++mask) {
            Word sub;
            for (auto c:w) if (mask>>c&1) sub.push_back(c);
            moments[mask]=moment(sub);
            auto k=moments[mask];
            const unsigned anchor=mask&-mask,rest=mask^anchor;
            if (rest) {
                unsigned s=(rest-1)&rest;
                while(true) {
                    const unsigned block=s|anchor;
                    k=subtract(k,multiply(cumulants[block],moments[mask^block]));
                    if(!s) break;
                    s=(s-1)&rest;
                }
            }
            cumulants[mask]=k;
        }
        return {ordinary,cumulants.back()};
    }
};

json::array residues(const Residues& v) {
    json::array a;for(auto x:v.value)a.push_back(x);return a;
}

int main(int argc,char** argv) {
    try {
        unsigned n=4,a=2,b=2,start=0,limit=5006,cap=500000;
        bool primitive=false,reverse=false,wide=false;
        std::array<unsigned,3> lane_n{};
        std::string file,word;
        for (int i=1;i<argc;++i) {
            std::string arg=argv[i];
            if(arg=="--primitive")primitive=true;
            else if(arg=="--reverse")reverse=true;
            else if(arg=="--wide-primes")wide=true;
            else {
                if(++i>=argc)throw std::invalid_argument("Missing option value");
                std::string value=argv[i];
                if(arg=="--N")n=std::stoul(value);
                else if(arg=="--width")a=std::stoul(value);
                else if(arg=="--height")b=std::stoul(value);
                else if(arg=="--start")start=std::stoul(value);
                else if(arg=="--limit")limit=std::stoul(value);
                else if(arg=="--memo-cap")cap=std::stoul(value);
                else if(arg=="--certificate")file=value;
                else if(arg=="--word")word=value;
                else if(arg=="--lane-Ns") {
                    std::istringstream in(value);
                    if(!(in>>lane_n[0]>>lane_n[1]>>lane_n[2]))throw std::invalid_argument("Expected three ranks");
                    std::string extra;if(in>>extra)throw std::invalid_argument("Expected exactly three ranks");
                }
                else throw std::invalid_argument("Unknown option "+arg);
            }
        }
        Probe probe(n,a,b,cap,wide,lane_n);
        const auto began=std::chrono::steady_clock::now();
        json::array entries;
        if(!word.empty()) {
            json::array w;std::istringstream in(word);unsigned x;
            while(in>>x){if(x>255)throw std::invalid_argument("Labels must fit in a byte");w.push_back(x);}
            if(!in.eof())throw std::invalid_argument("Invalid chord label");
            entries.push_back(json::object{{"id",0},{"word",w}});
            start=0;limit=1;
        } else {
            std::ifstream in(file);
            if(!in)throw std::invalid_argument("Provide --certificate FILE or --word WORD");
            std::ostringstream text;text<<in.rdbuf();
            entries=json::parse(text.str()).as_object().at("entries").as_array();
        }
        if(start>entries.size())throw std::invalid_argument("Start exceeds entry count");
        for(std::size_t i=start;i<entries.size() && i<std::size_t(start)+limit;++i) {
            const auto& entry=entries[i].as_object();Word w;
            for(const auto& x:entry.at("word").as_array()) {
                auto v=x.to_number<unsigned>();if(v>255)throw std::invalid_argument("Bad label");
                w.push_back(static_cast<std::uint8_t>(v));
            }
            if(reverse)std::reverse(w.begin(),w.end());
            auto [ordinary,k]=probe.evaluate(w,primitive);
            json::object result{{"id",entry.at("id")},{"N",n},{"width",a},{"height",b},
                                {"ordinary",residues(ordinary)}};
            if(wide)result["wide_primes"]=true;
            if(lane_n[0])result["lane_Ns"]=json::array{lane_n[0],lane_n[1],lane_n[2]};
            if(primitive)result["primitive"]=residues(k);
            std::cout<<json::serialize(result)<<'\n'<<std::flush;
            if((i-start+1)%20==0 || i+1==entries.size() || i+1==std::size_t(start)+limit) {
                auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
                std::cerr<<"processed="<<i-start+1<<" seconds="<<elapsed
                         <<" recurrence_nodes="<<probe.stats().recurrence_nodes
                         <<" memo="<<probe.stats().cache_entries<<" words="<<probe.word_count()<<'\n';
            }
        }
    } catch(const std::exception& e) {std::cerr<<"error: "<<e.what()<<'\n';return 1;}
}
