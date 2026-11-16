

#include "chord/symbolic_gl.hpp"
#include <chrono>
#include <iostream>
#include <map>
#include <sstream>

using Poly=chord::ExactPolynomial;
using Int=chord::ExactInt;
using Key=chord::MonomialKey;

chord::Word parse_word(const std::string& text){
    std::istringstream in(text);std::vector<int> labels;std::map<int,int> counts;int x;
    while(in>>x){if(x<0||x>255)throw std::invalid_argument("labels must be integers in 0..255");labels.push_back(x);++counts[x];}
    if(!in.eof()||labels.empty())throw std::invalid_argument("invalid or empty chord word");
    if(counts.size()>11)throw std::invalid_argument("at most 11 chords supported");
    std::map<int,int> ranks;int rank=0;
    for(auto [label,count]:counts){if(count!=2)throw std::invalid_argument("each label must occur exactly twice");ranks[label]=rank++;}
    chord::Word result;for(int label:labels)result.push_back(static_cast<std::uint8_t>(ranks.at(label)));return result;
}

Poly product(const Poly& a,const Poly& b){
    std::map<Key,Int> sums;
    for(const auto& x:a)for(const auto& y:b){
        Key key=0;
        for(unsigned k=0;k<23;++k){Key e=((x.monomial>>(5*k))&31)+((y.monomial>>(5*k))&31);
            if(e>31)throw std::overflow_error("monomial exponent overflow");key|=e<<(5*k);}
        Int c,total;
        if(__builtin_mul_overflow(x.coefficient,y.coefficient,&c)||__builtin_add_overflow(sums[key],c,&total))
            throw std::overflow_error("projection coefficient exceeds signed 128 bits");
        sums[key]=total;
    }
    Poly out;for(auto [m,c]:sums)if(c)out.push_back({m,c});return out;
}

std::uint64_t modular_product(std::uint64_t a,std::uint64_t b,std::uint64_t p){
    return static_cast<std::uint64_t>((static_cast<Key>(a)*b)%p);
}
std::uint64_t modular_power(std::uint64_t a,unsigned e,std::uint64_t p){
    std::uint64_t out=1;while(e){if(e&1)out=modular_product(out,a,p);a=modular_product(a,a,p);e>>=1;}return out;
}
chord::Residues evaluate_poly(const Poly& poly,const chord::GlEvaluator& evaluator){
    chord::Residues out;
    for(unsigned lane=0;lane<3;++lane){auto p=evaluator.primes()[lane];
        for(const auto& term:poly){Int remainder=term.coefficient%static_cast<Int>(p);if(remainder<0)remainder+=p;
            auto value=static_cast<std::uint64_t>(remainder);
            for(unsigned k=0;k<23;++k){unsigned e=static_cast<unsigned>((term.monomial>>(5*k))&31);
                if(e)value=modular_product(value,modular_power(evaluator.casimir_values()[lane][k],e,p),p);}
            out.value[lane]=(out.value[lane]+value)%p;
        }
    }
    return out;
}

int main(int argc,char** argv){try{
    std::string text;std::uint64_t seed=0x4c414e444fULL;bool show_primitive=false;
    for(int i=1;i<argc;++i){std::string a=argv[i];
        if(a=="--primitive-polynomial")show_primitive=true;
        else if(a=="--seed"){
            if(++i==argc)throw std::invalid_argument("missing seed");
            std::string s=argv[i];if(s.empty()||s.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("seed must be an unsigned integer");
            seed=std::stoull(s);
        }else if(a=="--help"){std::cout<<"Usage: gl_polynomial [--seed NUMBER] [--primitive-polynomial] [\"0 1 0 1\"]\n";return 0;}
        else if(text.empty())text=a;else throw std::invalid_argument("unexpected argument");
    }
    if(text.empty()){std::cerr<<"Enter chord word: ";std::getline(std::cin,text);}
    auto word=parse_word(text);unsigned n=static_cast<unsigned>(word.size()/2),end=1U<<n;
    auto start=std::chrono::steady_clock::now();
    Poly ordinary;
    {chord::SymbolicGlEvaluator symbolic;ordinary=symbolic.evaluate_chord_word(word);}
    std::cout<<"Wgl(D) = "<<chord::polynomial_string(ordinary)<<std::endl;
    std::vector<Poly> moments(end),primitive(end);moments[0]={{0,1}};moments.back()=ordinary;
    std::map<chord::Word,Poly> memo;
    for(unsigned mask=1;mask<end-1;++mask){chord::Word sub;for(auto label:word)if(mask&(1U<<label))sub.push_back(label);
        auto key=chord::canonical_oriented(sub);auto it=memo.find(key);
        if(it==memo.end()){chord::SymbolicGlEvaluator symbolic;it=memo.emplace(key,symbolic.evaluate_chord_word(key)).first;}
        moments[mask]=it->second;
    }
    
    for(unsigned mask=1;mask<end;++mask){auto value=moments[mask];unsigned anchor=mask&(~mask+1U);
        for(unsigned t=(mask-1)&mask;t;t=(t-1)&mask)if(t&anchor)
            value=chord::polynomial_sum(value,product(primitive[t],moments[mask^t]),-1);
        primitive[mask]=std::move(value);
    }
    chord::GlEvaluator numeric(seed^0x9e3779b97f4a7c15ULL,false);
    auto ordinary_values=evaluate_poly(ordinary,numeric);

    auto direct=numeric.evaluate_chord_word(word);
    if(!(direct==ordinary_values))throw std::runtime_error("numeric recurrence disagrees with exact polynomial");
    auto primitive_values=evaluate_poly(primitive.back(),numeric);
    
    std::vector<chord::Residues> m(end),k(end);m[0].value={1,1,1};
    for(unsigned mask=1;mask<end;++mask){m[mask]=evaluate_poly(moments[mask],numeric);k[mask]=m[mask];unsigned anchor=mask&(~mask+1U);
        for(unsigned t=(mask-1)&mask;t;t=(t-1)&mask)if(t&anchor)for(unsigned lane=0;lane<3;++lane){auto p=numeric.primes()[lane];
            k[mask].value[lane]=(k[mask].value[lane]+p-modular_product(k[t].value[lane],m[mask^t].value[lane],p))%p;}
        if(!(k[mask]==evaluate_poly(primitive[mask],numeric)))throw std::runtime_error("exact and modular projections disagree");
    }
    if(show_primitive)std::cout<<"Wgl(pi(D)) = "<<chord::polynomial_string(primitive.back())<<'\n';
    std::cout<<"primes: "<<chord::residues_to_string(chord::Residues{numeric.primes()})<<'\n';
    std::cout<<"Wgl(D): ";for(auto x:ordinary_values.value)std::cout<<x<<' ';std::cout<<'\n';
    std::cout<<"Wgl(pi(D)): ";for(auto x:primitive_values.value)std::cout<<x<<' ';std::cout<<'\n';

    std::cerr<<"CHECK PASS: direct numeric Wgl(D) equals polynomial substitution; all exact/modular projections agree\nseed="<<seed<<'\n';
    for(unsigned lane=0;lane<3;++lane){std::cerr<<"p="<<numeric.primes()[lane]<<" C0=N,C1,...,C22:";
        for(auto x:numeric.casimir_values()[lane])std::cerr<<' '<<x;std::cerr<<'\n';}
    std::cerr<<"seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}}
