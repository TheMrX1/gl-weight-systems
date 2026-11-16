#include "chord/gl_weight.hpp"
#include "chord/symbolic_gl.hpp"
#include <bit>
#include <chrono>
#include <iostream>
#include <sstream>
#include <map>

chord::Word parse(const std::string& text) {
    std::istringstream in(text); int x; std::vector<int> raw; std::map<int,int> counts;
    while(in>>x){if(x<0||x>255)throw std::invalid_argument("label outside 0..255");raw.push_back(x);++counts[x];}
    if(!in.eof()||raw.empty())throw std::invalid_argument("invalid or empty word");
    if(counts.size()>11)throw std::invalid_argument("at most 11 chords");
    std::map<int,int> rank;int k=0;
    for(auto [label,count]:counts){if(count!=2)throw std::invalid_argument("each label must occur twice");rank[label]=k++;}
    chord::Word w;for(int label:raw)w.push_back(static_cast<std::uint8_t>(rank.at(label)));return w;
}
chord::Word subset(const chord::Word& w,unsigned mask){chord::Word out;for(auto x:w)if(mask&(1U<<x))out.push_back(x);return out;}
void residues(const chord::Residues& r){std::cout<<'[';for(unsigned j=0;j<3;++j){if(j)std::cout<<',';std::cout<<r.value[j];}std::cout<<']';}
int main(int argc,char** argv){try{
    bool symbolic=false,sl=false,proper=false;std::string text;std::uint64_t seed=0x4c414e444fULL;
    for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--symbolic")symbolic=true;else if(a=="--sl")sl=true;else if(a=="--proper")proper=true;else if(a=="--seed"){if(++i==argc)throw std::invalid_argument("missing seed");seed=std::stoull(argv[i]);}else if(text.empty())text=a;else throw std::invalid_argument("unexpected argument");}
    auto w=parse(text);unsigned n=static_cast<unsigned>(w.size()/2),end=1U<<n;
    auto start=std::chrono::steady_clock::now();
    if(symbolic){
        if(sl)throw std::invalid_argument("symbolic mode is universal gl");
        std::map<chord::Word,chord::ExactPolynomial> memo;
        for(unsigned mask=1;mask<end-(proper?1U:0U);++mask){auto key=chord::canonical_oriented(subset(w,mask));auto it=memo.find(key);
            if(it==memo.end()){chord::SymbolicGlEvaluator e;it=memo.emplace(key,e.evaluate_chord_word(key)).first;}
            std::cout<<"{\"mask\":"<<mask<<",\"terms\":"<<chord::polynomial_json(it->second)<<"}\n";
            if(mask%128==0)std::cerr<<"symbolic masks="<<mask<<" unique="<<memo.size()<<'\n';
        }
        std::cerr<<"symbolic unique="<<memo.size()<<'\n';
    }else{
        chord::GlEvaluator e(seed^0x9e3779b97f4a7c15ULL,sl);
        std::vector<chord::Residues> m(end),k(end);m[0].value={1,1,1};
        std::map<chord::Word,chord::Residues> memo;

        for(unsigned mask=1;mask<end;++mask){auto key=chord::canonical_oriented(subset(w,mask));auto it=memo.find(key);
            if(it==memo.end()){
                if(mask==end-1||e.stats().cache_entries>1500000)e.clear_memo();
                it=memo.emplace(key,e.evaluate_chord_word(key)).first;
            }m[mask]=it->second;
            if(mask%256==0||mask==end-1)std::cerr<<"numeric masks="<<mask<<" unique="<<memo.size()<<'\n';
        }
        std::uint64_t operations=0;
        for(unsigned mask=1;mask<end;++mask){k[mask]=m[mask];unsigned anchor=mask&(~mask+1U);
            for(unsigned t=(mask-1)&mask;t;t=(t-1)&mask)if(t&anchor){++operations;for(unsigned j=0;j<3;++j){auto p=e.primes()[j];auto product=static_cast<std::uint64_t>((static_cast<unsigned __int128>(k[t].value[j])*m[mask^t].value[j])%p);k[mask].value[j]=(k[mask].value[j]+p-product)%p;}}
        }
        std::cout<<"{\"n\":"<<n<<",\"seed\":"<<seed<<",\"specialization\":\""<<(sl?"sl":"gl")<<"\",\"primes\":[";
        for(unsigned j=0;j<3;++j){if(j)std::cout<<',';std::cout<<e.primes()[j];}std::cout<<"],\"casimirs\":[";
        for(unsigned j=0;j<3;++j){if(j)std::cout<<',';std::cout<<'[';for(unsigned c=0;c<23;++c){if(c)std::cout<<',';std::cout<<e.casimir_values()[j][c];}std::cout<<']';}
        std::cout<<"],\"moments\":[";for(unsigned mask=0;mask<end;++mask){if(mask)std::cout<<',';residues(m[mask]);}
        std::cout<<"],\"primitive\":[";for(unsigned mask=0;mask<end;++mask){if(mask)std::cout<<',';residues(k[mask]);}
        std::cout<<"],\"operations\":"<<operations<<",\"unique_words\":"<<memo.size()<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
    }
    return 0;
}catch(const std::exception& e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}}
