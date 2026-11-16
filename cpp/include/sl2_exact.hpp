#pragma once
#include "chord/chord_diagram.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>

namespace primitive_sl2 {
using Integer=boost::multiprecision::cpp_int;
using Poly=std::vector<Integer>;
using chord::Word;

inline void trim(Poly& p) {while(p.size()>1 && p.back()==0)p.pop_back();}
inline Poly add(Poly a,const Poly& b,int sign=1) {
    a.resize(std::max(a.size(),b.size()));
    for(std::size_t k=0;k<b.size();++k)a[k]+=sign*b[k];
    trim(a);return a;
}
inline Poly multiply(const Poly& a,const Poly& b) {
    Poly c(a.size()+b.size()-1);
    for(std::size_t i=0;i<a.size();++i)if(a[i]!=0)
        for(std::size_t j=0;j<b.size();++j)if(b[j]!=0)c[i+j]+=a[i]*b[j];
    trim(c);return c;
}
inline Word subset(const Word& w,unsigned mask) {
    Word sub;for(auto a:w)if(mask>>a&1)sub.push_back(a);
    return chord::relabel_by_first_occurrence(sub);
}
inline bool connected(const Word& w) {
    unsigned n=w.size()/2;
    if(n<2)return true;
    std::vector<unsigned> start(n),end(n);
    for(unsigned i=0;i<w.size();++i)end[w[i]]=i;
    for(unsigned i=w.size();i-->0;)start[w[i]]=i;
    unsigned reached=1,last=0;
    while(last!=reached) {
        last=reached;
        for(unsigned a=0;a<n;++a)if(reached>>a&1)
            for(unsigned b=0;b<n;++b)
                if((start[a]<start[b] && start[b]<end[a] && end[a]<end[b]) ||
                   (start[b]<start[a] && start[a]<end[b] && end[b]<end[a]))
                    reached|=1u<<b;
    }
    return reached==(1u<<n)-1;
}
class Evaluator {
    std::unordered_map<std::string,Poly> moments_,primitives_;
public:
    std::size_t moment_entries()const{return moments_.size();}
    std::size_t primitive_entries()const{return primitives_.size();}
    Poly weight(const Word& input) {
        auto w=chord::canonical_unoriented(input);
        std::string key(w.begin(),w.end());
        if(auto it=moments_.find(key);it!=moments_.end())return it->second;
        if(w.empty())return {1};
        Word d=w;auto a=d[0];
        unsigned left=std::find(d.begin()+1,d.end(),a)-d.begin()-1;
        d.erase(std::remove(d.begin(),d.end(),a),d.end());
        std::vector<std::pair<unsigned,unsigned>> neighbours;
        for(unsigned i=0;i<left;++i) {
            auto end=std::find(d.begin()+i+1,d.end(),d[i]);
            if(end!=d.end() && unsigned(end-d.begin())>=left)
                neighbours.push_back({i,unsigned(end-d.begin())});
        }
        auto base=weight(d);Poly result(base.size()+1);
        for(unsigned k=0;k<base.size();++k) {
            result[k]-=4*neighbours.size()*base[k];result[k+1]+=base[k];
        }
        for(unsigned i=0;i<neighbours.size();++i)
            for(unsigned j=i+1;j<neighbours.size();++j) {
                auto [i0,i1]=neighbours[i];auto [j0,j1]=neighbours[j];
                auto ai=d[i0],aj=d[j0];
                d[i0]=aj;d[j1]=ai;auto positive=weight(d);
                d[j0]=ai;d[j1]=aj;auto negative=weight(d);
                d[i0]=ai;d[j0]=aj;
                auto diff=add(positive,negative,-1);
                result.resize(std::max(result.size(),diff.size()));
                for(unsigned k=0;k<diff.size();++k)result[k]+=4*diff[k];
            }
        trim(result);moments_.emplace(std::move(key),result);return result;
    }
    Poly primitive(const Word& input) {
        auto w=chord::canonical_unoriented(input);
        if(w.empty())return {0};
        std::string key(w.begin(),w.end());
        if(auto it=primitives_.find(key);it!=primitives_.end())return it->second;
        if(!connected(w)) {primitives_.emplace(std::move(key),Poly{0});return {0};}
        unsigned n=w.size()/2;
        if(n>=25)throw std::runtime_error("At most 24 chords in subset recursion");
        unsigned full=(1u<<n)-1,rest=full^1;
        auto result=weight(w);
        if(rest) {
            unsigned s=(rest-1)&rest;
            while(true) {
                auto block=s|1;
                result=add(std::move(result),multiply(primitive(subset(w,block)),
                                                    weight(subset(w,full^block))),-1);
                if(!s)break;s=(s-1)&rest;
            }
        }
        trim(result);primitives_.emplace(std::move(key),result);return result;
    }
};
inline Word leaf(const Word& input,unsigned root=0) {
    auto w=chord::relabel_by_first_occurrence(input);
    if(w.empty())throw std::runtime_error("Leaf needs a nonempty diagram");
    unsigned fresh=w.size()/2;Word out;bool inserted=false;
    for(auto a:w) {
        if(a==root && !inserted) {out.push_back(fresh);out.push_back(a);out.push_back(fresh);inserted=true;}
        else out.push_back(a);
    }
    return chord::relabel_by_first_occurrence(out);
}
}
