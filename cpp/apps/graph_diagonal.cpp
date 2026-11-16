#include <boost/json/src.hpp>
#include <array>
#include <bit>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
namespace json=boost::json;
using Poly=std::array<long long,14>;
Poly primitive(unsigned n,const std::vector<unsigned>& adj){
    std::vector<Poly> m(1u<<n),k(1u<<n);m[0][0]=1;
    for(unsigned mask=1;mask<(1u<<n);++mask){unsigned diag=mask,order=std::popcount(mask);
        while(true){unsigned basis[16]{},rank=0;
            for(unsigned v=0;v<n;++v)if(mask>>v&1){auto r=adj[v]&mask;if(diag>>v&1)r^=1u<<v;
                while(r){unsigned p=31-std::countl_zero(r);if(!basis[p]){basis[p]=r;++rank;break;}r^=basis[p];}}
            m[mask][order-rank]+=(std::popcount(diag)%2?-1:1);
            if(!diag)break;diag=(diag-1)&mask;
        }
        k[mask]=m[mask];unsigned anchor=mask&-mask,rest=mask^anchor;
        if(rest){unsigned s=(rest-1)&rest;while(true){unsigned block=s|anchor,other=mask^block;
            unsigned left=std::popcount(block),right=std::popcount(other);
            for(unsigned i=0;i<=left;++i)if(k[block][i])for(unsigned j=0;j<=right;++j)if(m[other][j])k[mask][i+j]-=k[block][i]*m[other][j];
            if(!s)break;s=(s-1)&rest;}}
    }return k.back();
}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("Expected JSON array of graphs");
 std::ifstream f(argv[1]);std::ostringstream t;t<<f.rdbuf();auto d=json::parse(t.str());
 for(const auto& item:d.as_array()){
  const auto& g=item.as_object();unsigned n=g.at("n").to_number<unsigned>();
  if(n>13)throw std::runtime_error("At most 13 vertices");
  std::vector<unsigned> adj(n);for(const auto& e:g.at("edges").as_array()){const auto& a=e.as_array();unsigned u=a[0].to_number<unsigned>(),v=a[1].to_number<unsigned>();adj[u]|=1u<<v;adj[v]|=1u<<u;}
  auto p=primitive(n,adj);json::array c;for(unsigned k=0;k<=n;++k)c.push_back(p[k]);
  std::cout<<json::serialize(json::object{{"id",g.at("id")},{"coefficients",c}})<<'\n'<<std::flush;
 }
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
