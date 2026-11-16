#include <boost/json/src.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <bit>
#include <array>
namespace json=boost::json;
std::array<long long,14> full(unsigned n,const std::vector<unsigned>& adj) {
    std::vector<unsigned> ranks(1u<<n);std::vector<std::array<long long,14>> k(1u<<n);
    for(unsigned mask=1;mask<(1u<<n);++mask){unsigned basis[16]{},rank=0;
        for(unsigned v=0;v<n;++v)if(mask>>v&1){auto row=adj[v]&mask;while(row){unsigned p=31-std::countl_zero(row);if(!basis[p]){basis[p]=row;++rank;break;}row^=basis[p];}}
        ranks[mask]=std::popcount(mask)-rank;k[mask][ranks[mask]]=1;
        unsigned anchor=mask&-mask,rest=mask^anchor;
        if(rest){auto s=(rest-1)&rest;while(true){unsigned block=s|anchor,shift=ranks[mask^block];
            for(unsigned j=0;j+shift<=n;++j)k[mask][j+shift]-=k[block][j];if(!s)break;s=(s-1)&rest;}}
    }return k.back();
}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("Expected JSON array of graphs");
 std::ifstream f(argv[1]);std::ostringstream t;t<<f.rdbuf();auto d=json::parse(t.str());
 for(const auto& item:d.as_array()){
  const auto& g=item.as_object();unsigned n=g.at("n").to_number<unsigned>();
  if(n>13)throw std::runtime_error("At most 13 vertices");
  std::vector<unsigned> adj(n);for(const auto& e:g.at("edges").as_array()){const auto& a=e.as_array();unsigned u=a[0].to_number<unsigned>(),v=a[1].to_number<unsigned>();adj[u]|=1u<<v;adj[v]|=1u<<u;}
  auto p=full(n,adj);json::array c;for(unsigned k=0;k<=n;++k)c.push_back(p[k]);
  std::cout<<json::serialize(json::object{{"id",g.at("id")},{"coefficients",c}})<<'\n'<<std::flush;
 }
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
