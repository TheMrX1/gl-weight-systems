#include <boost/json/src.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
namespace json=boost::json;
using Big=boost::multiprecision::cpp_int;
using Wide=__int128;
Big cumulant(unsigned n,const std::vector<Wide>& m){
 std::vector<Big> k(1u<<n);
 for(unsigned mask=1;mask<(1u<<n);++mask){k[mask]=m[mask];unsigned anchor=mask&-mask,rest=mask^anchor;
  if(rest){auto s=(rest-1)&rest;while(true){unsigned block=s|anchor;k[mask]-=k[block]*m[mask^block];if(!s)break;s=(s-1)&rest;}}
 }return k.back();
}
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("Expected input JSON graphs");
 std::ifstream f(argv[1]);std::ostringstream t;t<<f.rdbuf();auto data=json::parse(t.str());
 auto begin=std::chrono::steady_clock::now();
 for(auto& item:data.as_array()){
  auto& g=item.as_object();unsigned n=g.at("n").to_number<unsigned>();if(n>13)throw std::runtime_error("n<=13");
  std::vector<unsigned> adj(n);for(auto&e:g.at("edges").as_array()){auto&a=e.as_array();unsigned u=a[0].to_number<unsigned>(),v=a[1].to_number<unsigned>();adj[u]|=1u<<v;adj[v]|=1u<<u;}
  unsigned size=1;std::vector<unsigned> powers(n);for(unsigned i=0;i<n;++i){powers[i]=size;size*=3;}
  std::vector<unsigned> code(1u<<n);std::vector<int> sign(1u<<n,1);
  for(unsigned s=1;s<(1u<<n);++s){unsigned bit=s&-s,v=__builtin_ctz(bit),rest=s^bit;code[s]=code[rest]+powers[v];sign[s]=sign[rest]*(__builtin_parity(adj[v]&rest)?-1:1);}
  std::vector<unsigned> positive(size);
  for(unsigned state=0;state<size;++state){unsigned r=state;for(unsigned v=0;v<n;++v){if(r%3)positive[state]|=1u<<v;r/=3;}}
  std::vector<Wide> dp(size),next(size);
  for(unsigned s=0;s<(1u<<n);++s)dp[code[s]]=sign[s];
  auto moments=[&](){std::vector<Wide> m(1u<<n);for(unsigned mask=0;mask<(1u<<n);++mask){unsigned state=2*code[mask],s=mask;Wide sum=0;while(true){sum+=sign[s]*dp[state-code[s]];if(!s)break;s=(s-1)&mask;}m[mask]=sum;}return m;};
  std::vector<Big> values(17); 
  for(unsigned colors=2;colors<=16;++colors){
   values[colors]=cumulant(n,moments());
   if(colors<16){for(unsigned state=0;state<size;++state){unsigned s=positive[state];Wide sum=0;while(true){sum+=sign[s]*dp[state-code[s]];if(!s)break;s=(s-1)&positive[state];}next[state]=sum;}dp.swap(next);}
  }
  Big expected=0,power=1;for(auto& a:g.at("primitive_scaled").as_array()){expected+=Big(std::string(a.as_string()))*power;power*=3;}
  if(values[3]!=expected)throw std::runtime_error("S3 calibration failed");
  if(values[4]!=(values[3]<<n))throw std::runtime_error("S4 identity failed");
  json::array encoded;for(auto&v:values)encoded.push_back(json::value(v.str()));
  std::cout<<json::serialize(json::object{{"id",g.at("id")},{"primitive_S_values_N0_to16",encoded},{"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()}})<<'\n'<<std::flush;
 }
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
