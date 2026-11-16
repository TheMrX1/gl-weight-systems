#include <boost/json/src.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using Int=boost::multiprecision::cpp_int;
namespace json=boost::json;

unsigned number(const std::string& text) {
    if(text.empty() || text.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("Expected a nonnegative integer");
    auto value=std::stoull(text);
    if(value>1000000000ULL)throw std::invalid_argument("Parameter exceeds 10^9");
    return unsigned(value);
}

Int power(Int x,unsigned e) {
    Int result=1;
    while(e) {if(e&1)result*=x;e>>=1;if(e)x*=x;}
    return result;
}

std::vector<Int> pp(unsigned n,unsigned a,unsigned b) {
    
    std::vector<Int> series(8);series[0]=1;
    for(unsigned i=0;i<n;++i) {
        const unsigned row=i==0?a:(i<b?1:0);
        const Int level=Int(row)+n-i-1;
        auto old=series;
        for(unsigned k=1;k<series.size();++k)
            series[k]=old[k]-(level+1)*old[k-1]+level*series[k-1];
    }
    if(series[1]!=-Int(n))throw std::logic_error("PP normalization failed");
    std::vector<Int> c(7);c[0]=n;
    for(unsigned m=1;m<=6;++m)c[m]=-series[m+1];
    return c;
}

json::object compute(unsigned n,unsigned a,unsigned b,const json::array& delta) {
    if(!n || n>128 || !a || !b || b>n)
        throw std::invalid_argument("Require 1<=N<=128, a>=1, 1<=b<=N");
    auto c=pp(n,a,b);
    Int difference=0;
    for(const auto& item:delta) {
        const auto& t=item.as_object();Int term(std::string(t.at("coefficient").as_string()));
        const auto& powers=t.at("powers").as_array();
        for(unsigned m=0;m<powers.size();++m) {
            unsigned e=powers[m].to_number<unsigned>();
            if(e) {
                if(m>=c.size())throw std::logic_error("Unexpected higher Casimir");
                term*=power(c[m],e);
            }
        }
        difference+=term;
    }
    const Int N=n,A=a,B=b,S=A+B-1;
    const Int h=N*(N+A-B)-3*S;
    const Int factored=-20*N*(N-B)*(N+A)*(A-1)*(B-1)*S*S*h*h;
    if(difference!=factored)throw std::logic_error("General hook identity failed");
    std::string reason;
    if(a==1 || b==1)reason="row_or_column";
    else if(b==n)reason="full_height_hook";
    else if(n<=3)reason="rank_at_most_three_on_this_pair";
    else if(a==2 && b==n-1)reason="adjoint_hook";
    else reason="distinguishes_this_pair";
    const bool zero=reason!="distinguishes_this_pair";
    if(zero!=(difference==0) || difference>0)throw std::logic_error("Classification failed");
    return {{"N",n},{"first_row",a},{"first_column",b},
            {"ordinary_difference",difference.str()},{"primitive_difference",difference.str()},
            {"reason",reason},{"PP_and_factorization_agree",true},
            {"scope","Fixed order-11 mutant pair D,Dprime"}};
}

int main(int argc,char** argv) {
    try {
        const auto file=std::filesystem::absolute(argv[0]).parent_path()/"universal_delta.json";
        std::ifstream f(file);if(!f)throw std::runtime_error("Cannot read universal_delta.json next to executable");
        std::ostringstream buffer;buffer<<f.rdbuf();auto delta=json::parse(buffer.str()).as_array();
        if(argc==2 && std::string(argv[1])=="--batch") {
            std::string line;
            while(std::getline(std::cin,line)) {
                std::istringstream in(line);std::string n,a,b,extra;
                if(!(in>>n>>a>>b) || (in>>extra))throw std::invalid_argument("Expected N a b per line");
                std::cout<<json::serialize(compute(number(n),number(a),number(b),delta))<<'\n';
            }
        } else if(argc==4) {
            std::cout<<json::serialize(compute(number(argv[1]),number(argv[2]),number(argv[3]),delta))<<'\n';
        } else throw std::invalid_argument("Usage: hook_difference N a b, or hook_difference --batch");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
