#include <boost/json/src.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using Int=boost::multiprecision::cpp_int;
namespace json=boost::json;

unsigned number(const std::string& text) {
    if(text.empty()||text.find_first_not_of("0123456789")!=std::string::npos)
        throw std::invalid_argument("expected a nonnegative integer");
    auto x=std::stoull(text);
    if(x>1000000000ULL)throw std::invalid_argument("parameter exceeds 10^9");
    return static_cast<unsigned>(x);
}
std::vector<unsigned> partition(const std::string& text) {
    std::istringstream stream(text);std::string token;std::vector<unsigned> rows;
    while(stream>>token)rows.push_back(number(token));
    if(!std::is_sorted(rows.begin(),rows.end(),std::greater<unsigned>()))
        throw std::invalid_argument("Young row lengths must be nonincreasing");
    while(!rows.empty()&&!rows.back())rows.pop_back();
    return rows;
}
void rational_factor(std::vector<Int>& series,const Int& numerator,const Int& denominator) {
    auto old=series;
    for(std::size_t k=1;k<series.size();++k)
        series[k]=old[k]-numerator*old[k-1]+denominator*series[k-1];
}
std::vector<Int> casimirs(unsigned n,const std::vector<unsigned>& rows,unsigned degree,bool corners) {
    std::vector<Int> series(degree+2);series[0]=1;
    if(corners) {
        rational_factor(series,Int(n)-rows.size(),0);
        for(std::size_t start=0;start<rows.size();) {
            std::size_t end=start+1;while(end<rows.size()&&rows[end]==rows[start])++end;
            rational_factor(series,Int(n)+rows[start]-start,Int(n)+rows[start]-end);
            start=end;
        }
    } else {
        for(unsigned i=0;i<n;++i) {
            Int l=Int(i<rows.size()?rows[i]:0)+n-i-1;
            rational_factor(series,l+1,l);
        }
    }
    if(series[1]!=-Int(n))throw std::logic_error("PP normalization failed");
    std::vector<Int> result(degree+1);result[0]=n;
    for(unsigned k=1;k<=degree;++k)result[k]=-series[k+1];
    return result;
}
Int power(Int x,unsigned exponent) {
    Int result=1;
    for(;exponent;exponent>>=1,x*=x)if(exponent&1)result*=x;
    return result;
}
Int dimension(unsigned n,const std::vector<unsigned>& rows) {
    Int numerator=1,denominator=1;
    for(unsigned i=0;i<n;++i)for(unsigned j=i+1;j<n;++j) {
        numerator*=Int(i<rows.size()?rows[i]:0)-Int(j<rows.size()?rows[j]:0)+j-i;
        denominator*=j-i;
    }
    if(numerator%denominator!=0)throw std::logic_error("Weyl dimension is not integral");
    return numerator/denominator;
}
Int evaluate(const std::filesystem::path& path,const std::vector<Int>& cs) {
    std::ifstream stream(path);if(!stream)throw std::runtime_error("cannot open polynomial: "+path.string());
    std::string text((std::istreambuf_iterator<char>(stream)),{});
    auto root=json::parse(text);if(!root.is_array())throw std::invalid_argument("polynomial must be a JSON array");
    Int total=0;
    for(const auto& item:root.as_array()) {
        const auto& term=item.as_object();const auto& powers=term.at("powers").as_array();
        Int value(std::string(term.at("coefficient").as_string()));
        for(std::size_t k=0;k<powers.size();++k) {
            if(!powers[k].is_int64()||powers[k].as_int64()<0||powers[k].as_int64()>1000000)
                throw std::invalid_argument("invalid monomial exponent");
            auto exponent=static_cast<unsigned>(powers[k].as_int64());
            if(!exponent)continue;
            if(k>=cs.size())throw std::invalid_argument("increase --degree for this polynomial");
            value*=power(cs[k],exponent);
        }
        total+=value;
    }
    return total;
}
std::string decimal(const Int& value){return value.convert_to<std::string>();}

int main(int argc,char** argv) {
    try {
        unsigned n=0,degree=22,width=0,height=0;bool has_n=false,has_shape=false,rectangle=false;
        std::vector<unsigned> rows;std::vector<std::filesystem::path> paths;
        auto next=[&](int& i){if(++i>=argc)throw std::invalid_argument("missing option value");return std::string(argv[i]);};
        for(int i=1;i<argc;++i) {
            std::string option=argv[i];
            if(option=="--N"){n=number(next(i));has_n=true;}
            else if(option=="--degree")degree=number(next(i));
            else if(option=="--young") {
                if(has_shape)throw std::invalid_argument("supply exactly one shape");
                rows=partition(next(i));has_shape=true;
            } else if(option=="--rectangle") {
                if(has_shape)throw std::invalid_argument("supply exactly one shape");
                width=number(next(i));height=number(next(i));
                if(!width||!height||height>128)throw std::invalid_argument("rectangle must be nonempty, height <=128");
                rows.assign(height,width);has_shape=true;rectangle=true;
            } else if(option=="--polynomial")paths.emplace_back(next(i));
            else if(option=="--help") {
                std::cout<<"representation_weights --N N (--young '2 1' | --rectangle WIDTH HEIGHT) [--degree 22] [--polynomial FILE.json ...]\n";
                return 0;
            } else throw std::invalid_argument("unknown option: "+option);
        }
        if(!has_n||!has_shape||n<1||n>128||degree<2||degree>128||rows.size()>n)
            throw std::invalid_argument("require 1<=N<=128, 2<=degree<=128, and at most N rows");
        const auto cs=casimirs(n,rows,degree,false),boundary=casimirs(n,rows,degree,true);
        if(cs!=boundary)throw std::logic_error("row and corner formulas disagree");
        Int c1=0,c2=0;
        for(std::size_t i=0;i<rows.size();++i){c1+=rows[i];c2+=Int(rows[i])*(Int(rows[i])+n-1-2*Int(i));}
        if(cs[1]!=c1||cs[2]!=c2)throw std::logic_error("C1/C2 checks failed");
        if(rectangle)for(unsigned k=1;k<=degree;++k)
            if(cs[k]!=Int(width)*height*power(Int(n)+width-height,k-1))
                throw std::logic_error("rectangular formula failed");
        json::object out;
        out["N"]=n;out["degree"]=degree;
        json::array young;for(auto x:rows)young.push_back(x);out["young_rows"]=young;
        out["representation_dimension"]=decimal(dimension(n,rows));
        json::array casimir_values;for(const auto& x:cs)casimir_values.emplace_back(decimal(x));
        out["casimirs"]=casimir_values;
        out["rows_equal_corner_formula"]=true;
        json::array weights;
        for(const auto& path:paths)weights.push_back(json::object{{"polynomial",path.string()},{"scalar",decimal(evaluate(path,cs))}});
        out["weights"]=weights;
        std::cout<<json::serialize(out)<<'\n';
    }catch(const std::exception& error){std::cerr<<"error: "<<error.what()<<'\n';return 1;}
}
