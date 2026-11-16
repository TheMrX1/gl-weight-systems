

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace {
constexpr unsigned N=4, DIM=20;
using Matrix=std::array<std::int64_t,DIM*DIM>;
using Ambient=std::array<unsigned,3>;  
struct Entry { unsigned row,column; std::int64_t value; };
std::uint32_t prime=65521;
std::int64_t mod(__int128 x){
    if(prime){x%=prime;if(x<0)x+=prime;}
    if(x<std::numeric_limits<std::int64_t>::min()||x>std::numeric_limits<std::int64_t>::max())
        throw std::overflow_error("direct integer contraction exceeded 64 bits");
    return static_cast<std::int64_t>(x);
}
Matrix multiply(const Matrix& a,const Matrix& b){
    Matrix r{};
    for(unsigned i=0;i<DIM;++i)for(unsigned k=0;k<DIM;++k)if(a[i*DIM+k])
        for(unsigned j=0;j<DIM;++j)r[i*DIM+j]+=a[i*DIM+k]*b[k*DIM+j];
    return r;
}
std::array<Matrix,16> representation(){
    std::vector<Ambient> ambient;
    std::map<Ambient,unsigned> ambient_index;
    std::map<Ambient,std::vector<unsigned>> fibers;
    for(unsigned a=0;a<N;++a)for(unsigned b=a;b<N;++b)for(unsigned c=0;c<N;++c){
        Ambient t{a,b,c};unsigned k=static_cast<unsigned>(ambient.size());ambient.push_back(t);ambient_index[t]=k;
        auto sorted=t;std::sort(sorted.begin(),sorted.end());fibers[sorted].push_back(k);
    }
    std::vector<std::pair<unsigned,unsigned>> basis;
    for(const auto& [weight,fiber]:fibers)
        for(unsigned k=1;k<fiber.size();++k)basis.push_back({fiber[k],fiber[0]});
    if(basis.size()!=DIM)throw std::logic_error("wrong representation dimension");
    std::array<Matrix,16> matrices{};
    for(unsigned i=0;i<N;++i)for(unsigned j=0;j<N;++j)for(unsigned column=0;column<DIM;++column){
        std::vector<std::int64_t> image(ambient.size());
        for(auto [source,sign]:{std::pair{basis[column].first,1},std::pair{basis[column].second,-1}}){
            auto t=ambient[source];
            for(unsigned slot=0;slot<3;++slot)if(t[slot]==j){
                auto u=t;u[slot]=i;if(u[0]>u[1])std::swap(u[0],u[1]);
                image[ambient_index.at(u)]+=sign;
            }
        }
        for(const auto& [weight,fiber]:fibers){
            std::int64_t sum=0;for(auto k:fiber)sum+=image[k];
            if(sum)throw std::logic_error("generator does not preserve kernel");
        }
        for(unsigned row=0;row<DIM;++row)matrices[i*N+j][row*DIM+column]=image[basis[row].first];
    }
    
    for(unsigned i=0;i<N;++i)for(unsigned j=0;j<N;++j)
        for(unsigned k=0;k<N;++k)for(unsigned l=0;l<N;++l){
            const auto ab=multiply(matrices[i*N+j],matrices[k*N+l]);
            const auto ba=multiply(matrices[k*N+l],matrices[i*N+j]);
            for(unsigned x=0;x<DIM*DIM;++x)
                if(ab[x]-ba[x]!=(j==k?matrices[i*N+l][x]:0)-(i==l?matrices[k*N+j][x]:0))
                    throw std::logic_error("gl4 commutator failed");
        }
    return matrices;
}
std::vector<std::int64_t> casimirs(const std::array<Matrix,16>& matrices){
    
    constexpr unsigned size=N*DIM;
    std::vector<std::int64_t> block(size*size),power(size*size);
    for(unsigned i=0;i<N;++i)for(unsigned j=0;j<N;++j)
        for(unsigned a=0;a<DIM;++a)for(unsigned b=0;b<DIM;++b)
            block[(i*DIM+a)*size+j*DIM+b]=mod(matrices[i*N+j][a*DIM+b]);
    for(unsigned i=0;i<size;++i)power[i*size+i]=1;
    std::vector<std::int64_t> result{N};
    for(unsigned degree=1;degree<=22;++degree){
        std::vector<std::int64_t> next(size*size);
        for(unsigned i=0;i<size;++i)for(unsigned k=0;k<size;++k)if(power[i*size+k])
            for(unsigned j=0;j<size;++j)if(block[k*size+j])
                next[i*size+j]=mod(next[i*size+j]+__int128(power[i*size+k])*block[k*size+j]);
        power.swap(next);std::int64_t scalar=0;
        for(unsigned i=0;i<N;++i)scalar=mod(scalar+__int128(power[(i*DIM)*size+i*DIM]));
        for(unsigned a=0;a<DIM;++a)for(unsigned b=0;b<DIM;++b){
            std::int64_t value=0;
            for(unsigned i=0;i<N;++i)value=mod(value+__int128(power[(i*DIM+a)*size+i*DIM+b]));
            if(value!=(a==b?scalar:0))throw std::logic_error("Casimir not scalar");
        }
        result.push_back(scalar);
    }
    return result;
}
using Word=std::vector<unsigned>;
std::pair<Word,unsigned> choose_cut(Word word){
    Word best;unsigned best_width=100;
    for(unsigned cut=0;cut<word.size();++cut){
        std::rotate(word.begin(),word.begin()+1,word.end());
        std::set<unsigned> open;unsigned width=0;
        for(auto c:word){if(open.count(c))open.erase(c);else open.insert(c);width=std::max(width,unsigned(open.size()));}
        if(width<best_width){best=word;best_width=width;}
    }
    
    std::reverse(best.begin(),best.end());return {best,best_width};
}
std::int64_t direct_weight(const Word& original,const std::array<Matrix,16>& matrices,unsigned input){
    auto [word,width]=choose_cut(original);
    if(width>6)throw std::runtime_error("direct check width exceeds 6");
    std::array<std::vector<Entry>,16> sparse;
    for(unsigned g=0;g<16;++g)for(unsigned i=0;i<DIM;++i)for(unsigned j=0;j<DIM;++j)
        if(matrices[g][i*DIM+j])sparse[g].push_back({i,j,matrices[g][i*DIM+j]});
    std::vector<std::int64_t> values(DIM);values[input]=1;
    std::vector<std::uint32_t> active{0};std::vector<unsigned> open;
    std::size_t assignments=1;
    for(unsigned step=0;step<word.size();++step){
        auto c=word[step];auto found=std::find(open.begin(),open.end(),c);
        bool opening=found==open.end();unsigned slot=static_cast<unsigned>(found-open.begin());
        std::size_t next_assignments=opening?assignments*16:assignments/16;
        std::vector<std::int64_t> next(next_assignments*DIM);
        std::vector<std::uint32_t> next_active;
        std::vector<bool> touched(opening?0:next_assignments);
        for(auto key:active){
            if(opening){
                for(unsigned g=0;g<16;++g){
                    auto new_key=key+g*assignments;auto offset=new_key*DIM;
                    for(const auto& e:sparse[g])next[offset+e.row]=mod(next[offset+e.row]+__int128(e.value)*values[std::size_t(key)*DIM+e.column]);
                    if(std::any_of(next.data()+offset,next.data()+offset+DIM,[](auto x){return x!=0;}))next_active.push_back(static_cast<std::uint32_t>(new_key));
                }
            }else{
                unsigned shift=4*slot;unsigned g=(key>>shift)&15;
                unsigned dual=(g%N)*N+g/N;
                std::uint32_t new_key=(key&((std::uint32_t(1)<<shift)-1))|((key>>(shift+4))<<shift);
                auto offset=std::size_t(new_key)*DIM;
                if(!touched[new_key]){next_active.push_back(new_key);touched[new_key]=true;}
                for(const auto& e:sparse[dual])next[offset+e.row]=mod(next[offset+e.row]+__int128(e.value)*values[std::size_t(key)*DIM+e.column]);
            }
        }
        if(!opening){
            next_active.erase(std::remove_if(next_active.begin(),next_active.end(),[&](auto key){
                auto offset=std::size_t(key)*DIM;
                return std::all_of(next.data()+offset,next.data()+offset+DIM,[](auto x){return x==0;});
            }),next_active.end());
        }
        if(opening)open.push_back(c);else open.erase(found);
        assignments=next_assignments;values.swap(next);active.swap(next_active);
        std::cerr<<"direct step="<<step+1<<" open="<<open.size()<<" active="<<active.size()<<" width="<<width<<'\n';
    }
    if(!open.empty()||values.size()!=DIM)throw std::logic_error("unclosed chords");
    for(unsigned i=0;i<DIM;++i)if(i!=input&&values[i])throw std::logic_error("weight did not act as a scalar on test vector");
    return values[input];
}
}
int main(int argc,char** argv){
    try{
        if(argc>1){auto requested=std::stoul(argv[1]);if(requested>1000000007)throw std::invalid_argument("unsupported modulus");prime=static_cast<std::uint32_t>(requested);}
        if(prime && prime<100)throw std::invalid_argument("unsupported modulus");
        unsigned input=argc>2?static_cast<unsigned>(std::stoul(argv[2])):0;
        if(input>=DIM)throw std::invalid_argument("input column out of range");
        const auto matrices=representation();const auto cs=casimirs(matrices);
        const Word words[2]={
            {0,1,2,0,3,1,2,4,5,6,7,8,5,9,4,8,10,3,7,9,10,6},
            {3,4,5,6,7,8,5,9,4,8,10,0,1,2,0,3,1,2,7,9,10,6}};
        
        if(direct_weight({0,0},matrices,0)!=cs[2])throw std::logic_error("one-chord check failed");
        if(direct_weight({0,1,1,0},matrices,0)!=mod(__int128(cs[2])*cs[2]))throw std::logic_error("noncrossing check failed");
        if(direct_weight({0,1,0,1},matrices,0)!=mod(__int128(cs[2])*cs[2]+__int128(cs[1])*cs[1]-N*cs[2]))throw std::logic_error("crossing check failed");
        const auto first=direct_weight(words[0],matrices,input),second=direct_weight(words[1],matrices,input);
        std::cout<<"{\"method\":\"direct matrix-unit contraction in kernel(Sym2(V) tensor V -> Sym3(V))\",\"N\":4,\"dimension\":20,\"prime\":"<<prime<<",\"casimirs\":[";
        for(unsigned k=0;k<cs.size();++k){if(k)std::cout<<',';std::cout<<cs[k];}
        std::cout<<"],\"values\":["<<first<<','<<second<<"],\"difference\":"<<mod(std::int64_t(first)-second)<<"}\n";
        return first==second?2:0;
    }catch(const std::exception& e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}
}
