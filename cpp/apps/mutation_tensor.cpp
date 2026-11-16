

#define main original_fixed_witness_main
#include "audit_direct_gl4.cpp"
#undef main
#include <bit>
#include <chrono>
#include <sstream>
#include <string>
#include <unordered_map>

namespace full {
using U=std::uint32_t;
constexpr unsigned D=400;
using Vector=std::array<U,D>;
using Small=std::vector<std::vector<U>>;
using StateMap=std::unordered_map<U,Vector>;
U add(U a,U b){auto s=a+b;return s>=prime?s-prime:s;}
U sub(U a,U b){return a>=b?a-b:prime-(b-a);}
U mul(U a,U b){return U(std::uint64_t(a)*b%prime);}
U power(U a,U n){U r=1;for(;n;n>>=1,a=mul(a,a))if(n&1)r=mul(r,a);return r;}

struct Channel {
    std::array<int,4> weight{};
    unsigned dimension{};
    std::vector<Vector> basis;
    std::vector<unsigned> coordinates;
    Small flip;
};

void audit_spanning(const std::vector<Channel>& channels,const std::array<Matrix,16>& generators){
    auto lowering=[&](const Vector& input,unsigned k){
        Vector output{};const auto& g=generators[(k+1)*4+k];
        for(unsigned a=0;a<DIM;++a)for(unsigned b=0;b<DIM;++b)if(auto value=input[a*DIM+b]){
            for(unsigned x=0;x<DIM;++x){
                output[x*DIM+b]=add(output[x*DIM+b],mul(U(mod(g[x*DIM+a])),value));
                output[a*DIM+x]=add(output[a*DIM+x],mul(U(mod(g[x*DIM+b])),value));
            }
        }
        return output;
    };
    auto insert=[](Vector v,std::map<unsigned,Vector>& echelon){
        for(unsigned pos=0;pos<D;++pos){
            if(!v[pos])continue;auto found=echelon.find(pos);
            if(found!=echelon.end()){
                auto factor=v[pos];for(unsigned j=pos;j<D;++j)v[j]=sub(v[j],mul(factor,found->second[j]));
            }else{
                auto inverse=power(v[pos],prime-2);for(auto& x:v)x=mul(x,inverse);
                echelon.emplace(pos,v);return true;
            }
        }
        return false;
    };
    std::map<unsigned,Vector> total;
    for(const auto& channel:channels)for(const auto& highest:channel.basis){
        std::map<unsigned,Vector> local;std::vector<Vector> queue;
        insert(highest,local);queue.push_back(highest);
        for(unsigned cursor=0;cursor<queue.size();++cursor)for(unsigned k=0;k<3;++k){
            auto next=lowering(queue[cursor],k);
            if(insert(next,local))queue.push_back(next);
        }
        if(queue.size()!=channel.dimension)throw std::logic_error("lowering module dimension failed");
        for(const auto& v:queue)if(!insert(v,total))throw std::logic_error("dependent tensor component");
    }
    if(total.size()!=D)throw std::logic_error("full tensor basis rank is not 400");
}

std::vector<Channel> channels(const std::array<Matrix,16>& generators) {
    std::map<std::array<int,4>,std::vector<unsigned>> groups;
    for(unsigned a=0;a<DIM;++a)for(unsigned b=0;b<DIM;++b){
        std::array<int,4> weight{};
        for(unsigned i=0;i<4;++i)weight[i]=int(generators[i*4+i][a*DIM+a]+generators[i*4+i][b*DIM+b]);
        groups[weight].push_back(a*DIM+b);
    }
    std::vector<Channel> result;
    unsigned dimension_check=0;
    for(const auto& [weight,positions]:groups){
        std::map<unsigned,std::vector<U>> row_map;
        for(unsigned raising=0;raising<3;++raising)for(unsigned col=0;col<positions.size();++col){
            const auto a=positions[col]/DIM,b=positions[col]%DIM;
            const auto& g=generators[raising*4+raising+1];
            for(unsigned x=0;x<DIM;++x){
                for(auto [target,value]:{std::pair{x*DIM+b,g[x*DIM+a]},std::pair{a*DIM+x,g[x*DIM+b]}}){
                    if(!value)continue;
                    auto& row=row_map[raising*D+target];if(row.empty())row.resize(positions.size());
                    row[col]=add(row[col],U(mod(value)));
                }
            }
        }
        Small rows;for(auto& [key,row]:row_map)rows.push_back(row);
        std::vector<unsigned> pivots;
        for(unsigned col=0,r=0;col<positions.size()&&r<rows.size();++col){
            unsigned pivot=r;while(pivot<rows.size()&&!rows[pivot][col])++pivot;
            if(pivot==rows.size())continue;
            std::swap(rows[r],rows[pivot]);auto inv=power(rows[r][col],prime-2);
            for(auto& x:rows[r])x=mul(x,inv);
            for(unsigned k=0;k<rows.size();++k)if(k!=r&&rows[k][col]){
                auto factor=rows[k][col];
                for(unsigned j=0;j<positions.size();++j)rows[k][j]=sub(rows[k][j],mul(factor,rows[r][j]));
            }
            pivots.push_back(col);++r;
        }
        Channel channel;channel.weight=weight;
        for(unsigned col=0;col<positions.size();++col){
            if(std::find(pivots.begin(),pivots.end(),col)!=pivots.end())continue;
            Vector v{};v[positions[col]]=1;
            for(unsigned r=0;r<pivots.size();++r)v[positions[pivots[r]]]=sub(0,rows[r][col]);
            channel.basis.push_back(v);channel.coordinates.push_back(positions[col]);
        }
        if(channel.basis.empty())continue;
        std::int64_t numerator=1,denominator=1;
        for(unsigned i=0;i<4;++i)for(unsigned j=i+1;j<4;++j){
            if(weight[i]<weight[j])throw std::logic_error("nondominant highest weight");
            numerator*=weight[i]-weight[j]+j-i;denominator*=j-i;
        }
        if(numerator%denominator)throw std::logic_error("nonintegral Weyl dimension");
        channel.dimension=unsigned(numerator/denominator);
        dimension_check+=channel.dimension*unsigned(channel.basis.size());
        unsigned m=unsigned(channel.basis.size());channel.flip.assign(m,std::vector<U>(m));
        for(unsigned col=0;col<m;++col){
            for(unsigned r=0;r<m;++r){auto pos=channel.coordinates[r];channel.flip[r][col]=channel.basis[col][(pos%DIM)*DIM+pos/DIM];}
            for(unsigned pos=0;pos<D;++pos){
                U reconstructed=0;for(unsigned r=0;r<m;++r)reconstructed=add(reconstructed,mul(channel.basis[r][pos],channel.flip[r][col]));
                if(reconstructed!=channel.basis[col][(pos%DIM)*DIM+pos/DIM])throw std::logic_error("flip projection failed");
            }
        }
        result.push_back(channel);
    }
    if(dimension_check!=D)throw std::logic_error("tensor decomposition dimension mismatch");
    return result;
}

struct Endpoint {unsigned label;bool left;};
std::vector<Endpoint> schedule(Word a,Word c,unsigned& width){
    std::reverse(a.begin(),a.end());std::reverse(c.begin(),c.end());
    unsigned rows=unsigned(a.size()+1),cols=unsigned(c.size()+1);
    std::vector<U> ma(rows),mc(cols);
    for(unsigned i=0;i<a.size();++i)ma[i+1]=ma[i]^(U(1)<<a[i]);
    for(unsigned j=0;j<c.size();++j)mc[j+1]=mc[j]^(U(1)<<c[j]);
    std::vector<unsigned> costs(rows*cols,99);std::vector<bool> left(rows*cols);
    costs[0]=0;
    for(unsigned i=0;i<rows;++i)for(unsigned j=0;j<cols;++j){
        if(!i&&!j)continue;unsigned active=unsigned(std::popcount(ma[i]^mc[j]));
        unsigned x=i?std::max(costs[(i-1)*cols+j],active):99;
        unsigned y=j?std::max(costs[i*cols+j-1],active):99;
        left[i*cols+j]=x<=y;costs[i*cols+j]=std::min(x,y);
    }
    width=costs.back();if(width>5)throw std::runtime_error("independent audit width exceeds 5");
    std::vector<Endpoint> plan;unsigned i=rows-1,j=cols-1;
    while(i||j){if(left[i*cols+j]){plan.push_back({a[--i],true});}else{plan.push_back({c[--j],false});}}
    std::reverse(plan.begin(),plan.end());return plan;
}

Vector apply(const Vector& input,const std::vector<Entry>& sparse,bool left){
    Vector output{};
    for(const auto& entry:sparse){U coefficient=U(mod(entry.value));
        if(left){for(unsigned b=0;b<DIM;++b)if(auto x=input[entry.column*DIM+b])
            output[entry.row*DIM+b]=add(output[entry.row*DIM+b],mul(coefficient,x));}
        else{for(unsigned a=0;a<DIM;++a)if(auto x=input[a*DIM+entry.column])
            output[a*DIM+entry.row]=add(output[a*DIM+entry.row],mul(coefficient,x));}
    }
    return output;
}
bool nonzero(const Vector& v){return std::any_of(v.begin(),v.end(),[](auto x){return x!=0;});}
Vector contract(const Word& a,const Word& c,const Vector& initial,const std::array<std::vector<Entry>,16>& generators){
    unsigned width;const auto plan=schedule(a,c,width);
    StateMap states;states[0]=initial;std::vector<unsigned> open;
    for(const auto endpoint:plan){
        auto found=std::find(open.begin(),open.end(),endpoint.label);
        bool opening=found==open.end();unsigned slot=unsigned(found-open.begin());
        StateMap next;
        for(const auto& [key,value]:states){
            if(opening){
                for(unsigned g=0;g<16;++g){
                    auto image=apply(value,generators[g],endpoint.left);
                    if(nonzero(image))next.emplace(key|(g<<(4*open.size())),std::move(image));
                }
            }else{
                unsigned shift=4*slot,g=(key>>shift)&15,dual=(g%4)*4+g/4;
                auto image=apply(value,generators[dual],endpoint.left);if(!nonzero(image))continue;
                U new_key=(key&((U(1)<<shift)-1))|((key>>(shift+4))<<shift);
                auto& output=next[new_key];for(unsigned i=0;i<D;++i)output[i]=add(output[i],image[i]);
            }
        }
        if(opening)open.push_back(endpoint.label);else open.erase(found);
        if(next.size()>300000)throw std::runtime_error("independent audit state cap exceeded");
        states.swap(next);
    }
    if(!open.empty()||states.size()>1)throw std::logic_error("unclosed contraction");
    return states.empty()?Vector{}:states.at(0);
}

std::vector<Small> operator_blocks(const Word& a,const Word& c,const std::vector<Channel>& ch,
                                  const std::array<std::vector<Entry>,16>& generators){
    std::vector<Small> result;
    for(unsigned h=0;h<ch.size();++h){
        const auto& channel=ch[h];unsigned m=unsigned(channel.basis.size());
        Small matrix(m,std::vector<U>(m));
        for(unsigned col=0;col<m;++col){
            const auto output=contract(a,c,channel.basis[col],generators);
            for(unsigned row=0;row<m;++row)matrix[row][col]=output[channel.coordinates[row]];
            for(unsigned pos=0;pos<D;++pos){
                U reconstructed=0;
                for(unsigned row=0;row<m;++row)reconstructed=add(reconstructed,mul(channel.basis[row][pos],matrix[row][col]));
                if(reconstructed!=output[pos])throw std::logic_error("full-channel projection failed");
            }
        }
        result.push_back(matrix);
    }
    return result;
}
U sew(const std::vector<Small>& a,const std::vector<Small>& b,const std::vector<Channel>& ch){
    U total=0;
    for(unsigned h=0;h<ch.size();++h){
        unsigned m=unsigned(ch[h].basis.size());U trace=0;
        for(unsigned i=0;i<m;++i)for(unsigned j=0;j<m;++j)for(unsigned k=0;k<m;++k)
            trace=add(trace,mul(ch[h].flip[i][j],mul(a[h][j][k],b[h][k][i])));
        total=add(total,mul(ch[h].dimension,trace));
    }
    return mul(total,power(DIM,prime-2));
}
Word parse(const char* s){Word word;std::istringstream stream(s);unsigned x;while(stream>>x)word.push_back(x);if(!stream.eof())throw std::invalid_argument("invalid word");return word;}
}

int main(int argc,char** argv){
    try{
        if(argc!=5&&argc!=6)throw std::invalid_argument("A B C E [prime]");
        prime=argc==6?std::uint32_t(std::stoul(argv[5])):65521;
        if(prime<1000)throw std::invalid_argument("use an independently known prime greater than 1000");
        auto a=full::parse(argv[1]),b=full::parse(argv[2]),c=full::parse(argv[3]),e=full::parse(argv[4]);
        std::map<unsigned,unsigned> counts;
        for(auto w:{a,c})for(auto x:w)++counts[x];
        for(auto [x,n]:counts)if(n!=2)throw std::invalid_argument("invalid first share");
        std::set<unsigned> labels;for(auto [x,n]:counts)labels.insert(x);
        counts.clear();for(auto w:{b,e})for(auto x:w){if(labels.count(x))throw std::invalid_argument("shares overlap");++counts[x];}
        for(auto [x,n]:counts)if(n!=2)throw std::invalid_argument("invalid second share");
        
        auto relabel=[](Word& first,Word& second){std::map<unsigned,unsigned> names;for(auto* w:{&first,&second})for(auto& x:*w){auto [it,inserted]=names.emplace(x,unsigned(names.size()));x=it->second;}};
        relabel(a,c);relabel(b,e);
        auto started=std::chrono::steady_clock::now();
        const auto matrices=representation();const auto channels=full::channels(matrices);
        full::audit_spanning(channels,matrices);

        std::array<std::uint64_t,16> norms{};
        for(unsigned g=0;g<16;++g)for(unsigned row=0;row<DIM;++row){
            std::uint64_t row_sum=0;
            for(unsigned col=0;col<DIM;++col)row_sum+=std::uint64_t(std::abs(matrices[g][row*DIM+col]));
            norms[g]=std::max(norms[g],row_sum);
        }
        std::uint64_t norm_pair_bound=0;
        for(unsigned g=0;g<16;++g)norm_pair_bound+=norms[g]*norms[(g%4)*4+g/4];
        std::array<std::vector<Entry>,16> generators;
        for(unsigned g=0;g<16;++g)for(unsigned i=0;i<DIM;++i)for(unsigned j=0;j<DIM;++j)
            if(matrices[g][i*DIM+j])generators[g].push_back({i,j,matrices[g][i*DIM+j]});
        std::cerr<<"channels="<<channels.size()<<" total dimension=400\n";
        for(const auto& ch:channels){for(auto x:ch.weight)std::cerr<<x<<',';std::cerr<<" dim="<<ch.dimension<<" mult="<<ch.basis.size()<<'\n';}
        const auto complement=full::operator_blocks(b,e,channels,generators);
        std::vector<full::U> values;
        for(unsigned kind=0;kind<4;++kind){
            auto first=(kind&2)?c:a,second=(kind&2)?a:c;
            if(kind&1){std::reverse(first.begin(),first.end());std::reverse(second.begin(),second.end());}
            const auto blocks=full::operator_blocks(first,second,channels,generators);
            values.push_back(full::sew(blocks,complement,channels));
            std::cerr<<"kind="<<kind<<" full_weight="<<values.back()<<std::endl;
        }
        std::cout<<"{\"prime\":"<<prime<<",\"N\":4,\"highest_weight\":[2,1,0,0],\"dimension\":20,\"full_lowering_basis_rank\":400,\"operator_norm_pair_bound\":"<<norm_pair_bound<<",\"method\":\"independent full tensor-square decomposition and matrix-unit contraction\",\"order\":\"identity,R,S,RS\",\"values\":[";
        for(unsigned k=0;k<4;++k){if(k)std::cout<<',';std::cout<<values[k];}
        std::cout<<"],\"channels\":[";
        for(unsigned h=0;h<channels.size();++h){if(h)std::cout<<',';const auto& ch=channels[h];std::cout<<"{\"highest_weight\":[";for(unsigned k=0;k<4;++k){if(k)std::cout<<',';std::cout<<ch.weight[k];}std::cout<<"],\"dimension\":"<<ch.dimension<<",\"multiplicity\":"<<ch.basis.size()<<'}';}
        std::cout<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
    }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}
}
