#include "chord/symbolic_gl.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace chord {
namespace {
void rotate_cycle_to_minimum(Cycle& cycle) {
    if (cycle.empty()) {
        return;
    }
    const auto minimum = std::min_element(cycle.begin(), cycle.end());
    std::rotate(cycle.begin(), minimum, cycle.end());
}

Cycle rotate_left(const Cycle& cycle, std::size_t amount) {
    if (cycle.empty()) {
        return {};
    }
    Cycle result = cycle;
    amount %= result.size();
    std::rotate(result.begin(), result.begin() + static_cast<std::ptrdiff_t>(amount), result.end());
    return result;
}

std::string serialize(const Cycles& cycles) {
    std::string key;
    std::size_t size = 1U;
    for (const auto& cycle : cycles) {
        size += cycle.size() + 1U;
    }
    key.reserve(size);
    key.push_back(static_cast<char>(cycles.size()));
    for (const auto& cycle : cycles) {
        key.push_back(static_cast<char>(cycle.size()));
        for (const auto label : cycle) {
            key.push_back(static_cast<char>(label));
        }
    }
    return key;
}

bool is_increasing_identity_cycle(const Cycle& cycle) {
    for (std::size_t index = 0; index < cycle.size(); ++index) {
        if (cycle[index] != index) {
            return false;
        }
    }
    return true;
}

ExactInt checked_add(ExactInt a, ExactInt b) {
    ExactInt r; if (__builtin_add_overflow(a,b,&r)) throw std::overflow_error("coefficient addition exceeds signed 128 bits"); return r;
}
ExactInt checked_mul(ExactInt a, ExactInt b) {
    ExactInt r; if (__builtin_mul_overflow(a,b,&r)) throw std::overflow_error("coefficient product exceeds signed 128 bits"); return r;
}
MonomialKey monomial_product(MonomialKey a, MonomialKey b) {
    MonomialKey result=0;
    for(unsigned k=0;k<23;++k) {
        auto e=((a>>(5*k))&31)+((b>>(5*k))&31);
        if(e>31) throw std::overflow_error("monomial exponent exceeds 31");
        result |= e<<(5*k);
    }
    return result;
}
ExactPolynomial polynomial_product(const ExactPolynomial& a,const ExactPolynomial& b) {
    if(a.empty() || b.empty()) return {};
    ExactPolynomial r; r.reserve(a.size()*b.size());
    for(const auto& x:a) for(const auto& y:b)
        r.push_back({monomial_product(x.monomial,y.monomial),checked_mul(x.coefficient,y.coefficient)});
    std::sort(r.begin(),r.end(),[](const auto& x,const auto& y){return x.monomial<y.monomial;});
    std::size_t out=0;
    for(std::size_t i=0;i<r.size();) {
        auto key=r[i].monomial;ExactInt c=0;
        do {c=checked_add(c,r[i++].coefficient);} while(i<r.size()&&r[i].monomial==key);
        if(c) r[out++]={key,c};
    }
    r.resize(out);return r;
}
}
ExactPolynomial polynomial_sum(const ExactPolynomial& a,const ExactPolynomial& b,int sign) {
    ExactPolynomial r;r.reserve(a.size()+b.size()); std::size_t i=0,j=0;
    while(i<a.size()||j<b.size()) {
        ExactTerm t;
        if(j==b.size()||(i<a.size()&&a[i].monomial<b[j].monomial)) t=a[i++];
        else if(i==a.size()||b[j].monomial<a[i].monomial) {t=b[j++];t.coefficient=checked_mul(t.coefficient,sign);}
        else {t={a[i].monomial,checked_add(a[i].coefficient,checked_mul(b[j].coefficient,sign))};++i;++j;}
        if(t.coefficient)r.push_back(t);
    }
    return r;
}
std::string exact_integer_string(ExactInt x) {
    if(!x)return "0";
    bool negative=x<0; auto a=negative?MonomialKey(-(x+1))+1:MonomialKey(x);
    std::string s;while(a){s.push_back(char('0'+a%10));a/=10;}
    if(negative)s.push_back('-');std::reverse(s.begin(),s.end());return s;
}
std::string polynomial_string(const ExactPolynomial& p) {
    if(p.empty())return "0";std::ostringstream s;bool first=true;
    for(const auto& t:p) {
        auto c=exact_integer_string(t.coefficient);
        if(!first&&t.coefficient>0)s<<"+";
        s<<c;
        for(unsigned k=0;k<23;++k){unsigned e=unsigned((t.monomial>>(5*k))&31);if(!e)continue;
            s<<"*"<<(k==0?"N":"C"+std::to_string(k));if(e!=1)s<<"^"<<e;
        }
        first=false;
    }
    return s.str();
}
std::string polynomial_json(const ExactPolynomial& p) {
    std::ostringstream s;s<<"[";bool first=true;
    for(const auto& t:p){if(!first)s<<",";first=false;s<<"{\"coefficient\":\""<<exact_integer_string(t.coefficient)<<"\",\"powers\":[";
        for(unsigned k=0;k<23;++k){if(k)s<<",";s<<unsigned((t.monomial>>(5*k))&31);}s<<"]}";}
    s<<"]";return s.str();
}
SymbolicGlEvaluator::SymbolicGlEvaluator() {
    nodes_.push_back({Op::Zero});nodes_.push_back({Op::One});
    for(unsigned k=0;k<23;++k)leaves_[k]=node(Op::Casimir,k,0);
    cache_.reserve(1U<<18);
}
SymbolicGlEvaluator::NodeId SymbolicGlEvaluator::node(Op op,NodeId a,NodeId b){
    if(nodes_.size()>=std::numeric_limits<NodeId>::max())throw std::overflow_error("too many DAG nodes");
    nodes_.push_back({op,a,b});return NodeId(nodes_.size()-1);
}
SymbolicGlEvaluator::NodeId SymbolicGlEvaluator::add(NodeId a,NodeId b){return !a?b:!b?a:node(Op::Add,a,b);}
SymbolicGlEvaluator::NodeId SymbolicGlEvaluator::subtract(NodeId a,NodeId b){return a==b?0:!b?a:node(Op::Subtract,a,b);}
SymbolicGlEvaluator::NodeId SymbolicGlEvaluator::multiply(NodeId a,NodeId b){return !a||!b?0:a==1?b:b==1?a:node(Op::Multiply,a,b);}
SymbolicGlEvaluator::NodeId SymbolicGlEvaluator::power(NodeId a,std::size_t e){
    NodeId r=1;while(e){if(e&1)r=multiply(r,a);e>>=1;if(e)a=multiply(a,a);}return r;
}
ExactPolynomial SymbolicGlEvaluator::evaluate_chord_word(const Word& w) {
    if(w.size()>22)throw std::invalid_argument("symbolic evaluator supports at most 11 chords");
    const auto root=evaluate_impl(word_to_involution_cycles(w),0);
    std::cerr<<"DAG ready: nodes="<<nodes_.size()<<" cache="<<cache_.size()<<std::endl;
    cache_.clear();cache_.rehash(0);
    return expand(root);
}
ExactPolynomial SymbolicGlEvaluator::expand(NodeId root) {
#ifdef CHORD_AUDIT_WIDE_ARITHMETIC

    std::vector<bool> reachable(nodes_.size());std::vector<NodeId> pending{root};
    while(!pending.empty()){
        auto i=pending.back();pending.pop_back();if(reachable[i])continue;
        reachable[i]=true;const auto n=nodes_[i];
        if(n.op==Op::Add||n.op==Op::Subtract||n.op==Op::Multiply){pending.push_back(n.left);pending.push_back(n.right);}
    }
    std::vector<std::uint32_t> users(nodes_.size());
    for(NodeId i=0;i<nodes_.size();++i)if(reachable[i]){
        const auto n=nodes_[i];
        if(n.op==Op::Add||n.op==Op::Subtract||n.op==Op::Multiply){++users[n.left];++users[n.right];}
    }
    ++users[root];
    std::vector<std::unique_ptr<ExactPolynomial>> values(nodes_.size());
    std::uint64_t live_terms=0;auto tick=std::chrono::steady_clock::now();
    for(NodeId i=0;i<nodes_.size();++i){
        if(!users[i])continue;
        const auto n=nodes_[i];ExactPolynomial p;
        switch(n.op){
        case Op::Zero:break;
        case Op::One:p={{0,1}};break;
        case Op::Casimir:p={{MonomialKey(1)<<(5*n.left),1}};break;
        case Op::Add:case Op::Subtract:case Op::Multiply:
            p=n.op==Op::Multiply?polynomial_product(*values[n.left],*values[n.right]):
                polynomial_sum(*values[n.left],*values[n.right],n.op==Op::Add?1:-1);
            for(auto child:{n.left,n.right})if(!--users[child]){live_terms-=values[child]->size();values[child].reset();}
            break;
        }
        p.shrink_to_fit();live_terms+=p.size();values[i]=std::make_unique<ExactPolynomial>(std::move(p));
        if(live_terms>400000000ULL)throw std::runtime_error("wide audit exceeded 400 million live terms");
        if(std::chrono::steady_clock::now()-tick>std::chrono::seconds(10)){
            std::cerr<<"wide audit node="<<i<<"/"<<nodes_.size()<<" live_terms="<<live_terms<<std::endl;
            tick=std::chrono::steady_clock::now();
        }
    }
    return *values[root];
#else
    struct CompactTerm { std::uint64_t key; std::int64_t coefficient; };
    using CompactPolynomial=std::vector<CompactTerm>;

    std::array<unsigned,23> bits{},shift{};
    unsigned total=0;
    for(unsigned k=0;k<23;++k){
        unsigned bound=k?22/k:22;
        unsigned b=0;while(bound){++b;bound>>=1;}
        bits[k]=b;shift[k]=total;total+=b;
    }
    auto narrow=[](ExactInt x)->std::int64_t{
        if(x<std::numeric_limits<std::int64_t>::min()||x>std::numeric_limits<std::int64_t>::max())
            throw std::overflow_error("compact intermediate coefficient exceeds signed 64 bits");
        return static_cast<std::int64_t>(x);
    };
    auto sum=[&](const CompactPolynomial& a,const CompactPolynomial& b,int sign){
        CompactPolynomial p;p.reserve(a.size()+b.size());std::size_t i=0,j=0;
        while(i<a.size()||j<b.size()){
            CompactTerm t;
            if(j==b.size()||(i<a.size()&&a[i].key<b[j].key))t=a[i++];
            else if(i==a.size()||b[j].key<a[i].key){t=b[j++];t.coefficient=narrow(ExactInt(t.coefficient)*sign);}
            else {t={a[i].key,narrow(ExactInt(a[i].coefficient)+ExactInt(b[j].coefficient)*sign)};++i;++j;}
            if(t.coefficient)p.push_back(t);
        }
        p.shrink_to_fit();return p;
    };
    auto product=[&](const CompactPolynomial& a,const CompactPolynomial& b){
        CompactPolynomial p;p.reserve(a.size()*b.size());
        for(const auto& x:a)for(const auto& y:b){
            std::uint64_t key=0;
            for(unsigned k=0;k<23;++k){
                auto mask=(std::uint64_t(1)<<bits[k])-1;
                auto e=((x.key>>shift[k])&mask)+((y.key>>shift[k])&mask);
                if(e>mask)throw std::overflow_error("compact monomial exponent overflow");
                key|=e<<shift[k];
            }
            p.push_back({key,narrow(ExactInt(x.coefficient)*y.coefficient)});
        }
        std::sort(p.begin(),p.end(),[](const auto& x,const auto& y){return x.key<y.key;});
        std::size_t out=0;
        for(std::size_t i=0;i<p.size();){auto key=p[i].key;ExactInt c=0;
            do{c=checked_add(c,p[i++].coefficient);}while(i<p.size()&&p[i].key==key);
            if(c)p[out++]={key,narrow(c)};
        }
        p.resize(out);p.shrink_to_fit();return p;
    };
    std::vector<bool> reachable(nodes_.size());
    std::vector<NodeId> pending{root};std::uint64_t reachable_count=0;
    while(!pending.empty()){
        const auto i=pending.back();pending.pop_back();if(reachable[i])continue;
        reachable[i]=true;++reachable_count;const auto n=nodes_[i];
        if(n.op==Op::Add||n.op==Op::Subtract||n.op==Op::Multiply){pending.push_back(n.left);pending.push_back(n.right);}
    }
    std::cerr<<"reachable DAG nodes="<<reachable_count<<std::endl;
    std::vector<std::uint32_t> users(nodes_.size());
    for(NodeId i=0;i<nodes_.size();++i)if(reachable[i]){
        const auto n=nodes_[i];
        if(n.op==Op::Add||n.op==Op::Subtract||n.op==Op::Multiply){++users[n.left];++users[n.right];}
    }
    ++users[root];
    std::vector<std::unique_ptr<CompactPolynomial>> values(nodes_.size());
    std::uint64_t live_terms=0;
    auto tick=std::chrono::steady_clock::now();
    for(NodeId i=0;i<nodes_.size();++i){
        if(!users[i])continue;
        const auto n=nodes_[i];CompactPolynomial p;
        switch(n.op){
        case Op::Zero:break;
        case Op::One:p={{0,1}};break;
        case Op::Casimir:p={{std::uint64_t(1)<<shift[n.left],1}};break;
        case Op::Add:case Op::Subtract:case Op::Multiply:
            p=n.op==Op::Multiply?product(*values[n.left],*values[n.right]):sum(*values[n.left],*values[n.right],n.op==Op::Add?1:-1);
            for(auto child:{n.left,n.right})if(!--users[child]){live_terms-=values[child]->size();values[child].reset();}
            break;
        }
        live_terms+=p.size();values[i]=std::make_unique<CompactPolynomial>(std::move(p));
        if(live_terms>600000000ULL)throw std::runtime_error("symbolic live storage exceeded 600 million compact terms");
        if(std::chrono::steady_clock::now()-tick>std::chrono::seconds(10)){
            std::cerr<<"compact expansion node="<<i<<"/"<<nodes_.size()<<" live_terms="<<live_terms<<std::endl;
            tick=std::chrono::steady_clock::now();
        }
    }
    ExactPolynomial result;result.reserve(values[root]->size());
    for(const auto& t:*values[root]){MonomialKey key=0;
        for(unsigned k=0;k<23;++k)key|=MonomialKey((t.key>>shift[k])&((std::uint64_t(1)<<bits[k])-1))<<(5*k);
        result.push_back({key,t.coefficient});
    }
    std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.monomial<b.monomial;});
    return result;
#endif
}

SymbolicGlEvaluator::NormalizedState SymbolicGlEvaluator::normalize(const Cycles& input) const {
    NormalizedState result;
    Cycles cycles;
    std::vector<std::uint8_t> labels;
    for (const auto& cycle : input) {
        if (cycle.empty()) {
            ++result.empty_cycles;
            continue;
        }
        cycles.push_back(cycle);
        labels.insert(labels.end(), cycle.begin(), cycle.end());
    }
    if (cycles.empty()) {
        return result;
    }

    std::sort(labels.begin(), labels.end());
    if (std::adjacent_find(labels.begin(), labels.end()) != labels.end()) {
        throw std::invalid_argument("permutation cycles must use every label exactly once");
    }

    std::array<std::uint8_t, 256> compression;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        compression[labels[index]] = static_cast<std::uint8_t>(index);
    }
    for (auto& cycle : cycles) {
        for (auto& label : cycle) {
            label = static_cast<std::uint8_t>(compression[label]);
        }
    }

    const auto vertex_count = labels.size();
    Cycles best;
    bool initialized = false;
    auto candidate = cycles;

    for (std::size_t shift = 0; shift < vertex_count; ++shift) {
        if (shift != 0U) {
            for (auto& cycle : candidate) {
                for (auto& label : cycle) {
                    label = static_cast<std::uint8_t>(
                        label == 0U ? vertex_count - 1U
                                    : static_cast<std::size_t>(label) - 1U);
                }
            }
        }
        for (auto& cycle : candidate) {
            rotate_cycle_to_minimum(cycle);
        }
        std::sort(candidate.begin(), candidate.end());

        std::size_t prefix_vertices = 0U;
        std::size_t prefix_max = 0U;
        for (std::size_t cycle_index = 0; cycle_index + 1U < candidate.size(); ++cycle_index) {
            prefix_vertices += candidate[cycle_index].size();
            for (const auto label : candidate[cycle_index]) {
                prefix_max = std::max(prefix_max, static_cast<std::size_t>(label));
            }
            if (prefix_max + 1U == prefix_vertices) {
                result.decomposable = true;
                result.left.assign(candidate.begin(),
                                   candidate.begin() + static_cast<std::ptrdiff_t>(cycle_index + 1U));
                result.right.assign(candidate.begin() + static_cast<std::ptrdiff_t>(cycle_index + 1U),
                                    candidate.end());
                for (auto& cycle : result.right) {
                    for (auto& label : cycle) {
                        label = static_cast<std::uint8_t>(
                            static_cast<std::size_t>(label) - prefix_vertices);
                    }
                }
                return result;
            }
        }

        if (!initialized || candidate < best) {
            best = candidate;
            initialized = true;
        }
    }

    result.cycles = std::move(best);
    result.key = serialize(result.cycles);
    return result;
}

SymbolicGlEvaluator::NodeId SymbolicGlEvaluator::evaluate_impl(const Cycles& cycles, std::size_t depth) {
    ++stats_.calls;
    if (stats_.calls % 1000000 == 0) std::cerr << "recurrence calls=" << stats_.calls << " cache=" << cache_.size() << " dag=" << nodes_.size() << std::endl;
    stats_.max_depth = std::max(stats_.max_depth, depth);
    const auto state = normalize(cycles);
    const auto restore_empty_cycles = [this, count = state.empty_cycles](SymbolicGlEvaluator::NodeId value) {
        return count == 0U ? value : multiply(power(casimir(0U), count), value);
    };
    if (state.decomposable) {
        ++stats_.decompositions;
        const auto left_value = evaluate_impl(state.left, depth + 1U);
        if (is_zero(left_value)) {
            return zero();
        }
        return restore_empty_cycles(
            multiply(left_value, evaluate_impl(state.right, depth + 1U)));
    }
    if (state.cycles.empty()) {
        return restore_empty_cycles(one());
    }

    if (const auto found = cache_.find(state.key); found != cache_.end()) {
        ++stats_.cache_hits;
        return restore_empty_cycles(found->second);
    }

    if (state.cycles.size() == 1U && is_increasing_identity_cycle(state.cycles.front())) {
        const auto result = casimir(state.cycles.front().size());
        cache_.emplace(state.key, result);
        return restore_empty_cycles(result);
    }

    const auto& first = state.cycles.front();
    std::size_t k = 0U;
    while (k < first.size() && first[k] == k) {
        ++k;
    }
    if (k == first.size() || first[k] <= k) {
        throw std::logic_error("normalized irreducible state has no admissible recurrence position");
    }
    const auto p = first[k];
    const auto predecessor = static_cast<std::uint8_t>(p - 1U);

    std::size_t predecessor_cycle = 0U;
    std::size_t predecessor_position = 0U;
    bool located = false;
    for (std::size_t cycle_index = 0; cycle_index < state.cycles.size() && !located; ++cycle_index) {
        const auto& cycle = state.cycles[cycle_index];
        const auto iterator = std::find(cycle.begin(), cycle.end(), predecessor);
        if (iterator != cycle.end()) {
            predecessor_cycle = cycle_index;
            predecessor_position = static_cast<std::size_t>(iterator - cycle.begin());
            located = true;
        }
    }
    if (!located) {
        throw std::logic_error("predecessor label is missing from permutation");
    }

    auto swapped = state.cycles;
    for (auto& cycle : swapped) {
        for (auto& label : cycle) {
            if (label == p) {
                label = predecessor;
            } else if (label == predecessor) {
                label = p;
            }
        }
    }

    Cycles positive;
    Cycles negative;
    if (predecessor_cycle != 0U) {
        const auto alpha = rotate_left(first, k);
        const auto beta = rotate_left(state.cycles[predecessor_cycle], predecessor_position);
        for (std::size_t index = 1U; index < state.cycles.size(); ++index) {
            if (index != predecessor_cycle) {
                positive.push_back(state.cycles[index]);
                negative.push_back(state.cycles[index]);
            }
        }
        Cycle joined_positive = alpha;
        joined_positive.insert(joined_positive.end(), beta.begin() + 1, beta.end());
        Cycle joined_negative(alpha.begin() + 1, alpha.end());
        joined_negative.insert(joined_negative.end(), beta.begin(), beta.end());
        positive.insert(positive.begin(), std::move(joined_positive));
        negative.insert(negative.begin(), std::move(joined_negative));
    } else {
        if (predecessor_position <= k) {
            throw std::logic_error("recurrence predecessor must occur after the pivot");
        }
        const auto rotated = rotate_left(first, k);
        const auto distance = predecessor_position - k;
        Cycle alpha(rotated.begin(), rotated.begin() + static_cast<std::ptrdiff_t>(distance));
        Cycle beta(rotated.begin() + static_cast<std::ptrdiff_t>(distance), rotated.end());

        positive.push_back(alpha);
        positive.emplace_back(beta.begin() + 1, beta.end());
        negative.emplace_back(alpha.begin() + 1, alpha.end());
        negative.push_back(beta);
        positive.insert(positive.end(), state.cycles.begin() + 1, state.cycles.end());
        negative.insert(negative.end(), state.cycles.begin() + 1, state.cycles.end());
    }

    if (stats_.recurrence_nodes >= recurrence_deadline_) {
        throw EvaluationBudgetExceeded();
    }
    ++stats_.recurrence_nodes;
    auto result = add(evaluate_impl(swapped, depth + 1U),
                      subtract(evaluate_impl(positive, depth + 1U),
                               evaluate_impl(negative, depth + 1U)));
    cache_.emplace(state.key, result);
    return restore_empty_cycles(result);
}

} 
