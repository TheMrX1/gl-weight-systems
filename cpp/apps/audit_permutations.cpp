
#include "chord/gl_weight.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
int main(){
    chord::GlEvaluator evaluator(0x4155444954ULL);
    std::cout<<"{\"primes\":[";
    for(unsigned lane=0;lane<3;++lane){if(lane)std::cout<<',';std::cout<<evaluator.primes()[lane];}
    std::cout<<"],\"casimirs\":[";
    for(unsigned lane=0;lane<3;++lane){if(lane)std::cout<<',';std::cout<<'[';
        for(unsigned k=0;k<23;++k){if(k)std::cout<<',';std::cout<<evaluator.casimir_values()[lane][k];}std::cout<<']';}
    std::cout<<"]}\n";
    for(unsigned n=0;n<=7;++n){
        std::vector<unsigned> sigma(n);std::iota(sigma.begin(),sigma.end(),0);
        do{
            std::vector<bool> visited(n);chord::Cycles cycles;
            for(unsigned i=0;i<n;++i)if(!visited[i]){
                chord::Cycle c;unsigned x=i;
                do{c.push_back(static_cast<std::uint8_t>(x));visited[x]=true;x=sigma[x];}while(x!=i);
                cycles.push_back(c);
            }
            const auto value=evaluator.evaluate_cycles(cycles);
            std::cout<<"{\"permutation\":[";
            for(unsigned i=0;i<n;++i){if(i)std::cout<<',';std::cout<<sigma[i];}
            std::cout<<"],\"values\":[";
            for(unsigned lane=0;lane<3;++lane){if(lane)std::cout<<',';std::cout<<value.value[lane];}
            std::cout<<"]}\n";
        }while(std::next_permutation(sigma.begin(),sigma.end()));
    }
}
