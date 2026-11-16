#include "chord/symbolic_gl.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

chord::Word parse(const std::string& s) {
    std::istringstream in(s);unsigned x;chord::Word w;
    while(in>>x){if(x>255)throw std::invalid_argument("label exceeds 255");w.push_back(x);}
    if(!in.eof())throw std::invalid_argument("word contains a non-integer token");
    if(w.empty())throw std::invalid_argument("empty word");return w;
}
void save(const std::filesystem::path& dir,const std::string& name,const chord::ExactPolynomial& p) {
    std::filesystem::create_directories(dir);
    for(const auto& ext:{"txt","json"}) {
        auto path=dir/(name+"."+ext);auto temporary=path;temporary+=".partial";
        if(std::filesystem::exists(path))throw std::runtime_error("refusing to overwrite "+path.string());
        std::ofstream f(temporary);f<<(std::string(ext)=="txt"?chord::polynomial_string(p):chord::polynomial_json(p))<<"\n";
        f.close();if(!f)throw std::runtime_error("output write failed");
        std::filesystem::rename(temporary,path);
    }
}
int main(int argc,char** argv) {
    try {
        if(argc>2)throw std::invalid_argument("usage: symbolic_witness [output_directory | --word]");
        if(argc==2&&(std::string(argv[1])=="--word"||std::string(argv[1])=="--text")){
            std::string s;std::getline(std::cin,s);chord::SymbolicGlEvaluator e;
            auto poly=e.evaluate_chord_word(parse(s));
            std::cout<<(std::string(argv[1])=="--text"?chord::polynomial_string(poly):chord::polynomial_json(poly))<<"\n";return 0;
        }
        std::filesystem::path dir=argc>1?argv[1]:"data/symbolic_witness";
        for(const auto& name:{"D","D_prime","difference"})for(const auto& ext:{"txt","json"}){
            auto path=dir/(std::string(name)+"."+ext);
            if(std::filesystem::exists(path))throw std::runtime_error("refusing to overwrite "+path.string());
        }
        const std::string words[2]={
            "0 1 2 0 3 1 2 4 5 6 7 8 5 9 4 8 10 3 7 9 10 6",
            "3 4 5 6 7 8 5 9 4 8 10 0 1 2 0 3 1 2 7 9 10 6"};
        chord::ExactPolynomial results[2];
        for(unsigned i=0;i<2;++i){
            auto begin=std::chrono::steady_clock::now();chord::SymbolicGlEvaluator e;
            std::cerr<<"Starting "<<(i==0?"D":"D_prime")<<std::endl;
            results[i]=e.evaluate_chord_word(parse(words[i]));
            save(dir,i==0?"D":"D_prime",results[i]);
            std::cerr<<"Saved "<<i<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<std::endl;
        }
        auto difference=chord::polynomial_sum(results[0],results[1],-1);save(dir,"difference",difference);
        std::cout<<"D terms="<<results[0].size()<<" D_prime terms="<<results[1].size()<<" difference terms="<<difference.size()<<std::endl;
        return difference.empty()?2:0;
    }catch(const std::exception& e){std::cerr<<"error: "<<e.what()<<std::endl;return 1;}
}
