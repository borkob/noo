#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <limits>
#include <thread>

using namespace std;
using namespace std::chrono;

class LABS{
public:
    enum value { p=+1, n=-1};
    LABS(const size_t L): L(L), seq(L,p), c(L,0),
          e(numeric_limits<int>::max()), psl(numeric_limits<int>::max()) {};
    LABS(const LABS & l): L(l.L), seq(l.seq), c(l.c), e(l.e), psl(l.psl) {}
    LABS& operator=(const LABS & l);
    inline double get_mf() const { return (L*L)/(2.0*e); }
    inline int get_e() const { return e; }
    inline int get_psl() const { return psl; }
    void random(mt19937 & rand);
    void evaluate_e();
    void evaluate_psl();
    static LABS search_e(const size_t seed, const size_t n, const size_t L, const size_t num_threads);
    static LABS search_psl(const size_t seed, const size_t n, const size_t L, const size_t num_threads);

private:
    const size_t L;
    vector<value> seq;
    vector<int> c;
    int e, psl;
};

LABS& LABS::operator=(const LABS & l){
    if(L != l.L) throw std::string("Sequences with diferent length!");
    seq=l.seq; c=l.c; e=l.e; psl=l.psl; return *this;
}

void LABS::random(mt19937 & rand){
    for(size_t i=0; i<L; i++){
        if(rand()%2) seq[i] = p;
        else seq[i] = n;
    }
}

void LABS::evaluate_e(){
    e = 0;
    for (size_t k=1; k<L; k++) {
        c[k]=0;
        for (size_t i=0; i<=L-k-1; i++) c[k] += seq[i]*seq[i+k];
        e += c[k]*c[k];
    }
}

void LABS::evaluate_psl(){
    psl = 0;
    for (size_t k=1; k<L; k++) {
        c[k]=0;
        for (size_t i=0; i<=L-k-1; i++) c[k] += seq[i]*seq[i+k];
        if(abs(c[k]) > psl) psl = abs(c[k]);
    }
}

LABS LABS::search_e(const size_t seed, const size_t n, const size_t L, const size_t num_threads){
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_e();
    for(size_t i=0; i<n; i++){
        current.random(rand);
        current.evaluate_e();
        if(current.get_e() < best.get_e()) best = current;
    }
    return best;
}

LABS LABS::search_psl(const size_t seed, const size_t n, const size_t L, const size_t num_threads){
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_psl();
    for(size_t i=0; i<n; i++){
        current.random(rand);
        current.evaluate_psl();
        if(current.get_psl() < best.get_psl()) best = current;
    }
    return best;
}

int main(int argc, char *argv[]){
    try{
        if(argc < 4) throw string("Four arguments are required!");

        const size_t seed =atoi(argv[1]), n = atoi(argv[2]), D = atoi(argv[3], num_threads = atoi(argv[4]);

        cout<<"Searching ..."<<endl;
        auto start = system_clock::now();

        size_t num_threads;
        std::vector<std::thread> niti;

        for (int i = 0; i < num_threads; i++) {
            niti.emplace_back(LABS::search_e, seed+i, n, D, num_threads);
        }

        LABS best = LABS::search_e(seed, n/num_threads, D, num_threads);
        auto end = system_clock::now();
        auto elapsed = duration_cast<milliseconds>(end - start);
        cout<<"E: "<<best.get_e()<<" F: "<<best.get_mf();
        cout<<" speed: "<<n/(elapsed.count()/1000.0)<<" eval/sec"<<endl;

        for (std::thread& t : niti) {
            t.join();

        cout<<"Searching ..."<<endl;
        start = system_clock::now();

        size_t num_threads; // Poljubno število niti
        std::vector<std::thread> niti;

        for (int i = 0; i < num_threads; i++) {
            niti.emplace_back(LABS::search_psl, seed + i, n, D, num_threads);
        }

        best = LABS::search_psl(seed, n / num_threads, D, num_threads);
        end = system_clock::now();
        elapsed = duration_cast<milliseconds>(end - start);
        cout<<"PSL: "<<best.get_psl();
        cout<<" speed: "<<n/(elapsed.count()/1000.0)<<" eval/sec"<<endl;

        for (std::thread& t : niti) {
            t.join();
        }
        }
    }
    catch (string err) {
        cerr<<err<<std::endl;
        return 1;
    }
    return 0;
}
