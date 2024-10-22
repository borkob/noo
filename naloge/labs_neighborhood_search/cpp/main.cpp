#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <limits>

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
    int neighbor_e(const size_t n) const;
    int neighbor_psl(const size_t n) const;
    void update_e(const size_t n, const int e);
    void update_psl(const size_t n, const int psl);
    static LABS neighborhood_search_e(const size_t seed, const size_t n, const size_t L);
    static LABS neighborhood_search_psl(const size_t seed, const size_t n, const size_t L);

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

int LABS::neighbor_e(const size_t i) const{
    int e = 0, ck;
    const size_t lmt = std::max(L-i,i+1);
    size_t k=1;
    for(;k<lmt; k++) {
        ck = c[k];
        if(i+k<L) ck -= 2*seq[i]*seq[k+i];
        if(k<=i) ck -= 2*seq[i-k]*seq[i];
        e += ck*ck;
    }
    for (; k<L; k++) e += c[k]*c[k];
    return e;
}

int LABS::neighbor_psl(const size_t i) const{
    int psl = 0, ck;
    const size_t lmt = max(L-i,i+1);
    size_t k=1;
    for(; k<lmt; k++) {
        ck = c[k];
        if(i+k<L) ck -= 2*seq[i]*seq[k+i];
        if(k<=i) ck -= 2*seq[i-k]*seq[i];
        if(abs(ck) > psl) psl = abs(ck);
    }
    for (size_t k=1; k<L; k++){
        if(abs(c[k]) > psl) psl = abs(c[k]);
    }
    return psl;
}

void LABS::update_e(const size_t i, const int e){
    const size_t lmt = max(L-i,i+1);
    size_t k=1;
    for (;k<lmt; k++){
        int ck = c[k];
        if(i+k<L) ck -= 2*seq[i]*seq[k+i];
        if(k<=i) ck -= 2*seq[i-k]*seq[i];
        c[k] = ck;
    }
    this->e = e;
    seq[i] = (value)(-seq[i]);
    #ifndef NDEBUG
    int update_e = e;
    evaluate_e();
    if(e != update_e) throw string("Wrong E!");
    #endif
}

void LABS::update_psl(const size_t i, const int psl){
    const size_t lmt = max(L-i,i+1);
    size_t k=1;
    int ck;
    for (; k<lmt; k++){
        ck = c[k];
        if(i+k<L) ck -= 2*seq[i]*seq[k+i];
        if(k<=i) ck -= 2*seq[i-k]*seq[i];
        c[k] = ck;
    }
    this->psl = psl;
    seq[i] = (value)(-seq[i]);
    #ifndef NDEBUG
    int update_psl = psl;
    evaluate_psl();
    if(psl != update_psl) throw string("Wrong PSL!");
    #endif
}

LABS LABS::neighborhood_search_e(const size_t seed, const size_t n, const size_t L){
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_e();
    current = best;
    size_t nfes=0, best_neighbor;
    int best_neighbor_e, e;
    while(nfes < n){
        best_neighbor_e = numeric_limits<int>::max();
        for(size_t i=0; i<L; i++){
            e = current.neighbor_e(i);
            if(e < best_neighbor_e){
                best_neighbor = i;
                best_neighbor_e = e;
            }
        }
        nfes+=L;
        if(best_neighbor_e >= current.get_e()){
            current.random(rand);
            current.evaluate_e();
            nfes++;
        }
        else{
            current.update_e(best_neighbor,best_neighbor_e);
        }
        if(current.get_e() < best.get_e()) best = current;
    }
    return best;
}

LABS LABS::neighborhood_search_psl(const size_t seed, const size_t n, const size_t L){
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_psl();
    current = best;
    size_t nfes=0, best_neighbor;
    int best_neighbor_psl, psl;
    while(nfes < n){
        best_neighbor_psl = numeric_limits<int>::max();
        for(size_t i=0; i<L; i++){
            psl = current.neighbor_psl(i);
            if(psl < best_neighbor_psl){
                best_neighbor = i;
                best_neighbor_psl = psl;
            }
        }
        nfes+=L;
        if(best_neighbor_psl >= current.get_psl()){
            current.random(rand);
            current.evaluate_psl();
            nfes++;
        }
        else{
            current.update_psl(best_neighbor,best_neighbor_psl);
        }
        if(current.get_psl() < best.get_psl()) best = current;
    }
    return best;
}

int main(int argc, char *argv[]){
    try{
        if(argc < 4) throw string("Three arguments are required!");

        const size_t seed =atoi(argv[1]), n = atoi(argv[2]), D = atoi(argv[3]);
        cout<<"Searching ..."<<endl;
        auto start = system_clock::now();
        LABS best = LABS::neighborhood_search_e(seed,n,D);
        auto end = system_clock::now();
        auto elapsed = duration_cast<milliseconds>(end - start);
        cout<<"E: "<<best.get_e()<<" F: "<<best.get_mf() << endl;
		  best.evaluate_e();
		  cout << "evaluate_e(): " << best.get_e() << endl;
        cout<<" speed: "<<n/(elapsed.count()/1000.0)<<" eval/sec"<<endl;

        cout<<"Searching ..."<<endl;
        start = system_clock::now();
        best = LABS::neighborhood_search_psl(seed,n,D);
        end = system_clock::now();
        elapsed = duration_cast<milliseconds>(end - start);
        cout<<"PSL: "<<best.get_psl();
        cout<<" speed: "<<n/(elapsed.count()/1000.0)<<" eval/sec"<<endl;
    }
    catch (string err) {
        cerr<<err<<std::endl;
        return 1;
    }

    return 0;
}
