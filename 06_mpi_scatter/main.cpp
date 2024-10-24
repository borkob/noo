#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <limits>
#include "mpi.h"

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

void slave(const int rank, size_t seed, const size_t NFEs, const size_t D){
    const int tag_e=1, tag_psl=2;
    auto start = system_clock::now();
    LABS best = LABS::neighborhood_search_e(seed,NFEs,D);
    auto end = system_clock::now();
    auto elapsed = duration_cast<milliseconds>(end - start);
    int best_e = best.get_e();
    double speed = NFEs/(elapsed.count()/1000.0);
    MPI_Send(&best_e,1, MPI_INT, 0, tag_e, MPI_COMM_WORLD); // Posljemo E z oznako 1
    MPI_Send(&speed,1, MPI_DOUBLE, 0, tag_e, MPI_COMM_WORLD); // Posljemo hitrost z oznako 1

    MPI_Barrier(MPI_COMM_WORLD); // Cakamo na ostale procese

    start = system_clock::now();
    best = LABS::neighborhood_search_psl(seed,NFEs,D);
    end = system_clock::now();
    elapsed = duration_cast<milliseconds>(end - start);
    int best_psl = best.get_psl();
    speed = NFEs/(elapsed.count()/1000.0);
    MPI_Send(&best_psl,1, MPI_INT, 0, tag_psl, MPI_COMM_WORLD); // Posljemo PSL z oznako 2
    MPI_Send(&speed,1, MPI_DOUBLE, 0, tag_psl, MPI_COMM_WORLD); // Posljemo hitrost z oznako 2
}

void master(const size_t size){
    MPI_Status status;
    int best_e, best_psl;
    double speed;
    const int tag_e=1, tag_psl=2;
    std::cout<<"F"<<std::endl;
    for(size_t i=1; i<size; i++){
        // Prejmemo E z oznako 1
        MPI_Recv(&best_e, 1, MPI_INT, i, tag_e, MPI_COMM_WORLD, &status);
        // Prejmemo hitrost z oznako 1
        MPI_Recv(&speed, 1, MPI_DOUBLE, i, tag_e, MPI_COMM_WORLD, &status);
        std::cout<<"Slave "<<i<<" E: "<<best_e<<" speed: "<<speed<<" eval/sec"<<std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD); // Cakamo na ostale procese
    std::cout<<"PSL"<<std::endl;

    for(size_t i=1; i<size; i++){
        MPI_Recv(&best_psl, 1, MPI_INT, i, tag_psl, MPI_COMM_WORLD, &status);
        // Prejmemo PSL z oznako 2
        MPI_Recv(&speed, 1, MPI_DOUBLE, i, tag_psl, MPI_COMM_WORLD, &status);
        // Prejmemo hitrost z oznako 2
        std::cout<<"Slave "<<i<<" PSL: "<<best_psl<<" speed: "<<speed<<" eval/sec"<<std::endl;
    }
}

int main(int argc, char *argv[]){
    int size, rank;
    size_t NFEs, D, my_seed;
    std::vector<size_t> seed;
    MPI_Init(&argc, &argv);  // Inicializacija okolja MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Stevilo procesov
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Rank procesa
    seed.resize(size);
    try{
        if(rank == 0){
            if(argc < 4) throw string("Three arguments are required: seed NFEs D!");
            seed[0] =atoi(argv[1]);
            NFEs = atoi(argv[2]);
            D = atoi(argv[3]);
            for(size_t i=1; i<size; i++) seed[i] = seed[i-1]+1;
        }
        MPI_Scatter(&seed[0],1,MPI_UNSIGNED_LONG,&my_seed,1,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
        MPI_Bcast(&NFEs, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&D, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        if(rank == 0) master(size); // Gospodar - zbiratelj informacij
        else slave(rank,my_seed,NFEs,D); // Suznji
    }
    catch (string err) {
            cerr<<err<<std::endl;
            return 1;
    }
    MPI_Finalize(); // Koncamo okolje MPI
    return 0;
}
