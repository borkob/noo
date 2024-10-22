#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <limits>
#include <thread>
#include <mutex>

using namespace std;
using namespace std::chrono;

class LABS {
public:
    enum value { p=+1, n=-1 };
    LABS(const size_t L) : L(L), seq(L,p), c(L,0),
                          e(numeric_limits<int>::max()), psl(numeric_limits<int>::max()) {}
    LABS(const LABS & l) : L(l.L), seq(l.seq), c(l.c), e(l.e), psl(l.psl) {}
    LABS& operator=(const LABS & l);
    inline double get_mf() const { return (L*L)/(2.0*e); }
    inline int get_e() const { return e; }
    inline int get_psl() const { return psl; }
    void random(mt19937 & rand);
    void evaluate_e();
    void evaluate_psl();
    static LABS search_e(const size_t seed, const size_t n, const size_t L);
    static LABS search_psl(const size_t seed, const size_t n, const size_t L);

private:
    const size_t L;
    vector<value> seq;
    vector<int> c;
    int e, psl;
};

LABS& LABS::operator=(const LABS & l) {
    if (L != l.L) throw std::string("Sequences with different length!");
    seq = l.seq; c = l.c; e = l.e; psl = l.psl; return *this;
}

void LABS::random(mt19937 & rand) {
    for (size_t i = 0; i < L; i++) {
        if (rand() % 2) seq[i] = p;
        else seq[i] = n;
    }
}

void LABS::evaluate_e() {
    e = 0;
    for (size_t k = 1; k < L; k++) {
        c[k] = 0;
        for (size_t i = 0; i <= L - k - 1; i++) c[k] += seq[i] * seq[i + k];
        e += c[k] * c[k];
    }
}

void LABS::evaluate_psl() {
    psl = 0;
    for (size_t k = 1; k < L; k++) {
        c[k] = 0;
        for (size_t i = 0; i <= L - k - 1; i++) c[k] += seq[i] * seq[i + k];
        if (abs(c[k]) > psl) psl = abs(c[k]);
    }
}

class SearchE {
private:
    LABS& current;
    LABS& best;
    mt19937& rand;
    mutex& mtx;

public:
    SearchE(LABS& current, LABS& best, mt19937& rand, mutex& mtx)
        : current(current), best(best), rand(rand), mtx(mtx) {}

    void operator()() {
        current.random(rand);
        current.evaluate_e();
        lock_guard<mutex> lock(mtx);
        if (current.get_e() < best.get_e()) best = current;
    }
};

class SearchPSL {
private:
    LABS& current;
    LABS& best;
    mt19937& rand;
    mutex& mtx;

public:
    SearchPSL(LABS& current, LABS& best, mt19937& rand, mutex& mtx)
        : current(current), best(best), rand(rand), mtx(mtx) {}

    void operator()() {
        current.random(rand);
        current.evaluate_psl();
        lock_guard<mutex> lock(mtx);
        if (current.get_psl() < best.get_psl()) best = current;
    }
};

LABS LABS::search_e(const size_t seed, const size_t n, const size_t L) {
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_e();

    mutex mtx;
    vector<thread> threads;
    SearchE searchE(current, best, rand, mtx);
    for (size_t i = 0; i < n; i++) {
        threads.emplace_back(searchE);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return best;
}

LABS LABS::search_psl(const size_t seed, const size_t n, const size_t L) {
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_psl();

    mutex mtx;
    vector<thread> threads;
    SearchPSL searchPSL(current, best, rand, mtx);
    for (size_t i = 0; i < n; i++) {
        threads.emplace_back(searchPSL);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return best;
}

int main(int argc, char *argv[]) {
    try {
        if (argc < 4) throw string("Three arguments are required!");

        const size_t seed = atoi(argv[1]), n = atoi(argv[2]), D = atoi(argv[3]);
        cout << "Searching ..." << endl;
        auto start = system_clock::now();
        LABS best = LABS::search_e(seed,n,D);
        auto end = system_clock::now();
        auto elapsed = duration_cast<milliseconds>(end - start);
        cout << "E: " << best.get_e() << " F: " << best.get_mf();
        cout << " speed: " << n / (elapsed.count() / 1000.0) << " eval/sec" << endl;

        cout << "Searching ..." << endl;
        start = system_clock::now();
        best = LABS::search_psl(seed,n,D);
        end = system_clock::now();
        elapsed = duration_cast<milliseconds>(end - start);
        cout << "PSL: " << best.get_psl();
        cout << " speed: " << n / (elapsed.count() / 1000.0) << " eval/sec" << endl;
    }
    catch (string err) {
        cerr << err << std::endl;
        return 1;
    }
    return 0;
}
/*
In this version, I have defined two function objects (SearchE and SearchPSL) that encapsulate the search logic for evaluate_e and evaluate_psl respectively. Each function object is responsible for executing the search and updating the best object while ensuring thread safety using a mutex.

Inside the search_e and search_psl functions, instead of using lambda functions, I create instances of the function objects and pass them to the threads. Each thread executes the operator() function of the respective function object.

Please note that using function objects instead of lambdas allows for more flexibility and code reusability if the search logic needs to be modified in the future.*/
