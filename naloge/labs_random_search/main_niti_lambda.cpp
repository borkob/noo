#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <limits>
#include <thread>

using namespace std;
using namespace std::chrono;

class LABS {
public:
    enum value { p = +1, n = -1 };
    LABS(const size_t L) : L(L), seq(L, p), c(L, 0),
                          e(numeric_limits<int>::max()), psl(numeric_limits<int>::max()) {}
    LABS(const LABS& l) : L(l.L), seq(l.seq), c(l.c), e(l.e), psl(l.psl) {}
    LABS& operator=(const LABS& l);
    inline double get_mf() const { return (L * L) / (2.0 * e); }
    inline int get_e() const { return e; }
    inline int get_psl() const { return psl; }
    void random(mt19937& rand);
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

LABS& LABS::operator=(const LABS& l) {
    if (L != l.L) throw std::string("Sequences with different length!");
    seq = l.seq;
    c = l.c;
    e = l.e;
    psl = l.psl;
    return *this;
}

void LABS::random(mt19937& rand) {
    for (size_t i = 0; i < L; i++) {
        if (rand() % 2)
            seq[i] = p;
        else
            seq[i] = n;
    }
}

void LABS::evaluate_e() {
    e = 0;
    for (size_t k = 1; k < L; k++) {
        c[k] = 0;
        for (size_t i = 0; i <= L - k - 1; i++)
            c[k] += seq[i] * seq[i + k];
        e += c[k] * c[k];
    }
}

void LABS::evaluate_psl() {
    psl = 0;
    for (size_t k = 1; k < L; k++) {
        c[k] = 0;
        for (size_t i = 0; i <= L - k - 1; i++)
            c[k] += seq[i] * seq[i + k];
        if (abs(c[k]) > psl)
            psl = abs(c[k]);
    }
}

LABS LABS::search_e(const size_t seed, const size_t n, const size_t L) {
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_e();
    for (size_t i = 0; i < n; i++) {
        current.random(rand);
        current.evaluate_e();
        if (current.get_e() < best.get_e())
            best = current;
    }
    return best;
}

LABS LABS::search_psl(const size_t seed, const size_t n, const size_t L) {
    LABS current(L), best(L);
    mt19937 rand(seed);
    best.random(rand);
    best.evaluate_psl();
    for (size_t i = 0; i < n; i++) {
        current.random(rand);
        current.evaluate_psl();
        if (current.get_psl() < best.get_psl())
            best = current;
    }
    return best;
}

void performSearch(const size_t seed, const size_t n, const size_t L, LABS (*searchFunction)(const size_t, const size_t, const size_t)) {
    LABS best(L);
    mt19937 rand(seed);
    best.random(rand);
    auto start = system_clock::now();
    vector<thread> threads;
    for (size_t i = 0; i < n; i++) {
        threads.emplace_back([&best, &rand, searchFunction]() {
            LABS current(best);
            current.random(rand);
            if (searchFunction(current.get_e(), best.get_e(), best.get_psl(), current.get_psl()))
                best = current;
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    auto end = system_clock::now();
    auto elapsed = duration_cast<milliseconds>(end - start);
    if (searchFunction == &LABS::search_e) {
        cout << "E: " << best.get_e() << " F: " << best.get_mf();
    } else if (searchFunction == &LABS::search_psl) {
        cout << "PSL: " << best.get_psl();
    }
    cout << " speed: " << n / (elapsed.count() / 1000.0) << " eval/sec" << endl;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 4)
            throw string("Three arguments are required!");

        const size_t seed = atoi(argv[1]), n = atoi(argv[2]), D = atoi(argv[3]);
        cout << "Searching ..." << endl;
        performSearch(seed, n, D, &LABS::search_e);

        cout << "Searching ..." << endl;
        performSearch(seed, n, D, &LABS::search_psl);
    } catch (string err) {
        cerr << err << std::endl;
        return 1;
    }
    return 0;
}


/*
Added a performSearch function that takes a function pointer to the search function (search_e or search_psl) and performs the multithreaded search.
performSearch initializes the best variable outside the thread loop to avoid data races and passes it by reference to the threads.
Inside each thread, a copy of best is made (current) to perform the search.
The appropriate search function (search_e or search_psl) is called based on the function pointer passed to performSearch.
The thread objects are stored in a vector to ensure all threads are joined before moving on.
The timing and result output are done after all the threads have finished executing.*/
