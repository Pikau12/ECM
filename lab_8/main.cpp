#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <immintrin.h>
#include <limits>
#include <algorithm>
#include <cstdint>

const int K = 10;
const int N_MIN = 256;                 //1Kb
const int N_MAX = 4 * 1024 * 1024;    //16Mb
const int COUNT = 100;

void forward(std::vector<int>& arr) {
    for (size_t i = 0; i < arr.size() - 1; ++i) {
        arr[i] = i + 1;
    }
    arr.back() = 0;
}

void reverse(std::vector<int>& arr) {
    for (size_t i = arr.size() - 1; i > 0; --i) {
        arr[i] = i - 1;
    }
    arr[0] = arr.size() - 1;
}

void swap(int& a, int& b) {
    int c = a;
    a = b;
    b = c;
}

void random_(std::vector<int>& arr) {
    for (size_t i = 0; i < arr.size(); ++i) {
        arr[i] = i;
    }

    for (size_t i = arr.size() - 1; i > 0; --i) {
        swap(arr[i], arr[rand() % (i + 1)]);
    }

    int current = arr[0];
    for (size_t i = 1; i < arr.size(); ++i) {
        int next = arr[i];
        arr[current] = next;
        current = next;
    }
    arr[current] = arr[0];
}

long double tacts(const std::vector<int>& arr) {
    long double min_avg = std::numeric_limits<long double>::max();

    for (int it = 0; it < COUNT; ++it) {
        volatile int x = 0;
        for (size_t i = 0; i < arr.size(); ++i) {
            x = arr[x];
        }

        uint64_t start, end;
        start = __rdtsc();

        for (int i = 0; i < arr.size() * K; ++i) {
            x = arr[x];
        }

        end = __rdtsc();

        long double avg = (long double)(end - start) / (arr.size() * K);
        min_avg = std::min(min_avg, avg);
    }
    return min_avg;
}

int main() {
    srand(time(0));

    std::ofstream out("resultO1.csv");
    out << "Size;Forward;Reverse;Random\n";

    for (int N = N_MIN; N <= N_MAX; N = static_cast<int>(N * 1.2)) {
        std::vector<int> arr(N);

        forward(arr);
        long double time_forward = tacts(arr);

        reverse(arr);
        long double time_reverse = tacts(arr);

        random_(arr);
        long double time_random = tacts(arr);

        out << N * sizeof(int) << ";" << time_forward << ";" << time_reverse << ";" << time_random << "\n";
    }
    out.close();
    return 0;
}
