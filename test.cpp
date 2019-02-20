#include <iostream>
#include <iomanip>
#include <random>

#include "multi_modal.hpp"

void test_distributions();
void test_multi_modal();


int main(int argc, char * argv[]) {
    // test_distributions();
    test_multi_modal();

    return 0;
}


std::ostream & operator<<(std::ostream & os, std::vector<double> const & o) {
    os << "[ ";
    for (auto const & x : o) {
        os << x << " ";
    }
    return os << "]";
}

void test_multi_modal() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d11{-5, 1.5};
    std::normal_distribution<> d21{3, 1.5};
    std::normal_distribution<> d31{8, 1};
    std::normal_distribution<> d12{2, 1.5};
    std::normal_distribution<> d22{0, 1.5};
    std::normal_distribution<> d32{-1, 1};

    // multi_modal<float> mm(10);

    // for(int n = 0; n < 300; n++) {
    //     switch (gen() % 3) {
    //     case 0: mm.insert(d1(gen)); break;
    //     case 1: mm.insert(d2(gen)); break;
    //     case 2: mm.insert(d3(gen)); break;
    //     }
    // }

    multi_modal<std::vector<double>> mm(100);

    std::vector<double> v(2);
    for(int n = 0; n < 30000; n++) {
        switch (gen() % 3) {
        case 0:
            v[0] = d11(gen);
            v[1] = d12(gen);
            break;
        case 1:
            v[0] = d21(gen);
            v[1] = d22(gen);
            break;
        case 2:
            v[0] = d31(gen);
            v[1] = d32(gen);
            break;
        }
        mm.insert(v);
    }


    // mm.visit_children([](auto const & parent, auto const & left, auto const & right, size_t depth) {
    //     for(int d = (int)depth; d > 0; d-=1) std::cout << " ";

    //     std::cout << parent.mean << " / " << parent.standard_deviation() << " # " << parent.count;

    //     auto err = mixture_error(left, right);
    //     std::cout << " error: " << err << "\n";
    //     return true;
    // });

    // mm.visit_children([](auto const & parent, auto const & left, auto const & right, size_t depth) {
    //     for(int d = (int)depth; d > 0; d-=1) std::cout << " ";

    //     std::cout << parent.mean << " / " << parent.standard_deviation() << " # " << parent.count;

    //     auto err = mixture_error(left, right);
    //     std::cout << " error: " << err << "\n";
    //     return true;
    // });

    auto peaks = mm.extract_peaks();

    for(auto const & d : peaks) {
        std::cout << d.mean << " / " << d.standard_deviation() << " # " << d.count << "\n";
    }
}



void test_distributions() {
    auto a = distribution<float>::from_standard_deviation(0.0, 1.0, 1000);
    auto b = distribution<float>::from_standard_deviation(0.0, 2.0, 1000);
    auto c = distribution<float>::from_standard_deviation(1.0, 1.0, 1000);
    auto d = distribution<float>::from_standard_deviation(1.0, 1.0, 2000);
    auto e = distribution<float>::from_standard_deviation(4., 1.0, 1000);
    auto f = distribution<float>::from_standard_deviation(-4., 1.0, 1000);

    auto ab = mix(a, b), ac = mix(a, c), ad = mix(a, d), aa = mix(a,a), ae = mix(a,e), ef = mix(e, f);

    auto aba = unmix(ab, b);
    auto aca = unmix(ac, c);

    std::cout << std::setw(5) << "a:" << a << "\n";
    std::cout << std::setw(5) << "b:"  << b << "\n";
    std::cout << std::setw(5) << "c:"  << c << "\n";
    std::cout << std::setw(5) << "d:"  << d << "\n";
    std::cout << std::setw(5) << "e:"  << e << "\n";
    std::cout << std::setw(5) << "f:"  << f << "\n";
    std::cout << std::setw(5) << "ab:"  << ab << "\n";
    std::cout << std::setw(5) << "ac:" << ac << "\n";
    std::cout << std::setw(5) << "ad:" << ad << "\n";
    std::cout << std::setw(5) << "aa:" << aa << "\n";
    std::cout << std::setw(5) << "aba:" << aba << "\n";
    std::cout << std::setw(5) << "aca:" << aca << "\n";
    std::cout << std::setw(5) << "ae:" << ae << "\n";
    std::cout << std::setw(5) << "ef:" << ef << "\n";

    std::cout << std::setw(10) << "err ab:" << mixture_error2(a,b) << "\n";
    std::cout << std::setw(10) << "err ac:" << mixture_error2(a,c) << "\n";
    std::cout << std::setw(10) << "err ad:" << mixture_error2(a,d) << "\n";
    std::cout << std::setw(10) << "err bc:" << mixture_error2(b,c) << "\n";
    std::cout << std::setw(10) << "err bd:" << mixture_error2(b,d) << "\n";
    std::cout << std::setw(10) << "err cd:" << mixture_error2(c,d) << "\n";
    std::cout << std::setw(10) << "err aa:" << mixture_error2(a,a) << "\n";
    std::cout << std::setw(10) << "err ae:" << mixture_error2(a,e) << "\n";
    std::cout << std::setw(10) << "err ef:" << mixture_error2(e,f) << "\n";
}