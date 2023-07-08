#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <list>
#include <numeric>
#include <iterator>

#include "multi_modal.hpp"

using std::pair;

void test_distributions();
void test_multi_modal();
void test_many_distributions();


int main(int argc, char * argv[]) {
    // test_distributions();
    // test_multi_modal();

    test_many_distributions();

    return 0;
}


std::ostream & operator<<(std::ostream & os, std::vector<double> const & o) {
    os << "[ ";
    for (auto const & x : o) {
        os << x << " ";
    }
    return os << "]";
}

template<typename X>
double dist(X const & left, X const & right) {
    return std::abs(left - right);
}

template<typename T>
double dist(distribution<T> const & left, distribution<T> const & right) {
    return dist(left.mean, right.mean);
}

template<typename T>
double dist(T const & left_mean, distribution<T> const & right) {
    return dist(left_mean, right.mean);
}


template<typename T>
double dist(distribution<T> const & left, T const & right_mean) {
    return dist(left.mean, right_mean);
}



template<typename T>
double dist(std::vector<T> const & left, std::vector<T> const & right) {
    if (left.size() != right.size()) 
        throw std::logic_error("sizes must be equal");

    T sum = 0.;
    for(int i = 0; i < left.size(); i++) {
        T d = left[i] - right[i];
        sum += d * d;
    }

    return (double)std::sqrt(sum);
}

template<typename X>
double matched_distance(std::vector<distribution<X>> left, std::vector<distribution<X>> right) {
    // brute force for now
    
    double total = 0.;
    std::vector<distribution<X>> const * smaller = &left, * larger = &right;
    if(smaller->size() > larger->size()) {
        std::swap(smaller, larger);
    }

    int i = 0;
    for(; i < larger->size(); i++) {
        double min_distance = std::numeric_limits<double>::max();
        int min_index = -1;

        for(int j = 0; j < smaller->size(); j++) {
            double d = dist((*smaller)[i], (*larger)[j]);

            if(d < min_distance) {
                min_index = j;
                min_distance = d;
            }
        }
        // don't remove this one for now-- i dont' know how do deal with the two vectors being different sizes yet.

        total += min_distance;
    }

    return total;
}

void test_many_distributions() {
    const int iterations = 10;
    const int to_insert = 1000000;
    const int output_every = 10000;
    const int max_peaks = 4;
    const int peaks = 4;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<> means(-1, 1);
    std::uniform_real_distribution<> stds(0.000001, 0.0001);
    std::uniform_real_distribution<> weights(0.1, 1);
    
    for(int i = 0; i < iterations; i++) {
        multi_modal<double> mm(max_peaks);

        std::list<std::pair<double, std::normal_distribution<>>> dists;
        std::vector<double> my_means;
        double total_weight = 0.;
        for(int j = 0; j < peaks; j++) {
            auto w = weights(gen);
            // double w = 1.;
            total_weight += w;
            auto m = means(gen);
            auto s = stds(gen);

            dists.push_back({w, std::normal_distribution<>{m, s}});
            my_means.push_back(m);
        }
        // normalize
        for(auto & p : dists) {
            p.first /= total_weight;
        }


        for(int k = 0; k < to_insert; k++) {
            // choose a random peak:
            auto w = weights(gen);

            auto l = dists.begin();
            for(int m = 0; m < dists.size()-1; m++, l++) {
                w -= l->first;
                if(w < 0) {
                    break;
                }
            }

            // generate random number from peak
            auto val = l->second(gen);

            // insert into multi_modal
            mm.insert(val);

            if(k % output_every == 0) {
                auto peaks = mm.extract_peaks();

                std::vector<distribution<double>> p;
                std::vector<distribution<double>> q;

                std::transform(my_means.begin(), my_means.end(), std::back_inserter(q), [](auto const & tup) {
                    distribution<double> r;
                    r.mean = tup;
                    r.m2 = 1;
                    r.count = 1;
                    return r;
                });

                std::transform(peaks.begin(), peaks.end(), std::back_inserter(p), [](auto const & tup) {
                    return tup.second;
                });

                double dist = matched_distance(p, q);

                std::cout << i << "\t" << k << "\t" << dist << "\n";
            }
        }


        auto peaks = mm.extract_peaks();

        std::vector<decltype(peaks)::value_type> vp(peaks.begin(), peaks.end());
        std::vector<decltype(dists)::value_type> vd(dists.begin(), dists.end());

        std::cout << "dists\t";
        std::sort(vd.begin(), vd.end(), [](auto const & l, auto const & r) -> bool {
            return l.second.mean() < r.second.mean();
        });
        std::sort(vp.begin(), vp.end(), [](auto const & l, auto const & r) -> bool {
            return l.second.mean < r.second.mean;
        });

        std::for_each(vd.begin(), vd.end(), [](auto const & p) {
            std::cout << "\n\t(" << p.second.mean() << ";" << p.second.stddev() << ")*" << p.first;
        });
        std::cout << "\n"
                  << "found\t";

        double total = 0.;
        std::for_each(peaks.begin(), peaks.end(), [&total](auto const & p) { total += p.second.count; });

        std::for_each(vp.begin(), vp.end(), [total](auto const & p) {
            std::cout << "\n\t(" << p.second.mean << ";" << p.second.standard_deviation() << ")*" << p.second.count/total;
        });
        std::cout << "\n\n";
    }
}

void test_multi_modal() {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d11{-5, 1.5}, d12{2, 1.5};
    std::normal_distribution<> d21{0, 1.5}, d22{8, 1.5};
    std::normal_distribution<> d31{1.5, 3}, d32{-7, 3};
    std::normal_distribution<> d41{7, 2}, d42{4, 2};

    // multi_modal<float> mm(10);

    // for(int n = 0; n < 300; n++) {
    //     switch (gen() % 3) {
    //     case 0: mm.insert(d1(gen)); break;
    //     case 1: mm.insert(d2(gen)); break;
    //     case 2: mm.insert(d3(gen)); break;
    //     }
    // }

    multi_modal<std::vector<double>> mm(10);

    std::vector<double> v(2);
    std::string col = "red";

    int perframe = 100;

    std::vector<pair<std::vector<double>, std::string>> points;
  
    for(int n = 0; n < 6001; n++) {
        switch (gen() % 6) {
        case 0:
            col = "red";
            v[0] = d11(gen);
            v[1] = d12(gen);
            break;
        // make this point twice as likely
        case 1:
        case 2:
            col = "pink";
            v[0] = d21(gen);
            v[1] = d22(gen);
            break;
        case 3:
            col = "purple";
            v[0] = d31(gen);
            v[1] = d32(gen);
            break;
        case 4:
        case 5:
            col = "blue";
            v[0] = d41(gen);
            v[1] = d42(gen);
            break;
        }
        mm.insert(v);
        points.push_back({v, col});

        if (n > 5950) {
            auto found = mm.find_peak(v);
            std::cout << col << "[" << v[0] << " " << v[1] << "] found in " << found.first << "[" << found.second.mean[0] << " " << found.second.mean[1] << " / " << found.second.standard_deviation() << "]\n";
        }

        if (n > 0 && n % perframe == 0) {
            std::stringstream filename;
            filename << "example_output/frame" << std::setw(4) << std::setfill('0') << (n / perframe) << ".svg";

            std::fstream stream;
            stream.open(filename.str(), std::ios_base::out);
            stream << "<svg xmlns=\"http://www.w3.org/2000/svg\" height=\"1200\" width=\"1200\" viewbox=\"0 0 1200 1200\">\n";
            stream << "<g class=\"frame frame" << (n / perframe) << "\">\n";

            for (auto const & p : points) {
                stream << "<circle cx=\"" << p.first[0] * 50 + 500 << "\" cy=\"" << p.first[1] * 50 + 500 << "\" r=\"2\" stroke=\"" << p.second << "\" fill=\"" << p.second << "\" stroke-width=\"2\" />\n"; 
            }


            auto peaks = mm.extract_peaks();

            for(auto const & p : peaks) {
                auto d = p.second;
                auto x = d.mean[0]*50+500;
                auto y = d.mean[1]*50+500;

                stream << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << 2*d.standard_deviation()*50 << "\" fill=\"rgba(0.7,0.7,0.7,0.1)\" />\n";
                stream << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << d.standard_deviation()*50 << "\" fill=\"rgba(0.5,0.5,0.5,0.1)\" />\n";
                stream << "<text text-anchor=\"middle\" stroke-width=\"1\" alignment-baseline=\"middle\" style=\"font: normal 40px sans-serif; stroke: black; fill: white;\" x=\"" << x << "\" y=\"" << y << "\">id: " << p.first << " sig: " << d.standard_deviation() << "</text>\n";
            }

            stream << "</g>\n";
            stream << "</svg>\n";
        }
    }

    auto peaks = mm.extract_peaks();
    for(auto const & p : peaks) {
        auto d = p.second;

        std::cout << "peak " << p.first << " (" << d.mean[0] << " " << d.mean[1] << " / " << d.standard_deviation() << ")\n";
    }

    std::fstream fs;
    fs.open("serialize.multimodal", std::ios_base::trunc | std::ios_base::binary | std::ios_base::out);

    mm.serialize(fs);
    fs.close();

    fs.open("serialize.multimodal", std::ios_base::in | std::ios_base::binary);
    mm.deserialize(fs);
    fs.close();

    peaks = mm.extract_peaks();
    for(auto const & p : peaks) {
        auto d = p.second;

        std::cout << "peak " << p.first << " (" << d.mean[0] << " " << d.mean[1] << " / " << d.standard_deviation() << ")\n";
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

    // auto peaks = mm.extract_peaks();

    // for(auto const & d : peaks) {
    //     std::cout << d.mean << " / " << d.standard_deviation() << " # " << d.count << "\n";
    // }
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