/* 
 * File:   random_collection.h
 * Author: fujii
 *
 * Created on 2021/2/12, 16:12
*/

#ifndef RANDOM_COLLECTION_H
#define RANDOM_COLLECTION_H

#include <random>
#include <functional>
#include <complex>
#include <iostream>

template <typename T>
    class distribution {
    public:

        T operator()() const {
            return dist_();
        }
    protected:
        std::function<T() > dist_;
    };

    template <typename T = int, typename E = std::mt19937>
    class uniform_int_distribution : public distribution<T> {
    public:

        void init(T min, T max, unsigned long int seed) {
            distribution<T>::dist_ = std::bind(std::uniform_int_distribution<T>(min, max), E(seed));
        }
    };

    template <typename T = double, typename E = std::mt19937>
    class uniform_real_distribution : public distribution<T> {
    public:

        void init(T min, T max, unsigned long int seed) {
            distribution<T>::dist_ = std::bind(std::uniform_real_distribution<T>(min, max), E(seed));
        }
    };

    template <typename T = double, typename E = std::mt19937>
    class normal_distribution : public distribution<T> {
    public:

        void init(T mean, T sd, unsigned long int seed) {
            distribution<T>::dist_ = std::bind(std::normal_distribution<T>(mean, sd), E(seed));
        }
    };

    template <typename T = double, typename E = std::mt19937>
    class cnormal_distribution : public distribution<T> {
    public:
        void init(T mean, T sd, unsigned long int seed) {
            distribution<T>::dist_ = std::bind(std::normal_distribution<T>(mean, sd), E(seed));
        }

        std::complex<T> operator()() const {
            return std::complex<T>(distribution<T>::dist_(), distribution<T>::dist_());
        }
    };

#endif /* RANDOM_COLLECTION_H */