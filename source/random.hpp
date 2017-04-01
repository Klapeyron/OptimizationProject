#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

class RandomGenerator
{
    std::random_device rd;
    std::default_random_engine generator;
public:
    RandomGenerator();
    double generateDouble(double a, double b);
};

#endif // RANDOM_HPP
