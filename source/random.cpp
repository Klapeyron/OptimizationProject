#include "random.hpp"

RandomGenerator::RandomGenerator()
    :rd(),
     generator(rd())
{
}

double RandomGenerator::generateDouble(double a, double b)
{
  std::uniform_real_distribution<double> distribution(a,b);
  return distribution(generator);
}
