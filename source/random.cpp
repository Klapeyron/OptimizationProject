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

unsigned RandomGenerator::generateInt(unsigned a, unsigned b)
{
  std::uniform_int_distribution<unsigned> distribution(a,b);
  return distribution(generator);
}
