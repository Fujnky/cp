//CP Blatt 2, Aufgabe 1. Dag-Björn Hering und Lars Funke
//Linear Congruential Generator & XORShift

#include <iostream>
#include <fstream>
#include <limits>

//Random Number Generator class for testing
class RNG
{
  public:
    virtual double getRndDouble() = 0;
};


//XORShift generator (can use any int type to work with)
template <typename T>
class XORShift : public RNG
{
public:
  //constructor, initialize values
  XORShift(T x, T a, T b, T c) : x(x), a(a), b(b), c(c) {};

  //XORShift implementation
  T getRnd()
  {
    x ^= (x << a);
    x ^= (x >> b);
    x ^= (x << c);
    return x;
  }

  //convert random number to double in interval [0,1]
  double getRndDouble()
  {
    return getRnd() / (double) std::numeric_limits<T>::max();
  }

private:
  //state
  T x, a, b, c;
};

//Linear congruential generator
class LCG : public RNG
{
public:
  //constructor, initialize values
  LCG(uint32_t r0, uint32_t a, uint32_t c, uint32_t m) : r(r0), a(a), c(c), m(m) {};

  //LCG implementation
  uint32_t getRnd()
  {
    r = (a * r + c) % m;
    return r;
  };

  //convert random number to double in interval [0,1]
  double getRndDouble()
  {
    return getRnd() / (double) m;
  };

private:
  //state
  uint32_t r, a, c, m;
};

//find out period of given XORShift generator
template <typename T>
int testPeriod(XORShift<T>& rng)
{
  T first = rng.getRnd();
  T i = 1;
  //try until the same number occurs again
  for(; rng.getRnd() != first; i++);
  return i;
}


//sample RNG
void testRNG(RNG& rng, size_t samples, std::string filename)
{
  std::ofstream f;  //file to save results to
	f.open(filename); //file to save results to

  for(size_t i = 0; i < samples; i++)
    f << rng.getRndDouble() << std::endl;
  f.close();
}

int main()
{
  //initialize RNGs with different parameteres
  LCG lcg1(1234,      20,    120, 6075);
  LCG lcg2(1234,      137,   187, 256);
  LCG lcg3(123456789, 65539, 0,   2147483648);
  LCG lcg4(1234,      16807, 0,   2147483647);
  XORShift<uint16_t> xs1(123, 11, 1, 7);
  XORShift<uint16_t> xs2(123, 11, 4, 7);

  //save samples
  testRNG(lcg1, 1e5, "build/1.txt");
  testRNG(lcg2, 1e5, "build/2.txt");
  testRNG(lcg3, 1e5, "build/3.txt");
  testRNG(lcg4, 1e5, "build/4.txt");
  testRNG(xs1,  1e5, "build/5.txt");
  testRNG(xs2,  1e5, "build/6.txt");

  //test XORShift periods
  std::ofstream f;            //file to save results to
	f.open("build/period.txt"); //file to save results to
  for(int b = 1; b <= 15; b++)
  {
    for(int c = 1; c <= 15; c++)
    {
      //for every combination of b & c in interval [1,15] save period
      XORShift<uint16_t> xs(123, 11, b, c);
      f << testPeriod<uint16_t>(xs) << " ";
    }
    f << std::endl;
  }
  f.close();

  return EXIT_SUCCESS;
}
