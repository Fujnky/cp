#include <iostream>
#include <fstream>
#include <limits>

class RNG
{
  public:
    virtual double getRndDouble() = 0;
};


class XORShift : public RNG
{
public:
  XORShift(uint16_t x, uint16_t a, uint16_t b, uint16_t c) : x(x), a(a), b(b), c(c) {};

  uint16_t getRnd()
  {
    x ^= (x << a);
    x ^= (x >> b);
    x ^= (x << c);
    return x;
  }

  double getRndDouble()
  {
    return getRnd() / (double) std::numeric_limits<uint16_t>::max();
  }

private:
  uint16_t x, a, b, c;
};

class LCG : public RNG
{
public:
  LCG(uint32_t r0, uint32_t a, uint32_t c, uint32_t m) : r(r0), a(a), c(c), m(m) {};

  uint32_t getRnd()
  {
    r = (a * r + c) % m;
    return r;
  };

  double getRndDouble()
  {
    return getRnd() / (double) m;
  };

private:
  uint32_t r, a, c, m;
};

int testPeriod(XORShift& rng)
{
  uint16_t first = rng.getRnd();
  int i = 1;
  for(; rng.getRnd() != first; i++);
  return i;
}

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
  LCG lcg1(1234,      20,    120, 6075);
  LCG lcg2(1234,      137,   187, 256);
  LCG lcg3(123456789, 65539, 0,   2147483648);
  LCG lcg4(1234,      16807, 0,   2147483647);
  XORShift xs1(123, 11, 1, 7);
  XORShift xs2(123, 11, 4, 7);
  testRNG(lcg1, 1e5, "build/1.txt");
  testRNG(lcg2, 1e5, "build/2.txt");
  testRNG(lcg3, 1e5, "build/3.txt");
  testRNG(lcg4, 1e5, "build/4.txt");
  testRNG(xs1,  1e5, "build/5.txt");
  testRNG(xs2,  1e5, "build/6.txt");

  std::ofstream f;  //file to save results to
	f.open("build/period.txt"); //file to save results to
  for(int b = 1; b <= 15; b++)
  {
    for(int c = 1; c <= 15; c++)
    {
      XORShift xs(123, 11, b, c);
      //f << "b" << b << " c" << c << " " << testPeriod(xs) << " | ";
      f << testPeriod(xs) << " ";
    }
    f << std::endl;
  }
  f.close();

  return EXIT_SUCCESS;
}
