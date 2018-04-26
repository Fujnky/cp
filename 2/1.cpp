#include <iostream>
#include <fstream>

//#include <strin
#define RNGTYPE uint32_t


template <typename T>
class LCG
{
public:
  LCG(T r0, T a, T c, T m) : r(r0), a(a), c(c), m(m) {};

  T getRnd()
  {
    r = (a * r + c) % m;
    return r;
  };

  double getRndDouble()
  {
    return getRnd() / (double) m;
  };

private:
  T r, a, c, m;
};

template <typename T>
void testLCG(LCG<T>& lcg, size_t samples, std::string filename)
{
  std::ofstream f;  //file to save results to
	f.open(filename); //file to save results to
  for(size_t i = 0; i < samples; i++)
    f << lcg.getRndDouble() << std::endl;
  f.close();
}

int main()
{
  LCG<RNGTYPE> lcg1(1234,      20,    120, 6075);
  LCG<RNGTYPE> lcg2(1234,      137,   187, 256);
  LCG<RNGTYPE> lcg3(123456789, 65539, 0,   2147483648);
  LCG<RNGTYPE> lcg4(1234,      16807, 0,   2147483647);
  testLCG<RNGTYPE>(lcg1, 1e5, "build/1.txt");
  testLCG<RNGTYPE>(lcg2, 1e5, "build/2.txt");
  testLCG<RNGTYPE>(lcg3, 1e5, "build/3.txt");
  testLCG<RNGTYPE>(lcg4, 1e5, "build/4.txt");
  return EXIT_SUCCESS;
}
