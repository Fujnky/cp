#include <iostream>
#include <bitset>
#include <cmath>
#include <eigen>

using std::bitset;

template<size_t N>
void heisenberg(int S_z)
{
  bitset<N> spins[1 << N];
  //std::cout << spins[0].count() << std::endl;

  //todo: generalize for arbitrary system sizes > 64bit
  for(uint64_t i = 0; i < (1 << N); i++)
  {
    spins[i] = bitset<N>(i);
    std::cout << spins[i] << std::endl;
  }
}



int main()
{
  heisenberg<4>(1);
  //heisenberg(4);
}
