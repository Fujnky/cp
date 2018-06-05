#include <iostream>
#include <bitset>

using std::bitset;

template<size_t N>
void heisenberg(bitset<N>& spins, int S_z)
{
  std::cout << spins.count() << std::endl;
}

int main()
{
  bitset<10> foo;
  heisenberg(foo, 1);
  //heisenberg(4);
}
