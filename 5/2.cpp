#include <iostream>
#include <bitset>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <chrono>

using std::bitset;
using std::cout;
using Eigen::MatrixXd;


int modulo(int arg1, int arg2)
{
  return (arg1 % arg2 + arg2) % arg2;
}


template<size_t N>
MatrixXd heisenberg(std::vector<std::pair<uint64_t, uint64_t>>& states, double J)
{
  MatrixXd Untermatrix = MatrixXd::Identity(states.size(), states.size())* (-J/4) * N;
  for(int j = 0; j < states.size(); j++)
  {
    bitset<N> state(states[j].first);
    //cout << state << std::endl;
    for(int i = 0; i < N; i++)
    {
      // coole bits operationen
      auto new_state(state); // kann man auch besser machen
      new_state[modulo(i+1, N)] = state[i];
      new_state[i] = state[modulo(i+1, N)];
      for(int n = 0; n < states.size(); n++)
      {
        if(new_state.to_ullong() == states[n].first)
        Untermatrix(n, j)+= J/2;
      }
    }
  }
  /*if(states[0].second == N/2)
    cout << Untermatrix << std::endl;*/
  return Untermatrix;
}





template<size_t N>
void shizzle()
{
  auto start = std::chrono::high_resolution_clock::now();
  static_assert(N % 2 == 0, "Odd N does not work");
  MatrixXd hamiltonian (1 << N, 1 << N);
  std::vector<std::pair<uint64_t, uint64_t>> m(1 << N);
  //Matrix<uint64_t, 1 << N, 2> m;
  for(uint64_t i = 0; i < (1 << N); i++)
  {
    m[i].first = i;
    m[i].second = bitset<N>(i).count();
  }

  std::sort(m.begin(), m.end(),
        [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) -> bool
        { return a.second < b.second; });

  uint64_t previous = 0;
  size_t previous_index = 0;
  for(uint64_t i = 0; i <= (1 << N); i++)
  {
    //cout << m[i].second << " " << previous << std::endl;
    if(i == (1 << N) || m[i].second != previous)
    {
      //cout << "from" << previous_index << " to " << i << std::endl;
      std::vector<std::pair<uint64_t, uint64_t>> subvec(m.begin() + previous_index, m.begin() + i);


      hamiltonian.block(previous_index, previous_index, i-previous_index, i-previous_index) = heisenberg<N>(subvec, 1);

      if(i != (1<<N)) {
        previous_index = i;
        previous = m[i].second;
      }
    }
  }

  //cout << hamiltonian << std::endl;

  //cout << hamiltonian.eigenvalues() << std::endl;

  //Eigen::Matrix2d mmmmm;
  //mmmmm << 0, 1, -1, 0;
  //cout << mmmmm.eigenvalues() << std::endl;

  //for(uint64_t i = 0; i < (1 << N); i++)
  //  cout << i << " " << m[i].first << " " << m[i].second << std::endl;

  cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << "µs" << std::endl;
}



int main()
{
  //heisenberg<4>(1);
  //heisenberg(4);


  //shizzle<2>();
  //shizzle<4>();
  //shizzle<6>();
  //shizzle<8>();
  //shizzle<10>();
  //shizzle<12>();
  shizzle<14>();
  shizzle<16>();
  //Eigen::MatrixXd(65535ull, 65535ull> lol;
  //shizzle<18>();

}
