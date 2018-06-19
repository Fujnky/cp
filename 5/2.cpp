//CP Blatt 5, Aufgabe 1. Dag-Björn Hering und Lars Funke: Heisenberg
#include <bitset>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <chrono>
#include <fstream>

using std::bitset;
using Eigen::MatrixXd;

int modulo(int arg1, int arg2)
{
  return (arg1 % arg2 + arg2) % arg2;
}

template<size_t N>
MatrixXd generate_submatrix(std::vector<std::pair<uint64_t, uint64_t>>& states, double J)
{
  //generate diagonal matrix
  MatrixXd submatrix = MatrixXd::Identity(states.size(), states.size())* (-J/4) * N;

  //iterate states
  for(int j = 0; j < states.size(); j++)
  {
    //convert uint64_t to std::bitset
    bitset<N> state(states[j].first);
    for(int i = 0; i < N; i++)
    {
      //permute
      auto new_state(state);
      new_state[modulo(i+1, N)] = state[i];
      new_state[i] = state[modulo(i+1, N)];

      //construct subspace
      for(int n = 0; n < states.size(); n++)
      {
        if(new_state.to_ullong() == states[n].first)
          submatrix(n, j) += J/2;
      }
    }
  }

  return submatrix;
}

template<size_t N>
void heisenberg(std::ofstream& runtimes)
{
  static_assert(N <= 14, "Only systems with N <= 14 work.");
  auto start = std::chrono::high_resolution_clock::now();

  MatrixXd example_submatrix;
  MatrixXd hamiltonian (1 << N, 1 << N);
  std::vector<std::pair<uint64_t, uint64_t>> m(1 << N); // 1<<N == 2^N

  for(uint64_t i = 0; i < (1 << N); i++)
  {
    m[i].first = i; //all permutations
    m[i].second = bitset<N>(i).count(); //digit sum
  }

  //sort by digit sum
  std::sort(m.begin(), m.end(),
        [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) -> bool
        { return a.second < b.second; });

  //find subspaces (with same digit sum) and construct hamiltonian
  uint64_t previous = 0;
  size_t previous_index = 0;
  for(uint64_t i = 0; i <= (1 << N); i++)
  {
    if(i == (1 << N) || m[i].second != previous)
    {
      //get permutations with same digit sum
      std::vector<std::pair<uint64_t, uint64_t>> subvec(m.begin() + previous_index, m.begin() + i);
      auto subm = generate_submatrix<N>(subvec, 1);
      if(m[i].second == N/2 && N == 10)
        example_submatrix = subm;

      //insert subspace
      hamiltonian.block(previous_index, previous_index, i-previous_index, i-previous_index) = subm;

      if(i != 1 << N) {
        previous_index = i;
        previous = m[i].second;
      }
    }
  }

  runtimes << N << " " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << std::endl;

  if(N == 10)
  {
    std::ofstream f("build/submatrix");
    f << example_submatrix << std::endl;
    auto ev = hamiltonian.eigenvalues();

    std::ofstream g("build/ev");
    for(size_t i = 0; i < ev.rows(); i++)
    {
      g << ev(i).real() << std::endl;
    }
  }
}

int main()
{
  std::ofstream f("build/runtimes");
  f << "#N t/µs" << std::endl;
  //can't loop because of compile time constant template argument
  heisenberg<2>(f);
  heisenberg<4>(f);
  heisenberg<6>(f);
  heisenberg<8>(f);
  heisenberg<10>(f);
  heisenberg<12>(f);
  heisenberg<14>(f);
}
