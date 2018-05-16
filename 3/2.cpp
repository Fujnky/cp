//CP Blatt 3, Aufgabe 1. Dag-Björn Hering und Lars Funke: Single spin Metropolis
#include <fstream>
#include <random>
#include <algorithm>
#include <iostream>
#include <armadillo>



using std::cout;

double energy(arma::Mat<int>& spins)
{
  double energy = 0;
  double J = 1;
  for(int i = 0; i < spins.n_rows; i++)
  for(int j = 0; j < spins.n_cols; j++)
  {
    /*energy +=  -J* spins.at((i+1) % spins.n_cols, (j  ) % spins.n_rows) * spins.at((i-1) % spins.n_cols, (j  ) % spins.n_rows);
    energy +=  -J* spins.at((i+1) % spins.n_cols, (j  ) % spins.n_rows) * spins.at((i  ) % spins.n_cols, (j-1) % spins.n_rows);
    energy +=  -J* spins.at((i+1) % spins.n_cols, (j  ) % spins.n_rows) * spins.at((i  ) % spins.n_cols, (j+1) % spins.n_rows);
    energy +=  -J* spins.at((i-1) % spins.n_cols, (j  ) % spins.n_rows) * spins.at((i  ) % spins.n_cols, (j-1) % spins.n_rows);
    energy +=  -J* spins.at((i-1) % spins.n_cols, (j  ) % spins.n_rows) * spins.at((i  ) % spins.n_cols, (j+1) % spins.n_rows);
    energy +=  -J* spins.at((i  ) % spins.n_cols, (j-1) % spins.n_rows) * spins.at((i  ) % spins.n_cols, (j+1) % spins.n_rows);*/
    energy += J * spins.at(i, j) * (spins.at(i+1, j) + spins.at(i-1, j) + spins.at(i, j+1) + spins.at(i, j-1));
  }
  return energy;
}

arma::Mat<int> ising2d(std::mt19937& rng, unsigned int samples, size_t size, double kBT, int init)
{
  std::uniform_int_distribution<> dist_spin_flip (0, 100);
  std::uniform_real_distribution<> dist_uniform (0, 1);
  int accumulator = 0;
  arma::Mat<int> spins(size, size);

  for(auto& spin : spins)
    spin = init == 0? (dist_uniform(rng) < 0.5 ? -1 : 1) : init;

  cout << spins.at(-1, -1) << std::endl;
  double E = energy(spins);

  for(size_t k = 0; k < samples; k++){
    int i = dist_spin_flip(rng);
    int j = dist_spin_flip(rng);

    spins.at(i, j) *= -1;
    double delta_E = energy(spins) - E;
    //std::cout << "delta_E="<<delta_E << std::endl;
    //keep spin flipped if ΔE<0 or random condition is met
    spins.at(i, j) *= delta_E < 0 || dist_uniform(rng) < exp(-delta_E / kBT)? 1 : -1;
  }

  return spins;
}


int main()
{
  std::mt19937 rng;
  rng.seed(std::random_device()());
  std::ofstream f("build/ising2d.txt");

  ising2d(rng, 1e4, 100, 0.1, 0).save("build/momentaufnahme.mat", arma::arma_ascii);
  return EXIT_SUCCESS;
}
