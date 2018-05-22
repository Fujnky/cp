//CP Blatt 3, Aufgabe 1. Dag-Björn Hering und Lars Funke: 2D Ising fun
#include <fstream>
#include <random>
#include <algorithm>
#include <armadillo>
#include <sstream>
#include <iomanip>

//differs from c++ variant for negative numbers. needed for periodicity
int modulo(int arg1, int arg2)
{
  return (arg1 % arg2 + arg2) % arg2;
}

//energy of next neighbors
double proximity_energy(arma::Mat<double>& spins, int i, int j)
{
  constexpr double J = 1;
  return J * spins.at(i, j) * (spins.at(modulo(i+1, spins.n_cols), j)
  + spins.at(modulo(i-1, spins.n_cols), j)
  + spins.at(i, modulo(j+1, spins.n_rows))
  + spins.at(i, modulo(j-1, spins.n_rows)));
}

//simulate 2D Ising model
arma::Mat<double> ising2d(std::mt19937& rng, unsigned int sweeps, size_t size, double kBT, int init, arma::Mat<double>* history, arma::subview_row<double>* last)
{
  std::uniform_int_distribution<> dist_spin_flip (0, size-1);
  std::uniform_real_distribution<> dist_uniform (0, 1);
  int accumulator = 0;
  //spin field, double for easier calculations
  arma::Mat<double> spins(size, size);

  //initialize randomly if parameter is zero
  for(auto& spin : spins)
    spin = init == 0? (dist_uniform(rng) < 0.5 ? -1 : 1) : init;

  //calculate initial energy
  double E = 0;
  for(size_t i = 0; i < size; i++)
  for(size_t j = 0; j < size; j++)
    E += proximity_energy(spins, i, j);

  for(size_t t = 0; t < sweeps; t++)
  {
    //save history if wanted (energy and magnetization)
    if(history != nullptr)
    {
      history->at(t, 0) = (E / pow(size, 2));
      history->at(t, 1) = arma::mean(arma::mean(spins));
      history->at(t, 2) = fabs(history->at(t, 1));
    }

    //do sweep
    for(size_t k = 0; k < pow(size, 2); k++){
      int i = dist_spin_flip(rng);
      int j = dist_spin_flip(rng);
      double delta_E = 2 * proximity_energy(spins, i, j);

      //flip spin if ΔE<0 or random condition is met
      if(delta_E < 0 || dist_uniform(rng) < exp(-delta_E / kBT)) {
        spins.at(i, j) *= -1;
        E += delta_E;
      }
    }
  }

  //save last state if wanted (energy and magnetization)
  if(last != nullptr)
  {
    (*last)[1] = (E / pow(size, 2));;
    (*last)[2] = arma::mean(arma::mean(spins));
    (*last)[3] = fabs((*last)[2]);
  }

  return spins;
}

int main()
{
  std::mt19937 rng;
  rng.seed(std::random_device()());
  uint schweeps = 1e3; //equlibrium is reached after ~400 sweeps

  arma::Mat<double> history(schweeps, 3);
  std::vector<double> kBT({1, 1.66, 2.33, 3});

  //iterate different parameters (initialization mode & temperature)
  for(int init = 0; init <= 1; init++)
  for(const auto& temp : kBT)
  {
    std::stringstream file;
    file << "_" << std::fixed << std::setprecision(2) << temp << "_" << init;
    ising2d(rng, schweeps, 100, temp, init, &history, nullptr).save("build/momentaufnahme"+file.str()+".mat", arma::arma_ascii);
    history.save("build/history"+file.str()+".mat", arma::arma_ascii);
  }

  //iterate temperature steps to plot m(T) etc.
  uint T_steps = 80;
  uint T_schweeps = 400;
  arma::Mat<double> results(T_steps, 5);
  arma::Mat<double> T_history(T_schweeps, 3);
  for(size_t i = 0; i < T_steps; i++)
  {
    double T = 1 + 2 / (double) T_steps * i;
    auto row = results.row(i);
    ising2d(rng, T_schweeps, 100, T, 0, &T_history, &row);
    row[0] = T;
    row[4] = arma::mean(T_history.tail_rows(100).col(1)*10000) / 100;
  }
  results.save("build/temperature.mat", arma::arma_ascii);
  return EXIT_SUCCESS;
}
