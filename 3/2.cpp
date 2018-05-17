//CP Blatt 3, Aufgabe 1. Dag-Björn Hering und Lars Funke: Single spin Metropolis
#include <fstream>
#include <random>
#include <algorithm>
#include <iostream>
#include <armadillo>



using std::cout;

int geildulo(int arg1, int arg2)
{
  return (arg1 % arg2 + arg2) % arg2;
}

double proximity_energy(arma::Mat<double>& spins, int i, int j)
{
  constexpr double J = 1;
  if(geildulo(i-1, spins.n_cols) < 0 )
    cout << "go kill me" << geildulo(i+1, spins.n_cols);


  return J * spins.at(i, j) * (spins.at(geildulo(i+1, spins.n_cols), j) + spins.at(geildulo(i-1, spins.n_cols), j) + spins.at(i, geildulo(j+1, spins.n_rows)) + spins.at(i, geildulo(j-1, spins.n_rows)));
}

arma::Mat<double> ising2d(std::mt19937& rng, unsigned int samples, size_t size, double kBT, int init, arma::Mat<double>* history, arma::subview_row<double>* last)
{
  std::uniform_int_distribution<> dist_spin_flip (0, size-1);
  std::uniform_real_distribution<> dist_uniform (0, 1);
  int accumulator = 0;
  arma::Mat<double> spins(size, size);

  for(auto& spin : spins)
    spin = init == 0? (dist_uniform(rng) < 0.5 ? -1 : 1) : init;

  //cout << spins.at(-1, -1) << std::endl;
  double E = 0;
  double m = arma::mean(arma::mean(spins));
  for(size_t i = 0; i < size; i++)
  for(size_t j = 0; j < size; j++)
    E += proximity_energy(spins, i, j);

  for(size_t t = 0; t < samples; t++)
  {
    if(history != nullptr)
    {
      history->at(t, 0) = (E / pow(size, 2));
      history->at(t, 1) = arma::mean(arma::mean(spins));
      history->at(t, 2) = fabs(history->at(t, 1));
    }

    for(size_t k = 0; k < pow(size, 2); k++){
      int i = dist_spin_flip(rng);
      int j = dist_spin_flip(rng);
      double delta_E = 2 * proximity_energy(spins, i, j);
      //spin flipped if ΔE<0 or random condition is met
      if(delta_E < 0 || dist_uniform(rng) < exp(-delta_E / kBT)) {
        spins.at(i, j) *= -1;
        E += delta_E;
        m += 2 * spins.at(i, j) / pow(size, 2);
      }
    }
  }

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
  uint schweeps = 1e3;

  arma::Mat<double> history(schweeps, 3);
//  Anfangsbedingung alle Spins Zufällig gezogen
  ising2d(rng, schweeps, 100, 1, 0, &history, nullptr).save("build/momentaufnahme_1_0.mat", arma::arma_ascii);
  history.save("build/history_1_0.mat", arma::arma_ascii);
  // ising2d(rng, schweeps, 100, 1.66, 0, &history, nullptr).save("build/momentaufnahme_1,66_0.mat", arma::arma_ascii);
  // history.save("build/history_1,66_0.mat", arma::arma_ascii);
  // ising2d(rng, schweeps, 100, 2.33, 0, &history, nullptr).save("build/momentaufnahme_2,33_0.mat", arma::arma_ascii);
  // history.save("build/history_2,33_0.mat", arma::arma_ascii);
  // ising2d(rng, schweeps, 100, 3, 0, &history, nullptr).save("build/momentaufnahme.mat_3_0", arma::arma_ascii);
  // history.save("build/history_3_0.mat", arma::arma_ascii);

// Anfangsbedingung alle Spins in eine Richtung
history.save("build/history_1_0.mat", arma::arma_ascii);
// ising2d(rng, schweeps, 100, 1.66, 1, &history, nullptr).save("build/momentaufnahme_1,66_1.mat", arma::arma_ascii);
// history.save("build/history_1,66_1.mat", arma::arma_ascii);
// ising2d(rng, schweeps, 100, 2.33, 1, &history, nullptr).save("build/momentaufnahme_2,33_1.mat", arma::arma_ascii);
// history.save("build/history_2,33_1.mat", arma::arma_ascii);
// ising2d(rng, schweeps, 100, 3, 1, &history, nullptr).save("build/momentaufnahme.mat_3_1", arma::arma_ascii);
// history.save("build/history_3_1.mat", arma::arma_ascii);



  uint T_steps = 100;
  arma::Mat<double> results(T_steps, 4);
  for(size_t i = 0; i < T_steps; i++)
  {
    double T = 1 + 2 / (double) T_steps * i;
    auto row = results.row(i);
    row[0] = T;
    ising2d(rng, 300, 100, T, 0, nullptr, &row);
  }

  results.save("build/temperature.mat", arma::arma_ascii);

  cout << results;

  return EXIT_SUCCESS;
}
