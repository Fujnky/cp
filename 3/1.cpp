//CP Blatt 2, Aufgabe 1. Dag-Björn Hering und Lars Funke
//Single spin Metropolis

#include <iostream>
#include <fstream>
#include <random>

//simulate single spin in magnetic field H, return magnetization = <s>
double ising0d(std::mt19937& rng, unsigned int samples, double H)
{
  int accumulator = 0;
  std::discrete_distribution<> dist_initial_spin ({-1, 1});
  std::uniform_real_distribution<> dist (0, 1);

  //initialize spin randomly
  int s = dist_initial_spin(rng);

  for(size_t i = 0; i < samples; i++){
    accumulator += s; //add up current spin for averaging
    double delta_E = s * 2 * H;

    //flip spin if ΔE≤0 or random condition is met
    s *= delta_E <= 0 || dist(rng) < exp(-delta_E) ? -1 : 1;
  }

  return (double) accumulator / samples; //calculate mean
}

int main()
{
  std::mt19937 rng;
	rng.seed(std::random_device()());
  std::ofstream f;
  f.open("build/ising0d.txt");

  //iterate through different magnetic field strenghts
  for(double H = -5; H <= 5; H += (double) 10 / 1e3) //1e3 steps are enough
    f << H << " " << ising0d(rng, 1e4, H) << std::endl; //1e4 samples are enough
}
