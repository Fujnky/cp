//CP Blatt 3, Aufgabe 1. Dag-Björn Hering und Lars Funke: Single spin Metropolis
#include <fstream>
#include <random>

//simulate single spin in magnetic field at β = 1, return magnetization m = <s>
double ising0d(std::mt19937& rng, unsigned int samples, double magnetic_field)
{
  std::discrete_distribution<> dist_initial_spin ({-1, 1});
  std::uniform_real_distribution<> dist_uniform (0, 1);
  int accumulator = 0;
  int spin = dist_initial_spin(rng); //initialize spin randomly

  for(size_t i = 0; i < samples; i++){
    accumulator += spin; //add up every round's spin for averaging
    double delta_E = 2 * spin * magnetic_field;

    //flip spin if ΔE<0 or random condition is met
    spin *= delta_E < 0 || dist_uniform(rng) < exp(-delta_E)? -1 : 1;
  }

  return (double) accumulator / samples; //calculate mean
}

int main()
{
  std::mt19937 rng;
  std::ofstream f("build/ising0d.txt");

  //iterate through different magnetic field strengths
  for(double H = -5; H <= 5; H += 10.0 / 1e3) //1e3 steps are enough
    f << H << " " << ising0d(rng, 1e4, H) << std::endl; //1e4 samples are enough
  return EXIT_SUCCESS;
}
