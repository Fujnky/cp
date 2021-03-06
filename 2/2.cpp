//CP Blatt 2, Aufgabe 1. Dag-Björn Hering und Lars Funke
//Random number transformations

#include <iostream>
#include <fstream>
#include <random>
#include <functional>


//transform uniform distribution to dist_to and save to file
void sample_distribution(std::ofstream& file,
                         int samples,
                         std::uniform_real_distribution<> dist_from,
                         std::function<double(double)>& dist_to,
                         std::mt19937& rng) {
  for(size_t i = 0; i < samples; i++)
  {
    file << dist_to(dist_from(rng)) << " ";
  }
}

//implement box muller algorithm
void box_muller(std::ofstream& file, int samples, std::mt19937& rng, double sigma, double mu)
{
  std::uniform_real_distribution<> dist(0, 1);
  for(size_t i = 0; i < samples; i++)
  {
    double u1 = dist(rng);
    double u2 = dist(rng);
    file << sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) + mu << " ";
    file << sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) + mu << " ";
  }
}



int main()
{
  std::ofstream f_1;               //file to save results to
  std::ofstream f_2;               //file to save results to
	f_1.open ("build/output_1.txt"); //file to save results to
	f_2.open ("build/output_2.txt"); //file to save results to
  std::mt19937 rng;
  int samples = 1e5;
	rng.seed(std::random_device()());

  //define transformation functions
  std::function<double(double)> function1 = [](double x){return asin(x) ;};
  std::function<double(double)> function2 = [](double x){return asin((0.5)*(x+1)) ;};

  //apply transformation
  sample_distribution(f_1, samples, std::uniform_real_distribution<>(0, 1), function1, rng);
  sample_distribution(f_2, samples, std::uniform_real_distribution<>(-1, 1), function2, rng);

  std::ofstream f_3;
  f_3.open("build/output_3.txt");
  box_muller(f_3, samples, rng, 2, 3);

  return EXIT_SUCCESS;
}
