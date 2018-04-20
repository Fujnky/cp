#include <iostream>
#include <random>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <vector>
#include <math.h>
#include <functional>

void unit_sphere()
{
  //initialize PRNG
  std::mt19937 rng;
	rng.seed(std::random_device()());
	std::uniform_real_distribution<> dist(-1, 1);

  int points = 1e6;
  int accumulator = 0;

  //iterate over samples
  for(int i = 0; i < points; i++)
  {
    //progress
    if(i % 10000 == 0)
      std::cout << std::setprecision(2) << "\r" << (double) (i+1) * 100 / points << "%   ";
    //add if vector in unit sphere
    accumulator += (pow(dist(rng), 2) + pow(dist(rng), 2) + pow(dist(rng), 2)) < 1;
  }

  //compute volume and pi from samples
  double vol = (double) 8 * accumulator / points;
  double t_vol =  M_PI * 4 / 3;

  std::cout << std::setprecision(12) << "\runit sphere volume – MC: " << vol << " theory: " << t_vol << std::endl;
  std::cout << "pi – MC: " << vol * 3 / 4 << " theory: " << M_PI << std::endl;
  std::cout << "Relative error: " << 100 * abs(1 - (vol * 3 / 4 ) / M_PI) << "%" << std::endl;
}

double mc_integration(double a, double b, int samples, std::function<double(double)> func, std::mt19937& rng)
{
  std::uniform_real_distribution<> dist(a, b);
  double h = (b - a) / samples;
  double accumulator = 0;
  for(size_t i = 0; i < samples; i++)
  {
    accumulator += func(dist(rng));
  }
  return h * accumulator;
}


int main()
{
  std::mt19937 rng;
	rng.seed(std::random_device()());
  unit_sphere();


  std::function<double(double)> function1 = [](double x){return -1 / (sqrt(M_PI) * pow(x, 2)) * exp(-1 / pow(x, 2));};
  std::function<double(double)> function2 = [](double x){return -1 / (sqrt(M_PI) * (pow(x, 2) + 1)) * exp(-pow(atan(x), 2));};

  std::cout << mc_integration(0, -1, 1e7, function1, rng) << std::endl;
  std::cout << mc_integration(0, 1/1.1631, 1e7, function1, rng) << std::endl;
  std::cout << mc_integration(-M_PI / 2, M_PI / 2, 1e7, function2, rng) << std::endl;
	return EXIT_SUCCESS;
}
