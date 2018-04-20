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

  int points = 1e9;
  int accumulator = 0;

  //iterate over samples
  for(int i = 0; i < points; i++)
  {
    //progress
    if(i % 10000 == 0)
      std::cout << std::setprecision(2) << "\r" << (double) (i+1) * 100 / points << "%   ";
    //add if vector ends inside unit sphere
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
  //generate uniform random numbers between limits
  std::uniform_real_distribution<> dist(a, b);
  double h = (b - a) / samples;
  double accumulator = 0;

  for(size_t i = 0; i < samples; i++)
  {
    //add function value
    accumulator += func(dist(rng));
  }
  return h * accumulator;
}

double simpson_integration(double a, double b, int nodes, std::function<double(double)> func ){
  //compute simpson integration from a to b with nodes nodes of function func.

  double h = (b - a) / nodes;
  double accumulator = 0;
  for (int i = 0; i < nodes; i++){
    accumulator += func(a + h * i) + 4 * func(((a + h * i) + h / 2)) + func(a + h * i);
  }
  return h * accumulator / 6 ;
}


int main()
{
  std::mt19937 rng;
	rng.seed(std::random_device()());

  //a)
  unit_sphere();

  //b)

  //function to integrate
  std::function<double(double)> function2 = [](double x){return 1 / (sqrt(M_PI) * pow(cos(x), 2)) * exp(-pow(tan(x), 2));};

  std::cout << "Monte Carlo integration:" << std::endl;
  std::cout << std::setprecision(8) << mc_integration(-M_PI / 2, atan(-1),     1e9, function2, rng) << std::endl;
  std::cout << std::setprecision(8) << mc_integration(-M_PI / 2, atan(1.1631), 1e9, function2, rng) << std::endl;
  std::cout << std::setprecision(8) << mc_integration(-M_PI / 2, M_PI / 2,     1e9, function2, rng) << std::endl;

  std::cout << "Simpson rule:" << std::endl;
  std::cout << std::setprecision(8) << simpson_integration(-M_PI / 2, atan(-1),     1e6, function2) << std::endl;
  std::cout << std::setprecision(8) << simpson_integration(-M_PI / 2, atan(1.1631), 1e6, function2) << std::endl;
  std::cout << std::setprecision(8) << simpson_integration(-M_PI / 2, M_PI / 2,     1e6, function2) << std::endl;

	return EXIT_SUCCESS;
}
