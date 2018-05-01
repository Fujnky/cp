#include <iostream>
#include <fstream>
#include <random>
#include <functional>

void sample_distribution(std::ofstream& file,
                         int samples,
                         std::uniform_real_distribution<> dist_from,
                         std::function<double(double)> dist_to,
                         std::mt19937& rng) {
  for(size_t i = 0; i < samples; i++)
  {
    file << dist_to(dist_from(rng)) << " ";
  }
}



int main()
{
  std::ofstream f_1;             //file to save results to
  std::ofstream f_2;             //file to save results to
	f_1.open ("build/output_1.txt"); //file to save results to
	f_2.open ("build/output_2.txt"); //file to save results to
  std::mt19937 rng;
  int samples = 1e4;
	rng.seed(std::random_device()());
  std::function<double(double)> function1 = [](double x){return asin(x) ;};
  std::function<double(double)> function2 = [](double x){return asin((0.5)*(x+1)) ;};

  sample_distribution(f_1, samples, std::uniform_real_distribution<>(0, 1), function1, rng);
  sample_distribution(f_2, samples, std::uniform_real_distribution<>(-1, 1), function2, rng);

  return EXIT_SUCCESS;
}
