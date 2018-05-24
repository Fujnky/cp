//CP Blatt 3, Aufgabe 1. Dag-Björn Hering und Lars Funke: Wolff
#include <fstream>
#include <random>
#include <algorithm>
#include <armadillo>
#include <sstream>
#include <iomanip>
#include <queue>
#include <math.h>

//differs from c++ variant for negative numbers. needed for periodicity
int modulo(int arg1, int arg2)
{
  return (arg1 % arg2 + arg2) % arg2;
}

struct Position {
  int x;
  int y;
  size_t size;

  inline Position operator+(Position a) {
      return {modulo(a.x + x, size), modulo(a.y + y, size), size};
  }
};


//simulate 2D Ising model
std::vector<double> wolff(std::mt19937& rng, unsigned int steps, size_t size, double beta, int init, arma::Cube<uint8_t>* history, double J, uint q, uint save_every_nth_frame, bool compute_m_and_E, uint ignore_equilibration_steps)
{
  double p_c = 1 - exp(-beta*J);
  std::uniform_int_distribution<> dist_position (0, size-1);
  std::uniform_real_distribution<> dist_uniform (0, 1);
  std::uniform_int_distribution<> dist_value (0, q-1);
  std::uniform_int_distribution<> dist_value_flip (1, q-1);
  const Position neighbors[] = {{0, 1, size}, {0, -1, size}, {1, 0, size}, {-1, 0, size}};
  //spin field, double for easier calculations
  arma::Mat<uint8_t> spins(size, size);
  //initialize randomly if parameter is zero
  for(auto& spin : spins)
    spin = init < 0 ? dist_value(rng) : init;
  double m_ges;
  double E;
  double theta = 2 * M_PI / q;
  std::queue<Position> unvisited;

  for(size_t t = 0; t < steps; t++)
  {
    if(history != nullptr && t % save_every_nth_frame == 0)
      history->slice(t/save_every_nth_frame) = spins;

    Position cur = {dist_position(rng), dist_position(rng), size};
    unvisited.push(cur);
    int first_spin = spins.at(cur.x, cur.y);
    int motherflip = (first_spin + dist_value_flip(rng)) % q;
    // std::cerr << "s(" << cur.x << ", " << cur.y << ") = " << first_spin << " => " << motherflip << std::endl;
    while(!unvisited.empty())
    {
      cur = unvisited.front();
      unvisited.pop();
      if(spins.at(cur.x, cur.y) != first_spin) continue;
      //std::cout << cur.x << " " << cur.y << std::endl;
      for(const auto& n : neighbors){
        auto neww = cur + n;
        // std::cerr << cur.x << ", " << cur.y << " N.N" << neww.x << ", " << neww.y << std::endl;
        if(first_spin == spins.at(neww.x, neww.y) && dist_uniform(rng) < p_c)
        {
          unvisited.push(cur + n);
        }
      }
      spins.at(cur.x, cur.y) = motherflip;
    }

    if(!compute_m_and_E || t < ignore_equilibration_steps) continue;
    //herzlich willkommen zur berechnung vong magnetizirung und enrschie vong desch gitters her mein damen und herrn treten sie ein tadaaa
    double m_x = 0;
    double m_y = 0;
    m_ges = 0;
    E = 0;
    for(int x = 0; x < spins.n_cols; x++)
      for(int y = 0; y < spins.n_rows; y++)
      {
        Position cur = {x, y, size};
        for(const auto& n : neighbors)
        {
          auto neww = cur + n;
          E += spins.at(x,y) == spins.at(neww.x, neww.y)? -J : 0;
          m_x += cos((spins.at(x,y) + 1) * theta);
          m_y += sin((spins.at(x,y) + 1) * theta);
        }
      }
     m_ges += sqrt(m_x * m_x + m_y * m_y);
  }

  return std::vector<double>({m_ges/(steps-ignore_equilibration_steps), E/(steps-ignore_equilibration_steps)});
}

int main()
{
  std::mt19937 rng;
  rng.seed(std::random_device()());
  uint steps = 6000; //equlibrium is reached after ~400 sweeps
  uint size = 100;
  uint save_every_nth_frame = 10;
  int J = 1;
  uint q = 3;

  double beta_c = 1 / J * log (1 + sqrt(q));

  double beta = 1.165;
  arma::Cube<uint8_t> history(size, size, steps/save_every_nth_frame);

  //iterate different parameters (initialization mode & temperature)
  //for(int init = 0; init <= 1; init++)
  //for(const auto& temp : kBT)
  //{
    //std::stringstream file;
    //file << "_" << std::fixed << std::setprecision(2) << temp << "_" << init;


    //wolff(rng, steps, size, beta, -1, &history, J, q, save_every_nth_frame, false);

    //history.save(std::cout, arma::arma_ascii);
  //}

  //iterate temperature steps to plot m(T) etc.
  uint T_steps = 100;
  uint T_model_steps = 2500;
  arma::Mat<double> results(T_steps, 3);
  #pragma omp parallel for
  for(size_t i = 0; i < T_steps; i++)
  {
    //std::cerr << '\r' << i+1 << "/" << T_steps;
    double beta_ = 0.95 + beta_c * 0.15 / (double) T_steps * i;
    auto res = wolff(rng, T_model_steps, size, beta_, -1, nullptr, J, q, 0, true, 1500);
    results.at(i, 0) = beta_ ;
    results.at(i, 1) = res[0];
    results.at(i, 2) = res[1];
  }
  std::cerr << std::endl;
  results.save("build/beta.mat", arma::arma_ascii);
  return EXIT_SUCCESS;
}
