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

  //periodic boundary conditions implemented within addition
  inline Position operator+(Position a) {
      return {modulo(a.x + x, size), modulo(a.y + y, size), size};
  }
};


//simulate 2D Ising model
std::vector<double> wolff(std::mt19937& rng,
  unsigned int steps, size_t size,
  double beta, int init,
  arma::Cube<uint8_t>* history,
  double J,
  uint q,
  uint save_every_nth_step,
  bool compute_m_and_E,
  uint ignore_equilibration_steps)
{
  //cluster inclusion probability
  double p_c = 1 - exp(-beta*J);

  std::uniform_int_distribution<> dist_position (0, size-1);
  std::uniform_real_distribution<> dist_uniform (0, 1);
  std::uniform_int_distribution<> dist_value (0, q-1);
  std::uniform_int_distribution<> dist_value_flip (1, q-1);

  const Position neighbors[] = {{0, 1, size}, {0, -1, size},
    {1, 0, size}, {-1, 0, size}};

  //spin lattice
  arma::Mat<uint8_t> spins(size, size);

  //initialize randomly if parameter is -1
  for(auto& spin : spins)
    spin = init < 0 ? dist_value(rng) : init;

  double m_ges;
  double E;
  double theta = 2 * M_PI / q;
  std::queue<Position> unvisited;

  for(size_t t = 0; t < steps; t++)
  {
    //save history if wanted and only every save_every_nth_step steps
    if(history != nullptr && t % save_every_nth_step == 0)
      history->slice(t/save_every_nth_step) = spins;

    //roll the dice for start position
    Position cur = {dist_position(rng), dist_position(rng), size};
    unvisited.push(cur);
    int first_spin = spins.at(cur.x, cur.y); //spin at start
    int motherflip = (first_spin + dist_value_flip(rng)) % q; //spin to flip cluster to
    //process queue
    while(!unvisited.empty())
    {
      cur = unvisited.front();
      unvisited.pop();
      //ignore duplicate spins (already flipped)
      if(spins.at(cur.x, cur.y) != first_spin) continue;

      for(const auto& n : neighbors){ //iterate next neighbors
        auto neww = cur + n;//periodic boundary conditions implemented in struct
        if(first_spin == spins.at(neww.x, neww.y) && dist_uniform(rng) < p_c)
        {
          unvisited.push(cur + n);
        }
      }
      spins.at(cur.x, cur.y) = motherflip; //I'd flip that
    }

    //stop if m and E are not wanted or still warming up
    if(!compute_m_and_E || t < ignore_equilibration_steps) continue;
    double m_x = 0;
    double m_y = 0;
    m_ges = 0;
    E = 0;
    //iterate lattice and compute energy and magnetization
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
  return std::vector<double>({m_ges/(steps-ignore_equilibration_steps),
    E/(steps-ignore_equilibration_steps)});
}

int main(int argc, char* argv[])
{
  std::mt19937 rng;
  rng.seed(std::random_device()());
  uint size = 100;
  int J = 1;

  if(argc == 3 && strcmp(argv[1], "--history") == 0)
  {
    uint q = atoi(argv[2]);
    double beta_c = 1 / J * log (1 + sqrt(q));

    uint steps = pow(q, 2) * 500; //guessed
    uint nth_step = 10;
    arma::Cube<uint8_t> history(size, size, steps/nth_step);
    wolff(rng, steps, size, beta_c, -1, &history, J, q, nth_step, false, 0);
    std::cout << nth_step << " " << q << std::endl;
    history.save(std::cout, arma::arma_ascii);
  }

  if(argc == 3 && strcmp(argv[1], "--beta") == 0)
  {
    //iterate temperature steps to plot m(T) etc.
    uint q = atoi(argv[2]);
    double beta_c = 1 / J * log (1 + sqrt(q));
    uint beta_steps = 40;
    uint steps = 2500;
    arma::Mat<double> results(beta_steps, 3);
    #pragma omp parallel for
    for(size_t i = 0; i < beta_steps; i++)
    {
      double beta = 0 + beta_c * 2 / (double) beta_steps * i;
      auto res = wolff(rng, steps, size, beta, -1, nullptr, J, q, 0, true, 2250);
      results.at(i, 0) = beta ;
      results.at(i, 1) = res[0];
      results.at(i, 2) = res[1];
    }
    std::cout << q << std::endl;
    results.save(std::cout, arma::arma_ascii);
  }

  return EXIT_SUCCESS;
}
