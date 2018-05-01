//Simulate self-avoiding random walk in 2D and save mean squared start-to-end length for step counts 10...60
#include <iostream>
#include <random>
#include <cstring>
#include <fstream>

int saw2d(int steps, std::mt19937& rng)
{
	//create integer distribution to get int from [0,3]
	std::uniform_int_distribution<> dist(0, 3);
	//create boolean grid big enough for step count
	bool lattice[2 * steps + 1][2 * steps + 1];
	//initialize with false (= not visited)
	std::memset(lattice, false, sizeof(lattice));
	//set origin at center of lattice
	int x = steps;
	int y = steps;
	lattice[x][y] = true;

	//iterate SAW steps
	for(int i = 0; i < steps; i++)
	{
		int xnew, ynew;

		//terminate loop if no further steps are legal
		if(lattice[x+1][y] && lattice[x-1][y] && lattice[x][y+1] && lattice[x][y-1])
			break;

		//bad implementation, worst case time complexity is infinite
		//maybe remove illegal steps from random pool
		do
		{
			xnew = x;
			ynew = y;
			//generate random number and choose step direction
			switch(dist(rng))
			{
				case 0:	xnew = x + 1; break;
				case 1:	xnew = x - 1; break;
				case 2:	ynew = y + 1; break;
				case 3:	ynew = y - 1; break;
			}
		}
		//repeat when illegal
		while(lattice[xnew][ynew]);

		//new legal step found -> go and save
		x = xnew;
		y = ynew;
		lattice[x][y] = true;
	}
	//compute squared distance from origin
	return pow(x - steps, 2) + pow(y - steps, 2);
}

int main()
{
	//initialize PRNG
	std::mt19937 rng;
	rng.seed(std::random_device()());
	size_t size = 1e3; //sample size
	long result; //temporary output array
	std::ofstream f;             //file to save results to
	f.open ("build/output.txt"); //file to save results to
	//iterate step counts
	for(int j = 5; j <= 1445; j += 5)
	{
		result = 0;
		std::cout << j << "... " << std::flush;
		//iterate samples and save to temp array
		for(size_t i = 0; i < size; i++)
			result += saw2d(j, rng);
		//compute average R_N^2 and save to file
		f << j << " " << result / (double) size << std::endl;
	}
	std::cout << std::endl;
	return EXIT_SUCCESS;
}
