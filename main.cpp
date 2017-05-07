#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>

#include "cimg.h"
#include "sibson.h"

using namespace cimg_library;

void usage()
{
	std::cout << "Usage: sibson n input.bmp output.bmp" << std::endl;
	std::cout << "n: number of interpolation points (integer, >= 1)" << std::endl;
	std::cout << "input: input file in bmp format (24 bits)" << std::endl;
	std::cout << "output: output file in bmp format (24 bits)" << std::endl;
}

int main(int argc, char** argv)
{
	if (argc < 4)
	{ 
		usage();
		return 1;
	}

	// scan parameters
	int nb_points = 1;
	sscanf(argv[1], "%i", &nb_points);

	std::cout << yellow << "SIBSON INTERPOLATION" << white << std::endl;
	std::cout << "number of points: " << nb_points << std::endl;

	// Sibson interpolation
	Sibson sibson;
	sibson.run(nb_points, argv[2], argv[3]);

  return 0;
}
