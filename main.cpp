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
	std::cout << "n: number of interpolation points (integer, >= 1)(useless in mode 2, 3)" << std::endl;
  std::cout << "maximum error: the tolerance of error in the final image(between 0 and 255)(useless in mode 1,4)" << std::endl;
  std::cout << "mode: 1 - fix number of points" << std::endl; 
  std::cout << "      2 - add points one by one to meet the requirement of error" << std::endl;
  std::cout << "      3 - use a hierarchic method to meet the requirement of error" << std::endl;
  std::cout << "      4 - minimize the error when the number of point is fix." << std::endl;
	std::cout << "input: input file in bmp format (24 bits)" << std::endl;
	std::cout << "output: output file in bmp format (24 bits)" << std::endl;
}

int main(int argc, char** argv)
{
	if (argc < 6)
	{ 
		usage();
		return 1;
	}

	// scan parameters
  int nb_points = 1;
  sscanf(argv[1], "%i", &nb_points);
	float max_error = 50.0;
	sscanf(argv[2], "%f", &max_error);
  int mode = 1;
  sscanf(argv[3], "%i", &mode);
  if(mode < 1 || mode > 4){
    std::cout << "The mode is illegal> It should be 1, 2, 3 or 4." << std::endl;
    return 1;
  }

	std::cout << yellow << "SIBSON INTERPOLATION" << white << std::endl;
  std::cout << "MODE: " << mode << std::endl;

  Sibson sibson;
  if(mode == 1){
    std::cout << "Number of points: " << nb_points << std::endl;
    sibson.run(nb_points, mode, argv[4], argv[5]);
  }
  else if(mode == 4){
    std::cout << "Number of points: " << nb_points << std::endl;
    sibson.run(nb_points, argv[4], argv[5]);
  }
  else{
	  std::cout << "maximum error: " << max_error << std::endl;
    sibson.run(max_error, mode, argv[4], argv[5]);
  }

  return 0;
}
