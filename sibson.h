#include <list>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>

#include "cimg.h"
#include "point.h"
#include "random.h"
#include "console.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
typedef std::vector< std::pair< K::Point_2, K::FT  > >
                                                      Point_coordinate_vector;

using namespace cimg_library;

class Sibson
{
public:
	Sibson() {}
	~Sibson() {}

public:
	bool run(const int nb_points,
		const char *pInput,
		const char *pOutput)
	{
		// read input image
		std::cout << "read input image...";
		CImg<unsigned char> input(pInput);
		
		// verify image size
		const int width = input.width();
		const int height = input.height();
		if (width * height == 0)
		{
			std::cerr << red << "failed" << white << std::endl;
			return false;
		}
		std::cout << "done" << std::endl;

		// generate initial points to interpolate
		// at the four image corners
		std::list<Point> points;
		points.push_back(Point(0.0, 0.0));
		points.push_back(Point(width - 1.0, 0.0));
		points.push_back(Point(0.0, height - 1.0));
		points.push_back(Point(width - 1.0, height - 1.0));

		// TODO: add random points
		// tip: use eg function ::random_int(0, width - 1) to get a random x coordinate value 
		std::cout << "adding " << nb_points << " random points...";
		for (int i = 0; i < nb_points; i++)
		{
      int x = random_int(0, width - 1);
      int y = random_int(0, height - 1);
      
      points.push_back(Point(x, y));
		}
    points.unique();
		std::cout << "done" << std::endl;

		// TODO: incremental addition
		// tip: evaluate error in RGB or grey level space


    Delaunay_triangulation dt;
    // do a delaunay triangulation 
    for (std::list<Point>::iterator it = points.begin(); it != points.end(); it++)
      dt.insert(K::Point_2(it->x(),it->y()));

		// interpolation
		std::cout << "interpolating...";
		CImg<unsigned char> output(width, height, 1, 3, 0); // color image with pixels set to 
    
    // Create an error image
    CImg<unsigned char> diff(width, height, 1, 3, 0); 

		for (int x = 0; x < width; x++)
		{
			// simple ASCII progress bar
			std::cout << "."; 

			for (int y = 0; y < height; y++)
			{
        const Point query(x, y);

        double r = 0.0;
        double g = 0.0;
        double b = 0.0;

        std::list<Point>::iterator itc = points.end();
        for(std::list<Point>::iterator it = points.begin(); 
                                       it != points.end(); 
                                       it++){
          if(query.coincide(*it))
            itc = it;
        }

        if(itc != points.end()){
          r = (double)input(x, y, 0);
          g = (double)input(x, y, 1);
          b = (double)input(x, y, 2);
        }

        else{

           K::Point_2 p(x, y);
           Point_coordinate_vector coords;
           CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>,
                        K::FT, bool> result =
                        CGAL::natural_neighbor_coordinates_2(dt, p,
					              std::back_inserter(coords));

           for(int i = 0; i < coords.size(); i++){
             const double w = coords[i].second / result.second;
             const int xp = coords[i].first.x();
             const int yp = coords[i].first.y();

             r += w * (double)input(xp, yp, 0);
             g += w * (double)input(xp, yp, 1);
             b += w * (double)input(xp, yp, 2);
           }
        }
				// write pixel color
				output.atXY(x, y, 0) = (unsigned char)r;
				output.atXY(x, y, 1) = (unsigned char)g;
				output.atXY(x, y, 2) = (unsigned char)b;

        // write error image
        diff.atXY(x, y, 0) = (unsigned char)abs(r - input(x, y, 0));
        diff.atXY(x, y, 1) = (unsigned char)abs(g - input(x, y, 1));
        diff.atXY(x, y, 2) = (unsigned char)abs(b - input(x, y, 2));
			}
		}
		std::cout << "done" << std::endl;

		std::cout << "write image...";
		output.save_bmp(pOutput);
    diff.save_bmp("diff.bmp");
		std::cout << "done" << std::endl;

		return true;
	}
};

