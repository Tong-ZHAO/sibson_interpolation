#include <list>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iterator>
#include <Python.h>
#include <string>

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


bool pair_cmp(const std::pair<std::pair<int, int>, int>& firstElem, const std::pair<std::pair<int, int>, int>& secondElem) {
  return firstElem.second > secondElem.second;
}

class Sibson
{
public:
	Sibson() {}
	~Sibson() {}

public:
	bool run(const int nb_points,
    const int mode,
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

    if(mode == 1 || mode == 4){
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
    }

    int test_num = width * height - (int)points.size();
    std::cout << test_num << " points to calculate..."<< std::endl;

    Delaunay_triangulation dt;
    // do a delaunay triangulation 
    for (std::list<Point>::iterator it = points.begin(); it != points.end(); it++)
      dt.insert(K::Point_2(it->x(),it->y()));

    CImg<unsigned char> diff(width, height, 1, 3, 0);
    CImg<unsigned char> output(width, height, 1, 3, 0);

    float error[3] = {0.0};
      
	  // interpolation
	  std::cout << "interpolating...";
    
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

        error[0] += abs(r - input(x, y, 0));
        error[1] += abs(g - input(x, y, 1));
        error[2] += abs(b - input(x, y, 2));

        diff.atXY(x, y, 0) = (unsigned char)abs(r - input(x, y, 0));
        diff.atXY(x, y, 1) = (unsigned char)abs(g - input(x, y, 1));
        diff.atXY(x, y, 2) = (unsigned char)abs(b - input(x, y, 2));
		  }
		}
    std::cout << "done" << std::endl;

    error[0] = error[0] / test_num;
    error[1] = error[1] / test_num;
    error[2] = error[2] / test_num;

    std::cout << "The average error in R channel: " << error[0] << std::endl;
    std::cout << "The average error in G channel: " << error[1] << std::endl;
    std::cout << "The average error in B channel: " << error[2] << std::endl;


		std::cout << "write image...";
		output.save_bmp(pOutput);
    diff.save_bmp("diff.bmp");
		std::cout << "done" << std::endl << std::endl;

		return true;
	}
	
  bool run(const float max_error,
    const int mode,
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

    Delaunay_triangulation dt;
    // do a delaunay triangulation 
    for (std::list<Point>::iterator it = points.begin(); it != points.end(); it++)
      dt.insert(K::Point_2(it->x(),it->y()));
    
    float error[3] = {256.0};
    // Create an error image
    CImg<unsigned char> diff(width, height, 1, 3, 0);
    CImg<unsigned char> output(width, height, 1, 3, 0); // color image with pixels set to 

    while(error[0] > max_error || error[1] > max_error || error[2] > max_error){

      if(mode == 2){
        int x = random_int(0, width - 1);
        int y = random_int(0, height - 1);    
        points.push_back(Point(x, y));
        dt.insert(K::Point_2(x, y));
        std::cout << "add one point..." << std::endl;
      }
      else{
        int add_num = (int)points.size();
        for(int i = 0; i < add_num; i++){
          int x = random_int(0, width - 1);
          int y = random_int(0, height - 1);
          points.push_back(Point(x, y));
          dt.insert(K::Point_2(x, y));
        }

        std::cout << "add "<< add_num << " points..." << std::endl;
      }

		  // interpolation
      error[0] = 0.0;
      error[1] = 0.0;
      error[2] = 0.0;

      int test_num = width * height - (int)points.size();
		  std::cout << "interpolating...";
    
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

          error[0] += abs(r - input(x, y, 0));
          error[1] += abs(g - input(x, y, 1));
          error[2] += abs(b - input(x, y, 2));

          diff.atXY(x, y, 0) = (unsigned char)abs(r - input(x, y, 0));
          diff.atXY(x, y, 1) = (unsigned char)abs(g - input(x, y, 1));
          diff.atXY(x, y, 2) = (unsigned char)abs(b - input(x, y, 2));
			  }
		  }
      std::cout << "done" << std::endl;

      error[0] = error[0] / test_num;
      error[1] = error[1] / test_num;
      error[2] = error[2] / test_num;

      std::cout << "The average error in R channel: " << error[0] << std::endl;
      std::cout << "The average error in G channel: " << error[1] << std::endl;
      std::cout << "The average error in B channel: " << error[2] << std::endl;

    }

		std::cout << "write image...";
		output.save_bmp(pOutput);
    diff.save_bmp("diff.bmp");
		std::cout << "done" << std::endl << std::endl;

		return true;
	}


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
    std::ofstream myfile;

    // change it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    myfile.open("/home/zt/Maillage/sibson_interpolation/errors.txt");

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

    int test_num = width * height - (int)points.size();
    std::cout << test_num << " points to calculate..."<< std::endl;

    float p_error = 256.0 * 3;

    CImg<unsigned char> diff(width, height, 1, 3, 0);
    CImg<unsigned char> output(width, height, 1, 3, 0);
 
    float min_error = 256.0 * 3;
    std::list<Point> best_result;
    bool flag = false;
    int counter = 0;

    while(true){

      Delaunay_triangulation dt;
      // do a delaunay triangulation 
      for (std::list<Point>::iterator it = points.begin(); it != points.end(); it++)
        dt.insert(K::Point_2(it->x(),it->y()));

      float error[3] = {0.0};
      std::vector<std::pair<std::pair<int, int>, int> > pt_error;
      
	    // interpolation
	    std::cout << "interpolating...";
    
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

          int err_r = abs(r - input(x, y, 0));
          int err_g = abs(g - input(x, y, 1));
          int err_b = abs(b - input(x, y, 2));

          pt_error.push_back(std::make_pair(std::make_pair(x, y), err_r + err_g + err_b));

          error[0] += err_r;
          error[1] += err_g;
          error[2] += err_b;

          diff.atXY(x, y, 0) = (unsigned char)err_r;
          diff.atXY(x, y, 1) = (unsigned char)err_g;
          diff.atXY(x, y, 2) = (unsigned char)err_b;
		    }
		  }
    
      std::cout << "done" << std::endl;

      error[0] = error[0] / test_num;
      error[1] = error[1] / test_num;
      error[2] = error[2] / test_num;

      std::cout << "The average error in R channel: " << error[0] << std::endl;
      std::cout << "The average error in G channel: " << error[1] << std::endl;
      std::cout << "The average error in B channel: " << error[2] << std::endl;

      myfile << error[0] << " " << error[1] << " " << error[2] << " ";

      if(flag == true) break;

      // early stop
      if(error[0] + error[1] + error[2] < min_error){
        counter = 0;
        min_error = error[0] + error[1] + error[2];
        best_result = points;

        std::sort(pt_error.begin(), pt_error.end(), pair_cmp);
        p_error = error[0] + error[1] + error[2];

        for(int i = 0; i < nb_points / 100.0; i++){
          points.pop_front();
          points.push_back(Point(pt_error[i].first.first, pt_error[i].first.second));
        }

        points.push_front(Point(0.0, 0.0));
        points.push_front(Point(width - 1.0, 0.0));
        points.push_front(Point(0.0, height - 1.0));
        points.push_front(Point(width - 1.0, height - 1.0));
      }
      else if(counter < 10){
        counter++;
        std::sort(pt_error.begin(), pt_error.end(), pair_cmp);
        p_error = error[0] + error[1] + error[2];

        for(int i = 0; i < nb_points / 100.0; i++){
          points.pop_front();
          points.push_back(Point(pt_error[i].first.first, pt_error[i].first.second));
        }

        points.push_front(Point(0.0, 0.0));
        points.push_front(Point(width - 1.0, 0.0));
        points.push_front(Point(0.0, height - 1.0));
        points.push_front(Point(width - 1.0, height - 1.0));
      }
      else{
        points = best_result;
        flag = true;
      }
    }

    myfile.close();
		std::cout << "write image...";
		output.save_bmp(pOutput);
    diff.save_bmp("diff.bmp");
		std::cout << "done" << std::endl << std::endl;

    std::cout << "use python to draw figure..." << std::endl;
    FILE * fp = NULL;
    fp = fopen("/home/zt/Maillage/sibson_interpolation/draw_figure.py", "r");
    Py_Initialize();
    PyRun_SimpleFile(fp, "/home/zt/Maillage/sibson_interpolation/draw_figure.py");
    Py_Finalize();

		return true;
	}

};

