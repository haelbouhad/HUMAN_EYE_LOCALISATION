#ifndef FILTER_HPP
#define FILTER_HPP

#define TAILLE 256

#include <string>
#include <iostream>
#include "CImg.h"
#include <math.h>

using namespace cimg_library;

typedef struct
{
	float x;
	float y;
	float r;
	float val;
    float * A;
}tab_pix;


int comparator(const tab_pix* a, const tab_pix* b);
int max(int a, int b);
int min(int a, int b);


class Filter {

    public : 
        Filter();
        Filter(std::string);
        ~Filter();
        CImg<float> GaussianBlur(CImg<unsigned char>,double);
        CImg<double> GaussianKernel(double = 1.0);
        int GetAngle(int, int);
        CImg<unsigned char> CannyEdges(CImg<float>, double, double, bool = false);
         CImg<unsigned char> modified_hough_cercles(CImg<unsigned char>, int = 10, int = 50);
         CImg<float> modified_hough_cercles2(CImg<unsigned char>, int = 1, int = 100);
        CImg<unsigned char> hough_cercles(CImg<unsigned char>, int);
        CImg<unsigned char> hough_creer_cercles(CImg<unsigned char>, float, int);
        void compute();

    private :
        std::string filename;
        double sigma;
        double lowerThreshold;
        double upperThreshold;
        int numberCircle;
        float radius;
        tab_pix* pixels;
};


#endif