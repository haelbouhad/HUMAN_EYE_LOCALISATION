#ifndef FILTER_HPP
#define FILTER_HPP

#define TAILLE 256

#include <string>
#include <iostream>
#include "CImg.h"
#include <math.h>
#include <map>

using namespace cimg_library;

typedef struct
{
	float x;
	float y;
	float r;
	float val;
    std::map<int, int> Acc;     // Accumuler les rayons pour chaque point de centre : <rayon, occurence>
    std::map<int, int> AccLid;  // Accumuler les eyelids  pour les points choisis   : <y-cord, occurence>
}tab_pix;


int comparator(const tab_pix* a, const tab_pix* b);
int comparator_rayon(const tab_pix* a, const tab_pix* b);
int max_array(std::map<int, int>);
int max_array_key(std::map<int, int> &);



class Filter {

    public : 
        Filter();
        Filter(std::string);
        ~Filter();
        int GetAngle(int, int);
        CImg<unsigned char> CannyEdges(CImg<float>, double, double, bool = false);
         CImg<unsigned char> modified_hough_cercles(CImg<unsigned char>, int = 50, int = 70);
         CImg<unsigned char> modified_hough_creer_cercles(CImg<unsigned char>, float);
         void applyBresenham(CImg<unsigned char> &, int , int, int);
         void applyBresenhamForEyeLid(CImg<unsigned char> & , tab_pix &);
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