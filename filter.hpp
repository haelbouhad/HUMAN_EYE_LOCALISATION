#ifndef FILTER_HPP
#define FILTER_HPP

#define TAILLE 256

#include <string>
#include <iostream>
#include "CImg.h"


class Filter {

    public : 
        Filter();
        Filter(std::string);
        ~Filter();
        std::string getFilename();
        void compute();

    private :
        std::string filename;

};


#endif