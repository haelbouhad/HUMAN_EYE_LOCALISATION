#include "filter.hpp"

using namespace std;
using namespace cimg_library;



int main(int argc, char** argv){
    
    
    string filename = (argc == 2) ? argv[1] : "img/hough-orginal.jpg";

    Filter f(filename);

    f.compute();
    
}