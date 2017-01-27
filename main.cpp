#include "filter.hpp"

using namespace std;
using namespace cimg_library;



int main(int argc, char** argv){
    
    
    string filename = (argc == 2) ? argv[1] : "img/hough-orginal.jpg";

    Filter f(filename);

    f.compute();

    /*
    CImg<unsigned char> image(filename.c_str());
    CImgDisplay input(image, "Image de depart");
    CImg<unsigned char> image2 = f.applyBresenham(image,  90, 160, 160);
    CImgDisplay output(image2, "Image de depart");
    while(!input.is_closed()){
        input.wait();
    }
    */
    
}