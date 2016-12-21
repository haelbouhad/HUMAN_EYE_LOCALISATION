#include "filter.hpp"

using namespace std;
using namespace cimg_library;



int main(int argc, char** argv){
    
    
    string filename = (argc == 2) ? argv[1] : "img/cercle.tif";

    Filter f(filename);

    f.compute();
    
    /*
    CImg<unsigned char> image(filename.c_str());
    CImgDisplay input(image, "Image de depart");

    CImg<float> im_gauss = ApplyGaussian(image,4);
    CImgDisplay gauss(im_gauss, "Image apres le filtre gaussien");

    CImg<unsigned char> im_canny = CannyEdges(im_gauss, 2.5, 7.5, true);
    CImgDisplay canny(im_canny, "Image apres Canny"); 
    //im_canny.normalize(0, 255);
    //im_canny.save("img/out_canny.bmp");

    CImg<unsigned char> im_hough = hough_cercles(im_canny, 60);
    CImgDisplay hough_fct(im_hough, "Image resultante de Hough");

    im_hough = hough_creer_cercles(im_hough, 600, 60);
    CImgDisplay hough_cercles(im_hough, "Image de cercles");



    while(!canny.is_closed()){
        canny.wait();
    }
*/
}