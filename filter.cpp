#include "filter.hpp"

using namespace std;
using namespace cimg_library;

Filter::Filter(){

}

Filter::Filter(string pfilename){

    filename = (pfilename.size() ? pfilename : "img/vue.pgm" );

}

string Filter::getFilename(){
    return filename;
}

Filter::~Filter() {}


void Filter::compute(){

    CImg<unsigned char> image(filename.c_str());
    CImgDisplay input(image, "Image de depart");
    while(!input.is_closed()){
        input.wait();
    }

}