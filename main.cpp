#include "filter.hpp"

using namespace std;
using namespace cimg_library;

CImg<double> GaussianKernel(double sigma = 1.0)
{
  CImg<double> resultIm(5,5,1,1,0);
  int midX = 3, midY = 3;
  cimg_forXY(resultIm,X,Y) {
    resultIm(X,Y) = std::ceil(256.0*(std::exp(-(midX*midX + midY*midY)/(2*sigma*sigma)))/(2*cimg::PI*sigma*sigma));
  }
  return resultIm;
}


CImg<float> ApplyGaussian(CImg<unsigned char> im,double sigma)
{
  CImg<float> smoothIm(im.width(),im.height(),1,1,0);

  //make gaussian kernel
  CImg<float> gk = GaussianKernel(sigma);
  //apply gaussian

  CImg_5x5(I,int);
  cimg_for5x5(im,X,Y,0,0,I,int) {
    float sum = 0;
    sum += gk(0,0)*Ibb + gk(0,1)*Ibp + gk(0,2)*Ibc + gk(0,3)*Ibn + gk(0,4)*Iba;
    sum += gk(1,0)*Ipb + gk(1,1)*Ipp + gk(1,2)*Ipc + gk(1,3)*Ipn + gk(1,4)*Ipa;
    sum += gk(2,0)*Icb + gk(2,1)*Icp + gk(2,2)*Icc + gk(2,3)*Icn + gk(2,4)*Ica;
    sum += gk(3,0)*Inb + gk(3,1)*Inp + gk(3,2)*Inc + gk(3,3)*Inn + gk(3,4)*Ina;
    sum += gk(4,0)*Iab + gk(4,1)*Iap + gk(4,2)*Iac + gk(4,3)*Ian + gk(4,4)*Iaa;
    smoothIm(X,Y) = sum/256;
  }
  return smoothIm;
}


/**
 * PURPOSE: aux. function used by CannyEdges to quantize an angle theta given by gradients, dx and dy
 *  into 0 - 7
 * PARAM: dx,dy - gradient magnitudes
 * RETURN int value between 0 and 7
 **/
int GetAngle(int dy,int dx)
{
  double angle = cimg::abs(std::atan2((double)dy,(double)dx));
  if ((angle >= -cimg::PI/8)&&(angle <= cimg::PI/8))//-pi/8 to pi/8 => 0
    return 0;
  else if ((angle >= cimg::PI/8)&&(angle <= 3*cimg::PI/8))//pi/8 to 3pi/8 => pi/4
    return 1;
  else if ((angle > 3*cimg::PI/8)&&(angle <= 5*cimg::PI/8))//3pi/8 to 5pi/8 => pi/2
    return 2;
  else if ((angle > 5*cimg::PI/8)&&(angle <= 7*cimg::PI/8))//5pi/8 to 7pi/8 => 3pi/4
    return 3;
  else if (((angle > 7*cimg::PI/8) && (angle <= cimg::PI)) ||
           ((angle <= -7*cimg::PI/8)&&(angle >= -cimg::PI))) //-7pi/8 to -pi OR 7pi/8 to pi => pi
    return 4;
  else return 0;
}

/**
 * PURPOSE: create an edge map of the given image with hysteresis using thresholds T1 and T2
 * PARAMS: CImg<float> im the image to perform edge detection on
 * 		   T1 lower threshold
 *         T2 upper threshold
 * RETURN CImg<unsigned char> edge map
 **/
CImg<unsigned char> CannyEdges(CImg<float> im, double T1, double T2, bool doHysteresis=false)
{
  CImg<unsigned char> edges(im);
  CImg<float> secDerivs(im);
  secDerivs.fill(0);
  edges.fill(0);
  CImgList<float> gradients = im.get_gradient("xy",1);
  int image_width = im.width();
  int image_height = im.height();

  cimg_forXY(im,X,Y) {
    double Gr = std::sqrt(std::pow((double)gradients[0](X,Y),2.0) + std::pow((double)gradients[1](X,Y),2.0));
    double theta = GetAngle(Y,X);
    //if Gradient magnitude is positive and X,Y within the image
    //take the 2nd deriv in the appropriate direction
    if ((Gr > 0)&&(X < image_width - 2)&&(Y < image_height - 2)) {
      if (theta == 0)
        secDerivs(X,Y) = im(X + 2,Y) - 2*im(X + 1,Y) + im(X,Y);
      else if (theta == 1)
        secDerivs(X,Y) = im(X + 2,Y + 2) - 2*im(X + 1,Y + 1) + im(X,Y);
      else if (theta == 2)
        secDerivs(X,Y) = im(X,Y + 2) - 2*im(X,Y + 1) + im(X,Y);
      else if (theta == 3)
        secDerivs(X,Y) = im(X + 2,Y + 2) - 2*im(X + 1,Y + 1) + im(X,Y);
      else if (theta == 4)
        secDerivs(X,Y) = im(X + 2,Y) - 2*im(X + 1,Y) + im(X,Y);
    }
  }
  //for each 2nd deriv that crosses a zero point and magnitude passes the upper threshold.
  //Perform hysteresis in the direction of the gradient, rechecking the gradient
  //angle for each pixel that meets the threshold requirement.  Stop checking when
  //the lower threshold is not reached.
  CImg_5x5(I,float);
  cimg_for5x5(secDerivs,X,Y,0,0,I,float) {
    if (   (Ipp*Ibb < 0) ||
           (Ipc*Ibc < 0)||
           (Icp*Icb < 0)   ) {
      double Gr = std::sqrt(std::pow((double)gradients[0](X,Y),2.0) + std::pow((double)gradients[1](X,Y),2.0));
      int dir = GetAngle(Y,X);
      int Xt = X, Yt = Y, delta_x = 0, delta_y=0;
      double GRt = Gr;
      if (Gr >= T2)
        edges(X,Y) = 255;
      //work along the gradient in one direction
      if (doHysteresis) {
        while ((Xt > 0) && (Xt < image_width - 1) && (Yt > 0) && (Yt < image_height - 1)) {
          switch (dir){
          case 0 : delta_x=0;delta_y=1;break;
          case 1 : delta_x=1;delta_y=1;break;
          case 2 : delta_x=1;delta_y=0;break;
          case 3 : delta_x=1;delta_y=-1;break;
          case 4 : delta_x=0;delta_y=1;break;
          }
          Xt += delta_x;
          Yt += delta_y;
          GRt = std::sqrt(std::pow((double)gradients[0](Xt,Yt),2.0) + std::pow((double)gradients[1](Xt,Yt),2.0));
          dir = GetAngle(Yt,Xt);
          if (GRt >= T1)
            edges(Xt,Yt) = 255;
        }
        //work along gradient in other direction
        Xt = X; Yt = Y;
        while ((Xt > 0) && (Xt < image_width - 1) && (Yt > 0) && (Yt < image_height - 1)) {
          switch (dir){
          case 0 : delta_x=0;delta_y=1;break;
          case 1 : delta_x=1;delta_y=1;break;
          case 2 : delta_x=1;delta_y=0;break;
          case 3 : delta_x=1;delta_y=-1;break;
          case 4 : delta_x=0;delta_y=1;break;
          }
          Xt -= delta_x;
          Yt -= delta_y;
          GRt = std::sqrt(std::pow((double)gradients[0](Xt,Yt),2.0) + std::pow((double)gradients[1](Xt,Yt),2.0));
          dir = GetAngle(Yt,Xt);
          if (GRt >= T1)
            edges(Xt,Yt) = 255;
        }
      }
    }
  }
  return edges;
}

int main(int argc, char** argv){
    
    
    string filename = (argc == 2) ? argv[1] : "img/vue.pgm";

    /*
    Filter f(filename);
    
    cout << f.getFilename() << endl;

    f.compute();
    */

    CImg<unsigned char> image(filename.c_str());
    CImgDisplay input(image, "Image de depart");

    CImg<float> im_gauss = ApplyGaussian(image,2);
    CImgDisplay gauss(im_gauss, "Image apres le filtre gaussien");

    CImg<unsigned char> im_canny = CannyEdges(im_gauss, 2.5, 7.5, true);
    CImgDisplay canny(im_canny, "Image apres Canny"); 



    while(!canny.is_closed()){
        canny.wait();
    }

}