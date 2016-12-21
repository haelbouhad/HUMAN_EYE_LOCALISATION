#include "filter.hpp"

using namespace std;
using namespace cimg_library;

Filter::Filter(){

}

Filter::Filter(string pfilename):filename(pfilename){}

Filter::~Filter() {
    delete pixels;
}

CImg<double> Filter::GaussianKernel(double sigma)
{
  CImg<double> resultIm(5,5,1,1,0);
  int midX = 3, midY = 3;
  cimg_forXY(resultIm,X,Y) {
    resultIm(X,Y) = std::ceil(256.0*(std::exp(-(midX*midX + midY*midY)/(2*sigma*sigma)))/(2*cimg::PI*sigma*sigma));
  }
  return resultIm;
}


CImg<float> Filter::GaussianBlur(CImg<unsigned char> im,double sigma)
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

int Filter::GetAngle(int dy,int dx)
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

CImg<unsigned char> Filter::CannyEdges(CImg<float> im, double T1, double T2, bool doHysteresis)
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

CImg<unsigned char> Filter::hough_cercles(CImg<unsigned char> img, int rayon1)
{
  CImg<unsigned char> tab_cercles = CImg<unsigned char>(img.width(), img.height());
  tab_cercles.fill(0);
  pixels = new tab_pix[img.width()*img.height()];
  int cpt, ray_cur;

  rayon1 *= rayon1;

  for (int i = 0; i < img.width(); i++)
    for (int j = 0; j < img.height(); j++)
      pixels[i+j*img.width()].val = 0;

  for (int a = 0; a < img.width(); a++)
    for (int b = 0; b < img.height(); b++)
      if (img(a, b) != 0)
      {
        for (int i = 0; i < img.width(); i++)
	 				for (int j = 0; j < img.height(); j++)
          {
            ray_cur = (i-a)*(i-a) + (j-b)*(j-b);

            cpt = i+j*img.width();
            if (ray_cur >= rayon1-50 && ray_cur <= rayon1+50 && pixels[cpt].val < 255)
            {
              pixels[cpt].val++;
              pixels[cpt].x = i;
              pixels[cpt].y = j;
							pixels[cpt].r = rayon1;
            }
					/*	if (ray_cur >= ray2-50 && ray_cur <= ray2+50 && pixels[cpt].val < 255)
            {
              pixels[cpt].val++;
              pixels[cpt].x = i;
              pixels[cpt].y = j;
							pixels[cpt].r = ray2;
            }*/
          }
        }

  for (int i = 0; i < img.width(); i++)
    for (int j = 0; j < img.height(); j++)
			tab_cercles(i,j) = pixels[i+j*img.width()].val;

	return tab_cercles;
 // tab_cercles.save_jpeg("img/hough.jpg");
}


CImg<unsigned char> Filter::hough_creer_cercles(CImg<unsigned char> img, float seuil, int rayon1)
{
	CImg<unsigned char> tab_cercles = CImg<unsigned char>(img.width(), img.height());
  tab_cercles.fill(0);
  qsort(pixels, img.width()*img.height(), sizeof(tab_pix), (int (*)(const void*, const void*))comparator);
  seuil = img.width()*img.height()*(float)seuil/100000;
	rayon1 *= rayon1;

	for (int s = 0; s < (int)seuil; s++)
		for (int i = 0; i < img.width(); i++)
			for (int j = 0; j < img.height(); j++)
			{
				if (pixels[s].r == rayon1 && (i-pixels[s].x)*(i-pixels[s].x)+(j-pixels[s].y)*(j-pixels[s].y) >= rayon1-100 && (i-pixels[s].x)*(i-pixels[s].x)+(j-pixels[s].y)*(j-pixels[s].y) <= rayon1+100)
					tab_cercles(i, j) = 255;
			/*	if (pixels[s].r == ray2 && (i-pixels[s].x)*(i-pixels[s].x)+(j-pixels[s].y)*(j-pixels[s].y) >= ray2-100 && (i-pixels[s].x)*(i-pixels[s].x)+(j-pixels[s].y)*(j-pixels[s].y) <= ray2+100)
					tab_cercles(i, j) = 255;*/
			}

  cout << "debug..." << endl;

	return tab_cercles;
  //tab_cercles.save_jpeg("img/finale.jpg");
}


void Filter::compute(){

    CImg<unsigned char> image(filename.c_str());
    CImgDisplay input(image, "Image de depart");

    std::cout << "\n\n **************** Welcome to HUMAN_EYE_LOCATION application ****************\n" << std::endl ;

    std::cout << "\n 1. Gaussian blur " << std::endl;
    std::cout << "    > Insert sigma value (gaussian equation parameter) : " ;
    std::cin >> sigma;
    CImg<float> im_gauss = GaussianBlur(image,4);
    CImgDisplay gauss(im_gauss, "Image apres le filtre gaussien");

    std::cout << " \n 2. Canny Edge detection method (with hysteresis) using thresholds T1 and T2 " << std::endl;
    std::cout << "    > Insert T1 value (lower threshold) : "; std::cin >> lowerThreshold;
    std::cout << "    > Insert T2 value (upper threshold) : "; std::cin >> upperThreshold;
    CImg<unsigned char> im_canny = CannyEdges(im_gauss, lowerThreshold, upperThreshold, true);
    CImgDisplay canny(im_canny, "Image apres Canny"); 

    std::cout << " \n 3. Hough Transform " << std::endl;
    std::cout << "    > Insert Number of circles (NC)  : "; std::cin >> numberCircle;
    std::cout << "    > Insert radius value : "; std::cin >> radius;
    CImg<unsigned char> im_hough = hough_cercles(im_canny, radius);
    CImgDisplay hough_fct(im_hough, "Image resultante de Hough");

    im_hough = hough_creer_cercles(im_hough, numberCircle, radius);
    CImgDisplay hough_cercles(im_hough, "Image de cercles");
    

    while(!input.is_closed()){
        input.wait();
    }

}

int comparator(const tab_pix* a, const tab_pix* b)
{
	return ((*b).val - (*a).val);
}