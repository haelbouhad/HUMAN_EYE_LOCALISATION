#include "filter.hpp"

using namespace std;
using namespace cimg_library;

Filter::Filter(){

}

Filter::Filter(string pfilename):filename(pfilename){}

Filter::~Filter() {
    delete pixels;
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

CImg<unsigned char> Filter::modified_hough_cercles(CImg<unsigned char> img, int rayon1, int rayon2){
  CImgList<float> grad = img.get_gradient("xy",1);
  CImg<float> angles   = CImg<float>(img.width(), img.height());
  CImg<unsigned char> tab_cercles = CImg<unsigned char>(img.width(), img.height());
  
  angles = grad[1].get_atan2(grad[0]); // Compute atan2(y,x) for each pixel value.
  //int rayonpsi[2];
  int cpt;
  // Define pixels
  pixels = new tab_pix[img.width()*img.height()];
  for (int i = 0; i < img.width(); i++)
    for (int j = 0; j < img.height(); j++)
      (pixels[i+j*img.width()].Acc).insert( pair<int,int>(0,0) );
  int count = 0;
  
  for (int a = 0; a < img.width(); a++)
  for (int b = 0; b < img.height(); b++)
    if (img(a, b) != 0)
    {

          for(int r = rayon1; r <= rayon2; ++r){
              for(int sign = -1; sign <= 1; sign = sign+2 ){
                  int xc = a + r*cos(sign*M_PI/2 + angles(a,b));
                  int yc = b + r*sin(sign*M_PI/2 + angles(a,b));
                  cpt = xc+yc*img.width();
                  if( 0 <= xc && xc < img.width() && 0 <= yc && yc <= img.height() && cpt <= img.width()*img.height())
                  {
                    pixels[cpt].val++;
                    pixels[cpt].x = xc;
                    pixels[cpt].y = yc;
                    pixels[cpt].r = r;
                    count++;
                  }
              }
          }

    }

    for (int i = 0; i < img.width(); i++)
        for (int j = 0; j < img.height(); j++)
                tab_cercles(i,j) = pixels[i+j*img.width()].val;


     //std::cout << "COUNT = " << count << std::endl;

    qsort(pixels, img.width()*img.height(), sizeof(tab_pix), (int (*)(const void*, const void*))comparator);

    /*std::cout << "COUNT = " << count << std::endl;
    for(int i = 0; i < 10 ; ++i)
        std::cout << "MAX   = " << pixels[i].x << " - " <<  pixels[i].y << " - " << pixels[i].val << std::endl;
      */  
    
    
    int seuil = img.width()*img.height();
    int ray_cur;
    for (int s = 0; s < (int)seuil; s++)
        if(pixels[s].val > 1)
            for (int i = 0; i < img.width(); i++)
                for (int j = 0; j < img.height(); j++)
                {
                    if(img(i,j) != 0){
                        ray_cur = (i-pixels[s].x)*(i-pixels[s].x) + (j-pixels[s].y)*(j-pixels[s].y);
                        if (ray_cur >= rayon1*rayon1 && ray_cur <= rayon2*rayon2)
                        {
                            
                            if ( (pixels[s].Acc).find((int)sqrt(ray_cur)) == (pixels[s].Acc).end() ) {
                                (pixels[s].Acc).insert( std::pair<int,int>((int)sqrt(ray_cur),1));
                            } else {
                                pixels[s].Acc[(int)sqrt(ray_cur)]++;
                            } 
                            //cout << (int)sqrt(ray_cur) << endl;                
                        }
                    }
                }

    for (int s = 0; s < (int)seuil; s++)
      pixels[s].r = max_array_key(pixels[s].Acc);

    qsort(pixels, img.width()*img.height(), sizeof(tab_pix), (int (*)(const void*, const void*))comparator_rayon);

    /*for(int i = 0; i < 100 ; ++i){
        std::cout << "MAX   = x :" << pixels[i].x << " - y : " <<  pixels[i].y << " - val : " << pixels[i].val << endl;
        cout << " - " << " r : " << max_array(pixels[i].Acc) << std::endl;
    }*/

    // Eyelid detection
    int seuil2 = 2;
	int count1 = 0;
	int ray_eye, ray_eyelid;
    for(int s = 0 ; s < seuil2 ; s++){
      for(int other = seuil2; other < img.width()*img.height()/100; other++ ){
		  if(pixels[other].y >= pixels[s].y + 1*pixels[s].r && pixels[other].y <= pixels[s].y + 3*pixels[s].r ){

			  // Voisinage de l'oeil
			  for(int i = pixels[s].x - pixels[s].r  ; i <= pixels[s].x + pixels[s].r  ; i++)
			  	for(int j = pixels[s].y - pixels[s].r  ; j <= pixels[s].r ; j++){
				  ray_eyelid = (i - pixels[other].x)*(i - pixels[other].x) + (j - pixels[other].y)*(j - pixels[other].y);
				  // Appartenance à la première cercle (eyelid)
				  if( ray_eyelid  == pixels[other].r*pixels[other].r ){
					  ray_eye = (i - pixels[s].x)*(i - pixels[s].x) + (j - pixels[s].y)*(j - pixels[s].y);
					  if( ray_eye == pixels[s].r*pixels[s].r ){
						  if ( (pixels[s].AccLid).find((int)sqrt(ray_eyelid)) == (pixels[s].AccLid).end() ) {
                                (pixels[s].AccLid).insert( std::pair<int,int>( pixels[other].y - pixels[s].y ,1));
								cout << pixels[other].y - pixels[s].y  << endl;
                            } else {
                                pixels[s].AccLid[pixels[other].y - pixels[s].y]++;
								cout << pixels[other].y - pixels[s].y  << endl;
                         }  
					  }
				  }
				}

		  }
      }
    }




	return tab_cercles;
}


void Filter::applyBresenham(CImg<unsigned char> &  img, int r, int xc, int yc){

    int x,y,m;
    x = 0;
    y = r;
    m = 5 - 4*r;
    int c=0;
    while( x <= y){
        img( x+xc, y+yc ) = 255 ;
        img( y+xc, x+yc ) = 255 ;
        img( -x+xc, y+yc ) = 255 ;
        img( -y+xc, x+yc ) = 255 ;
        img( x+xc, -y+yc ) = 255 ;
        img( y+xc, -x+yc ) = 255 ;
        img( -x+xc, -y+yc ) = 255 ;
        img( -y+xc, -x+yc ) = 255 ;
        if(m > 0){
            y -= 1;
            m -= 8*y;
        }
        x += 1;
        m += 8*x + 4;
        c++;
    }

   // cout << "count = " << c << endl;

}

void Filter::applyBresenhamForEyeLid(CImg<unsigned char> &  img, tab_pix & pixel){

	int yvalue = max_array_key(pixel.AccLid);
	cout << yvalue << endl;
	applyBresenham(img, yvalue + pixel.r, pixel.x, yvalue + pixel.y);

}


CImg<unsigned char> Filter::modified_hough_creer_cercles(CImg<unsigned char> img, float seuil)
{
    cout << "Call draw circles" << endl;
	CImg<unsigned char> tab_cercles = CImg<unsigned char>(img.width(), img.height());
    tab_cercles.fill(0);
    
	for (int s = 0; s < (int)seuil; s++){
        pixels[s].r = max_array_key(pixels[s].Acc);
        applyBresenham(tab_cercles, pixels[s].r, pixels[s].x, pixels[s].y);
	}

	return tab_cercles;
}

CImg<float> Filter::modified_hough_cercles2(CImg<unsigned char> img, int rayon1, int rayon2){
  CImgList<float> grad = img.get_gradient("xy",1);
  CImg<float> angles   = CImg<float>(img.width(), img.height());

  return angles;
}

CImg<unsigned char> Filter::hough_cercles(CImg<unsigned char> img, int rayon1)
{
  // 
  CImgList<float> grad = img.get_gradient("xy",1);
  CImg<unsigned char> tab_cercles = CImg<unsigned char>(img.width(), img.height());
  CImg<float> angles = CImg<float>(img.width(), img.height());
  
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
           //if (img(i, j) != 0)
          {

            ray_cur = (i-a)*(i-a) + (j-b)*(j-b);

            cpt = i+j*img.width();
           // if (ray_cur >= rayon1-50 && ray_cur <= rayon1+50 && pixels[cpt].val < 255)
           if (ray_cur >= rayon1-2 && ray_cur <= rayon1+2 /*&& pixels[cpt].val < 255*/)
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

    //CImg<float> im_gauss = image.get_RGBtoYCbCr().get_channel(0);
    //im_gauss.blur(sigma);
    CImg<float> im_gauss = image.get_blur(sigma);
    CImgDisplay gauss(im_gauss, "Image apres le filtre gaussien");
    

    std::cout << " \n 2. Canny Edge detection method (with hysteresis) using thresholds T1 and T2 " << std::endl;
    std::cout << "    > Insert T1 value (lower threshold) : "; std::cin >> lowerThreshold;
    std::cout << "    > Insert T2 value (upper threshold) : "; std::cin >> upperThreshold;
    CImg<unsigned char> im_canny = CannyEdges(im_gauss, lowerThreshold, upperThreshold, true);
    CImgDisplay canny(im_canny, "Image apres Canny"); 

    std::cout << " \n 3. Hough Transform " << std::endl;
    /*std::cout << "    > Insert Number of circles (NC)  : "; std::cin >> numberCircle;
    std::cout << "    > Insert radius value : "; std::cin >> radius;*/
    CImg<float> im_hough = modified_hough_cercles(im_canny, 8, 12);
    CImgDisplay hough_fct(im_hough, "Image resultante de Hough");

    int seuil;
    std::cout << "    > Insert seuil value : "; std::cin >> seuil;
    im_hough = modified_hough_creer_cercles(im_hough, seuil);
    CImgDisplay hough_cercles(im_hough, "Image de cercles");
    


    while(!input.is_closed()){
        input.wait();
    }

}


int max_array(std::map<int, int> Acc){
    int max = 0;
    map<int, int>::iterator it;

    for ( it = Acc.begin(); it != Acc.end(); it++ )
    {
        if(it->second > max)
            max = it->second;
    }

    return max;
}

int max_array_key(std::map<int, int> & Acc){
	//cout << "Call max_array_key" << endl;
    int maxvalue = max_array(Acc);
    map<int, int>::iterator it;
    bool stop = false;
    int maxkey;

    for ( it = Acc.begin(); it != Acc.end() && !stop; it++ )
    {
        if(it->second == maxvalue){
            maxkey = it->first;
            stop = true;
        }
    }

    return maxkey;
}

int comparator(const tab_pix* a, const tab_pix* b)
{
	return ((*b).val - (*a).val);
}

int comparator_rayon(const tab_pix* a, const tab_pix* b)
{
	return (max_array((*b).Acc) - max_array((*a).Acc));
}


