
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <cv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#define scast(x) saturate_cast<uchar>(x)

#define PI 3.14159265359

using namespace cv;
using namespace std;

Mat image, result ;

double A[500][500], B[500][500], temp[500][500], tmpa[500][500], tmpb[500][500] ;

int rows, cols ; 


int dft(int u0, int v0)
{	
	int i, j, k, l, fa, fb ;
	
	double a, b, rad, tmp1, tmp2 ; // a + ib
	double shft ;
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			a = b = 0 ; 
		
			for(k=0;k<rows;k++)
			{
				for(l=0;l<cols;l++)
				{
					tmp1 = (double) (i*k) / (double) (rows) ;
					tmp2 = (double) (j*l) / (double) (cols) ;
				
					rad = 2 * PI * (tmp1 + tmp2) ; 
					rad = -rad ;
					
					fa = image.at<uchar>(k, l) ;
					fb = 0 ; 
					
					shft = 2 * PI * ((double) (u0*k) / (double) (rows) + (double) (v0*l) / (double) (cols)) ;
					
					fa = fa * cos(shft) ;
					fb = fa * sin(shft) ;
					
					a+= (fa * cos(rad)) - (fb * sin(rad)) ;
					b+= (fa * sin(rad)) + (fb * cos(rad)) ; 
				}
			}
			
			A[i][j] = a ;
			B[i][j] = b ;
		}
	}
}



int dftshift(int u0, int v0)
{	
	int i, j ; 
	
	while(u0 < 0)
	{
		u0+= rows ; 
	}
	
	while(v0 < 0)
	{
		v0+= cols ; 
	}    

	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{	
			tmpa[i][j] = A[(i+u0)%rows][(j+v0)%cols] ; 
			tmpb[i][j] = B[(i+u0)%rows][(j+v0)%cols] ; 
		}
	}
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{	
			A[i][j] = tmpa[i][j] ; 
			B[i][j] = tmpb[i][j] ; 
		}
	}
}



int show_fourier_spectrum()
{
	Mat res ;

	res.create(rows, cols, CV_8U);
    
	int i, j ;
	
	double la, lb, mod, mi, ma, sca ;
	
	mi = 99999999 ;
	ma = -99999999 ;
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			la = A[i][j] ; 
			lb = B[i][j] ;
			
			mod = sqrt(la*la + lb*lb) ;
		 
		 	temp[i][j] = mod ;
		 	
		 	if(mi > mod)
		 		mi = mod ;
		 		
		 	if(ma < mod)
		 		ma = mod ;
		}
	}
	
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			if(mi != ma)
			{
				sca = (temp[i][j] - mi) / (ma - mi) ; 
				sca = 5 * log(1 + sca) ;
				sca*= 255 ;
				
				sca = scast(sca) ; 
				
				res.at<uchar>(i,j) = sca ;
			}
		}
	}

	namedWindow( "Display", WINDOW_NORMAL );
    imshow( "Display", res );                                   
    
	imwrite("res.png", res);
    
    waitKey(0);                                      
} 


int show_power_spectrum()
{
	Mat res ;

	res.create(rows, cols, CV_8U);
    
	int i, j ;
	
	double la, lb, mod, mi, ma, sca ;
	
	mi = 99999999 ;
	ma = -99999999 ;
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			la = A[i][j] ; 
			lb = B[i][j] ;
			
			mod = la*la + lb*lb ;
		 
		 	temp[i][j] = mod ;
		 	
		 	if(mi > mod)
		 		mi = mod ;
		 		
		 	if(ma < mod)
		 		ma = mod ;	
		}
	}
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			if(mi != ma)
			{
				sca = (temp[i][j] - mi) / (ma - mi) ; 
				sca = 10 * log(1 + sca) ;
				sca*= 255 ;
				
				sca = scast(sca) ; 
				
				res.at<uchar>(i,j) = sca ;
			}
		}
	}

	namedWindow( "Display", WINDOW_NORMAL );
    imshow( "Display", res );                                   
    
	imwrite("res.png", res);
    
    waitKey(0);                                      
} 


int show_phase()
{
	Mat res ;

	res.create(rows, cols, CV_8U);
    
	int i, j ;
	
	double la, lb, mod, mi, ma, sca ;
	
	mi = 99999999 ;
	ma = -99999999 ;
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			la = A[i][j] ; 
			lb = B[i][j] ;
			
			if(la != 0)
			{
				mod = atan(lb / la) ; 
			}  
			else
			{
				if(lb > 0)
					mod = PI / 2 ;
				
				if(lb < 0) 
					mod = -PI / 2 ; 
			}
		 
		 	temp[i][j] = mod ;
		 	
		 	if(mi > mod)
		 		mi = mod ;
		 		
		 	if(ma < mod)
		 		ma = mod ;	
		}
	}
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			if(mi != ma)
			{
				sca = (temp[i][j] - mi) / (ma - mi) ; 
				
				sca*= 255 ;
				
				sca = scast(sca) ; 
				
				res.at<uchar>(i,j) = sca ;
			}
		}
	}

	namedWindow( "Display", WINDOW_NORMAL );
    imshow( "Display", res );                                   
    
    
	imwrite("res.png", res);
	 	
    waitKey(0);                                      
} 


int idft()
{
	result.create(rows, cols, CV_8U);
	
	int i, j, k, l ;
	
	uchar finval ;
	
	double la, lb, a, b, rad, tmp1, tmp2, mn ; // a + ib
	
	mn = rows*cols ;  
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			a = b = 0 ; 
		
			for(k=0;k<rows;k++)
			{
				for(l=0;l<cols;l++)
				{
					tmp1 = (double) (i*k) / (double) (rows) ;
					tmp2 = (double) (j*l) / (double) (cols) ;
				
					rad = 6.28318530718 * (tmp1 + tmp2) ;
					
					la = A[k][l] ; 
					lb = B[k][l] ; 
					
					a+= (la * cos(rad)) - (lb * sin(rad)) ;
					b+= (la * sin(rad)) + (lb * cos(rad)) ;
				}
			}
			
			a/= mn ;
			b/= mn ; 
			
			finval = scast(a) ;  
			
			result.at<uchar>(i,j) = finval ;	
		}
	}
}



int gaussian_low(double cutoff)
{
	dftshift(rows/2, cols/2) ; 

	int i, j ; 
	
	double dist, factor ; 
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			dist = (rows/2 - i)*(rows/2 - i) + (cols/2 - j)*(cols/2 - j) ; 
			
			factor = exp(-((dist*dist)/(2*cutoff*cutoff))) ; 
			
			tmpa[i][j] = A[i][j] ;
			tmpb[i][j] = B[i][j] ; 
			
			A[i][j]*= factor ; 
			B[i][j]*= factor ; 
		}
	}
	
	dftshift(-rows/2, -cols/2) ; 
}




int gaussian_high(double cutoff)
{
	dftshift(rows/2, cols/2) ; 

	int i, j ; 
	
	double dist, factor ; 
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			dist = (rows/2 - i)*(rows/2 - i) + (cols/2 - j)*(cols/2 - j) ; 
			
			factor = 1 - exp(-((dist*dist)/(2*cutoff*cutoff))) ; 
			factor+= 0.5 ;
			
			tmpa[i][j] = A[i][j] ;
			tmpb[i][j] = B[i][j] ; 
			
			A[i][j]*= factor ; 
			B[i][j]*= factor ; 
		}
	}
	
	dftshift(-rows/2, -cols/2) ; 
}



int butterworth_low(double cutoff, int order)
{
	dftshift(rows/2, cols/2) ; 

	int i, j ; 
	
	double dist, factor ; 
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			dist = (rows/2 - i)*(rows/2 - i) + (cols/2 - j)*(cols/2 - j) ; 
			
			factor = dist / cutoff ; 
			factor = 1 + pow(factor, 2*order) ; 
			factor = 1 / factor ;  
			
			tmpa[i][j] = A[i][j] ;
			tmpb[i][j] = B[i][j] ; 
			
			A[i][j]*= factor ; 
			B[i][j]*= factor ; 
		}
	}
	
	dftshift(-rows/2, -cols/2) ; 
}


int butterworth_high(double cutoff, int order)
{
	dftshift(rows/2, cols/2) ; 

	int i, j ; 
	
	double dist, factor ; 
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			dist = (rows/2 - i)*(rows/2 - i) + (cols/2 - j)*(cols/2 - j) ; 
			
			factor = (dist*dist) / (dist*dist + cutoff*cutoff) ; 
			factor+= 0.5 ;
			
			tmpa[i][j] = A[i][j] ;
			tmpb[i][j] = B[i][j] ; 
			
			A[i][j]*= factor ; 
			B[i][j]*= factor ; 
		}
	}
	
	dftshift(-rows/2, -cols/2) ; 
}


int restore()
{
	int i, j ;
	
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			A[i][j] = tmpa[i][j] ; 
			B[i][j] = tmpb[i][j] ;
		}
	}
}


int main(int argc, char** argv)
{
    if(argc < 2)
    {
    	cout << "Error : Not enough arguments" << endl ;
    	return -1;
    }

    image = imread(argv[1], 0) ;

    if(!image.data)
    {
        cout << "Error : Could not open or find the image" << endl ;
        return -1;
    }
    
    rows = image.rows ;
    cols = image.cols ;  
 
 
 	// Menu 
 	
 	int cont = 1, ch = 0, ord = 1 ;
 	
 	double cut = 0 ;
 	
 	dft(0, 0) ;
 	
 	while(cont > 0)  
 	{
	 	cout << "Please select from the options given below : " << endl << endl ; 
	 	cout << "1 : Compute Discrete Fourier Transform" << endl ; 
	 	cout << "2 : Display Fourier Spectrum" << endl ; 
	 	cout << "3 : Display Power Spectrum" << endl ; 
	 	cout << "4 : Display Phase Matrix" << endl ; 
	 	cout << "5 : Shift Fourier Transform" << endl ; 
	 	cout << "6 : Compute Inverse Discrete Fourier Transform" << endl ;
	 	cout << "7 : Apply Gaussian Low Pass Filter" << endl ;
	 	cout << "8 : Apply Gaussian High Pass Filter" << endl ;
	 	cout << "9 : Apply Butterworth Low Pass Filter" << endl ;
	 	cout << "10 : Apply Butterworth High Pass Filter" << endl ;
	 	cout << "11 : Show original image" << endl ;
	 	cout << "0 : Exit" << endl << endl ;
	 	
	 	cout << "Enter your choice : " ; 
	 	cin >> ch ; 
	 	
	 	switch(ch) 
	 	{
	 		case 1:
	    	dft(0, 0) ;	  // Some things here
	 		break ;
	 		
	 		case 2:
	 		dftshift(rows/2, cols/2) ; 
	 		show_fourier_spectrum() ;
	 		dftshift(-rows/2, -cols/2) ; 
	 		break ;
	 		
	 		case 3:
	 		dftshift(rows/2, cols/2) ;
	 		show_power_spectrum() ;
	 		dftshift(-rows/2, -cols/2) ;
	 		break ;
	 		
	 		case 4:
	 		dftshift(rows/2, cols/2) ;
	 		show_phase() ; 
	 		dftshift(-rows/2, -cols/2) ;
	 		break ;
	 		
	 		case 5:
	 		int u0, v0 ;
	 		
	 		cout << "How much to shift in Y direction : " ;
	 		cin >> u0 ;
	 		
	 		cout << "How much to shift in X direction : " ;
	 		cin >> v0 ;
	 		
    		dftshift(u0, v0) ; 
    		show_fourier_spectrum() ; 	
    		dftshift(-u0, -v0) ; 	
	 		break ;
	 		
	 		case 6:
	 		idft() ;
	 		
	 		namedWindow( "IDFT", WINDOW_NORMAL );
   			imshow( "IDFT", result );
   			
			imwrite("res.png", result);
   			
   			waitKey(0) ;
	 		break ;
	 		
	 		case 7:
	 		
	 		cout << "Enter cutoff frequency : " ;
	 		cin >> cut ;
	 		
	 		gaussian_low(cut) ;
	 		
	 		idft() ;
	 		
	 		namedWindow( "IDFT", WINDOW_NORMAL );
   			imshow( "IDFT", result );
   			
   			imwrite("res.png", result);
   			
   			waitKey(0) ;
   			
   			restore() ;
   			
	 		break ;
	 		
	 		case 8:
	 		
	 		cout << "Enter cutoff frequency : " ;
	 		cin >> cut ;
	 		
	 		gaussian_high(cut) ;
	 		
	 		idft() ;
	 		
	 		namedWindow( "IDFT", WINDOW_NORMAL );
   			imshow( "IDFT", result );
   			
   			imwrite("res.png", result);
   			
   			waitKey(0) ;
   			
   			restore() ;
   			
	 		break ;
	 		
	 		case 9: 
	 		
	 		cout << "Enter cutoff frequency : " ;
	 		cin >> cut ;
	 		cout << "Enter order of butterworth filter (n = ?) : " ;
	 		cin >> ord ;
	 		
	 		butterworth_low(cut, ord) ;
	 		
	 		idft() ;
	 		
	 		namedWindow( "IDFT", WINDOW_NORMAL );
   			imshow( "IDFT", result );
   			
   			imwrite("res.png", result);
   			
   			waitKey(0) ;
   			
   			restore() ;
   			
	 		break ;
	 		
	 		
	 		case 10:
	 		
	 		cout << "Enter cutoff frequency : " ;
	 		cin >> cut ;
	 		cout << "Enter order of butterworth filter (n = ?) : " ;
	 		cin >> ord ;
	 		
	 		butterworth_high(cut, ord) ;
	 		
	 		idft() ;
	 		
	 		namedWindow( "IDFT", WINDOW_NORMAL );
   			imshow( "IDFT", result );
   			
   			imwrite("res.png", result);
   			
   			waitKey(0) ;
   			
   			restore() ;
   			
	 		break ;
	 		
	 		case 11 :
	 		
	 		namedWindow( "Original", WINDOW_NORMAL );
   			imshow( "Original", image );
   			
   			imwrite("res.png", image);
   			
   			waitKey(0) ;
	 		break ;
	 		
	 		case 0:
	 		cont = 0 ;
	 		break ;		
	 		
	 		default:
	 		cout << "Error - Invalid option selected" << endl ; 
	 	}
	 	
	 	
	 	destroyAllWindows() ; 
	 	
	 	
	 	cout << endl << endl ; 
    }
                           
    
    return 0;
    
}









