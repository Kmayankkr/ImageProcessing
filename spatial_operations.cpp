

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <stdlib.h>
#include <algorithm> 
#include <math.h>

#define scast(x) saturate_cast<uchar>(x)

using namespace cv;
using namespace std;


Mat *image, *negative, *hist_eq, *scale, *filter, *rota, *trans, *shr, *ad_hist_eq, *bitplane, *hist_match, *sample, *spatial, *piece, *power, *loga   ;

double b_hist[255], g_hist[255], r_hist[255] ; 


Vec3b _bilinear(double rc, double cc)
{
	int r1, r2, c1, c2 ;

	r1 = floor(rc) ;
	r2 = ceil(rc) ;
	
	c1 = floor(cc) ;
	c2 = ceil(cc) ; 
	
	Vec3b i11, i12, i21, i22, ifin ;
	
	double iup[3], idown[3], icen[3] ;
	
	i11 = image->at<Vec3b>(r1, c1);
	i12 = image->at<Vec3b>(r1, c2);
	i21 = image->at<Vec3b>(r2, c1);
	i22 = image->at<Vec3b>(r2, c2); 
	
	if(c1 != c2)
	{
	iup[0] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i11[0] + ((double) (cc - c1) / (double) (c2 - c1)) * i12[0]) ; 
	iup[1] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i11[1] + ((double) (cc - c1) / (double) (c2 - c1)) * i12[1]) ;
	iup[2] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i11[2] + ((double) (cc - c1) / (double) (c2 - c1)) * i12[2]) ;
	
	idown[0] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i21[0] + ((double) (cc - c1) / (double) (c2 - c1)) * i22[0]) ; 
	idown[1] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i21[1] + ((double) (cc - c1) / (double) (c2 - c1)) * i22[1]) ;
	idown[2] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i21[2] + ((double) (cc - c1) / (double) (c2 - c1)) * i22[2]) ;
	}
	else
	{
	iup[0] = i11[0] ;
	iup[1] = i11[1] ;
	iup[2] = i11[2] ;
	
	idown[0] = i21[0] ;
	idown[1] = i21[1] ;
	idown[2] = i21[2] ;
	}
	
	if(r1 != r2)
	{
	icen[0] = (double)(((double) (r2 - rc) / (double) (r2 - r1)) * iup[0] + ((double) (rc - r1) / (double) (r2 - r1)) * idown[0]) ; 
	icen[1] = (double)(((double) (r2 - rc) / (double) (r2 - r1)) * iup[1] + ((double) (rc - r1) / (double) (r2 - r1)) * idown[1]) ;
	icen[2] = (double)(((double) (r2 - rc) / (double) (r2 - r1)) * iup[2] + ((double) (rc - r1) / (double) (r2 - r1)) * idown[2]) ;
	}
	else
	{
	icen[0] = iup[0] ;
	icen[1] = iup[1] ;
	icen[2] = iup[2] ;
	}
	
	ifin[0] = icen[0] ;
	ifin[1] = icen[1] ;
	ifin[2] = icen[2] ;  
			
	return ifin ;
}


void _negative()
{
	int i, j ;
    
    negative = new Mat(image->rows, image->cols, CV_8UC3, CV_RGB(0,0,0));
    
    for(i=0;i<image->rows;i++)
    {
    	for(j=0;j<image->cols;j++)
    	{
    		Vec3b intensity = image->at<Vec3b>(i, j);
			intensity.val[0] = 255 - intensity.val[0];
			intensity.val[1] = 255 - intensity.val[1];
			intensity.val[2] = 255 - intensity.val[2];
			negative->at<Vec3b>(i, j) = intensity; 
     	}
    }
}

// power, log, same as above 


void _histogram(Mat *ptr)
{	
	int i, j ;

	memset (b_hist, 0, sizeof(b_hist)) ;
	memset (g_hist, 0, sizeof(g_hist)) ;
	memset (r_hist, 0, sizeof(r_hist)) ;

	for(i=0;i<ptr->rows;i++)
    {
    	for(j=0;j<ptr->cols;j++)
    	{
    		Vec3b intensity = ptr->at<Vec3b>(i, j);
			b_hist[intensity.val[0]]++ ;
			g_hist[intensity.val[1]]++ ;
			r_hist[intensity.val[2]]++ ;
     	}
    }
}


void _hist_eq()
{
	_histogram(image) ; 
	
	int i, j ; 
	
	int total = (image->rows) * (image->cols) ;
	
	for(i=1;i<256;i++)
	{
		b_hist[i]+= b_hist[i-1] ;
		g_hist[i]+= g_hist[i-1] ;
		r_hist[i]+= r_hist[i-1] ;
	}
	 
 	for(i=0;i<256;i++)
	{
		b_hist[i] = (b_hist[i] * 255) / total ;
		g_hist[i] = (g_hist[i] * 255) / total ;
		r_hist[i] = (r_hist[i] * 255) / total ;
	}
	
	hist_eq = new Mat(image->rows, image->cols, CV_8UC3, CV_RGB(0,0,0));
    
    for(i=0;i<image->rows;i++)
    {
    	for(j=0;j<image->cols;j++)
    	{
    		Vec3b intensity = image->at<Vec3b>(i, j);
			intensity.val[0] = b_hist[intensity.val[0]];
			intensity.val[1] = g_hist[intensity.val[1]];
			intensity.val[2] = r_hist[intensity.val[2]];
			hist_eq->at<Vec3b>(i, j) = intensity; 
     	}
    }
}


int downsize(double rf, double cf, int interpolate)    // 0 - nearest neighbour 1 - bilinear
{
	int nrow = image->rows * rf ;
	int ncol = image->cols * cf ; 
	
	scale = new Mat(nrow, ncol, CV_8UC3, CV_RGB(0,0,0));
	
	int i, j ;
	
	
	
	if(interpolate == 1)
	{
	double rc, cc ;
	
	int r1, r2, c1, c2 ;
	
	for(i=0;i<nrow;i++)
	{
		for(j=0;j<ncol;j++)
		{
			rc = (double) i / rf ;
			cc = (double) j / cf ; 
			
			scale->at<Vec3b>(i, j) = _bilinear(rc, cc) ;  
			
			/*
			
			r1 = floor(rc) ;
			r2 = ceil(rc) ;
			
			c1 = floor(cc) ;
			c2 = ceil(cc) ; 
			
			Vec3b i11, i12, i21, i22, ifin ;
			
			double iup[3], idown[3], icen[3] ;
			
			i11 = image->at<Vec3b>(r1, c1);
			i12 = image->at<Vec3b>(r1, c2);
			i21 = image->at<Vec3b>(r2, c1);
			i22 = image->at<Vec3b>(r2, c2); 
			
			if(c1 != c2)
			{
			iup[0] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i11[0] + ((double) (cc - c1) / (double) (c2 - c1)) * i12[0]) ; 
			iup[1] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i11[1] + ((double) (cc - c1) / (double) (c2 - c1)) * i12[1]) ;
			iup[2] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i11[2] + ((double) (cc - c1) / (double) (c2 - c1)) * i12[2]) ;
			
			idown[0] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i21[0] + ((double) (cc - c1) / (double) (c2 - c1)) * i22[0]) ; 
			idown[1] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i21[1] + ((double) (cc - c1) / (double) (c2 - c1)) * i22[1]) ;
			idown[2] = (double)(((double) (c2 - cc) / (double) (c2 - c1)) * i21[2] + ((double) (cc - c1) / (double) (c2 - c1)) * i22[2]) ;
			}
			else
			{
			iup[0] = i11[0] ;
			iup[1] = i11[1] ;
			iup[2] = i11[2] ;
			
			idown[0] = i21[0] ;
			idown[1] = i21[1] ;
			idown[2] = i21[2] ;
			}
			
			if(r1 != r2)
			{
			icen[0] = (double)(((double) (r2 - rc) / (double) (r2 - r1)) * iup[0] + ((double) (rc - r1) / (double) (r2 - r1)) * idown[0]) ; 
			icen[1] = (double)(((double) (r2 - rc) / (double) (r2 - r1)) * iup[1] + ((double) (rc - r1) / (double) (r2 - r1)) * idown[1]) ;
			icen[2] = (double)(((double) (r2 - rc) / (double) (r2 - r1)) * iup[2] + ((double) (rc - r1) / (double) (r2 - r1)) * idown[2]) ;
			}
			else
			{
			icen[0] = iup[0] ;
			icen[1] = iup[1] ;
			icen[2] = iup[2] ;
			}
			
			ifin[0] = icen[0] ;
			ifin[1] = icen[1] ;
			ifin[2] = icen[2] ;  
			
			scale->at<Vec3b>(i, j) = ifin ;  
			
			*/ 
		}
	}
	
	}
	else
	{
	double rc, cc ;
	
	double r1, r2, c1, c2 ;
	
	int rfin, cfin ;
	
	for(i=0;i<nrow;i++)
	{
		for(j=0;j<ncol;j++)
		{
			rc = (double) i / rf ;
			cc = (double) j / cf ; 
			
			r1 = floor(rc) ;
			r2 = ceil(rc) ;
			
			c1 = floor(cc) ;
			c2 = ceil(cc) ; 
			
			if((rc-r1) < (r2-rc))
				rfin = r1 ;
			else
				rfin = r2 ; 
			
			if((cc-c1) < (c2-cc))
				cfin = c1 ;
			else
				cfin = c2 ;
				
			scale->at<Vec3b>(i, j) = image->at<Vec3b>(rfin, cfin) ;  
		}
	}    
	}
} 


int _avg_filter(int win)
{
	int orows, ocols ;

	int i, j, u, d, m, n ; 
	
	double ba, ga, ra ; 
	
	u = (win - 1) / 2;
	d = win - 1 - u ;
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	int cco = 0 ;
	
	filter = new Mat(image->rows, image->cols, CV_8UC3, CV_RGB(0,0,0)); 
	
	for(i=0;i<image->rows;i++)
	{
		for(j=0;j<image->cols;j++)
		{
		
			ba = ga = ra = 0 ;
			
			cco = 0 ;
			
			for(m=i;m>=i-u;m--)
			{
				for(n=j;n>=j-u;n--)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						ba+= intensity[0] ;
						ga+= intensity[1] ;
						ra+= intensity[2] ;
						cco++ ;
					}
				}
				for(n=j+1;n<=j+d;n++)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						ba+= intensity[0] ;
						ga+= intensity[1] ;
						ra+= intensity[2] ;
						cco++ ;
					}
				}
			}
			
			for(m=i+1;m<=i+d;m++)
			{
				for(n=j;n>=j-u;n--)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						ba+= intensity[0] ;
						ga+= intensity[1] ;
						ra+= intensity[2] ;
						cco++ ;
					}
				}
				for(n=j+1;n<=j+d;n++)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						ba+= intensity[0] ;
						ga+= intensity[1] ;
						ra+= intensity[2] ;
						cco++ ;
					}
				}
			}
			
			Vec3b ifin ;
			
			ba/= cco ;
			ga/= cco ;
			ra/= cco ;
			
			ifin[0] = ba ;
			ifin[1] = ga ;
			ifin[2] = ra ;
			
			filter->at<Vec3b>(i, j) = ifin ; 
		}
	}
}


int _median_filter(int win)
{
	int orows, ocols ;

	int i, j, u, d, m, n ; 
	
	double med ; 
	
	u = (win - 1) / 2;
	d = win - 1 - u ;
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	int rarr[win*win], garr[win*win], barr[win*win], cr, cb, cg ; 
	
	filter = new Mat(image->rows, image->cols, CV_8UC3, CV_RGB(0,0,0)); 
	
	for(i=0;i<image->rows;i++)
	{
		for(j=0;j<image->cols;j++)
		{
		
			cb = cg = cr = 0 ;
			
			for(m=i;m>=i-u;m--)
			{
				for(n=j;n>=j-u;n--)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						barr[cb++] = intensity[0] ;
						garr[cg++] = intensity[1] ;
						rarr[cr++] = intensity[2] ;
					}
				}
				for(n=j+1;n<=j+d;n++)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						barr[cb++] = intensity[0] ;
						garr[cg++] = intensity[1] ;
						rarr[cr++] = intensity[2] ;
					}
				}
			}
			
			for(m=i+1;m<=i+d;m++)
			{
				for(n=j;n>=j-u;n--)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						barr[cb++] = intensity[0] ;
						garr[cg++] = intensity[1] ;
						rarr[cr++] = intensity[2] ;
					}
				}
				for(n=j+1;n<=j+d;n++)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						barr[cb++] = intensity[0] ;
						garr[cg++] = intensity[1] ;
						rarr[cr++] = intensity[2] ;
					}
				}
			}
			
			
			sort(barr, barr+cb) ;
			sort(garr, garr+cg) ;
			sort(rarr, rarr+cr) ;
			
			
			Vec3b ifin ;
			
			
			if(cb%2 == 0)
			{
				med = (barr[cb/2] + barr[(cb+2)/2]) / 2 ; 
			}
			else
			{
				med = barr[cb/2] ; 
			}
			
			ifin[0] = med ; 
			
			
			if(cg%2 == 0)
			{
				med = (garr[cg/2] + garr[(cg+2)/2]) / 2 ; 
			}
			else
			{
				med = garr[cg/2] ; 
			}
			
			ifin[1] = med ; 
			
			
			if(cr%2 == 0)
			{
				med = (rarr[cr/2] + rarr[(cr+2)/2]) / 2 ; 
			}
			else
			{
				med = rarr[cr/2] ; 
			}
			
			ifin[2] = med ; 
			
			
			filter->at<Vec3b>(i, j) = ifin ; 
			
		}
	}
}



int _adapt_hist_eq(int win)
{
	int orows, ocols ;

	int i, j, u, d, m, n ,iter ; 
	
	u = (win - 1) / 2;
	d = win - 1 - u ;
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	int total ; //= orows * ocols ;
	
	int cco = 0 ;
	
	ad_hist_eq = new Mat(image->rows, image->cols, CV_8UC3, CV_RGB(0,0,0)); 
	
	for(i=0;i<image->rows;i++)
	{
		for(j=0;j<image->cols;j++)
		{
		
			cco = 0 ;
			
			for(iter=0;iter<256;iter++)
			{
				b_hist[iter] = 0 ;
				g_hist[iter] = 0 ;
				r_hist[iter] = 0 ;
			}
			
			for(m=i;m>=i-u;m--)
			{
				for(n=j;n>=j-u;n--)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						b_hist[intensity[0]]++ ;
						g_hist[intensity[1]]++ ;
						r_hist[intensity[2]]++ ;	
						cco++ ;
					}
				}
				for(n=j+1;n<=j+d;n++)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						b_hist[intensity[0]]++ ;
						g_hist[intensity[1]]++ ;
						r_hist[intensity[2]]++ ;	
						cco++ ;
					}
				}
			}
			
			for(m=i+1;m<=i+d;m++)
			{
				for(n=j;n>=j-u;n--)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						b_hist[intensity[0]]++ ;
						g_hist[intensity[1]]++ ;
						r_hist[intensity[2]]++ ;	
						cco++ ;
					}
				}
				for(n=j+1;n<=j+d;n++)
				{
					if(m>=0 && m<orows && n>=0 && n<ocols)
					{
						Vec3b intensity = image->at<Vec3b>(m, n);
						b_hist[intensity[0]]++ ;
						g_hist[intensity[1]]++ ;
						r_hist[intensity[2]]++ ;	
						cco++ ;  
					}
				}
			}
			
			total = cco ;
			
			for(iter=1;iter<256;iter++)
			{
				b_hist[iter]+= b_hist[iter-1] ;
				g_hist[iter]+= g_hist[iter-1] ;
				r_hist[iter]+= r_hist[iter-1] ;
			}
			
			for(iter=0;iter<256;iter++)
			{
				b_hist[iter] = (b_hist[iter] * 255) / total ;
				g_hist[iter] = (g_hist[iter] * 255) / total ;
				r_hist[iter] = (r_hist[iter] * 255) / total ;
			}
			
			Vec3b ifin = image->at<Vec3b>(i, j);
			ifin.val[0] = b_hist[ifin.val[0]];
			ifin.val[1] = g_hist[ifin.val[1]];
			ifin.val[2] = r_hist[ifin.val[2]];
			
			ad_hist_eq->at<Vec3b>(i, j) = ifin ; 
		}
	}
}




int _rotate(double deg, int interpolate)    // 0 - nearest neighbour 1 - bilinear
{
	deg = deg * 0.01745329251 ;

	int orows, ocols ;

	int i, j ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	double diag = orows * orows + ocols * ocols ;
	
	diag = sqrt(diag) ; 
	
	int size =  ceil(diag) ;
		
	rota = new Mat(size, size, CV_8UC3, CV_RGB(0,0,0)); 
	
	int ncr = size / 2, ncc = size / 2 ;
	int ocr = orows / 2, occ = ocols / 2 ; 
	
	double tcr, tcc, tmpr, tmpc ;
	
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			tmpr = ncr - i ;  // y
			tmpc = j - ncc ;  // x 
			
			tcr = -sin(deg) * tmpc + cos(deg) * tmpr ; 
			tcc = cos(deg) * tmpc + sin(deg) * tmpr ;
						
			tcc = tcc + occ ;
			tcr = ocr - tcr ;
			
			if(interpolate == 0)
			{
				if((tcc-floor(tcc)) < (ceil(tcc)-tcc))
					tcc = floor(tcc) ;
				else
					tcc = ceil(tcc) ; 
			
				if((tcr-floor(tcr)) < (ceil(tcr)-tcr))
					tcr = floor(tcr) ;
				else
					tcr = ceil(tcr) ;
			
				if(tcr>=0 && tcr<orows && tcc>=0 && tcc<ocols)
					rota->at<Vec3b>(i, j) = image->at<Vec3b>((int)tcr, (int)tcc); 
			}
			else
			{
				if(tcr>=0 && tcr<orows && tcc>=0 && tcc<ocols)
					rota->at<Vec3b>(i, j) = _bilinear(tcr, tcc) ; 
			}
		}
	}		
}


int _translate(double rv, double cv, int interpolate)     // 0 - nearest neighbour 1 - bilinear
{
	int orows, ocols ;

	int i, j ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	trans = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0)); 
	
	double tcr, tcc ;
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			tcr = i - rv ;  // y
			tcc = j - cv ;  // x 
			
			if(interpolate == 0)
			{
				if((tcc-floor(tcc)) < (ceil(tcc)-tcc))
					tcc = floor(tcc) ;
				else
					tcc = ceil(tcc) ; 
			
				if((tcr-floor(tcr)) < (ceil(tcr)-tcr))
					tcr = floor(tcr) ;
				else
					tcr = ceil(tcr) ;
			
				if(tcr>=0 && tcr<orows && tcc>=0 && tcc<ocols)
					trans->at<Vec3b>(i, j) = image->at<Vec3b>((int)tcr, (int)tcc); 
			}
			else
			{
				if(tcr>=0 && tcr<orows && tcc>=0 && tcc<ocols)
					trans->at<Vec3b>(i, j) = _bilinear(tcr, tcc) ; 
			}
		}
	}		
}



int _shear(double rv, double cv, int interpolate)     // 0 - nearest neighbour 1 - bilinear
{
	int orows, ocols ;

	int i, j ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	int size = 500 ;				// Please add variable sizes man !!!
		
	shr = new Mat(size, size, CV_8UC3, CV_RGB(0,0,0)); 
	
	int ncr = size / 2, ncc = size / 2 ;
	int ocr = orows / 2, occ = ocols / 2 ; 
	
	double tcr, tcc, tmpr, tmpc ;
	
	if(rv*cv == 1)
		return 0 ;
	
	double denom = 1 - rv * cv ; 
	
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			tmpr = ncr - i ;  // y
			tmpc = j - ncc ;  // x 
			
			tcr = (tmpr - rv*tmpc) / denom ; 
			tcc = (tmpc - cv*tmpr) / denom ; 
						
			tcc = tcc + occ ;
			tcr = ocr - tcr ;
			
			if(interpolate == 0)
			{
				if((tcc-floor(tcc)) < (ceil(tcc)-tcc))
					tcc = floor(tcc) ;
				else
					tcc = ceil(tcc) ; 
			
				if((tcr-floor(tcr)) < (ceil(tcr)-tcr))
					tcr = floor(tcr) ;
				else
					tcr = ceil(tcr) ;
			
				if(tcr>=0 && tcr<orows && tcc>=0 && tcc<ocols)
					shr->at<Vec3b>(i, j) = image->at<Vec3b>((int)tcr, (int)tcc); 
			}
			else
			{
				if(tcr>=0 && tcr<orows && tcc>=0 && tcc<ocols)
					shr->at<Vec3b>(i, j) = _bilinear(tcr, tcc) ; 
			}
		}
	}		
}


int _bit_plane(int mask)
{
	int orows, ocols ;

	int i, j ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	bitplane = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			Vec3b intensity = image->at<Vec3b>(i, j);
			intensity[0] = intensity[0] & mask ;
			intensity[1] = intensity[1] & mask ;
			intensity[2] = intensity[2] & mask ;
			bitplane->at<Vec3b>(i, j) = intensity; 
		}
	}
}



int _hist_match()
{
	hist_match = new Mat(image->rows, image->cols, CV_8UC3, CV_RGB(0,0,0));

	_histogram(image) ;
	
	int newb[256], newg[256], newr[256] ;
	
	memset (newb, 0, sizeof(newb)) ;
	memset (newg, 0, sizeof(newg)) ;
	memset (newr, 0, sizeof(newr)) ;
	
	int i, j ; 
	
	for(i=0;i<sample->rows;i++)
    {
    	for(j=0;j<sample->cols;j++)
    	{
    		Vec3b intensity = sample->at<Vec3b>(i, j);
			newb[intensity.val[0]]++ ;
			newg[intensity.val[1]]++ ;
			newr[intensity.val[2]]++ ;
     	}
    }
    
    
    for(i=1;i<256;i++)
    {
    	b_hist[i]+= b_hist[i-1] ;
    	g_hist[i]+= g_hist[i-1] ;
    	r_hist[i]+= r_hist[i-1] ;
    	
    	newb[i]+= newb[i-1] ;
    	newg[i]+= newg[i-1] ;
    	newr[i]+= newr[i-1] ;
    }

	int kb, ob, kg, og, kr, ora, temp ;
	
	int mapb[256], mapg[256], mapr[256] ; 

	for(i=0;i<256;i++)
	{
		kb = kg = kr = 2147483647 ;
		ob = og = ora = -1 ;
		
		for(j=0;j<256;j++)
		{
			temp = abs(b_hist[i] - newb[j]) ; 
			if(temp < kb)
			{
				kb = temp ;
				ob = j ; 
			}
			
			temp = abs(g_hist[i] - newg[j]) ; 
			if(temp < kg)
			{
				kg = temp ;
				og = j ; 
			}
			
			temp = abs(r_hist[i] - newr[j]) ; 
			if(temp < kr)
			{
				kr = temp ;
				ora = j ; 
			}
		}
		
		mapb[i] = ob ;
		mapg[i] = og ; 
		mapr[i] = ora ; 
	}

	for(i=0;i<image->rows;i++)
    {
    	for(j=0;j<image->cols;j++)
    	{
    		Vec3b intensity = image->at<Vec3b>(i, j);
			intensity[0] = mapb[intensity[0]] ;
			intensity[1] = mapg[intensity[1]] ;
			intensity[2] = mapr[intensity[2]] ;
			hist_match->at<Vec3b>(i, j) = intensity ; 
     	}
    }
}


int _laplacian()
{
	int orows, ocols ;

	int i, j, ba, ga, ra, cco ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	spatial = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			cco = 1 ;
			ba = ga = ra = 0 ; 
			
			if(i-1>=0)
			{
				cco++ ; 
				Vec3b intensity = image->at<Vec3b>(i-1, j);
				ba-= intensity[0] ;
				ga-= intensity[1] ;
				ra-= intensity[2] ;	
			}
			
			if(i+1<orows)
			{
				cco++ ; 
				Vec3b intensity = image->at<Vec3b>(i+1, j);
				ba-= intensity[0] ;
				ga-= intensity[1] ;
				ra-= intensity[2] ;	
			}
			
			if(j-1>=0)
			{
				cco++ ; 
				Vec3b intensity = image->at<Vec3b>(i, j-1);
				ba-= intensity[0] ;
				ga-= intensity[1] ;
				ra-= intensity[2] ;	
			}
			
			if(j+1<ocols)
			{
				cco++ ; 
				Vec3b intensity = image->at<Vec3b>(i, j+1);
				ba-= intensity[0] ;
				ga-= intensity[1] ;
				ra-= intensity[2] ;	
			}
			
			Vec3b intensity = image->at<Vec3b>(i, j);
			ba+= cco*intensity[0] ;
			ga+= cco*intensity[1] ;
			ra+= cco*intensity[2] ;  
			
			
			
			
			if(ba < 0)
				ba = 0 ;
			else if(ba > 255)
				ba = 255 ;
			
			if(ga < 0)
				ga = 0 ;
			else if(ga > 255)
				ga = 255 ;
			
			if(ra < 0)
				ra = 0 ;
			else if(ra > 255)
				ra = 255 ;
			
			intensity[0] = ba ;
			intensity[1] = ga ;
			intensity[2] = ra ;
				
			
			spatial->at<Vec3b>(i, j) = intensity ;  
		}
	}
}


int _highboost(double factor) 
{
	int orows, ocols ;

	int i, j  ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	spatial = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));
	
	_avg_filter(5) ;
	
	double value ; 
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			Vec3b orig = image->at<Vec3b>(i, j);
			Vec3b mean = filter->at<Vec3b>(i, j);
			Vec3b ifin ;
			
			value = factor * (orig[0] - mean[0]) + orig[0] ;
			ifin[0] = value ; 
			
			value = factor * (orig[1] - mean[1]) + orig[1] ;
			ifin[1] = value ;
			
			value = factor * (orig[2] - mean[2]) + orig[2] ;
			ifin[2] = value ;
			
			spatial->at<Vec3b>(i, j) = ifin ;  
		}
	}
}

int _threshold(double thres)
{
	int orows, ocols ;

	int i, j  ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	piece = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));
	
	double value ; 
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			Vec3b intensity = image->at<Vec3b>(i, j) ; 
			
			value = (double) intensity[0] / 255 ;
			
			if(value >= thres)
			{
				intensity[0] = 255 ; 
			}
			else
			{
				intensity[0] = 0 ;
			}
			
			value = (double) intensity[1] / 255 ;
			
			if(value >= thres)
			{
				intensity[1] = 255 ; 
			}
			else
			{
				intensity[1] = 0 ;
			}
			
			value = (double) intensity[2] / 255 ;
			
			if(value >= thres)
			{
				intensity[2] = 255 ; 
			}
			else
			{
				intensity[2] = 0 ;
			}
			
			piece->at<Vec3b>(i, j) = intensity ;  
		}
	}
}


int _contrast(double x1, double y1, double x2, double y2)
{
	int orows, ocols ;

	int i, j  ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	piece = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));

	double value ; 
	
	double s1 = y1 / x1 ;
	double s2 = (y2 - y1) / (x2 - x1) ; 
	double s3 = (255 - y2) / (255 - x2) ; 
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			Vec3b inte = image->at<Vec3b>(i, j) ; 
				
			if(inte[0] <= x1)
			{
				value = (double) (inte[0] * s1) ; 
			}
			else if(inte[0] > x1 && inte[0] <= x2)
			{
				value = (double) (inte[0] - x1) * s2 + y1 ; 
			} 
			else
			{
				value = (double) (inte[0] - x2) * s3 + y2 ;  
			}
			
			inte[0] = value ; 
			
			if(inte[1] <= x1)
			{
				value = (double) (inte[1] * s1) ; 
			}
			else if(inte[1] > x1 && inte[1] <= x2)
			{
				value = (double) (inte[1] - x1) * s2 + y1 ; 
			} 
			else
			{
				value = (double) (inte[1] - x2) * s3 + y2 ;  
			}
			
			inte[1] = value ; 
			
			if(inte[2] <= x1)
			{
				value = (double) (inte[2] * s1) ; 
			}
			else if(inte[2] > x1 && inte[2] <= x2)
			{
				value = (double) (inte[2] - x1) * s2 + y1 ; 
			} 
			else
			{
				value = (double) (inte[2] - x2) * s3 + y2 ;  
			}
			
			inte[2] = value ; 
			
			piece->at<Vec3b>(i, j) = inte ;  
		}
	}
}


int _power(double con, double gamma)
{
	int orows, ocols ;

	int i, j  ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	power = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));

	double value ; 
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			Vec3b inte = image->at<Vec3b>(i, j) ;
			
			value = (double) inte[0] / 255 ; 
			value = con * pow(value, gamma) ; 
			value = value * 255 ; 
			inte[0] = value ; 
			
			value = (double) inte[1] / 255 ; 
			value = con * pow(value, gamma) ; 
			value = value * 255 ; 
			inte[1] = value ; 
			
			value = (double) inte[2] / 255 ; 
			value = con * pow(value, gamma) ; 
			value = value * 255 ; 
			inte[2] = value ; 
			
			power->at<Vec3b>(i, j) = inte ;  
		}
	}
}



int _logarithm(double con)
{
	int orows, ocols ;

	int i, j  ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	loga = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));

	double value ; 
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			Vec3b inte = image->at<Vec3b>(i, j) ;
			
			value = (double) inte[0] / 255 ; 
			value = con * log(1 + value) ; 
			value = value * 255 ; 
			inte[0] = value ; 
			
			value = (double) inte[1] / 255 ; 
			value = con * log(1 + value) ; 
			value = value * 255 ; 
			inte[1] = value ; 
			
			value = (double) inte[2] / 255 ; 
			value = con * log(1 + value) ;  
			value = value * 255 ; 
			inte[2] = value ; 
			
			loga->at<Vec3b>(i, j) = inte ;  
		}
	}
}


int _inverse_logarithm(double con)
{
	int orows, ocols ;

	int i, j  ; 
	
	orows = image->rows ;
	ocols = image->cols ;  
	
	cout << orows << " " << ocols << endl ;;
	
	loga = new Mat(orows, ocols, CV_8UC3, CV_RGB(0,0,0));

	double value ; 
	
	for(i=0;i<orows;i++)
	{
		for(j=0;j<ocols;j++)
		{
			Vec3b inte = image->at<Vec3b>(i, j) ;
			
			value = (double) inte[0] / 255 ; 
			value = value / con ; 
			value = pow(2.71828, value) - 1 ;  
			value = value * 255 ; 
			inte[0] = value ; 
			
			value = (double) inte[1] / 255 ; 
			value = value / con ; 
			value = pow(2.71828, value) - 1 ; 
			value = value * 255 ; 
			inte[1] = value ; 
			
			value = (double) inte[2] / 255 ; 
			value = value / con ; 
			value = pow(2.71828, value) - 1 ; 
			value = value * 255 ; 
			inte[2] = value ; 
			
			loga->at<Vec3b>(i, j) = inte ;  
		}
	}
}


int main(int argc, char** argv)
{
    if( argc < 2)
    {
    	cout << "Error : Not enough arguments" << endl ;
    	return -1;
    }


    Mat input = imread(argv[1], CV_LOAD_IMAGE_COLOR), sam ;
    
    image = &input ; 
    
    if(argc > 2)
    {
    	sam = imread(argv[2], CV_LOAD_IMAGE_COLOR) ; 
    	sample = &sam ; 
    	
    }
    
    

    
    //cout << image->rows << " " << image->cols << endl ;

    if(!image->data)
    {
        cout << "Error : Could not open or find the image" << endl ;
        return -1;
    }
    
    // Code begins
    
    //   Mat M(100,100, CV_8UC3, Scalar(0,0,255));
    
    cout << "Menu : " << endl << endl ;
    cout << "1 - Translation" << endl ;
    cout << "2 - Rotation" << endl ;
    cout << "3 - Scaling" << endl ;
    cout << "4 - Shear" << endl ;
    cout << "5 - Negative" << endl ;
    cout << "6 - Log Transform" << endl ;
    cout << "7 - Power Law Tranform" << endl ;
    cout << "8 - Thresholding" << endl ;
    cout << "9 - Bit Plane Slicing" << endl ;
    cout << "10 - Contrast Stretching" << endl ;
    cout << "11 - Histogram Equalization" << endl ;
    cout << "12 - Adaptive Histogram Equalization" << endl ;
    cout << "13 - Histogram Matching" << endl ;
    cout << "14 - Mean Filter" << endl ;
    cout << "15 - Median Filter" << endl ;
    cout << "16 - Laplacian Filter" << endl ;
    cout << "17 - Unsharp Filter" << endl ;
    cout << "18 - High-boost Filter" << endl ;
    
    int ch ;
    
    cout << "\nEnter your choice : " ;
    
    cin >> ch ; 
    
	int in, win ;
	double rv, cv ;
	
	Mat *result  ;
    
    switch(ch)
    {
    	case 1: 
    	cout << "Enter Delta X : " ;
    	cin >> cv ;
    	cout << "Enter Delta Y : " ;
    	cin >> rv ;
    	cout << "Enter 0 for nearest neighbour and 1 for bilinear : " ;
    	cin >> in ;
    	rv = -rv ; 
    	_translate(rv, cv, in) ;
    	result = trans ;
    	break ;
    	
    	case 2:
    	double angle ;
    	cout << "Enter angle : " ;
    	cin >> angle ;
    	cout << "Enter 0 for nearest neighbour and 1 for bilinear : " ;
    	cin >> in ;
    	_rotate(angle, in) ;
    	result = rota ;
    	break ;
    	
    	case 3:
    	cout << "Enter Sx : " ;
    	cin >> cv ;
    	cout << "Enter Sy : " ;
    	cin >> rv ;
    	cout << "Enter 0 for nearest neighbour and 1 for bilinear : " ;
    	cin >> in ;
		downsize(rv, cv, in) ;
		result = scale ;
    	break ;
    	
    	case 4:
    	cout << "Enter X shear factor : " ;
    	cin >> cv ;
    	cout << "Enter Y shear factor : " ;
    	cin >> rv ;
    	cout << "Enter 0 for nearest neighbour and 1 for bilinear : " ;
    	cin >> in ;
		_shear(rv, cv, in) ;
		result = shr ;
    	break ;
    	
    	case 5:
		_negative() ;
		result = negative ;
    	break ;
    	
    	case 6:
    	cout << "Enter value of constant : " ;
    	cin >> cv ;
		_logarithm(cv) ;
		result = loga ;
    	break ;
    	
    	case 7:
    	cout << "Enter value of constant : " ;
    	cin >> cv ;
    	cout << "Enter value of gamma : " ;
    	cin >> rv ;
		_power(cv, rv) ;
		result = power ;
    	break ;
    	
    	case 8:
    	cout << "Enter value of threshold : " ;
    	cin >> cv ;
		_threshold(cv) ;
		result = piece ;
    	break ;
    	
    	case 9:
    	int mask ;
    	cout << "Enter value of mask : " ;
    	cin >> mask ;
		_bit_plane(mask) ;
		result = bitplane ;
    	break ;
    	
    	case 10:
    	double x1, y1, x2, y2 ;
    	cout << "Enter value of x1 : " ;
    	cin >> x1 ;
    	cout << "Enter value of y1 : " ;
    	cin >> y1 ;
    	cout << "Enter value of x2 : " ;
    	cin >> x2 ;
    	cout << "Enter value of y2 : " ;
    	cin >> y2 ;
		_contrast(x1, y1, x2, y2) ;
		result = piece ;
    	break ;
    	
    	case 11:
    	_hist_eq() ;
		result = hist_eq ;
    	break ;
    	
    	case 12:
    	cout << "Enter window size : " ;
    	cin >> win ;
    	_adapt_hist_eq(win) ;
		result = ad_hist_eq ;
    	break ;
    	
    	case 13:
    	_hist_match() ;
		result = hist_match ;
    	break ;
    	
    	case 14:
    	cout << "Enter window size : " ;
    	cin >> win ;
    	_avg_filter(win) ;
		result = filter ;
    	break ;
    	
    	case 15:
    	cout << "Enter window size : " ;
    	cin >> win ;
    	_median_filter(win) ;
		result = filter ;
    	break ;
    	
    	case 16:
    	_laplacian() ;
		result = spatial ;
    	break ;
    	
    	case 17:
    	_highboost(1) ;
		result = spatial ;
    	break ;
    	
    	case 18:
    	cout << "Enter the value of k : " ;
    	cin >> cv ;
    	_highboost(cv) ;
		result = spatial ;
    	break ;
    }
    
    //_negative() ;
    //_hist_eq() ;
    //downsize(1.5,1.5,1) ; 
    //_median_filter(5) ;
    //_avg_filter(10) ;
   	//_rotate(-23, 0) ;
 
 	//_translate(+150.24, -160.67, 0) ;
 	
	//_shear(0, 1, 0) ;
	//_bit_plane(96) ;  
	//_hist_match() ;
	
	//_adapt_hist_eq(100) ;  // Do something here
    
    
    //_hist_match() ;
    
    //_laplacian() ;
    //_highboost(3) ;
    
	//_threshold(0.4) ;
	//_contrast(40,50,100,230) ;
    
    //_power(1, 0.4) ;
    //_logarithm(1.45) ; 
	//_inverse_logarithm(4) ; 
    
      
    imwrite("result.png", *result);
    
    namedWindow( "Display window", WINDOW_AUTOSIZE );
    imshow( "Display window", *result ); 
    
          
    
    waitKey(0);                                      
    
    return 0;
}


// Adaptive hist needs check 
// Variable size in shear 
//


