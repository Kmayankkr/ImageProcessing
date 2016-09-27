
#include <stdio.h>
#include <queue> 
#include <iostream> 
#include <time.h> 
#include <map>
#include <list>
#include <cstdlib>
#include <cstring>
#include <algorithm> 
//#include <opencv2/opencv.hpp>

//using namespace cv;
using namespace std;  

int arr[1000][1000], visit[1000][1000], labelsize[1000001], seg, r, c, nbr, wrap, choice ;

int dfsdepth, bfsdepth, maxdfsdepth, maxbfsdepth ; 

double prob ; 

int seqgrp[1000001] ;
list<int> grp[1000001] ;


int dfs4(int x, int y)
{

	dfsdepth++ ; 
	
	if(maxdfsdepth < dfsdepth)
		maxdfsdepth = dfsdepth ; 

	visit[x][y] = seg ; 
	labelsize[seg]++ ; 

	if(x>0 && arr[x-1][y]==arr[x][y] && visit[x-1][y]==0)
	{
		dfs4(x-1,y) ; 
	}
	
	if(x<r-1 && arr[x+1][y]==arr[x][y] && visit[x+1][y]==0)
	{
		dfs4(x+1,y) ; 
	}
	
	if(y>0 && arr[x][y-1]==arr[x][y] && visit[x][y-1]==0)
	{
		dfs4(x,y-1) ; 
	}

	if(y<c-1 && arr[x][y+1]==arr[x][y] && visit[x][y+1]==0)
	{
		dfs4(x,y+1) ; 
	}
	
	if(nbr == 8)
	{
	
		if(x>0 && y>0 && arr[x-1][y-1]==arr[x][y] && visit[x-1][y-1]==0)
		{
			dfs4(x-1,y-1) ; 
		}
	
		if(x<r-1 && y>0 && arr[x+1][y-1]==arr[x][y] && visit[x+1][y-1]==0)
		{
			dfs4(x+1,y-1) ; 
		}
	
		if(x>0 && y<c-1 && arr[x-1][y+1]==arr[x][y] && visit[x-1][y+1]==0)
		{
			dfs4(x-1,y+1) ; 
		}

		if(x<r-1 && y<c-1 && arr[x+1][y+1]==arr[x][y] && visit[x+1][y+1]==0)
		{
			dfs4(x+1,y+1) ;
		}
	
	}
	
	dfsdepth-- ; 
	
}


int wdfs4(int x, int y)
{

	dfsdepth++ ; 
	
	if(maxdfsdepth < dfsdepth)
		maxdfsdepth = dfsdepth ; 

	visit[x][y] = seg ; 
	labelsize[seg]++ ; 

	int nx, ny ;
	
	nx = (x-1+r)%r ;
	ny = y ;
	if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
	{
		wdfs4(nx,ny) ; 
	}
	
	nx = (x+1)%r ;
	ny = y ;
	if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
	{
		wdfs4(nx,ny) ; 
	}
	
	nx = x ;
	ny = (y-1+c)%c ;
	if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
	{
		wdfs4(nx,ny) ; 
	}

	nx = x ;
	ny = (y+1)%c ;
	if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
	{
		wdfs4(nx,ny) ; 
	}
	
	if(nbr == 8)
	{
	
		nx = (x-1+r)%r ;
		ny = (y-1+c)%c ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			wdfs4(nx,ny) ; 
		}
	
		nx = (x+1)%r ;
		ny = (y-1+c)%c ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			wdfs4(nx,ny) ; 
		}
	
		nx = (x-1+r)%r ;
		ny = (y+1)%c ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			wdfs4(nx,ny) ; 
		}

		nx = (x+1)%r ;
		ny = (y+1)%c ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			wdfs4(nx,ny) ; 
		}
	
	}
	
	dfsdepth-- ; 
	
}


int bfs4(queue<int> qx, queue<int> qy)
{

	if(qx.size() > maxbfsdepth)
		maxbfsdepth = qx.size() ;

	int x, y ;

	queue<int> tx, ty ;  

	while(!qx.empty())
	{
	
		x = qx.front() ;
		y = qy.front() ; 
		
		labelsize[visit[x][y]]++ ;
		
		qx.pop() ;
		qy.pop() ; 
		

		if(x>0 && arr[x-1][y]==arr[x][y] && visit[x-1][y]==0)
		{
			tx.push(x-1) ;
			ty.push(y) ; 
			visit[x-1][y] = seg ;
		}
	
		if(x<r-1 && arr[x+1][y]==arr[x][y] && visit[x+1][y]==0)
		{
			tx.push(x+1) ;
			ty.push(y) ; 
			visit[x+1][y] = seg ; 
		}
	
		if(y>0 && arr[x][y-1]==arr[x][y] && visit[x][y-1]==0)
		{
			tx.push(x) ;
			ty.push(y-1) ; 
			visit[x][y-1] = seg ;
		}

		if(y<c-1 && arr[x][y+1]==arr[x][y] && visit[x][y+1]==0)
		{
			tx.push(x) ;
			ty.push(y+1) ;
			visit[x][y+1] = seg ;  
		}
		
		if(nbr == 8)
		{
			if(x>0 && y>0 && arr[x-1][y-1]==arr[x][y] && visit[x-1][y-1]==0)
			{
				tx.push(x-1) ;
				ty.push(y-1) ;
				visit[x-1][y-1] = seg ;  
			}
	
			if(x<r-1 && y>0 && arr[x+1][y-1]==arr[x][y] && visit[x+1][y-1]==0)
			{
				tx.push(x+1) ;
				ty.push(y-1) ;
				visit[x+1][y-1] = seg ;  
			}
	
			if(x>0 && y<c-1 && arr[x-1][y+1]==arr[x][y] && visit[x-1][y+1]==0)
			{
				tx.push(x-1) ;
				ty.push(y+1) ;
				visit[x-1][y+1] = seg ;   
			}

			if(x<r-1 && y<c-1 && arr[x+1][y+1]==arr[x][y] && visit[x+1][y+1]==0)
			{
				tx.push(x+1) ;
				ty.push(y+1) ;
				visit[x+1][y+1] = seg ;  
			}
		}
		
	}
	
	if(tx.size())
		bfs4(tx, ty) ;  

}



int wbfs4(queue<int> qx, queue<int> qy)
{

	if(qx.size() > maxbfsdepth)
		maxbfsdepth = qx.size() ; 

	int x, y, nx, ny ;

	queue<int> tx, ty ;  

	while(!qx.empty())
	{
	
		x = qx.front() ;
		y = qy.front() ; 
		
		labelsize[visit[x][y]]++ ;
		
		qx.pop() ;
		qy.pop() ; 
		
		nx = (x-1+r)%r ; 
		ny = y ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			tx.push(nx) ;
			ty.push(ny) ; 
			visit[nx][ny] = seg ;
		}
	
		nx = (x+1)%r ; 
		ny = y ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			tx.push(nx) ;
			ty.push(ny) ; 
			visit[nx][ny] = seg ;
		}
	
		nx = x ; 
		ny = (y-1+c)%c ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			tx.push(nx) ;
			ty.push(ny) ; 
			visit[nx][ny] = seg ;
		}

		nx = x ; 
		ny = (y+1)%c ;
		if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
		{
			tx.push(nx) ;
			ty.push(ny) ; 
			visit[nx][ny] = seg ;
		}
		
		if(nbr == 8)
		{
			nx = (x-1+r)%r ; 
			ny = (y-1+c)%c ;
			if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
			{
				tx.push(nx) ;
				ty.push(ny) ; 
				visit[nx][ny] = seg ;
			}
	
			nx = (x+1)%r ; 
			ny = (y-1+c)%c ;
			if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
			{
				tx.push(nx) ;
				ty.push(ny) ; 
				visit[nx][ny] = seg ;
			}
	
			nx = (x-1+r)%r ; 
			ny = (y+1)%c ;
			if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
			{
				tx.push(nx) ;
				ty.push(ny) ; 
				visit[nx][ny] = seg ;
			}

			nx = (x+1)%r ; 
			ny = (y+1)%c ;
			if(arr[nx][ny]==arr[x][y] && visit[nx][ny]==0)
			{
				tx.push(nx) ;
				ty.push(ny) ; 
				visit[nx][ny] = seg ;
			}
		}
		
	}
	
	
	
	if(tx.size())
		wbfs4(tx, ty) ;  

}




int seq4()
{

	int i, j ;
	
	seg = 0 ; 

	for(i=0;i<r;i++)
	{	
		for(j=0;j<c;j++)
		{
			if(arr[i][j]==0)
				continue ;
		
			if(j==0 || arr[i][j-1]==0)
			{
				seg++ ;
				visit[i][j] = seg ;
				labelsize[seg]++ ;
			}
			
			if(j<c-1 && arr[i][j+1]==arr[i][j])
			{
				if(visit[i][j+1]==0)
				{
					visit[i][j+1] = seg ;
					labelsize[seg]++ ;
				}
				else
				{
					grp[visit[i][j+1]].push_back(visit[i][j]) ;
				    grp[visit[i][j]].push_back(visit[i][j+1]) ;
				}
			}
			
			if(i>0)
			{
				if(arr[i][j]==arr[i-1][j])
				{
					grp[visit[i][j]].push_back(visit[i-1][j]) ;
				    grp[visit[i-1][j]].push_back(visit[i][j]) ;
				} 
				
				if(nbr == 8)
				{
					if(j<c-1 && arr[i][j]==arr[i-1][j+1])
					{
						grp[visit[i][j]].push_back(visit[i-1][j+1]) ;
						grp[visit[i-1][j+1]].push_back(visit[i][j]) ;	
					}
					
					if(j>0 && arr[i][j]==arr[i-1][j-1])
					{	
						grp[visit[i][j]].push_back(visit[i-1][j-1]) ;
						grp[visit[i-1][j-1]].push_back(visit[i][j]) ;
					}
				}
			}
		}	
	}
}


int wseq4()
{

	int i, j, ni, nj ;
	
	seg = 0 ; 

	for(i=0;i<r;i++)
	{	
		for(j=0;j<c;j++)
		{
		
			if(arr[i][j]==0)
				continue ;
		
			if(j==0 || arr[i][j-1]==0)
			{
				seg++ ;
				visit[i][j] = seg ;
				labelsize[seg]++ ;
			}
			
			ni = i ;
			nj = (j+1)%c ;
			if(arr[ni][nj]==arr[i][j])
			{
				if(visit[ni][nj]==0)
				{
					visit[ni][nj] = seg ;
					labelsize[seg]++ ;
				}
				else
				{
					grp[visit[ni][nj]].push_back(visit[i][j]) ;
				    grp[visit[i][j]].push_back(visit[ni][nj]) ;
				}
			}
			
			
			
			ni = (i-1+r)%r ;
			nj = j ;
			if(arr[i][j]==arr[ni][nj] && visit[ni][nj]>0)
			{
				grp[visit[i][j]].push_back(visit[ni][nj]) ;
			    grp[visit[ni][nj]].push_back(visit[i][j]) ;
			} 
			
			if(nbr == 8)
			{
				ni = (i-1+r)%r ;
				nj = (j+1)%c ;
				if(arr[i][j]==arr[ni][nj] && visit[ni][nj]>0)
				{
					grp[visit[i][j]].push_back(visit[ni][nj]) ;
					grp[visit[ni][nj]].push_back(visit[i][j]) ;	
				}
				
				ni = (i-1+r)%r ;
				nj = (j-1+c)%c ; 
				if(arr[i][j]==arr[ni][nj] && visit[ni][nj]>0)
				{	
					grp[visit[i][j]].push_back(visit[ni][nj]) ;
					grp[visit[ni][nj]].push_back(visit[i][j]) ;
				}
			}
		}	
	}
	
	i = 0 ;
	
	for(j=0;j<c;j++)
	{ 
		ni = (i-1+r)%r ;
		nj = j ;
		if(arr[i][j]==arr[ni][nj] && visit[ni][nj]>0)
		{
			grp[visit[i][j]].push_back(visit[ni][nj]) ;
		    grp[visit[ni][nj]].push_back(visit[i][j]) ;
		} 
		
		if(nbr == 8)
		{
			ni = (i-1+r)%r ;
			nj = (j+1)%c ;
			if(arr[i][j]==arr[ni][nj] && visit[ni][nj]>0)
			{
				grp[visit[i][j]].push_back(visit[ni][nj]) ;
				grp[visit[ni][nj]].push_back(visit[i][j]) ;	
			}
			
			ni = (i-1+r)%r ;
			nj = (j-1+c)%c ; 
			if(arr[i][j]==arr[ni][nj] && visit[ni][nj]>0)
			{	
				grp[visit[i][j]].push_back(visit[ni][nj]) ;
				grp[visit[ni][nj]].push_back(visit[i][j]) ;
			}
		}
	}
}



int seqdfs(int v)
{
	seqgrp[v] = 1 ;
	
	int child ;
	
	while(!grp[v].empty())
	{
		child = grp[v].front() ;
		grp[v].pop_front() ;
		
		if(seqgrp[child]==0)
			labelsize[v]+= seqdfs(child) ;
	}
	
	int val = labelsize[v] ;
	labelsize[v] = 0 ;
	
	return val ;
	
}


int generate()
{

	double ran, ma = RAND_MAX ;
	int i, j ;

	srand(time(NULL)) ;
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			ran = (rand() / ma) ;
			
			if(ran <= prob)
				arr[i][j] = 1 ;
			else
				arr[i][j] = 0 ; 
		}
	}

}


int print()
{
	int i, j, max = 0 ;
	
	double sum = 0, total ; 

	freopen ("output.txt","w",stdout);

	cout << "Output Matrix : \n" ;
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			cout << arr[i][j] << " " ;
		}
		cout << endl ; 
	}

	cout << "\nTotal number of clusters : " << seg << endl << endl ;

	for(i=1;i<=seg;i++)
	{
		sum+= labelsize[i] ;
		if(labelsize[i] > max)
			max = labelsize[i] ;
	}
	
	total = seg ;
	
	double avg ;
	avg = sum / seg ;
	
	cout << "Average size of clusters : " << avg << endl << endl ;
	cout << "Size of biggest cluster : " << max << endl << endl ;
	
	if(choice == 0)
		cout << "Maximum size of queue : " << maxbfsdepth << endl << endl ;
	else if(choice == 1)
		cout << "Maximum size of stack : " << maxdfsdepth << endl << endl ;
		
	fclose(stdout) ;

}


int main(int argc, char** argv )
{
	

	int i, j ;
	
	cout << "Enter no. of rows : " ;
	cin >> r ;
	
	cout << "\nEnter no. of columns : " ;
	cin >> c ;
	
	cout << "\nWhich neighbourhood - 4 or 8 : " ;
	cin >> nbr ;  
	
	cout << "\nEnter probability of occurence of 1 : " ;
	cin >> prob ;
	
	cout << "\nDo you need wrap around ? : " ;
	cin >> wrap ;
	
	cout << "\nEnter 0 for BFS, 1 for DFS, 2 for sequential " ;
	cin >> choice ;
	
	dfsdepth = bfsdepth = 0 ;
	maxdfsdepth = maxbfsdepth = 0 ;    
	
	
	generate() ;
	
	
	if(choice == 0)
	{
	
	memset(visit, 0, sizeof(visit[0][0])*1000000) ;  
	memset(labelsize, 0, sizeof(labelsize[0])*1000001) ;  
	
	seg = 0 ;
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			if(visit[i][j]==0 && arr[i][j]==1)
			{
				seg++ ; 
				visit[i][j] = seg ;
				
				queue<int> tx, ty ;
				tx.push(i) ;
				ty.push(j) ;
				
				if(wrap==1)
					wbfs4(tx, ty) ;
				else
					bfs4(tx, ty) ;  
			}
		}
	}
    
	cout << seg << " " << maxbfsdepth << endl ;
	
	labelsize[0] = -1 ;
	sort(labelsize, labelsize+seg+1) ;
	//print(seg) ;
	
	}
	
	if(choice == 1)
	{
	memset(visit, 0, sizeof(visit[0][0])*1000000) ;  
	memset(labelsize, 0, sizeof(labelsize[0])*1000001) ;  
	
	seg = 0 ;
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			if(visit[i][j]==0 && arr[i][j]==1)
			{
				seg++ ; 
				visit[i][j] = seg ;
				
				if(wrap==1)
					wdfs4(i, j) ;
				else
					dfs4(i, j) ;   
			}
		}
	}
    
	cout << seg << " " << maxdfsdepth << endl ;
	
	
	labelsize[0] = -1 ;
	sort(labelsize, labelsize+seg+1) ;
	//print(seg) ;
	}
	
	if(choice == 2)
	{
	memset(visit, 0, sizeof(visit[0][0])*1000000) ; 
	memset(labelsize, 0, sizeof(labelsize[0])*1000001) ; 
	memset(seqgrp, 0, sizeof(seqgrp[0])*1000001) ;  
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			if(visit[i][j]==0 && arr[i][j]==1)
			{
				seg = 1 ; 
				visit[i][j] = seg ;
				
				if(wrap==1)
					wseq4() ;
				else
					seq4() ;  
			}
		}
	}
	
	//cout << seg << endl ;
	

    
    int tot = 0 ;
    
    for(i=1;i<=seg;i++)
    {
    	if(seqgrp[i]==0)
    	{
    		tot++ ; 
    		labelsize[i]+= seqdfs(i) ;
    		
    	}
    }
    
    seg = tot ;
    
	cout << tot << endl ;
	
	labelsize[0] = -1 ;
	sort(labelsize, labelsize+seg+1) ;
	//print(seg) ;
	}

	print() ;

    return 0;
    
}




