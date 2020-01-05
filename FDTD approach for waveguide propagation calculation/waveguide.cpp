#include <stdio.h>
#include <math.h>
#include <time.h>


struct point{
	double Ez;
	double Hx;
	double Hy;
}p[3000][3000];       /*use structure to store electromagnetic field, TE wave*/

int main()
{
     /*define variables*/
	
	int Imax,Imin,Jmax,Jmin;
	int Wl,Wl1;
	int i,j;
	int step,Timestop;
	int Alter;                                 /*select variables for corners*/
	int Jsource,Ismin,Ismax;
	double t;
	double n0;
	double m[2]={1.0,3.0};	
	double EAB[4][3000];                              /*stroe boundaries*/
	double EAC[4][2][2];                             /*store corners*/
	double Ei;                                        /*defince source*/
	double T1,T2;
	double Eps;
	double ia,jb;
	double wavelength;                                /*define wavelength*/
	double r;                                              //waveguide radius  
	double d;                                              //waveguide width
	double temp;
    double PI=3.1415926535895;
    double u=PI*0.0000004;                                 /*permeability*/  
    double Eps0=8.85e-12;                            /*permittivity*/
    double c=299792457.4;                               /*speed of light*/
    double coeff_E;                            /*coefficient of electric field*/
    double coeff_H=1/c/u/2;                    /*coefficient of electric field*/ 
    
	FILE *fp;
	if((fp=fopen("waveguide.txt","w"))==NULL)
		printf("Can not open the file!\n");

	/*initialize variables*/
	Wl=10;
	Wl1=2;
	Imin=0;
	Imax=2999;
	Jmin=0;
	Jmax=2999;
	Timestop=13000;
	Alter=0;
	Jsource=1;
	Ismin=95;
	Ismax=105;
	r=0.0004495;
    d=0.00000155;
    wavelength=1.55e-6;

     time_t time1;
     time_t time2;
	 time1=time(NULL);
	for(i=Imin;i<=Imax;i++)
		for(j=Jmin;j<=Jmax;j++)
		{
			p[i][j].Ez=0.0;
			p[i][j].Hx=0.0;
			p[i][j].Hy=0.0;
		}
    
	for(i=0;i<4;i++)
		for(j=0;j<=2999;j++)
			EAB[i][j]=0.0;

	for(i=0;i<4;i++)
		for(j=0;j<2;j++)
		{
			EAC[i][j][0]=0.0;
			EAC[i][j][1]=0.0;
		}
  
	/*time based loop*/
    for(step=1;step<=Timestop;step++)
	{
		/*set source*/
		Ei=sin(step*2*PI/Wl/Wl1);
	
		for(i=Ismin;i<=Ismax;i++)
		{
			p[i][Jsource].Ez=p[i][Jsource].Ez+Ei;}


		   printf("Step %d is processing!\n",step);
		

		/*calcaulate Ez using FDTD except boundaries*/
		for(i=Imin+1;i<=Imax-1;i++)
			for(j=Jmin+1;j<=Jmax-1;j++)
			{   ia=(double)(2999-i);
			    jb=(double)j;
				temp=sqrt(ia*ia+jb*jb); 
				temp=temp*wavelength/Wl;
				if((temp>=r-d/2)&&(temp<=r+d/2))
			                  n0=m[1];
			      else 
					          n0=m[0]; 
				  
				   Eps=Eps0*sqrt(n0);
                   coeff_E=1/c/Eps/2;                          
                  p[i][j].Ez=p[i][j].Ez+coeff_E*(p[i][j].Hy-p[i][j-1].Hy-p[i-1][j].Hx+p[i][j].Hx);
			}
               

			
				
    	/*second-order Mur absorption boundary*/
			for(j=Jmin+1;j<=Jmax-1;j++)
			{
				/*bottom boundary*/
            T1=p[Imax][j].Hy-p[Imax][j-1].Hy+p[Imax-1][j].Hy-p[Imax-1][j-1].Hy;
            T2=p[Imax-1][j].Ez-p[Imax][j].Ez;
			p[Imax][j].Ez=EAB[0][j]+((1-Wl1)/(1+Wl1))*T2+c*u/(2+2*Wl1)*T1;
			EAB[0][j]=p[Imax-1][j].Ez; 
			
			    /*top boundary*/
			T1=p[Imin][j].Hy-p[Imin][j-1].Hy+p[Imin+1][j].Hy-p[Imin+1][j-1].Hy;
            T2=p[Imin+1][j].Ez-p[Imin][j].Ez;
			p[Imin][j].Ez=EAB[2][j]+((1-Wl1)/(1+Wl1))*T2+c*u/(2+2*Wl1)*T1;
			EAB[2][j]=p[Imin+1][j].Ez; 
			 
			}

			for(i=Imin+1;i<=Imax-1;i++)
			{
			    /*left boundary*/
			T1=p[i-1][Jmin].Hx-p[i][Jmin].Hx+p[i-1][Jmin+1].Hx-p[i][Jmin+1].Hx;
            T2=p[i][Jmin+1].Ez-p[i][Jmin].Ez;
			p[i][Jmin].Ez=EAB[3][i]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[3][i]=p[i][Jmin+1].Ez; 
			 
			    
			   /*right boundary*/
			T1=p[i-1][Jmax].Hx-p[i][Jmax].Hx+p[i-1][Jmax-1].Hx-p[i][Jmax-1].Hx;
            T2=p[i][Jmax-1].Ez-p[i][Jmax].Ez;
			p[i][Jmax].Ez=EAB[1][i]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[1][i]=p[i][Jmax-1].Ez;  
            
			}
			
		/*corner values*/
			Alter=1-Alter;

			/*bottom left*/
			p[Imax][Jmin].Ez=0.2928923219*EAC[0][0][Alter]+0.7071067812*EAC[0][1][Alter];
			EAC[0][0][Alter]=p[Imax][Jmin].Ez;
			EAC[0][1][Alter]=p[Imax-1][Jmin+1].Ez;

			
			/*bottom right*/
			p[Imax][Jmax].Ez=0.2928923219*EAC[1][0][Alter]+0.7071067812*EAC[1][1][Alter];
			EAC[1][0][Alter]=p[Imax][Jmax].Ez;
			EAC[1][1][Alter]=p[Imax-1][Jmax-1].Ez;

			
			/*top right*/
			p[Imin][Jmax].Ez=0.2928923219*EAC[2][0][Alter]+0.7071067812*EAC[2][1][Alter];
			EAC[2][0][Alter]=p[Imin][Jmax].Ez;
			EAC[2][1][Alter]=p[Imin+1][Jmax-1].Ez;

			
			/*top left*/
			p[Imin][Jmin].Ez=0.2928923219*EAC[3][0][Alter]+0.7071067812*EAC[3][1][Alter];
			EAC[3][0][Alter]=p[Imin][Jmin].Ez;
			EAC[3][1][Alter]=p[Imin+1][Jmin+1].Ez;

		/*calculate Hx,Hy using FDTD*/
           for(i=Imin;i<Imax;i++)
			   for(j=Jmin;j<=Jmax;j++)
			   p[i][j].Hx=p[i][j].Hx-coeff_H*(p[i][j].Ez-p[i+1][j].Ez);
			   
			   
              
		   for(i=Imin;i<=Imax;i++)
			   for(j=Jmin;j<Jmax;j++)
			   p[i][j].Hy=p[i][j].Hy+coeff_H*(p[i][j+1].Ez-p[i][j].Ez);
			   
				
			
	}
	     time2=time(NULL);
        t=difftime(time2,time1);
		

		//export result to file

		for(i=Imin;i<=Imax;i=i+5)
			for(j=Jmin;j<=Jmax;j=j+5)		
			fprintf(fp,"%d    %d          %f\n",i,j,p[i][j].Ez);
			fprintf(fp,"%f\n",t);		

 
		if(fclose(fp))
		{
              printf("Can not close file a.txt!\n"); 
              
		}
		return 0;
}
