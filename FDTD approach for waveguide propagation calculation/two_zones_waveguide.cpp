#include <stdio.h>
#include <math.h>
#include  <time.h>


struct point{
	double Ez;
	double Hx;
	double Hy;
}p1[1000][3000],p2[2000][1000];                                                   /*�Գ������ýṹ���洢*/

int main()
{
	/*�������*/
	 time_t time1;
     time_t time2;          
	int Imax1,Imin1,Jmax1,Jmin1,Imax2,Imin2,Jmax2,Jmin2;
	int Wl,Wl1;
	int i,j;
	int step,Timestop;
	int Alter;                                       /*�ڽǵ㴦��������ѡ�����*/
    int Jsource,Ismin,Ismax;
	double n0;
	double m[2]={1.0,3.0};
	double Eps;
	double EAB[4][3000];                              /*�����洢�߽�*/
	double EAC[5][2][2];                             /*�����洢�ǵ�*/
	double Ei;                                        /*���弤��Դ*/
	double T1,T2;
    double ia,jb;
	double wavelength;    /*���岨��*/
	double t;
   double r;                                              //���������뾶  
   double d;                                              //��������
   double temp;
   double PI=3.1415926535895;
   double u=PI*0.0000004;                                 /*���ʽ��ϵ��*/  
   double Eps0=8.85e-12;                            /*�ŵ�ϵ��*/
   double c=299792457.4;                               /*����*/
   double coeff_E;                          /*�糡�����ϵ��*/
   double coeff_H=1/c/u/2;                           /*�ų������ϵ��*/ 
    
	FILE *fp;
	if((fp=fopen("two_zones_waveguide.txt","w"))==NULL)
		printf("Can not open the file!\n");

	/*������ʼ����ֵ*/
	Wl=10;
	Wl1=2;
	Imin1=0;
	Imax1=999;
	Jmin1=0;
	Jmax1=2999;
	Imin2=0;
	Imax2=1999;
	Jmin2=0;
	Jmax2=999;
	Timestop=13000;
	Alter=0;
	Jsource=1;
	Ismin=95;
	Ismax=105;
	r=0.0004495;
    d=0.00000155;
    wavelength=1.55e-6;
	 
	 time1=time(NULL);

	for(i=Imin1;i<=Imax1;i++)             //����1��ʼ��
		for(j=Jmin1;j<=Jmax1;j++)
		{
			p1[i][j].Ez=0.0;
			p1[i][j].Hx=0.0;
			p1[i][j].Hy=0.0;
		}

    for(i=Imin2;i<=Imax2;i++)              //����2��ʼ��
		for(j=Jmin2;j<=Jmax2;j++)
		{
			p2[i][j].Ez=0.0;
			p2[i][j].Hx=0.0;
			p2[i][j].Hy=0.0;
		}
    
	for(i=0;i<4;i++)
		for(j=0;j<=2999;j++)
			EAB[i][j]=0.0;

	for(i=0;i<5;i++)
		for(j=0;j<2;j++)
		{
			EAC[i][j][0]=0.0;
			EAC[i][j][1]=0.0;
		}
  
	/*ʱ�䲽��ѭ��*/

    for(step=1;step<=Timestop;step++)
	{
		/*����Դ����*/
		Ei=sin(step*2*PI/Wl/Wl1);
		for(i=Ismin;i<=Ismax;i++)
		{
			p1[i][Jsource].Ez=p1[i][Jsource].Ez+Ei; }

		printf("Step %d is processing!\n",step);
		

		/*Ez��ȫ����FDTD�����������ձ߽��⣩*/
		for(i=Imin1+1;i<=Imax1-1;i++)
			for(j=Jmin1+1;j<=Jmax1-1;j++)
			{
				ia=(double)(2999-i);
			    jb=(double)j;
				temp=sqrt(ia*ia+jb*jb); 
				temp=temp*wavelength/Wl;
				if((temp>=r-d/2)&&(temp<=r+d/2))
			                  n0=m[1];
			      else 
					          n0=m[0]; 
				  
				   Eps=Eps0*sqrt(n0);
                   coeff_E=1/c/Eps/2;
				p1[i][j].Ez=p1[i][j].Ez+coeff_E*(p1[i][j].Hy-p1[i][j-1].Hy-p1[i-1][j].Hx+p1[i][j].Hx);
			}
		
		for(j=2000;j<=Jmax1-1;j++)
		{
			    ia=2000.0;
			    jb=(double)j;
				temp=sqrt(ia*ia+jb*jb); 
			temp=temp*wavelength/Wl;
				if((temp>=r-d/2)&&(temp<=r+d/2))
			                  n0=m[1];
			      else 
					          n0=m[0]; 
				  
				   Eps=Eps0*sqrt(n0);
                   coeff_E=1/c/Eps/2;
                p1[Imax1][j].Ez=p1[Imax1][j].Ez+coeff_E*(p1[Imax1][j].Hy-p1[Imax1][j-1].Hy-p1[Imax1-1][j].Hx+p1[Imax1][j].Hx);
		}

        for(j=Jmin2+1;j<Jmax2-1;j++)
		{       ia=1999.0;
			    jb=(double)(j+2000);
				temp=sqrt(ia*ia+jb*jb); 
				temp=temp*wavelength/Wl;
				if((temp>=r-d/2)&&(temp<=r+d/2))
			                  n0=m[1];
			      else 
					          n0=m[0]; 
				  
				   Eps=Eps0*sqrt(n0);
                   coeff_E=1/c/Eps/2;
			p2[Imin2][j].Ez=p2[Imin2][j].Ez+coeff_E*(p2[Imin2][j].Hy-p2[Imin2][j-1].Hy-p1[Imax1][j+2000].Hx+p2[Imin2][j].Hx);
		}

		for(i=Imin2+1;i<=Imax2-1;i++)
			for(j=Jmin2+1;j<=Jmax2-1;j++)
			{   ia=(double)(1999-i);
			    jb=(double)(j+2000);
				temp=sqrt(ia*ia+jb*jb); 
				temp=temp*wavelength/Wl;
				if((temp>=r-d/2)&&(temp<=r+d/2))
			                  n0=m[1];
			      else 
					          n0=m[0]; 
				  
				   Eps=Eps0*sqrt(n0);
                   coeff_E=1/c/Eps/2;
				p2[i][j].Ez=p2[i][j].Ez+coeff_E*(p2[i][j].Hy-p2[i][j-1].Hy-p2[i-1][j].Hx+p2[i][j].Hx);
			}

	
				
    	/*����Mur���ձ߽磬�ıߴ���*/
			for(j=Jmin1+1;j<=Jmax1-1;j++)
			{
				
			    /*���߽�*/
			T1=p1[Imin1][j].Hy-p1[Imin1][j-1].Hy+p1[Imin1+1][j].Hy-p1[Imin1+1][j-1].Hy;
            T2=p1[Imin1+1][j].Ez-p1[Imin1][j].Ez;
			p1[Imin1][j].Ez=EAB[2][j]+((1-Wl1)/(1+Wl1))*T2+c*u/(2+2*Wl1)*T1;
			EAB[2][j]=p1[Imin1+1][j].Ez; 
			 
			}
   
            /*�ױ߽�*/
			for(j=Jmin1+1;j<=1999;j++)
			{
				
            T1=p1[Imax1][j].Hy-p1[Imax1][j-1].Hy+p1[Imax1-1][j].Hy-p1[Imax1-1][j-1].Hy;
            T2=p1[Imax1-1][j].Ez-p1[Imax1][j].Ez;
			p1[Imax1][j].Ez=EAB[0][j]+((1-Wl1)/(1+Wl1))*T2+c*u/(2+2*Wl1)*T1;
			EAB[0][j]=p1[Imax1-1][j].Ez; 
			
			}

			for(j=Jmin2+1;j<=Jmax2-1;j++)
			{
            T1=p2[Imax2][j].Hy-p2[Imax2][j-1].Hy+p2[Imax2-1][j].Hy-p2[Imax2-1][j-1].Hy;
            T2=p2[Imax2-1][j].Ez-p2[Imax2][j].Ez;
			p2[Imax2][j].Ez=EAB[0][j+2000]+((1-Wl1)/(1+Wl1))*T2+c*u/(2+2*Wl1)*T1;
			EAB[0][j+2000]=p2[Imax2-1][j].Ez;
			}
             
			     
			/*��߽�*/
			for(i=Imin1+1;i<=Imax1-1;i++)
			{
			    
			T1=p1[i-1][Jmin1].Hx-p1[i][Jmin1].Hx+p1[i-1][Jmin1+1].Hx-p1[i][Jmin1+1].Hx;
            T2=p1[i][Jmin1+1].Ez-p1[i][Jmin1].Ez;
			p1[i][Jmin1].Ez=EAB[3][i]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[3][i]=p1[i][Jmin1+1].Ez; 			 			   
            
			}   
			
            T1=p1[Imax1][Jmin2+2000].Hx-p2[Imin2][Jmin2].Hx+p1[Imax1][Jmin2+2001].Hx-p2[Imin2][Jmin2+1].Hx;
            T2=p2[Imin2][Jmin2+1].Ez-p2[Imin2][Jmin2].Ez;
			p2[Imin2][Jmin2].Ez=EAB[3][i+1000]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[3][i+1000]=p2[Imin2][Jmin2+1].Ez; 

			for(i=Imin2+1;i<=Imax2-1;i++)
			{
            T1=p2[i-1][Jmin2].Hx-p2[i][Jmin2].Hx+p2[i-1][Jmin2+1].Hx-p2[i][Jmin2+1].Hx;
            T2=p2[i][Jmin2+1].Ez-p2[i][Jmin2].Ez;
			p2[i][Jmin2].Ez=EAB[3][i+1000]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[3][i+1000]=p2[i][Jmin2+1].Ez; 
			}
			     
           /*�ұ߽�*/
            for(i=Imin1+1;i<=Imax1;i++)
			{
				 			  
			T1=p1[i-1][Jmax1].Hx-p1[i][Jmax1].Hx+p1[i-1][Jmax1-1].Hx-p1[i][Jmax1-1].Hx;
            T2=p1[i][Jmax1-1].Ez-p1[i][Jmax1].Ez;
			p1[i][Jmax1].Ez=EAB[1][i]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[1][i]=p1[i][Jmax1-1].Ez;  
			}

           	T1=p1[Imax1][Jmax1].Hx-p2[Imin2][Jmax2].Hx+p1[Imax1][Jmax1-1].Hx-p2[Imin2][Jmax2-1].Hx;
            T2=p2[Imin2][Jmax2-1].Ez-p2[Imin2][Jmax2].Ez;
			p2[Imin2][Jmax2].Ez=EAB[1][1000]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[1][1000]=p2[i][Jmax2-1].Ez;

			for(i=Imin2+1;i<=Imax2-1;i++)
			{
				 			  
			T1=p2[i-1][Jmax2].Hx-p2[i][Jmax2].Hx+p2[i-1][Jmax2-1].Hx-p2[i][Jmax2-1].Hx;
            T2=p2[i][Jmax2-1].Ez-p2[i][Jmax2].Ez;
			p2[i][Jmax2].Ez=EAB[1][i+1000]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
			EAB[1][i+1000]=p2[i][Jmax2-1].Ez;  
			}
			    
		/*�ǵ㴦��*/
			Alter=1-Alter;

			/*����*/
			p1[Imax1][Jmin1].Ez=0.2928923219*EAC[0][0][Alter]+0.7071067812*EAC[0][1][Alter];
			EAC[0][0][Alter]=p1[Imax1][Jmin1].Ez;
			EAC[0][1][Alter]=p1[Imax1-1][Jmin1+1].Ez;
		
			/*����*/
            p2[Imax2][Jmin2].Ez=0.2928923219*EAC[4][0][Alter]+0.7071067812*EAC[4][1][Alter];
			EAC[4][0][Alter]=p2[Imax2][Jmin2].Ez;
			EAC[4][1][Alter]=p2[Imax2-1][Jmin2+1].Ez;
			
			/*����*/
			p2[Imax2][Jmax2].Ez=0.2928923219*EAC[1][0][Alter]+0.7071067812*EAC[1][1][Alter];
			EAC[1][0][Alter]=p2[Imax2][Jmax2].Ez;
			EAC[1][1][Alter]=p2[Imax2-1][Jmax2-1].Ez;

			
			/*����*/
			p1[Imin1][Jmax1].Ez=0.2928923219*EAC[2][0][Alter]+0.7071067812*EAC[2][1][Alter];
			EAC[2][0][Alter]=p1[Imin1][Jmax1].Ez;
			EAC[2][1][Alter]=p1[Imin1+1][Jmax1-1].Ez;

			
			/*����*/
			p1[Imin1][Jmin1].Ez=0.2928923219*EAC[3][0][Alter]+0.7071067812*EAC[3][1][Alter];
			EAC[3][0][Alter]=p1[Imin1][Jmin1].Ez;
			EAC[3][1][Alter]=p1[Imin1+1][Jmin1+1].Ez;

		/*Hx��Hy��ȫ����FDTD����*/
			//Hx
           for(i=Imin1;i<Imax1;i++)
			   for(j=Jmin1;j<=Jmax1;j++)
			   p1[i][j].Hx=p1[i][j].Hx-coeff_H*(p1[i][j].Ez-p1[i+1][j].Ez);
  
           for(j=2000;j<Jmax1;j++)
               p1[Imax1][j].Hx=p1[Imax1][j].Hx-coeff_H*(p1[Imax1][j].Ez-p2[Imin2][j-2000].Ez);
		 
		   for(i=Imin2;i<Imax2;i++)
			   for(j=Jmin2;j<=Jmax2;j++)
			   p2[i][j].Hx=p2[i][j].Hx-coeff_H*(p2[i][j].Ez-p2[i+1][j].Ez);	   
		     	   
          //Hy    
		   for(i=Imin1;i<=Imax1;i++)
			   for(j=Jmin1;j<Jmax1;j++)
			   p1[i][j].Hy=p1[i][j].Hy+coeff_H*(p1[i][j+1].Ez-p1[i][j].Ez);
			  
		   for(i=Imin2;i<=Imax2;i++)
			   for(j=Jmin2;j<Jmax2;j++)
			   p2[i][j].Hy=p2[i][j].Hy+coeff_H*(p2[i][j+1].Ez-p2[i][j].Ez);
				
		   	
	}
	     
			time2=time(NULL);
        t=difftime(time2,time1);

		//���������ļ�

		for(i=Imin1;i<=Imax1;i=i+5)
			for(j=Jmin1;j<=Jmax1;j=j+5)		
			fprintf(fp,"%d    %d          %f\n",i,j,p1[i][j].Ez);
		
		for(i=Imin2;i<=Imax2;i=i+5)
		{	for(j=0;j<2000;j=j+5)
		        fprintf(fp,"%d    %d          %f\n",i+1000,j,0.000000);
			for(j=Jmin2;j<=Jmax2;j=j+5)		
				fprintf(fp,"%d    %d          %f\n",i+1000,j+2000,p2[i][j].Ez);}
			

		fprintf(fp,"%f\n",t);	
        printf("End!\n");
		if(fclose(fp))
		{
              printf("Can not close file a.txt!\n"); 
              
		}
		return 0;
}
