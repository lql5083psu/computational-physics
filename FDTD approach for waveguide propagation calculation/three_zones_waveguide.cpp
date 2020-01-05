#include<stdio.h>
#include<math.h>
#include<time.h>

#define row1 500
#define list1 1550
#define row2  1000
#define list2 1150
#define row3  1500
#define list3 600
#define side  3000

struct point {
	double Ez;
	double Hx;
	double Hy;
}p1[row1][list1], p2[row2][list2], p3[row3][list3];

double EAB[4][side];                              //store boundary values


int main()
{
	int Imax1, Imin1, Jmax1, Jmin1, Imax2, Imin2, Jmax2, Jmin2, Imax3, Imin3, Jmax3, Jmin3;
	int Wl, Wl1;
	int i, j;
	int step, Timestop;
	int Alter;                          //select alternative variale for corner
	int Jsource, Ismin, Ismax;
	double n0;
	double m[2] = { 1.0,3.0 };
	double Eps;
	double EAC[8][2][2];                             //store corner values
	double Ei;                                        //define source
	double T1, T2;
	double ia, jb;
	double wavelength;    //define wavelength
	double t;
	double r;                                              //waveguide length  
	double d;                                              //waveguide width
	double temp;
	double PI = 3.1415926535895;
	double u = PI * 0.0000004;                                 /*permeability*/ 
	double Eps0 = 8.85e-12;                            /*permittivity*/
	double c = 299792457.4;                               /*speed of light*/
	double coeff_E;                          /*coefficient of electric field*/
	double coeff_H = 1/c/u/2;           /*coefficient of electric field*/

	FILE *fp;


	//initialize variables
	Wl = 10;
	Wl1 = 2;
	Imin1 = 0;
	Imax1 = row1 - 1;
	Jmin1 = 0;
	Jmax1 = list1 - 1;
	Imin2 = 0;
	Imax2 = row2 - 1;
	Jmin2 = 0;
	Jmax2 = list2 - 1;
	Imin3 = 0;
	Imax3 = row3 - 1;
	Jmin3 = 0;
	Jmax3 = list3 - 1;
	Timestop = 13000;
	Alter = 0;
	Jsource = 1;
	r = 0.0004495;
	d = 0.00000155;
	wavelength = 1.55e-6;
	Ismin = 95;
	Ismax = 105;

	time_t time1;
	time_t time2;
	time1 = time(NULL);

	for (i = Imin1; i <= Imax1; i++)             //initialize zone 1
		for (j = Jmin1; j <= Jmax1; j++)
		{
			p1[i][j].Ez = 0.0;
			p1[i][j].Hx = 0.0;
			p1[i][j].Hy = 0.0;
		}

	for (i = Imin2; i <= Imax2; i++)              //initialize zone 2
		for (j = Jmin2; j <= Jmax2; j++)
		{
			p2[i][j].Ez = 0.0;
			p2[i][j].Hx = 0.0;
			p2[i][j].Hy = 0.0;
		}

	for (i = Imin3; i <= Imax3; i++)              //initialize zone 3
		for (j = Jmin3; j <= Jmax3; j++)
		{
			p3[i][j].Ez = 0.0;
			p3[i][j].Hx = 0.0;
			p3[i][j].Hy = 0.0;
		}

	for (i = 0; i < 4; i++)
		for (j = 0; j <= side - 1; j++)
			EAB[i][j] = 0.0;

	for (i = 0; i < 8; i++)
		for (j = 0; j < 2; j++)
		{
			EAC[i][j][0] = 0.0;
			EAC[i][j][1] = 0.0;
		}

	//time based loop

	for (step = 1; step <= Timestop; step++)
	{
		//set up source
		Ei = sin(step * 2 * PI / Wl / Wl1);
		for (i = Ismin; i <= Ismax; i++)
		{
			p1[i][Jsource].Ez = p1[i][Jsource].Ez + Ei;
		}

		printf("Step %d is processing!\n", step);


		//calcaulate Ez using FDTD except boundaries
		//zone 1
		for (i = Imin1 + 1; i <= Imax1 - 1; i++)
			for (j = Jmin1 + 1; j <= Jmax1 - 1; j++)
			{
				ia = (double)(side - i);
				jb = (double)j;
				temp = sqrt(ia*ia + jb * jb);
				temp = temp * wavelength / Wl;
				if ((temp >= r - d / 2) && (temp <= r + d / 2))
					n0 = m[1];
				else
					n0 = m[0];

				Eps = Eps0 * sqrt(n0);
				coeff_E = 1 / c / Eps / 2;
				p1[i][j].Ez = p1[i][j].Ez + coeff_E * (p1[i][j].Hy - p1[i][j - 1].Hy - p1[i - 1][j].Hx + p1[i][j].Hx);
			}

		for (j = 1400; j <= Jmax1 - 1; j++)
		{
			ia = (double)(side - row1);
			jb = (double)j;
			temp = sqrt(ia*ia + jb * jb);
			temp = temp * wavelength / Wl;
			if ((temp >= r - d / 2) && (temp <= r + d / 2))
				n0 = m[1];
			else
				n0 = m[0];

			Eps = Eps0 * sqrt(n0);
			coeff_E = 1 / c / Eps / 2;
			p1[Imax1][j].Ez = p1[Imax1][j].Ez + coeff_E * (p1[Imax1][j].Hy - p1[Imax1][j - 1].Hy - p1[Imax1 - 1][j].Hx + p1[Imax1][j].Hx);
		}

		//zone 2
		for (j = Jmin2 + 1; j < 150; j++)
		{
			ia = (double)(side - row1);
			jb = (double)(j + 1400);
			temp = sqrt(ia*ia + jb * jb);
			temp = temp * wavelength / Wl;
			if ((temp >= r - d / 2) && (temp <= r + d / 2))
				n0 = m[1];
			else
				n0 = m[0];

			Eps = Eps0 * sqrt(n0);
			coeff_E = 1 / c / Eps / 2;
			p2[Imin2][j].Ez = p2[Imin2][j].Ez + coeff_E * (p2[Imin2][j].Hy - p2[Imin2][j - 1].Hy - p1[Imax1][j + 1400].Hx + p2[Imin2][j].Hx);
		}

		for (i = Imin2 + 1; i <= Imax2 - 1; i++)
			for (j = Jmin2 + 1; j <= Jmax2 - 1; j++)
			{
				ia = (double)(2500 - i);
				jb = (double)(j + 1400);
				temp = sqrt(ia*ia + jb * jb);
				temp = temp * wavelength / Wl;
				if ((temp >= r - d / 2) && (temp <= r + d / 2))
					n0 = m[1];
				else
					n0 = m[0];

				Eps = Eps0 * sqrt(n0);
				coeff_E = 1 / c / Eps / 2;
				p2[i][j].Ez = p2[i][j].Ez + coeff_E * (p2[i][j].Hy - p2[i][j - 1].Hy - p2[i - 1][j].Hx + p2[i][j].Hx);
			}


		for (j = 1000; j <= Jmax2 - 1; j++)
		{
			ia = 1500.0;
			jb = (double)(j + 1400);
			temp = sqrt(ia*ia + jb * jb);
			temp = temp * wavelength / Wl;
			if ((temp >= r - d / 2) && (temp <= r + d / 2))
				n0 = m[1];
			else
				n0 = m[0];

			Eps = Eps0 * sqrt(n0);
			coeff_E = 1 / c / Eps / 2;
			p2[Imax2][j].Ez = p2[Imax2][j].Ez + coeff_E * (p2[Imax2][j].Hy - p2[Imax2][j - 1].Hy - p2[Imax2 - 1][j].Hx + p2[Imax2][j].Hx);
		}

		//zone 3
		for (j = Jmin3 + 1; j < 150; j++)
		{
			ia = (double)(side - row1 - row2);
			jb = (double)(j + 2400);
			temp = sqrt(ia*ia + jb * jb);
			temp = temp * wavelength / Wl;
			if ((temp >= r - d / 2) && (temp <= r + d / 2))
				n0 = m[1];
			else
				n0 = m[0];

			Eps = Eps0 * sqrt(n0);
			coeff_E = 1 / c / Eps / 2;
			p3[Imin3][j].Ez = p3[Imin3][j].Ez + coeff_E * (p3[Imin3][j].Hy - p3[Imin3][j - 1].Hy - p2[Imax2][j + 1000].Hx + p3[Imin3][j].Hx);
		}

		for (i = Imin3 + 1; i <= Imax3 - 1; i++)
			for (j = Jmin3 + 1; j <= Jmax3 - 1; j++)
			{
				ia = (double)(1500 - i);
				jb = (double)(j + 2400);
				temp = sqrt(ia*ia + jb * jb);
				temp = temp * wavelength / Wl;
				if ((temp >= r - d / 2) && (temp <= r + d / 2))
					n0 = m[1];
				else
					n0 = m[0];

				Eps = Eps0 * sqrt(n0);
				coeff_E = 1 / c / Eps / 2;
				p3[i][j].Ez = p3[i][j].Ez + coeff_E * (p3[i][j].Hy - p3[i][j - 1].Hy - p3[i - 1][j].Hx + p3[i][j].Hx);
			}


		//second-order Mur absorption boundary			
		//���߽�
		//zone 1
		for (j = Jmin1 + 1; j <= Jmax1 - 1; j++)
		{
			T1 = p1[Imin1][j].Hy - p1[Imin1][j - 1].Hy + p1[Imin1 + 1][j].Hy - p1[Imin1 + 1][j - 1].Hy;
			T2 = p1[Imin1 + 1][j].Ez - p1[Imin1][j].Ez;
			p1[Imin1][j].Ez = EAB[2][j] + ((1 - Wl1) / (1 + Wl1))*T2 + c * u / (2 + 2 * Wl1)*T1;
			EAB[2][j] = p1[Imin1 + 1][j].Ez;

		}

		//zone 2
		for (j = 150; j <= Jmax2 - 1; j++)
		{
			T1 = p2[Imin2][j].Hy - p2[Imin2][j - 1].Hy + p2[Imin2 + 1][j].Hy - p2[Imin2 + 1][j - 1].Hy;
			T2 = p2[Imin2 + 1][j].Ez - p2[Imin2][j].Ez;
			p2[Imin2][j].Ez = EAB[2][j + 1400] + ((1 - Wl1) / (1 + Wl1))*T2 + c * u / (2 + 2 * Wl1)*T1;
			EAB[2][j + 1400] = p2[Imin2 + 1][j].Ez;

		}

		//zone 3
		for (j = 150; j <= Jmax3 - 1; j++)
		{
			T1 = p3[Imin3][j].Hy - p3[Imin3][j - 1].Hy + p3[Imin3 + 1][j].Hy - p3[Imin3 + 1][j - 1].Hy;
			T2 = p3[Imin3 + 1][j].Ez - p3[Imin3][j].Ez;
			p3[Imin3][j].Ez = EAB[2][j + 2400] + ((1 - Wl1) / (1 + Wl1))*T2 + c * u / (2 + 2 * Wl1)*T1;
			EAB[2][j + 2400] = p3[Imin3 + 1][j].Ez;

		}

		//�ױ߽�
		//zone 1
		for (j = Jmin1 + 1; j <= 1399; j++)
		{

			T1 = p1[Imax1][j].Hy - p1[Imax1][j - 1].Hy + p1[Imax1 - 1][j].Hy - p1[Imax1 - 1][j - 1].Hy;
			T2 = p1[Imax1 - 1][j].Ez - p1[Imax1][j].Ez;
			p1[Imax1][j].Ez = EAB[0][j] + ((1 - Wl1) / (1 + Wl1))*T2 + c * u / (2 + 2 * Wl1)*T1;
			EAB[0][j] = p1[Imax1 - 1][j].Ez;

		}

		//zone 2
		for (j = Jmin2 + 1; j <= 999; j++)
		{
			T1 = p2[Imax2][j].Hy - p2[Imax2][j - 1].Hy + p2[Imax2 - 1][j].Hy - p2[Imax2 - 1][j - 1].Hy;
			T2 = p2[Imax2 - 1][j].Ez - p2[Imax2][j].Ez;
			p2[Imax2][j].Ez = EAB[0][j + 1400] + ((1 - Wl1) / (1 + Wl1))*T2 + c * u / (2 + 2 * Wl1)*T1;
			EAB[0][j + 1400] = p2[Imax2 - 1][j].Ez;
		}

		//zone 3
		for (j = Jmin3 + 1; j <= Jmax3 - 1; j++)
		{
			T1 = p3[Imax3][j].Hy - p3[Imax3][j - 1].Hy + p3[Imax3 - 1][j].Hy - p3[Imax3 - 1][j - 1].Hy;
			T2 = p3[Imax3 - 1][j].Ez - p3[Imax3][j].Ez;
			p3[Imax3][j].Ez = EAB[0][j + 2400] + ((1 - Wl1) / (1 + Wl1))*T2 + c * u / (2 + 2 * Wl1)*T1;
			EAB[0][j + 2400] = p3[Imax3 - 1][j].Ez;
		}

		//��߽�
		//zone 1
		for (i = Imin1 + 1; i <= Imax1 - 1; i++)
		{

			T1 = p1[i - 1][Jmin1].Hx - p1[i][Jmin1].Hx + p1[i - 1][Jmin1 + 1].Hx - p1[i][Jmin1 + 1].Hx;
			T2 = p1[i][Jmin1 + 1].Ez - p1[i][Jmin1].Ez;
			p1[i][Jmin1].Ez = EAB[3][i] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
			EAB[3][i] = p1[i][Jmin1 + 1].Ez;

		}

		T1 = p1[Imax1][Jmin2 + 1400].Hx - p2[Imin2][Jmin2].Hx + p1[Imax1][Jmin2 + 1401].Hx - p2[Imin2][Jmin2 + 1].Hx;
		T2 = p2[Imin2][Jmin2 + 1].Ez - p2[Imin2][Jmin2].Ez;
		p2[Imin2][Jmin2].Ez = EAB[3][i + 500] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
		EAB[3][i + 500] = p2[Imin2][Jmin2 + 1].Ez;

		//zone 2
		for (i = Imin2 + 1; i <= Imax2 - 1; i++)
		{
			T1 = p2[i - 1][Jmin2].Hx - p2[i][Jmin2].Hx + p2[i - 1][Jmin2 + 1].Hx - p2[i][Jmin2 + 1].Hx;
			T2 = p2[i][Jmin2 + 1].Ez - p2[i][Jmin2].Ez;
			p2[i][Jmin2].Ez = EAB[3][i + 500] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
			EAB[3][i + 500] = p2[i][Jmin2 + 1].Ez;
		}

		T1 = p2[Imax2][Jmin3 + 1000].Hx - p3[Imin3][Jmin3].Hx + p2[Imax2][Jmin3 + 1001].Hx - p3[Imin3][Jmin3 + 1].Hx;
		T2 = p3[Imin3][Jmin3 + 1].Ez - p3[Imin3][Jmin3].Ez;
		p3[Imin3][Jmin3].Ez = EAB[3][i + 1500] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
		EAB[3][i + 1500] = p3[Imin3][Jmin3 + 1].Ez;

		//zone 3
		for (i = Imin3 + 1; i <= Imax3 - 1; i++)
		{
			T1 = p3[i - 1][Jmin3].Hx - p3[i][Jmin3].Hx + p3[i - 1][Jmin3 + 1].Hx - p3[i][Jmin3 + 1].Hx;
			T2 = p3[i][Jmin3 + 1].Ez - p3[i][Jmin3].Ez;
			p3[i][Jmin3].Ez = EAB[3][i + 1500] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
			EAB[3][i + 1500] = p3[i][Jmin3 + 1].Ez;
		}

		//�ұ߽�
		 //zone 1
		for (i = Imin1 + 1; i <= Imax1; i++)
		{

			T1 = p1[i - 1][Jmax1].Hx - p1[i][Jmax1].Hx + p1[i - 1][Jmax1 - 1].Hx - p1[i][Jmax1 - 1].Hx;
			T2 = p1[i][Jmax1 - 1].Ez - p1[i][Jmax1].Ez;
			p1[i][Jmax1].Ez = EAB[1][i] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
			EAB[1][i] = p1[i][Jmax1 - 1].Ez;
		}


		//zone 2
		for (i = Imin2 + 1; i <= Imax2; i++)
		{

			T1 = p2[i - 1][Jmax2].Hx - p2[i][Jmax2].Hx + p2[i - 1][Jmax2 - 1].Hx - p2[i][Jmax2 - 1].Hx;
			T2 = p2[i][Jmax2 - 1].Ez - p2[i][Jmax2].Ez;
			p2[i][Jmax2].Ez = EAB[1][i + 500] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
			EAB[1][i + 500] = p2[i][Jmax2 - 1].Ez;
		}

		/*T1=p1[Imax1][Jmax1].Hx-p2[Imin2][Jmax2].Hx+p1[Imax1][Jmax1-1].Hx-p2[Imin2][Jmax2-1].Hx;
		T2=p2[Imin2][Jmax2-1].Ez-p2[Imin2][Jmax2].Ez;
		p2[Imin2][Jmax2].Ez=EAB[1][5000]+((1-Wl1)/(1+Wl1))*T2-c*u/(2+2*Wl1)*T1;
		EAB[1][5000]=p2[i][Jmax2-1].Ez;*/

		//zone 3
		for (i = Imin3 + 1; i <= Imax3 - 1; i++)
		{

			T1 = p3[i - 1][Jmax3].Hx - p3[i][Jmax3].Hx + p3[i - 1][Jmax3 - 1].Hx - p3[i][Jmax3 - 1].Hx;
			T2 = p3[i][Jmax3 - 1].Ez - p3[i][Jmax3].Ez;
			p3[i][Jmax3].Ez = EAB[1][i + 1500] + ((1 - Wl1) / (1 + Wl1))*T2 - c * u / (2 + 2 * Wl1)*T1;
			EAB[1][i + 1500] = p3[i][Jmax3 - 1].Ez;
		}

		//corner values
		Alter = 1 - Alter;

		//����1�Ľǵ�
		//����
		p1[Imax1][Jmin1].Ez = 0.2928923219*EAC[0][0][Alter] + 0.7071067812*EAC[0][1][Alter];
		EAC[0][0][Alter] = p1[Imax1][Jmin1].Ez;
		EAC[0][1][Alter] = p1[Imax1 - 1][Jmin1 + 1].Ez;

		//����
		p1[Imin1][Jmax1].Ez = 0.2928923219*EAC[2][0][Alter] + 0.7071067812*EAC[2][1][Alter];
		EAC[2][0][Alter] = p1[Imin1][Jmax1].Ez;
		EAC[2][1][Alter] = p1[Imin1 + 1][Jmax1 - 1].Ez;


		//����
		p1[Imin1][Jmin1].Ez = 0.2928923219*EAC[3][0][Alter] + 0.7071067812*EAC[3][1][Alter];
		EAC[3][0][Alter] = p1[Imin1][Jmin1].Ez;
		EAC[3][1][Alter] = p1[Imin1 + 1][Jmin1 + 1].Ez;

		//����2�Ľǵ�
	   //����
		p2[Imax2][Jmin2].Ez = 0.2928923219*EAC[4][0][Alter] + 0.7071067812*EAC[4][1][Alter];
		EAC[4][0][Alter] = p2[Imax2][Jmin2].Ez;
		EAC[4][1][Alter] = p2[Imax2 - 1][Jmin2 + 1].Ez;

		//����
		p2[Imin2][Jmax2].Ez = 0.2928923219*EAC[5][0][Alter] + 0.7071067812*EAC[5][1][Alter];
		EAC[5][0][Alter] = p2[Imin2][Jmax2].Ez;
		EAC[5][1][Alter] = p2[Imin2 + 1][Jmax2 - 1].Ez;

		//����3�Ľǵ�
		//����
		p3[Imax3][Jmax3].Ez = 0.2928923219*EAC[1][0][Alter] + 0.7071067812*EAC[1][1][Alter];
		EAC[1][0][Alter] = p3[Imax3][Jmax3].Ez;
		EAC[1][1][Alter] = p3[Imax3 - 1][Jmax3 - 1].Ez;

		//����
		p3[Imax3][Jmin3].Ez = 0.2928923219*EAC[6][0][Alter] + 0.7071067812*EAC[6][1][Alter];
		EAC[6][0][Alter] = p3[Imax3][Jmin3].Ez;
		EAC[6][1][Alter] = p3[Imax3 - 1][Jmin3 + 1].Ez;

		//����
		p3[Imin3][Jmax3].Ez = 0.2928923219*EAC[7][0][Alter] + 0.7071067812*EAC[7][1][Alter];
		EAC[7][0][Alter] = p3[Imin3][Jmax3].Ez;
		EAC[7][1][Alter] = p3[Imin3 + 1][Jmax3 - 1].Ez;




		/*calculate Hx,Hy using FDTD*/
			//Hx
		for (i = Imin1; i < Imax1; i++)
			for (j = Jmin1; j <= Jmax1; j++)
				p1[i][j].Hx = p1[i][j].Hx - coeff_H * (p1[i][j].Ez - p1[i + 1][j].Ez);

		for (j = 1400; j < Jmax1; j++)
			p1[Imax1][j].Hx = p1[Imax1][j].Hx - coeff_H * (p1[Imax1][j].Ez - p2[Imin2][j - 1400].Ez);

		for (i = Imin2; i < Imax2; i++)
			for (j = Jmin2; j <= Jmax2; j++)
				p2[i][j].Hx = p2[i][j].Hx - coeff_H * (p2[i][j].Ez - p2[i + 1][j].Ez);


		for (j = 1000; j < Jmax2; j++)
			p2[Imax2][j].Hx = p2[Imax2][j].Hx - coeff_H * (p2[Imax2][j].Ez - p3[Imin3][j - 1000].Ez);

		for (i = Imin3; i < Imax3; i++)
			for (j = Jmin3; j <= Jmax3; j++)
				p3[i][j].Hx = p3[i][j].Hx - coeff_H * (p3[i][j].Ez - p3[i + 1][j].Ez);

		//Hy    
		for (i = Imin1; i <= Imax1; i++)
			for (j = Jmin1; j < Jmax1; j++)
				p1[i][j].Hy = p1[i][j].Hy + coeff_H * (p1[i][j + 1].Ez - p1[i][j].Ez);

		for (i = Imin2; i <= Imax2; i++)
			for (j = Jmin2; j < Jmax2; j++)
				p2[i][j].Hy = p2[i][j].Hy + coeff_H * (p2[i][j + 1].Ez - p2[i][j].Ez);

		for (i = Imin3; i <= Imax3; i++)
			for (j = Jmin3; j < Jmax3; j++)
				p3[i][j].Hy = p3[i][j].Hy + coeff_H * (p3[i][j + 1].Ez - p3[i][j].Ez);



	}
	time2 = time(NULL);
	t = difftime(time2, time1);

	//errno_t err;

	//if ((err = fopen_s(&fp, "trip1000.txt", "w")) != 0) printf("Can not open ///the file!\n");
	if ((fp = fopen("three_zones_waveguide.txt", "w")) == NULL)
		printf("Can not open the file!\n");

	for (i = Imin1; i <= Imax1; i = i + 5)
	{
		for (j = Jmin1; j <= Jmax1; j = j + 5)
			fprintf(fp, "%d    %d          %f\n", i, j, p1[i][j].Ez);
		for (j = Jmax1 + 1; j < side; j = j + 5)
			fprintf(fp, "%d    %d          %f\n", i, j, 0.000000);
	}

	for (i = Imin2; i <= Imax2; i = i + 5)
	{
		for (j = 0; j < 1400; j = j + 5)
			fprintf(fp, "%d    %d          %f\n", i + 500, j, 0.000000);
		for (j = Jmin2; j <= Jmax2; j = j + 5)
			fprintf(fp, "%d    %d          %f\n", i + 500, j + 1400, p2[i][j].Ez);
		for (j = Jmax2 + 1 + 1400; j < side; j = j + 5)
			fprintf(fp, "%d    %d          %f\n", i + 500, j, 0.000000);
	}

	for (i = Imin3; i <= Imax3; i = i + 5)
	{
		for (j = 0; j < 2400; j = j + 5)
			fprintf(fp, "%d    %d          %f\n", i + 1500, j, 0.000000);
		for (j = Jmin3; j <= Jmax3; j = j + 5)
			fprintf(fp, "%d    %d          %f\n", i + 1500, j + 2400, p3[i][j].Ez);
	}


	fprintf(fp, "%f\n", t);

	if (fclose(fp))
	{
		printf("Can not close file a.txt!\n");

	}

	//���������ļ�




	printf("End!\n");
	return 0;

}