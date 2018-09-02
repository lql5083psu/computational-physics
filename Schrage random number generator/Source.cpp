#include<stdio.h>
#include<time.h>
#include <iostream>
//Schrage methond to generate random number and check the independence
int main()
{
	//identify all the parameters
	long a, q, r, m, j, i, I[2], N[10][10], k;
	int x, y;
	double x1, x2, c, sum1, sum2, sum3, Cl;
	time_t t;
	struct tm *pt;
	FILE *fp;
	errno_t err;

	//Store value into file D01 
	if ((err = fopen_s(&fp, "D01.txt", "w")) != 0)
	{
		printf("File open error!\n");
	}

	//intialize all parameters, x1 and x2 are coordinate of a point 
	a = 16807;
	m = 2147483647;
	r = m % a;
	q = m / a;
	x1 = 0;
	x2 = 0;
	k = 100;
	sum1 = 0;
	sum2 = 0;
	sum3 = 0;
	Cl = 0;
	c = 0;

	for (i = 0; i < 10; i++)
	{
		for (j = 0; j < 10; j++)
			N[i][j] = 0;
	}


	//obtain system time as seed of random number 
	pt = new tm();
	localtime_s(pt, &t);
	I[0] = pt->tm_year + 70 * (pt->tm_mon + 12 * (pt->tm_mday + 31 * (pt->tm_hour + 23 * (pt->tm_min + 59 * pt->tm_sec))));

	//function to obtain random number 
	for (j = 0; j < 5000; j++)
	{
		I[1] = a * (I[0] % q) - r * (I[0] / q);
		if (I[1] < 0)
			I[1] = I[1] + m;
		I[1] = I[1] % m;
		I[0] = I[1];                                  //I[0]is a number bweteen 0~m-1
		x2 = (double)I[0] / (double)m;                  //x2 is a number bwteen 0~1

														//use adjacent random numbers as coordinate
		if (j % 2 == 0)
			fprintf(fp, "%f      %f\n", x1, x2);

		x = (int)(x1 * 10);
		y = (int)(x2 * 10);
		N[x][y] = N[x][y] + 1;

		sum1 = sum1 + x2;
		sum2 = sum2 + x2 * x2;
		sum3 = sum3 + x1 * x2;

		x1 = x2;
	}

	for (i = 0; i < 10; i++)                          //calculate the independence 
	{
		for (j = 0; j < 10; j++)
			c = c + (N[i][j] - k)*(N[i][j] - k) / k;

	}

	sum3 = sum3 / 10000;
	sum2 = sum2 / 10000;
	sum1 = sum1 / 10000;
	Cl = (sum3 - sum1 * sum1) / (sum2 - sum1 * sum1);                 //calculate the linear correlation

	printf("c=%f,Cl=%f\n", c, Cl);

	if (fclose(fp))
	{
		printf("Can not close file a.txt!\n");
	}
	std::cin.ignore();
	return 0;
}
