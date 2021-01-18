#include <SDKDDKVer.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <direct.h>
#include <sys/stat.h>
#include <stdio.h>
#include <tchar.h>
#include <omp.h>
#include "Glegendre.h"
using namespace std;

#define MYU				1.2566370614e-6
#define PI				3.14159265358979323846264338327950288419716939937511
#define R0				0.05
#define TURNS			1
#define Div_C			4
#define Width			0.004
#define THICKNESS		0.002
#define MP_MAXNUM	    64

typedef struct _tagIntegralInterval
{
	double t_upper;
	double t_lower;
	double z_upper;
	double z_lower;
}Intergral_Inv;

double M(double ri, double rj, double ti, double zi, double tj, double zj);
double CalculateM(double ri, double rj, double ti_upp, double ti_low, double tj_upp, double tj_low, double* ft, double* fz, double* fw, int nNum);
double L(double a, double b, double c, double x, double y1, double y2, double myu, double d);
double y11(double x);
double y22(double x);

void InductanceCalculation(char* Filename);

//self inductance functions:
double y11(double x)
{
	double yut = 0;
	yut = PI * x / 3 - log(1 + x * x) / (12 * x * x) - x * x * log(1 + 1 / (x * x)) / 12 - 2 * (x - 1 / x) * (atan(x)) / 3 - 1 / 12;
	return yut;
}

double y22(double x)
{
	double yur = 0;
	yur = (69 / 20 + 221 / (60 * x * x) - log(1 + x * x) / (10 * x * x * x * x) + x * x * log(1 + 1 / (x * x)) / 2 - 8 * PI * x / 5 + 16 * x * atan(x) / 5) / 6;
	return yur;
}

double L(double a, double b, double c, double x, double y1, double y2, double myu, double d)
{
	double Lut = 0;
	Lut = myu * a * ((1 + (3 * b * b + c * c) / (96 * a * a)) * log(8 * a / d) - y1 + b * b * y2 / (16 * a * a));
	return Lut;
}

//mutual inductance functions:
double M(double ri, double rj, double ti, double zi, double tj, double zj)
{
	double Mut = 0;

	Mut = ri * rj * cos(tj - ti) / sqrt(pow(ri, 2) + pow(rj, 2) - 2 * ri * rj * cos(tj - ti) + pow(zj - zi, 2));

	return Mut;
}
//integrate function:
double CalculateM(double ri, double rj, Intergral_Inv ti, Intergral_Inv tj, Intergral_Inv zi, Intergral_Inv zj, double* fti, double* ftj, double* fzi, double* fzj, double* fw, int nNum)
{
	int p = 0;
	int q = 0;
	int k = 0;
	int l = 0;
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	double sum = 0;
	double Coff = 0;

	Coff = 1e-7 / (Width * Width);

	for (p = 0; p < nNum; ++p)
	{
		sum1 = 0;

		for (q = 0; q < nNum; ++q)
		{
			sum2 = 0;

			for (k = 0; k < nNum; ++k)
			{
				sum3 = 0;

				for (l = 0; l < nNum; ++l)
				{
#pragma omp atomic
					sum3 += (zj.z_upper - zj.z_lower) / 2 * fw[l + 1] * M(ri, rj, (ti.t_upper + ti.t_lower) / 2 + (ti.t_upper - ti.t_lower) * (fti[p + 1] / 2), (zi.z_upper + zi.z_lower) / 2 + (zi.z_upper - zi.z_lower) * (fzi[q + 1] / 2), (tj.t_upper + tj.t_lower) / 2 + (tj.t_upper - tj.t_lower) * (ftj[k + 1] / 2), (zj.z_upper + zj.z_lower) / 2 + (zj.z_upper - zj.z_lower) * (fzj[l + 1] / 2));
				}

				sum2 += (tj.t_upper - tj.t_lower) / 2 * fw[k + 1] * sum3;
			}

			sum1 += (zi.z_upper - zi.z_lower) / 2 * fw[q + 1] * sum2;
		}

		sum += (ti.t_upper - ti.t_lower) / 2 * fw[p + 1] * sum1;
	}

	return Coff * sum;
}

void InductanceCalculation(char* Filename)
{
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	double ri = 0;
	double rj = 0;
	int nRetCode = 0;
	double fpoint[MP_MAXNUM + 2] = { 0 };
	double fweight[MP_MAXNUM + 2] = { 0 };

	Intergral_Inv ti;
	Intergral_Inv tj;
	Intergral_Inv zi;
	Intergral_Inv zj;

	double** Inductance_Matrix = new double* [Div_C * TURNS];//&取地址运算符 *取值运算符

	FILE* fpFile = NULL;

	cout << "-----------------Self and Mutual Inductance Calculation Start----------------\n" << endl;

	for (i = 0; i < Div_C * TURNS; ++i)
	{
		Inductance_Matrix[i] = new double[Div_C * TURNS];
		//分配一个存放double型一维数组的内存空间，该数组有???个元素，并将该数组的地址赋给Inductance_Matrix[i] delete/delete[]释放由new分配的内存空间
	}

	nRetCode = gauss_legendre(MP_MAXNUM, fpoint, fweight);

	for (i = 0; i < Div_C * TURNS; ++i)
	{
		ri = R0 + THICKNESS * (i / Div_C + 0.5);
		ti.t_lower = 2*PI / Div_C * i;
		ti.t_upper = 2*PI / Div_C * (static_cast<unsigned __int64>(i) + 1);

#pragma omp parallel for num_threads(4) private(rj,tj,zj)
		for (j = 0; j < Div_C * TURNS; ++j)
		{
			rj = R0 + THICKNESS * (j / Div_C + 0.5);
			tj.t_lower = 2*PI / Div_C * j;
			tj.t_upper = 2*PI / Div_C * (static_cast<unsigned __int64>(j) + 1);
			zi.z_lower = 0;
			zi.z_upper = Width;
			zj.z_lower = 0;
			zj.z_upper = Width;

			Inductance_Matrix[i][j] = CalculateM(ri, rj, ti, tj, zi, zj, fpoint, fpoint, fpoint, fpoint, fweight, MP_MAXNUM);
			printf_s("Caculation Process: %d | %d................%f \r\n", i, j, Inductance_Matrix[i][j]);
			printf_s(" %f %f %f %f \r \n", ti.t_lower, ti.t_upper, tj.t_lower, tj.t_upper);
		}
	}

	//计算自感值填入原矩阵
	double a = 0;
	double b = 0;
	double c = 0;
	double x = 0;
	double y1 = 0;
	double y2 = 0;
	double d = 0;
	double L_circle = 0;
	double sum = 0;
	double L_equal = 0;

	for (k = 0; k < Div_C * TURNS; k += Div_C)
	{
		a = (R0 + k / Div_C * THICKNESS + THICKNESS / 2);
		b = Width;
		c = THICKNESS;
		x = b / c;
		d = sqrt(b * b + c * c);
		y1 = y11(x);
		y2 = y22(x);
		L_circle = L(a, b, c, x, y1, y2, MYU, d);
		for (i = 1; i < Div_C; i++)
		{
			for (j = 0; j < i; j++)
			{
				sum += Inductance_Matrix[k + i][k + j];
			}
		}

		L_equal = (L_circle - 2*sum) / Div_C;
		for (i = 0; i < Div_C; i++)
		{
			Inductance_Matrix[k + i][k + i] = L_equal;
		}
	}

	fopen_s(&fpFile, Filename, "w+");

	for (i = 0; i < Div_C * TURNS; ++i)
	{
		for (j = 0; j < Div_C * TURNS; ++j)
		{
			fprintf_s(fpFile, "%3.18lf\t ", Inductance_Matrix[i][j]);
		}

		fprintf_s(fpFile, "\n");
	}

	fclose(fpFile);
}