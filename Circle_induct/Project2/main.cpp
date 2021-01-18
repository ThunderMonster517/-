#include <stdio.h>
#include <stdlib.h>
#include "Glegendre.h"
#include <math.h>
#include <iostream>
#include <time.h>
using namespace std;

extern void InductanceCalculation(char* Filename);

int main()
{
	clock_t begin, end;
	double cost = 0;
	begin = clock();

	char Output_InductanceFile[50] = { 0 };
	sprintf_s(Output_InductanceFile, ".\\OUTPUT\\Results\\Inductance_Matrix.dat");

	InductanceCalculation(Output_InductanceFile);

	end = clock();
	cost = (double)(end - begin) / CLOCKS_PER_SEC / 60;
	printf("The total time used: %lf mins\n", cost);
	system("pause");

	return 0;
}