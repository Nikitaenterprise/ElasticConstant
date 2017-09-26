// ElasticConstantC++.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

void CreatingCrystall(double *X, double *Y, double *Z, size_t size, double latticeParameter);
double F(double r);
double E(double *X, double *Y, double *Z, int e);
double Energy(double *X, double *Y, double *Z);
void SettingLatticeParameter(double *X, double *Y, double *Z, size_t size);
void PrintingMassive(double *X, double *Y, double *Z);

int main()
{
	double *X, *Y, *Z;
	size_t size = 10;
	X = new double[static_cast <size_t> (pow(size, 3))];
	Y = new double[static_cast <size_t> (pow(size, 3))];
	Z = new double[static_cast <size_t> (pow(size, 3))];
	std::cout << static_cast <size_t> (pow(size, 3)) << " " << sizeof(X) << std::endl;
	SettingLatticeParameter(X, Y, Z, size);
	PrintingMassive(X, Y, Z);
	system("pause");
    return 0;
}

void CreatingCrystall(double *X, double *Y, double *Z, size_t size, double latticeParameter)
{
	for (unsigned int i = 0; i < size; i++)
	{
		for (unsigned int j = 0; j < size; j++)
		{
			for (unsigned int k = 0; k < size; k++)
			{
				X[i + j + k] = latticeParameter*i;
				Y[i + j + k] = latticeParameter*j;
				Z[i + j + k] = latticeParameter*k;
			}
		}
	}
	PrintingMassive(X, Y, Z);
}

void SettingLatticeParameter(double *X, double *Y, double *Z, size_t size)
{
	double a0 = 2.845, h = a0 / 20;
	double latticeParameter = a0;
	int Nc = static_cast <int> (pow(size,3) / 2 + pow(size,2) / 2 + size / 2);
	CreatingCrystall(X, Y, Z, size, latticeParameter);
	double U0 = E(X, Y, Z, Nc);
	for (int i = 0; i < 1000; i++)
	{
		latticeParameter += h;
		CreatingCrystall(X, Y, Z, size, latticeParameter);
		if (E(X, Y, Z, Nc) == U0) break;
		if (E(X, Y, Z, Nc) >= U0)
		{
			latticeParameter -= 2 * h;
			CreatingCrystall(X, Y, Z, size, latticeParameter);
			U0 = E(X, Y, Z, Nc);
			h /= 8;
		}
		if (E(X, Y, Z, Nc) < U0) U0 = E(X, Y, Z, Nc);
	}

}

void PrintingMassive(double *X, double *Y, double *Z)
{
	for (int i = 0; i < sizeof(X); i++) std::cout << "i = " << i << "	X = " << X[i] << std::endl;	
}

double F(double r)
{
	double F = 0, E = 0.4174, R0 = 2.845, l = 1.3885;
	F = E*(exp(-2 * l*(r - R0)) - 2 * exp(-l*(r - R0)));
	return F;
}

double E(double *X, double *Y, double *Z, int e)
{
	double r = 0, E = 0;
	for (unsigned int i = 0; i < sizeof(X); i++)
	{
		r = sqrt(pow((X[e] - X[i]), 2) + pow((Y[e] - Y[i]), 2) + pow((Z[e] - Z[i]), 2));
		if (i != e)
		{
			E += F(r);
		}
		else
		{
			E = E;
		}
	}
	return E;
}

double Energy(double *X, double *Y, double *Z)
{
	double U = 0;
	for (unsigned int i = 0; i < sizeof(X); i++)
	{
		U += E(X, Y ,Z , i);
	}
	U /= 2;
	return U;
}

