// ElasticConstantC++.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

void CreatingCrystall(double *X, double *Y, double *Z, size_t size, double latticeParameter);
void CreatingVacancy(double *X, double *Y, double *Z, size_t size);
double F(double r);
double E(double *X, double *Y, double *Z, size_t size, int e);
double Energy(double *X, double *Y, double *Z, size_t size);
void SettingLatticeParameter(double *X, double *Y, double *Z, size_t size);
void PrintingMassive(double *X, double *Y, double *Z, size_t size);

int main()
{
	size_t parameter = 10;
	size_t size = static_cast <size_t> (pow(parameter, 3));
	double	*X = new double[size],
			*Y = new double[size],
			*Z = new double[size];

	//std::vector <double> arr1;
	//std::vector <double*> arr2;

	//std::vector <std::vector<double>> arr3;
	//arr3.push_back(std::vector <double>());
	//arr3[arr3.size() - 1].push_back(100);


	SettingLatticeParameter(X, Y, Z, size);

	double energyWithoutDefect = Energy(X, Y, Z, size);
	CreatingVacancy(X, Y, Z, size);
	double energyWithDefect = Energy(X, Y, Z, size);

	std::cout << energyWithoutDefect << "\t" << energyWithDefect << std::endl;
	system("pause");
    return 0;
}

void CreatingCrystall(double *X, double *Y, double *Z, size_t size, double latticeParameter)
{
	int index = 0;
	for (unsigned int i = 0; i < (pow(size, 1.0 / 3.0)); i++)
	{
		for (unsigned int j = 0; j < (pow(size, 1.0 / 3.0)); j++)
		{
			for (unsigned int k = 0; k < (pow(size, 1.0 / 3.0)); k++)
			{
				X[index] = latticeParameter*i;
				Y[index] = latticeParameter*j;
				Z[index] = latticeParameter*k;
				index += 1;
			}
		}
	}
}

void SettingLatticeParameter(double *X, double *Y, double *Z, size_t size)
{
	double a0 = 2.845, h = a0 / 20;
	double latticeParameter = a0;
	int Nc = static_cast <int> (pow(size, 1.0 / 3.0) / 2 + pow(size, 2.0 / 3.0) / 2 + size / 2);
	CreatingCrystall(X, Y, Z, size, latticeParameter);
	double U0 = E(X, Y, Z, size, Nc);
	for (int i = 0; i < 1000; i++)
	{
		latticeParameter += h;
		CreatingCrystall(X, Y, Z, size, latticeParameter);
		if (E(X, Y, Z, size, Nc) == U0) break;
		if (E(X, Y, Z, size, Nc) >= U0)
		{
			latticeParameter -= 2 * h;
			CreatingCrystall(X, Y, Z, size, latticeParameter);
			U0 = E(X, Y, Z, size, Nc);
			h /= 8;
		}
		if (E(X, Y, Z, size, Nc) < U0) U0 = E(X, Y, Z, size, Nc);
	}

}

void CreatingVacancy(double *X, double *Y, double *Z, size_t size)
{
	int Nc = static_cast <int> (pow(size, 1.0 / 3.0) / 2 + pow(size, 2.0 / 3.0) / 2 + size / 2);
	X[Nc] = 1000000;
	Y[Nc] = 1000000;
	Z[Nc] = 1000000;
}

void PrintingMassive(double *X, double *Y, double *Z, size_t size)
{
	for (unsigned int i = 0; i < size; i++) std::cout << "i = " << i << "	X = " << X[i] << std::endl;
}

double F(double r)
{
	double F = 0, E = 0.4174, R0 = 2.845, l = 1.3885;
	F = E*(exp(-2 * l*(r - R0)) - 2 * exp(-l*(r - R0)));
	return F;
}

double E(double *X, double *Y, double *Z, size_t size, int e)
{
	double r = 0, E = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		if (i != e)
		{
			r = sqrt(pow((X[e] - X[i]), 2) + pow((Y[e] - Y[i]), 2) + pow((Z[e] - Z[i]), 2));
			E += F(r);
		}
		else
		{
			E = E;
		}
	}
	return E;
}

double Energy(double *X, double *Y, double *Z, size_t size)
{
	double U = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		U += E(X, Y, Z, size, i);
	}
	U /= 2;
	return U;
}

