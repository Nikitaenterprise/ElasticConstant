// ElasticConstantC++.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

void CreatingCrystall(double *X, double *Y, double *Z, size_t size, double latticeParameter);
void CreatingVacancy(double *X, double *Y, double *Z, size_t size, std::vector<double> &vacancy);
double F(double r);
double E(double *X, double *Y, double *Z, size_t size, int e);
double Energy(double *X, double *Y, double *Z, size_t size);
double SettingLatticeParameter(double *X, double *Y, double *Z, size_t size);
void PrintingMassive(double *X, double *Y, double *Z, size_t size);
int *NearestToVacancy(double *X, double *Y, double *Z, size_t size, double Xcentral, double Ycentral, double Zcentral, double latticeParameter, size_t &sizeOfIndex);
void MakingSphere(double *X, double *Y, double *Z, size_t &size, std::vector<double> &vacancy, double latticeParameter);

int main()
{
	size_t parameter = 14;
	size_t size = static_cast<size_t> (pow(parameter, 3));
	double	*X = new double[size],
			*Y = new double[size],
			*Z = new double[size];

	std::vector<double> vacancy;
	vacancy.push_back(0);
	vacancy.push_back(0);
	vacancy.push_back(0);
	
	double latticeParameter = SettingLatticeParameter(X, Y, Z, size);

	double energyWithoutDefect = Energy(X, Y, Z, size);
	CreatingVacancy(X, Y, Z, size, vacancy);
	double energyWithDefect = Energy(X, Y, Z, size);

	std::cout << "size = " << size << std::endl;
	MakingSphere(X, Y, Z, size, vacancy, latticeParameter);
	std::cout << "size = " << size << std::endl;
	
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

double SettingLatticeParameter(double *X, double *Y, double *Z, size_t size)
{
	double a0 = 2.845, h = a0 / 20;
	double latticeParameter = a0;
	int central = static_cast<int> (pow(size, 1.0 / 3.0) / 2 + pow(size, 2.0 / 3.0) / 2 + size / 2);
	CreatingCrystall(X, Y, Z, size, latticeParameter);
	double U0 = E(X, Y, Z, size, central);
	for (int i = 0; i < 1000; i++)
	{
		latticeParameter += h;
		CreatingCrystall(X, Y, Z, size, latticeParameter);
		if (E(X, Y, Z, size, central) == U0) break;
		if (E(X, Y, Z, size, central) >= U0)
		{
			latticeParameter -= 2 * h;
			CreatingCrystall(X, Y, Z, size, latticeParameter);
			U0 = E(X, Y, Z, size, central);
			h /= 8;
		}
		if (E(X, Y, Z, size, central) < U0) U0 = E(X, Y, Z, size, central);
	}
	return latticeParameter;
}

void CreatingVacancy(double *X, double *Y, double *Z, size_t size, std::vector<double> &vacancy)
{
	int central = static_cast <int> (pow(size, 1.0 / 3.0) / 2 + pow(size, 2.0 / 3.0) / 2 + size / 2);
	vacancy[0] = X[central];
	vacancy[1] = Y[central];
	vacancy[2] = Z[central];
	X[central] = 1000000;
	Y[central] = 1000000;
	Z[central] = 1000000;
}

void PrintingMassive(double *X, double *Y, double *Z, size_t size)
{
	for (unsigned int i = 0; i < size; i++)	std::cout << "i = " << i << "\tX = " << X[i] << "\tY = " << Y[i] << "\tZ = " << Z[i] << std::endl;
}

void MakingSphere(double *X, double *Y, double *Z, size_t &size, std::vector<double> &vacancy, double latticeParameter)
{
	double	*tempX = new double[size],
			*tempY = new double[size],
			*tempZ = new double[size];
	double Xcentral = vacancy[0], Ycentral = vacancy[1], Zcentral = vacancy[2];
	std::cout << Xcentral << "\t" << Ycentral << "\t" << Zcentral << std::endl;
	int counter = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * pow(size, 1.0 / 3.0) / 2 + 0.001)
		{
			tempX[counter] = X[i];
			tempY[counter] = Y[i];
			tempZ[counter] = Z[i];
			counter += 1;
		}
	}

	delete[] X, Y, Z;
	X = new double[counter];
	Y = new double[counter];
	Z = new double[counter];

	for (int i = 0; i < counter; i++)
	{
		X[i] = tempX[i];
		Y[i] = tempY[i];
		Z[i] = tempZ[i];
	}
	
	delete[] tempX, tempY, tempZ;
	tempX = new double[counter];
	tempY = new double[counter];
	tempZ = new double[counter];

	int newCounter = 0;

	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * 5)
		{
			tempX[newCounter] = X[i];
			tempY[newCounter] = Y[i];
			tempZ[newCounter] = Z[i];
			newCounter += 1;
		}
	}
	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * pow(size, 1.0 / 3.0) / 2 + 0.001 && R > latticeParameter * 5)
		{
			tempX[newCounter] = X[i];
			tempY[newCounter] = Y[i];
			tempZ[newCounter] = Z[i];
			newCounter += 1;
		}
	}
	std::cout << "counters = " << counter << "\t" << newCounter << std::endl;
	delete[] tempX, tempY, tempZ;
	size = newCounter;
	size_t sizeOfIndex;
	int *index = NearestToVacancy(X, Y, Z, size, Xcentral, Ycentral, Zcentral, latticeParameter, sizeOfIndex);
}

int *NearestToVacancy(double *X, double *Y, double *Z, size_t size, double Xcentral, double Ycentral, double Zcentral, double latticeParameter, size_t &sizeOfIndex)
{
	int counter = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R >= latticeParameter * 2 && R <= latticeParameter * 3)	counter += 1;
	}
	int *index = new int[counter];
	counter = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R >= latticeParameter * 2 && R <= latticeParameter * 3)
		{
			index[counter] = i;
			counter += 1;
		}
	}
	sizeOfIndex = counter;
	return index;
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
	}
	return E;
}

double Energy(double *X, double *Y, double *Z, size_t size)
{
	double U = 0;
	for (unsigned int i = 0; i < size; i++) U += E(X, Y, Z, size, i);
	U /= 2;
	return U;
}

