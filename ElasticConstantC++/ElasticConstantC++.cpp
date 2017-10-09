#include "stdafx.h"

void CreatingCrystall(size_t, double);
void CreatingVacancy(size_t, std::vector<double> &);
double F(double );
double E(size_t, int);
double E(size_t, int, double **);
double Energy(size_t);
double ParameterC();
double SettingLatticeParameter(size_t);
void PrintingMassive(size_t);
void Relaxation(double);
int *NearestToVacancy(size_t, double, double, double, double, size_t &);
void MakingSphere(size_t &, std::vector<double> &, double, int &);

double *X = NULL, *Y = NULL, *Z = NULL;
double ***MyPrettyMassive = NULL;
double **BeforeRelax = NULL;
int counterForRFrom2UpTo3, 
	counterForRFrom5UpTo7, 
	counterForRLessThan5, 
	counterForRLessThan7;


int main()
{
	size_t parameter = 14;
	size_t size = static_cast<size_t> (pow(parameter, 3));
	double	*_X = new double[size],
			*_Y = new double[size],
			*_Z = new double[size];
	X = _X;
	Y = _Y;
	Z = _Z;

	BeforeRelax = new double*[3];
	MyPrettyMassive = new double**[3];
	for (int i = 0; i < 4; i++) MyPrettyMassive[i] = new double*[3];
	
	std::vector<double> vacancy;
	vacancy.push_back(0);
	vacancy.push_back(0);
	vacancy.push_back(0);
	
	double latticeParameter = SettingLatticeParameter(size);

	double energyWithoutDefect = Energy(size);
	CreatingVacancy(size, vacancy);
	double energyWithDefect = Energy(size);
	std::cout << "size = " << size << std::endl;
	MakingSphere(size, vacancy, latticeParameter, counterForRLessThan5);
	std::cout << "size = " << size << std::endl;
	for (int i = 0; i < counterForRFrom2UpTo3; i++)
	{
		std::cout << "i = " << MyPrettyMassive[2][0][i] << " X = " << MyPrettyMassive[2][1][i] << " Y = " << MyPrettyMassive[2][2][i] << " Z = " << MyPrettyMassive[2][3][i] << std::endl;
	}
	Relaxation(latticeParameter);
	std::ofstream out("out.txt");
	double C = 0;
	for (int i = 1; i < 10; i++)
	{
		Relaxation(latticeParameter);
		C = ParameterC();
		std::cout << i << std::endl;
		out << i << "\t" << C << "\n";
	}
	out.close();
	std::cout << energyWithoutDefect << "\t" << energyWithDefect << std::endl;
	system("pause");
    return 0;
}

void CreatingCrystall(size_t size, double latticeParameter)
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

double SettingLatticeParameter(size_t size)
{
	double a0 = 2.845, h = a0 / 20;
	double latticeParameter = a0;
	int central = static_cast<int> (pow(size, 1.0 / 3.0) / 2 + pow(size, 2.0 / 3.0) / 2 + size / 2);
	CreatingCrystall(size, latticeParameter);
	double U0 = E(size, central);
	for (int i = 0; i < 1000; i++)
	{
		latticeParameter += h;
		CreatingCrystall(size, latticeParameter);
		if (E(size, central) == U0) break;
		if (E(size, central) >= U0)
		{
			latticeParameter -= 2 * h;
			CreatingCrystall(size, latticeParameter);
			U0 = E(size, central);
			h /= 8;
		}
		if (E(size, central) < U0) U0 = E(size, central);
	}
	std::cout << "Lattice Parameter = " << latticeParameter << std::endl;
	return latticeParameter;
}

void CreatingVacancy(size_t size, std::vector<double> &vacancy)
{
	int central = static_cast <int> (pow(size, 1.0 / 3.0) / 2 + pow(size, 2.0 / 3.0) / 2 + size / 2);
	vacancy[0] = X[central];
	vacancy[1] = Y[central];
	vacancy[2] = Z[central];
	X[central] = 1000000;
	Y[central] = 1000000;
	Z[central] = 1000000;
}

void PrintingMassive(size_t size)
{
	for (unsigned int i = 0; i < size; i++)	std::cout << "i = " << i << "\tX = " << X[i] << "\tY = " << Y[i] << "\tZ = " << Z[i] << std::endl;
}

void MakingSphere(size_t &size, std::vector<double> &vacancy, double latticeParameter, int &_counterForRLessThan5)
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
			newCounter++;
		}
	}
	counterForRLessThan5 = newCounter;
	
	newCounter = 0;
	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R > latticeParameter * 5 && R <= latticeParameter * pow(size, 1.0 / 3.0) / 2 + 0.001)
		{
			tempX[newCounter] = X[i];
			tempY[newCounter] = Y[i];
			tempZ[newCounter] = Z[i];
			newCounter++;
		}
	}
	counterForRFrom5UpTo7 = newCounter;

	counterForRFrom2UpTo3 = 0;
	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R >= latticeParameter * 2 && R <= latticeParameter * 3) counterForRFrom2UpTo3++;
	}

	counterForRLessThan7 = 0;
	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * 7) counterForRLessThan7++;
	}

	for (int i = 0; i < 4; i++) MyPrettyMassive[0][i] = new double[counterForRLessThan5];
	for (int i = 0; i < 4; i++) MyPrettyMassive[1][i] = new double[counterForRLessThan7];
	for (int i = 0; i < 4; i++) MyPrettyMassive[2][i] = new double[counterForRFrom2UpTo3];
	for (int i = 0; i < 4; i++) MyPrettyMassive[3][i] = new double[counterForRFrom5UpTo7];
	for (int i = 0; i < 4; i++) BeforeRelax[i] = new double[counterForRFrom2UpTo3];

	int newCounter0 = 0, newCounter1 = 0, newCounter2 = 0, newCounter3 = 0;
	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * 5)
		{
			MyPrettyMassive[0][0][newCounter0] = i;
			MyPrettyMassive[0][1][newCounter0] = X[i];
			MyPrettyMassive[0][2][newCounter0] = Y[i];
			MyPrettyMassive[0][3][newCounter0] = Z[i];
			newCounter0++;
		}
		if (R <= latticeParameter * 7)
		{
			MyPrettyMassive[1][0][newCounter1] = i;
			MyPrettyMassive[1][1][newCounter1] = X[i];
			MyPrettyMassive[1][2][newCounter1] = Y[i];
			MyPrettyMassive[1][3][newCounter1] = Z[i];
			newCounter1++;
		}
		if (R >= latticeParameter * 2 && R <= latticeParameter * 3)
		{
			MyPrettyMassive[2][0][newCounter2] = i;
			MyPrettyMassive[2][1][newCounter2] = X[i];
			MyPrettyMassive[2][2][newCounter2] = Y[i];
			MyPrettyMassive[2][3][newCounter2] = Z[i];

			BeforeRelax[0][newCounter2] = i;
			BeforeRelax[1][newCounter2] = X[i];
			BeforeRelax[2][newCounter2] = Y[i];
			BeforeRelax[3][newCounter2] = Z[i];
			newCounter2++;
		}
		if (R > latticeParameter * 5 && R <= latticeParameter * pow(size, 1.0 / 3.0) / 2 + 0.001)
		{
			MyPrettyMassive[3][0][newCounter3] = i;
			MyPrettyMassive[3][1][newCounter3] = X[i];
			MyPrettyMassive[3][2][newCounter3] = Y[i];
			MyPrettyMassive[3][3][newCounter3] = Z[i];
			newCounter3++;
		}
	}

	for (int i = 0; i < counter; i++)
	{
		X[i] = tempX[i];
		Y[i] = tempY[i];
		Z[i] = tempZ[i];
	}

	size = newCounter;
	size_t sizeOfIndex;
	int *index = NearestToVacancy(size, Xcentral, Ycentral, Zcentral, latticeParameter, sizeOfIndex);
}

int *NearestToVacancy(size_t size, double Xcentral, double Ycentral, double Zcentral, double latticeParameter, size_t &sizeOfIndex)
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

void Relaxation(double latticeParameter)
{
	//std::cout << "I`m relaxation" << std::endl;
	double h = 0.005 * latticeParameter;
	for (int i = 0; i < counterForRLessThan5; i++)
	{
		int counter = 0;
		for (unsigned int j = 1; j < 1000; j++)
		{
			double Estart = E(counterForRLessThan5, i, MyPrettyMassive[0]);
			MyPrettyMassive[0][1][i] += h;
			double Ex = E(counterForRLessThan5, i, MyPrettyMassive[0]);
			double dEx = Ex - Estart;
			MyPrettyMassive[0][1][i] -= h;
			MyPrettyMassive[0][2][i] += h;
			double Ey = E(counterForRLessThan5, i, MyPrettyMassive[0]);
			double dEy = Ey - Estart;
			MyPrettyMassive[0][2][i] += h;
			MyPrettyMassive[0][3][i] -= h;
			double Ez = E(counterForRLessThan5, i, MyPrettyMassive[0]);
			double dEz = Ez - Estart;
			MyPrettyMassive[0][3][i] += h;
			double b = sqrt(pow(dEx, 2) + pow(dEy, 2) + pow(dEz, 2));
			if (b == 0) break;
			double hx = -h*dEx / b, hy = -h*dEy / b, hz = -h*dEz / b;
			MyPrettyMassive[0][1][i] += hx;
			MyPrettyMassive[0][2][i] += hy;
			MyPrettyMassive[0][3][i] += hz;
			double dE = E(counterForRLessThan5, i, MyPrettyMassive[0]);
			if (dE >= Estart)
			{
				h /= (j * 8);
				MyPrettyMassive[0][1][i] -= hx;
				MyPrettyMassive[0][2][i] -= hy;
				MyPrettyMassive[0][3][i] -= hz;
				counter++;
			}
			if (counter >= 10) break;
		}
	}
}

double ParameterC()
{
	std::cout << "I`m parameter" << std::endl;
	double Cx = 0, Cy = 0, Cz = 0, C1 = 0;
	for (int i = 0; i < counterForRFrom2UpTo3; i++)
	{
		double R = sqrt(pow((MyPrettyMassive[2][1][i] - BeforeRelax[1][i]), 2) + pow((MyPrettyMassive[2][2][i] - BeforeRelax[2][i]), 2) + pow((MyPrettyMassive[2][3][i] - BeforeRelax[3][i]), 2));
		std::cout << MyPrettyMassive[2][1][i] << "\t" << MyPrettyMassive[2][2][i] << "\t" << MyPrettyMassive[2][3][i] << std::endl;
		if (int(MyPrettyMassive[2][1][i]) != 0) Cx = (MyPrettyMassive[2][1][i] - BeforeRelax[1][i])*pow(R, 3) / MyPrettyMassive[2][1][i];
		if (int(MyPrettyMassive[2][2][i]) != 0) Cy = (MyPrettyMassive[2][2][i] - BeforeRelax[2][i])*pow(R, 3) / MyPrettyMassive[2][2][i];
		if (int(MyPrettyMassive[2][3][i]) != 0) Cz = (MyPrettyMassive[2][3][i] - BeforeRelax[3][i])*pow(R, 3) / MyPrettyMassive[2][3][i];
		C1 += Cx + Cy + Cz;
		std::cout << "C = " << C1 << "\t" << Cx << "\t" << Cy << "\t" << Cz << "\t" << R << std::endl;
	}
	std::cout << "C = " << C1 << std::endl;
	std::cout << "counterForRFrom2UpTo3 = " << counterForRFrom2UpTo3 << std::endl;
	C1 /= counterForRFrom2UpTo3*3;
	std::cout << "C = " << C1 << std::endl;
	return C1;
}

double F(double r)
{
	double F = 0, E = 0.4174, R0 = 2.845, l = 1.3885;
	F = E*(exp(-2 * l*(r - R0)) - 2 * exp(-l*(r - R0)));
	return F;
}

double E(size_t size, int e, double **Massive)
{
	double r = 0, E = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		if (i != e)
		{
			r = sqrt(pow((Massive[1][e] - Massive[1][i]), 2) + pow((Massive[2][e] - Massive[2][i]), 2) + pow((Massive[3][e] - Massive[3][i]), 2));
			E += F(r);
		}
	}
	return E;
}

double E(size_t size, int e)
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

double Energy(size_t size)
{
	double U = 0;
	for (unsigned int i = 0; i < size; i++) U += E(size, i);
	U /= 2;
	return U;
}

