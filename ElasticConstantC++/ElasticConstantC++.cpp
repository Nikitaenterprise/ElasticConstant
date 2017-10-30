#include "stdafx.h"

void CreatingCrystall(double);
void CreatingVacancy(size_t, std::vector<double> &);
double F(double );
double Fsharp(double );
double E(size_t, int);
double E(size_t, int, double **);
double EForAt(size_t, int, double **);
double Energy(size_t);
double ParameterC();
double SettingLatticeParameter(size_t);
void PrintingMassive(size_t);
void Relaxation();
void Relaxation(int, double);
void MakingSphere(size_t &, std::vector<double> &);
void RelaxationOuterSphere(const std::vector<double> &);
void CopyData(double **, int, double **, int);


double *X = NULL, *Y = NULL, *Z = NULL;
double ***MyPrettyMassive = NULL;
std::vector<std::vector<double>> vecNearestTo;
std::vector<double> vacancy;
double **NearestTo = NULL;
double **BeforeRelax = NULL;
size_t parameter = 16;
long double latticeParameter;
int counterForRFrom2UpTo3, 
	counterForRFrom5UpTo7, 
	counterForRLessThan5, 
	counterForRLessThan7;
double C = 0;

int main()
{
	int startClock = clock();
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
	
	vacancy.push_back(0);
	vacancy.push_back(0);
	vacancy.push_back(0);
	
	latticeParameter = SettingLatticeParameter(size);

	double energyWithoutDefect = Energy(size);
	CreatingVacancy(size, vacancy);
	double energyWithDefect = Energy(size);
	MakingSphere(size, vacancy);

	double VfE1 = 0;
	for (int i = 0; i < counterForRLessThan7; i++)
	{
		for (int j = 0; j < counterForRLessThan7; j++)
		{
			if (i != j)
			{
				double R = sqrt(pow((MyPrettyMassive[1][1][i] - MyPrettyMassive[1][1][j]), 2) + pow((MyPrettyMassive[1][2][i] - MyPrettyMassive[1][2][j]), 2) + pow((MyPrettyMassive[1][3][i] - MyPrettyMassive[1][3][j]), 2));
				VfE1 += R*Fsharp(R);
			}
		}
	}
	VfE1 /= 2;
	std::ofstream out0("nonrelaxed.txt");
	for (int i = 0; i < counterForRLessThan7; i++)
	{
		out0 << MyPrettyMassive[1][1][i] << "\t" << MyPrettyMassive[1][2][i] << "\t" << MyPrettyMassive[1][3][i] << std::endl;
	}
	out0.close();

	std::ofstream out("out.txt");
	if (C != 0)
	{
		RelaxationOuterSphere(vacancy);
		Relaxation();
		//C = ParameterC();
		std::cout << "I`m in main alyo bliat\n";
		out << 0 << "\t" << C << "\n";
	}
	for (int i = 1; i < 10; i++)
	{
		RelaxationOuterSphere(vacancy);
		Relaxation();
		C = ParameterC();
		std::cout << "\t" << i << std::endl;
		out << i << "\t" << C << "\n";
	}
	out.close();

	double dVf = 4 * 3.14*C;
	double VfE2 = 0;
	std::ofstream out1("relaxed.txt");
	for (int i = 0; i < counterForRLessThan5; i++)
	{
		out1 << MyPrettyMassive[0][1][i] << "\t" << MyPrettyMassive[0][2][i] << "\t" << MyPrettyMassive[0][3][i] << std::endl;
	}
	out1.close();
	for (int i = 0; i < counterForRLessThan7; i++)
	{
		for (int j = 0; j < counterForRLessThan7; j++)
		{
			if (i != j)
			{
				double R = sqrt(pow((MyPrettyMassive[1][1][i] - MyPrettyMassive[1][1][j]), 2) + pow((MyPrettyMassive[1][2][i] - MyPrettyMassive[1][2][j]), 2) + pow((MyPrettyMassive[1][3][i] - MyPrettyMassive[1][3][j]), 2));
				VfE2 += R*Fsharp(R);
			}
		}
	}
	VfE2 /= 2;

	double VfE = VfE2 - VfE1;
	VfE /= -6;

	double Vf = VfE + dVf;
	double Vf1 = Vf / (4 / 3 * 3.14*pow(parameter / 2, 3));

	std::cout << "VfE1 = " << VfE1 << " VfE2 = " << VfE2 << std::endl;
	std::cout << "VfE2 - VfE1 = " << VfE2 - VfE1 << std::endl;
	std::cout << "VfE = " << VfE << std::endl;
	std::cout << "dVf = " << dVf << std::endl;
	std::cout << "Vf = " << Vf << std::endl;
	std::cout << "Vf1 = " << Vf1 << std::endl;
	std::cout << energyWithoutDefect << "\t" << energyWithDefect << std::endl;
	
	int index;
	double tempR = 100;
	for (int i = 0; i < counterForRLessThan5; i++)
	{
		double R = sqrt(pow(vacancy[0] - MyPrettyMassive[0][1][i], 2) + pow(vacancy[1] - MyPrettyMassive[0][2][i], 2) + pow(vacancy[2] - MyPrettyMassive[0][3][i], 2));
		if (R < tempR && R != 0)
		{
			tempR = R;
			index = i;
			std::cout << "tempR = " << tempR << " index = " << index << std::endl;
		}
	}
	std::cout << "index is = " << index << std::endl;
	std::ofstream out2("traj.txt");
	double dx = abs(vacancy[0] - MyPrettyMassive[0][1][index]) / 30;
	dx = latticeParameter / 30;
	if (vacancy[0] > MyPrettyMassive[0][1][index]) dx *= 1;
	else if (vacancy[0] < MyPrettyMassive[0][1][index]) dx*=-1;
	std::cout << "vac = " << vacancy[0] << " x = " << MyPrettyMassive[0][1][index] << std::endl;
	for (int i = 0; i < 60; i++)
	{
		std::cout << i << std::endl;
		Relaxation(index, dx);
		double En = E(counterForRLessThan5, index, MyPrettyMassive[0]);
		out2 << En << "\t" << MyPrettyMassive[0][1][index] << "\t" << MyPrettyMassive[0][2][index] << "\t" << MyPrettyMassive[0][3][index] << std::endl;
	}
	out2.close();
	std::cout << "\a" << std::endl;
	int endClock = clock();
	std::cout << "Time of calculation = " << (endClock - startClock)/CLOCKS_PER_SEC << std::endl;
	//std::cout << "\a" << std::endl;
	system("pause");
    return 0;
}

void CreatingCrystall(double latticeParameter)
{
	int counter = 0;
	for (unsigned int i = 0; i < parameter; i++)
	{
		for (unsigned int j = 0; j < parameter; j++)
		{
			for (unsigned int k = 0; k < parameter; k++)
			{
				X[counter] = latticeParameter*i;
				Y[counter] = latticeParameter*j;
				Z[counter] = latticeParameter*k;
				counter ++;
				/*X[counter] = latticeParameter*i + latticeParameter / 2;
				Y[counter] = latticeParameter*j + latticeParameter / 2;
				Z[counter] = latticeParameter*k + latticeParameter / 2;
				counter++;*/
			}
		}
	}
}

double SettingLatticeParameter(size_t size)
{
	double a0 = 2.845, h = a0 / 20;
	latticeParameter = a0;
	int central = static_cast <int> (pow(parameter, 3) / 2 + pow(parameter, 2) / 2 + parameter / 2);
	std::cout << X[central] << "\t" << Y[central] << "\t" << Z[central] << std::endl;
	CreatingCrystall(latticeParameter);
	double U0 = E(size, central);
	for (int i = 0; i < 1000000; i++)
	{
		latticeParameter += h;
		CreatingCrystall(latticeParameter);
		if (E(size, central) == U0) break;
		if (E(size, central) >= U0)
		{
			latticeParameter -= 8 * h;
			CreatingCrystall(latticeParameter);
			U0 = E(size, central);
			h /= 10;
		}
		if (E(size, central) < U0) U0 = E(size, central);
	}
	std::cout << "Lattice Parameter = " << latticeParameter << std::endl;
	return latticeParameter;
}

void CreatingVacancy(size_t size, std::vector<double> &vacancy)
{
	int central = static_cast <int> (pow(parameter, 3) / 2 + pow(parameter, 2) / 2 + parameter / 2);
	std::cout << central << std::endl;
	std::cout << X[central] << "\t" << Y[central] << "\t" << Z[central] << std::endl;
	vacancy[0] = X[central];
	vacancy[1] = Y[central];
	vacancy[2] = Z[central];
	X[central] = 1000000;
	Y[central] = 1000000;
	Z[central] = 1000000;
}

void MakingSphere(size_t &size, std::vector<double> &vacancy)
{
	double	*tempX = new double[size],
			*tempY = new double[size],
			*tempZ = new double[size];
	double Xcentral = vacancy[0], Ycentral = vacancy[1], Zcentral = vacancy[2];
	int counter = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * parameter / 2)
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
	
	/*::VfE1 = 0;
	for (int i = 0 ; i < counter; i++)
	{
		for (int j = 0; j < counter; j++)
		{
			double R = sqrt(pow(X[i] - X[j], 2) + pow(Y[i] - Y[j], 2) + pow(Y[i] - Y[j], 2));
			::VfE1 += R*Fsharp(R);
		}
	}
	::VfE1 /= 2;
	std::cout << "VfE1 = " << ::VfE1 << std::endl;*/

	delete[] tempX, tempY, tempZ;
	tempX = new double[counter];
	tempY = new double[counter];
	tempZ = new double[counter];

	int newCounter = 0;
	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * (parameter / 2 - 5))
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
		if (R > latticeParameter * (parameter / 2 - 5) && R <= latticeParameter * parameter / 2)
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
		if (R >= latticeParameter * ((parameter / 2 - 5) / 2 - 1) && R <= latticeParameter * ((parameter / 2 - 5) / 2 + 1)) counterForRFrom2UpTo3++;
	}

	counterForRLessThan7 = 0;
	for (int i = 0; i < counter; i++)
	{
		double R = sqrt(pow((Xcentral - X[i]), 2) + pow((Ycentral - Y[i]), 2) + pow((Zcentral - Z[i]), 2));
		if (R <= latticeParameter * parameter / 2) counterForRLessThan7++;
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
		if (R <= latticeParameter * (parameter / 2 - 5))
		{
			MyPrettyMassive[0][0][newCounter0] = i;
			MyPrettyMassive[0][1][newCounter0] = X[i];
			MyPrettyMassive[0][2][newCounter0] = Y[i];
			MyPrettyMassive[0][3][newCounter0] = Z[i];
			newCounter0++;
		}
		if (R <= latticeParameter * parameter / 2)
		{
			MyPrettyMassive[1][0][newCounter1] = i;
			MyPrettyMassive[1][1][newCounter1] = X[i];
			MyPrettyMassive[1][2][newCounter1] = Y[i];
			MyPrettyMassive[1][3][newCounter1] = Z[i];
			newCounter1++;
		}
		if (R >= latticeParameter * ((parameter / 2 - 5) / 2 - 1) && R <= latticeParameter * ((parameter / 2 - 5) / 2 + 1))
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
		if (R > latticeParameter * ((parameter / 2) - 5) && R <= latticeParameter * parameter / 2)
		{
			MyPrettyMassive[3][0][newCounter3] = i;
			MyPrettyMassive[3][1][newCounter3] = X[i];
			MyPrettyMassive[3][2][newCounter3] = Y[i];
			MyPrettyMassive[3][3][newCounter3] = Z[i];
			newCounter3++;
		}
	}

	NearestTo = new double*[counterForRLessThan5];
	for (int i = 0; i < counterForRLessThan5; i++)
	{
		int counter = 0;
		for (int j = 0; j < counterForRLessThan5; j++)
		{
			double R = sqrt(pow((MyPrettyMassive[0][1][i] - MyPrettyMassive[0][1][j]), 2)
						+ pow((MyPrettyMassive[0][2][i] - MyPrettyMassive[0][2][j]), 2)
						+ pow((MyPrettyMassive[0][3][i] - MyPrettyMassive[0][3][i]), 2));

			if (R<=5) counter++;
		}
		vecNearestTo.push_back(std::vector<double>());
		vecNearestTo[i].push_back(MyPrettyMassive[0][0][i]);
		//std::cout << counter << std::endl;
		NearestTo[i] = new double[counter];
		counter = 0;
		for (int j = 0; j < counterForRLessThan5; j++)
		{
			double R = sqrt(pow((MyPrettyMassive[0][1][i] - MyPrettyMassive[0][1][j]), 2)
				+ pow((MyPrettyMassive[0][2][i] - MyPrettyMassive[0][2][j]), 2)
				+ pow((MyPrettyMassive[0][3][i] - MyPrettyMassive[0][3][i]), 2));

			if (R <= 5)
			{
				vecNearestTo[i].push_back(MyPrettyMassive[0][0][j]);
				NearestTo[i][counter] = MyPrettyMassive[0][0][j];
				counter++;
			}
		}
	}

	for (int i = 0; i < vecNearestTo.size(); i++)
	{
		//for (int j = 0; j < vecNearestTo[i].size(); j++)
		//{
			std::cout << vecNearestTo[i][0] << "\t";
		//}
		std::cout << std::endl;
	}

	for (int i = 0; i < counter; i++)
	{
		X[i] = tempX[i];
		Y[i] = tempY[i];
		Z[i] = tempZ[i];
	}
	size = counter;
}

void Relaxation()
{
	std::cout << "C from relax = " << C << std::endl;
	double *tempX = new double[counterForRFrom2UpTo3], *tempY = new double[counterForRFrom2UpTo3], *tempZ = new double[counterForRFrom2UpTo3];
	for (int i = 0; i < counterForRFrom2UpTo3; i++)
	{
		tempX[i] = MyPrettyMassive[2][1][i];
		tempY[i] = MyPrettyMassive[2][2][i];
		tempZ[i] = MyPrettyMassive[2][3][i];
	}
	//std::cout << "I`m relaxation" << std::endl;
	for (int j = 0; j < 20; j++)
	{
		std::cout << j << std::endl;
		//std::cout << MyPrettyMassive[0][1][20] << "\t" << MyPrettyMassive[0][2][20] << "\t" << MyPrettyMassive[0][3][20] << std::endl;
		for (int i = 0; i < counterForRLessThan5; i++)
		{
			double h = latticeParameter / 100;
			int counter = 1;
			while (counter < 15)
			{
				double Estart = E(counterForRLessThan5, i, MyPrettyMassive[0]);
				//std::cout << "Estart = " << Estart << std::endl;
				MyPrettyMassive[0][1][i] += h;
				double Ex = E(counterForRLessThan5, i, MyPrettyMassive[0]);
				double dEx = Ex - Estart;
				MyPrettyMassive[0][1][i] -= h;
				MyPrettyMassive[0][2][i] += h;
				double Ey = E(counterForRLessThan5, i, MyPrettyMassive[0]);
				double dEy = Ey - Estart;
				MyPrettyMassive[0][2][i] -= h;
				MyPrettyMassive[0][3][i] += h;
				double Ez = E(counterForRLessThan5, i, MyPrettyMassive[0]);
				double dEz = Ez - Estart;
				MyPrettyMassive[0][3][i] -= h;
				double b = sqrt(pow(dEx, 2) + pow(dEy, 2) + pow(dEz, 2));
				if (b == 0) break;
				double hx = -h*dEx / b, hy = -h*dEy / b, hz = -h*dEz / b;
				MyPrettyMassive[0][1][i] += hx;
				MyPrettyMassive[0][2][i] += hy;
				MyPrettyMassive[0][3][i] += hz;
				double dE = E(counterForRLessThan5, i, MyPrettyMassive[0]);
				if (dE >= Estart)
				{
					h /= (counter * 8);
					MyPrettyMassive[0][1][i] -= hx;
					MyPrettyMassive[0][2][i] -= hy;
					MyPrettyMassive[0][3][i] -= hz;
					counter++;
				}
			}
		}
	}
	CopyData(MyPrettyMassive[2], counterForRFrom2UpTo3, MyPrettyMassive[0], counterForRLessThan5);
	/*for (int i = 0; i < counterForRFrom2UpTo3; i++)
	{
		std::cout << "i = " << MyPrettyMassive[2][0][i] << " dX = " << MyPrettyMassive[2][1][i] - tempX[i] << " dY = " << MyPrettyMassive[2][2][i] - tempY[i] << " dZ = " << MyPrettyMassive[2][3][i] - tempZ[i] << std::endl;
	}*/
}

void Relaxation(int a, double dx)
{
	/*double *tempX = new double[counterForRFrom2UpTo3], *tempY = new double[counterForRFrom2UpTo3], *tempZ = new double[counterForRFrom2UpTo3];
	for (int i = 0; i < counterForRFrom2UpTo3; i++)
	{
		tempX[i] = MyPrettyMassive[2][1][i];
		tempY[i] = MyPrettyMassive[2][2][i];
		tempZ[i] = MyPrettyMassive[2][3][i];
	}*/
	MyPrettyMassive[0][1][a] += dx;
	for (int j = 0; j < 10; j++)
	{
		double h = latticeParameter / 100;
		int counter = 1;
		while (counter < 10)
		{
			double Estart = E(counterForRLessThan5, a, MyPrettyMassive[0]);
			MyPrettyMassive[0][2][a] += h;
			double Ey = E(counterForRLessThan5, a, MyPrettyMassive[0]);
			double dEy = Ey - Estart;
			MyPrettyMassive[0][2][a] -= h;
			MyPrettyMassive[0][3][a] += h;
			double Ez = E(counterForRLessThan5, a, MyPrettyMassive[0]);
			double dEz = Ez - Estart;
			MyPrettyMassive[0][3][a] -= h;
			double b = sqrt(pow(dEy, 2) + pow(dEz, 2))+0.0000001;
			//std::cout << "Estart = " << Estart << " Ey = " << Ey << " Ez = " << Ez << std::endl;
			//std::cout << "b = " << b << std::endl;
			if (b == 0) break;
			double hy = -h*dEy / b, hz = -h*dEz / b;
			MyPrettyMassive[0][2][a] += hy;
			MyPrettyMassive[0][3][a] += hz;
			double dE = E(counterForRLessThan5, a, MyPrettyMassive[0]);
			if (dE >= Estart)
			{
				h /= (counter * 8);
				MyPrettyMassive[0][2][a] -= hy;
				MyPrettyMassive[0][3][a] -= hz;
				counter++;
			}
			//std::cout << "hy = " << hy << " hz = " << hz << std::endl;
		}
		for (int i = 0; i < counterForRLessThan5; i++)
		{
			if (i != a)
			{
				h = latticeParameter / 100;
				counter = 1;
				while (counter < 10)
				{
					double Estart = E(counterForRLessThan5, i, MyPrettyMassive[0]);
					MyPrettyMassive[0][1][i] += h;
					double Ex = E(counterForRLessThan5, i, MyPrettyMassive[0]);
					double dEx = Ex - Estart;
					MyPrettyMassive[0][1][i] -= h;
					MyPrettyMassive[0][2][i] += h;
					double Ey = E(counterForRLessThan5, i, MyPrettyMassive[0]);
					double dEy = Ey - Estart;
					MyPrettyMassive[0][2][i] -= h;
					MyPrettyMassive[0][3][i] += h;
					double Ez = E(counterForRLessThan5, i, MyPrettyMassive[0]);
					double dEz = Ez - Estart;
					MyPrettyMassive[0][3][i] -= h;
					double b = sqrt(pow(dEx, 2) + pow(dEy, 2) + pow(dEz, 2));
					if (b == 0) break;
					double hx = -h*dEx / b, hy = -h*dEy / b, hz = -h*dEz / b;
					MyPrettyMassive[0][1][i] += hx;
					MyPrettyMassive[0][2][i] += hy;
					MyPrettyMassive[0][3][i] += hz;
					double dE = E(counterForRLessThan5, i, MyPrettyMassive[0]);
					if (dE >= Estart)
					{
						h /= (counter * 8);
						MyPrettyMassive[0][1][i] -= hx;
						MyPrettyMassive[0][2][i] -= hy;
						MyPrettyMassive[0][3][i] -= hz;
						counter++;
					}
				}
			}
		}
	}
	CopyData(MyPrettyMassive[2], counterForRFrom2UpTo3, MyPrettyMassive[0], counterForRLessThan5);
	CopyData(MyPrettyMassive[1], counterForRLessThan7, MyPrettyMassive[0], counterForRLessThan5);
	CopyData(MyPrettyMassive[3], counterForRFrom5UpTo7, MyPrettyMassive[1], counterForRLessThan7);
}

double ParameterC()
{
	std::cout << "I`m parameter" << std::endl;
	std::cout << "C = " << C << std::endl;
	double Cx = 0, Cy = 0, Cz = 0, C1 = 0;
	int counter = 3;
	for (int i = 0; i < counterForRFrom2UpTo3; i++)
	{
		double R = sqrt(pow((MyPrettyMassive[2][1][i] - BeforeRelax[1][i]), 2) + pow((MyPrettyMassive[2][2][i] - BeforeRelax[2][i]), 2) + pow((MyPrettyMassive[2][3][i] - BeforeRelax[3][i]), 2));
		if (int(MyPrettyMassive[2][1][i]) != 0) Cx = (MyPrettyMassive[2][1][i] - BeforeRelax[1][i])*pow(R, 3) / MyPrettyMassive[2][1][i];
		else counter--;
		if (int(MyPrettyMassive[2][2][i]) != 0) Cy = (MyPrettyMassive[2][2][i] - BeforeRelax[2][i])*pow(R, 3) / MyPrettyMassive[2][2][i];
		else counter--;
		if (int(MyPrettyMassive[2][3][i]) != 0) Cz = (MyPrettyMassive[2][3][i] - BeforeRelax[3][i])*pow(R, 3) / MyPrettyMassive[2][3][i];
		else counter--;
		C1 += Cx + Cy + Cz;
		//std::cout << "C1 = " << C1 << std::endl;
	}
	C1 /= (counterForRFrom2UpTo3*counter);
	//std::cout << "C1 = " << C1 << std::endl;
	CopyData(MyPrettyMassive[0], counterForRLessThan5, MyPrettyMassive[2], counterForRFrom2UpTo3);
	CopyData(MyPrettyMassive[1], counterForRLessThan7, MyPrettyMassive[0], counterForRLessThan5);
	CopyData(MyPrettyMassive[3], counterForRFrom5UpTo7, MyPrettyMassive[1], counterForRLessThan7);
	return C1;
}

void CopyData(double **MassiveTo, int sizeTo, double **MassiveFrom, int sizeFrom)
{
	//int counter = 0;
	for (int i = 0; i < sizeFrom; i++)
	{
		for (int j = 0; j < sizeTo; j++)
		{
			if (MassiveFrom[0][i] == MassiveTo[0][j])
			{
				MassiveTo[1][j] = MassiveFrom[1][i];
				MassiveTo[2][j] = MassiveFrom[2][i];
				MassiveTo[3][j] = MassiveFrom[3][i];
				//counter++;
			}
		}
	}
	//std::cout << "counter = " << counter << " sizeFrom = " << sizeFrom << " sizeTo = " << sizeTo << std::endl;
}

void RelaxationOuterSphere(const std::vector<double> &vacancy)
{
	std::cout << "I`m outer sphere" << std::endl;
	std::cout << "C from outer shphere = " << C << std::endl;
	double R, dX, dY, dZ;
	for (int i = 0; i < counterForRFrom5UpTo7; i++)
	{
		R = sqrt(pow((MyPrettyMassive[1][1][i] - vacancy[0]), 2) + pow((MyPrettyMassive[1][2][i] - vacancy[1]), 2) + pow((MyPrettyMassive[1][3][i] - vacancy[2]), 2));
		dX = C*MyPrettyMassive[1][1][i] / (pow(R, 3));
		dY = C*MyPrettyMassive[1][2][i] / (pow(R, 3));
		dZ = C*MyPrettyMassive[1][3][i] / (pow(R, 3));
		//std::cout << "dX = " << dX << " dY = " << dY << " dZ = " << dZ << std::endl;
		MyPrettyMassive[1][1][i] += dX;
		MyPrettyMassive[1][2][i] += dY;
		MyPrettyMassive[1][3][i] += dZ;
	}
}

double F(double r)
{
	double F, E = 0.4174, R0 = 2.845, l = 1.3885;
	F = E*(exp(-2 * l*(r - R0)) - 2 * exp(-l*(r - R0)));
	return F;
}

double Fsharp(double r)
{
	double Fsharp, E = 0.4174, R0 = 2.845, l = 1.3885;
	Fsharp = E*(-2 * l*exp(-2 * l*(r - R0)) + 2 * l*exp(-l*(r - R0)));
	return Fsharp;
}

double E(size_t size, int e, double **Massive)
{
	double R = 0, E = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		if (i != e)
		{
			R = sqrt(pow((Massive[1][e] - Massive[1][i]), 2) + pow((Massive[2][e] - Massive[2][i]), 2) + pow((Massive[3][e] - Massive[3][i]), 2));
			E += F(R);
		}
	}
	return E;
}

double EForAt(size_t size, int index, double **Massive)
{
	double R = 0, E = 0;
	for (unsigned int i = 0; i < size; i++)
	{
		if (Massive[0][i] != Massive[0][index])
		{
			//std::cout << "index1 = " << Massive[0][i] << " index2 = " << Massive[0][index] << std::endl;
			R = sqrt(pow((Massive[1][index] - Massive[1][i]), 2) + pow((Massive[2][index] - Massive[2][i]), 2) + pow((Massive[3][index] - Massive[3][i]), 2));
			E += F(R);
			//std::cout << "R = " << R << " E = " << E << std::endl;
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

