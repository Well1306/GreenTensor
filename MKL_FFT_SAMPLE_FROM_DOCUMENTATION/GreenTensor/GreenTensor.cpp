
#include <iostream>
#include "Hankel-transform.cpp"
#include <stdlib.h>
#include <complex>

const double lambda_w = 500;
const double mu0 = 4 * PI * pow(10, -7);
const std::complex<double> eps0 = (pow(10, 7) / (4 * PI), 0);
const std::complex<double> eps_up = 1.4; //esp в верхнем полупространстве
const std::complex<double> eps1 = (-0.33564, -2.9236);
const std::complex<double> eps2 = 1.2;
const std::complex<double> eps3 = 1.5;
const std::complex<double> epsp1 = 2.0;
const std::complex<double> epsp2 = (-0.33564, -2.9236);
std::complex<double> k0 = 2 * PI / lambda_w; //k в верхнем полупространстве
std::complex<double> k1 = k0 * sqrt(eps1 / eps_up);
std::complex<double> k2 = k0 * sqrt(eps2 / eps_up);
std::complex<double> k3 = k0 * sqrt(eps3 / eps_up);
std::complex<double> kp1 = k0 * sqrt(epsp1 / eps_up);
std::complex<double> kp2 = k0 * sqrt(epsp2 / eps_up);

double z0 = -2 * lambda_w;
double z1 = -3.5 * lambda_w;
double z2 = -5 * lambda_w;



std::complex<double> eta(std::string m, double l)
{
	std::complex<double> r;
	if (m == "0")
		r = sqrt(pow(l, 2) - pow(k0, 2));
	else if (m == "1")
		r = sqrt(pow(l, 2) - pow(k1, 2));
	else if(m == "2")
		r = sqrt(pow(l, 2) - pow(k2, 2));
	else if(m =="3")
		r = sqrt(pow(l, 2) - pow(k3, 2));
	else if (m == "p1")
		r = sqrt(pow(l, 2) - pow(kp1, 2));
	else if (m == "p2")
		r = sqrt(pow(l, 2) - pow(kp2, 2));
	return r;
}


std::complex<double> r2(double l)
{
	return exp(eta("2", l) * z2) * (eta("2", l) - eta("3", l)) / (eta("2", l) + eta("3", l));
}

std::complex<double> r1(double l)
{
	return exp(eta("0", l) * z1) * ((eta("1", l) + eta("2", l)) * r2(l) * exp(eta("2", l) * z2)
		+ (eta("1", l) - eta("2", l)) * exp(2.0 * eta("2", l) * z1)) / ((eta("1", l) - eta("2", l)) * r2(l) * exp(eta("2", l) * z2) + (eta("1", l) + eta("2", l)) * exp(2.0 * eta("2", l) * z1));
}

std::complex<double> u(double l, double z)
{
	return exp(-eta("0", l) * (abs(z - z0))) / eta("0", l) - exp(-eta("0", l) * (z + z0)) / eta("0", l) + exp(-eta("0", l) * (z + z0)) * (2.0 * (r1(l) * exp(eta("1", l) * z1) + 1.0))
		/ (eta("0", l) * (r1(l) * exp(eta("1", l) * z1) + 1.0) - eta("1", l) * (r1(l) * exp(eta("1", l) * z1) - 1.0));

}

struct dot2D
{
	double z, r;
};

struct dot3D
{
	double z, r, phi;
};


double GetRandomNumberFloat(double min, double max, int precision)
{

	double value;

	// получить случайное число как целое число с порядком precision
	value = rand() % (int)pow(10, precision);

	// получить вещественное число
	value = min + (value / pow(10, precision)) * (max - min);

	return value;
}

int main() {

	// Частица

	double PartSize = 3 * lambda_w;
	std::vector<dot2D> RotGraph;
	RotGraph.resize(30);
	for (int i = 1; i < 29; i++)
	{
		RotGraph[i].z = PartSize / 30 * (i + 1);
		RotGraph[i].r = GetRandomNumberFloat(0, lambda_w / 4, 2);
	}

	RotGraph[0].z = 0;
	RotGraph[0].r = 0;
	RotGraph[29].z = 3 * lambda_w;
	RotGraph[29].r = 0;

	//
	// Сетка для частицы

	dot3D PartGrid[30][16];

	for(int i = 0; i < 16; i++)
		for(int j = 0 ; j < 30; j++)
		{
			PartGrid[i][j].r = RotGraph[j].r;
			PartGrid[i][j].z = RotGraph[j].z;
			PartGrid[i][j].phi = 2 * PI / 16 * (i + 1);
		}


	/*for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 30; j++)
		{
			printf("[%3.2f, %3.2f, %3.2f]", PartGrid[i][j].r, PartGrid[i][j].z, PartGrid[i][j].phi);
		}
		printf("\n");
	}*/
	
	//
	//


	std::cout << k1 << "\n";
	std::cout << k2 << "\n";
	std::cout << k3 << "\n";

	std::cout << eta("2", 1) * z2 << "\n";

	std::cout << u(0.01, 0.01);

}

