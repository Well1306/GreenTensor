
#include <iostream>
#include "Hankel-transform.cpp"
#include <stdlib.h>
#include <complex>

const std::complex<long double> mu0 = 1.25663706 * 1.0E3;
const std::complex<long double> eps0 = 8.8541878e-6;
const std::complex<long double> lambda_w = 500;
const std::complex<long double> eps_up = (long double)1.4 * eps0; //esp в верхнем полупространстве
const std::complex<long double> eps1 = ((long double)-0.33564 * eps0, (long double)-2.9236 * eps0);
const std::complex<long double> eps2 = (long double)1.2 * eps0;
const std::complex<long double> eps3 = (long double)3.0 * eps0;
const std::complex<long double> epsp1 = (long double)2.0 * eps0;
const std::complex<long double> epsp2 = ((long double)-0.33564 * eps0, (long double)-2.9236 * eps0);
std::complex<long double> k0 = (long double)2 * PI / lambda_w; //k в верхнем полупространстве
std::complex<long double> k1 = k0 * sqrt(eps1 / eps_up);
std::complex<long double> k2 = k0 * sqrt(eps2 / eps_up);
std::complex<long double> k3 = k0 * sqrt(eps3 / eps_up);
std::complex<long double> kp1 = k0 * sqrt(epsp1 / eps_up);
std::complex<long double> kp2 = k0 * sqrt(epsp2 / eps_up);

std::complex<long double> z0 = (long double)2 * lambda_w;
std::complex<long double> z1 = (long double)-3.5 * lambda_w;
std::complex<long double> z2 = (long double)-5 * lambda_w;



std::complex<long double> eta(std::string m, long double l)
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


std::complex<long double> r2(long double l)
{
	return exp(eta("2", l) * z2) * (eta("2", l) - eta("3", l)) / (eta("2", l) + eta("3", l));
}

std::complex<long double> r1(long double l)
{
	return exp(eta("0", l) * z1) * ((eta("1", l) + eta("2", l)) * r2(l) * exp(eta("2", l) * z2)
		+ (eta("1", l) - eta("2", l)) * exp((long double)2.0 * eta("2", l) * z1)) / ((eta("1", l) - eta("2", l)) * r2(l) * exp(eta("2", l) * z2) + (eta("1", l) + eta("2", l)) * exp((long double)2.0 * eta("2", l) * z1));
}

std::complex<long double> u(long double l, long double z)
{
	return exp(-eta("0", l) * (abs(z - z0))) / eta("0", l) - exp(-eta("0", l) * (z + z0)) / eta("0", l) + exp(-eta("0", l) * (z + z0)) * ((long double)2.0 * (r1(l) * exp(eta("1", l) * z1) + (long double)1.0))
		/ (eta("0", l) * (r1(l) * exp(eta("1", l) * z1) + (long double)1.0) - eta("1", l) * (r1(l) * exp(eta("1", l) * z1) - (long double)1.0));
}

struct dot2D
{
	long double z, r;
};

struct dot3D
{
	long double z, r, phi;
};


double GetRandomNumberFloat(long double min, long double max, int precision)
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

	long double PartSize = (long double)3 * real(lambda_w);
	std::vector<dot2D> RotGraph;
	RotGraph.resize(30);
	for (int i = 1; i < 29; i++)
	{
		RotGraph[i].z = PartSize / 30 * (i + 1);
		RotGraph[i].r = GetRandomNumberFloat(1, real(lambda_w) / (long double)4, 2);
	}

	RotGraph[0].z = 0;
	RotGraph[0].r = 0;
	RotGraph[29].z = (long double)3 * real(lambda_w);
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


	//std::cout << k1 << "\n";
	//std::cout << k2 << "\n";
	//std::cout << k3 << "\n\n\n\n";

	//std::cout << k0 << "\n";
	//std::cout << eta("2", 11) - eta("3", 11) << "\n";

	std::cout << u(1, 1) << "\n";
	std::cout << u(0, 12312) << "\n";

	std::cout << u(0, 1356666666) << "\n";

	//std::cout << u(6, 1);

}

