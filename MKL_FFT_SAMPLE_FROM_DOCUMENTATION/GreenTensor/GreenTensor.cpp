
#include <iostream>
#include "Hankel-transform.cpp"
#include <stdlib.h>
#include <complex>
#include <cmath>


const double cLight = 299792458;
const double sqrt_eps0mu0 = (long double)1.0 / cLight;
const double mu0 = 4 * PI * 1.0E-7;
const double eps0 = 8.8541878e-12;
const double lambda_w = 500 * 1.0E-9;
//const std::complex<long double> eps_up = (long double)1.4 * eps0; //esp в верхнем полупространстве
//const std::complex<long double> eps1 = ((long double)-0.33564 * eps0, (long double)-2.9236 * eps0);
//const std::complex<long double> eps2 = (long double)1.2 * eps0;
//const std::complex<long double> eps3 = (long double)3.0 * eps0;
//const std::complex<long double> epsp1 = (long double)2.0 * eps0;
//const std::complex<long double> epsp2 = ((long double)-0.33564 * eps0, (long double)-2.9236 * eps0);
const double eps_up = 1.4; //esp в верхнем полупространстве
const std::complex<double> eps1 = (-0.33564, -2.9236);
const double eps2 = (long double)1.2;
const double eps3 = (long double)1.5;
const std::complex<double> epsp1 = 2.0;
const std::complex<double> epsp2 = (-0.33564, -2.9236);
std::complex<long double> k0 = (long double)2 * PI / lambda_w; //k в верхнем полупространстве


double k1;    // верхнее полупространстве подложки
std::complex<double> k2;   // металлический слой
std::complex<double> k3;   // второй слой
std::complex<double> k4;   // нижнее полупространство
std::complex<double> kp1;  // верхнее полупространство частицы 
std::complex<double> kp2;  // нижнее  
double k1n;    // верхнее полупространстве подложки
std::complex<double> k2n;   // металлический слой
std::complex<double> k3n;   // второй слой
std::complex<double> k4n;   // нижнее полупространство
std::complex<double> kp1n;  // верхнее полупространство частицы 
std::complex<double> kp2n;  // нижнее полупространство частицы 

double ro;
double z0;
double z1;
const int nro = 50;
const double roMin = 0;
const double roMax = 10;

std::complex<double> veps[6] =
{ std::complex<double>(eps_up, 0),
  std::complex<double>(-0.33564, 2.9236),
  std::complex<double>(eps_up * 1.2, 0),
  std::complex<double>(eps_up * 1.5, 0),
  std::complex<double>(eps_up * 2.0, 0),            // epsp1
  std::complex<double>(-0.33564, 2.9236)        // epsp2
};

void U_hs(double tau,  double& resRe, double& resIm)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (z1 > 0)
		{
			e1 = ctmp * std::exp(-eta1 * (z0 + z1));
		}
		else
		{
			e1 = ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * z1);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = ctmp * std::exp(-eta1 * z1) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = ctmp * std::exp(eta2 * z0) * std::exp(eta2 * z1);
		}
	}


	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
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

	long double value;

	// получить случайное число как целое число с порядком precision
	value = rand() % (int)pow(10, precision);

	// получить вещественное число
	value = min + (value / pow(10, precision)) * (max - min);

	return value;
}

int main() {

	std::complex<double>* kw = new std::complex<double>[6];

	double omega = 2.0 * PI * cLight / lambda_w / sqrt(eps_up);

	for (int j = 0; j < 6; j++)
		kw[j] = omega * sqrt_eps0mu0 * sqrt(veps[j]);

	k1 = real(kw[0]);
	k2 = kw[1];
	k3 = kw[2];
	k4 = kw[3];
	kp1 = kw[4];
	kp2 = kw[5];

	k1n = 1;
	k2n = k2 / k1;
	k3n = k3 / k1;
	k4n = k4 / k1;
	kp1n = kp1 / k1;
	kp2n = kp2 / k1;


	z0 = 0.01;
	z1 = 0.1;

	double* ro = new double[nro];
	double hro = (roMax - roMin) / nro;

	for (int j = 0; j < nro; j++)
	{
		ro[j] = roMin + j * hro;
	}


	double* resRe = new double[nro];
	double* resIm = new double[nro];
	double* resMod = new double[nro];


	double re = 0;
	double im = 0;

	for (int j = 0; j < nro; j++)
	{
		U_hs(ro[j], re, im);
		resRe[j] = re;
		resIm[j] = im;
		resMod[j] = sqrt(re * re + im * im);
	}

	std::vector<double> filter = MakeFilter(hro, roMax, 0);


	for (int j = 0; j < nro; j++)
	{
		std::cout << j << ": " << filter[j] << "\n";
	}

	std::vector<double> hankel = HankelTransform(filter, hro, roMax);





	std::ofstream out1;          // поток для записи
	out1.open("test.txt");

	out1 << "en\n###\nЗависимость электромагнитного поля от коэффицента z\nРасстояние\nЭлектромагнитное поле\n100\n";
	for (int i = 0; i < nro; i++)
	{
		out1 << ro[i] << '\n';
	}
	out1 << "real" << '\n';
	for (int i = 0; i < nro; i++)
	{
		out1 << resRe[i] << '\n';
	}
	out1 << "imag" << '\n';
	for (int i = 0; i < nro; i++)
	{
		out1 << resIm[i] << '\n';
	}
	out1 << "mod" << '\n';
	for (int i = 0; i < nro; i++)
	{
		out1 << resMod[i] << '\n';
	}

	// Частица

	long double PartSize = (long double)3 * lambda_w;
	std::vector<dot2D> RotGraph;
	RotGraph.resize(30);
	for (int i = 1; i < 29; i++)
	{
		RotGraph[i].z = PartSize / 30 * (i + 1);
		RotGraph[i].r = GetRandomNumberFloat(1, lambda_w / (long double)4, 2);
	}

	RotGraph[0].z = 0;
	RotGraph[0].r = 0;
	RotGraph[29].z = (long double)3 * lambda_w;
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


	//std::cout << z << "\n";
	//std::cout << k2 << "\n";
    //std::cout << k3 << "\n\n\n\n";

	//std::cout << k0 << "\n";
	//std::cout << eta("2", 11) - eta("3", 11) << "\n";

	std::ofstream out;          // поток для записи
	out.open("hello.txt");

	//std::cout << u(123213, 213123) << "\n";
	//std::cout << u(0, 0) << "\n\n\n\n!!!\n";
	std::complex<long double> s = 0;
	out << "\n\n\n";
	for (int i = 0; i < 10; i++)
	{
		out << i << '\n';
	}
	out << "real\n";
	for (int i = 0; i < 10; i++)
	{
		//out << real(u(i, 0.5)) << '\n';
	}
	out << "imag\n";
	for (int i = 0; i < 10; i++)
	{
		//out << imag(u(i, 0.5)) << '\n';
	}
	

	//std::cout << u(6, 1);

}

