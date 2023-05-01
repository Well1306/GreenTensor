
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


std::complex<double> k1;    // верхнее полупространстве подложки
std::complex<double> k2;   // металлический слой
std::complex<double> k3;   // второй слой
std::complex<double> k4;   // нижнее полупространство
std::complex<double> kp1;  // верхнее полупространство частицы 
std::complex<double> kp2;  // нижнее  
std::complex<double> k1n;    // верхнее полупространстве подложки
std::complex<double> k2n;   // металлический слой
std::complex<double> k3n;   // второй слой
std::complex<double> k4n;   // нижнее полупространство
std::complex<double> kp1n;  // верхнее полупространство частицы 
std::complex<double> kp2n;  // нижнее полупространство частицы 

double ro;
double z0;
double _z;
double z1 = -0.15;
double z2 = -1.65;
const int nro = 100;
const double roMin = 0;
const double roMax = 0.1;

std::complex<double> veps[6] =
{ std::complex<double>(eps_up, 0),
  std::complex<double>(-0.33564, 2.9236),
  std::complex<double>(eps_up * 1.2, 0),
  std::complex<double>(eps_up * 1.5, 0),
  std::complex<double>(eps_up * 2.0, 0),            // epsp1
  std::complex<double>(-0.33564, 2.9236)        // epsp2
};



// два полупространства

void U_hs(double tau, double& resRe, double& resIm)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (_z > 0)
		{
			e1 = ctmp * std::exp(-eta1 * (z0 + _z));
		}
		else
		{
			e1 = ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * _z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = ctmp * std::exp(-eta1 * _z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = ctmp * std::exp(eta2 * z0) * std::exp(eta2 * _z);
		}
	}


	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}


double U_hs_real(double tau)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (_z > 0)
		{
			e1 = ctmp * std::exp(-eta1 * (z0 + _z));
		}
		else
		{
			e1 = ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * _z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = ctmp * std::exp(-eta1 * _z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = ctmp * std::exp(eta2 * z0) * std::exp(eta2 * _z);
		}
	}

	return real(e1);
}


double U_hs_imag(double tau)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (_z > 0)
		{
			e1 = ctmp * std::exp(-eta1 * (z0 + _z));
		}
		else
		{
			e1 = ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * _z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = ctmp * std::exp(-eta1 * _z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = ctmp * std::exp(eta2 * z0) * std::exp(eta2 * _z);
		}
	}

	return imag(e1);
}



void W_hs(double tau, double& resRe, double& resIm)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (z > 0)
		{
			e1 = ctmp * epsp2 * std::exp(-eta1 * (z0 + z));
		}
		else
		{
			e1 = ctmp * epsp2 * std::exp(-eta1 * z0) * std::exp(eta2 * z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = ctmp * epsp1 * std::exp(-eta1 * z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = ctmp * epsp1 * std::exp(eta2 * z0) * std::exp(eta2 * z);
		}
	}


	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}

double W_hs_real(double tau)
{

	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (z > 0)
		{
			e1 = ctmp * epsp2 * std::exp(-eta1 * (z0 + z));
		}
		else
		{
			e1 = ctmp * epsp2 * std::exp(-eta1 * z0) * std::exp(eta2 * z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = ctmp * epsp1 * std::exp(-eta1 * z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = ctmp * epsp1 * std::exp(eta2 * z0) * std::exp(eta2 * z);
		}
	}

	return real(e1);
}

double W_hs_imag(double tau)
{

	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (z > 0)
		{
			e1 = ctmp * epsp2 * std::exp(-eta1 * (z0 + z));
		}
		else
		{
			e1 = ctmp * epsp2 * std::exp(-eta1 * z0) * std::exp(eta2 * z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = ctmp * epsp1 * std::exp(-eta1 * z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = ctmp * epsp1 * std::exp(eta2 * z0) * std::exp(eta2 * z);
		}
	}

	return imag(e1);
}



void der_U_hs(double tau, double& resRe, double& resIm)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (_z > 0)
		{
			e1 = -eta1 * ctmp * std::exp(-eta1 * (z0 + _z));
		}
		else
		{
			e1 = eta2 * ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * _z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = -eta1 * ctmp * std::exp(-eta1 * _z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = eta2 * ctmp * std::exp(eta2 * z0) * std::exp(eta2 * _z);
		}
	}


	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}


double der_U_hs_real(double tau)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (_z > 0)
		{
			e1 = -eta1 * ctmp * std::exp(-eta1 * (z0 + _z));
		}
		else
		{
			e1 = eta2 * ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * _z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = -eta1 * ctmp * std::exp(-eta1 * _z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = eta2 * ctmp * std::exp(eta2 * z0) * std::exp(eta2 * _z);
		}
	}

	return real(e1);
}


double der_U_hs_imag(double tau)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 2.0 / (eta1 + eta2);
	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (_z > 0)
		{
			e1 = -eta1 * ctmp * std::exp(-eta1 * (z0 + _z));
		}
		else
		{
			e1 = eta2 * ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * _z);
		}

	}
	else
	{
		if (z > 0)
		{
			e1 = -eta1 * ctmp * std::exp(-eta1 * _z) * std::exp(eta2 * z0);
		}
		else
		{
			e1 = eta2 * ctmp * std::exp(eta2 * z0) * std::exp(eta2 * _z);
		}
	}

	return imag(e1);
}

void g_hs(double tau, double& resRe, double& resIm)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 1.0 / (eta1 * epsp2 + eta2 * epsp1);
	std::complex<double> ctmp1 = (eta1 * epsp2 - eta2 * epsp1);

	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (z > z0)
		{
			e1 = ctmp * ctmp1 * std::exp(-eta1 * (z0 + z)) - std::exp(eta1 * (z0 - z));
		}
		else if (z > 0)
		{
			e1 = ctmp * ctmp1 * std::exp(-eta1 * (z0 + z)) + std::exp(-eta1 * (z0 - z));
		}
		else
		{
			e1 = 2.0 * eta1 * epsp2 * ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * z);
		}
	}
	else
	{
		if (z > 0)
		{
			e1 = -2.0 * eta2 * epsp1 * ctmp * std::exp(-eta1 * z) * std::exp(eta2 * z0);
		}
		else if (z > z0)
		{
			e1 = ctmp * ctmp1 * std::exp(eta2 * (z0 + z)) - std::exp(eta2 * (z0 - z));
		}
		else
		{
			e1 = ctmp * ctmp1 * std::exp(eta2 * (z0 + z)) + std::exp(-eta2 * (z0 - z));
		}
	}
	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}


double g_hs_real(double tau)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 1.0 / (eta1 * epsp2 + eta2 * epsp1);
	std::complex<double> ctmp1 = (eta1 * epsp2 - eta2 * epsp1);

	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (z > z0)
		{
			e1 = ctmp * ctmp1 * std::exp(-eta1 * (z0 + z)) - std::exp(eta1 * (z0 - z));
		}
		else if (z > 0)
		{
			e1 = ctmp * ctmp1 * std::exp(-eta1 * (z0 + z)) + std::exp(-eta1 * (z0 - z));
		}
		else
		{
			e1 = 2.0 * eta1 * epsp2 * ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * z);
		}
	}
	else
	{
		if (z > 0)
		{
			e1 = -2.0 * eta2 * epsp1 * ctmp * std::exp(-eta1 * z) * std::exp(eta2 * z0);
		}
		else if (z > z0)
		{
			e1 = ctmp * ctmp1 * std::exp(eta2 * (z0 + z)) - std::exp(eta2 * (z0 - z));
		}
		else
		{
			e1 = ctmp * ctmp1 * std::exp(eta2 * (z0 + z)) + std::exp(-eta2 * (z0 - z));
		}
	}
	return real(e1);
}

double g_hs_imag(double tau)
{
	std::complex<double> tmp = tau * tau - kp1n * kp1n;
	std::complex<double> eta1 = sqrt(tmp);

	std::complex<double> ctmp = tau * tau - kp2n * kp2n;
	std::complex<double> eta2 = sqrt(ctmp);
	ctmp = 1.0 / (eta1 * epsp2 + eta2 * epsp1);
	std::complex<double> ctmp1 = (eta1 * epsp2 - eta2 * epsp1);

	std::complex<double> e1 = 0;

	if (z0 > 0)
	{
		if (z > z0)
		{
			e1 = ctmp * ctmp1 * std::exp(-eta1 * (z0 + z)) - std::exp(eta1 * (z0 - z));
		}
		else if (z > 0)
		{
			e1 = ctmp * ctmp1 * std::exp(-eta1 * (z0 + z)) + std::exp(-eta1 * (z0 - z));
		}
		else
		{
			e1 = 2.0 * eta1 * epsp2 * ctmp * std::exp(-eta1 * z0) * std::exp(eta2 * z);
		}
	}
	else
	{
		if (z > 0)
		{
			e1 = -2.0 * eta2 * epsp1 * ctmp * std::exp(-eta1 * z) * std::exp(eta2 * z0);
		}
		else if (z > z0)
		{
			e1 = ctmp * ctmp1 * std::exp(eta2 * (z0 + z)) - std::exp(eta2 * (z0 - z));
		}
		else
		{
			e1 = ctmp * ctmp1 * std::exp(eta2 * (z0 + z)) + std::exp(-eta2 * (z0 - z));
		}
	}
	return imag(e1);
}

void for_g2_hs(double tau, double& resRe, double& resIm)
{
	std::complex<double> e1 = 0;

	double resRe_g;
	double resIm_g;
	double resRe_der_u;
	double resIm_der_u;

	g_hs(tau, resRe_g, resIm_g);
	der_U_hs(tau, resRe_der_u, resIm_der_u);


	resRe = resRe_g - resRe_der_u;       // Re
	resIm = resIm_g - resIm_der_u;   // Im
}


double for_g2_hs_real(double tau)
{
	std::complex<double> e1 = 0;

	double resRe_g;
	double resIm_g;
	double resRe_der_u;
	double resIm_der_u;

	g_hs(tau, resRe_g, resIm_g);
	der_U_hs(tau, resRe_der_u, resIm_der_u);


	return resRe_g - resRe_der_u;
}

double for_g2_hs_imag(double tau)
{
	std::complex<double> e1 = 0;

	double resRe_g;
	double resIm_g;
	double resRe_der_u;
	double resIm_der_u;

	g_hs(tau, resRe_g, resIm_g);
	der_U_hs(tau, resRe_der_u, resIm_der_u);


	return resIm_g - resIm_der_u;
}




// два полупространства, разделенные слоями


void U_le(double tau, double& resRe, double& resIm)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//r2
	std::complex<double> r2 = std::exp(eta2 * z2) * (eta2 - eta3) / (eta2 - eta3);

	//r1
	std::complex<double> eta_s = eta1 + eta2;
	std::complex<double> eta_d = eta1 - eta2;
	std::complex<double> r1 = std::exp(eta1 * z1) * (eta_s * r2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * r2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));



	e1 = std::exp(-eta0 * abs(_z - z0)) / eta0 - std::exp(-eta0 * (_z + z0)) / eta0 + std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));

	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}



double U_le_real(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//r2
	std::complex<double> r2 = std::exp(eta2 * z2) * (eta2 - eta3) / (eta2 - eta3);

	//r1
	std::complex<double> eta_s = eta1 + eta2;
	std::complex<double> eta_d = eta1 - eta2;
	std::complex<double> r1 = std::exp(eta1 * z1) * (eta_s * r2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * r2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));



	e1 = std::exp(-eta0 * abs(_z - z0)) / eta0 - std::exp(-eta0 * (_z + z0)) / eta0 + std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));


	return real(e1);    // Im
}

double U_le_imag(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//r2
	std::complex<double> r2 = std::exp(eta2 * z2) * (eta2 - eta3) / (eta2 + eta3);

	//r1
	std::complex<double> eta_s = eta1 + eta2;
	std::complex<double> eta_d = eta1 - eta2;
	std::complex<double> r1 = std::exp(eta1 * z1) * (eta_s * r2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * r2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));



	e1 = std::exp(-eta0 * abs(_z - z0)) / eta0 - std::exp(-eta0 * (_z + z0)) / eta0 + std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));


	return imag(e1);    // Im
}




void W_le(double tau, double& resRe, double& resIm)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//R2
	std::complex<double> R2 = std::exp(eta2 * z2) * (veps[3] * eta2 - veps[2] * eta3) / (veps[3] * eta2 + veps[2] * eta3);

	//R1
	std::complex<double> eta_s = veps[2] * eta1 + veps[1] * eta2;
	std::complex<double> eta_d = veps[2] * eta1 - veps[1] * eta2;
	std::complex<double> R1 = std::exp(eta1 * z1) * (eta_s * R2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * R2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));



	e1 = std::exp(-eta0 * abs(_z - z0)) / eta0 - std::exp(-eta0 * (_z + z0)) / eta0 + std::exp(-eta0 * (_z + z0)) * (2.0 * veps[1] * (R1 * std::exp(eta1 * z1) + 1.0)) /
		(veps[1] * eta0 * (R1 * std::exp(eta1 * z1) + 1.0) - veps[0] * eta1 * (R1 * std::exp(eta1 * z1) - 1.0));

	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}


double W_le_real(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//R2
	std::complex<double> R2 = std::exp(eta2 * z2) * (veps[3] * eta2 - veps[2] * eta3) / (veps[3] * eta2 + veps[2] * eta3);

	//R1
	std::complex<double> eta_s = veps[2] * eta1 + veps[1] * eta2;
	std::complex<double> eta_d = veps[2] * eta1 - veps[1] * eta2;
	std::complex<double> R1 = std::exp(eta1 * z1) * (eta_s * R2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * R2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));



	e1 = std::exp(-eta0 * abs(_z - z0)) / eta0 - std::exp(-eta0 * (_z + z0)) / eta0 + std::exp(-eta0 * (_z + z0)) * (2.0 * veps[1] * (R1 * std::exp(eta1 * z1) + 1.0)) /
		(veps[1] * eta0 * (R1 * std::exp(eta1 * z1) + 1.0) - veps[0] * eta1 * (R1 * std::exp(eta1 * z1) - 1.0));

	return real(e1);        // Re
}


double W_le_imag(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//R2
	std::complex<double> R2 = std::exp(eta2 * z2) * (veps[3] * eta2 - veps[2] * eta3) / (veps[3] * eta2 + veps[2] * eta3);

	//R1
	std::complex<double> eta_s = veps[2] * eta1 + veps[1] * eta2;
	std::complex<double> eta_d = veps[2] * eta1 - veps[1] * eta2;
	std::complex<double> R1 = std::exp(eta1 * z1) * (eta_s * R2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * R2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));



	e1 = std::exp(-eta0 * abs(_z - z0)) / eta0 - std::exp(-eta0 * (_z + z0)) / eta0 + std::exp(-eta0 * (_z + z0)) * (2.0 * veps[1] * (R1 * std::exp(eta1 * z1) + 1.0)) /
		(veps[1] * eta0 * (R1 * std::exp(eta1 * z1) + 1.0) - veps[0] * eta1 * (R1 * std::exp(eta1 * z1) - 1.0));

	return imag(e1);        // Re
}



void g_le(double tau, double& resRe, double& resIm)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//R2
	std::complex<double> R2 = std::exp(eta2 * z2) * (veps[3] * eta2 - veps[2] * eta3) / (veps[3] * eta2 + veps[2] * eta3);

	//R1
	std::complex<double> eta_s = veps[2] * eta1 + veps[1] * eta2;
	std::complex<double> eta_d = veps[2] * eta1 - veps[1] * eta2;
	std::complex<double> R1 = std::exp(eta1 * z1) * (eta_s * R2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * R2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));

	//R0
	eta_s = veps[1] * eta0 + veps[0] * eta1;
	eta_d = veps[1] * eta0 - veps[0] * eta1;
	std::complex<double> R0 = (R1 * std::exp(eta1 * z1) * eta_s + eta_d) / (R1 * std::exp(eta1 * z1) * eta_d + eta_s);


	if (_z > z0)
		e1 = R0 * std::exp(-eta0 * (_z + z0)) - std::exp(-eta0 * (_z - z0));
	else
		e1 = R0 * std::exp(-eta0 * (_z + z0)) + std::exp(-eta0 * (_z - z0));


	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}


double g_le_real(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//R2
	std::complex<double> R2 = std::exp(eta2 * z2) * (veps[3] * eta2 - veps[2] * eta3) / (veps[3] * eta2 + veps[2] * eta3);

	//R1
	std::complex<double> eta_s = veps[2] * eta1 + veps[1] * eta2;
	std::complex<double> eta_d = veps[2] * eta1 - veps[1] * eta2;
	std::complex<double> R1 = std::exp(eta1 * z1) * (eta_s * R2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * R2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));

	//R0
	eta_s = veps[1] * eta0 + veps[0] * eta1;
	eta_d = veps[1] * eta0 - veps[0] * eta1;
	std::complex<double> R0 = (R1 * std::exp(eta1 * z1) * eta_s + eta_d) / (R1 * std::exp(eta1 * z1) * eta_d + eta_s);


	if (_z > z0)
		e1 = R0 * std::exp(-eta0 * (_z + z0)) - std::exp(-eta0 * (_z - z0));
	else
		e1 = R0 * std::exp(-eta0 * (_z + z0)) + std::exp(-eta0 * (_z - z0));


	return real(e1);        // Re
}

double g_le_imag(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//R2
	std::complex<double> R2 = std::exp(eta2 * z2) * (veps[3] * eta2 - veps[2] * eta3) / (veps[3] * eta2 + veps[2] * eta3);

	//R1
	std::complex<double> eta_s = veps[2] * eta1 + veps[1] * eta2;
	std::complex<double> eta_d = veps[2] * eta1 - veps[1] * eta2;
	std::complex<double> R1 = std::exp(eta1 * z1) * (eta_s * R2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * R2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));

	//R0
	eta_s = veps[1] * eta0 + veps[0] * eta1;
	eta_d = veps[1] * eta0 - veps[0] * eta1;
	std::complex<double> R0 = (R1 * std::exp(eta1 * z1) * eta_s + eta_d) / (R1 * std::exp(eta1 * z1) * eta_d + eta_s);


	if (_z > z0)
		e1 = R0 * std::exp(-eta0 * (_z + z0)) - std::exp(-eta0 * (_z - z0));
	else
		e1 = R0 * std::exp(-eta0 * (_z + z0)) + std::exp(-eta0 * (_z - z0));


	return imag(e1);        // Re
}



void der_U_le(double tau, double& resRe, double& resIm)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//r2
	std::complex<double> r2 = std::exp(eta2 * z2) * (eta2 - eta3) / (eta2 - eta3);

	//r1
	std::complex<double> eta_s = eta1 + eta2;
	std::complex<double> eta_d = eta1 - eta2;
	std::complex<double> r1 = std::exp(eta1 * z1) * (eta_s * r2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * r2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));


	if (_z > z0)
		e1 = -eta0 * std::exp(-eta0 * abs(_z - z0)) / eta0 + eta0 * std::exp(-eta0 * (_z + z0)) / eta0 - eta0 * std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));
	else
		e1 = eta0 * std::exp(-eta0 * abs(_z - z0)) / eta0 + eta0 * std::exp(-eta0 * (_z + z0)) / eta0 - eta0 * std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));



	resRe = tau * real(e1);        // Re
	resIm = tau * imag(e1);   // Im
}


double der_U_le_real(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1 * k1);
	std::complex<double> eta1 = sqrt(tau * tau - k2 * k2);
	std::complex<double> eta2 = sqrt(tau * tau - k3 * k3);
	std::complex<double> eta3 = sqrt(tau * tau - k4 * k4);

	std::complex<double> e1 = 0;

	//r2
	std::complex<double> r2 = std::exp(eta2 * z2) * (eta2 - eta3) / (eta2 - eta3);

	//r1
	std::complex<double> eta_s = eta1 + eta2;
	std::complex<double> eta_d = eta1 - eta2;
	std::complex<double> r1 = std::exp(eta1 * z1) * (eta_s * r2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * r2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));


	if (_z > z0)
		e1 = -eta0 * std::exp(-eta0 * abs(_z - z0)) / eta0 + eta0 * std::exp(-eta0 * (_z + z0)) / eta0 - eta0 * std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));
	else
		e1 = eta0 * std::exp(-eta0 * abs(_z - z0)) / eta0 + eta0 * std::exp(-eta0 * (_z + z0)) / eta0 - eta0 * std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));



	return real(e1);        // Re
}

double der_U_le_imag(double tau)
{
	std::complex<double> eta0 = sqrt(tau * tau - k1n * k1n);
	std::complex<double> eta1 = sqrt(tau * tau - k2n * k2n);
	std::complex<double> eta2 = sqrt(tau * tau - k3n * k3n);
	std::complex<double> eta3 = sqrt(tau * tau - k4n * k4n);

	std::complex<double> e1 = 0;

	//r2
	std::complex<double> r2 = std::exp(eta2 * z2) * (eta2 - eta3) / (eta2 - eta3);

	//r1
	std::complex<double> eta_s = eta1 + eta2;
	std::complex<double> eta_d = eta1 - eta2;
	std::complex<double> r1 = std::exp(eta1 * z1) * (eta_s * r2 * std::exp(eta2 * z2) + eta_d * std::exp(2.0 * eta2 * z2)) / (eta_d * r2 * std::exp(eta2 * z2) + eta_s * std::exp(2.0 * eta2 * z2));


	if (_z > z0)
		e1 = -eta0 * std::exp(-eta0 * abs(_z - z0)) / eta0 + eta0 * std::exp(-eta0 * (_z + z0)) / eta0 - eta0 * std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));
	else
		e1 = eta0 * std::exp(-eta0 * abs(_z - z0)) / eta0 + eta0 * std::exp(-eta0 * (_z + z0)) / eta0 - eta0 * std::exp(-eta0 * (_z + z0)) * (2.0 * (r1 * std::exp(eta1 * z1) + 1.0)) /
		(eta0 * (r1 * std::exp(eta1 * z1) + 1.0) - eta1 * (r1 * std::exp(eta1 * z1) - 1.0));



	return imag(e1);        // Re
}


void for_g2_le(double tau, double& resRe, double& resIm)
{
	std::complex<double> e1 = 0;

	double resRe_g;
	double resIm_g;
	double resRe_der_u;
	double resIm_der_u;

	g_le(tau, resRe_g, resIm_g);
	der_U_le(tau, resRe_der_u, resIm_der_u);


	resRe = resRe_g - resRe_der_u;       // Re
	resIm = resIm_g - resIm_der_u;   // Im
}

double for_g2_le_real(double tau)
{
	std::complex<double> e1 = 0;

	double resRe_g;
	double resIm_g;
	double resRe_der_u;
	double resIm_der_u;

	g_le(tau, resRe_g, resIm_g);
	der_U_le(tau, resRe_der_u, resIm_der_u);


	return resRe_g - resRe_der_u;       // Re
}


double for_g2_le_imag(double tau)
{
	std::complex<double> e1 = 0;

	double resRe_g;
	double resIm_g;
	double resRe_der_u;
	double resIm_der_u;

	g_le(tau, resRe_g, resIm_g);
	der_U_le(tau, resRe_der_u, resIm_der_u);


	return resIm_g - resIm_der_u;       // Re
}


void calc(std::string mode,double* ro, double re, double im, double* resRe, double* resIm, double* resMod, double hro, double calc_z0, double calc_z, std::string name)
{
	z0 = calc_z0;
	_z = calc_z;
	for (int j = 0; j < nro; j++)
	{
		if (mode == "U_hs")
			U_hs(ro[j], re, im);
		else if(mode == "W_hs")
			W_hs(ro[j], re, im);
		else if(mode == "for_g2_hs")
			for_g2_hs(ro[j], re, im);
		else if(mode == "U_le")
			U_le(ro[j], re, im);
		else if (mode == "W_le")
			W_le(ro[j], re, im);
		else if (mode == "for_g2_le")
			for_g2_le(ro[j], re, im);
		resRe[j] = re;
		resIm[j] = im;
		resMod[j] = sqrt(re * re + im * im);
	}
	std::vector<double> filter = MakeFilter(hro, roMax, 0);
	std::ofstream out2;          // поток для записи
	out2.open("real.txt");
	out2 << nro << '\n';
	for (int i = 0; i < nro; i++)
	{
		out2 << ro[i] << '\n';
	}
	std::vector<double> G1_real;
	if (mode == "U_hs")
		G1_real = HankelTransform(filter, hro, roMax, U_hs_real);
	else if (mode == "W_hs")
		G1_real = HankelTransform(filter, hro, roMax, W_hs_real);
	else if (mode == "for_g2_hs")
		G1_real = HankelTransform(filter, hro, roMax, for_g2_hs_real);
	else if (mode == "U_le")
		G1_real = HankelTransform(filter, hro, roMax, U_le_real);
	else if (mode == "W_le")
		G1_real = HankelTransform(filter, hro, roMax, W_le_real);
	else if (mode == "for_g2_le")
		G1_real = HankelTransform(filter, hro, roMax, for_g2_le_real);
	out2 << name << '\n';
	for (int i = 0; i < nro; i++)
	{
		out2 << G1_real[i] << '\n';
	}
	std::ofstream out3;          // поток для записи
	out3.open("imag.txt");
	out3 << nro << '\n';
	for (int i = 0; i < nro; i++)
	{
		out3 << ro[i] << '\n';
	}
	std::vector<double> G1_imag = HankelTransform(filter, hro, roMax, for_g2_le_imag);
	if (mode == "U_hs")
		G1_imag = HankelTransform(filter, hro, roMax, U_hs_imag);
	else if (mode == "W_hs")
		G1_imag = HankelTransform(filter, hro, roMax, W_hs_imag);
	else if (mode == "for_g2_hs")
		G1_imag = HankelTransform(filter, hro, roMax, for_g2_hs_imag);
	else if (mode == "U_le")
		G1_imag = HankelTransform(filter, hro, roMax, U_le_imag);
	else if (mode == "W_le")
		G1_imag = HankelTransform(filter, hro, roMax, W_le_imag);
	else if (mode == "for_g2_le")
		G1_imag = HankelTransform(filter, hro, roMax, for_g2_le_imag);
	out3 << name << '\n';
	for (int i = 0; i < nro; i++)
	{
		out3 << G1_imag[i] << '\n';
	}
	std::ofstream out4;          // поток для записи
	out4.open("mod.txt");
	out4 << nro << '\n';
	for (int i = 0; i < nro; i++)
	{
		out4 << ro[i] << '\n';
	}
	out4 << name << '\n';
	for (int i = 0; i < nro; i++)
	{
		out4 << sqrt(G1_imag[i] * G1_imag[i] + G1_real[i] * G1_real[i]) << '\n';
	}
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




	calc("for_g2_le", ro, re, im, resRe, resIm, resMod, hro, 0.1, 1, "h1 = 0.15, h2 = 1.5, z > z0");


	z1 = -1.5;
	z2 = -1.65;
	calc("for_g2_le", ro, re, im, resRe, resIm, resMod, hro, 1, 0.1, "h1 = 1.5, h2 = 0.15, z < z0");
	



}

