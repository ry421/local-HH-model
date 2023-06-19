//============================================================================
// Name        : LocalHH.cpp
// Author      : Ryan Nesselrodt
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

//This is a code for calculating the local atomic green's function for the holstien-Hubbard model with time dependent couplings g(t),U(t)

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;
typedef complex<double> cx;

// Definitions
const cx one(1, 0);
const cx I(0, 1);
double pi = 4 * atan(1);

//define physical parameters
double mu = 0.25;
double beta = 1;
double m = 2;

//freqency
double w = 0.5;
double fw = ((double)1) / ((double)(1 - exp(-beta * w)));

//COMPUTATIONAL PARAMS
//define number of rectangles in integrals

double dt = 0.1;

double tmax = 250;
int ntot = round((double)tmax / (double)dt);

vector<cx> CC(2 * ntot + 1);
vector<cx> intf(2 * ntot + 1);
vector<cx> cosi(2 * ntot + 1);
vector<cx> Uint(2 * ntot + 1);

//Define all the necessary functions to calculate the local nonequilibrium Holstien GF

//e-ph coupling
double g(double t)
{
	//return 1-0.005*exp(-0.1*t);
	//return 0.4-0.1*sin(0.001*t)*exp(-0*t);
	return 1-0.1*sin(0.1*t)*exp(-0.01*t);
	//return 1-0.1*sin(0.1*t);
	//return 1+t*t;
}

double dgdt(double t, double dt)
{
	return ((double)1. / (double)(2 * dt)) * (g(t + dt) - g(t - dt));
}

double d2gdt(double t, double dt)
{
	return ((double)1 / (double)(dt * dt)) * (g(t + dt) - 2 * g(t) + g(t - dt));
}

//more definitions
double g0 = g(0);

//electron-electron coupling
double U(double t)
{
	//return 1-0.01*t*t;
	return 1+0.01*sin(0.1*t);
	//return 1;
}
double U0 = U(0);

double d2Udt(double t, double dt)
{
	return ((double)1 / (double)(dt * dt)) * (U(t + dt) - 2 * U(t) + U(t - dt));
}

//build partition function

double factor1 = exp(beta * mu + ((double)beta * g0 * g0) / ((double)2 * m * w * w));
double factor2 = exp(beta * (2. * mu - U0) + ((double)2 * beta * g0 * g0) / ((double)m * w * w));
double Z = fw * (1. + 2. * factor1 + factor2);
double Zm = (1. + 2. * factor1 + factor2);
double overZ = ((double)1 / (double)Z);
double overZm = ((double)1 / (double)Zm);

//build ave nf

double avenf = overZm * (2. * factor1 + 2. * factor2);

double avenupndown = overZm * exp(beta * (2 * mu - U0) + ((double)2 * beta * g0 * g0) / ((double)m * w * w));

double avencube = avenf + 6. * avenupndown;

//FUNCTIONS NEEDED TO GET MOMENTS

//function for expectation value of x as function of time
double xave(double t)
{

	int num = round((double)(t) / (double)(dt));
	cx sum(0, 0);

	for (int i = 0; i <= num; i++)
	{
		double time = i * dt;
		sum += exp(I * w * time) * g(time) * dt;
	}

	double cons = ((double)1 / (double)(m * w));

	cx number = I * exp(-I * w * t) * sum * cons;
	double realio = number.real();

	double facy = ((double)g0 / ((double)(m * w * w)));

	return -avenf * (facy * cos(w * t) + realio);
}

double dxavedt(double t, double dt)
{
	return ((double)(1.) / (double)(2 * dt)) * (xave(t + dt) - xave(t - dt));
}

//average <n_bsigma x(t)>
double avenx(double t)
{
	int num = round((double)(t) / (double)(dt));
	cx sum(0, 0);

	for (int i = 0; i <= num; i++)
	{
		double time = i * dt;
		sum += exp(I * w * time) * g(time) * dt;
	}

	double cons = ((double)1 / (double)(m * w));

	cx number = I * exp(-I * w * t) * sum * cons;
	double realio = number.real();

	double facy = ((double)g0 / ((double)(m * w * w)));

	return -(0.5 * avenf + avenupndown) * (facy * cos(w * t) + realio);
}

//average <x^2(t)>
double avex2(double t)
{
	int num = round((double)(t) / (double)(dt));
	cx sum(0, 0);

	for (int i = 0; i <= num; i++)
	{
		double time = i * dt;
		sum += exp(I * w * time) * g(time) * dt;
	}

	double cons = ((double)1 / (double)(m * w));

	cx number = I * exp(-I * w * t) * sum * cons;
	double realio = number.real();

	double facy = ((double)g0 / ((double)(m * w * w)));

	return (avenf + 2 * avenupndown) * (facy * cos(w * t) + realio) * (facy * cos(w * t) + realio) + ((double)1 / (double)(2 * m * w * tanh(0.5 * beta * w)));
}

//average <n_bsig x^2(t)>
double avenx2(double t)
{
	int num = round((double)(t) / (double)(dt));
	cx sum(0, 0);

	for (int i = 0; i <= num; i++)
	{
		double time = i * dt;
		sum += exp(I * w * time) * g(time) * dt;
	}

	double cons = ((double)1 / (double)(m * w));

	cx number = I * exp(-I * w * t) * sum * cons;
	double realio = number.real();

	double facy = ((double)g0 / ((double)(m * w * w)));

	return (0.5 * avenf + 3 * avenupndown) * (facy * cos(w * t) + realio) * (facy * cos(w * t) + realio) + ((double)1 / (double)(2 * m * w * tanh(0.5 * beta * w))) * 0.5 * avenf;
}

//variance of x operator, <x^2>-<x>^2
double varx(double t)
{
	double consty = avenf - avenf * avenf + 2 * avenupndown;

	int num = round((double)(t) / (double)(dt));
	cx sum(0, 0);

	for (int i = 0; i <= num; i++)
	{
		double time = i * dt;
		sum += exp(I * w * time) * g(time) * dt;
	}

	double cons = ((double)1 / (double)(m * w));

	cx number = I * exp(-I * w * t) * sum * cons;
	double realio = number.real();

	double facy = ((double)g0 / ((double)(m * w * w)));

	double inter = (facy * cos(w * t) + realio);

	return consty * inter * inter + ((double)1 / (double)(2 * m * w * tanh(0.5 * beta * w)));
}

//return variance of x^3 operator, <x^3>-<x>^3

double varxcube(double t)
{
	//build the real integral
	int num = round((double)(t) / (double)(dt));
	cx sum(0, 0);

	for (int i = 0; i <= num; i++)
	{
		double time = i * dt;
		sum += exp(I * w * time) * g(time) * dt;
	}

	double cons = ((double)1 / (double)(m * w));

	cx number = I * exp(-I * w * t) * sum * cons;
	double realio = number.real();

	/*	//build other constants
	double c1= ((double)(g0*g0*g0))/((double)(4*m*m*m*w*w*w*w*w*w));
	double c2 =((double)(3*g0*g0))/((double)(2*m*m*w*w*w*w));
	double c3 = ((double)(3*g0))/((double)(2*m*w*w));

	//term proportional to <n^3>
	double term1 = -avencube*(c1*cos(3*w*t)+3.*c1*cos(w*t)+c2*realio*(cos(2*w*t)+1)+c3*realio*realio*cos(w*t)+realio*realio*realio);*/

	//other constants
	double a1 = ((double)(g0) / (double)(m * m * w * w * w));
	double a2 = ((double)1 / (double)(m * w));

	double a3 = (double)g0 / (double)(m * w * w);

	double cothb = ((double)1 / (double)tanh(0.5 * beta * w));

	//term proportional to <n>
	double term2 = -avenf * cothb * 1.5 * a2 * (a3 * cos(w * t) + realio);
	//double term2=1.5*a2*cothb*xave(t);

	//term proportional to <n^3>-<n>^3
	double term3 = -(avenf * (1. - avenf * avenf) + 6 * avenupndown) * (a3 * cos(w * t) + realio) * (a3 * cos(w * t) + realio) * (a3 * cos(w * t) + realio);
	return term2 + term3;
}

//FUNCTIONS TO COMPUTE MOMENTS OF GF

//return first moment of the GF
double firstmom(double t)
{
	return -mu + g(t) * xave(t) + U(t) * 0.5 * avenf;
}
//return second moment of the GF

double secondmom(double t)
{
	return (-mu + g(t) * xave(t) + 0.5 * U(t) * avenf) * (-mu + g(t) * xave(t) + 0.5 * U(t) * avenf) + g(t) * g(t) * (varx(t)) + U(t) * U(t) * (0.5 * avenf - 0.25 * avenf * avenf) + 2. * U(t) * g(t) * (avenx(t) - xave(t) * 0.5 * avenf);
}

//extra terms for third moment
double extraterms(double t)
{
	return U(t) * U(t) * U(t) * (0.5 * avenf - 0.125 * avenf * avenf * avenf) - 3 * mu * U(t) * U(t) * (0.5 * avenf - 0.25 * avenf * avenf) + 3 * U(t) * U(t) * g(t) * (avenx(t) - 0.25 * avenf * avenf * xave(t)) - 6 * mu * U(t) * g(t) * (avenx(t) - 0.5 * avenf * xave(t)) + 3 * U(t) * g(t) * g(t) * (avenx2(t) - 0.5 * avenf * xave(t) * xave(t));
}

//third moment of the GF
double thirdmom(double t)
{
	//	return firstmom(t)*firstmom(t)*firstmom(t)+0.25*g(t)*w*w*xave(t)+g(t)*g(t)*g(t)*varxcube(t)-3.*mu*g(t)*g(t)*varx(t)+0.5*g(t)*g(t)*((double) 1/(double) m)*(1.-0.5*avenf);
	//return firstmom(t)*firstmom(t)*firstmom(t)+g(t)*w*w*xave(t)+g(t)*g(t)*g(t)*varxcube(t)-3.*mu*g(t)*g(t)*varx(t);
	return firstmom(t) * firstmom(t) * firstmom(t) + 1 * g(t) * w * w * xave(t) + g(t) * g(t) * g(t) * varxcube(t) - 3. * mu * g(t) * g(t) * varx(t) + g(t) * g(t) * 0.5 * ((double)1 / (double)m) * (0.5 * avenf + 1.) - 0.25 * (d2gdt(t, dt) * xave(t) + d2Udt(t, dt)) + extraterms(t) - 0.5 * dxavedt(t, dt) * dgdt(t, dt);
}

double thirdmomwrong(double t)
{
	//moment as analytic third deriv
	return firstmom(t) * firstmom(t) * firstmom(t) + 0.25 * g(t) * w * w * xave(t) + g(t) * g(t) * g(t) * varxcube(t) - 3. * mu * g(t) * g(t) * varx(t) + g(t) * g(t) * 0.5 * ((double)1 / (double)m) * (-0.5 * avenf + 1.) - 0.25 * (d2gdt(t, dt) * xave(t) + d2Udt(t, dt)) + extraterms(t) - 0.5 * dgdt(t, dt) * dxavedt(t, dt);
}

//FUCTIONS TO COMPUTE THE GF

//define f(t) the prefactor with an integral over sine
double f(double t)
{

	int num = round((double)t / (double)dt);

	double sum = 0;

	if (num >= 0)
	{
		for (int i = 0; i <= num; i++)
		{
			double tp = i * dt;
			sum += dt * g(tp) * sin(w * (t - tp));
		}
	}

	else
	{
		for (int i = 0; i >= num; i -= 1)
		{
			double tp = i * dt;
			sum += dt * g(tp) * sin(w * (t - tp));
		}
	}
	return g(t) * sum;
}

//Define the function C that appears frequently in our expressions
//has one internal integration
std::complex<double> C(double t)
{
	//make coefficient
	std::complex<double> coeff((double)1 / (double)(std::sqrt(2 * m * w)), 0);

	int num = round((double)t / (double)dt);
	std::complex<double> sum(0, 0);
	if (num >= 0)
	{
		for (int i = 0; i <= num; i++)
		{
			double tp = i * dt;
			sum += g(tp) * exp(std::complex<double>(0, w * tp)) * dt;
		}
	}
	else
	{
		for (int i = 0; i >= num; i -= 1)
		{
			double tp = i * dt;
			sum += g(tp) * exp(std::complex<double>(0, w * tp)) * dt;
		}
	}

	return coeff * sum;
}

double cosint(double t)
{

	int num = round((double)t / (double)dt);
	double sum = 0;
	if (num >= 0)
	{
		for (int i = 0; i <= num; i++)
		{
			double tp = i * dt;
			sum += g(tp) * cos(w * tp) * dt;
		}
	}
	else
	{
		for (int i = 0; i >= num; i -= 1)
		{
			double tp = i * dt;
			sum += g(tp) * cos(w * tp) * dt;
		}
	}
	double coeff = ((double)g0) / ((double)m * w * w);

	return coeff * sum;
}

//integral over U
std::complex<double> Uin(double t)
{

	int num = round((double)t / (double)dt);

	std::complex<double> sum(0, 0);
	if (num >= 0)
	{
		for (int i = 0; i <= num; i++)
		{
			double tp = i * dt;
			sum += U(tp) * dt;
		}
	}
	else
	{
		for (int i = 0; i >= num; i -= 1)
		{
			double tp = i * dt;
			sum += U(tp) * dt;
		}
	}

	return -I * sum;
}

//define greater and lesser green's functions by accessing elements from vectors we've built

cx Ggreater(double t1, double t2)
{

	int nt1 = round((double)t1 / (double)dt) + ntot;
	int nt2 = round((double)t2 / (double)dt) + ntot;

	//define quantites that go in exponents
	double c = 0.5 - fw;
	cx mag = CC[nt1] * conj(CC[nt1]) + CC[nt2] * conj(CC[nt2]);

	cx imy = -1. * (CC[nt1] * conj(CC[nt2]));
	double fac = 2. * fw * (CC[nt1] * conj(CC[nt2])).real();

	cx magi = abs(CC[nt1] - CC[nt2]);
	std::complex<double> coeff(0, (double)1 / (double)(2 * m * w));

	cx Qew = coeff * (intf[nt1] - intf[nt2]);

	//define prefactor
	cx prefactor = -I * overZm * exp(I * mu * (t1 - t2) + one * 0.5 * (mag) + imy - fw * magi * magi);

	//second term in sum
	cx factor2 = exp(one * beta * mu + one * ((double)(g0 * g0 * beta) / ((double)2 * m * w * w)) + one * 3. * Qew + cosi[nt1] - cosi[nt2] + Uint[nt1] - Uint[nt2]);

	return prefactor * (exp(Qew) + factor2);
}

/*
cx Ggreater(double t1,double t2){

	int nt1= round((double)t1/(double)dt)+ntot;
	int nt2=round((double)t2/(double)dt)+ntot;



	//define quantites that go in exponents
	double c=0.5-fw;
	cx  mag=CC[nt1]*conj(CC[nt1])+CC[nt2]*conj(CC[nt2]);

	cx imy= -I*(CC[nt1]*conj(CC[nt2])).imag();
	double fac=2.*fw*(CC[nt1]*conj(CC[nt2])).real();


	cx magi =abs(CC[nt1]-CC[nt2]);
	std::complex<double> coeff(0,(double)1/(double)(2*m*w));

	cx Qew= coeff*(intf[nt1]-intf[nt2]);


	//define prefactor
	cx prefactor=-I*overZm*exp(I*mu*(t1-t2)+one*c*(magi*magi)+imy);

	//second term in sum
	cx factor2=exp(one*beta*mu+one*((double)(g0*g0*beta)/((double)2*m*w*w))+one*3.*Qew +cosi[nt1]-cosi[nt2]);




	return prefactor*(exp(Qew)+factor2);
}

*/

cx Glesser(double t1, double t2)
{

	int nt1 = round((double)t1 / (double)dt) + ntot;
	int nt2 = round((double)t2 / (double)dt) + ntot;

	std::complex<double> coeff(0, (double)1 / (double)(2 * m * w));

	cx Qew = coeff * (intf[nt1] - intf[nt2]);

	cx magi = abs(CC[nt1] - CC[nt2]);

	//define quantites that go in exponents
	double c = 0.5 - fw;
	cx mag = CC[nt1] * conj(CC[nt1]) + CC[nt2] * conj(CC[nt2]);

	cx imy = -1. * (CC[nt2] * conj(CC[nt1]));
	double fac = 2. * fw * (CC[nt1] * conj(CC[nt2])).real();

	//define prefactor
	cx prefactor = I * overZm * exp(I * mu * (t1 - t2) + one * 0.5 * (mag) + imy - fw * magi * magi);

	//first term in sum
	cx factor1 = exp(one * beta * mu + one * ((double)(g0 * g0 * beta) / ((double)2 * m * w * w)) + one * Qew + cosi[nt1] - cosi[nt2]);

	//second term in sum

	cx factor2 = exp(one * 2. * beta * mu - one * beta * U0 + 2. * one * ((double)(g0 * g0 * beta) / ((double)m * w * w)) + 3. * one * Qew + 2. * (cosi[nt1] - cosi[nt2]) + Uint[nt1] - Uint[nt2]);

	return prefactor * (factor1 + factor2);
}

/*
cx Glesser(double t1,double t2){

	int nt1= round((double)t1/(double)dt)+ntot;
	int nt2=round((double)t2/(double)dt)+ntot;

	std::complex<double> coeff(0,(double)1/(double)(2*m*w));

	cx Qew= coeff*(intf[nt1]-intf[nt2]);

	cx magi =abs(CC[nt1]-CC[nt2]);


	//define quantites that go in exponents
	double c=0.5-fw;
	cx mag=CC[nt1]*conj(CC[nt1])+CC[nt2]*conj(CC[nt2]);

	cx imy= I*(CC[nt1]*conj(CC[nt2])).imag();
	double fac=2.*fw*(CC[nt1]*conj(CC[nt2])).real();


	//define prefactor
	cx prefactor=I*overZm*exp(I*mu*(t1-t2)+one*c*(magi*magi)+imy);

	//first term in sum
	cx factor1=exp(one*beta*mu+one*((double)(g0*g0*beta)/((double)2*m*w*w))+one*Qew+cosi[nt1]-cosi[nt2]);

	//second term in sum

	cx factor2=exp(one*2.*beta*mu+2.*one*((double)(g0*g0*beta)/((double)m*w*w))+3.*one*Qew+2.*(cosi[nt1]-cosi[nt2]));








	return prefactor*(factor1+factor2);
}

*/

double theta1(double t1, double t2)
{
	if (t1 >= t2)
	{
		return 1;
	}
	else
		return 0;
}

double theta0(double t1, double t2)
{
	if (t1 > t2)
	{
		return 1;
	}
	else
		return 0;
}
cx gR(double t1, double t2)
{
	//define quantities that go in exponents

	return (theta1(t1, t2) * Ggreater(t1, t2) - theta1(t1, t2) * Glesser(t1, t2));
}

cx g(double t1, double t2)
{
	return Ggreater(t1, t2) - Glesser(t1, t2);
}

cx gA(double t1, double t2)
{
	return theta0(t2, t1) * (Glesser(t1, t2) - Ggreater(t1, t2));
}

//FUNCTIONS TO COMPUTE NUMERICAL DERIVATIVES FOR MOMENTS

//first deriv
cx deriv1(double t1, double t2, double tstep)
{

	double realy = ((double)(g(t1 + tstep, t2).real() - g(t1, t2).real()) / ((double)tstep));
	double imagi = ((double)(g(t1 + tstep, t2).imag() - g(t1, t2).imag()) / ((double)tstep));
	return one * realy + I * imagi;
}

cx deriv2(double t1, double t2, double tstep)
{
	double realy = ((double)(g(t1, t2).real() - g(t1, t2 - tstep).real()) / ((double)tstep));
	double imagi = ((double)(g(t1, t2).imag() - g(t1, t2 - tstep).imag()) / ((double)tstep));
	return one * realy + I * imagi;
}
cx totalderiv(double t1, double t2, double tstep)
{
	return 0.5 * (deriv1(t1, t2, tstep) - deriv2(t1, t2, tstep));
}

//functions for calculating second derivative

//pieces needed
cx d2t1(double t1, double t2, double tstep)
{
	return ((cx)(g(t1 + tstep, t2) - 2. * g(t1, t2) + g(t1 - tstep, t2)) / ((cx)tstep * tstep));
}

cx d2t2(double t1, double t2, double tstep)
{
	return ((cx)(g(t1, t2 + tstep) - 2. * g(t1, t2) + g(t1, t2 - tstep)) / ((cx)tstep * tstep));
}

cx d2t1t2(double t1, double t2, double tstep)
{
	return 0.25 * ((cx)(g(t1 + tstep, t2 + tstep) - g(t1 + tstep, t2 - tstep) + g(t1 - tstep, t2 - tstep) - g(t1 - tstep, t2 + tstep)) / ((cx)tstep * tstep));
}

//total second deriv wrt t2
cx d2dtr2(double t1, double t2, double tstep)
{
	return 0.25 * (d2t1(t1, t2, tstep) + d2t2(t1, t2, tstep) - 2. * d2t1t2(t1, t2, tstep));
}

//functions for calculating third deriv numerically

//components required

cx d3t1(double t1, double t2, double tstep)
{
	return ((cx)one / (cx)(tstep * tstep * tstep)) * (-0.5 * g(t1 - 2 * tstep, t2) + g(t1 - tstep, t2) - g(t1 + tstep, t2) + 0.5 * g(t1 + 2 * tstep, t2));
}

cx d3t2(double t1, double t2, double tstep)
{
	return ((cx)one / (cx)(tstep * tstep * tstep)) * (-0.5 * g(t1, t2 - 2 * tstep) + g(t1, t2 - tstep) - g(t1, t2 + tstep) + 0.5 * g(t1, t2 + 2 * tstep));
}

cx d3dt1(double t1, double t2, double tstep)
{
	return 0.125 * ((cx)one / (cx)(tstep * tstep * tstep)) * (g(t1 + 3 * tstep, t2) - 3. * g(t1 + tstep, t2) + 3. * g(t1 - tstep, t2) - g(t1 - 3 * tstep, t2));
}

cx d3dt2(double t1, double t2, double tstep)
{
	return 0.125 * ((cx)one / (cx)(tstep * tstep * tstep)) * (g(t1, t2 + 3 * tstep) - 3. * g(t1, t2 + tstep) + 3. * g(t1, t2 - tstep) - g(t1, t2 - 3 * tstep));
}

cx d3t12t2(double t1, double t2, double tstep)
{
	return 0.125 * ((cx)one / (cx)(tstep * tstep * tstep)) * (g(t1 + 2 * tstep, t2 + tstep) - g(t1 + 2 * tstep, t2 - tstep) - 2. * g(t1, t2 + tstep) + 2. * g(t1, t2 - tstep) + g(t1 - 2 * tstep, t2 + tstep) - g(t1 - 2 * tstep, t2 - tstep));
}

cx d3t1t22(double t1, double t2, double tstep)
{
	return 0.125 * ((cx)one / (cx)(tstep * tstep * tstep)) * (g(t1 + tstep, t2 + 2 * tstep) - g(t1 - tstep, t2 + 2 * tstep) - 2. * g(t1 + tstep, t2) + 2. * g(t1 - tstep, t2) + g(t1 + tstep, t2 - 2 * tstep) - g(t1 - tstep, t2 - 2 * tstep));
}

//d^3g/dt_1^2dt2

cx d3dt11dt2(double t1, double t2, double tstep)
{
	return 0.5 * ((cx)one / (cx)(tstep * tstep * tstep)) * (g(t1 + tstep, t2 + tstep) - 2. * g(t1, t2 + tstep) + g(t1 - tstep, t2 + tstep) - g(t1 + tstep, t2 - tstep) + 2. * g(t1, t2 - tstep) - g(t1 - tstep, t2 - tstep));
}

//d^3g/dt1dt2^2

cx d3dt1dt22(double t1, double t2, double tstep)
{
	return 0.5 * ((cx)one / (cx)(tstep * tstep * tstep)) * (g(t1 + tstep, t2 + tstep) - 2. * g(t1 + tstep, t2) + g(t1 + tstep, t2 - tstep) - g(t1 - tstep, t2 + tstep) + 2. * g(t1 - tstep, t2) - g(t1 - tstep, t2 - tstep));
}

//full third deriv wrt tr

cx d3dtr(double t1, double t2, double tstep)
{
	return 0.125 * (d3t1(t1, t2, tstep) - d3t2(t1, t2, tstep) + 3. * (d3t1t22(t1, t2, tstep) - d3t12t2(t1, t2, tstep)));
}

int main()
{

	clock_t tStart = clock();

	ofstream test, gee;
	test.open("testplotOLD.txt");
	gee.open("gtimedataOLD.txt");

	//calculate all the integrals we need one time, store them in a vector (first if block)
	//also write them to a file so we can use them quick

	//if we've already calculated them, set 1 to 0 to execute second block for fixed dt, tmax
	if (1)
	{
		ofstream integrals;
		integrals.open("integralsOLD.txt");

		cx holder(0, 0);

		for (int i = 1; i <= ntot; i += 1)
		{
			double t = -i * dt;
			holder += f(t) * dt;
			intf[ntot - i] = holder;
			//	cout<<t<< " "<<ntot-i<<" "<<intf[ntot-i]<<endl;
		}

		cx holder1(0, 0);

		for (int i = 0; i <= 2 * ntot; i++)
		{
			double t = i * dt - tmax;

			if (ntot <= i && i <= 2 * ntot)
			{
				holder1 += dt * f(t);

				intf[i] = holder1;
				//cout<<t<<" "<<intf[i]<<endl;
			}
			CC[i] = C(t);
			Uint[i] = Uin(t);
			cosi[i] = cx(0, ((double)g0 * 2. / ((double)sqrt(2 * m * w) * w)) * (CC[i]).real());
			//	cosi[i]=I*cosint(t);
			integrals << t << " " << CC[i].real() << " " << CC[i].imag() << " " << intf[i].real() << " " << intf[i].imag() << " " << cosi[i].real() << " " << cosi[i].imag() << endl;
		}
		integrals.close();
	}
	else
	{
		ifstream integrals;
		integrals.open("g1dt01tmax10cf.txt");
		for (int i = 0; i <= ntot; i++)
		{
			double cr = CC[i].real();
			double ci = CC[i].imag();
			double intfr = intf[i].real();
			double intfi = intf[i].imag();
			double cosir = cosi[i].real();
			double cosii = cosi[i].imag();
			integrals >> cr >> ci >> intfr >> intfi >> cosir >> cosii;
		}

		integrals.close();
	}

	cout << "Integrals done in " << ((double)(clock() - tStart) / CLOCKS_PER_SEC) << endl;

	double ta = tmax * 0.5;

	cx holdy = 0;

	vector<cx> green(ntot + 1);

	cout << avenf << endl;
	double Omega = 0.000;

	double testy = 0;
	ofstream file2;
	file2.open("datadampedOLD.txt");

	int halftot = 0.5 * ntot;
	for (int i = 3; i <= halftot - 3; i++)
	{

		double t = i * 2 * dt;

		test << t << " " << d3dtr(t, t, dt).real() << " " << thirdmomwrong(t) << " " << thirdmomwrong(t)-d3dtr(t,t,dt).real() << endl;

		//		test<<t<< " "<<1*(d3dtr(t,t,dt)).real()+0*thirdmom(t)<< " "<<0*d3dtr(t,t,dt).real()+thirdmomwrong(t)<<endl;

		//test<<t<<" "<<(-totalderiv(t,t,dt)).real()<<" "<<firstmom(t)<<endl;

		//	test<<t<<" "<<d2dtr2(t,t,dt).imag()<<" "<<secondmom(t)<<endl;
		//	test<<t<<" "<<gR(ta+0.5*t, ta-0.5*t).imag()<<endl;
		//file2<<t<<" "<<exp(-Omega*t)*gR(ta+0.5*t, ta-0.5*t).imag()<<endl;
		file2 << t << " " << d3dtr(t, t, dt).real() - thirdmomwrong(t) << endl;
		gee << t << " " << g(t) << endl;
		green[i] = gR(ta + 0.5 * t, ta - 0.0 * t);
	}

	double dw = 0.1;
	double wmax = 5;
	double wmin = -5;
	int nw = ((wmax - wmin) / dw);

	ofstream freq;
	freq.open("freqold.txt");

	ofstream fourier;
	fourier.open("fourierold.txt");

	vector<cx> deltaweights(nw + 1);

	for (int i = 0; i <= nw; i++)
	{
		double w = i * dw + wmin;
		holdy = 0;

		cx c0(0, 0);
		for (int j = 0; j <= ntot; j++)
		{
			double t = j * dt;
			double index = i;
			c0 += cx(((double)1) / (double)(8 * pi), 0) * dt * green[j] * exp(-I * index * 0.25 * t);
			holdy += dt * green[j] * exp(I * w * t) * exp(-Omega * t);
		}
		deltaweights[i] = c0;
		//cout<<deltaweights[i]<<endl;
		//cout<<c0<<endl;
		freq << w << " " << (-1. / pi) * holdy.imag() << endl;
	}
	fourier.close();
	freq.close();

	cout << "ay all done in " << ((double)(clock() - tStart) / CLOCKS_PER_SEC) << endl;
	;

	test.close();

	return 0;
}
