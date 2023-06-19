//============================================================================
// Name        : LocalHH.cpp
// Author      : Ryan Nesselrodt
// Version     :


//Description: This is a code for calculating the local atomic green's function for the 
//Holstien-Hubbard model with arbitrary time dependent couplings g(t),U(t) from derived expression

//In addition in this code we verify the self-consistency of the expression by comparing the first
//three spectral moments with their derived moment sum rule expressions for the local case
//We calculate these moment sum rules via expectation values with exact values, and calculate the 
//moments via three numerical derivatives of the expression here
//============================================================================


#include<iostream>
#include<fstream>
#include<complex>
#include<vector>
#include<cmath>


using namespace std;
typedef complex<double> cx;


// Definitions
const cx one(1,0);
const cx I(0,1);
double pi = 4*atan(1);


//define physical parameters
double mu=0.25;
double beta=1;
double m=1;

//freqency
double w=1;
double fw= ((double)1)/((double)(1-exp(-beta*w)));

double nw=((double)1)/((double)(exp(beta*w)-1.));

//COMPUTATIONAL PARAMS
//define number of rectangles in integrals


double DT=0.005;
double tmax=30;
double tmin=0;

int ntot= round((double)(tmax-tmin)/(double)DT);

vector<cx> CC(ntot+1);
vector<cx> intf(ntot+1);
vector<cx>cosi(ntot+1);
vector<cx>Uint(ntot+1);


//Define all the necessary functions to calculate the local nonequilibrium Holstien GF

//e-ph coupling
double g(double t){
//return 1-0.005*exp(-0.1*t);
//return 1;
//return 1+0.1*sin(0.001*t)*exp(-0.0*t);
return 1-0.1*sin(0.1*t)*exp(-0*t);
//return 1+0.1*t;
}

double dgdt(double t, double dt){
	return ((double)1./(double)(2*dt))*(g(t+dt)-g(t-dt));
	
}

double d2gdt(double t, double dt){
	return ((double)1/(double)(dt*dt))*(g(t+dt)-2*g(t)+g(t-dt));
}



//more definitions
double g0=g(tmin);


//if 0, g,U constant, this evaluates x(t) integrals exactly. otherwise set to 1
int eqflag=1;



//electron-electron coupling
double U(double t){
	//return 1-0.01*t*t;
	//return 1;
	return 1;
}
double U0=U(tmin);

double d2Udt(double t, double dt){
	return ((double)1/(double)(dt*dt))*(U(t+dt)-2*U(t)+U(t-dt));
}


//build partition function

double factor1= exp(beta*mu+((double)beta*g0*g0)/((double)2*m*w*w));
double factor2= exp(beta*(2.*mu-U0)+((double)2*beta*g0*g0)/((double)m*w*w));
double Z=fw*(1.+2.*factor1+factor2);
double Zm=(1.+2.*factor1+factor2);
double overZ= ((double)1/(double)Z);
double overZm=((double)1/(double)Zm);


//build ave nf


double avenf= overZm*(2.*factor1+2.*factor2);

double avenupndown=overZm*exp(beta*(2*mu-U0)+((double)2*beta*g0*g0)/((double)m*w*w));

double avencube = avenf+6.*avenupndown;





//FUNCTIONS NEEDED TO GET MOMENTS

//Define the function C that appears frequently in our expressions
//has one internal integration
std::complex<double> C(double t){
	//make coefficient
	std::complex<double> coeff((double)1/(double)(std::sqrt(2*m*w)),0);


	int num=round((double)(t-tmin)/(double)DT);
	std::complex<double> sum(0,0);
	
	for(int i=0; i<=num; i++){
		double tp=i*DT+tmin;
		sum += g(tp)*exp(std::complex<double>(0,w*(tp-tmin)))*DT;
	}
	
return coeff*sum;
}




//function for expectation value of x as function of time
double xave(double t, int flag){

	if(flag ==0){
		return -((double)g0*avenf/(double)(m*w*w));
	}
	else{

	int num = round((double)(t-tmin)/(double)(DT));
	cx sum(0,0);


	for(int i=0; i<=num; i++){
		double time= i*DT+tmin;
		sum += exp(I*w*(time))*g(time)*DT;
	}

	double cons=((double)1/(double)(m*w));

	cx number= I*exp(-I*w*(t))*sum*cons;
	double realio =number.real();

	double facy = ((double)g0/((double)(m*w*w)));


	return -avenf*(facy*cos(w*(t-tmin))+realio);
	}
}

double dxavedt(double t,double dt){
return ((double)(1.)/(double)(2*dt))*(xave(t+dt,eqflag)-xave(t-dt,eqflag));
}

double dgdxaveexplicit(double t, double dt){
	double cons = dgdt(t,dt)*avenf*w;
	
double term1=sin(w*t)*((double)g0/(double)(m*w*w));
cx holder(0,0);
int number=round((double)(t-tmin)/(double)dt);
for(int i=0; i<=number; i++){
	double time =tmin+i*dt;
	holder+= dt*g(time)*exp(I*w*time)*((double)1/(double)(m*w));
}
double term2= (exp(-I*w*t)*holder).real();
return cons*(term1-term2);
}
//average <n_bsigma x(t)>
double avenx(double t){
	int num = round((double)(t-tmin)/(double)(DT));
	cx sum(0,0);


	for(int i=0; i<=num; i++){
		double time= i*DT+tmin;
		sum += exp(I*w*(time-tmin))*g(time)*DT;
	}

	double cons=((double)1/(double)(m*w));

	cx number= I*exp(-I*w*(t-tmin))*sum*cons;
	double realio =number.real();

	double facy = ((double)g0/((double)(m*w*w)));



	return -(0.5*avenf+avenupndown)*(facy*cos(w*(t-tmin))+realio);
}

//average <x^2(t)>
double avex2(double t){
	int num = round((double)(t-tmin)/(double)(DT));
	cx sum(0,0);


	for(int i=0; i<=num; i++){
		double time= i*DT+tmin;
		sum += exp(I*w*(time-tmin))*g(time)*DT;
	}

	double cons=((double)1/(double)(m*w));

	cx number= I*exp(-I*w*(t-tmin))*sum*cons;
	double realio =number.real();

	double facy = ((double)g0/((double)(m*w*w)));


	return (avenf+2*avenupndown)*(facy*cos(w*(t-tmin))+realio)*(facy*cos(w*(t-tmin))+realio)+((double)1/(double)(2*m*w*tanh(0.5*beta*w)));
}

//average <n_bsig x^2(t)>
double avenx2(double t){
	int num = round((double)(t-tmin)/(double)(DT));
	cx sum(0,0);


	for(int i=0; i<=num; i++){
		double time= i*DT+tmin;
		sum += exp(I*w*(time-tmin))*g(time)*DT;
	}

	double cons=((double)1/(double)(m*w));

	cx number= I*exp(-I*w*(t-tmin))*sum*cons;
	double realio =number.real();

	double facy = ((double)g0/((double)(m*w*w)));


	return (0.5*avenf+3*avenupndown)*(facy*cos(w*(t-tmin))+realio)*(facy*cos(w*(t-tmin))+realio)+((double)1/(double)(2*m*w*tanh(0.5*beta*w)))*0.5*avenf;
}

//variance of x operator, <x^2>-<x>^2
double varx(double t, int flag){
	double consty= avenf-avenf*avenf+2*avenupndown;
	if(flag ==0){
		return xave(t,flag)*xave(t,flag)*consty+((double)1/(double)(2*m*w*tanh(0.5*beta*w)));
	}
	else{
	int num = round((double)(t-tmin)/(double)(DT));
	cx sum(0,0);


	for(int i=0; i<=num; i++){
		double time= i*DT+tmin;
		sum += exp(I*w*(time))*g(time)*DT;
	}

	double cons=((double)1/(double)(m*w));

	cx number= I*exp(-I*w*(t))*sum*cons;
	double realio =number.real();

	double facy = ((double)g0/((double)(m*w*w)));

	double inter=(facy*cos(w*(t-tmin))+realio);

	return consty*inter*inter+((double)1/(double)(2*m*w*tanh(0.5*beta*w)));}

}


//return variance of x^3 operator, <x^3>-<x>^3

double varxcube(double t, int flag){
	double cothb =((double)1/(double)tanh(0.5*beta*w));

	if(flag==0){

		return (avenf*(1.-avenf*avenf)+6*avenupndown)*xave(t,flag)*xave(t,flag)*xave(t,flag)+((double)1/(double)(m*w))*1.5*cothb*xave(t,flag);
	}
	else{
	//build the real integral
	int num = round((double)(t-tmin)/(double)(DT));
	cx sum(0,0);


	for(int i=0; i<=num; i++){
		double time= i*DT+tmin;
		sum += exp(I*w*(time-tmin))*g(time)*DT;
	}

	double cons=((double)1/(double)(m*w));

	cx number= I*exp(-I*w*(t-tmin))*sum*cons;
	double realio =number.real();


	//other constants
	double a1=((double)(g0)/(double)(m*m*w*w*w));
	double a2= ((double)1/(double)(m*w));

	double a3= (double)g0/(double)(m*w*w);

	

	//term proportional to <n>
	double term2=-avenf*cothb*1.5*a2*(a3*cos(w*(t-tmin))+realio);
	//double term2=1.5*a2*cothb*xave(t);

	//term proportional to <n^3>-<n>^3
	double term3 =-(avenf*(1.-avenf*avenf)+6*avenupndown)*(a3*cos(w*(t-tmin))+realio)*(a3*cos(w*(t-tmin))+realio)*(a3*cos(w*(t-tmin))+realio);
	return term2+term3;}
}


//FUNCTIONS TO COMPUTE MOMENTS OF GF


//return first moment of the GF
double firstmom(double t){
	return -mu+g(t)*xave(t,eqflag)+U(t)*0.5*avenf;
}
//return second moment of the GF

double secondmom(double t){
	return (-mu+g(t)*xave(t,eqflag)+0.5*U(t)*avenf)*(-mu+g(t)*xave(t,eqflag)+0.5*U(t)*avenf)+g(t)*g(t)*(varx(t,eqflag))+U(t)*U(t)*(0.5*avenf-0.25*avenf*avenf)+2.*U(t)*g(t)*(avenx(t)-xave(t,eqflag)*0.5*avenf);
}

//extra terms for third moment
double extraterms(double t){
	return U(t)*U(t)*U(t)*(0.5*avenf-0.125*avenf*avenf*avenf)-3*mu*U(t)*U(t)*(0.5*avenf-0.25*avenf*avenf)+(3*U(t)*U(t)*g(t))*(avenx(t)-0.25*avenf*avenf*xave(t,eqflag))-6*mu*U(t)*g(t)*(avenx(t)-0.5*avenf*xave(t,eqflag))+3*U(t)*g(t)*g(t)*(avenx2(t)-0.5*avenf*xave(t,eqflag)*xave(t,eqflag));
}



//third moment of the GF
double thirdmom(double t){
	//moment as derived from commutators
	return firstmom(t)*firstmom(t)*firstmom(t)+g(t)*w*w*xave(t,eqflag)+g(t)*g(t)*g(t)*varxcube(t,eqflag)-3.*mu*g(t)*g(t)*varx(t,eqflag)+g(t)*g(t)*0.5*((double)1/(double)m)*(2*avenf+1)-0.25*(d2gdt(t,DT)*xave(t,eqflag)+d2Udt(t,DT))+extraterms(t)-0.5*dgdt(t,DT)*dxavedt(t,DT);
//moment as analytic third deriv
	//return firstmom(t)*firstmom(t)*firstmom(t)+0.25*g(t)*w*w*xave(t,eqflag)+g(t)*g(t)*g(t)*varxcube(t,eqflag)-3.*mu*g(t)*g(t)*varx(t,eqflag)+g(t)*g(t)*0.5*((double)1/(double)m)*(-0.5*avenf+1.)-0.25*(d2gdt(t,DT)*xave(t,eqflag)+d2Udt(t,DT))+extraterms(t)+0.5*dgdt(t,DT)*dxavedt(t,DT);

}

double thirdmomwrong(double t){
//moment as analytic third deriv
	return firstmom(t)*firstmom(t)*firstmom(t)+0.25*g(t)*w*w*xave(t,eqflag)+g(t)*g(t)*g(t)*varxcube(t,eqflag)-3.*mu*g(t)*g(t)*varx(t,eqflag)+g(t)*g(t)*0.5*((double)1/(double)m)*(-0.5*avenf+1.)-0.25*(d2gdt(t,DT)*xave(t,eqflag)+d2Udt(t,DT))+extraterms(t)+0.5*dgdt(t,DT)*dxavedt(t,DT);	
}




//FUCTIONS TO COMPUTE THE GF


//define f(t) the prefactor with an integral over sine
double f(double t){

	int num =round((double)(t-tmin)/(double)DT);

	double sum=0;

	for(int i=0; i<=num; i++){
		double tp=i*DT+tmin;
		sum+=DT*g(tp)*sin(w*(t-tp));
	}

	return g(t)*sum;
}






double cosint(double t){

	int num=round((double)(t-tmin)/(double)DT);
	double sum=0;
	
	for(int i=0; i<=num; i++){
		double tp=i*DT+tmin;
		sum += g(tp)*cos(w*(tp-tmin))*DT;
	}
	
	double coeff=((double)g0)/((double)m*w*w);

return coeff*sum;
}

//integral over U
std::complex<double> Uin(double t){

	int num=round((double)(t-tmin)/(double)DT);

	std::complex<double> sum(0,0);
	
	for(int i=0; i<=num; i++){
		double tp=i*DT+tmin;
		sum += U(tp)*DT;
	}
	
return -I*sum;
}



//define greater and lesser green's functions by accessing elements from vectors we've built



cx Ggreater(double t1,double t2){

	if(t1<tmin ||t2<tmin){
		cx fail = cx(10000.,10000.);
		return fail;
	}
	else{
	int nt1= round((double)(t1-tmin)/(double)DT);
	int nt2=round((double)(t2-tmin)/(double)DT);


	//define quantites that go in exponents
	double c=0.5-fw;
	cx  mag=CC[nt1]*conj(CC[nt1])+CC[nt2]*conj(CC[nt2]);

	cx imy= -1.*(CC[nt1]*conj(CC[nt2]));
	double fac=2.*fw*(CC[nt1]*conj(CC[nt2])).real();


	cx magi =abs(CC[nt1]-CC[nt2]);
	std::complex<double> coeff(0,(double)1/(double)(2*m*w));

	cx Qew= coeff*(intf[nt1]-intf[nt2]);


	//define prefactor
	cx prefactor=-I*overZm*exp(I*mu*(t1-t2)+one*0.5*(mag)+imy-fw*magi*magi);

	//second term in sum
	cx factor2=exp(one*beta*mu+one*((double)(g0*g0*beta)/((double)2*m*w*w))+one*3.*Qew +cosi[nt1]-cosi[nt2]+Uint[nt1]-Uint[nt2]);

	return prefactor*(exp(Qew)+factor2);}
}




cx Glesser(double t1,double t2){
if(t1<tmin ||t2<tmin){
	cx fail = cx(100000.,100000.);
		return fail;
	}
	else{
	int nt1= round((double)(t1-tmin)/(double)DT);
	int nt2=round((double)(t2-tmin)/(double)DT);

	std::complex<double> coeff(0,(double)1/(double)(2*m*w));

	cx Qew= coeff*(intf[nt1]-intf[nt2]);

	cx magi =abs(CC[nt1]-CC[nt2]);


	//define quantites that go in exponents
	double c=0.5-fw;
	cx mag=CC[nt1]*conj(CC[nt1])+CC[nt2]*conj(CC[nt2]);

	cx imy= -1.*(CC[nt2]*conj(CC[nt1]));
	double fac=2.*fw*(CC[nt1]*conj(CC[nt2])).real();


	//define prefactor
	cx prefactor=I*overZm*exp(I*mu*(t1-t2)+one*0.5*(mag)+imy-fw*magi*magi);

	//first term in sum
	cx factor1=exp(one*beta*mu+one*((double)(g0*g0*beta)/((double)2*m*w*w))+one*Qew+cosi[nt1]-cosi[nt2]);

	//second term in sum

	cx factor2=exp(one*2.*beta*mu-one*beta*U0+2.*one*((double)(g0*g0*beta)/((double)m*w*w))+3.*one*Qew+2.*(cosi[nt1]-cosi[nt2])+Uint[nt1]-Uint[nt2]);

	return prefactor*(factor1+factor2);}
}


double theta1(double t1, double t2){
	if(t1>=t2){
		return 1;
	}
	else
		return 0;
}

double theta0(double t1, double t2){
	if(t1>t2){
		return 1;
	}
	else
		return 0;
}
cx gR(double t1, double t2){
	//define quantities that go in exponents

	return (theta1(t1,t2)*Ggreater(t1,t2)-theta1(t1,t2)*Glesser(t1,t2));
}

cx g(double t1, double t2){
	return Ggreater(t1,t2)-Glesser(t1,t2);
}

cx gA(double t1,double t2){
	return theta0(t2,t1)*(Glesser(t1,t2)-Ggreater(t1,t2));
}



cx gequil(double t){
	double cons1 =((double) g0*g0)/(double)(2*m*w*w);
	double cons2 =((double) g0*g0)/(double)(m*w*w*w);
	//cx prefactor =-I*theta1(t,0)*overZm*exp(I*(mu+cons1)*t-cons2*(nw+0.5)*(1.-cos(w*t)));
	cx prefactor =-I*overZm*exp(I*(mu+cons1)*t-cons2*(nw+0.5)*(1.-cos(w*t)));
	cx term1=exp(-I*sin(w*t)*0.5*cons2)*(1.+exp(beta*(mu+cons1)+I*(cons1-U0)*t));
	cx term2=exp(I*0.5*cons2*sin(w*t))*(exp(beta*(mu+cons1))+exp(beta*(4.*cons1+2.*mu-U0)+I*(2.*cons1-U0)*t));
	return prefactor*(term1+term2);
}




//FUNCTIONS TO COMPUTE NUMERICAL DERIVATIVES FOR MOMENTS

//first deriv
cx deriv1(double t1, double t2,double tstep){

	double realy = ((double)(g(t1+tstep,t2).real()-g(t1,t2).real())/((double)tstep));
	double imagi =((double)(g(t1+tstep,t2).imag()-g(t1,t2).imag())/((double)tstep));
	return one*realy+I*imagi;
}

cx deriv2(double t1, double t2, double tstep){
	double realy = ((double)(g(t1,t2).real()-g(t1,t2-tstep).real())/((double)tstep));
		double imagi =((double)(g(t1,t2).imag()-g(t1,t2-tstep).imag())/((double)tstep));
		return one*realy+I*imagi;
	}
cx totalderiv(double t1, double t2, double tstep){
	return 0.5*(deriv1(t1,t2,tstep)-deriv2(t1,t2,tstep));
}



//functions for calculating second derivative

//pieces needed
cx d2t1(double t1,double t2, double tstep){
	return ((cx)(g(t1+tstep,t2)-2.*g(t1,t2)+g(t1-tstep,t2))/((cx)tstep*tstep));
}

cx d2t2(double t1,double t2, double tstep){
	return ((cx)(g(t1,t2+tstep)-2.*g(t1,t2)+g(t1,t2-tstep))/((cx)tstep*tstep));
}

cx d2t1t2(double t1, double t2, double tstep){
	return 0.25*((cx)(g(t1+tstep,t2+tstep)-g(t1+tstep,t2-tstep)+g(t1-tstep,t2-tstep)-g(t1-tstep,t2+tstep))/((cx)tstep*tstep));
}

//total second deriv wrt t2
cx d2dtr2(double t1, double t2, double tstep){
	return 0.25*(d2t1(t1,t2,tstep)+d2t2(t1,t2,tstep)-2.*d2t1t2(t1,t2,tstep));
}


//functions for calculating third deriv numerically

//components required

cx d3t1(double t1, double t2, double tstep){
	return ((cx)one/(cx)(tstep*tstep*tstep))*(-0.5*g(t1-2*tstep,t2)+g(t1-tstep,t2)-g(t1+tstep,t2)+0.5*g(t1+2*tstep,t2));
}

cx d3t2(double t1, double t2, double tstep){
	return ((cx)one/(cx)(tstep*tstep*tstep))*(-0.5*g(t1,t2-2*tstep)+g(t1,t2-tstep)-g(t1,t2+tstep)+0.5*g(t1,t2+2*tstep));
}


cx d3dt1(double t1, double t2, double tstep){
	return 0.125*((cx)one/(cx)(tstep*tstep*tstep))*(g(t1+3*tstep,t2)-3.*g(t1+tstep,t2)+3.*g(t1-tstep,t2)-g(t1-3*tstep,t2));
}

cx d3dt2(double t1, double t2, double tstep){
	return 0.125*((cx)one/(cx)(tstep*tstep*tstep))*(g(t1,t2+3*tstep)-3.*g(t1,t2+tstep)+3.*g(t1,t2-tstep)-g(t1,t2-3*tstep));
}





cx d3t12t2(double t1, double t2, double tstep){
	return 0.125*((cx)one/(cx)(tstep*tstep*tstep))*(g(t1+2*tstep,t2+tstep)-g(t1+2*tstep,t2-tstep)-2.*g(t1,t2+tstep)+2.*g(t1,t2-tstep)+g(t1-2*tstep,t2+tstep)
			-g(t1-2*tstep,t2-tstep));
}

cx d3t1t22(double t1, double t2, double tstep){
	return 0.125*((cx)one/(cx)(tstep*tstep*tstep))*(g(t1+tstep,t2+2*tstep)-g(t1-tstep,t2+2*tstep)-2.*g(t1+tstep,t2)+2.*g(t1-tstep,t2)+g(t1+tstep,t2-2*tstep)
			-g(t1-tstep,t2-2*tstep));
}

//d^3g/dt_1^2dt2

cx d3dt11dt2(double t1, double t2, double tstep){
	return 0.5*((cx)one/(cx)(tstep*tstep*tstep))*(g(t1+tstep,t2+tstep)-2.*g(t1,t2+tstep)+g(t1-tstep,t2+tstep)
			-g(t1+tstep,t2-tstep)+2.*g(t1,t2-tstep)-g(t1-tstep,t2-tstep));
}

//d^3g/dt1dt2^2

cx d3dt1dt22(double t1, double t2, double tstep){
	return 0.5*((cx)one/(cx)(tstep*tstep*tstep))*(g(t1+tstep,t2+tstep)-2.*g(t1+tstep,t2)+g(t1+tstep,t2-tstep)
			-g(t1-tstep,t2+tstep)+2.*g(t1-tstep,t2)-g(t1-tstep,t2-tstep));
}


//full third deriv wrt tr

cx d3dtr(double t1, double t2, double tstep){
	return 0.125*(d3t1(t1,t2,tstep)-d3t2(t1,t2,tstep)+3.*(d3t1t22(t1,t2,tstep)-d3t12t2(t1,t2,tstep)));
}



//Equilibrium third deriv
cx equild3dt(double t, double tstep){
	return 0.125*((cx)one/(cx)(tstep*tstep*tstep))*(gequil(t+3.*tstep)-3.*gequil(t+tstep)+3.*gequil(t-tstep)-gequil(t-3*tstep));
}

//series third deriv (small tr)
cx seriesthird(double t1, double t2, double tstep){
	double ta=0.5*(t1+t2);
	double tr=t1-t2;
	return ((double)6/((double)tstep*tstep*tstep))*(g(t1,t2)-(g(ta,ta)+totalderiv(ta,ta,tstep)*tr+0.5*tr*tr*d2dtr2(ta,ta,tstep)));
}

int main(){

	ofstream test,gee;
	test.open("testplot.txt");
	gee.open("gtimedata.txt");

//calculate all the integrals we need one time, store them in a vector (first if block)
//also write them to a file so we can use them quick

//if we've already calculated them, set 1 to 0 to execute second block for fixed dt, tmax
//calculate from tmin to tmax
if(1){
	ofstream integrals;
	integrals.open("integrals.txt");

	cx holder(0,0);

	for(int i=0; i<=ntot; i++){
		double t=i*DT+tmin;
		holder+=DT*f(t);
		intf[i]=holder;
			
		CC[i]=C(t);
		Uint[i]=Uin(t);
		cosi[i] =cx(0,((double)g0*2./((double)sqrt(2*m*w)*w))*(CC[i]).real());
		integrals<<t<<" "<<CC[i].real()<<" "<<CC[i].imag()<<" "<<intf[i].real()<<" "<<intf[i].imag()<<" "<<cosi[i].real()<<" "<<cosi[i].imag()<<endl;
	}
	integrals.close();

}
else{
	ifstream integrals;
	integrals.open("g1dt01tmax10cf.txt");
	for(int i=0; i<=ntot; i++){
		double cr=CC[i].real();
		double ci=CC[i].imag();
		double intfr=intf[i].real();
		double intfi=intf[i].imag();
		double cosir=cosi[i].real();
		double cosii=cosi[i].imag();
		integrals >> cr>>ci>>intfr>>intfi>>cosir>>cosii;}
	
	integrals.close();
}




double ta=0.5*tmax;



cx holdy=0;

vector<cx> green(ntot+1);


double Omega=0.01;

double testy=0;
ofstream file2;
file2.open("datadamped.txt");
ofstream x2;
x2.open("x2.txt");
ofstream moments;
moments.open("momentsright.txt");
//want start index at time =0, stop index at tave+0.5t=tmax;
int starti= ((double)(-tmin)/(double)DT);

cout<< "average n_f= "<< avenf<<endl;
int halftot= 0.5*ntot;

int stopi =-((double)(tmax-ta)/DT)+ntot;
int other=0.5*stopi;
cout<<other;
	for(int i=3; i<halftot;i++){
		double t=i*2*DT+tmin;
		//double t=i*DT+tmin;
		if (i==stopi){
			cout<<t<<" "<<tmax<<" "<<ta-0.5*t<<endl;
		}

	test<<t<<" "<<0.75*g(t)*w*w*xave(t,eqflag)+0.75*g(t)*g(t)*avenf*((double)1/(double)m)<<" "<<dgdxaveexplicit(t,DT)<< " "<<dgdt(t,DT)*dxavedt(t,DT)<<endl;
	file2<<t<<" "<<exp(-Omega*t)*gR(ta+0.5*t, ta-0.5*t).imag()<<endl;
	gee<<t<<" "<<g(t)<<endl;
	x2<<t<< " "<<avex2(t)<<endl;
	//moments<<t<< " "<<d3dtr(t,t,DT).real()<<" "<<thirdmom(t)<<endl;
	//moments<<t<<" "<<seriesthird(ta+0.5*t,ta-0.5*t,DT).real()<< " "<<thirdmom(t)<<endl;
	moments<<t<< " "<<d3dtr(t,t,DT).real()-1*thirdmom(t)<<" "<<0*d3dtr(t,t,DT).real()+1*thirdmom(t)<<endl;
	//moments<<t<<" "<<d3dtr(0.5*t+ta,-0.5*t+ta,dt).real()<<" "<<equild3dt(t,dt).real()<<endl;
	green[i]=gR(ta+0.5*t,ta-0.5*t);

	}

moments.close();
double dw=0.01;
double wmax=5;
double wmin=-5;
int nw= ((wmax-wmin)/dw);

ofstream freq;
freq.open("freq.txt");

ofstream fourier;
fourier.open("fourier.txt");

vector<cx> deltaweights(nw+1);


for(int i =0; i<=nw; i++){
	double w=i*dw+wmin;
	holdy=0;
	
	cx c0(0,0);
		for (int j=0; j<halftot; j++){
			double t=j*DT+tmin;
			double index =i;
			c0 += cx(((double)1)/(double)(8*pi),0)*DT*green[j]*exp(-I*index*0.25*t);
			holdy+=DT*green[j]*exp(I*w*t)*exp(-Omega*t);
		}
		deltaweights[i] = c0;
		//cout<<deltaweights[i]<<endl;
		//cout<<c0<<endl;
	freq<<w<< " "<<(-1./pi)*holdy.imag()<<endl;}
fourier.close();
freq.close();
gee.close();
x2.close();
file2.close();



	cout<<"ay"<<endl;

	test.close();

	return 0;
}


