#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdio>
#include <algorithm>
//���������� ������� �����-������ � ����������� ����������
//������� ����� �� ����������� ��������� ���.41,52-55
//� ��������� "������������� ���������� � ������ ����������"

//����� ������������� ��� ���������� invI05
double a[5]=
{
    0.23960888E+0,
    0.25551970E-1,
    0.16001700E-2,
    0.62234299E-4,
    0.19789713E-5
};

double alpha[4]=
{
	0.72771980E+3,
    -0.11754360E+4,
    0.32688108E+3,
    -0.22333565E+1
};

double b1=-0.26354443E-1;
double beta1=0.33069903E+3;
double beta2=0.10575554E+1;

//����� ������������� ��� ���������� I15invI05
double aa[5]=
{
	1,
	0.50625059,
	0.93656560E-1,
	0.74482877E-2,
	0.21519676E-3
};

double bbb[3]=
{
	1,
	0.10730909,
	0.33951152E-2
};

double aal[4]=
{
	0.58432930E+5,
	0.21851397E+5,
	0.18917921E+4,
	0.57806497E+2
};

double bbt[2]=
{
	0.11348583E+4,
	0.41356686E+2
};

const double invI05argmin = 1e-161;

double ComKsi1(double Y)
{
	Y = std::max(invI05argmin, Y);
	return log(4*Y*Y*(1+a[0]*pow(Y,2)+a[1]*pow(Y,4)+
		+a[2]*pow(Y,6)+a[3]*pow(Y,8)+a[4]*pow(Y,10))/(3.*pow(M_PI,0.5)*(1+b1*Y*Y)));
}

double ComKsi2(double Y)
{
  	return pow( ( (alpha[0] + alpha[1]*pow(Y,(8./3.)) +alpha[2]*pow(Y,(16./3.))+
            alpha[3]*pow(Y,(24./3.)) + pow(Y,(32./3.)) ) / ( beta1+
			+beta2*pow(Y,(8./3.))+pow(Y,(16./3.)) )),0.25);
}


double invI05(double z)
{
	double y=sqrt(1.5*z);
	if(y<3.02) return ComKsi1(y);
	else return ComKsi2(y);
}

// page 48
double I15low(double y2)
{
	double hs,zn;
	hs=aa[0]+y2*(aa[1]+y2*(aa[2]+y2*(aa[3]+y2*aa[4])));
	zn=bbb[0]+y2*(bbb[1]+y2*bbb[2]);
	return y2*pow(hs/zn,1/3.0);
}

// page 48
double I15hi(double y)
{
	double y8d3=pow(y,8/3.0);
	double hs,zn;
	hs=aal[0]+y8d3*(aal[1]+y8d3*(aal[2]+y8d3*(aal[3]+y8d3)));
	zn=bbt[0]+y8d3*(bbt[1]+y8d3);
	return 0.4*y*y*pow(hs/zn,1/4.0);
}

// page 48
double I15invI05(double z)
{
	double y=sqrt(1.5*z);
	if(y<3.75) return I15low(1.5*z);
	else return I15hi(y);
}

//page  48 I05(mu/t)
double I05mu_d_t(double t,double V,double xe)
{
	return M_PI*M_PI*xe/(sqrt(2.0)*pow(t,1.5)*V);
}

//page 48 I15(mu/t)
double I15mu_d_t(double t,double V,double xe)
{
	return I15invI05(I05mu_d_t(t,V,xe));
}

// page 42
double mu(double t,double V,double xe)
{
	return invI05(M_PI*M_PI*xe / (sqrt(2.0)*pow(t, 1.5)*V)) * t;
}
