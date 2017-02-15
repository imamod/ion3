#include "sahasolver.h"
#include "mu_fi.h"
#include "atom_ed.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <cstdio>
#include <stdexcept>

SahaSolver::SahaSolver(const TElement &element):_element(element), _xe(0)
{
    _x.resize(element.Z+1);
    _H0.resize(element.Z+1);
}

SahaPoint SahaSolver::Calculate_TVae(double T, double V)
{
    double a = -746, b = log(double(_element.Z)), c;
    double fa, fb, fc;
    fa = ff(exp(a), T, V);
    fb = ff(exp(b), T, V);
	if(((fa >= 0) && (fb >= 0)) || ((fa <= 0) && (fb <= 0)))
	{
		if (fabs(fa) < fabs(fb)) _xe = exp(a);
		else _xe = exp(b);
	}
	else
	{
		do
		{
			c = 0.5*(a + b);
			fc = ff(exp(c), T, V);
			if (((fa <= 0) && (fc >= 0)) || ((fa >= 0) && (fc <= 0)))
			{
				b = c; fb = fc;
			}
			else if (((fb <= 0) && (fc >= 0)) || ((fb >= 0) && (fc <= 0)))
			{
				a = c; fa = fc;
			}
			else
			{
				char message[256];
				sprintf(message, "ln(a) = %g ln(b) = %g ln(c) = %g fa = %g fb = %g fc = %g\n", a, b, c, fa, fb, fc);
				error("find root error", message, T, V);

			}
		}
		while (b - a > 1e-7);
		_xe = exp(0.5*(a + b));
	}

	if (!isfinite(_xe))
	{
		char message[256];
		sprintf(message, "xe=%g\n",_xe);
		error("invalid xe", message, T, V);
	}

    formX(T, V);

    double vFree = Vfree(V);
    SahaPoint result;
    double E = e(T,vFree);
    double S = s(T,vFree);
    double P = p(T,vFree);

    result.Z = _element.Z;
    result.T = T;
    result.V = V;
    result.E = E;
    result.P = P;
    result.S = S;
    result.Xe = _xe;
    result.M = mu(T, vFree, _xe);
    result.F = E - T * S;
    result.K = 1 - vFree / V;

    return result;
}

SahaPoint SahaSolver::Calculate_lgTeV_lgVae(double lgT, double lgV)
{
    return Calculate_TVae(pow(10.0,lgT) / eFi, pow(10.0,lgV));
}

double SahaSolver::ff(double xe, double T, double V)
{
    _xe = xe;
    double vFree = Vfree(V);
	if (vFree < 0)
	{
		return std::numeric_limits<double>::max();
	}

    double maxH0;
    formH0(mu(T, vFree, _xe), p(T, vFree), T, maxH0);

    double expTemp1, Asum = 0, Bsum = 0;
    for(unsigned int i = 0; i <= _element.Z; i++)
    {
        //Вычитая maxH0, выполянем деление на макс. слагаемое знаменателя, см. ниже в formH0
        expTemp1 = exp(_H0[i] - maxH0); // что за формулы в данном цикле, что высчитываем? - см. ниже
        Asum += i * expTemp1;
        Bsum += expTemp1;
    }

    return Asum / Bsum - xe; // похоже на формулу на стр 42, где считается хе, но почему здесь вычитаем хе
    // Потому что это уравнение относительно xe, и чтобы его решить надо все перегнать все в левую часть,
    // чтобы было F(xe) = 0, а потом решать.

}

void SahaSolver::error(const std::string & errorType, const std::string &message, double T, double V)
{
	char result[1024];
	sprintf(result, "Internal Saha Error\nError type:%s\nlgT(eV) = %g lgV(a.e.) = %g Z = %d\nInfo:%s\n",
		errorType.c_str(), log10(T) + log10(eFi), log10(V), _element.Z, message.c_str());
	throw std::runtime_error(result);
}
// не понятен смысл данной функции
// тут вся физика и сидит :)
void SahaSolver::formH0(double mu, double P, double T, double &maxH0)
{
    _H0[0] = 0;
    for(unsigned int i = 1; i <= _element.Z; i++)
    {
      // стр 2 из src/theorys
        _H0[i] =  (-_element.fi[i-1] - mu - P * (_element.v[i] - _element.v[i-1])) / T;
    }
	_H0[_element.Z] -= log(2.0); // Что такое H0? Это логарифм из формулы (7) на стр 2 теории? Почему вычитаем логарифм 2?
    //Лучше посмотрите на формулу (25) из реферата и ту, что прямо под ней
    //Тут _H0 это fk(xe)
    //Там в роли DFk берем P * (_element.v[i] - _element.v[i-1]), "круглого" DFk нет совсем
    //а Iinv1/2(...) - это и есть mu, при этом t=T
    //Все Gk равны 2, кроме последнего, который = 1, из-за этого приходится делать это вычитание

    double sum = 0;
    maxH0 = 0;
    for(unsigned int i = 0; i <= _element.Z; i++)
    {
        sum += _H0[i];_H0[i] = sum;
        if(sum > maxH0) maxH0 = sum;
    }
    //В этом цикле считаются куммулятивные суммы, см. опять же (25) из реферата.
    //Здесь _H0 уже куммулятивная сумма fk(xe)
    //Вычисление максимума нужно для борьбы с переполнением при расчете.
    //Для этого числитель и знаменатель дроби (25) делится на максимальное слагаемое из тех, что
    //в знаменателе, что происходит уже в функции ff
}

void SahaSolver::formX(double T, double V)
{
    double maxH0, vFree = Vfree(V);
    formH0(mu(T, vFree, _xe), p(T,vFree), T, maxH0);

    double expSum = 0;
    for(unsigned int i = 0; i <= _element.Z; i++)
    {
        expSum += exp(_H0[i] - maxH0);
    }

    double logX0 = -log(expSum) - maxH0;

    for(unsigned int i = 0; i <= _element.Z; i++)
    {
        _x[i] = exp(logX0 + _H0[i]); // что такое _x[i]?
        // это массив концентраций ионов различной кратности
    }
}

double SahaSolver::Vfree(double V)
{
    unsigned int i = floor(_xe); // _xe - это концентрация здесь? что показывает i?
    // _xe - электронная концентрация. Плохо, что это глобальная переменная,
    // надо бы убрать потом и сделать локальной
    // i - целая часть _xe, что и так, наверное, понятно
    // Содержимое этой функции нигде не описано. Суть происходящего - приближенное вычисление
    // среднего объема ионов, а затем свободного объема. Собственно из-за этой вот "приближенности" весь этот
    // расчет не вполне точен и надо использовать метод Ньютона, чем Вы, собственно, и займетесь.
    double fracXe = _xe - i;
    if(i < _element.Z)
    {
        return V - (_element.v[i] * (1-fracXe) + _element.v[i+1] * fracXe);
    }
    return V - _element.v[_element.Z];
}

double SahaSolver::Vion(double rCoeff)
{
	double V = 0;
	for (int i = 0; i < _element.Z; i++)
	{
		double vi = 4 / 3.0 * M_PI * pow(rCoeff * (i + 1) / _element.fi[i], 3.0);
		V += _x[i] * vi;
	}
	return V;
}

double SahaSolver::p(double T, double vFree)
{
    return 2*sqrt(2.0)/(3*M_PI*M_PI) * pow(T,2.5) * I15mu_d_t(T,vFree,_xe) + T / vFree;
}

double SahaSolver::s(double T, double vFree)
{
    double Se = 0;
    if(_xe > 0)
    {
        Se = sqrt(2.0) / (M_PI*M_PI) * pow(T, 1.5) * vFree * (5.0 / 3.0 * I15mu_d_t(T, vFree, _xe) - mu(T,vFree,_xe) / T * I05mu_d_t(T, vFree, _xe));
    }

    const double M = 1822.887 * _element.A;
    double Si = 2.5, logG = log(2.0);
    for (unsigned int i = 0; i <= _element.Z; i++)
    {
       if (_x[i] > 0)
       {
          if(i == _element.Z) logG = 0;
          Si += _x[i] * (logG + log(vFree) - log(_x[i]) + 1.5 * log(M * T / 2.0 / M_PI));
       }
    }
    return Si + Se;
}

double SahaSolver::e(double T, double vFree)
{
    double Ee = sqrt(2.0)/(M_PI*M_PI) * pow(T,2.5) * vFree * I15mu_d_t(T,vFree,_xe);

    double Efi = 0;
    for(unsigned int i = 1; i <= _element.Z; i++)
    {
        Efi += _element.cumFi[i-1] * _x[i];
    }
    return 1.5*T + Ee  + Efi;
}

void SahaSolver::GetX(std::vector<double> &x)
{
    x = _x;
}
