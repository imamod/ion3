#include "sahasolver.h"
#include "mu_fi.h"
#include "atom_ed.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <cstdio>

SahaSolver::SahaSolver(const TElement &element):_element(element), _xe(0)
{
    _x.resize(element.Z+1);
    _temp.resize(element.Z+1);
}

double SahaSolver::calculate_TVae(double T, double V)
{
    double a = -745, b = log(_element.Z), c;
    double fa, fb, fc;
    fa = ff(exp(a), T, V);
    fb = ff(exp(b), T, V);

    do
    {
        c = 0.5*(a + b);
        fc = ff(exp(c), T, V);
        if(((fa <= 0) && (fc >= 0)) || ((fa >= 0) && (fc <= 0)))
        {
            b = c;fb = fc;
        }
        else if(((fb <= 0) && (fc >= 0)) || ((fb >= 0) && (fc <= 0)))
        {
            a = c;fa = fc;
        }
    }
    while(b - a > 1e-6);

    _xe = exp(0.5*(a+b));
    return _xe;
}

double SahaSolver::calculate_lgTeV_lgVae(double lgT, double lgV)
{
    return calculate_TVae(pow(10.0,lgT) / eFi, pow(10.0,lgV));
}

double SahaSolver::ff(double xe, double T, double V)
{
    _xe = xe;
    double vFree = Vfree(V);
    if(vFree < 0) return std::numeric_limits<double>::max();

    double _mu = mu(T, vFree, _xe);
    double _P = p(T,vFree);

    _temp[0] = 0;
    for(int i = 1; i <= _element.Z; i++)
    {
        _temp[i] =  (-_element.fi[i-1] - _mu - _P * (_element.v[i] - _element.v[i-1])) / T;
    }
    _temp[_element.Z] -= log(2.0);

    double sum = 0, maxTemp = 0;
    for(int i = 0; i <= _element.Z; i++)
    {
        sum += _temp[i];_temp[i] = sum;
        if(sum > maxTemp) maxTemp = sum;
    }

    double expTemp1, Asum = 0, Bsum = 0;
    for(int i = 0; i <= _element.Z; i++)
    {
        expTemp1 = exp(_temp[i] - maxTemp);
        Asum += i * expTemp1;
        Bsum += expTemp1;
    }

    return Asum / Bsum - xe;

}

void SahaSolver::formX()
{

}

double SahaSolver::Vfree(double V)
{
    int i = floor(_xe);
    double fracXe = _xe - i;
    if(i < _element.Z)
    {
        return V - (_element.v[i] * (1-fracXe) + _element.v[i+1] * fracXe);
    }
    return V - _element.v[_element.Z];
}

double SahaSolver::P(double T, double V)
{
    return P(T, Vfree(V));
}

double SahaSolver::p(double T, double vFree)
{
    return 2*sqrt(2)/(3*M_PI*M_PI)*pow(T,2.5)*I15mu_d_t(T,vFree,_xe)+T/vFree;
}

double SahaSolver::E(double T, double V)
{
    return 0;
}

double SahaSolver::xe()
{
    return _xe;
}
