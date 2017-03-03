#include "sahasolver.h"
#include "mu_fi.h"
#include "atom_ed.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <cstdio>
#include <stdexcept>

SahaSolver::SahaSolver(const TElement &element):_element(element)
{
    _x.resize(element.Z+1);
    _H0.resize(element.Z+1);
}

double SahaSolver::findroot(double logA, double logB, const std::function<double (double)> &F, double eps, double T, double V)
{
    double a = logA, b = logB, c;
    double fa, fb, fc;
    double root = 0;

    int cc = 0;

    fa = F(exp(a));
    fb = F(exp(b));
    if(((fa >= 0) && (fb >= 0)) || ((fa <= 0) && (fb <= 0)))
    {
        if (fabs(fa) < fabs(fb)) root = exp(a);
        else root = exp(b);
    }
    else
    {
        do
        {
            c = 0.5*(a + b);
            fc = F(exp(c));
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
            cc++;
        }
        while ((b - a > eps) && (cc < 60));
        root = exp(0.5*(a + b));
    }

    double root2 = exp(a) - (exp(b) - exp(a))/(fb-fa)*fa;
    if(fabs(F(root2)) < fabs(F(root)))
    {
        return root2;
    }
    else return root;
}

void SahaSolver::calcCore1(double T, double V, calcCoreResult &result)
{
    double xe = findroot(-746, log(double(_element.Z)), [&](double x) {return ff(x, T, V);}, 1e-7, T, V);

    if (!isfinite(xe))
    {
        char message[256];
        sprintf(message, "xe=%g\n",xe);
        error("invalid xe", message, T, V);
    }

    result.vFree = Vfree(V, xe);
    formX(T, V, result.vFree, xe);
    result.vError = fabs((result.vFree + vion()) / V - 1);
    result.xe = xe;
}

void SahaSolver::calcCore2(double T, double V, calcCoreResult &result, double vEps)
{

    double xe = result.xe;
    double vFree = result.vFree;

    double dxe;
    double vError, vErrorOld;
    double vFreeOld, xeOld;

    vError = 1;

    for(int i = 0; i < 10; i++)
    {
        vFreeOld = vFree;
        vErrorOld = vError;

        vFree = findroot(log(V) - 30, log(V), [&](double vfree) {return vFun(xe, T, V, vfree);}, 1e-14, T, V);

        formX(T, V, vFree, xe);
        dxe = ffV(xe,T,V,vFree);

        vError = fabs((vFree + vion()) / V - 1);

        if(vError > 1e-1) break;

        xeOld = xe;
        xe = xe + dxe;

        if((log(fabs(xe-dxe)) - log(fabs(xe)) < 1e-7))
        {
            if(vError > vErrorOld)
            {
                vFree = vFreeOld;
                vError = vErrorOld;
                xe = xeOld;
                break;
            }

            if(vError < vEps) break;
        }

    }

    result.xe = xe;
    result.vFree = vFree;
    result.vError = vError;

}

double SahaSolver::ff(double xe, double T, double V)
{
    double vFree = Vfree(V, xe);

    if (vFree < 0)
    {
        return std::numeric_limits<double>::max();
    }

    return ffV(xe, T, V, vFree);
}

double SahaSolver::ffV(double xe, double T, double V, double vFree)
{
    double maxH0;
    formH0(mu(T, vFree, xe), p(T, vFree, xe), T, maxH0);

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

double SahaSolver::vFun(double xe, double T, double V, double vFree)
{
    formX(T, V, vFree, xe);
    return (vFree+vion()) / V - 1;
}

SahaPoint SahaSolver::Calculate_TVae(double T, double V)
{
    const double thresholdEps = 1e-4;
    double xe, vFree, vError;
    calcCoreResult res1, res2;

    calcCore1(T, V, res1);
    xe = res1.xe;
    vFree = res1.vFree;
    vError = res1.vError;

    if(res1.vError > thresholdEps)
    {
        res2 = res1;
        calcCore2(T, V, res2, thresholdEps);
        if(res2.vError < res1.vError)
        {
            xe = res2.xe;
            vFree = res2.vFree;
            vError = res2.vError;
        }
    }

    SahaPoint result;
    double E = e(T,vFree, xe);
    double S = s(T,vFree, xe);
    double P = p(T,vFree, xe);

    result.Z = _element.Z;
    result.T = T;
    result.V = V;
    result.E = E;
    result.P = P;
    result.S = S;
    result.Xe = xe;
    result.M = mu(T, vFree, xe);
    result.F = E - T * S;
    result.K = 1 - vFree / V;
    result.vError = vError;

    return result;
}

SahaPoint SahaSolver::Calculate_lgTeV_lgVae(double lgT, double lgV)
{
    return Calculate_TVae(pow(10.0,lgT) / eFi, pow(10.0,lgV)); // что ткое eFi?
    //eFi - это столько эВ в одной атомной единице температуры
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
    //Там в роли DFk берем P * (_element.v[i] - _element.v[i-1]), "круглого" DFk нет совсем // что такое DFk? не нахожу такого обозначения
    //DFk - это "дельта фи ik", индекс i можно выкинуть, он для смеси элементов
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

void SahaSolver::formX(double T, double V, double vFree, double xe)
{
    double maxH0;
    formH0(mu(T, vFree, xe), p(T,vFree, xe), T, maxH0);

    double expSum = 0;
    for(unsigned int i = 0; i <= _element.Z; i++)
    {
        expSum += exp(_H0[i] - maxH0);
    }

    double logX0 = -log(expSum) - maxH0; // что такое logX0 и почему оно считается так?
    //logX0 - это по смыслу log(_x[0]), см. формулу (26)
    //Манипуляции с maxH0 нужны для борьбы с переполнением, см. formH0

    for(unsigned int i = 0; i <= _element.Z; i++)
    {
        _x[i] = exp(logX0 + _H0[i]); // что такое _x[i]?
        // это массив концентраций ионов различной кратности
    }
}

double SahaSolver::Vfree(double V, double xe)
{
    unsigned int i = floor(xe); // _xe - это концентрация здесь? что показывает i?
    // xe - электронная концентрация
    // i - целая часть _xe, что и так, наверное, понятно
    // Содержимое этой функции нигде не описано. Суть происходящего - приближенное вычисление
    // среднего объема ионов, а затем свободного объема. Собственно из-за этой вот "приближенности" весь этот
    // расчет не вполне точен и надо использовать метод Ньютона, чем Вы, собственно, и займетесь.
    double fracXe = xe - i;
    if(i < _element.Z)
    {
        return V - (_element.v[i] * (1-fracXe) + _element.v[i+1] * fracXe);
    }
    return V - _element.v[_element.Z];
}

double SahaSolver::vion()
{
    double V = 0;
    for (int i = 0; i < _element.Z; i++)
    {
        V += _x[i] * _element.v[i];
    }
    return V;
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

void SahaSolver::SahaLeft(std::vector<double> &result)
{
    //Здесь считаем левую часть системы уравнений Саха и помещаем результат в result
    result.resize(_element.Z, -1);//Заглушка, ее надо будет убрать
}

void SahaSolver::vgraph(double lgT, double lgV, double xe)
{
    double T = pow(10.0, lgT) / eFi;
    double V = pow(10.0, lgV);

    FILE *gf = fopen("graph.m","wt+");
    fprintf(gf,"r=[");
    for(double mm = -6; mm < 0.0005; mm += 0.001)
    {
        double vfree = V * pow(10.0, mm);
        formX(T, V, vfree, xe);
        double vi = vion();

        fprintf(gf,"%g %g\n",log10(vfree+vi),ffV(xe,T,V,vfree) + xe);
    }
    fprintf(gf,"];plot(r(:,1),r(:,2),'r+-');\n");
    fclose(gf);

    double vf = findroot(log(V) - 30, log(V), [&](double vfree) {return vFun(xe, T, V, vfree);}, 1e-14, T, V);
    formX(T, V, vf, xe);
    xe = ffV(xe,T,V,vf) + xe;
    printf("vError = %g xe = %g\n",vf + vion() - V, xe);
}

double SahaSolver::p(double T, double vFree, double xe)
{
    return 2*sqrt(2.0)/(3*M_PI*M_PI) * pow(T,2.5) * I15mu_d_t(T,vFree,xe) + T / vFree;
}

double SahaSolver::s(double T, double vFree, double xe)
{
    double Se = 0;
    if(xe > 0)
    {
        Se = sqrt(2.0) / (M_PI*M_PI) * pow(T, 1.5) * vFree * (5.0 / 3.0 * I15mu_d_t(T, vFree, xe) - mu(T,vFree,xe) / T * I05mu_d_t(T, vFree, xe));
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

double SahaSolver::e(double T, double vFree, double xe)
{
    double Ee = sqrt(2.0)/(M_PI*M_PI) * pow(T,2.5) * vFree * I15mu_d_t(T,vFree,xe);

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
