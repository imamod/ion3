#ifndef SAHASOLVER
#define SAHASOLVER

#include "elements.h"
#include <vector>

struct SahaPoint
{
   unsigned int Z;  // атомный номер
   double lgT;		  // температура, а. е.
   double lgV;		  // объём атомной ячейки, а. е.
   double lgP;		  // давление, а. е.
   double lgE;		  // энергия, а. е.
   double lgS;		  // энтропия, а. е.
   double M;		  // химический потенциал, а. е.
   double F;        // свободная энергия
   double Xe;       // Ионизация

   double lgKappa;  // объёмная доля элекронных остовов в электронном газе, б/р
};

class SahaSolver
{
public:
    SahaSolver(const TElement &element);
    SahaPoint Calculate_TVae(double T, double V);
    SahaPoint Calculate_lgTeV_lgVae(double lgT, double lgV);
    void GetX(std::vector<double> &x);

private:
    void formH0(double mu, double P, double T, double &maxH0);
    double ff(double xe, double T, double V);
    double Vfree(double V);
    double e(double T, double vFree);
    double p(double T, double vFree);
    double s(double T, double vFree);
    void formX(double T, double V);
    const TElement &_element;
    std::vector<double> _x;
    std::vector<double> _H0;
    double _xe;
};

#endif // SAHASOLVER
