#ifndef SAHASOLVER
#define SAHASOLVER

#include "elements.h"
#include <vector>
#include <string>
#include <functional>

struct SahaPoint
{
   unsigned int Z;  // атомный номер
   double T;		// температура, а. е.
   double V;		// объём атомной ячейки, а. е.
   double P;		// давление, а. е.
   double E;		// энергия, а. е.
   double S;		// энтропия, а. е.
   double M;		// химический потенциал, а. е.
   double F;        // свободная энергия
   double Xe;       // Ионизация

   double K;  // объёмная доля элекронных остовов в электронном газе, б/р
   double vError;
};

class SahaSolver
{
public:
    SahaSolver(const TElement &element);
    SahaPoint Calculate_TVae(double T, double V);
    SahaPoint Calculate_lgTeV_lgVae(double lgT, double lgV);
    void GetX(std::vector<double> &x);
	double Vion(double rCoeff);
    void SahaLeft(std::vector<double> &result);
    void vgraph(double lgT, double lgV, double xe);

private:

    struct calcCoreResult
    {
        double xe;
        double vFree;
        double vError;
    };

    double findroot(double logA, double logB, const std::function<double(double)> &F, double eps, double T, double V);
    void calcCore1(double T, double V, calcCoreResult &result);
    void calcCore2(double T, double V, calcCoreResult &result, double vEps);
	void error(const std::string & errorType, const std::string & message, double T, double V);
    void formH0(double mu, double P, double T, double &maxH0);
    double ff(double xe, double T, double V);
    double ffV(double xe, double T, double V, double vFree);
    double vFun(double xe, double T, double V, double vFree);
    double Vfree(double V, double xe);
    double vion();

    double e(double T, double vFree, double xe);
    double p(double T, double vFree, double xe);
    double s(double T, double vFree, double xe);
    void formX(double T, double V, double vFree, double xe);
    const TElement &_element;
    std::vector<double> _x;
    std::vector<double> _H0;


};

#endif // SAHASOLVER
