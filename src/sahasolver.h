#ifndef PLASMACOMPONENT
#define PLASMACOMPONENT

#include "elements.h"
#include <vector>


class SahaSolver
{
public:
    SahaSolver(const TElement &element);
    double calculate_TVae(double T, double V);
    double calculate_lgTeV_lgVae(double lgT, double lgV);
    double P(double T, double V);
    double E(double T, double V);
    double xe();

private:
    double ff(double xe, double T, double V);
    double Vfree(double V);
    double p(double T, double vFree);
    void formX();
    const TElement &_element;
    std::vector<double> _x;
    std::vector<double> _temp;
    double _xe;
};

#endif // PLASMACOMPONENT
