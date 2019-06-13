#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <vector>

class TElement
{
public:
    TElement(unsigned int z, double rCoeff = 1, bool correctV0 = false);
    unsigned int Z;//Заряд элемента
    double A;//Атомный вес в атомных единицах массы
    double ro;//Кристаллографичесая плотность
    std::vector<double> fi;//Массив потенциалов ионизации в а.е.
    std::vector<double> cumFi;//Кумулятивная сумма потенциалов ионизации в а.е.
    std::vector<double> v;//Массив объемов ионных остовов в а.е.
};

namespace elements
{
    const TElement H(1);
    const TElement Fe(26);
    const TElement Cu(29);
    double GetA(unsigned int z);
    double GetRo(unsigned int z);
}

#endif
