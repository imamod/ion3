#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <vector>

struct TElement
{
    TElement(unsigned int z, double rCoeff = 1);
    unsigned int Z;//Заряд элемента
	double A;//Атомный вес в атомных единицах массы
    std::vector<double> fi;//Массив потенциалов ионизации в а.е.
    std::vector<double> cumFi;//Кумулятивная сумма потенциалов ионизации в а.е.
    std::vector<double> v;//Массив объемов ионных остовов в а.е.
};

namespace elements
{
    const TElement H(1);
    const TElement Fe(26);
    const TElement Cu(29);
}

#endif
