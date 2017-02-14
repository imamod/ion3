#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <vector>
// проблемы с кодировкой. Можно ли переписать комментарии?

struct TElement
{
    TElement(unsigned int z, double rCoeff = 1);
    unsigned int Z;//����� ��������
	double A;//������� ��� � ������� �������� �����
    std::vector<double> fi;//������ ����������� ��������� � �.�.
    std::vector<double> cumFi;//������������ ����� ����������� ��������� � �.�.
    std::vector<double> v;//������ ������� ������ ������� � �.�.
};

namespace elements
{
    const TElement H(1);
    const TElement Fe(26);
    const TElement Cu(29);
}

#endif
