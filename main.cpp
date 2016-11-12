#include "src/elements.h"
#include "src/sahasolver.h"
#include "src/atom_ed.h"
#include <math.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void outputArray(std::ostream& os, const std::string dataName, const std::vector<double> &data)
{
   os << dataName << " = [";
   for(std::vector<double>::const_iterator it = data.begin(); it!= data.end(); ++it) {os << *it << " ";};
   os << "];" << std::endl;
}

void outputTable(std::ostream& os, std::string tableName, const std::vector<std::vector<double>> &table)
{
    os << tableName << " = [" << std::endl;
    for (size_t i = 0; i < table.size(); ++i)
    {
        const std::vector<double>& line = table[i];
        for (size_t j = 0; j < line.size(); ++j)
        {
            os << (line[j]) << " ";
        }
        os << std::endl;
    }
    os << "];" << std::endl;
}

void CrashTest(double rCoeff, double lgVMin, double lgVMax, double lgVStep, double lgTMin, double lgTMax, double lgTStep)
{
	for (int Z = 1; Z <= 103; Z++)
	{
		const TElement elem(Z, rCoeff);
		SahaSolver solver(elem);

		printf("[%d] ",Z);

		for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep)
		{
			for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
			{
				SahaPoint res = solver.Calculate_lgTeV_lgVae(lgT, lgV);
			}
		}

	}
}

void calculator(unsigned int Z, double rCoeff, double lgVMin, double lgVMax, double lgVStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    const TElement elem(Z, rCoeff);
    SahaSolver solver(elem);

    std::vector<double> lgTPhys;
    std::vector<double> lgVa;
    std::vector<double> lgRho;

    std::vector<std::vector<double>> ionizationTable;
    std::vector<std::vector<double>> pTable;
    std::vector<std::vector<double>> eTable;

    for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep) lgTPhys.push_back(lgT);

    double roConst = log10(eRo*elem.A);
    for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
    {
        lgVa.push_back(lgV);
        lgRho.push_back(roConst - lgV);
    }

    for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep)
    {
        std::cout << "[" << lgT << "]";
        std::vector<double> ionizationLine;
        std::vector<double> pLine;
        std::vector<double> eLine;
        for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
        {
            SahaPoint res = solver.Calculate_lgTeV_lgVae(lgT,lgV);
            ionizationLine.push_back(res.Xe);
            pLine.push_back(pow(10,res.lgP));
            eLine.push_back(pow(10,res.lgE));
        }
        ionizationTable.push_back(ionizationLine);
        pTable.push_back(pLine);
        eTable.push_back(eLine);
    }

    std::fstream f(filename.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=" << Z << ";" << std::endl;
    outputArray(f, "lgT", lgTPhys);
    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", lgRho);
    outputTable(f, "xe_Saha", ionizationTable);
    outputTable(f, "P_Saha", pTable);
    outputTable(f, "E_Saha", eTable);
}

int main()
{
	CrashTest(0.6, -3, 6.01, 0.05, -0.51, 4.6, 0.05);
    //calculator(13,0.6,-3,6.01,0.05,0.49,4.6,0.05,"../../mion2/saha_Al.m");
    return 0;
}
