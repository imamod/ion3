#include "src/elements.h"
#include "src/sahasolver.h"
#include "src/atom_ed.h"
#include "saha.h"
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
            pLine.push_back(res.P);
            eLine.push_back(res.E);
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

void testSahaLeft()
{
    const TElement elem(29, 0.6); //Расчет для меди c Z=29
    SahaSolver solver(elem);

    SahaPoint res = solver.Calculate_lgTeV_lgVae(3,-0.5);
    printf("xe = %g\n",res.Xe);

    printf("Vector SahaLeft:");
    std::vector<double> result;
    solver.SahaLeft(result);

    for(auto &x : result) printf("%g ",x);
    printf("\n");
}

void vtest()
{
    const TElement elem(29, 0.6); //Расчет для меди c Z=29
    SahaSolver solver(elem);
    solver.vtest(3, -1.7, 19.0398, 0.01, 150);
}

int main()
{
	try
	{
        testSahaLeft();
        vtest();
        //testSahaLeft();
		//CrashTest(0.6, -3, 6.01, 0.05, -1.51, 4.6, 0.05);
        //calculator(82, 0.6, -3, 6.01, 0.05, -5.51, 4.6, 0.05, "saha_Pb.m");

        /*saha::Point ppp;
		ppp = saha::Calculate(26, 1.5, 2);

		printf("\n%d %g %g\n", ppp.Z, ppp.P, ppp.lgKappa);

		ppp = saha::Calculate(26, 1.5, 2.5);
        printf("%d %g %g\n", ppp.Z, ppp.P, ppp.lgKappa);*/
	}
	catch (std::exception& r)
	{
		printf("\n%s", r.what());
	}

    return 0;
}
