#include "saha.h"
#include "src/sahasolver.h"
#include "src/atom_ed.h"
#include <cmath>
#include <sstream>
#include <memory>

namespace saha
{
    static double externalAtomicWeight = 0;
    static double externalRho = 0;
    static double teta = 0;
    static bool correctV0 = false;
    void SetTeta(double value) {teta = value;}
    void SetCorrectV0(double value) {correctV0 = value;}

	namespace
	{
		const unsigned int c_maxZ = 103;
        double _ionRadiusCoeff = 0; //Нулевое значение - режим совместимости со старой Сахой
	
        std::shared_ptr<TElement> elem, elem0;
        std::shared_ptr<SahaSolver> solver, solver0;
	   
        SahaPoint calculate(unsigned int i_Z, double i_lgT, double i_lgV, double ionRadiusCoeff)
		{
			if ((i_Z > c_maxZ) || (i_Z < 1))
			{
				std::ostringstream oss;
				oss << "Invalid element number " << i_Z << ". Should be <= " << c_maxZ << " and >= 1." << std::endl;
				throw std::invalid_argument(oss.str());
			}

            //printf("[%g]",teta);

			if (!elem)
			{
                _ionRadiusCoeff = ionRadiusCoeff;
                elem.reset(new TElement(i_Z, _ionRadiusCoeff, correctV0));
                solver.reset(new SahaSolver(*elem, teta));
                //elem0.reset(new TElement(i_Z, 0, false));
                //solver0.reset(new SahaSolver(*elem0, teta));

                /*printf("volumes:");
                for(auto v : elem->v) printf("%g ",v);
                printf("\nro = %g\n",elem->ro);*/
			}
			else
			{
                if ((elem->Z != i_Z) || (_ionRadiusCoeff != ionRadiusCoeff))
				{
                    _ionRadiusCoeff = ionRadiusCoeff;
                    elem.reset(new TElement(i_Z, _ionRadiusCoeff, correctV0));
                    solver.reset(new SahaSolver(*elem, teta));
                    //elem0.reset(new TElement(i_Z, 0, false));
                    //solver0.reset(new SahaSolver(*elem0, teta));
                }
			}

			SahaPoint result = solver->Calculate_TVae(pow(10.0, i_lgT), pow(10.0, i_lgV));
            //SahaPoint result0 = solver0->Calculate_TVae(pow(10.0, i_lgT), pow(10.0, i_lgV));

            //double V = pow(10.0, i_lgV);
            //result.vFactor = exp(-solver0->Vion(_ionRadiusCoeff) / V);
            //result.Xe = result0.Xe;
            //printf("<%g>",result.vFactor);

			//В режиме совместимости с обычной Сахой считаем K по-старому
            //if (_ionRadiusCoeff == 0) result.K = solver->Vion(1.0) / pow(10.0, i_lgV);

            //Склеивающий параметр - химический потенциал
            //result.K = pow(result.IMu,1.5);

			return result;
      }
   }

    void SetExternalAtomicWeight(double A)
    {
        externalAtomicWeight = A;
    }

    void SetExternalRho(double rho)
    {
        externalRho = rho;
    }

    double GetA(unsigned int i_Z)
    {
        return (externalAtomicWeight > 0) ? externalAtomicWeight : elements::GetA(i_Z);
    }

    Point Calculate(unsigned int i_Z, double i_lgT, double i_lgV, double ionRadiusCoeff)
	{
      const SahaPoint sp = calculate(i_Z, i_lgT, i_lgV, ionRadiusCoeff);
      const double dArg = 0.01;

      double dLgKdLgV = 0;
      double dLgKdLgT = 0;
      double dLgPdLgT = 0;

      const SahaPoint sp_v_left = calculate(i_Z, i_lgT, i_lgV - dArg, ionRadiusCoeff);
      const SahaPoint sp_v_right = calculate(i_Z, i_lgT, i_lgV + dArg, ionRadiusCoeff);

      //Костыли для борьбы с переполнением...
      const double k_v_left = log10(std::max(sp_v_left.K, 1e-308));
      const double k_v_right = log10(std::max(sp_v_right.K, 1e-308));

      const SahaPoint sp_t_left = calculate(i_Z, i_lgT - dArg, i_lgV, ionRadiusCoeff);
      const SahaPoint sp_t_right = calculate(i_Z, i_lgT + dArg, i_lgV, ionRadiusCoeff);

      const double k_t_left = log10(sp_t_left.K);
      const double k_t_right = log10(sp_t_right.K);

      dLgPdLgT = (log10(sp_t_right.P) - log10(sp_t_left.P)) / 2.0 / dArg;

      {
        dLgKdLgV = (k_v_right - k_v_left) / 2.0 / dArg;
        dLgKdLgT = (k_t_right - k_t_left) / 2.0 / dArg;
      }

      const double V = pow(10.0, i_lgV);
      const double T = pow(10.0, i_lgT);

      Point pt;
      pt.Z = i_Z;
      pt.T = T;
      pt.V = V;
      pt.P = sp.P;
	  pt.DPQuip = 0;
      pt.E = sp.E;
      pt.S = sp.S;
      pt.M = sp.M;
      pt.lgKappa = log10(sp.K);
      pt.lgIMu = log10(sp.IMu);
      pt.F = sp.E - sp.S * T;
      pt.dLgKdLgV = dLgKdLgV;
      pt.dLgKdLgT = dLgKdLgT;
      pt.dLgPdLgT = dLgPdLgT;
      pt.Xe = sp.Xe;
      pt.vFactor = sp.vFactor;
      pt.x = sp.x;
      pt.zd = sp.zd;
      return pt;
    }

    double GetRo(unsigned int i_Z)
    {
        return (externalRho > 0) ? externalRho : elements::GetRo(i_Z);
    }

    void TestIonVolumes()
    {
        std::shared_ptr<TElement> elem;
        for(int i = 1; i <= 103; i++)
        {
            elem.reset(new TElement(i, 0.6, true));
            for(int j = 1; j < i;j++)
            {
                if(elem->v[j] > elem->v[j-1])
                {
                    printf("element %d Zion = %d: %g %g\n", i, j-1, elem->v[j-1], elem->v[j]);
                }
            }
        }
    }

    void FullSahaTable::Clear()
    {
        _data.clear();
    }

    void FullSahaTable::AddCalculatedPoint(unsigned int Z, double lgT_eV, double lgRho, double ionRadiusCoeff)
    {
        const double lgRoConst = log10(saha::GetA(Z)) + log10(eRo);
        const double log_eV = log10(eFi);

        Point point = Calculate(Z, lgT_eV - log_eV, lgRoConst - lgRho, ionRadiusCoeff);
        _data.push_back(point);
    }

    void FullSahaTable::CalculateTable(unsigned int Z, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep)
    {
        _data.clear();
        for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2; lgT -= lgTStep)
        {
            printf("[%g]",lgT);fflush(0);
            for(double lgRho = lgRhoMin; lgRho < lgRhoMax + lgRhoStep / 2; lgRho += lgRhoStep)
            {
                AddCalculatedPoint(Z, lgT, lgRho, rCoeff);
            }
        }
    }

    void FullSahaTable::PrintIonVector(const std::string &fileName)
    {
        FILE *f = fopen(fileName.c_str(), "wt+");

        const double log_eV = log10(eFi);
        for(Point &elem : _data)
        {
            double zh = 0, zd = 0;
            for(int i = 1; i <= elem.Z; i++)
            {
                zh += elem.x[i] * pow(i, 1.5);
                zd += elem.x[i] * i * i;
            }

            zh = pow(zh, 2/3.0);
            zd = pow(zd, 0.5);

            fprintf(f, "%6.3f %6.3f %20.15f %20.15f %20.15f ", log10(elem.T) + log_eV, log10(saha::GetA(elem.Z)) + log10(eRo) - log10(elem.V), elem.Xe, zd, zh);
            for(double xElem : elem.x)
            {
                fprintf(f, "%20.15f ", xElem);
            }
            fprintf(f, "\n");
        }

        fclose(f);
    }

} // namespace saha
