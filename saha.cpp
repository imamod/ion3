#include "saha.h"
#include "src/sahasolver.h"
#include <cmath>
#include <sstream>
#include <memory>

namespace saha
{

	namespace
	{
		const unsigned int c_maxZ = 103;
		//const double ionRadiusCoeff = 0; //Нулевое значение - режим совместимости со старой Сахой
		const double ionRadiusCoeff = 0.6;

		std::shared_ptr<TElement> elem;
		std::shared_ptr<SahaSolver> solver;

		SahaPoint calculate(unsigned int i_Z, double i_lgT, double i_lgV)
		{
			if ((i_Z > c_maxZ) || (i_Z < 1))
			{
				std::ostringstream oss;
				oss << "Invalid element number " << i_Z << ". Should be <= " << c_maxZ << " and >= 1." << std::endl;
				throw std::invalid_argument(oss.str());
			}

			if (!elem)
			{
				elem.reset(new TElement(i_Z, ionRadiusCoeff));
				solver.reset(new SahaSolver(*elem));
			}
			else
			{
				if (elem->Z != i_Z)
				{
					elem.reset(new TElement(i_Z, ionRadiusCoeff));
					solver.reset(new SahaSolver(*elem));
				}
			}

			SahaPoint result = solver->Calculate_TVae(pow(10.0, i_lgT), pow(10.0, i_lgV));

			//В режиме совместимости с обычной Сахой считаем K по-старому
			if (ionRadiusCoeff == 0) result.K = solver->Vion(1.0) / pow(10.0, i_lgV);

			return result;
      }
   }

	Point Calculate(unsigned int i_Z, double i_lgT, double i_lgV)
	{
      const SahaPoint sp = calculate(i_Z, i_lgT, i_lgV);
      const double dArg = 0.05;

      const double k_v_left = log10(calculate(i_Z, i_lgT, i_lgV - dArg).K);
      const double k_v_right = log10(calculate(i_Z, i_lgT, i_lgV + dArg).K);
      const double k_t_left = log10(calculate(i_Z, i_lgT - dArg, i_lgV).K);
      const double k_t_right = log10(calculate(i_Z, i_lgT + dArg, i_lgV).K);

      const double dLgKdLgV = (k_v_right - k_v_left) / 2.0 / dArg;
      const double dLgKdLgT = (k_t_right - k_t_left) / 2.0 / dArg;

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
      pt.F = sp.E - sp.S * T;
      pt.dLgKdLgV = dLgKdLgV;
      pt.dLgKdLgT = dLgKdLgT;
      pt.Xe = sp.Xe;
      return pt;
   }

} // namespace saha
