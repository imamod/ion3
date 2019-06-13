#ifndef SAHA
#define SAHA

#include <stdexcept>
#include <vector>

namespace saha
{

   // Точка расчета по модели Саха
   struct Point
   {
      unsigned int Z;  // атомный номер
      double T;      // температура, а. е.
      double V;      // объём атомной ячейки, а. е.
      double P;      // давление, а. е.
      double DPQuip;   // log P quip
      double E;      // энергия, а. е.
      double S;      // энтропия, а. е.
      double M;        // химический потенциал, а. е.
      double F;        // свободная энергия
      double Xe;       // Ионизация
      double zd;

      std::vector<double> x; //Вектор ионизаций

      double vFactor;
      double lgKappa;  // объёмная доля элекронных остовов в электронном газе, б/р
      double lgIMu; // lg(I05(Mu/T))
      double dLgKdLgV;
      double dLgKdLgT;
      double dLgPdLgT;
   };

   class FullSahaTable
   {
   public:
       void Clear();
       void AddCalculatedPoint(unsigned int Z, double lgT_eV, double lgRho, double ionRadiusCoeff);
       void CalculateTable(unsigned int Z, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep);
       void PrintIonVector(const std::string &fileName);

   protected:
       std::vector<Point> _data;
   };

   /// Функция расчета терммодинамических величин по модели Саха.
   /// Z - атомный номер элемента
   /// lgT - температура, а. е.
   /// lgV - объём электронной ячейки, а. е.
   ///
   /// В случае ошибки выбрасывается исключение std::exception.
   Point Calculate(unsigned int i_Z, double i_lgT, double i_lgV, double ionRadiusCoeff = 0);
   void SetExternalAtomicWeight(double A);
   void SetExternalRho(double rho);
   double GetA(unsigned int i_Z);
   double GetRo(unsigned int i_Z);
   void SetTeta(double value);
   void SetCorrectV0(double value);
   void TestIonVolumes();
} // namespace

#endif
