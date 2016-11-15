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

      double lgKappa;  // объёмная доля элекронных остовов в электронном газе, б/р
      double dLgKdLgV;
      double dLgKdLgT;
   };

   /// Функция расчета терммодинамических величин по модели Саха.
   /// Z - атомный номер элемента
   /// lgT - температура, а. е.
   /// lgV - объём электронной ячейки, а. е.
   ///
   /// В случае ошибки выбрасывается исключение std::exception.
   Point Calculate(unsigned int i_Z, double i_lgT, double i_lgV);

} // namespace

#endif
