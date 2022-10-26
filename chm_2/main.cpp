#include "datatypes.h"
#include "matrix_IO.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <format>

using namespace std;

vector<vector<real_t>> matrixA;  // Вектор диагоналей матрицы А (от верхней к нижней)
vector<real_t> vectorX;
vector<real_t> vectorB;
vector<real_t> nextVectorX;      // Вектор Х для следующих операций
vector<int64_t> diagsShift;

size_t n;                        // Размер матрицы
int64_t m;                       // Число нулевых диагоналей
size_t maxIterations;            // Максимальное число итераций
real_t maxDif;                   // Максимальная невязка
real_t w;                        // Коэф. релаксации

void ReadData()
{
   ifstream matrixParams("./iofiles/matrixParams.txt");
   matrixParams >> n >> m;
   matrixParams.close();

   ifstream solversParams("./iofiles/solversParams.txt");
   solversParams >> maxIterations >> maxDif >> w;
   solversParams.close();

   matrixA.resize(9);
   diagsShift = { 4 + m, 3 + m, 2 + m, 1, 0, -1, -2 - m, -3 - m, -4 - m };
   vectorB.resize(n);
   vectorX.resize(n);
   nextVectorX.resize(n);

   ifstream matrixAFile("./iofiles/matrixA.txt");

   // Данные читаются с верхней диагонали к нижней
   // причём нижние диагонали должны быть прижатыми в правый край,
   // а верхние диагонали в левый край. Пустые элементы заполнить нулями
   for (int i = 0; i < 9; i++)
   {
      matrixA[i] = GetVectorFromFile<real_t>(matrixAFile, n);
   }
   matrixAFile.close();

   vectorB = GetVectorFromFile<real_t>("./iofiles/vectorB.txt", n);

   vectorX = GetVectorFromFile<real_t>("./iofiles/initialX.txt", n);

   for (auto& elem : nextVectorX)
   {
      elem = 0;
   }
}

real_t Norm(vector<real_t> X)
{
   real_t norma = 0;
   for (size_t i = 0; i < n; i++)
   {
      norma += X[i] * X[i];
   }
   norma = sqrt(norma);
   return norma;
}

size_t Iterations(vector<real_t>& vectorXPrev, vector<real_t>& vectorXNext, real_t& resultDif, bool debugOutput = false)
{
   size_t k;
   accum_t sum = 0;
   real_t norm = Norm(vectorB);     // Норма
   accum_t dif = norm;                 // Невязка, делаем её по умолчанию равной norm, чтобы зашло в цикл

   // Итерационный цикл
   for (k = 1; (k <= maxIterations) && (dif > maxDif); k++)
   {
      dif = 0;
      // Пробегаемся по всем строкам матрицы А
      for (size_t i = 0; i < n; i++)
      {
         sum = 0;
         // Пробегаемся по всем столбцам матрицы А
         for (size_t j = 0; j < 9; j++)
         {
            int64_t indX = i + diagsShift[j];
            if (indX >= 0 && indX < n)
            {
               sum += vectorXPrev[indX] * matrixA[j][i];
            }
         }
         vectorXNext[i] = vectorXPrev[i] + (vectorB[i] - sum) * w / matrixA[4][i];
         dif += (vectorB[i] - sum) * (vectorB[i] - sum);       // В данный момент sum = скалярное произведение i-ой строки А на вектор Х
      }

      dif = sqrt(dif) / norm;                               // относительная невязка

      // Выводим на то же место, что и раньше (со сдвигом каретки)
      if (debugOutput)
      {
         cout << format("\rИтерация: {0:<10} относительная невязка: {1:<15.3e}", k, dif);
      }
      if (isinf(dif))
      {
         break;
      }
      swap(vectorXNext, vectorXPrev);
   }
   cout << endl;
   if (debugOutput)
   {
      if (isinf(dif))
      {
         cout << "Выход по переполнению метода" << endl << endl;
      }
      else if (k > maxIterations)
      {
         cout << "Выход по числу итераций" << endl << endl;
      }
      else
      {
         cout << "Выход по относительной невязке" << endl << endl;
      }
   }

   resultDif = dif;
   return k - 1;
}

void RelaxationTester() {
   double w1, w2, step;
   cout << "Тестирование на решение матрицы с разными коэф. релаксации" << endl << endl;
   cout << "Для тестирования метода Якоби диапазон значений должен быть таким:   0 < w1 <= w2 <= 1" << endl;
   cout << "Для тестирования метода Зейделя диапазон значений должен быть таким: 0 < w1 <= w2 < 2" << endl << endl;

   cout << "Введите диапазон коэф. релаксации:" << endl;
   cout << "  начало диапазона w1: ";
   cin >> w1;
   cout << "  конец диапазона w2:  ";
   cin >> w2;
   cout << "  шаг: ";
   cin >> step;

   if (w1 > w2 || w1 < 0 || w2 < 0 || w1 > 2 || w2 > 2)
   {
      cout << "Неправильно введены данные." << endl;
      return;
   }

   int methodNum = 0;
   cout << endl << "Выберите метод для исследования: " << endl;
   cout << "  1) Якоби" << endl;
   cout << "  2) Зейдель" << endl;

   cout << " > ";
   cin >> methodNum;

   auto vec1 = vectorX;
   auto vec2 = nextVectorX;
   string methodName;
   switch (methodNum)
   {
      case 1:
         methodName = "Якоби";
         break;
      case 2: 
         vec2 = vec1;
         methodName = "Зейделя";
         break;
      default:
         cout << "Неверно введено значение." << endl;
         return;
         break;
   }

   cout << endl << "Начало исследования для метода " << methodName << endl << endl;

   size_t minIter = maxIterations;
   real_t minIterRelax = INFINITY;
   real_t minDif = INFINITY;

   vector<real_t> bestVec;
   auto defaultVec1 = vector<real_t>(vec1);
   auto defaultVec2 = vector<real_t>(vec2);

   for (size_t i = 0; i <= (w2 - w1) / step; i++)
   {
      cout << endl << "****************" << endl << endl;

      w = w2 - i * step;
      if (w2 - i * step < w1) w = w1;

      cout << "Текущее начение коэф.релаксации: " << w << endl;

      size_t currIter = 0;
      real_t dif;
      currIter = Iterations(vec1, vec2, dif);
      if (isinf(dif) || isnan(dif))
      {
         cout << "Решение не было получено, метод разошёлся." << endl;
      }
      else
      {
         if (currIter < minIter)
         {
            minIter = currIter;
            minIterRelax = w;
            minDif = dif;
            bestVec = vector<real_t>(vec1);
         }
         if (currIter == minIter && dif < minDif)
         {
            minIterRelax = w;
            minDif = dif;
            bestVec = vector<real_t>(vec1);
         }
         cout << format("Число итераций: {0}\nПолученная относительная невязка: {1:<15.3e}\n", currIter, dif);
         cout << "Полученный вектор решения: " << endl;
         PrintArray(vec1.data(), vec1.size(), g_coutPrecision);
      }

      vec1 = vector<real_t>(defaultVec1);
      vec2 = vector<real_t>(defaultVec2);
   }

   cout << endl << "****************" << endl << endl;

   cout << "Исследование завершено. Получены следующие результаты: " << endl;
   cout << "Наименьшее значение итераций метода на этом диапазоне: " << minIter << endl;
   cout << "Значение релаксации при этом:   " << minIterRelax << endl;
   cout << "Значение относительной невязки: " << minDif << endl;
   cout << "Наиболее точный полученный вектор решения: " << endl;
   PrintArray(bestVec.data(), bestVec.size(), g_coutPrecision);
}

void main()
{
   setlocale(LC_ALL, "ru-RU");

   ReadData();

   int operationNum;
   real_t dif;
   cout << "Выберите метод решения СЛАУ:" << endl;
   cout << "  1) Якоби" << endl;
   cout << "  2) Зейдель" << endl;
   cout << "  4) Исследование на зависимость скорости сходимости от релаксации" << endl;
   cout << "> ";
   cin >> operationNum;
   switch (operationNum)
   {
      case 1:
         Iterations(vectorX, nextVectorX, dif, true);

         cout << "Полученный вектор решения системы: " << endl;
         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, cout);

         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, ofstream(g_outputFileName));
         break;
      case 2:
         Iterations(vectorX, vectorX, dif, true);

         cout << "Полученный вектор решения системы: " << endl;
         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, cout);

         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, ofstream(g_outputFileName));
         break;

      case 4:
         RelaxationTester();
         break;
    }

}
