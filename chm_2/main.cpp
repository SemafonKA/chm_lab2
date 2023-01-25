#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <format>

#include "matrix_IO.h"
#include "Chrono_Timer.h"

using namespace std;

vector<double> matrixA[9];  // Вектор диагоналей матрицы А (от верхней к нижней)
vector<double> vectorX;
vector<double> vectorB;
vector<double> nextVectorX;      // Вектор Х для следующих операций
vector<int64_t> diagsShift;

size_t n;                        // Размер матрицы
int64_t m;                       // Число нулевых диагоналей
size_t maxIterations;            // Максимальное число итераций
double maxDif;                   // Максимальная невязка
double w;                        // Коэф. релаксации

void ReadData()
{
   ifstream matrixParams("./iofiles/matrixParams.txt");
   matrixParams >> n >> m;
   matrixParams.close();

   ifstream solversParams("./iofiles/solversParams.txt");
   solversParams >> maxIterations >> maxDif >> w;
   solversParams.close();

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
      matrixA[i] = GetVectorFromFile<double>(matrixAFile, n);
   }
   matrixAFile.close();

   vectorB = GetVectorFromFile<double>("./iofiles/vectorB.txt", n);

   vectorX = GetVectorFromFile<double>("./iofiles/initialX.txt", n);

   for (auto& elem : nextVectorX)
   {
      elem = 0;
   }
}

double Norm(vector<double> X)
{
   double norma = 0;
   for (size_t i = 0; i < n; i++)
   {
      norma += X[i] * X[i];
   }
   norma = sqrt(norma);
   return norma;
}

size_t Iterations(vector<double>& vectorXPrev, vector<double>& vectorXNext, double& resultDif, bool debugOutput = false)
{
   size_t k;
   double sum = 0;
   double norm = Norm(vectorB);     // Норма
   double dif = norm;                  // Невязка, делаем её по умолчанию равной norm, чтобы зашло в цикл

   // Итерационный цикл
   for (k = 1; (k <= maxIterations) && (dif > maxDif); k++)
   {
      dif = 0;
      // Пробегаемся по всем строкам матрицы А
      for (size_t i = 0; i < n; i++)
      {
         sum = 0;
         int64_t indX;
         // Пробегаемся по всем столбцам матрицы А
         for (size_t j = 0; j < 4; j++)
         {
            indX = i + diagsShift[j];
            if (indX < n)
            {
               sum += vectorXPrev[indX] * matrixA[j][i];
            }
         }
         sum += vectorXPrev[i] * matrixA[4][i];
         for (size_t j = 5; j < 9; j++)
         {
            indX = i + diagsShift[j];
            if (indX >= 0)
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
   if (debugOutput)
   {
      cout << endl;
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

   auto& vec1 = vectorX;
   auto& vec2 = methodNum == 1 ? nextVectorX : vectorX;
   string methodName;
   switch (methodNum)
   {
      case 1:
         methodName = "Якоби";
         break;
      case 2:
         methodName = "Зейделя";
         break;
      default:
         cout << "Неверно введено значение." << endl;
         return;
         break;
   }

   cout << endl << "Начало исследования для метода " << methodName << endl << endl;

   size_t minIter = maxIterations;
   double minIterRelax = INFINITY;
   double minDif = INFINITY;

   vector<double> bestVec;
   auto defaultVec1 = vector<double>(vec1);
   auto defaultVec2 = vector<double>(vec2);

   for (size_t i = 0; i <= (w2 - w1) / step; i++)
   {
      cout << endl << "****************" << endl << endl;

      w = w2 - i * step;
      if (w2 - i * step < w1) w = w1;

      cout << "Текущее начение коэф.релаксации: " << w << endl;

      size_t currIter = 0;
      double dif;
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
            bestVec = vector<double>(vec1);
         }
         if (currIter == minIter && dif < minDif)
         {
            minIterRelax = w;
            minDif = dif;
            bestVec = vector<double>(vec1);
         }
         cout << format("Число итераций: {0}\nПолученная относительная невязка: {1:<15.3e}\n", currIter, dif);
         cout << "Полученный вектор решения: " << endl;
         PrintArray(vec1.data(), vec1.size(), 15);
      }

      vec1 = vector<double>(defaultVec1);
      vec2 = vector<double>(defaultVec2);
   }

   cout << endl << "****************" << endl << endl;

   cout << "Исследование завершено. Получены следующие результаты: " << endl;
   cout << "Наименьшее значение итераций метода на этом диапазоне: " << minIter << endl;
   cout << "Значение релаксации при этом:   " << minIterRelax << endl;
   cout << "Значение относительной невязки: " << minDif << endl;
   cout << "Наиболее точный полученный вектор решения: " << endl;
   PrintArray(bestVec.data(), bestVec.size(), 15);
}

void main()
{
   setlocale(LC_ALL, "ru-RU");

   ReadData();

   int operationNum;
   double dif;
   cout << "Выберите метод решения СЛАУ:" << endl;
   cout << "  1) Якоби" << endl;
   cout << "  2) Зейдель" << endl;
   cout << "  3) Исследование на зависимость скорости сходимости от релаксации" << endl;
   cout << "> ";
   cin >> operationNum;
   switch (operationNum)
   {
      case 1:
      {
         Timer timer;
         size_t it = Iterations(vectorX, nextVectorX, dif, false);
         timer.elapsed();
         cout << "Время решения: " << timer.elapsedValue * 100 << " мс" << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << dif << endl;

         cout << "Полученный вектор решения системы: " << endl;
         PrintArray<double>(vectorX.data(), vectorX.size(), 15, cout);

         PrintArray<double>(vectorX.data(), vectorX.size(), 15, ofstream("./iofiles/Output.txt"));
         break;
      }
      case 2:
      {
         Timer timer;
         size_t it = Iterations(vectorX, vectorX, dif, false);
         timer.elapsed();
         cout << "Время решения: " << timer.elapsedValue * 100 << " мс" << endl;
         cout << "Количество итераций: " << it << endl;
         cout << "Относительная невязка: " << dif << endl;

         cout << "Полученный вектор решения системы: " << endl;
         PrintArray<double>(vectorX.data(), vectorX.size(), 15, cout);

         PrintArray<double>(vectorX.data(), vectorX.size(), 15, ofstream("./iofiles/Output.txt"));
         break;
      }
      case 3:
         RelaxationTester();
         break;
   }

}
