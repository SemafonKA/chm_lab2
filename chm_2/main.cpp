#include "datatypes.h"
#include "matrix_IO.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

vector<vector<real_t>> matrixA;
vector<real_t> vectorX, vectorB;
vector<real_t> nextVectorX;   // Вектор Х для следующих операций
vector<int64_t> diagIndexes;

size_t n, m, maxIterations;   // размер матрицы, число нулевых диагоналей, максимальное число итераций
real_t maxDif, w;    // Максимальная невязка, коэф. релаксации

void ReadData()
{
   ifstream matrixParams("./iofiles/matrixParams.txt");
   matrixParams >> n >> m >> w;
   matrixParams.close();

   ifstream solversParams("./iofiles/solversParams.txt");
   solversParams >> maxIterations >> maxDif;
   solversParams.close();

   matrixA.resize(9);
   diagIndexes.resize(9);
   vectorB.resize(n);
   vectorX.resize(n);
   nextVectorX.resize(n);

   diagIndexes[3] = -1;
   diagIndexes[4] = 0;
   diagIndexes[5] = 1;
   for (size_t i = 0; i < 3; i++)
   {
      diagIndexes[6 + i] = 2 + m + i;
      diagIndexes[2 - i] = -diagIndexes[6 + i];
   }

   ifstream matrixAFile("./iofiles/matrixA.txt");

   // Данные читаются с нижней диагонали к верхней
   // причём нижние диагонали должны быть прижатыми в правый край,
   // а верхние диагонали в нижний край. Пустые элементы заполнить нулями
   for (int i = 0; i < 4; i++)
   {
      matrixA[i] = GetVectorFromFile<real_t>(matrixAFile, n);
   }
   matrixA[4] = GetVectorFromFile<real_t>(matrixAFile, n);
   for (int i = 5; i < 9; i++)
   {
      matrixA[i] = GetVectorFromFile<real_t>(matrixAFile, n);
   }
   matrixAFile.close();

   vectorB = GetVectorFromFile<real_t>("./iofiles/vectorB.txt", n);

   vectorX = GetVectorFromFile<real_t>("./iofiles/initialX.txt", n);

   for (auto& elem : nextVectorX) 
      elem = 0;
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

void Iterations(vector<double>& vectorXPrev, vector<double>& vectorXNext)
{
   size_t k;
   accum_t sum = 0;
   real_t norm = Norm(vectorB);     // Норма
   accum_t dif = norm;                 // Невязка, делаем её по умолчанию равной norm, чтобы зашло в цикл
   
   // Итерационный цикл
   for (k = 0; (k < maxIterations) && (dif > maxDif); k++)  
   {
      dif = 0;
      // Пробегаемся по всем строкам матрицы А
      for (int i = 0; i < n; i++)
      {
         // Пробегаемся по всем столбцам матрицы А
         for (size_t j = 0; j < 9; j++)
         {
            // Суммируем с диагональным элементом
            if (j == 4)
            {
               sum = vectorXPrev[i] * matrixA[4][i];
            }
            // Суммируем по нижним диагоналям
            if (j < 4 && diagIndexes[j] + i >= 0)
            {
               sum += vectorXPrev[diagIndexes[j] + i] * matrixA[j][i];
            }
            // Суммируем по верхним диагоналям
            if (j > 4 && diagIndexes[j] + i < n)
            {
               sum += vectorXPrev[diagIndexes[j] + i] * matrixA[j][i];
            }
         }
         vectorXNext[i] = vectorXPrev[i] + (vectorB[i] - sum) * w / matrixA[4][i];
         dif += (vectorB[i] - sum) * (vectorB[i] - sum);       // В данный момент sum = скалярное произведение i-ой строки А на вектор Х
      }

      dif = sqrt(dif) / norm;                               //относительная невязка

      cout << "Итерация: " << k << endl << "\t относительная невязка: " << dif << endl;

      swap(vectorXNext, vectorXPrev);
   }

   if (k == maxIterations)
   {
      cout << "Выход по числу итераций" << endl;
   }
   else
   {
      cout << "Выход по относительной невязке" << endl;
   }
}

void main()
{
   setlocale(LC_ALL, "ru-RU");

   ReadData();

   int operationNum;
   cout << "Выберите метод решения СЛАУ:" << endl << "1 - Якоби" << endl << "2 - Зейдель" << endl;
   cout << "> ";
   cin >> operationNum;
   switch (operationNum)
   {
      case 1:
         Iterations(vectorX, nextVectorX);
         break;
      case 2:
         Iterations(vectorX, vectorX);
         break;
   }

   cout << "Полученный вектор решения системы: " << endl;
   PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, cout);

   PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, ofstream(g_outputFileName));
}