#include "datatypes.h"
#include "matrix_IO.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <format>

using namespace std;

vector<vector<real_t>> matrixA;  // ������ ���������� ������� � (�� ������� � ������)
vector<real_t> vectorX;
vector<real_t> vectorB;
vector<real_t> nextVectorX;      // ������ � ��� ��������� ��������
vector<int64_t> diagsShift;

size_t n;                        // ������ �������
int64_t m;                       // ����� ������� ����������
size_t maxIterations;            // ������������ ����� ��������
real_t maxDif;                   // ������������ �������
real_t w;                        // ����. ����������

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

   // ������ �������� � ������� ��������� � ������
   // ������ ������ ��������� ������ ���� ��������� � ������ ����,
   // � ������� ��������� � ����� ����. ������ �������� ��������� ������
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
   real_t norm = Norm(vectorB);     // �����
   accum_t dif = norm;                 // �������, ������ � �� ��������� ������ norm, ����� ����� � ����

   // ������������ ����
   for (k = 1; (k <= maxIterations) && (dif > maxDif); k++)
   {
      dif = 0;
      // ����������� �� ���� ������� ������� �
      for (size_t i = 0; i < n; i++)
      {
         sum = 0;
         // ����������� �� ���� �������� ������� �
         for (size_t j = 0; j < 9; j++)
         {
            int64_t indX = i + diagsShift[j];
            if (indX >= 0 && indX < n)
            {
               sum += vectorXPrev[indX] * matrixA[j][i];
            }
         }
         vectorXNext[i] = vectorXPrev[i] + (vectorB[i] - sum) * w / matrixA[4][i];
         dif += (vectorB[i] - sum) * (vectorB[i] - sum);       // � ������ ������ sum = ��������� ������������ i-�� ������ � �� ������ �
      }

      dif = sqrt(dif) / norm;                               // ������������� �������

      // ������� �� �� �� �����, ��� � ������ (�� ������� �������)
      if (debugOutput)
      {
         cout << format("\r��������: {0:<10} ������������� �������: {1:<15.3e}", k, dif);
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
         cout << "����� �� ������������ ������" << endl << endl;
      }
      else if (k > maxIterations)
      {
         cout << "����� �� ����� ��������" << endl << endl;
      }
      else
      {
         cout << "����� �� ������������� �������" << endl << endl;
      }
   }

   resultDif = dif;
   return k - 1;
}

void RelaxationTester() {
   double w1, w2, step;
   cout << "������������ �� ������� ������� � ������� ����. ����������" << endl << endl;
   cout << "��� ������������ ������ ����� �������� �������� ������ ���� �����:   0 < w1 <= w2 <= 1" << endl;
   cout << "��� ������������ ������ ������� �������� �������� ������ ���� �����: 0 < w1 <= w2 < 2" << endl << endl;

   cout << "������� �������� ����. ����������:" << endl;
   cout << "  ������ ��������� w1: ";
   cin >> w1;
   cout << "  ����� ��������� w2:  ";
   cin >> w2;
   cout << "  ���: ";
   cin >> step;

   if (w1 > w2 || w1 < 0 || w2 < 0 || w1 > 2 || w2 > 2)
   {
      cout << "����������� ������� ������." << endl;
      return;
   }

   int methodNum = 0;
   cout << endl << "�������� ����� ��� ������������: " << endl;
   cout << "  1) �����" << endl;
   cout << "  2) �������" << endl;

   cout << " > ";
   cin >> methodNum;

   auto vec1 = vectorX;
   auto vec2 = nextVectorX;
   string methodName;
   switch (methodNum)
   {
      case 1:
         methodName = "�����";
         break;
      case 2: 
         vec2 = vec1;
         methodName = "�������";
         break;
      default:
         cout << "������� ������� ��������." << endl;
         return;
         break;
   }

   cout << endl << "������ ������������ ��� ������ " << methodName << endl << endl;

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

      cout << "������� ������� ����.����������: " << w << endl;

      size_t currIter = 0;
      real_t dif;
      currIter = Iterations(vec1, vec2, dif);
      if (isinf(dif) || isnan(dif))
      {
         cout << "������� �� ���� ��������, ����� ���������." << endl;
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
         cout << format("����� ��������: {0}\n���������� ������������� �������: {1:<15.3e}\n", currIter, dif);
         cout << "���������� ������ �������: " << endl;
         PrintArray(vec1.data(), vec1.size(), g_coutPrecision);
      }

      vec1 = vector<real_t>(defaultVec1);
      vec2 = vector<real_t>(defaultVec2);
   }

   cout << endl << "****************" << endl << endl;

   cout << "������������ ���������. �������� ��������� ����������: " << endl;
   cout << "���������� �������� �������� ������ �� ���� ���������: " << minIter << endl;
   cout << "�������� ���������� ��� ����:   " << minIterRelax << endl;
   cout << "�������� ������������� �������: " << minDif << endl;
   cout << "�������� ������ ���������� ������ �������: " << endl;
   PrintArray(bestVec.data(), bestVec.size(), g_coutPrecision);
}

void main()
{
   setlocale(LC_ALL, "ru-RU");

   ReadData();

   int operationNum;
   real_t dif;
   cout << "�������� ����� ������� ����:" << endl;
   cout << "  1) �����" << endl;
   cout << "  2) �������" << endl;
   cout << "  4) ������������ �� ����������� �������� ���������� �� ����������" << endl;
   cout << "> ";
   cin >> operationNum;
   switch (operationNum)
   {
      case 1:
         Iterations(vectorX, nextVectorX, dif, true);

         cout << "���������� ������ ������� �������: " << endl;
         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, cout);

         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, ofstream(g_outputFileName));
         break;
      case 2:
         Iterations(vectorX, vectorX, dif, true);

         cout << "���������� ������ ������� �������: " << endl;
         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, cout);

         PrintArray<real_t>(vectorX.data(), vectorX.size(), g_coutPrecision, ofstream(g_outputFileName));
         break;

      case 4:
         RelaxationTester();
         break;
    }

}
