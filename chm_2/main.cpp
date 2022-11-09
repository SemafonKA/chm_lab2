#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <format>

#include "matrix_IO.h"

using namespace std;

vector<double> matrixA[9];  // ������ ���������� ������� � (�� ������� � ������)
vector<double> vectorX;
vector<double> vectorB;
vector<double> nextVectorX;      // ������ � ��� ��������� ��������
vector<int64_t> diagsShift;

size_t n;                        // ������ �������
int64_t m;                       // ����� ������� ����������
size_t maxIterations;            // ������������ ����� ��������
double maxDif;                   // ������������ �������
double w;                        // ����. ����������

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

   // ������ �������� � ������� ��������� � ������
   // ������ ������ ��������� ������ ���� ��������� � ������ ����,
   // � ������� ��������� � ����� ����. ������ �������� ��������� ������
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
   double norm = Norm(vectorB);     // �����
   double dif = norm;                  // �������, ������ � �� ��������� ������ norm, ����� ����� � ����

   // ������������ ����
   for (k = 1; (k <= maxIterations) && (dif > maxDif); k++)
   {
      dif = 0;
      // ����������� �� ���� ������� ������� �
      for (size_t i = 0; i < n; i++)
      {
         sum = 0;
         int64_t indX;
         // ����������� �� ���� �������� ������� �
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
   if (debugOutput)
   {
      cout << endl;
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

   auto& vec1 = vectorX;
   auto& vec2 = methodNum == 1 ? nextVectorX : vectorX;
   string methodName;
   switch (methodNum)
   {
      case 1:
         methodName = "�����";
         break;
      case 2:
         methodName = "�������";
         break;
      default:
         cout << "������� ������� ��������." << endl;
         return;
         break;
   }

   cout << endl << "������ ������������ ��� ������ " << methodName << endl << endl;

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

      cout << "������� ������� ����.����������: " << w << endl;

      size_t currIter = 0;
      double dif;
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
            bestVec = vector<double>(vec1);
         }
         if (currIter == minIter && dif < minDif)
         {
            minIterRelax = w;
            minDif = dif;
            bestVec = vector<double>(vec1);
         }
         cout << format("����� ��������: {0}\n���������� ������������� �������: {1:<15.3e}\n", currIter, dif);
         cout << "���������� ������ �������: " << endl;
         PrintArray(vec1.data(), vec1.size(), 15);
      }

      vec1 = vector<double>(defaultVec1);
      vec2 = vector<double>(defaultVec2);
   }

   cout << endl << "****************" << endl << endl;
   
   cout << "������������ ���������. �������� ��������� ����������: " << endl;
   cout << "���������� �������� �������� ������ �� ���� ���������: " << minIter << endl;
   cout << "�������� ���������� ��� ����:   " << minIterRelax << endl;
   cout << "�������� ������������� �������: " << minDif << endl;
   cout << "�������� ������ ���������� ������ �������: " << endl;
   PrintArray(bestVec.data(), bestVec.size(), 15);
}

void main()
{
   setlocale(LC_ALL, "ru-RU");

   ReadData();

   int operationNum;
   double dif;
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
         PrintArray<double>(vectorX.data(), vectorX.size(), 15, cout);

         PrintArray<double>(vectorX.data(), vectorX.size(), 15, ofstream("./iofiles/allFloatOutput.txt"));
         break;
      case 2:
         Iterations(vectorX, vectorX, dif, true);

         cout << "���������� ������ ������� �������: " << endl;
         PrintArray<double>(vectorX.data(), vectorX.size(), 15, cout);

         PrintArray<double>(vectorX.data(), vectorX.size(), 15, ofstream("./iofiles/allFloatOutput.txt"));
         break;

      case 4:
         RelaxationTester();
         break;
    }

}
