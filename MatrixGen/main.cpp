#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>

using namespace std;

void MatrixGen(vector<vector<double>>& matrix, size_t zeroDiags) {
   auto rand = mt19937(time(nullptr));
   auto getRandNumber = [&]() -> double { return rand() % 100 / 10.0; };

   for (size_t i = 0; i < matrix.size(); i++)
   {
      if (i < matrix.size() - 1)
      {
         matrix[i + 1][i] = getRandNumber();
      }

      if (i >= 1)
      {
         matrix[i - 1][i] = getRandNumber();
      }

      for (size_t k = 1; k < 4; k++)
      {
         if (i < matrix.size() - zeroDiags - 1 - k)
         {
            matrix[i + zeroDiags + k + 1][i] = getRandNumber();
         }
         if (i >= 1 + zeroDiags + k)
         {
            matrix[i - 1 - zeroDiags - k][i] = getRandNumber();
         }
      }
   }

   for (size_t i = 0; i < matrix.size(); i++)
   {
      for (size_t k = 0; k < matrix.size(); k++)
      {
         if (k == i) continue;
         matrix[i][i] += matrix[i][k];
      }
   }
}

void GetVectorB(vector<vector<double>>& mas, vector<double>& x, vector<double>& b) {
   for (int i = 0; i < b.size(); i++)
   {
      for (int j = 0; j < x.size(); j++)
      {
         b[i] += mas[i][j] * x[j];
      }
   }
}

void MatrixOut(vector<vector<double>>& matrix, size_t zeroDiags, ostream&& out) {
   // ������� ��������� � ������� ����
   out.precision(15);
   out.setf(std::ios::scientific);

   long long lZeroDiags = static_cast<long long>(zeroDiags);
   vector<long long> diagsShift = {4 + lZeroDiags, 3 + lZeroDiags, 2 + lZeroDiags, 1, 0, -1, -2 - lZeroDiags, -3 - lZeroDiags, -4 - lZeroDiags};
   for (size_t k = 0; k < 9; k++)
   {
      long long i = 0, j = diagsShift[k];
      while (j < 0)
      {
         j++; i++;
      }
      if (i > 0)
      {
         for (int n = i; n > 0; n--)
         {
            out << 0. << " ";
         }
      }
      while (j < matrix.size() && i < matrix.size())
      {
         out << matrix[i][j] << " ";
         j++; i++;
      }
      while (i++ < matrix.size())
      {
         out << 0. << " ";
      }
      out << endl;
   }

}

void VectorOut(vector<double>& vec, ostream&& out) {
   out.precision(15);
   out.setf(ios::scientific);

   for (auto elem : vec)
   {
      out << elem << " ";
   }
   out << endl;
}

int main(int argc, char** argv) {
   setlocale(LC_ALL, "ru-RU");

   size_t matrixSize;
   size_t zeroDiags;
   if (argc >= 3)
   {
      stringstream ss;
      ss << argv[1];
      ss >> matrixSize;
      ss << argv[2];
      ss >> zeroDiags;
   }
   else
   {
      cout << "������� ������ ������������ �������: ";
      cin >> matrixSize;
      cout << "������� ����� ������� ����������: ";
      cin >> zeroDiags;
   }

   if (zeroDiags > matrixSize - 5)
   {
      cout << "������� ������� ����� ������� ����������." << endl;
      cerr << "Fatal error: wrong number of zero diags." << endl;
      return -1;
   }

   vector<vector<double>> matrixA(matrixSize);
   vector<double> vectorX(matrixSize);
   vector<double> vectorB(matrixSize);

   for (auto& elem : matrixA)
   {
      elem.resize(matrixSize);
   }

   for (int i = 0; i < vectorX.size(); i++)
   {
      vectorX[i] = i + 1;
   }

   MatrixGen(matrixA, zeroDiags);
   GetVectorB(matrixA, vectorX, vectorB);

   MatrixOut(matrixA, zeroDiags, ofstream("./iofiles/matrixA.txt"));
   VectorOut(vectorB, ofstream("./iofiles/vectorB.txt"));
   VectorOut(vectorX, ofstream("./iofiles/absoluteVectorX.txt"));

   auto matrixParamsFile = ofstream("./iofiles/matrixParams.txt");
   matrixParamsFile << matrixA.size() << " " << zeroDiags << endl;
   matrixParamsFile.close();
   return 0;
}