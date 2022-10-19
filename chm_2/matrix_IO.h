#pragma once
#include <fstream>
#include <iostream>
#include <vector>

/// <summary>
/// ������� ��������� ����� ���� size_t �� �����
/// </summary>
/// <param name="filePath"> - ���� �� ����� </param>
/// <returns> ����� �� ����� </returns>
inline size_t GetSizetFromFile(const char* filePath)
{
   auto fin = std::ifstream(filePath);
   size_t matrix_size = 0;
   fin >> matrix_size;
   fin.close();
   return matrix_size;
}

/// <summary>
/// ������� ���������� ������� ������ � �����
/// </summary>
/// <typeparam name="T"> - ��� ������, ������� ���������� ������� </typeparam>
/// <param name="filePath"> - ���� �� �����</param>
/// <param name="count"> - ����� ���������, ����������� ������� </param>
/// <returns> ��������� �� ������ ���������� ������� </returns>
template<typename T>
T* GetArrayFromFile(const char* filePath, const size_t count)
{
   T* arr = new T[count];
   auto fin = std::ifstream(filePath);
   for (int i = 0; i < count; i++)
   {
      fin >> arr[i];
   }
   fin.close();

   return arr;
}

/// <summary>
/// ������� ���������� ������� ������ � �����
/// </summary>
/// <typeparam name="T"> - ��� ������, ������� ���������� ������� </typeparam>
/// <param name="filePath"> - ���� �� �����</param>
/// <param name="count"> - ����� ���������, ����������� ������� </param>
/// <returns> ��������� ������ </returns>
template<typename T>
std::vector<T> GetVectorFromFile(const char* filePath, const size_t count)
{
   std::vector<T> arr(count);
   auto fin = std::ifstream(filePath);
   for (int i = 0; i < count; i++)
   {
      fin >> arr[i];
   }
   fin.close();

   return arr;
}

/// <summary>
/// ������� ���������� ������� ������ � �����
/// </summary>
/// <typeparam name="T"> - ��� ������, ������� ���������� ������� </typeparam>
/// <param name="file"> - �������������� �������� ����</param>
/// <param name="count"> - ����� ���������, ����������� ������� </param>
/// <returns> ��������� ������ </returns>
template<typename T>
std::vector<T> GetVectorFromFile(std::istream& file, const size_t count)
{
   std::vector<T> arr(count);
   for (int i = 0; i < count; i++)
   {
      file >> arr[i];
   }

   return arr;
}

/// <summary>
/// ����� ������� � ������������ �����
/// </summary>
/// <typeparam name="T"> - ��� ������ ��� ������</typeparam>
/// <param name="_arr"> - ��������� �� ������ ������� </param>
/// <param name="_size"> - ����� ������� </param>
/// <param name="_out">  - ����� ������ (cout �� ���������)</param>
template <typename T>
void PrintArray(const T* _arr, const size_t _size, const int _precision = 7, std::ostream& _out = std::cout)
{
   auto prevPrec = _out.precision();
   auto prevFlags = _out.flags();
   _out.precision(_precision);
   _out.setf(std::ios_base::fixed);

   for (int i = 0; i < _size; i++)
   {
      _out << _arr[i] << std::endl;
   }

   _out.precision(prevPrec);
   _out.flags(prevFlags);
}

/// <summary>
/// ����� ������� � ������������ �����
/// </summary>
/// <typeparam name="T"> - ��� ������ ��� ������</typeparam>
/// <param name="_arr"> - ��������� �� ������ ������� </param>
/// <param name="_size"> - ����� ������� </param>
/// <param name="_out">  - ����� ������ (cout �� ���������)</param>
template <typename T>
void PrintArray(const T* _arr, const size_t _size, const int _precision, std::ostream&& _out)
{
   auto prevPrec = _out.precision();
   auto prevFlags = _out.flags();
   _out.precision(_precision);
   _out.setf(std::ios_base::fixed);

   for (int i = 0; i < _size; i++)
   {
      _out << _arr[i] << std::endl;
   }

   _out.precision(prevPrec);
   _out.flags(prevFlags);
}
