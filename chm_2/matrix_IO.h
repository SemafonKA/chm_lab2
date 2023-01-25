#pragma once
#include <fstream>
#include <iostream>
#include <vector>

/// <summary>
/// Функция считывания вектора данных с файла
/// </summary>
/// <typeparam name="T"> - тип данных, который необходимо считать </typeparam>
/// <param name="filePath"> - путь до файла</param>
/// <param name="count"> - число элементов, необходимое считать </param>
/// <returns> считанный вектор </returns>
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
/// Функция считывания вектора данных с файла
/// </summary>
/// <typeparam name="T"> - тип данных, который необходимо считать </typeparam>
/// <param name="file"> - предварительно открытый файл</param>
/// <param name="count"> - число элементов, необходимое считать </param>
/// <returns> считанный вектор </returns>
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
/// Вывод массива в произвольный поток
/// </summary>
/// <typeparam name="T"> - тип данных для вывода</typeparam>
/// <param name="_arr"> - указатель на начало массива </param>
/// <param name="_size"> - длина массива </param>
/// <param name="_out">  - поток вывода (cout по умолчанию)</param>
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
/// Вывод массива в произвольный поток
/// </summary>
/// <typeparam name="T"> - тип данных для вывода</typeparam>
/// <param name="_arr"> - указатель на начало массива </param>
/// <param name="_size"> - длина массива </param>
/// <param name="_out">  - поток вывода (cout по умолчанию)</param>
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
