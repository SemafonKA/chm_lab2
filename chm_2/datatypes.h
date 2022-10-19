#include <cmath>
#include <string>

// Выбор точности работы программы:
// __ALL_FLOAT__     - все дробные числа и счётчики в float, 
// __SEMI_FLOAT__    - все дробные числа в float, счётчики в double, 
//__ALL_DOUBLE__     - все дробные числа и счётчики в double
#define __ALL_DOUBLE__

#ifdef __ALL_FLOAT__
   using real_t = float;
   using accum_t = float;
   constexpr real_t g_eps = FLT_EPSILON * 10.0;
   constexpr int g_coutPrecision = 7;
   const std::string g_outputFileName = "./iofiles/allFloatOutput.txt";
#else 
   #ifdef __SEMI_FLOAT__
      using real_t = float;
      using accum_t = double;
      constexpr real_t g_eps = FLT_EPSILON * 10.0;
      constexpr int g_coutPrecision = 7;
      const std::string g_outputFileName = "./iofiles/scalarDoubleOutput.txt";
    #else // __ALL_DOUBLE__
      using real_t = double;
      using accum_t = double;
      constexpr real_t g_eps = DBL_EPSILON * 10.0;
      constexpr int g_coutPrecision = 15;
      const std::string g_outputFileName = "./iofiles/allDoubleOutput.txt";
    #endif 
#endif

inline bool IsAlmostEq(real_t _first, real_t _second, real_t maxRelDif = g_eps) {
   real_t delta = abs(_first - _second);
   _first = abs(_first); _second = abs(_second);
   real_t maxNum = _first > _second ? _first : _second;
   return delta <= maxNum * maxRelDif;
}