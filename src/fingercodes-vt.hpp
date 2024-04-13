#ifndef FINGERCODES_VT_HPP
#define FINGERCODES_VT_HPP

#include "triple.hpp"
#include "mp.hpp"

using namespace std;

class FingerCodesVT {
 public:
  const static size_t dim = 640;
  const static int threshold = 20;
  const static size_t mod = 255 * 255 * dim;
  const static uint32_t BITLEN = 32;
  static uint32_t bitlen;
  const static string defaultDB;
//  const static auto circType = S_BOOL;
  // const static auto circType = S_YAO;

  ofstream log;

  size_t n = 100;
  size_t nc = 1;

  /*
  0: only final matching
  1: matching per face
  2: hamming distance per face
  */
  int MODE = 0;

  static uint64_t square(vector<uint8_t> x, size_t mod);
  static size_t euclidianDistance(vector<uint8_t> x, vector<uint8_t> y);
  static vector<uint8_t> Sample(size_t x);

  static void gen(size_t n);
};

class FingerCodesServerVT : FingerCodesVT {
 public:
  FingerCodesServerVT(string db);
};

class FingerCodesClientVT : FingerCodesVT {
 public:
  FingerCodesClientVT(string s);
};

#endif