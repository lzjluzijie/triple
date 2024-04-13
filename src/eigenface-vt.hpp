#ifndef EIGENFACE_VT_HPP
#define EIGENFACE_VT_HPP

#include "triple.hpp"
#include "mp.hpp"

using namespace std;

class EigenFaceVT {
 public:
  const static size_t spDim = 10304;
  const static size_t edDim = 12;
  const static int threshold = 200;
  const static size_t mod = -1;
  const static uint32_t BITLEN = 64;
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
  int MODE = 2;

  static uint64_t square(vector<uint8_t> x);
  static size_t euclidianDistance(vector<uint8_t> x, vector<uint8_t> y);
};

class EigenFaceServerVT : EigenFaceVT {
 public:
  EigenFaceServerVT(string db);
  static void gen(size_t n);
};

class EigenFaceClientVT : EigenFaceVT {
 public:
  EigenFaceClientVT(string s);
  static void gen(size_t n);
};

#endif