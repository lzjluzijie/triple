#include "fingercodes-vt.hpp"

#include <coproto/Socket/AsioSocket.h>
#include <cryptoTools/Common/CLP.h>
#include <libOTe/Base/BaseOT.h>
#include <libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h>

#include "triple.hpp"

using namespace osuCrypto;

const string FingerCodesVT::defaultDB = "fingercodes.txt";
uint32_t FingerCodesVT::bitlen = FingerCodesVT::BITLEN;
using UINT = uint64_t;
using INT = int64_t;

void FingerCodesVT::gen(size_t num) {
  ofstream out(defaultDB);
  out << num << endl;
  for (size_t i = 0; i < num; i++) {
    auto data = randomUint8(dim);
    auto s = base64_encode(data);
    out << s << endl;
  }
}

vector<uint8_t> FingerCodesVT::Sample(size_t x) {
  vector<uint8_t> data(dim, 0);
  if (x < 256) {
    data = vector<uint8_t>(dim, x);
  } else if (x == 256) {
    for (size_t i = 0; i < dim; i++) {
      data[i] = i % 256;
    }
  }
  return data;
}

size_t FingerCodesVT::euclidianDistance(vector<uint8_t> x, vector<uint8_t> y) {
  size_t ed = 0;
  for (size_t i = 0; i < dim; i++) {
    ed += (x[i] - y[i]) * (x[i] - y[i]);
  }
  return ed;
}

inline u64 FingerCodesVT::square(vector<uint8_t> x, size_t mod) {
  u64 sum = 0;
  for (size_t i = 0; i < dim; i++) {
    sum += x[i] * x[i];
    sum = sum % mod;
  }
  return sum;
}

FingerCodesServerVT::FingerCodesServerVT(string db) {
  log = ofstream("server-log.txt");
  if (db == "") {
    db = defaultDB;
  }
  PRNG prng(sysRandomSeed());

  vector<vector<uint8_t>> data;
  ifstream in(db);
  if (in.good()) {
    string line;
    in >> n;
    for (size_t j = 0; j < n; j++) {
      in >> line;
      data.push_back(base64_decode(line));
    }
  } else {
    for (size_t j = 0; j < n; j++) {
      data.push_back(randomUint8(dim));
    }
  }

  auto chl = cp::asioConnect("127.0.0.1:7700", true);
  cp::sync_wait(chl.send(n));
  cp::sync_wait(chl.recv(nc));

  vector<vector<u32>> x0(nc, vector<u32>(dim, 0));
  vector<vector<u32>> y0(dim, vector<u32>(n));
  for (size_t k = 0; k < dim; k++) {
    for (size_t i = 0; i < n; i++) {
      y0[k][i] = data[i][k];
    }
  }

  auto triple = fakeZ<u32>(nc, n, dim);

  vector<u32> R(n, 0);

  auto otStart = chrono::high_resolution_clock::now();
  auto mpShare = mp0<u32, u32>(chl, triple.first, x0, y0);
  auto otEnd = chrono::high_resolution_clock::now();

  u64 traffic = chl.bytesSent() + chl.bytesReceived();
  cout << "Matrix product time: "
       << chrono::duration_cast<chrono::microseconds>(otEnd - otStart).count()
       << "us with transferred "  << traffic << " "
       << (traffic >> 10) << "KiB" << endl;

  for (u64 j = 0; j < n; j++) {
    R[j] = square(data[j], mod) - 2 * mpShare[0][j];
  }

  vector<u32> uR(n);
  for (u64 j = 0; j < n; j++) {
    uR[j] = (R[j]);
  }
//  if (MODE == 0) {
//    ABY aby(mod, threshold, bitlen, circType);
//    auto res = aby.AddMatchingR(uR);
//    cout << "Final matching: " << (res ? "True" : "False") << endl;
//  } else if (MODE == 1) {
//    ABY aby(mod, threshold, bitlen, circType);
//    auto res = aby.EveryMatchingR(uR);
//    for (size_t i = 0; i < n; i++) {
//      cout << "Result of finger " << i << ": "
//          << (res[i] == 1 ? "True" : "False") << endl;
//    }
//  } else if (MODE == 2) {
//    ABY aby(mod, threshold, bitlen, circType);
//    auto res = aby.AddEveryDistanceR(uR);
//    for (size_t i = 0; i < n; i++) {
//      cout << "ED of finger " << i << ": " << res[i] << endl;
//    }
//  }
}

FingerCodesClientVT::FingerCodesClientVT(string db) {
  log = ofstream("client-log.txt");
  PRNG prng(sysRandomSeed());

  if (db == "") {
    db = defaultDB;
  }

  vector<vector<uint8_t>> data;
  ifstream in(db);
  if (in.good()) {
    string line;
    in >> nc;
    for (size_t j = 0; j < nc; j++) {
      in >> line;
      data.push_back(base64_decode(line));
    }
  } else {
    for (size_t j = 0; j < nc; j++) {
      data.push_back(randomUint8(dim));
    }
  }

  auto chl = cp::asioConnect("127.0.0.1:7700", false);
  IknpOtExtReceiver receiver;
  cp::sync_wait(chl.recv(n));
  cp::sync_wait(chl.send(nc));

  vector<u32> T(n);

  vector<vector<u32>> x0(nc, vector<u32>(dim));
  for (size_t k = 0; k < dim; k++) {
    for (size_t i = 0; i < nc; i++) {
      x0[i][k] = data[i][k];
    }
  }
  vector<vector<u32>> y0(dim, vector<u32>(n, 0));

  auto triple = fakeZ<u32>(nc, n, dim);

  auto mpShare = mp1<u32, u32>(chl, triple.second, x0, y0);

  cp::sync_wait(chl.flush());

  for (u64 j = 0; j < n; j++) {
    T[j] = square(data[0], mod) - 2 * mpShare[0][j];
  }

  vector<u32> uT(n);
  for (u64 j = 0; j < n; j++) {
    uT[j] = (T[j]);
  }
//  if (MODE == 0) {
//    ABY aby(mod, threshold, bitlen, circType);
//    auto res = aby.AddMatchingT(uT);
//    cout << "Final matching: " << (res ? "True" : "False") << endl;
//  } else if (MODE == 1) {
//    ABY aby(mod, threshold, bitlen, circType);
//    auto res = aby.EveryMatchingT(uT);
//    for (size_t i = 0; i < n; i++) {
//      cout << "Result of finger " << i << ": "
//          << (res[i] == 1 ? "True" : "False") << endl;
//    }
//  } else if (MODE == 2) {
//    ABY aby(mod, threshold, bitlen, circType);
//    auto res = aby.AddEveryDistanceT(uT);
//    for (size_t i = 0; i < n; i++) {
//      cout << "ED of finger " << i << ": " << res[i] << endl;
//    }
//  }
}
