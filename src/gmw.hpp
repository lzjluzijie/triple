#ifndef GMW_HPP
#define GMW_HPP

#include <libOTe/Base/BaseOT.h>

#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace std;
using namespace osuCrypto;

using Triple = array<vector<u8>, 3>;
using BV = BitVector;

class GMW {
 public:
  u64 bitlen = 16;

  template<typename D>
  vector<u8> MatchingR(coproto::Socket chl, vector<D> R, Triple triple) {
    size_t n = R.size();
    vector<u8> carry(n, 0);
    BV d0(n);
    BV d1(n);
    BV e0(n);
    BV e1(n);

    for (u64 k = 0; k < bitlen; k++) {
      for (u64 i = 0; i < n; i++) {
        u8 x1c = ((0 >> k) & 1) ^ carry[i];
        u8 y1c = ((R[i] >> k) & 1) ^ carry[i];
        d1[i] = x1c ^ triple[0][k * n + i];
        e1[i] = y1c ^ triple[1][k * n + i];
      }
      auto d0Task = chl.recv(d0);
      auto e0Task = chl.recv(e0);
      auto d1Task = chl.send(d1);
      auto e1Task = chl.send(e1);
      cp::sync_wait(d0Task);
      cp::sync_wait(e0Task);
      cp::sync_wait(d1Task);
      cp::sync_wait(e1Task);
      // auto res = macoro::sync_wait(
      //     macoro::when_all_ready(d0Task, e0Task, d1Task, e1Task));

      for (u64 i = 0; i < n; i++) {
        carry[i] = carry[i] ^ ((d0[i] ^ d1[i]) & triple[1][k * n + i]) ^
            ((e0[i] ^ e1[i]) & triple[0][k * n + i]) ^
            triple[2][k * n + i];
      }
    }

    return carry;
  }

  template<typename D>
  vector<u8> MatchingT(coproto::Socket chl, vector<D> T, Triple triple) {
    size_t n = T.size();
    vector<u8> carry(n, 0);
    BV d0(n);
    BV d1(n);
    BV e0(n);
    BV e1(n);

    auto abyStart = chrono::high_resolution_clock::now();

    for (u64 k = 0; k < bitlen; k++) {
      for (u64 i = 0; i < n; i++) {
        u8 x0c = ((T[i] >> k) & 1) ^ carry[i];
        u8 y0c = ((0 >> k) & 1) ^ carry[i];
        d0[i] = x0c ^ triple[0][k * n + i];
        e0[i] = y0c ^ triple[1][k * n + i];
      }

      auto d0Task = chl.send(d0);
      auto e0Task = chl.send(e0);
      auto d1Task = chl.recv(d1);
      auto e1Task = chl.recv(e1);
      cp::sync_wait(d0Task);
      cp::sync_wait(e0Task);
      cp::sync_wait(d1Task);
      cp::sync_wait(e1Task);
      // auto res = macoro::sync_wait(
      // macoro::when_all_ready(d0Task, e0Task, d1Task, e1Task));

      for (u64 i = 0; i < n; i++) {
        carry[i] = carry[i] ^ ((d0[i] ^ d1[i]) & (e0[i] ^ e1[i])) ^
            ((d0[i] ^ d1[i]) & triple[1][k * n + i]) ^
            ((e0[i] ^ e1[i]) & triple[0][k * n + i]) ^
            triple[2][k * n + i];
      }
    }

    auto abyEnd = chrono::high_resolution_clock::now();
    cerr << "GMW time: "
         << chrono::duration_cast<chrono::milliseconds>(abyEnd - abyStart)
             .count()
         << "ms" << endl;
    return carry;
  }
};

#endif
