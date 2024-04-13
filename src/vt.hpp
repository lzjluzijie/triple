#ifndef VT_HPP
#define VT_HPP

#include <iostream>
#include <coproto/Socket/AsioSocket.h>
#include "libOTe/Vole/Noisy/NoisyVoleSender.h"
#include "libOTe/Vole/Noisy/NoisyVoleReceiver.h"
#include "libOTe/Vole/Silent/SilentVoleSender.h"
#include "libOTe/Vole/Silent/SilentVoleReceiver.h"
#include "libOTe/Vole/Subfield/NoisyVoleSender.h"
#include "libOTe/Vole/Subfield/NoisyVoleReceiver.h"
#include "libOTe/Vole/Subfield/SilentVoleSender.h"
#include "libOTe/Vole/Subfield/SilentVoleReceiver.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h"
#include "cryptoTools/Network/Session.h"
#include "cryptoTools/Network/IOService.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Common/Range.h"
#include "cryptoTools/Common/TestCollection.h"
#include "coproto/Socket/BufferingSocket.h"
#include "libOTe/Tools/Subfield/Subfield.h"

using namespace osuCrypto;
using namespace osuCrypto::Subfield;
using namespace std;

const bool verifyVT = false;

template<typename T, size_t N>
void noisyBench(u64 m, u64 n) {
  u64 iterations = 5;
  Timer timer;
  timer.setTimePoint("start");
  block seed = block(0, 0);
//  block seed = sysRandomSeed();
  PRNG prng(seed);

  using TypeTrait = TypeTraitVec<u32, N>;
  u64 bitsF = TypeTrait::bitsF;
  using F = typename TypeTrait::F;
  using G = typename TypeTrait::G;

  std::vector<F> d0(m);
  prng.get(d0.data(), d0.size());
  std::vector<G> y1(m * n);
  prng.get(y1.data(), y1.size());
  std::vector<F> s0(m * n), r1(m * n);

  std::vector<F> d1(m);
  prng.get(d1.data(), d1.size());
  std::vector<G> y0(m * n);
  prng.get(y0.data(), y0.size());
  std::vector<F> s1(m * n), r0(m * n);

  u64 traffic = 0;
  vector<u64> times(iterations);

  AlignedVector<std::array<block, 2>> otSendMsg0(bitsF * m);
  BitVector choice0((u8 *) d0.data(), bitsF * m);
  AlignedVector<block> otRecvMsg0(bitsF * m);

  AlignedVector<std::array<block, 2>> otSendMsg1(bitsF * m);
  BitVector choice1((u8 *) d1.data(), bitsF * m);
  AlignedVector<block> otRecvMsg1(bitsF * m);

  for (u64 count = 0; count < iterations; count++) {
    auto start = std::chrono::high_resolution_clock::now();
    std::thread t0 = thread([&]() {
      auto chl0 = cp::asioConnect("localhost:7700", true);

      block seed0 = sysRandomSeed();
      block seed1 = sysRandomSeed();
      PRNG prng0(seed0);
      PRNG prng1(seed1);

      IknpOtExtReceiver otr0;
      IknpOtExtSender ots0;

      timer.setTimePoint("ot0");
      auto f0 = chl0.fork();
      auto f1 = chl0.fork();
      auto task0 = otr0.receive(choice0, otRecvMsg0, prng0, f0);
      auto task1 = ots0.send(otSendMsg0, prng1, f1);
      macoro::sync_wait(macoro::when_all_ready(std::move(task0), std::move(task1)));
//      macoro::sync_wait(std::move(task0));
//      macoro::sync_wait(std::move(task1));

      timer.setTimePoint("vole0");
      NoisySubfieldVoleReceiver<TypeTrait> recv;
      NoisySubfieldVoleSender<TypeTrait> send;
      vector<macoro::eager_task<>> tasks0;
      tasks0.reserve(2 * m);
      vector<cp::Socket> chls(2 * m);
      for (u64 i = 0; i < 2 * m; i++) {
        chls[i] = chl0.fork();
      }
      for (u64 k = 0; k < m; k++) {
//        cout << k << endl;
        {
          auto fork0 = chl0.fork();
          auto kRecvMsg = oc::span<block>(otRecvMsg0.data() + k * bitsF, bitsF);
          auto ks0 = oc::span<F>(s0.data() + k * n, n);
          auto p0 = send.send(d0[k], ks0, prng0, kRecvMsg, chls[2 * k]);
          tasks0.emplace_back(std::move(p0) | macoro::make_eager());
        }

        {
          auto fork1 = chl0.fork();
          auto kSendMsg = oc::span<std::array<block, 2>>(otSendMsg0.data() + k * bitsF, bitsF);
          auto ky = oc::span<G>(y1.data() + k * n, n);
          auto kr1 = oc::span<F>(r1.data() + k * n, n);
          auto p1 = recv.receive(ky, kr1, prng1, kSendMsg, chls[2 * k + 1]);
          tasks0.emplace_back(std::move(p1) | macoro::make_eager());
        }
      }

      for (u64 i = 0; i < tasks0.size(); i++) {
//        timer.setTimePoint("0 waiting " + std::to_string(i));
        macoro::sync_wait(tasks0[i]);
      }

      timer.setTimePoint("compute triple0");
      for (u64 k = 0; k < m; k++) {
        for (u64 i = 0; i < n; ++i) {
          s0[k * n + i] = d0[k] * y1[k * n + i] + s0[k * n + i] - r1[k * n + i];
        }
      }

      timer.setTimePoint("done0");

      traffic += chl0.bytesSent() + chl0.bytesReceived();
      macoro::sync_wait(chl0.flush());
    });

    std::thread t1 = thread([&]() {
      auto chl1 = cp::asioConnect("localhost:7700", false);
      block seed0 = sysRandomSeed();
      block seed1 = sysRandomSeed();
      PRNG prng0(seed0);
      PRNG prng1(seed1);

      IknpOtExtSender ots;
      IknpOtExtReceiver otr;

      timer.setTimePoint("ot1");
      auto f0 = chl1.fork();
      auto f1 = chl1.fork();
      auto task0 = ots.send(otSendMsg1, prng0, f0);
      auto task1 = otr.receive(choice1, otRecvMsg1, prng1, f1);
      macoro::sync_wait(macoro::when_all_ready(std::move(task0), std::move(task1)));
//      macoro::sync_wait(std::move(task0));
//      macoro::sync_wait(std::move(task1));

      timer.setTimePoint("vole1");
      NoisySubfieldVoleReceiver<TypeTrait> recv;
      NoisySubfieldVoleSender<TypeTrait> send;
      std::vector<macoro::eager_task<>> tasks0;
      tasks0.reserve(2 * m);
      vector<cp::Socket> chls(2 * m);
      for (u64 i = 0; i < 2 * m; i++) {
        chls[i] = chl1.fork();
      }
      for (u64 k = 0; k < m; k++) {
        {
          auto fork0 = chl1.fork();
          auto kSendMsg = oc::span<std::array<block, 2>>(otSendMsg1.data() + k * bitsF, bitsF);
          auto ky = oc::span<G>(y0.data() + k * n, n);
          auto kr0 = oc::span<F>(r0.data() + k * n, n);
          auto p0 = recv.receive(ky, kr0, prng0, kSendMsg, chls[2 * k]);
          tasks0.emplace_back(std::move(p0) | macoro::make_eager());
        }
        {
          auto fork1 = chl1.fork();
          auto kRecvMsg = oc::span<block>(otRecvMsg1.data() + k * bitsF, bitsF);
          auto ks1 = oc::span<F>(s1.data() + k * n, n);
          auto p1 = send.send(d1[k], ks1, prng1, kRecvMsg, chls[2 * k + 1]);
          tasks0.emplace_back(std::move(p1) | macoro::make_eager());
        }
      }

      for (u64 i = 0; i < tasks0.size(); i++) {
//        timer.setTimePoint("1 waiting " + std::to_string(i));
        macoro::sync_wait(tasks0[i]);
      }

      timer.setTimePoint("compute triple1");
      for (u64 k = 0; k < m; k++) {
        for (u64 i = 0; i < n; ++i) {
          s1[k * n + i] = d1[k] * y0[k * n + i] + s1[k * n + i] - r0[k * n + i];
        }
      }
      timer.setTimePoint("done1");

      macoro::sync_wait(chl1.flush());
    });

    t0.join();
    t1.join();
    auto end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    times[count] = time;
  }

  traffic /= iterations;
  cout << "time: ";
  for (auto t : times) {
    cout << t << " ";
  }
  sort(times.begin(), times.end());
  u64 medianTime = times[times.size() / 2];
  cout << " median " << medianTime << endl;
  cout << "transferred " << traffic
       << " " << ((double) traffic / (1 << 10)) << "KiB " << ((double) traffic / (1 << 20)) << "MiB " << endl;
//  std::cout << timer;

//  for (u64 k = 0; k < m; k++) {
//    for (u64 i = 0; i < n; ++i) {
//      for (u64 j = 0; j < N; j++) {
//        G left = d0[k][j] * y0[k * n + i];
//        G right = s0[k * n + i][j] - r0[k * n + i][j];
//        if (left != right) {
//          throw RTE_LOC;
//        }
//        left = d1[k][j] * y1[k * n + i];
//        right = s1[k * n + i][j] - r1[k * n + i][j];
//        if (left != right) {
//          throw RTE_LOC;
//        }
//      }
//    }
//    cout << "OK" << endl;
//  }

  if (verifyVT) {
    for (u64 k = 0; k < m; k++) {
      for (u64 i = 0; i < n; i++) {
        for (u64 j = 0; j < N; j++) {
          F left = (d0[k] + d1[k]) * (y0[k * n + i] + y1[k * n + i]);
          F right = s0[k * n + i] + s1[k * n + i];
          if (left != right) {
            throw RTE_LOC;
          }
        }
      }
    }
  }
}

//template<typename T, size_t N>
//void silentBench(u64 m, u64 n) {
//  u64 iterations = 5;
//  vector<u64> times(iterations);
//  u64 traffic = 0;
//
//  cout << "time: ";
//  for (auto t : times) {
//    cout << t << " ";
//  }
//  sort(times.begin(), times.end());
//  u64 medianTime = times[times.size() / 2];
//  cout << " median " << medianTime << endl;
//  cout << "transferred " << traffic
//       << " " << ((double) traffic / (1 << 10)) << "KiB " << ((double) traffic / (1 << 20)) << "MiB " << endl;
//}

template<typename T, size_t N>
void silentBench(u64 m, u64 n) {
  u64 iterations = 5;
  Timer timer;
  timer.setTimePoint("start");
  block seed = block(0, 0);
  PRNG prng(seed);

  using TypeTrait = TypeTraitVec<T, N>;
  using F = typename TypeTrait::F;
  using G = typename TypeTrait::G;

  std::vector<F> d0(m);
  prng.get(d0.data(), d0.size());
  std::vector<G> y1(m * n);
  prng.get(y1.data(), y1.size());
  std::vector<F> s0(m * n), r1(m * n);

  std::vector<F> d1(m);
  prng.get(d1.data(), d1.size());
  std::vector<G> y0(m * n);
  prng.get(y0.data(), y0.size());
  std::vector<F> s1(m * n), r0(m * n);

//  auto chls = cp::LocalAsyncSocket::makePair();
  u64 traffic = 0;
  vector<u64> times(iterations);

  for (u64 count = 0; count < iterations; count++) {
    auto start = std::chrono::high_resolution_clock::now();

    std::thread t0 = thread([&]() {
      vector<PRNG> prng0(m);
      for (u64 k = 0; k < m; k++) {
        prng0[k].SetSeed(block(00, k));
      }
      vector<PRNG> prng1(m);
      for (u64 k = 0; k < m; k++) {
        prng1[k].SetSeed(block(01, k));
      }

      timer.setTimePoint("vole0");

//      auto chl0 = cp::asioConnect("localhost:7700", true);
      vector<SilentSubfieldVoleReceiver<TypeTrait>> recv(m);
      vector<SilentSubfieldVoleSender<TypeTrait>> send(m);
      vector<cp::Socket> chls(2 * m);
      for (u64 i = 0; i < 2 * m; i++) {
//        chls[i] = chl0.fork();
      }

      std::vector<macoro::eager_task<>> tasks0;
      tasks0.reserve(2 * m);
      for (u64 k = 0; k < m; k++) {
        {
          chls[2 * k] = cp::asioConnect("localhost:" + to_string(7800 + 2 * k), true);
          auto ks0 = oc::span<F>(s0.data() + k * n, n);
          auto p0 = send[k].silentSend(d0[k], ks0, prng0[k], chls[2 * k]);
          tasks0.emplace_back(std::move(p0) | macoro::make_eager());
        }
        {
          chls[2 * k + 1] = cp::asioConnect("localhost:" + to_string(7800 + 2 * k + 1), true);
          auto ky = oc::span<G>(y1.data() + k * n, n);
          auto kr1 = oc::span<F>(r1.data() + k * n, n);
          auto p1 = recv[k].silentReceive(ky, kr1, prng1[k], chls[2 * k + 1]);
          tasks0.emplace_back(std::move(p1) | macoro::make_eager());
        }
      }
      for (u64 i = 0; i < tasks0.size(); i++) {
        macoro::sync_wait(tasks0[i]);
      }

//      traffic = chl0.bytesSent() + chl0.bytesReceived();
//      macoro::sync_wait(chl0.flush());
      for (u64 i = 0; i < 2 * m; i++) {
        macoro::sync_wait(chls[i].flush());
        traffic += chls[i].bytesSent() + chls[i].bytesReceived();
      }

      timer.setTimePoint("compute triple0");
      for (u64 k = 0; k < m; k++) {
        for (u64 i = 0; i < n; ++i) {
          s0[k * n + i] = d0[k] * y1[k * n + i] + r1[k * n + i] - s0[k * n + i];
        }
      }

      timer.setTimePoint("done0");
    });

    std::thread t1 = thread([&]() {
      vector<PRNG> prng0(m);
      for (u64 k = 0; k < m; k++) {
        prng0[k].SetSeed(block(10, k));
      }
      vector<PRNG> prng1(m);
      for (u64 k = 0; k < m; k++) {
        prng1[k].SetSeed(block(11, k));
      }

      timer.setTimePoint("vole1");
      vector<SilentSubfieldVoleReceiver<TypeTrait>> recv(m);
      vector<SilentSubfieldVoleSender<TypeTrait>> send(m);
//      auto chl1 = cp::asioConnect("localhost:7700", false);
      vector<cp::Socket> chls(2 * m);
      for (u64 i = 0; i < 2 * m; i++) {
//        chls[i] = chl1.fork();
//        chls[i] = cp::asioConnect("localhost:" + to_string(7800+i), false);
      }
      std::vector<macoro::eager_task<>> tasks0;
      tasks0.reserve(2 * m);
      for (u64 k = 0; k < m; k++) {
        {
          chls[2 * k] = cp::asioConnect("localhost:" + to_string(7800 + 2 * k), false);
          auto ky = oc::span<G>(y0.data() + k * n, n);
          auto kr0 = oc::span<F>(r0.data() + k * n, n);
          auto p0 = recv[k].silentReceive(ky, kr0, prng0[k], chls[2 * k]);
          tasks0.emplace_back(std::move(p0) | macoro::make_eager());
        }
        {
          chls[2 * k + 1] = cp::asioConnect("localhost:" + to_string(7800 + 2 * k + 1), false);
          auto ks1 = oc::span<F>(s1.data() + k * n, n);
          auto p1 = send[k].silentSend(d1[k], ks1, prng1[k], chls[2 * k + 1]);
          tasks0.emplace_back(std::move(p1) | macoro::make_eager());
        }
      }
      for (u64 i = 0; i < tasks0.size(); i++) {
        macoro::sync_wait(tasks0[i]);
      }

//      macoro::sync_wait(chl1.flush());
      for (u64 i = 0; i < 2 * m; i++) {
        macoro::sync_wait(chls[i].flush());
      }

      timer.setTimePoint("compute triple1");
      for (u64 k = 0; k < m; k++) {
        for (u64 i = 0; i < n; ++i) {
          s1[k * n + i] = d1[k] * y0[k * n + i] + r0[k * n + i] - s1[k * n + i];
        }
      }

      timer.setTimePoint("done1");
    });

    t0.join();
    t1.join();
    auto end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    times[count] = time;

//  std::cout << timer;
//  std::cout << "transferred " << traffic
//            << " " << (traffic >> 10) << "KiB " << (traffic >> 20) << "MiB " << std::endl;

//  for (u64 k = 0; k < m; k++) {
//    for (u64 i = 0; i < n; ++i) {
//      for (u64 j = 0; j < N; j++) {
//        G left = d0[k][j] * y0[k * n + i];
//        G right = r0[k * n + i][j] - s0[k * n + i][j];
//        if (left != right) {
//          throw RTE_LOC;
//        }
//        left = d1[k][j] * y1[k * n + i];
//        right = r1[k * n + i][j] - s1[k * n + i][j];
//        if (left != right) {
//          throw RTE_LOC;
//        }
//      }
//    }
//  }

    if (verifyVT) {
      for (u64 k = 0; k < m; k++) {
        for (u64 i = 0; i < n; i++) {
          for (u64 j = 0; j < N; j++) {
            F left = (d0[k] + d1[k]) * (y0[k * n + i] + y1[k * n + i]);
            F right = s0[k * n + i] + s1[k * n + i];
            if (left != right) {
              throw RTE_LOC;
            }
          }
        }
      }
    }
  }

  cout << "time: ";
  for (auto t : times) {
    cout << t << " ";
  }
  sort(times.begin(), times.end());
  u64 medianTime = times[times.size() / 2];
  cout << " median " << medianTime << endl;
  traffic /= iterations;
  cout << "transferred " << traffic
       << " " << ((double) traffic / (1 << 10)) << "KiB " << ((double) traffic / (1 << 20)) << "MiB " << endl;
}

template<typename T>
pair<array<vector<T>, 3>, array<vector<T>, 3>> fakeTriple(u64 n, block seed = ZeroBlock) {
  PRNG prng(seed);
  vector<T> a0(n);
  vector<T> a1(n);
  vector<T> b0(n);
  vector<T> b1(n);
  vector<T> c0(n);
  vector<T> c1(n);
  prng.get(a0.data(), a0.size());
  prng.get(a1.data(), a1.size());
  prng.get(b0.data(), b0.size());
  prng.get(b1.data(), b1.size());
  prng.get(c0.data(), c0.size());
  for (size_t i = 0; i < n; i++) {
    c1[i] = (a0[i] + a1[i]) * (b0[i] + b1[i]) - c0[i];
  }
  return pair<array<vector<T>, 3>, array<vector<T>, 3>>({a0, b0, c0}, {a1, b1, c1});
}

template<typename T>
vector<T> mul0(cp::Socket chl, array<vector<T>, 3> &t, vector<T> &x0) {
  u64 n = x0.size();
  vector<T> d0(n);
  vector<T> e0(n);
  vector<T> d(n);
  vector<T> e(n);
  for (size_t i = 0; i < n; i++) {
    d0[i] = x0[i] - t[0][i];
    e0[i] = -t[1][i];
  }
  auto td0 = chl.send(d0);
  auto te0 = chl.send(e0);
  auto td1 = chl.recv(d);
  auto te1 = chl.recv(e);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
  for (size_t i = 0; i < n; i++) {
    d[i] = d0[i] + d[i];
    e[i] = e0[i] + e[i];
  }
  vector<T> res(n);
  for (size_t i = 0; i < n; i++) {
    res[i] = d[i] * e[i] + d[i] * t[1][i] + t[0][i] * e[i] + t[2][i];
  }
  return res;
}

template<typename T>
vector<T> mul1(cp::Socket chl, array<vector<T>, 3> &t, vector<T> &x1) {
  u64 n = x1.size();
  vector<T> d(n);
  vector<T> e(n);
  vector<T> d1(n);
  vector<T> e1(n);

  for (size_t i = 0; i < n; i++) {
    d1[i] = -t[0][i];
    e1[i] = x1[i] - t[1][i];
  }
  auto td0 = chl.recv(d);
  auto te0 = chl.recv(e);
  auto td1 = chl.send(d1);
  auto te1 = chl.send(e1);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
  for (size_t i = 0; i < n; i++) {
    d[i] = d[i] + d1[i];
    e[i] = e[i] + e1[i];
  }
  vector<T> res(n);
  for (size_t i = 0; i < n; i++) {
    res[i] = d[i] * t[1][i] + t[0][i] * e[i] + t[2][i];
  }
  return res;
}

template<typename T>
vector<vector<T>> mul0(cp::Socket chl, array<vector<T>, 3> &t, vector<vector<T>> &x0) {
  u64 m = x0.size();
  u64 n = x0[0].size();
  vector<T> d0(m * n);
  vector<T> e0(m * n);
  vector<T> d(m * n);
  vector<T> e(m * n);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      d0[i * n + j] = x0[i][j] - t[0][i * n + j];
      e0[i * n + j] = -t[1][i * n + j];
    }
  }
  auto td0 = chl.send(d0);
  auto te0 = chl.send(e0);
  auto td1 = chl.recv(d);
  auto te1 = chl.recv(e);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      d[i * n + j] = d0[i * n + j] + d[i * n + j];
      e[i * n + j] = e0[i * n + j] + e[i * n + j];
    }
  }
  vector<vector<T>> res(m, vector<T>(n));
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      res[i][j] = d[i * n + j] * e[i * n + j] + d[i * n + j] * t[1][i * n + j] + t[0][i * n + j] * e[i * n + j]
          + t[2][i * n + j];
    }
  }
  return res;
}

template<typename T>
vector<vector<T>> mul1(cp::Socket chl, array<vector<T>, 3> &t, vector<vector<T>> &x1) {
  u64 m = x1.size();
  u64 n = x1[0].size();
  vector<T> d(m * n);
  vector<T> e(m * n);
  vector<T> d1(m * n);
  vector<T> e1(m * n);

  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      d1[i * n + j] = -t[0][i * n + j];
      e1[i * n + j] = x1[i][j] - t[1][i * n + j];
    }
  }
  auto td0 = chl.recv(d);
  auto te0 = chl.recv(e);
  auto td1 = chl.send(d1);
  auto te1 = chl.send(e1);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      d[i * n + j] = d[i * n + j] + d1[i * n + j];
      e[i * n + j] = e[i * n + j] + e1[i * n + j];
    }
  }
  vector<vector<T>> res(m, vector<T>(n));
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      res[i][j] = d[i * n + j] * t[1][i * n + j] + t[0][i * n + j] * e[i * n + j] + t[2][i * n + j];
    }
  }
  return res;
}

template<typename T>
void verifyMul(u64 n, u64 m) {
  auto triple = fakeTriple<T>(n * m);
  vector<vector<T>> x0(m, vector<T>(n));
  vector<vector<T>> x1(m, vector<T>(n));
  PRNG prng(ZeroBlock);
  for (size_t i = 0; i < m; i++) {
    prng.get(x0[i].data(), x0[i].size());
    prng.get(x1[i].data(), x1[i].size());
  }

  vector<vector<T>> res0, res1;
  thread t0 = thread([&]() {
    auto chl0 = cp::asioConnect("localhost:7700", true);
    res0 = mul0(chl0, triple.first, x0);
  });
  thread t1 = thread([&]() {
    auto chl1 = cp::asioConnect("localhost:7700", false);
    res1 = mul1(chl1, triple.second, x1);
  });
  t0.join();
  t1.join();

  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      T left = res0[i][j] + res1[i][j];
      T right = x0[i][j] * x1[i][j];
      if (left != right) {
        throw RTE_LOC;
      }
    }
  }
}

#endif
