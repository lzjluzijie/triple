#ifndef MP_HPP
#define MP_HPP

#include <coproto/Socket/AsioSocket.h>
#include <cryptoTools/Common/CLP.h>
#include <coproto/Socket/Socket.h>

#include <span>
#include <utility>
#include <cryptoTools/Crypto/PRNG.h>
#include <coproto/Socket/LocalAsyncSock.h>

using namespace osuCrypto;
using namespace std;
namespace cp = coproto;

template<typename T>
using slice = nonstd::span<T>;

template<typename T>
using Triples = vector<array<vector<T>, 3>>;

template<typename T>
inline void op(const vector<T> &a, const vector<T> &b, vector<vector<T>> &res) {
  size_t m = a.size();
  size_t n = b.size();
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      res[i][j] = a[i] * b[j];
    }
  }
}

template<typename T>
vector<vector<T>> mp(const vector<vector<T>> &X, const vector<vector<T>> &Y) {
  size_t a = X.size();
  size_t b = X[0].size();
  size_t c = Y[0].size();
  assert(b == Y.size());
  vector<vector<T>> res(a, vector<T>(c, 0));
  for (size_t i = 0; i < a; i++) {
    for (size_t j = 0; j < c; j++) {
      for (size_t k = 0; k < b; k++) {
        res[i][j] += X[i][k] * Y[k][j];
      }
    }
  }
  return res;
}

template<typename T, typename RES = T>
vector<vector<RES>> mp0(cp::Socket chl, vector<array<vector<T>, 3>> &t,
                        const vector<vector<T>> &x0, const vector<vector<T>> &y0) {
  u64 a = x0.size();
  u64 b = x0[0].size();
  u64 c = y0[0].size();
  assert(b == y0.size());
//  T *d0 = new T[a * b];
//  T *e0 = new T[a * c];
//  T *d = new T[a * b];
//  T *e = new T[a * c];
//  for (u64 k = 0; k < b; k++) {
//    for (u64 i = 0; i < a; i++) {
//      d[k * b + i] = x0[i][k] - t[k][0][i];
//    }
//    for (u64 i = 0; i < c; i++) {
//      e[k * c + i] = y0[k][i] - t[k][1][i];
//    }
//  }
//  auto td1 = chl.recv(slice<T>(d0, a * b));
//  auto te1 = chl.recv(slice<T>(e0, a * c));
//  auto td0 = chl.send(slice<T>(d, a * b));
//  auto te0 = chl.send(slice<T>(e, a * c));
  vector<T> d0(b * a);
  vector<T> e0(b * c);
  vector<T> d(b * a);
  vector<T> e(b * c);
  for (u64 k = 0; k < b; k++) {
    for (u64 i = 0; i < a; i++) {
      d[k * a + i] = x0[i][k] - t[k][0][i];
    }
    for (u64 i = 0; i < c; i++) {
      e[k * c + i] = y0[k][i] - t[k][1][i];
    }
  }
  auto td1 = chl.recv(d0);
  auto te1 = chl.recv(e0);
  auto td0 = chl.send(d);
  auto te0 = chl.send(e);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
//  cp::sync_wait(chl.recv(d0));
//  cp::sync_wait(chl.recv(e0));
//  cp::sync_wait(chl.send(d));
//  cp::sync_wait(chl.send(e));
  for (u64 i = 0; i < b * a; i++) {
    d[i] += d0[i];
  }
  for (u64 i = 0; i < b * c; i++) {
    e[i] += e0[i];
  }
  vector<vector<RES>> res(a, vector<RES>(c, 0));
  for (u64 k = 0; k < b; k++) {
    for (size_t i = 0; i < a; i++) {
      for (size_t j = 0; j < c; j++) {
        res[i][j] +=
            d[k * a + i] * e[k * c + j] + d[k * a + i] * t[k][1][j] + t[k][0][i] * e[k * c + j] + t[k][2][i * c + j];
      }
    }
  }
  return res;
}

template<typename T, typename RES = T>
vector<vector<RES>> mp1(cp::Socket chl, vector<array<vector<T>, 3>> &t,
                        const vector<vector<T>> &x1, const vector<vector<T>> &y1) {
  u64 a = x1.size();
  u64 b = x1[0].size();
  u64 c = y1[0].size();
  assert(b == y1.size());
//  T *d = new T[a * b];
//  T *e = new T[a * c];
//  T *d1 = new T[a * b];
//  T *e1 = new T[a * c];
//  for (u64 k = 0; k < b; k++) {
//    for (u64 i = 0; i < a; i++) {
//      d[k * b + i] = x1[i][k] - t[k][0][i];
//    }
//    for (u64 i = 0; i < c; i++) {
//      e[k * c + i] = y1[k][i] - t[k][1][i];
//    }
//  }
//  auto td0 = chl.send(slice<T>(d, a * b));
//  auto te0 = chl.send(slice<T>(e, a * c));
//  auto td1 = chl.recv(slice<T>(d1, a * b));
//  auto te1 = chl.recv(slice<T>(e1, a * c));
  vector<T> d(b * a);
  vector<T> e(b * c);
  vector<T> d1(b * a);
  vector<T> e1(b * c);
  for (u64 k = 0; k < b; k++) {
    for (u64 i = 0; i < a; i++) {
      d[k * a + i] = x1[i][k] - t[k][0][i];
    }
    for (u64 i = 0; i < c; i++) {
      e[k * c + i] = y1[k][i] - t[k][1][i];
    }
  }
  auto td0 = chl.send(d);
  auto te0 = chl.send(e);
  auto td1 = chl.recv(d1);
  auto te1 = chl.recv(e1);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
//  cp::sync_wait(chl.send(d));
//  cp::sync_wait(chl.send(e));
//  cp::sync_wait(chl.recv(d1));
//  cp::sync_wait(chl.recv(e1));
  for (u64 i = 0; i < b * a; i++) {
    d[i] += d1[i];
  }
  for (u64 i = 0; i < b * c; i++) {
    e[i] += e1[i];
  }
  vector<vector<RES>> res(a, vector<RES>(c, 0));
  for (u64 k = 0; k < b; k++) {
    for (size_t i = 0; i < a; i++) {
      for (size_t j = 0; j < c; j++) {
        res[i][j] += d[k * a + i] * t[k][1][j] + t[k][0][i] * e[k * c + j] + t[k][2][i * c + j];
      }
    }
  }
  return res;
}

template<typename T>
pair<Triples<T>, Triples<T>> fakeZ(u64 m, u64 n, u64 num, block seed = ZeroBlock) {
  PRNG prng(seed);
  Triples<T> r0(num);
  Triples<T> r1(num);
  for (u64 k = 0; k < num; ++k) {
    array<vector<T>, 3> t0;
    array<vector<T>, 3> t1;
    t0[0].resize(m);
    prng.get(t0[0].data(), t0[0].size());
    t0[1].resize(n);
    prng.get(t0[1].data(), t0[1].size());
    t0[2].resize(m * n);
    prng.get(t0[2].data(), t0[2].size());
    t1[0].resize(m);
    prng.get(t1[0].data(), t1[0].size());
    t1[1].resize(n);
    prng.get(t1[1].data(), t1[1].size());
    t1[2].resize(m * n);
    for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
        t1[2][i * n + j] =
            (t0[0][i] + t1[0][i]) * (t0[1][j] + t1[1][j]) - t0[2][i * n + j];
      }
    }
    r0[k] = t0;
    r1[k] = t1;
  }
  return pair<Triples<T>, Triples<T>>(r0, r1);
}

template<typename T>
void verifyZ(Triples<T> t0, Triples<T> t1) {
  u64 nTriples = t0.size();
  u64 a = t0[0][0].size();
  u64 c = t0[0][1].size();
  for (u64 k = 0; k < nTriples; ++k) {
    auto a0 = t0[k][0];
    auto b0 = t0[k][1];
    auto c0 = t0[k][2];
    auto a1 = t1[k][0];
    auto b1 = t1[k][1];
    auto c1 = t1[k][2];
    for (u64 i = 0; i < a; ++i) {
      for (u64 j = 0; j < c; ++j) {
        T left = (a0[i] + a1[i]) * (b0[j] + b1[j]);
        T right = c0[i * c + j] + c1[i * c + j];
        if (left != right) {
          cout << "failed " << i << " " << j << endl;
          throw RTE_LOC;
        }
      }
    }
  }
}

template<typename T>
void mpTest(u64 a, u64 b, u64 c, u64 num) {
  PRNG prng(sysRandomSeed());
  vector<vector<T>> x0(a, vector<T>(b));
  vector<vector<T>> x1(a, vector<T>(b));
  vector<vector<T>> y0(b, vector<T>(c));
  vector<vector<T>> y1(b, vector<T>(c));

//  vector<vector<T>> x(a, vector<T>(b));
//  vector<vector<T>> y(b, vector<T>(c));
//
//  for (size_t i = 0; i < a; i++) {
//    prng.get(x0[i].data(), x0[i].size());
//    prng.get(x1[i].data(), x1[i].size());
//    for (size_t j = 0; j < b; j++) {
//      x[i][j] = x0[i][j] + x1[i][j];
//    }
//  }
//
//  for (size_t i = 0; i < b; i++) {
//    prng.get(y0[i].data(), y0[i].size());
//    prng.get(y1[i].data(), y1[i].size());
//    for (size_t j = 0; j < c; j++) {
//      y[i][j] = y0[i][j] + y1[i][j];
//    }
//  }
//
//  vector<vector<T>> p = mp(x, y);
//  vector<vector<vector<T>>> p0(num);
//  vector<vector<vector<T>>> p1(num);

  auto chl = cp::LocalAsyncSocket::makePair();

  vector<array<vector<T>, 3>> triples0;
  vector<array<vector<T>, 3>> triples1;

  {
    auto t = fakeZ<T>(a, c, b);
    triples0 = t.first;
    triples1 = t.second;
  }

  verifyZ(triples0, triples1);

  size_t t0sent = 0;
  size_t t0recv = 0;
  size_t t1sent = 0;
  size_t t1recv = 0;

  vector<u64> times(num);

  {
    thread t0 = thread([&]() {
      auto chl0 = cp::asioConnect("127.0.0.1:7700", true);
      for (size_t i = 0; i < num; i++) {
//        p0[i] = mp0(chl0, triples0, x0, y0);
        auto start = chrono::steady_clock::now();
        mp0(chl0, triples0, x0, y0);
        auto end = chrono::steady_clock::now();
        auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
        times[i] = time;
      }
      t0sent = chl0.bytesSent();
      t0recv = chl0.bytesReceived();
    });
    thread t1 = thread([&]() {
      auto chl1 = cp::asioConnect("127.0.0.1:7700", false);
      vector<T> tp1(num);
      for (size_t i = 0; i < num; i++) {
//        p1[i] = mp1(chl1, triples1, x1, y1);
        mp1(chl1, triples1, x1, y1);
      }
      t1sent = chl1.bytesSent();
      t1recv = chl1.bytesReceived();
    });
    t0.join();
    t1.join();
    u64 traffic = t0sent + t0recv;
    traffic /= num;
    cout << "computing: " << num << " matrix products of size " << a << "x" << b << " " << b << "x" << c;
//    cout << " t0 sent " << t0sent << " recv " << t0recv;
//    cout << " t1 sent " << t1sent << " recv " << t1recv;
    cout << " traffic " << traffic << "   " << ((traffic) >> 10) << "KiB ";
    cout << (traffic >> 20) << "MiB" << endl;
//    for (u64 i = 0; i < num; i++) {
//      cout << times[i] << " ";
//    }
    sort(times.begin(), times.end());
    cout << "median " << times[num / 2] << endl;

//    for (u64 k = 0; k < num; k++) {
//      for (u64 i = 0; i < a; i++) {
//        for (u64 j = 0; j < c; j++) {
//          T v = p0[k][i][j] + p1[k][i][j];
//          if (p[i][j] != v) {
//            cout << "error at " << i << " " << j << endl;
//            exit(1);
//          }
//        }
//      }
//    }


  }
}

#endif
