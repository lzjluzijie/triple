#ifndef SP_HPP
#define SP_HPP

#include "vt.hpp"

using namespace std;

// scalar product
template<typename T>
inline T sp(const vector<T> &a, const vector<T> &b) {
  // if (a.size() != b.size()) {
  //   throw std::runtime_error("a.size() != b.size()");
  // }
  T res = 0;
  for (size_t i = 0; i < a.size(); i++) {
    res = res + (a[i] * b[i]);
  }
  return res;
}

// trace
template<typename T>
inline T tr(const vector<T> &a, u64 n) {
  // if (a.size() != n * n) {
  //   throw std::runtime_error("a.size() != n * n");
  // }
  T res = 0;
  for (size_t i = 0; i < n; i++) {
    res = res + a[i * n + i];
  }
  return res;
}

template<typename T>
T sp0(cp::Socket chl, const array<vector<T>, 3> &t, const vector<T> &x0,
      const vector<T> &y0) {
  u64 n = x0.size();
  vector<T> d(n);
  vector<T> e(n);
  vector<T> d1(n);
  vector<T> e1(n);
  for (u64 i = 0; i < n; i++) {
    d[i] = x0[i] - t[0][i];
    e[i] = y0[i] - t[1][i];
  }
  auto td0 = chl.send(d);
  auto te0 = chl.send(e);
  auto td1 = chl.recv(d1);
  auto te1 = chl.recv(e1);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
  for (u64 i = 0; i < n; i++) {
    d[i] += d1[i];
    e[i] += e1[i];
  }
  T res = sp(d, e) + sp(d, t[1]) + sp(t[0], e) + tr(t[2], n);
  return res;
}

template<typename T>
T sp1(cp::Socket chl, const array<vector<T>, 3> &t, const vector<T> &x1,
      const vector<T> &y1) {
  u64 n = x1.size();
  vector<T> d0(n);
  vector<T> e0(n);
  vector<T> d(n);
  vector<T> e(n);
  for (size_t i = 0; i < n; i++) {
    d[i] = x1[i] - t[0][i];
    e[i] = y1[i] - t[1][i];
  }
  auto td0 = chl.recv(d0);
  auto te0 = chl.recv(e0);
  auto td1 = chl.send(d);
  auto te1 = chl.send(e);
  cp::sync_wait(cp::when_all_ready(std::move(td0), std::move(te0), std::move(td1), std::move(te1)));
  for (size_t i = 0; i < n; i++) {
    d[i] += d0[i];
    e[i] += e0[i];
  }
  T res = sp(d, t[1]) + sp(t[0], e) + tr(t[2], n);
  return res;
}

template<typename T>
void spTest(u64 n, u64 num) {
  PRNG prng(sysRandomSeed());
  vector<T> x0(n);
  vector<T> x1(n);
  vector<T> x(n);
  vector<T> y0(n);
  vector<T> y1(n);
  vector<T> y(n);

  prng.get(x0.data(), x0.size());
  prng.get(x1.data(), x1.size());
  for (size_t i = 0; i < n; i++) {
    x[i] = x0[i] + x1[i];
  }
  prng.get(y0.data(), y0.size());
  prng.get(y1.data(), y1.size());
  for (size_t i = 0; i < n; i++) {
    y[i] = y0[i] + y1[i];
  }
  T p = sp(x, y);
  // print(p, "sp");

  auto chl = cp::LocalAsyncSocket::makePair();

  vector<array<vector<T>, 3>> triples0;
  vector<array<vector<T>, 3>> triples1;

  // {
  //   auto start = chrono::steady_clock::now();
  //   thread t0 = thread([&]() { triples0 = vtSenderZ(n, n, num, chl[0]); });
  //   thread t1 = thread([&]() { triples1 = vtReceiverZ(n, n, num, chl[1]); });
  //   t0.join();
  //   t1.join();
  //   auto end = chrono::steady_clock::now();
  //   auto ms = chrono::duration_cast<chrono::milliseconds>(end -
  //   start).count(); cout << "generating: " << num << " triples of size " << n
  //   << "x" << n
  //        << " takes " << ms << "ms" << endl;
  // }

  {
    auto t = fakeZ<T>(n, n, num);
    triples0 = t.first;
    triples1 = t.second;
  }

  // verifyZ(triples0, triples1);

  vector<T> p0(num);
  vector<T> p1(num);

  {
    auto start = chrono::steady_clock::now();
    thread t0 = thread([&]() {
      for (size_t i = 0; i < num; i++) {
        p0[i] = sp0(chl[0], triples0[i], x0, y0);
      }
    });
    thread t1 = thread([&]() {
      vector<T> tp1(num);
      for (size_t i = 0; i < num; i++) {
        p1[i] = sp1(chl[1], triples1[i], x1, y1);
      }
    });
    t0.join();
    t1.join();
    auto end = chrono::steady_clock::now();
    auto ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "computing: " << num << " scalar products of size " << n << "x" << n
         << " takes " << ms << "ms" << endl;

    for (u64 i = 0; i < num; i++) {
      if (p0[i] + p1[i] != p) {
        cout << "error at " << i << endl;
        print(p0[i], "p0");
        print(p1[i], "p1");
        print(p, "p");
        exit(1);
      }
    }
  }
}

void spBench();

#endif