#ifndef ED_HPP
#define ED_HPP

#include "mp.hpp"
#include <cryptoTools/Common/Timer.h>

template<typename T>
void edTest(u64 m, u64 n, u64 k) {
  Timer timer;
  timer.setTimePoint("start");
  PRNG prng(ZeroBlock);
  vector<vector<T>> x0(m, vector<T>(k));
  vector<vector<T>> x1(m, vector<T>(k));
  vector<vector<T>> y0(k, vector<T>(n));
  vector<vector<T>> y1(k, vector<T>(n));
  for (u64 i = 0; i < m; i++) {
    prng.get(x0[i].data(), k);
//    prng.get(x1[i].data(), k);
  }
  for (u64 i = 0; i < k; i++) {
//    prng.get(y0[i].data(), n);
    prng.get(y1[i].data(), n);
  }
  vector<vector<T>> ed(m, vector<T>(n, 0));
  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < n; j++) {
      for (u64 l = 0; l < k; l++) {
        T t = (x0[i][l] + x1[i][l]) - (y0[l][j] + y1[l][j]);
        ed[i][j] += t * t;
      }
    }
  }
  auto triples = fakeZ<T>(m, n, k);

  vector<vector<T>> ed0(m, vector<T>(n));
  vector<vector<T>> ed1(m, vector<T>(n));
  timer.setTimePoint("init");

  thread t0 = thread([&]() {
    vector<T> sqx0(m, 0);
    vector<T> sqy0(n, 0);
    for (u64 i = 0; i < m; i++) {
      for (u64 l = 0; l < k; l++) {
        sqx0[i] += x0[i][l] * x0[i][l];
      }
    }
    for (u64 i = 0; i < k; i++) {
      for (u64 j = 0; j < n; j++) {
        sqy0[j] += y0[i][j] * y0[i][j];
      }
    }
    timer.setTimePoint("t0 local");

    auto chl = cp::asioConnect("localhost:7700", true);
    auto res = mp0<T>(chl, triples.first, x0, y0);
    timer.setTimePoint("mp0 done");

    for (u64 i = 0; i < m; i++) {
      for (u64 j = 0; j < n; j++) {
        ed0[i][j] = sqx0[i] + sqy0[j] - 2 * res[i][j];
      }
    }
    timer.setTimePoint("t0 done");
  });

  thread t1 = thread([&]() {
    vector<T> sqx1(m, 0);
    vector<T> sqy1(n, 0);
    for (u64 i = 0; i < m; i++) {
      for (u64 l = 0; l < k; l++) {
        sqx1[i] += x1[i][l] * x1[i][l];
      }
    }
    for (u64 i = 0; i < k; i++) {
      for (u64 j = 0; j < n; j++) {
        sqy1[j] += y1[i][j] * y1[i][j];
      }
    }
    timer.setTimePoint("t1 local");

    auto chl = cp::asioConnect("localhost:7700", false);
    auto res = mp1<T>(chl, triples.second, x1, y1);
    timer.setTimePoint("mp1 done");

    for (u64 i = 0; i < m; i++) {
      for (u64 j = 0; j < n; j++) {
        ed1[i][j] = sqx1[i] + sqy1[j] - 2 * res[i][j];
      }
    }
    timer.setTimePoint("t1 done");
  });

  t0.join();
  t1.join();
  timer.setTimePoint("done");

  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < n; j++) {
      if (ed[i][j] != T(ed0[i][j] + ed1[i][j])) {
        cout << "ed[" << i << "][" << j << "] = " << ed[i][j] << endl;
        cout << "ed0[" << i << "][" << j << "] = " << ed0[i][j] << endl;
        cout << "ed1[" << i << "][" << j << "] = " << ed1[i][j] << endl;
        throw runtime_error("ed0 + ed1 != ed");
      }
    }
  }
  timer.setTimePoint("verify");
  cout << timer << endl;
}

template<typename T>
void edTest1(u64 m, u64 n, u64 k, u64 iterations) {
  Timer timer;
  timer.setTimePoint("start");
//  PRNG prng(ZeroBlock);
  PRNG prng(sysRandomSeed());
  vector<vector<T>> x0(m, vector<T>(k));
  vector<vector<T>> x1(m, vector<T>(k));
  vector<vector<T>> y0(k, vector<T>(n));
  vector<vector<T>> y1(k, vector<T>(n));
  for (u64 i = 0; i < m; i++) {
    prng.get(x0[i].data(), k);
//    prng.get(x1[i].data(), k);
  }
  for (u64 i = 0; i < k; i++) {
//    prng.get(y0[i].data(), n);
    prng.get(y1[i].data(), n);
  }
  vector<vector<T>> ed(m, vector<T>(n, 0));
  for (u64 i = 0; i < m; i++) {
    for (u64 j = 0; j < n; j++) {
      for (u64 l = 0; l < k; l++) {
        T t = (x0[i][l] + x1[i][l]) - (y0[l][j] + y1[l][j]);
        ed[i][j] += t * t;
      }
    }
  }
  auto triples = fakeZ<T>(m, n, k);

  Matrix<T> ed0(m, n);
  Matrix<T> ed1(m, n);
  Matrix<T> temp(m, n);
  timer.setTimePoint("init");

  vector<u64> times(iterations);
  u64 traffic = 0;

  for (u64 it = 0; it < iterations; it++) {
    auto start = chrono::high_resolution_clock::now();

    thread t0 = thread([&]() {
      vector<T> sqx0(m, 0);
      vector<T> sqy0(n, 0);
      for (u64 i = 0; i < m; i++) {
        for (u64 l = 0; l < k; l++) {
          sqx0[i] += x0[i][l] * x0[i][l];
        }
      }
      for (u64 i = 0; i < k; i++) {
        for (u64 j = 0; j < n; j++) {
          sqy0[j] += y0[i][j] * y0[i][j];
        }
      }
      timer.setTimePoint("t0 local");

      auto chl = cp::asioConnect("localhost:7700", true);
      auto res = mp0<T>(chl, triples.first, x0, y0);
      timer.setTimePoint("mp0 done");

      for (u64 i = 0; i < m; i++) {
        for (u64 j = 0; j < n; j++) {
          ed0[i][j] = sqx0[i] + sqy0[j] - 2 * res[i][j];
        }
      }

      cp::sync_wait(chl.send(ed0));

      traffic += chl.bytesReceived() + chl.bytesSent();
      timer.setTimePoint("t0 done");
    });

    thread t1 = thread([&]() {
      vector<T> sqx1(m, 0);
      vector<T> sqy1(n, 0);
      for (u64 i = 0; i < m; i++) {
        for (u64 l = 0; l < k; l++) {
          sqx1[i] += x1[i][l] * x1[i][l];
        }
      }
      for (u64 i = 0; i < k; i++) {
        for (u64 j = 0; j < n; j++) {
          sqy1[j] += y1[i][j] * y1[i][j];
        }
      }
      timer.setTimePoint("t1 local");

      auto chl = cp::asioConnect("localhost:7700", false);
      auto res = mp1<T>(chl, triples.second, x1, y1);
      timer.setTimePoint("mp1 done");

      for (u64 i = 0; i < m; i++) {
        for (u64 j = 0; j < n; j++) {
          ed1[i][j] = sqx1[i] + sqy1[j] - 2 * res[i][j];
        }
      }

      cp::sync_wait(chl.recv(temp));
      for (u64 i = 0; i < m; i++) {
        for (u64 j = 0; j < n; j++) {
          ed1[i][j] += temp[i][j];
        }
      }

      timer.setTimePoint("t1 done");
    });

    t0.join();
    t1.join();

    auto end = chrono::high_resolution_clock::now();
    times[it] = chrono::duration_cast<chrono::microseconds>(end - start).count();
    timer.setTimePoint("done");

    for (u64 i = 0; i < m; i++) {
      for (u64 j = 0; j < n; j++) {
//        if (ed[i][j] != T(ed0[i][j] + ed1[i][j])) {
        if (ed[i][j] != ed1[i][j]) {
          cout << "ed[" << i << "][" << j << "] = " << ed[i][j] << endl;
          cout << "ed0[" << i << "][" << j << "] = " << ed0[i][j] << endl;
          cout << "ed1[" << i << "][" << j << "] = " << ed1[i][j] << endl;
//          throw runtime_error("ed0 + ed1 != ed");
        }
      }
    }
    timer.setTimePoint("verify");
    //  cout << timer << endl;
  }

  traffic /= iterations;
  cout << "computing: " << iterations << " sqED of size " << m << "x" << k << " " << k << "x" << n;
//    cout << " t0 sent " << t0sent << " recv " << t0recv;
//    cout << " t1 sent " << t1sent << " recv " << t1recv;
  cout << " traffic " << traffic << "   " << ((traffic) >> 10) << "KiB ";
  cout << (traffic >> 20) << "MiB" << endl;
//    for (u64 i = 0; i < num; i++) {
//      cout << times[i] << " ";
//    }
  sort(times.begin(), times.end());
  cout << "median " << times[iterations / 2] << endl;
}

#endif
