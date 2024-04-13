#include <coproto/Socket/AsioSocket.h>
#include <coproto/Socket/LocalAsyncSock.h>
#include <cryptoTools/Common/CLP.h>
#include <libOTe/Base/BaseOT.h>

#include <array>

#include "fingercodes-vt.hpp"
#include "gmw.hpp"
#include "triple.hpp"
#include "mp.hpp"
#include "vt.hpp"
#include "ed.hpp"
#include "eigenface-vt.hpp"

using namespace osuCrypto;
using namespace std;

std::random_device rd;
std::default_random_engine gen(rd());
std::uniform_int_distribution<u8> r1(0, 1);
u8 randBool() { return r1(gen); }

void testBase64() {
  {
    vector<u8> data(900);
    for (size_t i = 0; i < 900; i++) {
      data[i] = i % 2;
    }
    auto s = boolsToBase64(data);
    // cout << s << endl;
    auto data2 = base64ToBools(s);
    for (size_t i = 0; i < 900; i++) {
      if (data[i] != data2[i]) {
        cout << "base64 error: " << i << " " << data[i] << " " << data2[i]
             << endl;
        return;
      }
    }
  }

  {
    string s =
        "8fDXOJWWy/Jlym+3HqGkFftJqw/AvTquGMx+q/"
        "G1P85tK6D5rug14yJJRhLQ8AaDEYhmkt9I5Tsr5a6ZBMhORFKWeon00000000004cU7Km+"
        "xcRtSKXSC6gj0w17sbgv0ZXPaGBcC/WGmQBevxgp5j==";
    auto data = base64ToBools(s);
    if (data.size() != 900) {
      cout << "size error" << endl;
      return;
    }
    auto ss = boolsToBase64(data);
    if (s != ss) {
      cout << "base64 error" << endl;
      return;
    }
  }

  {
    string s =
        "AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQ"
        "EBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB"
        "AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQ"
        "EBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB"
        "AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQ"
        "EBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB"
        "AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQ"
        "EBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB"
        "AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQ"
        "EBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB"
        "AQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQ"
        "EBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEBAQEB"
        "AQEBAQEBAQEBAQ==";
    auto data = base64_decode(s);
    // if (data.size() != 640) {
    //   cout << "size error: " << data.size() << endl;
    //   return;
    // }
    auto ss = base64_encode(data);
    if (s != ss) {
      cout << "base64 error" << endl;
      return;
    }
  }

  {
    vector<u8> data(640, 4);
    auto s = base64_encode(data);
    // cout << s << endl;
    auto data2 = base64_decode(s);
    // if (data.size() != data2.size()) {
    //   cout << "size error: " << data.size() << " " << data2.size() << endl;
    //   return;
    // }
    for (size_t i = 0; i < data.size(); i++) {
      if (data[i] != data2[i]) {
        cout << "base64 error: " << i << " " << data[i] << " " << data2[i]
             << endl;
        return;
      }
    }
  }
}

//void testABY() {
//  {
//    u64 n = 10000;
//    u64 mod = 900;
//    u64 threshold = 10;
//    u32 rand = rd() % mod;
//    ABY aby(mod, 10, 32);
//    vector<u32> T(n);
//    for (size_t i = 0; i < n; i++) {
//      T[i] = (rand + i * 11 + threshold + mod) % mod;
//    }
//    vector<u32> R(n);
//    for (size_t i = 0; i < n; i++) {
//      R[i] = (rand + i * 11 - threshold + mod) % mod;
//    }
//
//    {
//      thread t0 = thread([&]() {
//        bool res = aby.MatchingT(T);
//        if (res) {
//          cout << "ABY error" << endl;
//          exit(1);
//        }
//      });
//      thread t1 = thread([&]() {
//        bool res = aby.MatchingR(R);
//        if (res) {
//          cout << "ABY error" << endl;
//          exit(1);
//        }
//      });
//      t0.join();
//      t1.join();
//    }
//
//    u64 i = rd() % n;
//    R[i] = (rand + i * 11 + threshold + mod) % mod;
//
//    {
//      thread t0 = thread([&]() {
//        bool res = aby.MatchingT(T);
//        if (!res) {
//          cout << "ABY error" << endl;
//          exit(1);
//        }
//      });
//      thread t1 = thread([&]() {
//        bool res = aby.MatchingR(R);
//        if (!res) {
//          cout << "ABY error" << endl;
//          exit(1);
//        }
//      });
//      t0.join();
//      t1.join();
//    }
//  }
//}

void testGMW() {
  u64 n = 50000;
  u64 mod = 900;
  u64 threshold = 10;
  u32 rand = rd() % mod;
  GMW gmw;
  vector<u32> T(n);
  for (size_t i = 0; i < n; i++) {
    T[i] = (rand + i * 11 + mod);
  }
  vector<u32> R(n);
  for (size_t i = 0; i < n; i++) {
    R[i] = (rand + i * 11 - threshold + mod);
  }

  // auto chl = cp::LocalAsyncSocket::makePair();

  {
    u64 m = n * 16;
    vector<u8> a0(m);
    vector<u8> b0(m);
    vector<u8> c0(m);
    vector<u8> a1(m);
    vector<u8> b1(m);
    vector<u8> c1(m);

    for (u64 i = 0; i < m; i++) {
      a0[i] = randBool();
      b0[i] = randBool();
      c0[i] = randBool();
      a1[i] = randBool();
      b1[i] = randBool();
      c1[i] = c0[i] ^ ((a0[i] ^ a1[i]) & (b0[i] ^ b1[i]));
    }

    vector<u8> carry0;
    vector<u8> carry1;

    thread t0 = thread([&]() {
      auto c = cp::asioConnect("127.0.0.1:7700", true);
      carry0 = gmw.MatchingT(c, T, {a0, b0, c0});
    });
    thread t1 = thread([&]() {
      auto c = cp::asioConnect("127.0.0.1:7700", false);
      carry1 = gmw.MatchingR(c, R, {a1, b1, c1});
    });
    t0.join();
    t1.join();

    for (u64 i = 0; i < 100; i++) {
      if ((carry0[i] ^ carry1[i])) {
        cout << "GMW error: " << i << endl;
        return;
      }
    }
  }
}

void verify() {
  testBase64();
//  testABY();
  testGMW();
}

template<u64 m>
void noisyBenchHelper(u64 n) {
  std::cout << "*** " << m << "x" << n << " ***" << std::endl;
  noisyBench<u32, m>(1, n);
  noisyBench<u32, m>(8, n);
  noisyBench<u32, m>(64, n);
  std::cout << "\n\n" << std::endl;
}

void vtBenchNoisy() {
  using namespace oc::Subfield;

  std::cout << "*** noisy ***\n\n\n\n\n" << std::endl;
//  noisyBenchHelper<8>(8);
//  noisyBenchHelper<8>(32);
//  noisyBenchHelper<8>(256);
//  noisyBenchHelper<8>(1024);
//  noisyBenchHelper<8>(4096);
//  noisyBenchHelper<8>(16384);
//  noisyBenchHelper<32>(32);
//  noisyBenchHelper<32>(256);
//  noisyBenchHelper<32>(1024);
//  noisyBenchHelper<32>(4096);
//  noisyBenchHelper<32>(16384);
//  noisyBenchHelper<256>(256);
//  noisyBenchHelper<256>(1024);
//  noisyBenchHelper<256>(4096);
//  noisyBenchHelper<256>(16384);
//  noisyBenchHelper<1024>(1024);
//  noisyBenchHelper<1024>(4096);
//  noisyBenchHelper<1024>(16384);
}

template<u64 m>
void silentBenchHelper(u64 n) {
  std::cout << "*** " << m << "x" << n << " ***" << std::endl;
  silentBench<u32, m>(1, n);
  silentBench<u32, m>(32, n);
  silentBench<u32, m>(64, n);
  std::cout << "\n\n" << std::endl;
}

void vtBenchSilent() {
  std::cout << "*** silent ***\n\n\n\n\n" << std::endl;

//  silentBenchHelper<8>(8);
//  silentBenchHelper<8>(32);
//  silentBenchHelper<8>(256);
//  silentBenchHelper<8>(1024);
//  silentBenchHelper<8>(4096);
//  silentBenchHelper<8>(16384);
//  silentBenchHelper<32>(32);
//  silentBenchHelper<32>(256);
//  silentBenchHelper<32>(1024);
//  silentBenchHelper<32>(4096);
//  silentBenchHelper<32>(16384);
//  silentBenchHelper<256>(256);
//  silentBenchHelper<256>(1024);
//  silentBenchHelper<256>(4096);
//  silentBenchHelper<256>(16384);
//  silentBenchHelper<1024>(1024);
//  silentBenchHelper<1024>(4096);
//  silentBenchHelper<1024>(16384);
}

void forkTest() {
  u64 s = 123;
  u64 r;

  thread t0 = thread([&]() {
    auto chl0 = cp::asioConnect("localhost:7700", true);
    cp::sync_wait(chl0.send(s));
    auto f0 = chl0.fork();
    s = 312;
    cp::sync_wait(f0.send(s));
    auto f1 = f0.fork();
    s = 321;
    cp::sync_wait(f1.send(s));
  });

  thread t1 = thread([&]() {
    auto chl1 = cp::asioConnect("localhost:7700", false);
    cp::sync_wait(chl1.recv(r));
    auto f0 = chl1.fork();
    cp::sync_wait(f0.recv(r));
    auto f1 = f0.fork();
    cp::sync_wait(f1.recv(r));
  });

  t0.join();
  t1.join();
  assert(s == r);
  cout << "OK" << endl;
}

void forkTest2() {
  u64 n = 100;
  u64 m = 100;
  bool fork = true;

  thread t0 = thread([&]() {
    auto chl = cp::asioConnect("localhost:7700", true);
    vector<cp::Socket> forks(m);
    for (u64 i = 0; i < m; i++) {
      if (fork) {
        forks[i] = chl.fork();
      } else {
        forks[i] = cp::asioConnect("localhost:" + to_string(7800 + i), true);
      }
    }

    vector<SilentVoleSender> send(m);
    vector<PRNG> prngs(m);
    vector<block> d(m);
    vector<vector<block>> z0(m, vector<block>(n));
    vector<macoro::eager_task<>> tasks;
    tasks.reserve(m);

    for (u64 i = 0; i < m; i++) {
      prngs[i].SetSeed(ZeroBlock);
      auto task = send[i].silentSend(d[i], z0[i], prngs[i], forks[i]);
      tasks.emplace_back(std::move(task) | macoro::make_eager());
    }

    for (u64 i = 0; i < m; i++) {
      cp::sync_wait(tasks[i]);
    }

    cp::sync_wait(chl.flush());
    if (!fork) {
      for (u64 i = 0; i < m; i++) {
        cp::sync_wait(forks[i].flush());
      }
    }
  });

  thread t1 = thread([&]() {
    auto chl = cp::asioConnect("localhost:7700", false);
    vector<cp::Socket> forks(m);
    for (u64 i = 0; i < m; i++) {
      if (fork) {
        forks[i] = chl.fork();
      } else {
        forks[i] = cp::asioConnect("localhost:" + to_string(7800 + i), false);
      }
    }

    vector<SilentVoleReceiver> recv(m);
    vector<PRNG> prngs(m);
    vector<vector<block>> y(m, vector<block>(n));
    vector<vector<block>> z1(m, vector<block>(n));
    vector<macoro::eager_task<>> tasks;
    tasks.reserve(m);

    for (u64 i = 0; i < m; i++) {
      prngs[i].SetSeed(ZeroBlock);
      auto task = recv[i].silentReceive(y[i], z1[i], prngs[i], forks[i]);
      tasks.emplace_back(std::move(task) | macoro::make_eager());
    }

    for (u64 i = 0; i < m; i++) {
      cp::sync_wait(tasks[i]);
    }

    cp::sync_wait(chl.flush());
    if (!fork) {
      for (u64 i = 0; i < m; i++) {
        cp::sync_wait(forks[i].flush());
      }
    }
  });

  t0.join();
  t1.join();
  cout << "OK" << endl;
}


void test() {

//  {
//    thread t0 = thread([&]() {
//      FingerCodesServerVT("fingercodes-128.txt");
////      FingerCodesServerVT("fingercodes-100.txt");
//    });
//    thread t1 = thread([&]() {
//      FingerCodesClientVT("fingercodes-128.txt");
////      FingerCodesClientVT("fingercodes-1.txt");
//    });
//
//    t0.join();
//    t1.join();
//  }

//  forkTest2();

//  vtBenchNoisy();
//  vtBenchSilent();

  noisyBench<u32, 64>(1, 64);
  silentBench<u32, 64>(1, 64);
}
