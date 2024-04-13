#include "eigenface-vt.hpp"

#include <coproto/Socket/AsioSocket.h>
#include <cryptoTools/Common/CLP.h>
#include <libOTe/Base/BaseOT.h>
#include <libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h>

#include "triple.hpp"
#include "vt.hpp"

using namespace osuCrypto;

const string EigenFaceVT::defaultDB = "eigenface.txt";
uint32_t EigenFaceVT::bitlen = EigenFaceVT::BITLEN;

void EigenFaceServerVT::gen(size_t num) {
  ofstream out(defaultDB);
  out << num << endl;
  for (size_t i = 0; i < num; i++) {
    auto data = randomUint8(edDim);
    auto s = base64_encode(data);
    out << s << endl;
  }
}

void EigenFaceClientVT::gen(size_t num) {
  ofstream out(defaultDB);
  out << num << endl;
  for (size_t i = 0; i < num; i++) {
    auto data = randomUint8(spDim);
    auto s = base64_encode(data);
    out << s << endl;
  }
}

size_t EigenFaceVT::euclidianDistance(vector<uint8_t> x, vector<uint8_t> y) {
  size_t ed = 0;
  for (size_t i = 0; i < edDim; i++) {
    ed += (x[i] - y[i]) * (x[i] - y[i]);
  }
  return ed;
}

inline u64 EigenFaceVT::square(vector<uint8_t> x) {
  u64 sum = 0;
  for (size_t i = 0; i < edDim; i++) {
    sum += x[i] * x[i];
  }
  return sum;
}

EigenFaceServerVT::EigenFaceServerVT(string db) {
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
  } else if (db == "test1") {
    vector<u8> t0(edDim, 0);
    vector<u8> t1(edDim);
    for (size_t i = 0; i < edDim; i++) {
      t1[i] = i;
    }
    data.push_back(t0);
    data.push_back(t1);
    n = 2;
  } else {
    for (size_t j = 0; j < n; j++) {
      data.push_back(randomUint8(edDim));
    }
  }

  // generate eigenfaces
  vector<vector<u8>> eigenfaces(edDim, vector<u8>(spDim));
//  for (size_t i = 0; i < edDim; i++) {
//    prng.get(eigenfaces[i].data(), eigenfaces[i].size());
//  }
  for (size_t i = 0; i < edDim; i++) {
    for (size_t j = 0; j < spDim; j++) {
      eigenfaces[i][j] = (i == j);
    }
  }

  // generate average face
  vector<u8> averageFace(spDim);
//  prng.get(averageFace.data(), averageFace.size());
  for (size_t i = 0; i < spDim; i++) {
    averageFace[i] = i;
  }

  vector<u64> averageFaceNeg(spDim);
  for (size_t i = 0; i < spDim; i++) {
    averageFaceNeg[i] = -averageFace[i];
  }

  auto chl = cp::asioConnect("127.0.0.1:7700", true);
  cp::sync_wait(chl.send(n));
  cp::sync_wait(chl.recv(nc));

  vector<vector<u64>> spx0(nc, averageFaceNeg);
  vector<vector<u64>> spy0(spDim, vector<u64>(edDim));
  for (size_t k = 0; k < spDim; k++) {
    for (size_t i = 0; i < edDim; i++) {
      spy0[k][i] = eigenfaces[i][k];
    }
  }

  vector<vector<u64>> edy0(edDim, vector<u64>(n));
  for (size_t i = 0; i < edDim; i++) {
    for (size_t j = 0; j < n; j++) {
      edy0[i][j] = data[j][i];
    }
  }

  auto spTriple = fakeZ<u64>(nc, edDim, spDim);

  auto spStart = chrono::high_resolution_clock::now();
  auto spShare = mp0<u64, u64>(chl, spTriple.first, spx0, spy0);
  auto spEnd = chrono::high_resolution_clock::now();
  auto spTime = chrono::duration_cast<chrono::microseconds>(spEnd - spStart).count();
  u64 spTraffic = chl.bytesSent() + chl.bytesReceived();
  cout << "SP time: " << spTime
       << "us with transferred " << spTraffic << " "
       << (spTraffic >> 10) << "KiB" << endl;

  auto edTriple = fakeZ<u64>(nc, n, edDim);
  vector<vector<u64>> edx0(nc, vector<u64>(edDim));
  for (size_t i = 0; i < nc; i++) {
    for (size_t j = 0; j < edDim; j++) {
      edx0[i][j] = spShare[i][j];
    }
  }

  auto edStart = chrono::high_resolution_clock::now();
  auto edShare = mp0<u64>(chl, edTriple.first, edx0, edy0);
  auto edEnd = chrono::high_resolution_clock::now();
  auto edTime = chrono::duration_cast<chrono::microseconds>(edEnd - edStart).count();
  u64 edTraffic = chl.bytesSent() + chl.bytesReceived() - spTraffic;
  cout << "ED time: " << edTime
       << "us with transferred " << edTraffic << " "
       << (edTraffic >> 10) << "KiB" << endl;

  auto mulTriple = fakeTriple<u64>(nc * edDim);
  auto mulStart = chrono::high_resolution_clock::now();
  auto mulShare = mul0<u64>(chl, mulTriple.first, spShare);
  auto mulEnd = chrono::high_resolution_clock::now();
  auto mulTime = chrono::duration_cast<chrono::microseconds>(mulEnd - mulStart).count();
  u64 mulTraffic = chl.bytesSent() + chl.bytesReceived() - spTraffic - edTraffic;
  cout << "MUL time: " << mulTime
       << "us with transferred " << mulTraffic << " "
       << (mulTraffic >> 10) << "KiB" << endl;

  cp::sync_wait(chl.flush());

  auto localStart = chrono::high_resolution_clock::now();
  vector<u64> sqXshare(nc, 0);
  for (size_t i = 0; i < nc; i++) {
    for (size_t j = 0; j < edDim; j++) {
      sqXshare[i] += spShare[i][j] * spShare[i][j] + 2 * mulShare[i][j];
    }
  }
  vector<u64> sqY(n, 0);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < edDim; j++) {
      sqY[i] += edy0[j][i] * edy0[j][i];
    }
  }
  vector<vector<u64>> R(nc, vector<u64>(n));
  for (size_t i = 0; i < nc; i++) {
    for (size_t j = 0; j < n; j++) {
      R[i][j] = sqXshare[i] - 2 * edShare[i][j] + sqY[j];
    }
  }
  auto localEnd = chrono::high_resolution_clock::now();
  auto localTime = chrono::duration_cast<chrono::microseconds>(localEnd - localStart).count();
  cout << "Local time: " << localTime << "us" << endl;
}

EigenFaceClientVT::EigenFaceClientVT(string db) {
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
  } else if (db == "test1") {
    vector<u8> t1(spDim);
    for (size_t i = 0; i < spDim; i++) {
      t1[i] = i;
    }
    t1[0] = 1;
    data.push_back(t1);
    nc = 1;
  } else {
    for (size_t j = 0; j < nc; j++) {
      data.push_back(randomUint8(spDim));
    }
  }

  auto chl = cp::asioConnect("127.0.0.1:7700", false);
  IknpOtExtReceiver receiver;
  cp::sync_wait(chl.recv(n));
  cp::sync_wait(chl.send(nc));

  vector<vector<u64>> spx1(nc, vector<u64>(spDim));
  for (size_t k = 0; k < spDim; k++) {
    for (size_t i = 0; i < nc; i++) {
      spx1[i][k] = data[i][k];
    }
  }
  vector<vector<u64>> spy1(spDim, vector<u64>(edDim, 0));

  auto spTriple = fakeZ<u64>(nc, edDim, spDim);

  auto spShare = mp1<u64, u64>(chl, spTriple.second, spx1, spy1);

  vector<vector<u64>> edx1(nc, vector<u64>(edDim));
  for (size_t i = 0; i < nc; i++) {
    for (size_t j = 0; j < edDim; j++) {
      edx1[i][j] = spShare[i][j];
    }
  }
  vector<vector<u64>> edy1(edDim, vector<u64>(n, 0));

  auto edTriple = fakeZ<u64>(nc, n, edDim);
  auto edShare = mp1<u64>(chl, edTriple.second, edx1, edy1);

  auto mulTriple = fakeTriple<u64>(nc * edDim);
  auto mulShare = mul1<u64>(chl, mulTriple.second, spShare);

  cp::sync_wait(chl.flush());

  vector<u64> sqXshare(nc, 0);
  for (size_t i = 0; i < nc; i++) {
    for (size_t j = 0; j < edDim; j++) {
      sqXshare[i] += spShare[i][j] * spShare[i][j] + 2 * mulShare[i][j];
    }
  }
  vector<vector<u64>> T(nc, vector<u64>(n));
  for (size_t i = 0; i < nc; i++) {
    for (size_t j = 0; j < n; j++) {
      T[i][j] = sqXshare[i] - 2 * edShare[i][j];
    }
  }
}
