#ifndef TRIPLE_H
#define TRIPLE_H

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <thread>

using namespace std;

string boolsToBase64(const vector<uint8_t> &data);
vector<uint8_t> base64ToBools(const string &base64String);
string base64_encode(const vector<uint8_t> &data);
vector<uint8_t> base64_decode(const string &encoded);
vector<uint8_t> randomBools(size_t l);
vector<uint8_t> randomUint8(size_t l);
uint64_t hd(const vector<uint8_t> &arr1, const vector<uint8_t> &arr2);
void test();
#endif
