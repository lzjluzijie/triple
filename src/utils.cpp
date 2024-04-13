#include "triple.hpp"

const string base64Table =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// Function to encode a bool array to Base64
std::string boolsToBase64(const std::vector<uint8_t> &data) {
  auto size = data.size();
  std::string base64String;
  size_t i = 0;
  size_t j = 0;
  int value = 0;

  while (i < size) {
    value |= data[i++] << (5 - j++);
    if (j == 6 || i == size) {
      base64String += base64Table[value];
      j = 0;
      value = 0;
    }
  }

  while (base64String.length() % 4 != 0) {
    base64String += "=";  // Add padding if needed
  }

  return base64String;
}

// Function to decode a Base64 string to a bool array
vector<uint8_t> base64ToBools(const string &base64String) {
  vector<uint8_t> bools;

  for (char c : base64String) {
    if (c == '=') {
      break;  // Padding, stop decoding
    }
    size_t pos = base64Table.find(c);
    if (pos != std::string::npos) {
      for (int k = 5; k >= 0; --k) {
        bools.push_back((pos >> k) & 1);
      }
    }
  }
  return bools;
}

// Function to encode a uint8_t vector to base64
std::string base64_encode(const vector<uint8_t> &data) {
  std::string ret;
  size_t i = 0;
  uint32_t x = 0;
  int j = 0;
  for (; i + 3 < data.size(); i += 3) {
    x = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    for (j = 0; j < 4; j++) {
      ret.push_back(base64Table[(x >> (6 * (3 - j))) & 0x3f]);
    }
  }
  if (i < data.size()) {
    x = data[i] << 16;
    ret.push_back(base64Table[(x >> 18) & 0x3f]);
    ret.push_back(base64Table[(x >> 12) & 0x3f]);
    ret.push_back('=');
    ret.push_back('=');
  } else if (i + 1 < data.size()) {
    x = (data[i] << 16) | (data[i + 1] << 8);
    ret.push_back(base64Table[(x >> 18) & 0x3f]);
    ret.push_back(base64Table[(x >> 12) & 0x3f]);
    ret.push_back(base64Table[(x >> 6) & 0x3f]);
    ret.push_back('=');
  } else if (i + 2 < data.size()) {
    x = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    ret.push_back(base64Table[(x >> 18) & 0x3f]);
    ret.push_back(base64Table[(x >> 12) & 0x3f]);
    ret.push_back(base64Table[(x >> 6) & 0x3f]);
    ret.push_back(base64Table[x & 0x3f]);
  }
  return ret;
}

// Function to decode a base64 string to a uint8_t array
std::vector<uint8_t> base64_decode(const std::string &encoded) {
  std::vector<uint8_t> ret;
  size_t i = 0;
  uint32_t x = 0;
  int j = 0;
  for (; i + 4 < encoded.size(); i += 4) {
    x = (base64Table.find(encoded[i]) << 18) |
        (base64Table.find(encoded[i + 1]) << 12) |
        (base64Table.find(encoded[i + 2]) << 6) |
        base64Table.find(encoded[i + 3]);
    for (j = 0; j < 3; j++) {
      ret.push_back((x >> (8 * (2 - j))) & 0xff);
    }
  }
  if (i + 1 < encoded.size()) {
    x = (base64Table.find(encoded[i]) << 18);
    x |= (base64Table.find(encoded[i + 1]) << 12);
    ret.push_back((x >> 16) & 0xff);
  }
  if (i + 2 < encoded.size()) {
    x = (base64Table.find(encoded[i]) << 18);
    x |= (base64Table.find(encoded[i + 1]) << 12);
    x |= (base64Table.find(encoded[i + 2]) << 6);
    ret.push_back((x >> 16) & 0xff);
    ret.push_back((x >> 8) & 0xff);
  }
  return ret;
}

vector<uint8_t> randomBools(size_t l) {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_int_distribution<uint8_t> r(0, 1);
  vector<uint8_t> res(l);
  // cout << sizeof(res) << endl;
  for (size_t i = 0; i < l; i++) {
    res[i] = r(gen) == 1;
  }
  return res;
}

vector<uint8_t> randomUint8(size_t l) {
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_int_distribution<uint8_t> r(0, 255);
  vector<uint8_t> res(l);
  for (size_t i = 0; i < l; i++) {
    res[i] = r(gen);
  }
  return res;
}

uint64_t hd(const vector<uint8_t> &arr1, const vector<uint8_t> &arr2) {
  if (arr1.size() != arr2.size()) {
    throw std::invalid_argument("arr1.size() != arr2.size()");
  }
  uint64_t distance = 0;
  for (size_t i = 0; i < arr1.size(); ++i) {
    if (arr1[i] != arr2[i]) {
      ++distance;
    }
  }
  return distance;
}
