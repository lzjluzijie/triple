# Triple

https://eprint.iacr.org/2023/1863

This project is the implementation of the paper "Efficient Secure Multiparty Computation for Multidimensional Arithmetics and Its Application in Privacy-Preserving Biometric Identification".  The functionalities of this project include generating tensor(vector) triples, using tensor triples to compute scalar products, tensor(outer) products, and matrix multiplications. There are also examples of applications in biometric identification, including the FingerCodes and EigenFace.

Examples and benchmarks are located in `src/test.cpp`. The `noisy` functions stands for the Correlated OT(COT) based protocol in the paper, and the `silent` functions stands for the Silent OT(SOT). The words `noisy` and `silent` are originated from the two types of VOLE in `libOTe`, and we interchangeably use `vt`, `vector triple`, and `tensor triple` in the project and the paper.

## Build

This project depends on `boost` and `libOTe`. The author used `gcc 12.2.0`.

```shell
cd libOTe
python build.py --all --boost --sodium
cd .. 
mkdir build
cd build
cmake ..
make -j
```

## Contact

Contact `lzjluzijie@gmail.com` if you have any question.
