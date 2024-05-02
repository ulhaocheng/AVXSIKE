# AVXSIKE

AVXSIKE is a software library of AVX-512 implementations (using AVX-512IFMA
extension) of Supersingular Isogeny Key Encapsulation (SIKE).

It contains the low-latency implementation AVXSIKE-LL and the high-throughput
implementation AVXSIKE-HT.

## The Structure of AVXSIKE Software Library
<pre> AVXSIKE
      ├── HT (High-Throughput) 
      |      ├── SIKEp434 
      |      ├── SIKEp503     
      |      ├── SIKEp610     
      |      └── SIKEp751         
      |
      └── LL (Low-Latency) 
             ├── SIKEp434 
             ├── SIKEp503     
             ├── SIKEp610     
             └── SIKEp751  
</pre>

## Usage

### Compile the AVXSIKE implementation (need Intel Cannon Lake, Ice Lake or Tiger Lake machine!)

```bash
    $ cd AVXSIKE/AVXSIKE-[HT/LL]/SIKEp[434/503/610/751] 
    $ make 
    $ ./sike
```
### KAT test to verify the correctness of AVXSIKE-LL

```bash
    $ cd AVXSIKE/AVXSIKE-LL/SIKEp[434/503/610/751] 
    $ make kat
    $ ./kat
```

## Paper
A paper describing the various implementations in library has been published in
*IACR Transactions on Cryptographic Hardware and Embedded Systems, 2022(2),
41-68*.

The paper online link is
[here](https://tches.iacr.org/index.php/TCHES/article/view/9480/9021).

## Software Author
Hao Cheng (University of Luxembourg).

## Copyright
Copyright (C) 2021-2022 by University of Luxembourg.

## LICENSE
GPLv3 (see details in LICENSE file).
