# SafetyNets

SafetyNets provides a framework which allows a client to provably verify the correctness of deep neural network inference computations, outsourced to an untrusted cloud. By leveraging techniques in interactive proof systems, the cloud server can provide (in addition to the result) a *proof* of the correctess of the result. Read more in our [paper](http://papers.nips.cc/paper/7053-safetynets-verifiable-execution-of-deep-neural-networks-on-an-untrusted-cloud.pdf).


## Build
`safetynets` implements the interactive proof protocol, and measures the running time of the client (verifier) and the server (prover). To build and use the framework, run:
```shell
$ make
$ ./safetynets <arch filepath>
```
#### Example
```shell
$ ./safetynets.o timit_arch.txt
```

## Usage
`safetynets` takes as input a file `<arch filepath>` containing the input batch size and the fully connected network architecture. The network architecture is described as input size and the number of neurans in each layer. As an example, `timit_arch` describes an input batch size of 512, and a neural network with input size of 1845, 3 hidden layers of 2000 neurans each, and an output size of 183. Therefore, `timit_arch` contains: 
```txt
512
1845
2000
2000
2000
183
```
Please note that convolutional neural networks should be first converted to their fully connected equivalents.


