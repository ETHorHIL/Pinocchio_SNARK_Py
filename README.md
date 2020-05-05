# Pinocchio_Py_TH
My "for fun" Python Implementation of Pinocchio.

I took inspiration from these sources:
https://medium.com/@imolfar/why-and-how-zk-snark-works-8-zero-knowledge-computation-f120339c2c55

https://github.com/initc3/babySNARK

https://github.com/ethereum/research/blob/master/zksnark/code_to_r1cs.py

SNARK.py is missing the optional calculations to make it Zero Knowledge. zk-SNARK.py has those calculations and is a bit slower.

There are two trivial circuits that you can run. A NAND gate and check if a number is binary. The SNARK.py file also contains a function that allows you to compile your
own code to R1CS and prove an assignment to that. The compiler, code_to_r1cs.py
was written by Vitalik Buterin.

# Protocol:

![Image description](https://github.com/ETHorHIL/Pinocchio_Py_TH/blob/master/Proctocol_SNARK.png)
