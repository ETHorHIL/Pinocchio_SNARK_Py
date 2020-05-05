# Pinocchio_Py_TH
My "for fun" Python Implementation of Pinocchio.

I took inspiration from these sources:

* why-and-how-zk-snark-works(https://medium.com/@imolfar/why-and-how-zk-snark-works-8-zero-knowledge-computation-f120339c2c55)

* babySNARK(https://github.com/initc3/babySNARK)

* Vitaliks R1CS compiler(https://github.com/ethereum/research/blob/master/zksnark/code_to_r1cs.py)

SNARK.py is missing the optional calculations to make it Zero Knowledge. zk-SNARK.py has those calculations and is a bit slower. I included some trivial examples to see it running: A NAND gate and a check if a number is binary. I have also included Vitaliks code_to_r1cs.py in my SNARK.py file which allows you to compile any code down to R1CS and prove that its assignment checks out.

# Protocol:

![Image description](https://github.com/ETHorHIL/Pinocchio_Py_TH/blob/master/Proctocol_SNARK.png)
