This repo contains two short Python/Sage proof of concept implementations (about $50$ lines of code each) that enable **fast  multiplication of two polynomials** of **large degree** (up to $N=8192$) and **large integer coefficients** (up to $800$ bit integers). The idea is to use available implementations of Fast Fourier Transform (FFT) that are already optimized for floating-point numbers. Even though the floating-point FFTs produce only *approximate results*, they can be used to compute *exact results* in the case of arbitrary-precision polynomial multiplication.


The main tool that is leveraged to achieve this is [scipy's](https://pypi.org/project/scipy/) approximate Fast Fourier Transform (FFT) that makes use of complex numbers. The general idea is to break the two polynomials we want to multiply into small-coefficient polynomials of the same degree. $\texttt{scipy}$ can fast multiply these small-coefficient polynomials (FFT-based multiplication takes ${O}(N\log N)$ time), using floating-point arithmetic, such that we get a correct result when the floating-point errors are shaved off. In the end, all these partial computations are used to construct the integer polynomial that is the product of the two initial polynomials. To break the coefficients, two methods are used, each with a separate implementation:

**a.** **Base-decomposition** 

[bd_pmul.py](https://github.com/rtitiu/pmul-approx-ffts/blob/main/bd_pmul.py)


**b.** **Chinese Remainder Theorem (CRT)**     

[crt_pmul.sage](https://github.com/rtitiu/pmul-approx-ffts/blob/main/crt_pmul.sage)

:arrow_right: For more details check out the implementation_details.pdf file! 