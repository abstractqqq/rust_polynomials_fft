# Polynomials 

1. A vec representation of polynomials over real numbers. Obviously not the most general, as polynomials can be defined over any arbitrary rings. Obviously not the most space efficient, as a lot of space is wasted for polynomials like x^5, which is internally represented as vec![0,0,0,0,1].
2. Add, subtract, multiply, divide (long division)
3. It is built for fun, and eventually I will build a fft for multiplication for polynomials. 
4. I do not aim for the utmost speed, as this is more of a learning project.


## Learning Rust

0. Tring to implement FFT for polynomial multiplication (Or convolution)
2. FFT implementation for the polynomial struct (multiplication). 