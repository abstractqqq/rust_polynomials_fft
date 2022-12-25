# Polynomials and Fast Fourier Transform

A learning project in which I defined a polynomial struct and implemented FFT for polynomial multiplication, with everything written in Rust.

# Retro

This is my summary of things I did in this project.

1. Algorithms
    - Basic mathematical operations for polynomials, including long division.
    - Raise any polynomial to a power n, where n is a non-negative integer, by the use of the "squaring" algorithm. 
    - FFT for polynomial multiplication.
    - Pretty print a polynomial.

2. Rust
    - Unit testing.
    - Use of trait bounds to make general polynomials.
    - Other Rust syntax and concepts.
    - Getting more familiar with ndarray.
    
# Implementation Details

1. Polynomials\<T\>, polynomials over T, where T: Num + Clone + Copy + Display
    - T is like a mathematical field.
    - Mathematically speaking, T can be more abstract, say T can be any mathematical ring. But then divide_by will not make sense and has to be rewritten.
    - If we ditch the Div trait bound in Num(NumOps), then we will get something close to a mathematical ring.

2. Implemented plus, minus, multiply, divide_by (long division) for polynomials, and you may use p1 + p2, p1 - p2, and p1 * p2.

3. The normal multiply method on the Polynomial struct uses the O(n^2) way. FFT for multiplication follows the [Cooley–Tukey FFT algorithm.](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm). FFT is  implemented for Polynomial\<f64\>. It can be extended to Polynomial\<Complex64\> easily. If one wants to extend this to other more general fields, one has to define more complicated structures and define n-th roots of unity in those fields (e.g. Finite fields, etc...), which is very beyond this project.

4. FFT has float point precision issues, so I rounded everything to 5 decimal places.

5. It is built for fun and I am not aiming for the utmost speed, as this is more of a learning project.

# Resources:

1. [Learning FFT](https://www.youtube.com/watch?v=h7apO7q16V0&t=1265s)
2. [FFT Lecture Notes](http://www.cs.toronto.edu/~denisp/csc373/docs/tutorial3-adv-writeup.pdf)
3. Documentation of ndarray, num_traits, etc.
