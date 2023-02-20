use itertools::{EitherOrBoth::*, Itertools};
use ndarray::{Array1, ArrayView1, s};
use std::f64::consts::PI;
use std::thread;
use num_complex::*;
use num_traits::{Num, Zero, One};
use std::{cmp::PartialEq, fmt::Display, ops::{Add, Sub, Mul}};

//-------------------------------------------------------------------------------------------------------------
// Abstract Implementation of polynomials
#[derive(Debug)]
pub struct Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    coeffs: Array1<T>
}

impl <T> Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    #[inline]
    pub fn new(c: Array1<T>) -> Polynomial<T> {
        Polynomial::no_leading_zeros(c.to_vec())
    }

    #[inline]
    pub fn from_vec(c: Vec<T>) -> Polynomial<T> {
        Polynomial::no_leading_zeros(c)
    }

    /// If we get something like vec![0,0,0,1,0], then we should remove this last 0.
    /// The polynomial is vec![0,0,0,1]. We only keep 0 if the polynomial is vec![0], the zero polynomial.
    /// This method will consume c.
    /// 
    /// returns a polynomial.
    pub fn no_leading_zeros(mut c:Vec<T>) -> Polynomial<T> {
        if c.len() == 0 {
            panic!("Cannot generate polynomial from empty vec.")
        }
        while c.len() > 1 {
            let last = c.last().unwrap();
            if *last == T::zero() {
                c.pop();
            } else {
                break
            }
        }
        Polynomial{coeffs: Array1::from_vec(c)}
    }

    #[inline]
    pub fn copy(&self) -> Polynomial<T> {
        Polynomial{coeffs: self.coeffs.clone()}
    }

    // Utility
    pub fn eval(&self, x:T) -> T {
        if x.is_zero() {
            self.coeffs[0]
        } else {
            let mut x_pow = T::one();
            self.coeffs.iter().enumerate().fold(
                T::zero(), |acc, (power, coef)| {
                    if power == 0 {
                        *coef
                    } else {
                        x_pow = x_pow * x;
                        acc + *coef * x_pow
                    }
                }
            )
        }
    }

    pub fn deg(&self) -> usize {
        let n = self.coeffs.len();
        match n {
            0|1 => 0,
            _ => n-1 // complier knows this is > 0 and therefore always a usize!!!
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    #[inline]
    pub fn get_coeffs(&self) -> Array1<T> {
        self.coeffs.clone()
    }

    #[inline]
    pub fn get_coeffs_view(&self) -> ArrayView1<T> {
        self.coeffs.view()
    }

    #[inline]
    pub fn get_const(&self) -> &T {
        self.coeffs.get(0).unwrap()
    }

    /// Polynomial of degree n with constant coefficient c
    pub fn const_coef(c:T, n:usize) -> Polynomial<T> {
        if c == T::zero() {
            return Polynomial::zero()
        } else {
            match n {
                0 => Polynomial{coeffs: Array1::from_elem(1, T::zero())},
                _ => Polynomial{coeffs: Array1::from_elem(n, c)}
            }
        }
    }

    /// Creates the n basis polynomial with coefficient c.
    /// E.g. n = 0, c = 1, returns 1 as a polynomial
    /// n = 2, c = 2 returns 2x^2
    pub fn basis(c:T, n:usize) -> Polynomial<T> {
        if c == T::zero() {
            return Polynomial::zero()
        } else {
            match n {
                0 => Polynomial{coeffs: Array1::from_elem(1, c)},
                _ => {
                    let mut v = vec![T::zero(); n];
                    v.push(c);
                    Polynomial{coeffs: Array1::from_vec(v)}
                }
            }
        }
    }

    #[inline]
    pub fn highest_coeff(&self) -> T {
        *self.coeffs.last().unwrap()
    }

    // Arithmetic
    pub fn plus(&self, p:&Polynomial<T>) -> Polynomial<T> {
        let mut new_poly:Vec<T> = Vec::with_capacity(self.coeffs.len().max(p.coeffs.len()));
        for pair in self.coeffs.iter().zip_longest(p.coeffs.iter()) {
            match pair {
                Both(l,r) => new_poly.push(*l+*r),
                Left(l) => new_poly.push(*l),
                Right(r) => new_poly.push(*r)
            }
        }
        // it is possible that the leading term becomes 0 after adding/subtracting
        // so we use new to deal with that case
        Polynomial::no_leading_zeros(new_poly)

    }

    #[inline]
    pub fn minus(&self, p:&Polynomial<T>) -> Polynomial<T> {
        let p2 = Polynomial{coeffs: p.coeffs.iter().map(|x| T::zero()-*x).collect()};
        self.plus(&p2)
    }

    pub fn multiply(&self, p:&Polynomial<T>) -> Polynomial<T> {
        let mut new_poly:Vec<T> = vec![T::zero(); self.coeffs.len() + p.coeffs.len() - 1];
        for (i,a) in self.coeffs.iter().enumerate() {
            for (j,b) in p.coeffs.iter().enumerate() {
                let k = i + j;
                new_poly[k] = new_poly[k] + (*a)*(*b);
            }
        }
        Polynomial::no_leading_zeros(new_poly)
    }

    /// Long Division
    pub fn divide_by(&self, p:&Polynomial<T>) -> Option<(Polynomial<T>, Polynomial<T>)> {
        // (P1, P2) = (Quotient, Remainder)
        if p.is_zero(){
            return None
        }

        let dividee = self.copy();
        let dividee_deg = dividee.deg();
        let divider_deg = p.deg();
        if dividee_deg < divider_deg {
            return Some((Polynomial::zero(), Polynomial{coeffs: p.coeffs.clone()}))
        }
        let mut quotient = vec![T::zero(); dividee_deg - divider_deg + 1];
        let remainder = Polynomial::_long_div(dividee, p, &mut quotient);
        Some((Polynomial{coeffs: Array1::from_vec(quotient)}, remainder))
    }


    fn _long_div(
        dividee:Polynomial<T>
        , divider:&Polynomial<T>
        , quotient:&mut Vec<T>
    ) -> Polynomial<T> { // returns the remainder    
        let dividee_deg = dividee.deg();
        let divider_deg = divider.deg();
        if dividee_deg < divider_deg {// dividee becomes remainder
            return  dividee
        } else {
            // dividee_deg >= divider_deg.
            let new_term_deg = dividee_deg - divider_deg; // always >= 0
            let new_term_coeff = dividee.highest_coeff() / divider.highest_coeff();
            // modify the quotient
            if new_term_coeff == T::zero() || quotient[new_term_deg] != T::zero() {
                // This means that we did not reduce degree. 
                // This might happen when we are working with Polynomials over integers.
                // E.g x^2 + 1 divided by 2x^2 in Z[x], or 3x^2 + 1 divided by 2x^2
                return dividee // this is the remainder
            } else { // first time setting coeff for this deg
                quotient[new_term_deg] = new_term_coeff;
                // update dividee
                let basis = Polynomial::basis(new_term_coeff, new_term_deg);
                let new_dividee = dividee.minus(&basis.multiply(divider));
                // println!("Division steps: dividing {} by {}", new_dividee, divider);
                Self::_long_div(new_dividee, divider, quotient)
            }
        }
    }

    // also known as formal derivative
    pub fn ddx(&self) -> Polynomial<T> {
        let deg = self.deg();
        match deg {
            0 => Polynomial::zero(),
            _ => {
                Polynomial{
                    coeffs:self.coeffs.slice(s![1..=deg]).iter()
                                .enumerate()
                                .map(|(i, v)| Self::_fast_self_add(*v, i+1))
                                .collect()
                }
            }
        }
    }

    // add T times, e.g. _fast_self_add(1, 6) = 1 + 1 + 1 + 1 + 1 + 1 = 6
    // This might seem dumb, but this works for general T, even when T is not real or complex numbers.
    fn _fast_self_add(value:T, times:usize) -> T {
        if value == T::zero() {
            return value;
        } else {
            match times {
                0|1 => value,
                _ => {
                    let two_sum = value + value;
                    if times % 2 == 1 {
                        Self::_fast_self_add(two_sum, (times-1) >> 1) + value
                    } else {
                        Self::_fast_self_add(two_sum, times >> 1)
                    }
                }
            }
        }   
    }

    // raise a polynomial p to a deg.
    pub fn pow(&self, n:usize) -> Polynomial<T> {
        match n {
            0 => {
                    println!("DON'T DO THIS.");
                    Polynomial::one()
                },
            _ => {
                let cur = self.copy();
                Self::_power(cur, n)
            }
        }
    }

    fn _power(current:Polynomial<T>, deg:usize) -> Polynomial<T> {
        match deg {
            0|1 => current,
            _ => {
                let squared = current.multiply(&current);
                if deg % 2 == 1 {
                    Self::_power(squared, (deg-1) >> 1).multiply(&current)
                } else {
                    Self::_power(squared, deg >> 1)
                }
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------------------
// implementations for mathmatical traits

impl <T> PartialEq for Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs 
    }

    #[inline]
    fn ne(&self, other: &Self) -> bool {
        self.coeffs != other.coeffs
     }
}

impl <T> Add for Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    type Output = Self;
    fn add(self, p:Polynomial<T>) -> Polynomial<T> {
        self.plus(&p)
    }
}

impl <T> Add for &Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    type Output = Polynomial<T>;
    fn add(self, p:&Polynomial<T>) -> Polynomial<T> {
        self.plus(p)
    }
}


impl <T> Sub for Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    type Output = Self;
    fn sub(self, p:Polynomial<T>) -> Polynomial<T> {
        self.minus(&p)
    }
}

impl <T> Sub for &Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    type Output = Polynomial<T>;
    fn sub(self, p:&Polynomial<T>) -> Polynomial<T> {
        self.minus(&p)
    }
}

impl <T> Mul for Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    type Output = Self;
    fn mul(self, p:Polynomial<T>) -> Polynomial<T> {
        self.multiply(&p)
    }
}

impl <T> Mul for &Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    type Output = Polynomial<T>;
    fn mul(self, p:&Polynomial<T>) -> Polynomial<T> {
        self.multiply(&p)
    }
}

impl <T> Zero for Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    #[inline]
    fn zero() -> Self {
        Polynomial { coeffs: Array1::from_elem(1, T::zero())}
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.coeffs == Array1::from_elem(1, T::zero())
    }
}

impl <T> One for Polynomial<T> 
    where T: Num + Clone + Copy + Display
{
    #[inline]
    fn one() -> Self {
        Polynomial { coeffs: Array1::from_elem(1, T::one())}
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.coeffs == Array1::from_elem(1, T::one())
    }
}

//-------------------------------------------------------------------------------------------------------------
// This is pretty print for Polynomials
impl <T> std::fmt::Display for Polynomial<T>
    where T: Num + Clone + Copy + Display
{
    
    fn fmt(&self, f:&mut std::fmt::Formatter) -> std::fmt::Result {
        let l = self.coeffs.len();
        if l == 0 {
            return write!(f, "Empty Polynomial. (Empty coefficients)")
        }
        let mut p = String::new();
        for (i,coef) in self.coeffs.iter().enumerate().rev() {
            if coef.is_zero() { 
                if i == 0 && p.len() == 0 {
                    p.push_str(&coef.to_string());
                }    
                continue
            }
            let mut term = String::new();
            if i == l - 1 { 
                // leading term, no need to manually add signs
                if self.deg() == 0 {
                    // push everything if poly is constant
                    term.push_str(&coef.to_string());                
                } else { // poly is not constant and must have a x^k term, k > 0, don't print 1.
                    if !coef.is_one() {
                        term.push_str(&coef.to_string());
                    }
                }
            } else { // not leading term
                if i == 0 { // constant term
                    if !coef.is_zero() {
                        term.push_str(&" + ");
                        let coef_str = &coef.to_string();
                        if coef_str.starts_with("-") { // ad-hoc catch of 'negative' values of T
                            term.push('(');
                            term.push_str(coef_str);
                            term.push(')');
                        } else {
                            term.push_str(coef_str);
                        }
                    }
                } else { // non leading, non constant terms
                    term.push_str(&" + ");
                    if !coef.is_one() {
                        let coef_str = &coef.to_string();
                        if coef_str.starts_with("-") { // ad-hoc catch of 'negative' values of T
                            term.push('(');
                            term.push_str(coef_str);
                            term.push(')');
                        } else {
                            term.push_str(coef_str);
                        }
                    }
                }
            }
            if i > 0 {
                term.push('x');
            }
            if i > 1 {
                term.push('^');
                term.push_str(&i.to_string());
            }
            p.push_str(&term);
        }
        write!(f, "{}", p)
    }
}

//-------------------------------------------------------------------------------------------------------------
// Only for real polynomials. Can be made slightly more general. But I will stop here.
impl Polynomial<f64> {

    #[inline]
    pub fn get_value_repr(&self) -> Array1<Complex64> {
        Self::fft(self.coeffs.view())
    }

    fn fft(p:ArrayView1<f64>) -> Array1<Complex64> {
        // returns the value representation of p, 
        let n = p.len();
        if n == 1 {
            return Array1::from_elem(1, Complex::new(p[0], 0.))
        }
        let deg = 2.0 * PI / (n as f64);
        let w_n = Complex64::new(deg.cos(), deg.sin());
        let even = p.slice(s![..;2]);
        let odd = p.slice(s![1..p.len();2]);
        let y_e = Self::fft(even);
        let y_o = Self::fft(odd);
        let mut y = Array1::from_elem(n, Complex64::zero());
        let mut w = Complex64::one();
        let half = n/2;
        for j in 0..half {
            let odd_term = w * y_o[j];
            y[j] = y_e[j] + odd_term;
            y[j + half] = y_e[j] - odd_term;
            w *= w_n;
        }
        y
    }


    fn inverse_fft(p:ArrayView1<Complex64>) -> Array1<Complex64> {
    
        let n = p.len();
        if n == 1 { 
            return Array1::from_elem(1, p[0])
        }
    
        let deg = -2.0 * PI / (n as f64);
        let w_n = Complex64::new(deg.cos(), deg.sin());
        let even = p.slice(s![..;2]);
        let odd = p.slice(s![1..p.len();2]);
        let y_e = Self::inverse_fft(even);
        let y_o = Self::inverse_fft(odd);
        let mut y = Array1::from_elem(n, Complex64::zero());
        let mut w = Complex64::one();
        let half = n/2;
        for j in 0..half {
            let odd_term = w * y_o[j];
            y[j] = y_e[j] + odd_term;
            y[j + half] = y_e[j] - odd_term;
            w *= w_n;
        }
        y
    }

    fn smallest_pow_2(n:usize) -> usize {
        // smallest power of 2 that is >= n
        let mut test:usize = 1;
        while test < n {
            test <<= 1;
        }
        test
    }
    
    fn with_leading_zeros(&self, target_len:usize) -> Array1<f64>{
        let q = self.get_coeffs_view();
        let mut new_array = Array1::from_elem(target_len, 0.);
        for i in 0..q.len() {
            new_array[i] = q[i];
        }
        new_array
    }
    
    /// Performs polynomial multiplication using FFT
    /// 
    /// self
    /// q: &Polynomial<f64>
    /// decimal_places: how many decimal places should the output keep? FFT has some small numerical error. 
    /// Although it is small, it will often give 0.000000000000012312 instead of 0. 
    /// You can decide what precision you need for your output.
    /// 
    /// returns: product of self and q 
    pub fn fft_mul(&self, q:&Polynomial<f64>, decimal_places:usize) -> Polynomial<f64> {
        let q_deg = q.deg();
        let p_deg = self.deg();
        if q_deg == 0 || p_deg == 0 {
            return self.multiply(q)
        }
        // q.deg + p.deg + 1 = length of the output
        let target_len = Self::smallest_pow_2(p_deg + q_deg + 1);
        // 
        let new_p = self.with_leading_zeros(target_len);
        let new_q = q.with_leading_zeros(target_len);
        // pointwise multiplication, then apply inverse
        let product_ptwise = Self::fft(new_p.view()) * Self::fft(new_q.view());
        let new_coeff_array = Self::inverse_fft(product_ptwise.view())/(target_len as f64);
        // extract real parts and return, round real parts to 5 decimals.
        let rounding_factor = (10.).powi(decimal_places as i32);
        Polynomial::no_leading_zeros(new_coeff_array.map(|z| ((z.re()*rounding_factor).trunc())/rounding_factor).to_vec())
    }

    pub fn fft_mul_threaded(&self, q:&Polynomial<f64>, decimal_places:usize) -> Polynomial<f64> {
        let q_deg = q.deg();
        let p_deg = self.deg();
        if q_deg == 0 || p_deg == 0 {
            return self.multiply(q)
        }
        // q.deg + p.deg + 1 = length of the output
        let target_len = Self::smallest_pow_2(q_deg + p_deg + 1);
        // 
        let (fft_p, fft_q) = thread::scope(|s| {
            let first = s.spawn(|| 
                Self::fft(self.with_leading_zeros(target_len).view())
            );
            let second = s.spawn(|| 
                Self::fft(q.with_leading_zeros(target_len).view())
            );
            (first.join().unwrap(), second.join().unwrap())
        });
        // pointwise multiplication, then apply inverse
        let product_ptwise = fft_p * fft_q;
        let new_coeff_array = Self::inverse_fft(product_ptwise.view())/(target_len as f64);
        // extract real parts and return, round real parts to 5 decimals.
        let rounding_factor = (10.).powi(decimal_places as i32);
        Polynomial::no_leading_zeros(new_coeff_array.map(|z| ((z.re()*rounding_factor).trunc())/rounding_factor).to_vec())
    }


}
