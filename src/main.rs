mod polynomial;
use polynomial::Polynomial;

use std::f64::consts::PI;
use num_traits::{Zero, One};
use num_complex::*;
use ndarray::{Array1, s, ArrayView1};


fn fft(p:ArrayView1<f64>) -> Array1<Complex64> {

    let n = p.len();
    if n == 1 { 
        return Array1::from_elem(1, Complex::new(p[0], 0.))
    }
    let deg = 2.0 * PI / (n as f64);
    let w_n = Complex64::new(deg.cos(), deg.sin());
    let even = p.slice(s![..;2]);
    let odd = p.slice(s![1..p.len();2]);
    let y_e = fft(even);
    let y_o = fft(odd);
    let mut y = Array1::from_elem(n, Complex64::zero());
    let mut w = Complex64::one();
    for j in 0..(n/2) {
        let odd_term = w * y_o[j];
        y[j] = y_e[j] + odd_term;
        y[j + n/2] = y_e[j] - odd_term;
        w *= w_n;
    }
    return y
}

fn inverse_fft(p:ArrayView1<Complex64>) -> Array1<Complex64> {

    let n = p.len();
    if n == 1 { 
        return Array1::from_elem(1, p[0])
    }

    let divisor = n as f64;
    let deg = -2.0 * PI / divisor;
    let w_n = Complex64::new(deg.cos(), deg.sin());
    let even = p.slice(s![..;2]);
    let odd = p.slice(s![1..p.len();2]);
    let y_e = inverse_fft(even);
    let y_o = inverse_fft(odd);
    let mut y = Array1::from_elem(n, Complex64::zero());
    let mut w = Complex64::one();
    for j in 0..(n/2) {
        let odd_term = w * y_o[j];
        y[j] = y_e[j] + odd_term;
        y[j + n/2] = y_e[j] - odd_term;
        w *= w_n;
    }
    return y
}

fn smallest_pow_2(n:usize) -> usize {
    // smallest power of 2 that is >= n
    let mut test:usize = 1;
    while test < n {
        test <<= 1;
    }
    test
}

fn add_leading_zeros(p:&Polynomial<f64>, n:usize) -> Polynomial<f64>{
    let q = p.get_coeffs_view();
    let mut new_array = Array1::from_elem(n, 0.);
    for i in 0..q.len() {
        new_array[i] = q[i];
    }
    Polynomial::new(new_array)
}

fn value_repr(p:Polynomial<f64>) -> Array1<Complex64> {
    // it is assumed that the input to this function is already of len 2^n
    fft(p.get_coeffs_view())
}

fn fft_mul(p:&Polynomial<f64>, q:&Polynomial<f64>) -> Polynomial<f64> {
    // q.deg + p.deg + 1 = length of the output
    let target_len = smallest_pow_2(q.deg() + p.deg() + 1);
    // I could potentially launch two threads to do this part
    let new_p = add_leading_zeros(p, target_len);
    let new_q = add_leading_zeros(q, target_len);
    // pointwise multiplication, then apply inverse
    let product_ptwise = value_repr(new_p) * value_repr(new_q);
    let new_coeff_array = inverse_fft(product_ptwise.view())/(target_len as f64);
    println!("Product: {:?}", new_coeff_array);
    // extract real parts and return, round real parts to 5 decimals.
    Polynomial::no_leading_zeros(new_coeff_array.map(|z| (z.re()*100000.).trunc()/100000.).to_vec())
}




fn main() {
    // let p1 = vec![0.5,1.,-2.,-3.];
    // let p2 = vec![0.5,-1.,3.,4.];
    // //let p2 = vec![0., -1., 2., 3.];
    // let q1 = Polynomial {coeffs: p1};
    // let q2 = Polynomial{coeffs: p2};
    // let q3 = q1.add(&q2);
    // println!("{}", q3);
    // println!("Evaluating at 1 is: {}", q3.eval(1.));

    // let p3 = vec![0.,1.,1.];
    // let p4 = vec![0., 2.];
    // let q3 = Polynomial{coeffs: p3};
    // let q4 = Polynomial{coeffs: p4};
    // println!("Multiply {} by {} is:\n{}", q3, q4, q3.multiply(&q4));

    // // 1 + x + x^5 + 2x^6
    // let p1 = Polynomial::new(vec![1., 1., 0., 0., 0., 1., 2.]);
    // // 1 + x^2
    // let p2 = Polynomial::new(vec![1., 0., 1.]);
    // let (quotient, remainder) = p1.divide_by(&p2);

    // println!("The quotient is {}. \nThe remainder is {}.", quotient, remainder);

    // let p1 = Polynomial::new(vec![-1., 0., 0., 0., 0., 0., 0., 1.]);
    // let p2 = Polynomial::new(vec![-1., 1.]);
    // let (quotient, remainder) = p1.divide_by(&p2);

    // println!("The quotient is {:?}. \nThe remainder is {:?}.", quotient, remainder);


    // let p1 = Polynomial::new(vec![-1., 1.]);
    // let n = 7;
    // let p1_power = p1.to_power(n);

    // println!("Raising ({:?}) to the {:?}th power yields\n{:?}.", p1, n, p1_power);

    let p1 = Polynomial::from_vec(vec![-1., 1.]);
    let p2 = Polynomial::from_vec(vec![1., 1., 1., 1.]);
    println!("The product is: {:?}", fft_mul(&p1, &p2));
    println!("The product is: {:?}", p1.multiply(&p2));

}
