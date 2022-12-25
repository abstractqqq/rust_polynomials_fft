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

    let deg = -2.0 * PI / (n as f64);
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
    // extract real parts and return, round real parts to 5 decimals.
    Polynomial::no_leading_zeros(new_coeff_array.map(|z| ((z.re()*100000.).trunc())/100000.0).to_vec())
}




fn main() {

    let p1 = Polynomial::from_vec(vec![1,2,3]);
    let p2 = Polynomial::from_vec(vec![0,0,2]);
    let (quotient, remainder) = p1.divide_by(&p2).unwrap();
    println!("Dividing {} by {},", p1, p2);
    println!("The quotient is {},",quotient);
    println!("The remainder is {}.", remainder);
    println!("Hence {} = ({}) * ({}) + {}.\n\n", p1, quotient, p2, remainder);
    
    let q1 = Polynomial::from_vec(vec![-1.,1.]);
    let q2 = Polynomial::from_vec(vec![1.,1.,1.,1.,1.,1.,1.]);
    println!("Multiplying {} by {} ...", q1, q2);
    println!("The result is {}.",fft_mul(&q1, &q2));

}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_add_1() {
        let p1 = Polynomial::from_vec(vec![1, 2, 3, 4, 5]);
        let p2 = Polynomial::from_vec(vec![5, 4, 3, 2, -5]);
        assert_eq!(p1 + p2, Polynomial::from_vec(vec![6,6,6,6]));
    }

    // p1 + p2 internally calls plus. So the test is redudant. It is important to keep in mind that p1 + p2 will consume both.

    #[test]
    fn test_add_2() {
        let p1 = Polynomial::from_vec(vec![1, 2, 3, 4, 5]);
        let p2 = Polynomial::from_vec(vec![5, 4, 3, 2, -5]);
        assert_eq!(p1.plus(&p2), Polynomial::from_vec(vec![6,6,6,6]));
    }

    #[test]
    fn test_minus_1() {
        let p1 = Polynomial::from_vec(vec![1,2,3,4,5]);
        let p2 = Polynomial::from_vec(vec![2,1,3]);
        assert_eq!(p1.minus(&p2), Polynomial::from_vec(vec![-1,1,0,4,5]));
    }

    #[test]
    fn test_minus_2() {
        let p1 = Polynomial::from_vec(vec![1,2,3,4,5]);
        let p2 = Polynomial::from_vec(vec![1,2,3,4,5]);
        assert_eq!(p1.minus(&p2), Polynomial::from_vec(vec![0]));
    }

    #[test]
    fn test_multiply_1() {
        let p1 = Polynomial::from_vec(vec![-1,1]);
        let p2 = Polynomial::from_vec(vec![1,1,1,1,1,1,1]);
        assert_eq!(p1.multiply(&p2), Polynomial::from_vec(vec![-1,0,0,0,0,0,0,1]));
    }

    #[test]
    fn test_multiply_2() {
        let p1 = Polynomial::from_vec(vec![-1,1]);
        let p2 = Polynomial::from_vec(vec![0]);
        assert_eq!(p1.multiply(&p2), Polynomial::from_vec(vec![0]));
    }

    #[test]
    fn test_divide_1() {
        let p1 = Polynomial::from_vec(vec![-1,1]);
        let p2 = Polynomial::from_vec(vec![0]);
        println!("{:?}", p1.divide_by(&p2));
        assert_eq!(p1.divide_by(&p2), None);
    }

    // long division with 0 remainder
    #[test]
    fn test_divide_2() {
        let p1 = Polynomial::from_vec(vec![-1,0,0,0,0,0,0,1]);
        let p2 = Polynomial::from_vec(vec![-1,1]);
        assert_eq!(p1.divide_by(&p2).unwrap(), (Polynomial::from_vec(vec![1,1,1,1,1,1,1]), Polynomial::zero()));
    }

    // long division with a nonzero remainder
    #[test]
    fn test_divide_3() {
        let p1 = Polynomial::from_vec(vec![5,-8,0,6]);
        let p2 = Polynomial::from_vec(vec![-4,2]);
        assert_eq!(p1.divide_by(&p2).unwrap(), (Polynomial::from_vec(vec![8,6,3]), Polynomial::const_coef(37, 1)));
    }

    // Integer long division special case treatment 
    #[test]
    fn test_divide_4() {
        let p1 = Polynomial::from_vec(vec![1,2,3]);
        let p2 = Polynomial::from_vec(vec![0,0,2]);
        assert_eq!(p1.divide_by(&p2).unwrap(), (Polynomial::from_vec(vec![1]), Polynomial::from_vec(vec![1,2,1])));
    }

    #[test]
    fn test_fft_1() {
        let p1 = Polynomial::from_vec(vec![-1.,1.]);
        let p2 = Polynomial::from_vec(vec![1.,1.,1.,1.,1.,1.,1.]);
        assert_eq!(p1.multiply(&p2), fft_mul(&p1, &p2));
    }

    #[test]
    fn test_fft_2() {
        let p1 = Polynomial::from_vec(vec![-1.,1.]).pow(4);
        let p2 = Polynomial::from_vec(vec![-1.,1.]); 
        let fft = fft_mul(&p1, &p2);
        assert_eq!(p1.multiply(&p2), fft);
        assert_eq!(p2.pow(5), fft);
    }

    #[test]
    fn pow_1() {
        let p1 = Polynomial::from_vec(vec![-1.,1.]).pow(4);
        let p2 = Polynomial::from_vec(vec![-1.,1.]); 
        assert_eq!(p1.multiply(&p2), p2.pow(5));
    }

    #[test]
    fn pow_2() {
        let p2 = Polynomial::from_vec(vec![-1,1]); 
        assert_eq!(Polynomial::from_vec(vec![-1, 5,-10, 10, -5, 1]), p2.pow(5));
    }



}