use itertools::{EitherOrBoth::*, Itertools};
use ndarray::{Array1, ArrayView1};
use num_traits::Num;

#[derive(Debug)]
pub struct Polynomial<T> 
where T: Num + Clone + Copy
{
    coeffs: Array1<T>
}

impl <T> Polynomial<T> 
    where T: Num + Clone + Copy
{
    pub fn new(c: Array1<T>) -> Polynomial<T> {
        Polynomial {coeffs: c}
    }

    pub fn from_vec(c: Vec<T>) -> Polynomial<T> {
        Polynomial {coeffs: Array1::from_vec(c)}
    }

    pub fn no_leading_zeros(mut c:Vec<T>) -> Polynomial<T> {
        // For a polynomial, if we have  a coefficient which is 0, like in 0x^5, 
        // then we should remove this term, except when 0 is the constant term.        
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

    pub fn copy(&self) -> Polynomial<T> {
        Polynomial{coeffs: self.coeffs.clone()}
    }

    // Utility
    pub fn eval(&self, x:T) -> T {
        if x == T::zero() {
            self.coeffs[0]
        } else {
            let mut x_pow = T::zero();
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

    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    pub fn get_coeffs(&self) -> Array1<T> {
        self.coeffs.clone()
    }

    pub fn get_coeffs_view(&self) -> ArrayView1<T> {
        self.coeffs.view()
    }

    pub fn get_const(&self) -> &T {
        self.coeffs.get(0).unwrap()
    }

    /// Polynomial of degree n with constant coefficient c
    pub fn const_coef(c:T, n:usize) -> Polynomial<T> {
        match n {
            0 => Polynomial{coeffs: Array1::from_elem(1, T::zero())},
            _ => Polynomial{coeffs: Array1::from_elem(n, c)}
        }
    }

    /// Creates the n basis polynomial with coefficient c.
    /// E.g. n = 0, c = 1, returns 1 
    /// n = 2, c = 2 returns 2x^2
    pub fn basis(c:T, n:usize) -> Polynomial<T> {
        match n {
            0 => Polynomial{coeffs: Array1::from_elem(1, c)},
            _ => {
                let mut v = vec![T::zero(); n];
                v.push(c);
                Polynomial{coeffs: Array1::from_vec(v)}
            }
        }
    }

    pub fn highest_coeff(&self) -> T {
        *self.coeffs.last().unwrap()
    }

    // Arithmetic
    pub fn add(&self, p:&Polynomial<T>) -> Polynomial<T> {
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

    pub fn minus(&self, p:&Polynomial<T>) -> Polynomial<T> {
        let p2 = Polynomial{coeffs: p.coeffs.iter().map(|x| T::zero()-*x).collect()};
        self.add(&p2)
    }

    pub fn multiply(&self, p:&Polynomial<T>) -> Polynomial<T> {
        let mut new_poly:Vec<T> = vec![T::zero(); self.coeffs.len() + p.coeffs.len() - 1];
        for (i,a) in self.coeffs.iter().enumerate() {
            for (j,b) in p.coeffs.iter().enumerate() {
                let k = i + j;
                new_poly[k] = new_poly[k] + (*a)*(*b);
            }
        }
        Polynomial {coeffs: Array1::from_vec(new_poly)}
    }

    /// Long Division
    pub fn divide_by(&self, p:&Polynomial<T>) -> (Polynomial<T>, Polynomial<T>) {
        // (P1, P2) = (Quotient, Remainder)
        let dividee = self.copy();
        let dividee_deg = dividee.deg();
        let divider_deg = p.deg();
        if dividee_deg < divider_deg {
            return (Polynomial::const_coef(T::zero(), 1), Polynomial{coeffs: p.coeffs.clone()})
        }
        let mut quotient = vec![T::zero(); dividee_deg - divider_deg + 1];
        let remainder = Polynomial::_long_div(dividee, p, &mut quotient);
        (Polynomial{coeffs: Array1::from_vec(quotient)}, remainder)
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

    // raise a polynomial p to a deg.
    pub fn to_power(&self, n:usize) -> Polynomial<T> {
        match n {
            0 => {
                    println!("DON'T DO THIS.");
                    Polynomial { coeffs: Array1::from_elem(1, T::one()) }
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


// This is pretty print for Polynomials defined over f64. Hard to generalize this. So I will put this in comments.
// impl std::fmt::Display for Polynomial

// {
    
//     fn fmt(&self, f:&mut std::fmt::Formatter) -> std::fmt::Result {
//         let l = self.coeffs.len();
//         if l == 0 {
//             return write!(f, "Empty Polynomial. (Empty coefficients)")
//         }
//         let mut p = String::new();
//         for (i,coef) in self.coeffs.iter().enumerate().rev() {
//             if *coef == f32::default() { // if we have 0 constant, then only print it when the polynomial has no higher terms.
//                 if i == 0 && p.len() == 0 {
//                     p.push;
//                 }    
//                 continue
//             }
//             let mut signed = String::new();
//             if i == l - 1 { 
//                 // leading term, no need to manually add signs
//                 if self.deg() == 0 {
//                     // push everything if poly is constant
//                     signed.push_str(&coef.to_string());                   
//                 } else { // poly is not constant and must have a x^k term, k > 0, don't print 1.
//                     if *coef != 1. {
//                         signed.push_str(&coef.to_string());
//                     }
//                 }
//             } else { // not leading term, add sign and use abs
//                 if self.coeffs[i] < 0.0 {
//                     signed.push_str(&" - ");
//                 } else {
//                     signed.push_str(&" + ");
//                 }
//                 if i == 0 || *coef != 1. {
//                     // always push constants
//                     // No need to push 1 in front of x
//                     signed.push_str(&coef.abs().to_string());
//                 }
//             }
//             if i > 0 {
//                 signed.push('x');
//             }
//             if i > 1 {
//                 signed.push('^');
//                 signed.push_str(&i.to_string());
//             }
//             p.push_str(&signed);
//         }
//         write!(f, "{}", p)
//     }

// }