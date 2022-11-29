use itertools::{EitherOrBoth::*, Itertools};

struct Polynomial {
    coeffs: Vec<f64>
}

impl Polynomial {
    fn new(mut c:Vec<f64>) -> Polynomial {
        // do this later.
        // if c.is_empty() {
        //     panic!("Polynomial cannot be");
        // }

        // For a polynomial, if we have 0x^5, then we should remove this term, except
        // when 0 is the constant term, which happens when c.len() == 1.        
        while c.len() > 1 {
            let last = *c.last().unwrap();
            if last == 0. {
                c.pop();
            } else {
                break
            }
        }

        Polynomial{coeffs: c}
    }

    fn add(&self, p:&Polynomial) -> Polynomial {
        let mut new_poly:Vec<f64> = Vec::with_capacity(self.coeffs.len().max(p.coeffs.len()));
        for pair in self.coeffs.iter().zip_longest(p.coeffs.iter()) {
            match pair {
                Both(l,r) => new_poly.push(l+r),
                Left(l) => new_poly.push(*l),
                Right(r) => new_poly.push(*r)
            }
        }
        // in the minus case, it is possible that the leading term becomes 0
        // so we use new to deal with that case
        Polynomial::new(new_poly)

    }

    fn minus(&self, p:&Polynomial) -> Polynomial {
        let p2 = Polynomial{coeffs: p.coeffs.iter().map(|x| -x).collect()};
        self.add(&p2)
    }

}

impl std::fmt::Display for Polynomial {
    
    fn fmt(&self, f:&mut std::fmt::Formatter) -> std::fmt::Result {
        let l = self.coeffs.len();
        if l == 0 {
            return write!(f, "Empty Polynomial. (Empty coefficients)")
        }
        let mut p = String::new();
        for (i,coef) in self.coeffs.iter().enumerate().rev() {
            if *coef == 0.0 {
                if i == 0 && p.len() == 0 {
                    p.push('0');
                }    
                continue
            }
            let mut signed = String::new();
            if i == l - 1 { // leading term, don't add sign
                if *coef != 1. {
                    signed.push_str(&coef.to_string());
                }
            } else { // not leading term, add sign and use abs
                if self.coeffs[i] < 0.0 {
                    signed.push_str(&" - ");
                } else {
                    signed.push_str(&" + ");
                }
                if i == 0 || *coef != 1. {
                    // always push constants
                    // No need to push 1 in front of x
                    signed.push_str(&coef.abs().to_string());
                }
            }
            if i > 0 {
                signed.push('x');
            }
            if i > 1 {
                signed.push('^');
                signed.push_str(&i.to_string());
            }
            p.push_str(&signed);
        }
        write!(f, "{}", p)
    }

}


fn main() {
    let p1 = vec![0.5,1.,-2.,-3.];
    let p2 = vec![0.5,-1.,3.,4.];
    //let p2 = vec![0., -1., 2., 3.];
    let q1 = Polynomial {coeffs: p1};
    let q2 = Polynomial{coeffs: p2};
    println!("{}", q1.add(&q2));
}
