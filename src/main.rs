use itertools::{EitherOrBoth::*, Itertools};
use rand::distributions::{Distribution, Uniform};
struct Polynomial {
    coeffs: Vec<f64>
}

impl Polynomial {
    /// Owns c
    fn new(mut c:Vec<f64>) -> Polynomial {
        // do this later.
        // if c.is_empty() {
        //     panic!("Polynomial cannot be");
        // }

        // For a polynomial, if we have  a coefficient which is 0, like in 0x^5, 
        // then we should remove this term, except when 0 is the constant term.        
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

    /// Polynomial of degree n with constant coefficient c
    fn const_coef(c:f64, n:usize) -> Polynomial {
        match n {
            0 => Polynomial{coeffs: vec![0.; 1]},
            _ => Polynomial{coeffs: vec![c; n]}
        }
    }

    /// Creates the n basis polynomial with coefficient c.
    /// E.g. n = 0, c = 1, returns 1 
    /// n = 2, c = 2 returns 2x^2
    fn basis(c:f64, n:usize) -> Polynomial {
        match n {
            0 => Polynomial{coeffs: vec![1.; 1]},
            _ => {
                let mut v = vec![0.; n-1];
                v.push(c);
                Polynomial{coeffs: v}
            }
        }
    }

    /// Randomly generate a polynomial of deg n
    /// with uniformly distributed coefficients between lower and upper
    fn uniform_rand(lower:f64, upper:f64, n:usize) -> Polynomial {
        match n {
            0 => Polynomial{coeffs: vec![0.; 1]},
            _ => {
                let mut rng = rand::thread_rng();
                let u = Uniform::new(lower, upper);
                Polynomial::new((0..n).map(|_| u.sample(&mut rng)).collect()) // Well, maybe leading coeff is 0..
            }
        }
    }

    // Arithmetic

    fn add(&self, p:&Polynomial) -> Polynomial {
        let mut new_poly:Vec<f64> = Vec::with_capacity(self.coeffs.len().max(p.coeffs.len()));
        for pair in self.coeffs.iter().zip_longest(p.coeffs.iter()) {
            match pair {
                Both(l,r) => new_poly.push(l+r),
                Left(l) => new_poly.push(*l),
                Right(r) => new_poly.push(*r)
            }
        }
        // it is possible that the leading term becomes 0 after adding/subtracting
        // so we use new to deal with that case
        Polynomial::new(new_poly)

    }

    fn minus(&self, p:&Polynomial) -> Polynomial {
        let p2 = Polynomial{coeffs: p.coeffs.iter().map(|x| -x).collect()};
        self.add(&p2)
    }

    fn multiply(&self, p:&Polynomial) -> Polynomial {
        let mut new_poly:Vec<f64> = vec![0.; self.coeffs.len() + p.coeffs.len() - 1];
        for (i,a) in self.coeffs.iter().enumerate() {
            for (j,b) in p.coeffs.iter().enumerate() {
                new_poly[i + j] += a*b;
            }
        }
        Polynomial {coeffs: new_poly}
    }

    ///
    /// Long Division by Euclidean Algorithm
    /// 
    fn divide_by(&self, p:&Polynomial) -> (Polynomial, Polynomial) {
        // let mut p1 = Polynomial{coeffs:self.coeffs.clone()};
        // let mut p2 = Polynomial{coeffs:p.coeffs.clone()};
        // let mut diff = p1.deg() - p2.deg();
        // while diff >= 0 {
        //     let coef = p1.coeffs.last().unwrap() / p2.coeffs.last().unwrap();
        //     let multiplier = Polynomial::basis(coef, diff);
        //     let new_p = p.minus(& q.multiply(&Polynomial::basis(coef, d)));
        //     Polynomial::long_div(new_p, q);
        //     todo!()
        // }
        todo!()
    }

    fn long_div(p:Polynomial, q:&Polynomial) -> (Polynomial, Polynomial) {
        let mut output_q = vec![0.; 1];
        let p_deg = p.deg();
        let q_deg = q.deg();
        if p_deg < q_deg {
            (Polynomial{coeffs: vec![0.;1]}, Polynomial{coeffs: q.coeffs})
        } else { // p_deg >= q_deg

            let d = p_deg - q_deg;
            let coef = p.coeffs.last().unwrap() / q.coeffs.last().unwrap();
            let multiplier = Polynomial::basis(coef, d);
            let new_p = p.minus(& q.multiply(&Polynomial::basis(coef, d)));
            Polynomial::long_div(new_p, q);
            todo!()
        }

    }

    // Utility

    fn eval(&self, x:f64) -> f64 {
        if x == 0. {
            self.coeffs[0]
        } else {
            self.coeffs.iter()
            .enumerate()
            .fold(0., |acc, (idx, coef)| acc + coef * (x.powi(idx as i32)))
        }
    }

    fn deg(&self) -> usize {
        let n = self.coeffs.len();
        match n {
            0|1 => 0,
            _ => n-1 // complier knows this is > 0 and therefore always a usize!!!
        }
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
    let q3 = q1.add(&q2);
    println!("{}", q3);
    println!("Evaluating at 1 is: {}", q3.eval(1.));

    let p3 = vec![0.,1.,1.];
    let p4 = vec![0., 2.];
    let q3 = Polynomial{coeffs: p3};
    let q4 = Polynomial{coeffs: p4};
    println!("Multiply {} by {} is:\n{}", q3, q4, q3.multiply(&q4));

}
