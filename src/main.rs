use itertools::{EitherOrBoth::*, Itertools};
use rand::distributions::{Distribution, Uniform};
struct Polynomial {
    coeffs: Vec<f64>
}

impl Polynomial {
    fn new(mut c:Vec<f64>) -> Polynomial {
        // do this later.
        // if c.is_empty() {
        //     panic!("Polynomial cannot be");
        // }

        // For a polynomial, if we have  a coefficient which is 0, like in 0x^5, 
        // then we should remove this term, except when 0 is the constant term.        
        while c.len() > 1 {
            let last = c.last().unwrap();
            if *last == 0. {
                c.pop();
            } else {
                break
            }
        }
        Polynomial{coeffs: c}
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
            0 => Polynomial{coeffs: vec![c; 1]},
            _ => {
                let mut v = vec![0.; n];
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

    fn highest_coeff(&self) -> &f64 {
        self.coeffs.last().unwrap()
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

    /// Long Division
    fn divide_by(&self, p:&Polynomial) -> (Polynomial, Polynomial) {
        // (P1, P2) = (Quotient, Remainder)
        let dividee = Polynomial{coeffs: self.coeffs.clone()};
        let dividee_deg = dividee.deg();
        let divider_deg = p.deg();
        if dividee_deg < divider_deg {
            return (Polynomial::const_coef(0., 1), Polynomial{coeffs: p.coeffs.clone()})
        }
        let mut quotient = vec![0.; dividee_deg - divider_deg + 1];
        let remainder = Polynomial::_long_div(dividee, p, &mut quotient);
        (Polynomial::new(quotient), remainder)
    }

    fn _long_div(
        dividee:Polynomial
        , divider:&Polynomial
        , quotient:&mut Vec<f64>
    ) -> Polynomial { // returns the remainder    
        let dividee_deg = dividee.deg();
        let divider_deg = divider.deg();
        if dividee_deg < divider_deg {// dividee becomes remainder
            return  dividee
        } else {
            // dividee_deg >= divider_deg.
            let new_term_deg = dividee_deg - divider_deg; // always >= 0
            let new_term_coeff = dividee.highest_coeff() / divider.highest_coeff();
            // modify the quotient
            quotient[new_term_deg] = new_term_coeff;
            // update dividee
            let basis = Polynomial::basis(new_term_coeff, new_term_deg);
            let new_dividee = dividee.minus(&basis.multiply(divider));
            println!("Division steps: dividing {} by {}", new_dividee, divider);
            Self::_long_div(new_dividee, divider, quotient)
        }

    }

    // raise a polynomial p to a deg.
    fn to_power(&self, n:usize) -> Polynomial {
        match n {
            0 => {
                    println!("Raising a polynomial to the 0th power always returns the polynomial 1. This may not be the mathematical convention, nor the definition you choose to use.");
                    Polynomial { coeffs: vec![1.;1] }
                },
            _ => {
                let cur = Polynomial {coeffs : self.coeffs.clone()};
                Self::_power(cur, n)
            }
        }
    }

    fn _power(current:Polynomial, deg:usize) -> Polynomial {
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

impl std::fmt::Display for Polynomial {
    
    fn fmt(&self, f:&mut std::fmt::Formatter) -> std::fmt::Result {
        let l = self.coeffs.len();
        if l == 0 {
            return write!(f, "Empty Polynomial. (Empty coefficients)")
        }
        let mut p = String::new();
        for (i,coef) in self.coeffs.iter().enumerate().rev() {
            if *coef == 0.0 { // if we have 0 constant, then only print it when the polynomial has no higher terms.
                if i == 0 && p.len() == 0 {
                    p.push('0');
                }    
                continue
            }
            let mut signed = String::new();
            if i == l - 1 { 
                // leading term, no need to manually add signs
                if self.deg() == 0 {
                    // push everything if poly is constant
                    signed.push_str(&coef.to_string());                   
                } else { // poly is not constant and must have a x^k term, k > 0, don't print 1.
                    if *coef != 1. {
                        signed.push_str(&coef.to_string());
                    }
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

    // println!("The quotient is {}. \nThe remainder is {}.", quotient, remainder);


    let p1 = Polynomial::new(vec![-1., 1.]);
    let n = 7;
    let p1_power = p1.to_power(n);

    println!("Raising ({}) to the {}th power yields\n{}.", p1, n, p1_power);

}
