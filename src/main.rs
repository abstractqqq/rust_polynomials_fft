
struct Polynomial {
    coeffs: Vec<f64>,
    deg: usize
}

impl Polynomial {
    fn new(mut c:Vec<f64>) -> Polynomial {

        loop {
            match c.last(){
                Some(v) => {
                    if v == &0. {
                        c.pop();
                    } else {
                        break
                    }
                }
                _ => break
            }
        }
        let n = (c.len() - 1).max(0); // if c is empty, return 0
        Polynomial {coeffs: c, deg: n}
    }

}

impl std::fmt::Display for Polynomial {
    
    fn fmt(&self, f:&mut std::fmt::Formatter) -> std::fmt::Result {
        let l = self.coeffs.len();
        let mut p = String::new();
        for i in (0..l).rev() {
            if self.coeffs[i] == 0.0 {
                continue;
            }
            let mut signed = String::new();
            if self.coeffs[i] < 0.0 {
                signed.push_str(&" - ");
            } else {
                if i != l - 1 {
                    signed.push_str(&" + ");
                }
            }
            signed.push_str(&self.coeffs[i].to_string());
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
    let p1 = vec![0.5,1.,2.,3.];
    let p = Polynomial { coeffs: p1, deg: 3};
    println!("{}", p);
}
