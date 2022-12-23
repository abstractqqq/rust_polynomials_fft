mod polynomial;
use polynomial::Polynomial;

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


    let p1 = Polynomial::new(vec![-1., 0., 0., 0., 0., 0., 0., 1.]);
    let p2 = Polynomial::new(vec![-1., 1.]);
    let (quotient, remainder) = p1.divide_by(&p2);

    println!("The quotient is {:?}. \nThe remainder is {:?}.", quotient, remainder);


    // let p1 = Polynomial::new(vec![-1., 1.]);
    // let n = 7;
    // let p1_power = p1.to_power(n);

    // println!("Raising ({:?}) to the {:?}th power yields\n{:?}.", p1, n, p1_power);

}
