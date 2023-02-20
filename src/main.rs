mod polynomial;
use polynomial::Polynomial;
use polars::prelude::*;
use std::time::Instant;

fn run_performance_test(runs:usize) -> Result<DataFrame, PolarsError> {
    let mut regular_mul:Vec<f64> = Vec::with_capacity(runs);
    let mut fft_mul:Vec<f64> = Vec::with_capacity(runs);
    let mut fft_mul_threaded:Vec<f64> = Vec::with_capacity(runs);

    for i in 10..(runs+10) {
        let p1 = Polynomial::const_coef(3., i);
        let p2 = Polynomial::const_coef(5., i);
        let now = Instant::now();
        let _result = p1.multiply(&p2);
        regular_mul.push(now.elapsed().as_secs_f64())
    }

    for i in 10..(runs+10) {
        let p1 = Polynomial::const_coef(3., i);
        let p2 = Polynomial::const_coef(5., i);
        let now = Instant::now();
        let _result = p1.fft_mul(&p2, 10);
        fft_mul.push(now.elapsed().as_secs_f64())
    }

    for i in 10..(runs+10) {
        let p1 = Polynomial::const_coef(3., i);
        let p2 = Polynomial::const_coef(5., i);
        let now = Instant::now();
        let _result = p1.fft_mul_threaded(&p2, 10);
        fft_mul_threaded.push(now.elapsed().as_secs_f64())
    }

    let regular = Series::from_vec("regular", regular_mul);
    let fft = Series::from_vec("fft", fft_mul);
    let fft_threaded = Series::from_vec("fft_threaded", fft_mul_threaded);

    DataFrame::new(vec![regular, fft, fft_threaded])
}

fn main() {
    // // Trying to see when will FFT be faster than regular multiplication.
    // let args: Vec<String> = env::args().collect();
    // let n = args[1].parse::<usize>().unwrap();
    // let result = run_performance_test(n);
    // match result {
    //     Ok(mut df) => {
    //         use std::fs::File;
    //         println!("Finished measuring runtime. Saving...");
    //         let mut file = File::create("test_results.csv").expect("Could not create file");
    //         let write_result = CsvWriter::new(&mut file)
    //                                                     .has_header(true)
    //                                                     .with_delimiter(b',')
    //                                                     .finish(&mut df);
    //         if write_result.is_err() {
    //             println!("Error happened when writing to csv.");
    //         }
    //     }

    //     Err(e) => {
    //         println!("Some error occured during the test: {}", e);
    //     }
    // }


    let p1 = Polynomial::from_vec(vec![1,2,3]);
    let p2 = Polynomial::from_vec(vec![0,0,2]);
    let (quotient, remainder) = p1.divide_by(&p2).unwrap();
    println!("Dividing {} by {}.", p1, p2);
    println!("The quotient is {},",quotient);
    println!("The remainder is {}.", remainder);
    println!("Hence {} = ({}) * ({}) + {}.\n\n", p1, quotient, p2, remainder);

    
    let q1 = Polynomial::from_vec(vec![-1.,1.]);
    let q2 = Polynomial::from_vec(vec![1.,1.,1.,1.,1.,1.,1.]);
    println!("Multiplying {} by {}.", q1, q2);
    println!("The result is {}.\n",q1.fft_mul(&q2, 5)); // compute product using FFT, and keep 5 decimal places.

    let q3 = q1.multiply(&q2);
    let (quotient2, remainder2) = q3.divide_by(&q1).unwrap();
    println!("Dividing {} by {}.", q3, q1);
    println!("The quotient is {},",quotient2);
    println!("The remainder is {}.", remainder2);
    println!("Hence {} = ({}) * ({}) + {}.\n\n", q3, quotient2, q1, remainder2);
    
    let r1 = Polynomial::from_vec(vec![1,2,3,4,5]); 
    println!("The derivative of {} is:", r1);
    println!("{}", r1.ddx());

    let r2 = Polynomial::from_vec(vec![1.,1.,1.,1.,1.]); 
    println!("The derivative of {} is:", r2);
    println!("{}", r2.ddx());

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
        assert_eq!(p1.divide_by(&p2).unwrap(), (Polynomial::from_vec(vec![1,1,1,1,1,1,1]), Polynomial::from_vec(vec![0])));
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
        assert_eq!(p1.multiply(&p2), p1.fft_mul(&p2, 5));
    }

    #[test]
    fn test_fft_2() {
        let p1 = Polynomial::from_vec(vec![-1.,1.]).pow(4);
        let p2 = Polynomial::from_vec(vec![-1.,1.]); 
        let fft = p1.fft_mul(&p2, 10);
        assert_eq!(p1.multiply(&p2), fft);
        assert_eq!(p2.pow(5), fft);
    }

    #[test]
    fn test_fft_3() {
        let p1 = Polynomial::from_vec(vec![-1.,1.]).pow(4);
        let p2 = Polynomial::from_vec(vec![-1.,1.]); 
        let fft = p1.fft_mul_threaded(&p2, 10);
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
        let p1 = Polynomial::from_vec(vec![-1,1]); 
        assert_eq!(Polynomial::from_vec(vec![-1, 5,-10, 10, -5, 1]), p1.pow(5));
    }

    #[test]
    fn eval_1() {
        let p1 = Polynomial::from_vec(vec![-1,0,0,0,0,0,0,1]); 
        assert_eq!(p1.eval(0), -1);
    }

    #[test]
    fn eval_2() {
        let p1 = Polynomial::from_vec(vec![1,1]).pow(5); 
        assert_eq!(p1.eval(1), 32);
    }

    #[test]
    fn ddx_1() {
        let p1 = Polynomial::from_vec(vec![0]); 
        assert_eq!(p1.ddx(), Polynomial::from_vec(vec![0]));
    }

    #[test]
    fn ddx_2() {
        let p1 = Polynomial::from_vec(vec![1,2,3,4,5]); 
        assert_eq!(p1.ddx(), Polynomial::from_vec(vec![2,6,12,20]));
    }

}