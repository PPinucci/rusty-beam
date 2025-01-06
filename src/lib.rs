/*!
    This crate is a simple implementtion of a Timoshenko Beam in 3D
    
    The crate is built upon the work by Yunhua Luo "An Efficient 3D Timoshenko Beam Element with Consistent Shape Functions"
    Adv. Theor. Appl. Mech., Vol. 1, 2008, no. 3, 95 - 106
    
    The stiffness matrix is computed as per the Appendix A of the paper.

*/

pub mod beam;

/// Cross product of 2 3D vectors
/// (a b c)x(d e f) = (bf-ce,cd-af,ae-bd)
pub fn cross_prodct(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Scalar multiplication of a 3D vector
/// (a b c)*d = (ad,bd,cd)
pub fn scalar_mult(a: &[f64; 3], b: f64) -> [f64; 3] {
    [a[0] * b, a[1] * b, a[2] * b]
}

/// Modulus of a 3D vector
/// |a b c| = sqrt(a^2 + b^2 + c^2)
pub fn modulus(a: &[f64; 3]) -> f64 {
    (a[0].powi(2) + a[1].powi(2) + a[2].powi(2)).sqrt()
}
