/*!
    Contains the structural properties of a [BeamElement], the [BeamProperties] struct, along with the [BeamShape] enum.

    The [BeamProperties] struct contains all the necessry properties to compute the stiffness matrix of a beam element.

    The [BeamShape] enum contains the different shapes of beams that can be used to compute the shape factor.

*/

use ndarray::Array2;

/// This struct contains the properties of a beam element.
///
/// The properties are:
/// - the Young's modulus,
/// - the shear modulus,
/// - the cross-sectional area,
/// - the moment of inertia about the y-axis
/// - the moment of inertia about the z-axis.
/// - the shape of the beam which dictates the shape coefficient.
#[derive(Debug, Clone, PartialEq)]
pub struct BeamProperties {
    e: f64,
    g: f64,
    a: f64,
    iy: f64,
    iz: f64,
    shape: BeamShape,
}

impl BeamProperties {
    /// This function creates a new [BeamProperties] object with the given properties.
    ///
    /// The convention is that x is positive along the beam axis, y is where the y versor is defined in the [BeamElement] struct
    /// and z is the right-hand rule cross product of the two versors.
    /// 
    /// # Arguments
    /// - `e` - Young's modulus
    /// - `g` - Shear modulus
    /// - `a` - Cross-sectional area
    /// - `iy` - Moment of inertia about the y-axis
    /// - `iz` - Moment of inertia about the z-axis
    /// - `shape` - Shape of the beam (from the [BeamShape] enum)
    /// 
    pub fn new(e: f64, g: f64, a: f64, iy: f64, iz: f64, shape: BeamShape) -> Self {
        Self {
            e,
            g,
            a,
            iy,
            iz,
            shape,
        }
    }

    /// This function computes the stiffness matrix of the beam element.
    ///
    /// It is computed according to Appendiz A of the paper by Yunhua Luo.
    /// It is a 12x12 matrix, calculated for a beam of length `length` which is uniform and symmetric.
    /// In the loal coordinate it is [K] such that it soves the equation:
    /// [u1, u2, v1, v2, w1, w2, phi_x1, phi_x2, phi_y1, phi_y2, phi_z1, phi_z2] = [K] * [F]
    /// where [u, v, w] are the displacements along the x, y and z axis, and [phi_x, phi_y, phi_z] are the rotations about the x, y and z axis
    /// and the appendix marks the location of the node as 1 and 2 (start and end of the beam).
    /// [F] is the force vector.
    /// 
    pub(crate) fn compute_stiffness_matrix(&self, length: f64) -> Array2<f64> {
        let k = self.shape.get_shape_coefficient();
        let mut stiffness_matrix = Array2::zeros([12, 12]);
        let ea = self.e * self.a;
        let eiy = self.e * self.iy;
        let eiz = self.e * self.iz;
        let eiz12 = 12.0 * eiz;
        let eiy12 = 12.0 * eiy;
        let ga = self.g * self.a;
        let kgal2 = k * ga * (length.powi(2));
        let k11 = ea / length;
        let k22 = k * ga * eiy12 * (eiy12 + kgal2) / (length * (eiy12 - kgal2).powi(2));
        let k26 = k22 * length * 0.5;
        let k33 = k * ga * eiz12 * (eiz12 + kgal2) / (length * (eiz12 - kgal2).powi(2));
        let k35 = -k33 * length * 0.5;
        let k44 = self.g * (self.iy + self.iz) / length;
        let k55 = (4.0 * eiz * (kgal2.powi(2) + 3.0 * kgal2 * eiz + 36.0 * eiz.powi(2)))
            / (length * (eiz12 - kgal2).powi(2));
        let k59 = (6.0 * k * ga * eiz * (eiz12 + kgal2)) / ((eiz12 - kgal2).powi(2));
        let k511 = -(2.0 * eiz * (72.0 * (eiz.powi(2)) - kgal2.powi(2) - 30.0 * kgal2 * eiz))
            / (length * (eiz12 - kgal2).powi(2));
        let k66 = (4.0 * eiy * (kgal2.powi(2) + 3.0 * kgal2 * eiy + 36.0 * eiy.powi(2)))
            / (length * (eiy12 - kgal2).powi(2));
        let k68 = -(6.0 * k * ga * eiy * (eiy12 + kgal2)) / ((eiy12 - kgal2).powi(2));
        let k612 = -(2.0 * eiy * (-kgal2.powi(2) - 30.0 * kgal2 * eiy + 72.0 * (eiy.powi(2))))
            / (length * (eiy12 - kgal2).powi(2));
        let k88 = -k68 * 2.0 / length;
        let k812 = k68;
        let k99 = k59 * 2.0 / length;
        let k911 = k59;
        let k1111 = k55;
        let k1212 = k66;

        stiffness_matrix[[0, 0]] = k11;
        stiffness_matrix[[6, 6]] = k11;
        stiffness_matrix[[0, 6]] = -k11;
        stiffness_matrix[[1, 1]] = k22;
        stiffness_matrix[[1, 7]] = -k22;
        stiffness_matrix[[1, 5]] = k26;
        stiffness_matrix[[1, 11]] = k26;
        stiffness_matrix[[2, 2]] = k33;
        stiffness_matrix[[2, 8]] = k33;
        stiffness_matrix[[2, 4]] = k35;
        stiffness_matrix[[2, 10]] = k35;
        stiffness_matrix[[3, 3]] = k44;
        stiffness_matrix[[6, 6]] = k44;
        stiffness_matrix[[3, 9]] = -k44;
        stiffness_matrix[[4, 4]] = k55;
        stiffness_matrix[[4, 8]] = k59;
        stiffness_matrix[[4, 10]] = k511;
        stiffness_matrix[[5, 5]] = k66;
        stiffness_matrix[[5, 7]] = k68;
        stiffness_matrix[[5, 11]] = k612;
        stiffness_matrix[[7, 7]] = k88;
        stiffness_matrix[[7, 11]] = k812;
        stiffness_matrix[[8, 8]] = k99;
        stiffness_matrix[[8, 10]] = k911;
        stiffness_matrix[[10, 10]] = k1111;
        stiffness_matrix[[11, 11]] = k1212;

        stiffness_matrix
    }
}

/// This enum contains the different shapes of beams that can be used to compute the shape factor.
///
/// This is non-exaustive, meaning that it can be extended in the future.
///
/// A [BeamShape::Custom] can be used to define a custom shape factor.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum BeamShape {
    /// Rectangular beam, shape coefficient is 5/4
    Rectangular,
    /// Custom beam shape, shape coefficient value is given by the user
    Custom(f64),
}

const RECTANGULAR_SHAPE_COEFF: f64 = 5.0 / 4.0;

impl BeamShape {
    /// This function returns the shape coefficient of the beam.
    ///
    ///
    fn get_shape_coefficient(&self) -> f64 {
        match self {
            BeamShape::Rectangular => RECTANGULAR_SHAPE_COEFF,
            BeamShape::Custom(value) => *value,
            _ => RECTANGULAR_SHAPE_COEFF,
        }
    }
}
