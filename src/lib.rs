/*!
 * This crate is a simple implementtion of a Timoshenko Beam in 3D
 * 
 * The crate is built upon the work by Yunhua Luo "An Efficient 3D Timoshenko Beam Element with Consistent Shape Functions"
 * Adv. Theor. Appl. Mech., Vol. 1, 2008, no. 3, 95 - 106
 * The stiffness matrix is computed as per the Appendix A of the paper.

*/

use ndarray::Array2;


/// This is a beam, defined to be solved with the Timoshenko beam theory.
/// 
#[derive(Debug, Clone, PartialEq)]
pub struct Beam {

}

impl Beam {

}

/// This is the elementary element of a [Beam]. 
/// 
/// It discretizes the [Beam] into its singular elements.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct BeamElement {
    id: usize,
    length: f64,
    versor: Array2<f64>,
    properties: Option<BeamProperties>,
    stiffness_matrix: Array2<f64>,
    weight: Option<f64>,
}

impl BeamElement {

    /// This function creates a new [BeamElement] oriented in space
    /// 
    /// The local axis of the [BeamElement] is such that the x axis is the unit vector between the two given points [u1] and [u2], and the y axis is the vector [y], which can also be non-unitary.
    /// If the properties are given as a [BeamProperties] object, the stiffness matrix is computed.
    /// If the weight is computed, this is applied on the z axis.
    fn new(u1: &[f64;3], u2:[f64;3], id: usize, y: &[f64;3], properties: Option<BeamProperties>, weight: Option<f64>) -> Self {
        let mut x = [(u2[0] - u1[0]), (u2[1] - u1[1]), (u2[2] - u1[2])];
        let length = (x[0].powi(2) + x[1].powi(2) + x[2].powi(2)).sqrt();
        x = x*1.0/length;
        let length_y = (y[0].powi(2) + y[1].powi(2) + y[2].powi(2)).sqrt();
        let mut loc_y = y;
        if length_y != 1.0 {
             loc_y = loc_y*1.0/length_y;
        }
        let z = todo!("(a b c)x(d e f) = (bf-ce,cd-af,ae-bd)");
        let versor = array![x, loc_y, z];
        let stiffness_matrix;
        
        if let Some(properties) = &properties {
            stiffness_matrix = properties.compute_stiffness_matrix(length);
            //TODO: add change of base!!!
        } else {
            stiffness_matrix = Array2::zeros([12,12]);
        }
        
        Self {
            id,
            length,
            versor,
            properties,
            stiffness_matrix,
            weight,
        }   
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct BeamProperties {
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
    /// The convention is that x is positive along the beam axis, y is positive upwards and z is positive to the right.
    pub(crate) fn new(e: f64, g: f64, a: f64, iy: f64, iz: f64, shape: BeamShape) -> Self {
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
    fn compute_stiffness_matrix(&self, length: f64) -> Array2<f64> {
        let k = self.shape.get_shape_coefficient();
        let mut  stiffness_matrix = Array2::zeros([12,12]);
        let ea = self.e*self.a;
        let eiy = self.e*self.iy;
        let eiz = self.e*self.iz;
        let eiz12 = 12.0*eiz;
        let eiy12 = 12.0*eiy;
        let ga = self.g*self.a;
        let kgal2 = k*ga*(length.powi(2));
        let k11 = ea/length;
        let k22 = k*ga*eiy12*(eiy12+kgal2)/(length*(eiy12-kgal2).powi(2));
        let k26 = k22*length*0.5;
        let k33 = k*ga*eiz12*(eiz12+kgal2)/(length*(eiz12-kgal2).powi(2));
        let k35 = -k33*length*0.5;
        let k44 = self.g*(self.iy+self.iz)/length;
        let k55 = (4.0*eiz*(kgal2.powi(2)+3.0*kgal2*eiz+36.0*eiz.powi(2)))/(length*(eiz12-kgal2).powi(2));
        let k59 = (6.0*k*ga*eiz*(eiz12+kgal2))/((eiz12-kgal2).powi(2));
        let k511 = -(2.0*eiz*(72.0*(eiz.powi(2))-kgal2.powi(2)-30.0*kgal2*eiz))/(length*(eiz12-kgal2).powi(2));
        let k66 = (4.0*eiy*(kgal2.powi(2)+3.0*kgal2*eiy+36.0*eiy.powi(2)))/(length*(eiy12-kgal2).powi(2));
        let k68 = -(6.0*k*ga*eiy*(eiy12+kgal2))/((eiy12-kgal2).powi(2));
        let k612 = -(2.0*eiy*(-kgal2.powi(2)-30.0*kgal2*eiy+72.0*(eiy.powi(2))))/(length*(eiy12-kgal2).powi(2));
        let k88 = -k68*2.0/length;
        let k812 = k68;
        let k99 = k59*2.0/length;
        let k911 = k59;
        let k1111 = k55;
        let k1212 = k66;

        stiffness_matrix[[0,0]] =  k11;
        stiffness_matrix[[6,6]]= k11;
        stiffness_matrix[[0,6]] = -k11;
        stiffness_matrix[[1,1]] = k22;
        stiffness_matrix[[1,7]] = -k22;
        stiffness_matrix[[1,5]] = k26;
        stiffness_matrix[[1,11]] = k26;
        stiffness_matrix[[2,2]] = k33;
        stiffness_matrix[[2,8]]= k33;
        stiffness_matrix[[2,4]] = k35;
        stiffness_matrix[[2,10]] = k35;
        stiffness_matrix[[3,3]] = k44;
        stiffness_matrix[[6,6]] = k44;
        stiffness_matrix[[3,9]] = -k44;
        stiffness_matrix[[4,4]] = k55;
        stiffness_matrix[[4,8]] = k59;
        stiffness_matrix[[4,10]] = k511;
        stiffness_matrix[[5,5]] = k66;
        stiffness_matrix[[5,7]] = k68;
        stiffness_matrix[[5,11]] = k612;
        stiffness_matrix[[7,7]] = k88;
        stiffness_matrix[[7,11]] = k812;
        stiffness_matrix[[8,8]] = k99;
        stiffness_matrix[[8,10]] = k911;
        stiffness_matrix[[10,10]] = k1111;
        stiffness_matrix[[11,11]] = k1212;

        stiffness_matrix
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum BeamShape {
    Rectangular,
}

const RECTANGULAR_SHAPE_COEFF: f64 = 5.0/4.0;

impl BeamShape {
    fn get_shape_coefficient(&self) -> f64 {
        match self {
            BeamShape::Rectangular => RECTANGULAR_SHAPE_COEFF,
        }
    }
}