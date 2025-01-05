/*
    This module contains the elementary element of the beam, the BeamElement struct, and its properties, the BeamProperties struct.

    The BeamElement struct is the elementary element of the beam, and it is defined by two points in space, 
    the id of the element, the versor of the element, the properties of the element and optionally by the stiffness matrix of the element and the weight of the element.
*/

use ndarray::{array, Array2};
use crate::{beam::Beam, cross_prodct, modulus, scalar_mult};

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
    /// 
    /// # Arguments
    /// 
    /// * `u1` - A reference to a 3D array representing the first point of the element.
    /// * `u2` - A 3D array representing the second point of the element.
    /// * `id` - The id of the element.
    /// * `y` - A reference to a 3D array representing the y axis of the element (where the local Iyy is defined).
    /// * `properties` - An optional [BeamProperties] object defining the structural properties of the element.
    /// * `weight` - An optional f64 defining the weight of the element, always pointing in the negative global z axis.
    /// 
    /// # Example
    /// 
    /// '''
    /// use timoshenko_beam::beam::BeamElement;
    /// use timoshenko_beam::beam::BeamProperties;
    /// use timoshenko_beam::beam::BeamShape;
    /// let properties = BeamProperties::new(1.0, 1.0, 1.0, 1.0, 1.0, BeamShape::Rectangular);
    ///
    /// let beam_element = BeamElement::new(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 0, &[0.0, 1.0, 0.0], Some(properties), None);
    /// '''
    fn new(u1: &[f64;3], u2:[f64;3], id: usize, y: &[f64;3], properties: Option<BeamProperties>, weight: Option<f64>) -> Self {
        // compute the axis of the beam
        let mut x = [(u2[0] - u1[0]), (u2[1] - u1[1]), (u2[2] - u1[2])];
        // calculate the length
        let length = modulus(&x);
        // normalize the x axis
        x = scalar_mult(&x, 1.0/length);
        // calculate the modulus of y
        let length_y = modulus(&y);
        // copy y in order to normalize without declaring the input y to be mutable
        let mut loc_y = *y;
        // we only normalyze if the modulus is different from 1
        if length_y != 1.0 {
             loc_y = scalar_mult(&loc_y,1.0/length_y);
        }
        // calculate the z axis as the cross product of x and y
        //  since both are normalized oleady, this is normalized as well
        let z = cross_prodct(&x, &loc_y);
        // create the versor matrix
        let versor = array![x, loc_y, z];

        // allocate the stiffness matrix
        let stiffness_matrix;
        if let Some(properties) = &properties {
            // the stifness matrix in the local reference is computed if we have properties input
            stiffness_matrix = properties.compute_stiffness_matrix(length);
            //TODO: add change of base to the global reference!!!

        } else {
            // to avoid an Option, we just allocate a zero matrix
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

    /// Tests if properties are defined for the [BeamElement]
    /// 
    /// Returns true if the properties are defined, false otherwise.
    /// 
    /// # Example
    /// '''
    /// use timoshenko_beam::beam::BeamElement;
    /// use timoshenko_beam::beam::BeamProperties;
    /// use timoshenko_beam::beam::BeamShape;
    /// let properties = BeamProperties::new(1.0, 1.0, 1.0, 1.0, 1.0, BeamShape::Rectangular);
    /// 
    /// let beam_element = BeamElement::new(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 0, &[0.0, 1.0, 0.0], Some(properties), None);
    /// 
    /// assert_eq!(beam_element.has_properties(), true);
    /// 
    /// let beam_element = BeamElement::new(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 0, &[0.0, 1.0, 0.0], None, None);
    /// 
    /// 
    /// assert_eq!(beam_element.has_properties(), false);
    /// '''
    pub(crate) fn has_properties(&self) -> bool {
        self.properties.is_some()
    }

    /// Tests if weight is defined for the [BeamElement]
    pub(crate) fn has_weight(&self) -> bool {
        self.weight.is_some()
    }

    /// Gets the id of the [BeamElement]
    pub(crate) fn get_id(&self) -> usize {
        self.id
    }

    /// Gets the length of the [BeamElement]
    pub(crate) fn get_length(&self) -> f64 {
        self.length
    }

    /// Gets a reference to the versor of the [BeamElement]
    /// 
    /// The versor is a 3x3 matrix where the columns are the x, y and z axis of the local reference of the element.
    /// 
    /// # Example
    /// '''
    /// use timoshenko_beam::beam::BeamElement;
    /// use timoshenko_beam::beam::BeamProperties;
    /// use timoshenko_beam::beam::BeamShape;
    /// use ndarray::array;
    /// let properties = BeamProperties::new(1.0, 1.0, 1.0, 1.0, 1.0, BeamShape::Rectangular);
    /// 
    /// let beam_element = BeamElement::new(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 0, &[0.0, 1.0, 0.0], Some(properties), None);
    /// 
    /// let versor = beam_element.get_versor();
    /// let x = versor.column(0);
    /// let y = versor.column(1);
    /// let z = versor.column(2);
    /// 
    /// assert_eq!(x, array![1.0, 0.0, 0.0]);
    /// assert_eq!(y, array![0.0, 1.0, 0.0]);
    /// assert_eq!(z, array![0.0, 0.0, 1.0]);
    /// '''
    
    pub(crate) fn get_versor(&self) -> &Array2<f64> {
        &self.versor
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
    /// The convention is that x is positive along the beam axis, y is where the y versor is defined in the [BeamElement] struct
    /// and z is the right-hand rule cross product of the two versors.
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