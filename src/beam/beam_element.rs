/*!
     Contains the elementary element of the beam, the [BeamElement] struct.

    The BeamElement struct is the elementary element of the beam, and it is defined by two points in space,
    the id of the element, the versor of the element, the properties of the element and optionally by the stiffness matrix of the element and the weight of the element.
*/

use crate::{beam::beam_properties::BeamProperties, cross_prodct, modulus, scalar_mult};
use ndarray::{array, s, Array2};

/// This is the elementary element of a [Beam].
///
/// It discretizes the [Beam] into its singular elements.
#[derive(Debug, Clone, PartialEq)]
pub struct BeamElement {
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
    /// ```
    ///  use rusty_beam::beam::BeamElement;
    ///  use rusty_beam::beam::BeamProperties;
    ///  use rusty_beam::beam::BeamShape;
    ///
    ///  let properties = BeamProperties::new(
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         BeamShape::Rectangular
    ///                   );
    /// let beam_element = BeamElement::new(
    ///                         &[0.0, 0.0, 0.0],
    ///                         &[1.0, 0.0, 0.0],
    ///                         0,
    ///                         &[0.0, 1.0, 0.0],
    ///                         Some(properties),
    ///                          None
    ///                    );
    /// ```
    pub fn new(
        u1: &[f64; 3],
        u2: &[f64; 3],
        id: usize,
        y: &[f64; 3],
        properties: Option<BeamProperties>,
        weight: Option<f64>,
    ) -> Self {
        // compute the axis of the beam
        let mut x = [(u2[0] - u1[0]), (u2[1] - u1[1]), (u2[2] - u1[2])];
        // calculate the length
        let length = modulus(&x);
        // normalize the x axis
        x = scalar_mult(&x, 1.0 / length);
        // calculate the modulus of y
        let length_y = modulus(&y);
        // copy y in order to normalize without declaring the input y to be mutable
        let mut loc_y = *y;
        // we only normalyze if the modulus is different from 1
        if length_y != 1.0 {
            loc_y = scalar_mult(&loc_y, 1.0 / length_y);
        }
        // calculate the z axis as the cross product of x and y
        //  since both are normalized oleady, this is normalized as well
        let z = cross_prodct(&x, &loc_y);
        // create the versor matrix
        let versor = array![x, loc_y, z];

        // allocate the stiffness matrix
        let mut stiffness_matrix= Array2::zeros([12, 12]);
        if let Some(prop) = &properties {
           stiffness_matrix = Self::stiffness_from_properties(prop, length, &versor);
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

    fn stiffness_from_properties(properties: &BeamProperties, length: f64, versor: &Array2<f64>) -> Array2<f64> {
            // the stifness matrix in the local reference is computed if we have properties input
            let stiffness_matrix = properties.compute_stiffness_matrix(length);
            // TODO: add change of base to the global reference
            // The rotation matrix to apply is: 
            // [R 0 0 0]
            // [0 I 0 0]
            // [0 0 R 0]
            // [0 0 0 I]
            // to the left of the stiffness matrix.
            // R is the rotation matrix from the local reference to the global reference
            let identity = Array2::eye(3);
            let mut rotation = Array2::zeros([12,12]);
            rotation.slice_mut(s![0..3, 0..3]).assign(&versor);
            rotation.slice_mut(s![3..6, 3..6]).assign(&identity);
            rotation.slice_mut(s![6..9, 6..9]).assign(&versor);
            rotation.slice_mut(s![9..12, 9..12]).assign(&identity);
            // multiply the stiffness matrix by the rotation matrix
            rotation.dot(&stiffness_matrix)
    }

    /// Tests if properties are defined for the [BeamElement]
    ///
    /// Returns true if the properties are defined, false otherwise.
    ///
    /// # Example
    /// ```
    /// use rusty_beam::beam::BeamElement;
    /// use rusty_beam::beam::BeamProperties;
    /// use rusty_beam::beam::BeamShape;
    ///
    /// let properties = BeamProperties::new(
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         BeamShape::Rectangular
    ///                   );
    /// let beam_element = BeamElement::new(
    ///                         &[0.0, 0.0, 0.0],
    ///                         &[1.0, 0.0, 0.0],
    ///                         0,
    ///                         &[0.0, 1.0, 0.0],
    ///                         Some(properties),
    ///                          None
    ///                    );
    /// assert_eq!(beam_element.has_properties(), true);
    ///
    /// let beam_element = BeamElement::new(
    ///                         &[0.0, 0.0, 0.0],
    ///                         &[1.0, 0.0, 0.0],
    ///                         0,
    ///                         &[0.0, 1.0, 0.0],
    ///                         None,
    ///                         None
    ///                    );
    /// assert_eq!(beam_element.has_properties(), false);
    /// ```
    pub fn has_properties(&self) -> bool {
        self.properties.is_some()
    }

    /// Tests if weight is defined for the [BeamElement]
    pub fn has_weight(&self) -> bool {
        self.weight.is_some()
    }

    /// Gets the id of the [BeamElement]
    pub fn get_id(&self) -> usize {
        self.id
    }

    /// Gets the length of the [BeamElement]
    pub fn get_length(&self) -> f64 {
        self.length
    }

    /// Gets a reference to the versor of the [BeamElement]
    ///
    /// The versor is a 3x3 matrix where the columns are the x, y and z axis of the local reference of the element.
    ///
    /// # Example
    /// ```
    /// use rusty_beam::beam::BeamElement;
    /// use rusty_beam::beam::BeamProperties;
    /// use rusty_beam::beam::BeamShape;
    /// use ndarray::array;
    ///
    /// let properties = BeamProperties::new(
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         1.0,
    ///                         BeamShape::Rectangular
    ///                   );
    ///
    /// let beam_element = BeamElement::new(
    ///                         &[0.0, 0.0, 0.0],
    ///                         &[1.0, 0.0, 0.0],
    ///                         0,
    ///                         &[0.0, 1.0, 0.0],
    ///                         Some(properties),
    ///                          None
    ///                    );
    ///
    /// let versor = beam_element.get_versor();
    /// let x = versor.column(0);
    /// let y = versor.column(1);
    /// let z = versor.column(2);
    ///
    /// assert_eq!(x, array![1.0, 0.0, 0.0]);
    /// assert_eq!(y, array![0.0, 1.0, 0.0]);
    /// assert_eq!(z, array![0.0, 0.0, 1.0]);
    /// ```
    pub fn get_versor(&self) -> &Array2<f64> {
        &self.versor
    }
}
