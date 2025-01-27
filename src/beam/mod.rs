/*!
 *  Contains the [Beam] struct, which is made by various [BeamElement]s.


*/
pub mod beam_element;
pub mod beam_properties;

use std::collections::HashMap;

pub use beam_element::BeamElement;
pub use beam_properties::{BeamProperties, BeamShape};
use ndarray::{s, Array2};

/// This is the main struct of the crate, defined to be solved with the Timoshenko beam theory.
///
#[derive(Debug, Clone, PartialEq)]
pub struct Beam {
    /// id of the beam
    id: usize,
    /// optional name of the beam
    name: Option<String>,
    /// the topology of the beam, with nodes lists and elements
    topology: BeamTopology,
    /// external applied loads on the specific beam node
    /// the key is the node id, the value is the load
    loads:HashMap<usize,Load>,
    /// boundary conditions of the beam
    /// the key is the node id, the value is the boundary conditions
    conditions: HashMap<usize,BeamBoundaryConditions>,
    /// computed stiffness matrix of the beam
    stiffness_mmatrix: Array2<f64>,
}

impl Beam {

    pub fn new(id:usize, name:Option<String>, topology:BeamTopology, loads:HashMap<usize,Load>, conditions:HashMap<usize,BeamBoundaryConditions>) -> Self {
        let stiffness_mmatrix = Array2::zeros((topology.nodes.len()*6,topology.nodes.len()*6));
        Self {
            id,
            name,
            topology,
            loads,
            conditions,
            stiffness_mmatrix,
        }
    }

}

/// This struct contains the topology of the beam.
/// 
#[derive(Debug, Clone, PartialEq)]
pub struct BeamTopology {
    /// list of nodes in the beam
    nodes: Vec<usize>,
    /// list of elements in the beam with the relative id
    elements: HashMap<(usize,usize),BeamElement>,
}

impl BeamTopology {
    /// creates a new empty topology
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            elements: HashMap::new(),
        }
    }

    /// Adds a new node to the topology
    pub fn add_element(&mut self, start: usize, end: usize, element: BeamElement) {
        if !self.nodes.contains(&start) {
            self.nodes.push(start);
        }
        if !self.nodes.contains(&end) {
            self.nodes.push(end);
        }
        self.elements.insert((start,end),element);
    }

    /// Adds a straight segment to the topology
    /// 
    /// The nodes are given in order, and the elements are given in 
    pub fn add_straight_segment(&mut self, nodes: &[usize], elements: &[BeamElement]) {
        if nodes.len() != elements.len() + 1 {
            panic!("The number of nodes must be one more than the number of elements");
        }
        let mut start = nodes[0];
        for (id,node) in nodes.iter().enumerate().skip(1) {
            let element = elements[id].clone();
            self.add_element(start,*node,element);
            start = *node;
        }
    }

    /// This function checks if the topology is valid
    /// 
    /// each node must be connected to at least one element
    // Does this makes sense? It should be automatically valid by the add_element function
    pub fn check_topology(&self) -> bool {
        for id in self.nodes.iter() {
            if self.elements.iter().filter(|(key,_)| key.0 == *id || key.1 == *id).count() == 0 {
                return false;
            }
        }
        true
    }



    /// This function computes the stiffness matrix of the beam starting from its elements beams
    /// 
    fn get_stiffness_matrix(&self) -> Array2<f64> {
        let mut stiffness_matrix = Array2::zeros((self.nodes.len()*6,self.nodes.len()*6));
        for (key,element) in self.elements.iter() {
            let start = key.0;
            let end = key.1;
            let element_stiffness = element.get_stiffness_matrix();
            let start_index = start*6;
            let end_index = end*6;
            stiffness_matrix.slice_mut(s![start_index..start_index+6,start_index..start_index+6]).assign(&element_stiffness.slice(s![0..6,0..6]));
            stiffness_matrix.slice_mut(s![end_index..end_index+6,end_index..end_index+6]).assign(&element_stiffness.slice(s![6..12,6..12]));
            stiffness_matrix.slice_mut(s![start_index..start_index+6,end_index..end_index+6]).assign(&element_stiffness.slice(s![0..6,6..12]));
            stiffness_matrix.slice_mut(s![end_index..end_index+6,start_index..start_index+6]).assign(&element_stiffness.slice(s![6..12,0..6]));
        }
        stiffness_matrix
    }
}

/// This enum contains the different types of loads that can be applied to a beam node.
#[derive(Debug, Clone, PartialEq)]
pub enum Load {
    /// Point load defined on a beam node with a vectoir in the 3D space (global coordinates)
    PointLoad([f64;3]),
    /// Distributed load defined up to a beam node with a vector in the 3D space (global coordinates)
    DistributedLoad([f64;3],(usize,[f64;3])),
    /// A point moment defined on a beam node with a vector in the 3D space (global coordinates)
    Moment([f64;3]),
}

/// This struct contains the boundary conditions of a beam node.
/// 
/// The boundary conditions are defined by the fixed degrees of freedom and the prescribed displacements od a node id.
/// The fixed degrees of freedom are defined by a boolean array of 6 elements, where true means that the degree of freedom is fixed.
/// The prescribed displacement are defined wiht a [DisplacementCondition] struct.
#[derive(Debug, Clone, PartialEq)]  
pub struct BeamBoundaryConditions {
    fixed_dofs: HashMap<usize,[bool;6]>,
    prescribed_displacements : HashMap<usize,DisplacementCondition>,
}

impl BeamBoundaryConditions {
    pub fn new(id: usize) -> Self {
        Self {
            fixed_dofs: HashMap::new(),
            prescribed_displacements: HashMap::new(),
        }
    }
    pub fn add_fixed(&mut self,id: usize, fixed_dofs: [bool;6])  {
        self.fixed_dofs.insert(id,fixed_dofs);
    }

    pub fn add_prescribed(&mut self,id: usize, condition: DisplacementCondition) {
        self.prescribed_displacements.insert(id,condition);
    }
    
}

/// This is a prescribed displacement condition of a beam node.
/// 
/// The condition is determined by the axis of the displacement and the value of the displacement.
/// The axis is in the local coordinate system of the beam element, with x along the beam axis, y where the y versor is defined in the [BeamElement] struct, and z by consequence
/// The other 3 axis are the rotational axis around the x, y and z axis.
#[derive(Debug, Clone, PartialEq)]
pub struct DisplacementCondition {
    /// the three translations and rotations
    axis: [bool;6],
    /// the values of each prescribed displacement.
    values: Vec<f64>,
}

impl DisplacementCondition {
    pub fn new(axis: [bool;6], values: Vec<f64>) -> Self {
        Self {
            axis,
            values,
        }
    }
    
}

