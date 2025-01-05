/*!
 *  Contains the [Beam] struct, which is made by various [BeamElement]s.


*/

/// This is the main struct of the crate, defined to be solved with the Timoshenko beam theory.
///
#[derive(Debug, Clone, PartialEq)]
pub struct Beam {}

impl Beam {}

pub mod beam_element;
pub mod beam_properties;

pub use beam_element::BeamElement;
pub use beam_properties::{BeamProperties, BeamShape};
