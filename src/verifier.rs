use ark_ec::pairing::Pairing;
use ark_std::rand::Rng;
use crate::error::Error;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};

pub fn verify<P: Pairing, R: Rng + ?Sized>(
    pp: &PublicParameters<P>,
    statement: &Statement<P>,
) -> Result<(), Error> {
    let mut transcript = Transcript::<P::ScalarField>::new();
    transcript.append_elements(&[
        (Label::PublicParameters, pp.hash_representation.clone()),
        (Label::Statement, statement.hash_representation.clone()),
    ])?;
    
    Ok(())
}