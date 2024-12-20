use crate::error::Error;
use crate::kzg::Kzg;
use crate::public_parameters::PublicParameters;
use crate::witness::Witness;
use crate::COMPRESS_MOD;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_serialize::CanonicalSerialize;
use blake2::{Blake2b512, Digest};

pub struct Statement<P: Pairing> {
    pub(crate) g1_affine_left_values: P::G1Affine,
    pub(crate) g1_affine_right_values: P::G1Affine,
    pub(crate) hash_representation: Vec<u8>,
}

impl<P: Pairing> Witness<P> {
    pub fn generate_statement(&self, pp: &PublicParameters<P>) -> Result<Statement<P>, Error> {
        let g1_affine_srs = &pp.g1_affine_srs;
        let g1_affine_left_values = Kzg::<P::G1>::commit(g1_affine_srs, &self.poly_left_values).into_affine();
        let g1_affine_right_values = Kzg::<P::G1>::commit(g1_affine_srs, &self.poly_right_values).into_affine();

        let mut buf = Vec::new();
        g1_affine_left_values.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        g1_affine_right_values.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        let mut hasher = Blake2b512::new();
        hasher.update(&buf);
        let hash_representation = hasher.finalize().to_vec();

        Ok(Statement {
            g1_affine_left_values,
            g1_affine_right_values,
            hash_representation,
        })
    }
}