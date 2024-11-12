use ark_serialize::Compress;

pub mod kzg;
pub mod error;
pub mod public_parameters;
mod domain;
pub mod prover;
pub mod verifier;
pub mod witness;
pub mod statement;
mod transcript;

const COMPRESS_MOD: Compress = Compress::Yes;

#[cfg(test)]
mod tests {
    use crate::prover::prover;
    use crate::public_parameters::PublicParameters;
    use crate::verifier::verify;
    use crate::witness::Witness;
    use ark_bn254::{Bn254, Fr};
    use ark_std::{test_rng, UniformRand};
    use std::collections::BTreeMap;

    #[test]
    fn end_to_end() {
        let rng = &mut test_rng();
        let mut mappings = BTreeMap::new();
        mappings.insert(0, 0);
        mappings.insert(2, 4);
        mappings.insert(4, 8);
        mappings.insert(6, 12);

        let pp = PublicParameters::<Bn254>::builder()
            .left_element_size(8)
            .right_element_size(16)
            .positions_left(&[0, 2, 4, 6])
            .positions_right(&[0, 4, 8, 12])
            .position_mappings(&mappings)
            .build(rng).unwrap();
        
        let left_witness_values = (0..8).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let mut right_witness_values = (0..16).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        right_witness_values[0] = left_witness_values[0];
        right_witness_values[4] = left_witness_values[2];
        right_witness_values[8] = left_witness_values[4];
        right_witness_values[12] = left_witness_values[6];

        let witness = Witness::new(&pp, &left_witness_values, &right_witness_values).unwrap();
        let statement = witness.generate_statement(&pp).unwrap();

        let proof = prover::<Bn254>(&pp, &witness, &statement).unwrap();
        verify::<Bn254>(&pp, &statement, &proof).unwrap();
    }
}
