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

const COMPRESS_MOD: Compress = Compress::No;

#[cfg(test)]
mod tests {
    use crate::prover::prove;
    use crate::public_parameters::PublicParameters;
    use crate::verifier::verify;
    use crate::witness::Witness;
    use ark_bn254::{Bn254, Fr};
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn end_to_end() {
        let rng = &mut test_rng();

        let pp = PublicParameters::<Bn254>::builder()
            .size_left_values(8)
            .size_right_values(16)
            .size_positions(4)
            .build(rng).unwrap();

        // Correct verification.
        let left_witness_values = (0..8).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let mut right_witness_values = (0..16).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        right_witness_values[0] = left_witness_values[0];
        right_witness_values[4] = left_witness_values[2];
        right_witness_values[8] = left_witness_values[4];
        right_witness_values[12] = left_witness_values[6];

        let witness = Witness::new(&pp, &left_witness_values, &right_witness_values).unwrap();
        let statement = witness.generate_statement(&pp).unwrap();

        let proof = prove::<Bn254>(&pp, &witness, &statement).unwrap();
        verify::<Bn254>(&pp, &statement, &proof).unwrap();

        // Wrong common witness value.
        let mut left_witness_values = left_witness_values;
        left_witness_values[4] = Fr::from(42u64);
        right_witness_values[8] = Fr::from(12u64);

        let witness = Witness::new(&pp, &left_witness_values, &right_witness_values).unwrap();
        let statement = witness.generate_statement(&pp).unwrap();

        assert!(prove::<Bn254>(&pp, &witness, &statement).is_err());
    }
}
