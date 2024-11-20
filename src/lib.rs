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
    use ark_ec::pairing::Pairing;
    use ark_ff::FftField;
    use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
    use ark_std::{test_rng, UniformRand};
    use crate::domain::create_sub_domain;

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

        let proof = prove::<Bn254>(&pp, &witness, &statement).unwrap();
        assert!(verify::<Bn254>(&pp, &statement, &proof).is_err());
    }
    
    #[test]
    fn test_coset() {
        let domain_1 = Radix2EvaluationDomain::<Fr>::new(16).unwrap();
        let domain_2 = create_sub_domain::<Bn254>(&domain_1, 4, 4).unwrap();
        println!("domain_1 elements: {:?}", domain_1.elements().collect::<Vec<_>>());
        println!("domain_2 elements: {:?}", domain_2.elements().collect::<Vec<_>>());
        
        let domain_1_coset = domain_1.get_coset(Fr::GENERATOR).unwrap();
        let domain_2_coset = domain_2.get_coset(Fr::GENERATOR).unwrap();
        println!("domain_1 coset: {:?}", domain_1_coset.elements().collect::<Vec<_>>());
        println!("domain_2 coset: {:?}", domain_2_coset.elements().collect::<Vec<_>>());
    }
}
