use ark_bn254::{Bn254, Fr};
use ark_isep::prover::prove;
use ark_isep::public_parameters::PublicParameters;
use ark_isep::statement::Statement;
use ark_isep::verifier::verify;
use ark_isep::witness::Witness;
use ark_std::{test_rng, UniformRand};
use std::collections::BTreeMap;

fn generate_inputs(num_tx: usize, pow_seg: usize, pow_shared: usize) -> (
    PublicParameters<Bn254>,
    Witness<Bn254>,
    Statement<Bn254>,
) {
    let rng = &mut test_rng();
    let mappings = (0..num_tx).map(|i| (i << pow_seg, i << pow_shared)).collect::<BTreeMap<_, _>>();
    let num_left_values = (1 << pow_seg) * num_tx;
    let num_right_values = (1 << pow_shared) * num_tx;
    let curr_time = std::time::Instant::now();
    let pp = PublicParameters::<Bn254>::builder()
        .left_element_size(num_left_values)
        .right_element_size(num_right_values)
        .positions_left(&(0..num_tx).map(|i| i << pow_seg).collect::<Vec<_>>())
        .positions_right(&(0..num_tx).map(|i| i << pow_shared).collect::<Vec<_>>())
        .position_mappings(&mappings)
        .build(rng).unwrap();
    println!("setup time: {:?} ms", curr_time.elapsed().as_millis());

    let left_witness_values = (0..num_left_values).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let mut right_witness_values = (0..num_right_values).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    mappings.iter().for_each(|(k, v)| {
        right_witness_values[*v] = left_witness_values[*k];
    });

    let witness = Witness::new(&pp, &left_witness_values, &right_witness_values).unwrap();
    let statement = witness.generate_statement(&pp).unwrap();

    (pp, witness, statement)
}
fn main() {
    let (pp, witness, statement) = generate_inputs(1024, 6, 14);
    let curr_time = std::time::Instant::now();
    let proof = prove(&pp, &witness, &statement).unwrap();
    println!("prove time: {:?} ms", curr_time.elapsed().as_millis());
    let curr_time = std::time::Instant::now();
    verify(&pp, &statement, &proof).unwrap();
    println!("verify time: {:?} ms", curr_time.elapsed().as_millis());
}
