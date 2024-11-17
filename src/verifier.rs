use crate::error::Error;
use crate::prover::Proof;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_std::{One, Zero};
use std::ops::Mul;

pub fn verify<P: Pairing>(
    pp: &PublicParameters<P>,
    statement: &Statement<P>,
    proof: &Proof<P>,
) -> Result<(), Error> {
    let mut transcript = Transcript::<P::ScalarField>::new();
    transcript.append_elements(&[
        (Label::PublicParameters, pp.hash_representation.clone()),
        (Label::Statement, statement.hash_representation.clone()),
    ])?;
    transcript.append_elements(&[
        (Label::G1Sv, proof.g1_affine_sv),
        (Label::G1Qlv, proof.g1_affine_qlv),
        (Label::G1Qrv, proof.g1_affine_qrv),
    ])?;

    // Sample random beta, gamma.
    let beta = transcript.squeeze_challenge(Label::ChallengeBeta)?;
    let gamma = transcript.squeeze_challenge(Label::ChallengeGamma)?;
    

    let fr_one = P::ScalarField::one();
    let z_at_beta = beta.pow(&[pp.size_positions as u64]) - fr_one;
    let inv_z_at_beta = z_at_beta.inverse().ok_or(Error::FailedToInverseFieldElement)?;
    let check = proof.lv_at_beta - proof.sv_at_beta;
    let check = check * inv_z_at_beta;
    if check != proof.qlv_at_beta {
        return Err(Error::Check1Failed);
    }
    
    let check = proof.rv_at_beta - proof.sv_at_beta;
    let check = check * inv_z_at_beta;
    if check != proof.qrv_at_beta {
        return Err(Error::Check2Failed);
    }

    let g1_list = vec![
        statement.g1_affine_left_values,
        statement.g1_affine_right_values,
        proof.g1_affine_sv,
        proof.g1_affine_qlv,
        proof.g1_affine_qrv,
    ];

    let fr_list = vec![
        proof.lv_at_beta,
        proof.rv_at_beta,
        proof.sv_at_beta,
        proof.qlv_at_beta,
        proof.qrv_at_beta,
    ];

    let mut g1_batched = P::G1::zero();
    let mut fr_batched = P::ScalarField::zero();
    let mut fr_pow_gamma = fr_one;
    g1_list
        .iter()
        .zip(fr_list.iter())
        .for_each(|(g1, &fr)| {
            g1_batched += g1.mul(fr_pow_gamma);
            fr_batched += fr * fr_pow_gamma;
            fr_pow_gamma *= gamma;
        });

    let pairing_left = P::pairing(
        g1_batched - pp.g1_affine_srs[0].mul(fr_batched) + proof.batch_proof.mul(beta),
        pp.g2_affine_srs[0],
    );
    let pairing_right = P::pairing(
        proof.batch_proof,
        pp.g2_affine_srs[1],
    );
    if pairing_left != pairing_right {
        return Err(Error::PairingFailed);
    }

    Ok(())
}