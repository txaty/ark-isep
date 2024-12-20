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

    // Sample random beta, gamma.
    let beta = transcript.squeeze_challenge(Label::ChallengeBeta)?;
    let gamma = transcript.squeeze_challenge(Label::ChallengeGamma)?;

    transcript.append_elements(
        &[
            (Label::G1L, proof.g1_affine_l),
            (Label::G1R, proof.g1_affine_r),
            (Label::G1Ql, proof.g1_affine_ql),
            (Label::G1Qr, proof.g1_affine_qr),
        ]
    )?;

    let delta = transcript.squeeze_challenge(Label::ChallengeDelta)?;
    let epsilon = transcript.squeeze_challenge(Label::ChallengeEpsilon)?;

    transcript.append_element(Label::G1BatchProofAtRandPoint, &proof.batch_proof_at_rand_point)?;
    transcript.append_elements(
        &[
            (Label::FrLAtDelta, proof.l_at_delta),
            (Label::FrRAtDelta, proof.r_at_delta),
            (Label::FrLvAtDelta, proof.lv_at_delta),
            (Label::FrRvAtDelta, proof.rv_at_delta),
            (Label::FrPlAtDelta, proof.pl_at_delta),
            (Label::FrPrAtDelta, proof.pr_at_delta),
            (Label::FrPmAtDelta, proof.pm_at_delta),
            (Label::FrLAtZero, proof.l_at_zero),
            (Label::FrRAtZero, proof.r_at_zero),
        ]
    )?;

    let zeta = transcript.squeeze_challenge(Label::ChallengeZeta)?;

    // Pairing check of batch proof at random point.
    let fr_one = P::ScalarField::one();
    let fr_zl_at_delta = delta.pow(&[pp.size_left_values as u64]) - fr_one;
    let fr_inv_zl_at_delta = fr_zl_at_delta.inverse().ok_or(Error::FailedToInverseFieldElement)?;
    let fr_ql_at_delta = beta + proof.lv_at_delta + gamma * proof.pm_at_delta;
    let fr_ql_at_delta = fr_ql_at_delta * proof.l_at_delta;
    let fr_ql_at_delta = fr_ql_at_delta - proof.pl_at_delta;
    let fr_ql_at_delta = fr_ql_at_delta * fr_inv_zl_at_delta;

    let fr_zr_at_delta = delta.pow(&[pp.size_right_values as u64]) - fr_one;
    let fr_inv_zr_at_delta = fr_zr_at_delta.inverse().ok_or(Error::FailedToInverseFieldElement)?;
    let fr_qr_at_delta = beta + proof.rv_at_delta + gamma * delta;
    let fr_qr_at_delta = fr_qr_at_delta * proof.r_at_delta;
    let fr_qr_at_delta = fr_qr_at_delta - proof.pr_at_delta;
    let fr_qr_at_delta = fr_qr_at_delta * fr_inv_zr_at_delta;

    let g1_list = vec![
        proof.g1_affine_l,
        proof.g1_affine_r,
        proof.g1_affine_ql,
        proof.g1_affine_qr,
        statement.g1_affine_left_values,
        statement.g1_affine_right_values,
        pp.g1_affine_positions_left,
        pp.g1_affine_positions_right,
        pp.g1_affine_position_mappings,
    ];

    let fr_list = vec![
        proof.l_at_delta,
        proof.r_at_delta,
        fr_ql_at_delta,
        fr_qr_at_delta,
        proof.lv_at_delta,
        proof.rv_at_delta,
        proof.pl_at_delta,
        proof.pr_at_delta,
        proof.pm_at_delta,
    ];

    let mut g1_batched = P::G1::zero();
    let mut fr_batched = P::ScalarField::zero();
    let mut fr_pow_epsilon = fr_one;
    g1_list
        .iter()
        .zip(fr_list.iter())
        .for_each(|(g1, &fr)| {
            g1_batched += g1.mul(fr_pow_epsilon);
            fr_batched += fr * fr_pow_epsilon;
            fr_pow_epsilon = fr_pow_epsilon * epsilon;
        });

    let pairing_left = P::pairing(
        g1_batched - pp.g1_affine_srs[0].mul(fr_batched) + proof.batch_proof_at_rand_point.mul(delta),
        pp.g2_affine_srs[0],
    );
    let pairing_right = P::pairing(
        proof.batch_proof_at_rand_point,
        pp.g2_affine_srs[1],
    );
    if pairing_left != pairing_right {
        return Err(Error::Pairing1Failed);
    }

    // Pairing check of batch proof at zero.
    let tmp = proof.g1_affine_r.mul(zeta);
    let tmp = tmp + proof.g1_affine_l;
    let tmp = tmp - pp.g1_affine_srs[0].mul(proof.l_at_zero + proof.r_at_zero * zeta);
    let pairing_left = P::pairing(tmp, pp.g2_affine_srs[0]);
    let pairing_right = P::pairing(proof.batch_proof_at_zero, pp.g2_affine_srs[1]);
    if pairing_left != pairing_right {
        return Err(Error::Pairing2Failed);
    }

    // Sumcheck Lemma.
    if proof.l_at_zero * P::ScalarField::from(pp.size_left_values as u64) != proof.r_at_zero * P::ScalarField::from(pp.size_right_values as u64) {
        return Err(Error::EqualityCheckFailed);
    }

    Ok(())
}