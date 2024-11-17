use crate::error::Error;
use crate::prover::Proof;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};
use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_poly::Polynomial;
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
            (Label::G2L, proof.g1_affine_l),
            (Label::G2R, proof.g1_affine_r),
        ]
    )?;
    transcript.append_elements(
        &[
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
            (Label::FrLAtZero, proof.l_at_zero),
            (Label::FrRAtZero, proof.r_at_zero),
        ]
    )?;

    let zeta = transcript.squeeze_challenge(Label::ChallengeZeta)?;

    // Pairing check of batch proof at random point.
    let fr_one = P::ScalarField::one();
    let fr_zl_at_delta = delta.pow(&[pp.size_left_values as u64]) - fr_one;
    let fr_inv_zl_at_delta = fr_zl_at_delta.inverse().ok_or(Error::FailedToInverseFieldElement)?;
    let fr_pl_at_delta = pp.poly_positions_left.evaluate(&delta);
    let fr_pm_at_delta = pp.poly_position_mappings.evaluate(&delta);
    let fr_ql_at_delta = beta + proof.lv_at_delta + gamma * fr_pm_at_delta;
    let fr_ql_at_delta = fr_ql_at_delta * proof.l_at_delta;
    let fr_ql_at_delta = fr_ql_at_delta - fr_pl_at_delta;
    let fr_ql_at_delta = fr_ql_at_delta * fr_inv_zl_at_delta;

    let fr_zr_at_delta = delta.pow(&[pp.size_right_values as u64]) - fr_one;
    let fr_inv_zr_at_delta = fr_zr_at_delta.inverse().ok_or(Error::FailedToInverseFieldElement)?;
    let fr_pr_at_delta = pp.poly_positions_right.evaluate(&delta);
    let fr_qr_at_delta = beta + proof.rv_at_delta + gamma * delta;
    let fr_qr_at_delta = fr_qr_at_delta * proof.r_at_delta;
    let fr_qr_at_delta = fr_qr_at_delta - fr_pr_at_delta;
    let fr_qr_at_delta = fr_qr_at_delta * fr_inv_zr_at_delta;

    let g1_list = vec![
        proof.g1_affine_l,
        proof.g1_affine_r,
        proof.g1_affine_ql,
        proof.g1_affine_qr,
        statement.g1_affine_left_values,
        statement.g1_affine_right_values,
        pp.g1_affine_position_mappings,
    ];

    let fr_list = vec![
        proof.l_at_delta,
        proof.r_at_delta,
        fr_ql_at_delta,
        fr_qr_at_delta,
        proof.lv_at_delta,
        proof.rv_at_delta,
        fr_pm_at_delta,
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


    // Pairing check of left-side.
    // let g1_tmp = pp.g1_affine_position_mappings.mul(gamma) + statement.g1_affine_left_values;
    // let pairing_left = P::pairing(g1_tmp, proof.g1_affine_l);
    // let g1_one = pp.g1_affine_srs[0].into_group();
    // let g2_l = proof.g1_affine_l.into_group();
    // let g2_positions_left = pp.g1_affine_positions_left.into_group();
    // let g2_tmp = g2_l.mul(-beta) + g2_positions_left;
    // let g2_zl = pp.g2_affine_zl.into_group();
    // let pairing_right = P::multi_pairing(
    //     &[proof.g1_affine_ql.into_group(), g1_one],
    //     &[g2_zl, g2_tmp],
    // );
    // if pairing_left != pairing_right {
    //     return Err(Error::Pairing1Failed);
    // }

    // Pairing check of right-side.
    // let g1_tau = pp.g1_affine_srs[1].into_group();
    // let pairing_left = P::pairing(
    //     statement.g1_affine_right_values + g1_tau.mul(gamma),
    //     proof.g1_affine_r,
    // );
    // let tmp = proof.g1_affine_r.mul(-beta);
    // let tmp = tmp + pp.g1_affine_positions_right;
    // let pairing_right = P::multi_pairing(
    //     &[proof.g1_affine_qr.into_group(), g1_one],
    //     &[pp.g2_affine_zr, tmp.into_affine()],
    // );
    // if pairing_left != pairing_right {
    //     return Err(Error::Pairing2Failed);
    // }

    // transcript.append_elements(
    //     &[
    //         (Label::G2L, proof.g1_affine_l),
    //         (Label::G2R, proof.g1_affine_r),
    //     ]
    // )?;
    // transcript.append_elements(
    //     &[
    //         (Label::G1Ql, proof.g1_affine_ql),
    //         (Label::G1Qr, proof.g1_affine_qr),
    //     ]
    // )?;
    // transcript.append_elements(
    //     &[
    //         (Label::FrLAtZero, proof.l_at_zero),
    //         (Label::FrRAtZero, proof.r_at_zero),
    //     ]
    // )?;
    // // Sample random delta.
    // let delta = transcript.squeeze_challenge(Label::ChallengeDelta)?;

    // Check KZG batch proof at zero.
    // let batch_proof_at_zero = proof.batch_proof_at_zero;
    // let pairing_left = P::pairing(pp.g1_affine_srs[1], batch_proof_at_zero);
    // let tmp = proof.g1_affine_r.mul(delta);
    // let tmp = tmp + proof.g1_affine_l;
    // let tmp = tmp - pp.g2_affine_srs[0].mul(proof.l_at_zero + proof.r_at_zero * delta);
    // let pairing_right = P::pairing(pp.g1_affine_srs[0], tmp);
    // if pairing_left != pairing_right {
    //     return Err(Error::Pairing3Failed);
    // }

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