use crate::error::Error;
use crate::prover::Proof;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
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

    // Pairing check of left-side.
    let g1_tmp = pp.g1_affine_position_mappings.mul(gamma) + statement.g1_affine_left_elements;
    let pairing_left = P::pairing(g1_tmp, proof.g2_affine_l);
    let g1_one = pp.g1_affine_srs[0].into_group();
    let g2_l = proof.g2_affine_l.into_group();
    let g2_positions_left = pp.g2_affine_positions_left.into_group();
    let g2_tmp = g2_l.mul(-beta) + g2_positions_left;
    let g2_zl = pp.g2_affine_zl.into_group();
    let pairing_right = P::multi_pairing(
        &[proof.g1_affine_ql.into_group(), g1_one],
        &[g2_zl, g2_tmp],
    );
    if pairing_left != pairing_right {
        return Err(Error::Pairing1Failed);
    }

    // Pairing check of right-side.
    let g1_tau = pp.g1_affine_srs[1].into_group();
    let pairing_left = P::pairing(
        statement.g1_affine_right_elements + g1_tau.mul(gamma),
        proof.g2_affine_r,
    );
    let tmp = proof.g2_affine_r.mul(-beta);
    let tmp = tmp + pp.g2_affine_positions_right;
    let pairing_right = P::multi_pairing(
        &[proof.g1_affine_qr.into_group(), g1_one],
        &[pp.g2_affine_zr, tmp.into_affine()],
    );
    if pairing_left != pairing_right {
        return Err(Error::Pairing2Failed);
    }

    transcript.append_elements(
        &[
            (Label::G2L, proof.g2_affine_l),
            (Label::G2R, proof.g2_affine_r),
        ]
    )?;
    transcript.append_elements(
        &[
            (Label::G1Ql, proof.g1_affine_ql),
            (Label::G1Qr, proof.g1_affine_qr),
        ]
    )?;
    transcript.append_elements(
        &[
            (Label::FrLAtZero, proof.l_at_zero),
            (Label::FrRAtZero, proof.r_at_zero),
        ]
    )?;
    // Sample random delta.
    let delta = transcript.squeeze_challenge(Label::ChallengeDelta)?;

    // Check KZG batch proof.
    let batch_proof = proof.batch_proof;
    let pairing_left = P::pairing(pp.g1_affine_srs[1], batch_proof);
    let tmp = proof.g2_affine_r.mul(delta);
    let tmp = tmp + proof.g2_affine_l;
    let tmp = tmp - pp.g2_affine_srs[0].mul(proof.l_at_zero + proof.r_at_zero * delta);
    let pairing_right = P::pairing(pp.g1_affine_srs[0], tmp);
    if pairing_left != pairing_right {
        return Err(Error::Pairing3Failed);
    }


    // Sumcheck Lemma.
    if proof.l_at_zero * P::ScalarField::from(pp.size_left_values as u64) != proof.r_at_zero * P::ScalarField::from(pp.size_right_values as u64) {
        return Err(Error::EqualityCheckFailed);
    }

    Ok(())
}