use crate::error::Error;
use crate::prover::Proof;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_std::rand::Rng;
use std::ops::Mul;

pub fn verify<P: Pairing, R: Rng + ?Sized>(
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
    let beta = transcript.get_and_append_challenge(Label::ChallengeBeta)?;
    let gamma = transcript.get_and_append_challenge(Label::ChallengeGamma)?;

    // Pairing check of left-side.
    let tmp = pp.g1_affine_position_mappings.mul(gamma);
    let pairing_left = P::pairing(statement.left_commitment + tmp, proof.g2_affine_l);
    let g1_one = pp.g1_affine_srs[0];
    let tmp = proof.g2_affine_l.mul(-beta);
    let tmp = tmp + pp.g2_affine_positions_left;
    let pairing_right = P::multi_pairing(
        &[proof.g1_affine_ql, g1_one],
        &[pp.g2_affine_zl, tmp.into_affine()],
    );
    if pairing_left != pairing_right {
        return Err(Error::Pairing1Failed);
    }

    // Pairing check of right-side.
    let tmp = pp.g1_affine_srs[1].mul(gamma);
    let pairing_left = P::pairing(statement.right_commitment + tmp, proof.g2_affine_r);
    let tmp = proof.g2_affine_r.mul(-beta);
    let tmp = tmp + pp.g2_affine_positions_right;
    let pairing_right = P::multi_pairing(
        &[proof.g1_affine_qr, g1_one],
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
    let delta = transcript.get_and_append_challenge(Label::ChallengeDelta)?;

    // Check KZG batch proof.
    let batch_proof = proof.batch_proof;
    let pairing_left = P::pairing(pp.g1_affine_srs[1], batch_proof);
    let tmp = proof.g2_affine_r.mul(delta);
    let tmp = tmp + proof.g2_affine_l;
    let tmp = tmp - P::G2::generator().mul(-(proof.l_at_zero + proof.r_at_zero * delta));
    let pairing_right = P::pairing(pp.g1_affine_srs[0], tmp);
    if pairing_left != pairing_right {
        return Err(Error::Pairing3Failed);
    }
    

    // Sumcheck Lemma.
    if proof.l_at_zero * P::ScalarField::from(pp.left_element_size as u64) != proof.r_at_zero * P::ScalarField::from(pp.right_element_size as u64) {
        return Err(Error::EqualityCheckFailed);
    }

    Ok(())
}