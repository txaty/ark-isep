use crate::domain::{divide_by_vanishing_poly_on_coset_in_place, roots_of_unity};
use crate::error::Error;
use crate::kzg::Kzg;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};
use crate::witness::Witness;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::Zero;
use rayon::prelude::*;

pub struct Proof<P: Pairing> {
    pub(crate) g2_affine_l: P::G2Affine,
    pub(crate) g2_affine_r: P::G2Affine,
    pub(crate) g1_affine_ql: P::G1Affine,
    pub(crate) g1_affine_qr: P::G1Affine,
    pub(crate) batch_proof: P::G2Affine,
    pub(crate) l_at_zero: P::ScalarField,
    pub(crate) r_at_zero: P::ScalarField,
}

pub fn prove<P: Pairing>(
    pp: &PublicParameters<P>,
    witness: &Witness<P>,
    statement: &Statement<P>,
) -> Result<Proof<P>, Error> {
    let mut transcript = Transcript::<P::ScalarField>::new();
    transcript.append_elements(&[
        (Label::PublicParameters, pp.hash_representation.clone()),
        (Label::Statement, statement.hash_representation.clone()),
    ])?;

    // Sample random beta, gamma.
    let beta = transcript.squeeze_challenge(Label::ChallengeBeta)?;
    let gamma = transcript.squeeze_challenge(Label::ChallengeGamma)?;

    // Construct the polynomial representing the left half.
    let mut poly_eval_l = vec![P::ScalarField::zero(); pp.size_left_values];
    let non_zero_eval_list: Result<Vec<(usize, P::ScalarField)>, Error> = pp.positions_left
        .par_iter() // Parallel iterator
        .map(|&i| {
            let eval = beta + witness.left_elements[i] + gamma * pp.position_mappings[&i];
            let inv = eval.inverse().ok_or(Error::FailedToInverseFieldElement)?;

            Ok((i, inv))
        })
        .collect();
    let non_zero_eval_list = non_zero_eval_list?;
    non_zero_eval_list.iter().for_each(|(i, eval)| {
        poly_eval_l[*i] = *eval;
    });
    let poly_coeff_l = pp.domain_l.ifft(&poly_eval_l);
    let poly_l = DensePolynomial::from_coefficients_vec(poly_coeff_l);
    let g2_affine_l = Kzg::<P::G2>::commit(&pp.g2_affine_srs, &poly_l).into_affine();

    // Construct the quotient polynomial of the left half.
    let domain_coset_l = pp.domain_l
        .get_coset(P::ScalarField::GENERATOR)
        .ok_or(Error::FailedToCreateCosetOfEvaluationDomain)?;

    let poly_coset_eval_list_l = domain_coset_l.fft(&poly_l);
    let poly_coset_eval_list_left_elements = domain_coset_l.fft(&witness.poly_left_elements);
    let poly_coset_eval_list_positions_left = domain_coset_l.fft(&pp.poly_positions_left);
    let poly_coset_eval_position_mappings = domain_coset_l.fft(&pp.poly_position_mappings);
    let mut poly_coset_eval_list_ql: Vec<P::ScalarField> = poly_coset_eval_list_l
        .par_iter()
        .zip(poly_coset_eval_list_left_elements.par_iter())
        .zip(poly_coset_eval_list_positions_left.par_iter())
        .zip(poly_coset_eval_position_mappings.par_iter())
        .map(|(((&l, &e), &p), &m)| l * (beta + e + gamma * m) - p)
        .collect();
    domain_coset_l.ifft_in_place(&mut poly_coset_eval_list_ql);
    let mut poly_coset_coeff_list_ql = poly_coset_eval_list_ql;
    divide_by_vanishing_poly_on_coset_in_place::<P::G1>(&pp.domain_l, &mut
        poly_coset_coeff_list_ql)?;
    let coeff_ql = poly_coset_coeff_list_ql;
    let poly_ql = DensePolynomial::from_coefficients_vec(coeff_ql);
    let g1_affine_ql = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_ql).into_affine();

    // Construct the polynomial representing the right half.
    let mut poly_eval_r = vec![P::ScalarField::zero(); pp.size_right_values];
    let roots_of_unity_r = roots_of_unity::<P>(&pp.domain_r);
    let non_zero_eval_list: Result<Vec<(usize, P::ScalarField)>, Error> = pp.positions_right.par_iter().map(|&i| {
        let eval = beta + witness.right_elements[i] + gamma * roots_of_unity_r[i];
        let inv = eval.inverse().ok_or(Error::FailedToInverseFieldElement)?;

        Ok((i, inv))
    }).collect();
    let non_zero_eval_list = non_zero_eval_list?;
    non_zero_eval_list.iter().for_each(|(i, eval)| {
        poly_eval_r[*i] = *eval;
    });
    let poly_coeff_r = pp.domain_r.ifft(&poly_eval_r);
    let poly_r = DensePolynomial::from_coefficients_vec(poly_coeff_r);
    let g2_affine_r = Kzg::<P::G2>::commit(&pp.g2_affine_srs, &poly_r).into_affine();

    // Construct the quotient polynomial of the right half.
    let domain_coset_r = pp.domain_r
        .get_coset(P::ScalarField::GENERATOR)
        .ok_or(Error::FailedToCreateCosetOfEvaluationDomain)?;
    let poly_coset_eval_list_r = domain_coset_r.fft(&poly_r);
    let poly_coset_eval_list_right_elements = domain_coset_r.fft(&witness.poly_right_elements);
    let poly_coset_eval_list_positions_right = domain_coset_r.fft(&pp.poly_positions_right);
    let roots_of_unity_coset_r = roots_of_unity::<P>(&domain_coset_r);
    let mut poly_coset_eval_list_qr: Vec<P::ScalarField> = poly_coset_eval_list_r
        .par_iter()
        .zip(poly_coset_eval_list_right_elements.par_iter())
        .zip(poly_coset_eval_list_positions_right.par_iter())
        .zip(roots_of_unity_coset_r.par_iter())
        .map(|(((&r, &e), &p), &c)| r * (beta + e + gamma * c) - p)
        .collect();
    domain_coset_r.ifft_in_place(&mut poly_coset_eval_list_qr);
    let mut poly_coset_coeff_list_qr = poly_coset_eval_list_qr;
    divide_by_vanishing_poly_on_coset_in_place::<P::G1>(&pp.domain_r, &mut poly_coset_coeff_list_qr)?;
    let coeff_qr = poly_coset_coeff_list_qr;
    let poly_qr = DensePolynomial::from_coefficients_vec(coeff_qr);
    let g1_affine_qr = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_qr).into_affine();

    transcript.append_elements(
        &[
            (Label::G2L, g2_affine_l),
            (Label::G2R, g2_affine_r),
        ]
    )?;
    transcript.append_elements(
        &[
            (Label::G1Ql, g1_affine_ql),
            (Label::G1Qr, g1_affine_qr),
        ]
    )?;

    let fr_zero = P::ScalarField::zero();
    let l_at_zero = poly_l.evaluate(&fr_zero);
    let r_at_zero = poly_r.evaluate(&fr_zero);

    transcript.append_elements(
        &[
            (Label::FrLAtZero, l_at_zero),
            (Label::FrRAtZero, r_at_zero),
        ]
    )?;
    // Sample random delta.
    let delta = transcript.squeeze_challenge(Label::ChallengeDelta)?;

    let mut poly_batched = poly_l + &poly_r * delta;
    poly_batched.coeffs.drain(0..1);
    let batch_proof = Kzg::<P::G2>::commit(&pp.g2_affine_srs, &poly_batched).into_affine();

    Ok(Proof {
        g2_affine_l,
        g2_affine_r,
        g1_affine_ql,
        g1_affine_qr,
        batch_proof,
        l_at_zero,
        r_at_zero,
    })
}