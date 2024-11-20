use crate::domain::{divide_by_vanishing_poly_on_coset_in_place, roots_of_unity};
use crate::error::Error;
use crate::kzg::Kzg;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};
use crate::witness::Witness;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::Zero;
use rayon::prelude::*;

pub struct Proof<P: Pairing> {
    pub(crate) g1_affine_l: P::G1Affine,
    pub(crate) g1_affine_r: P::G1Affine,
    pub(crate) g1_affine_ql: P::G1Affine,
    pub(crate) g1_affine_qr: P::G1Affine,
    pub(crate) batch_proof_at_rand_point: P::G1Affine,
    pub(crate) batch_proof_at_zero: P::G1Affine,
    pub(crate) l_at_delta: P::ScalarField,
    pub(crate) r_at_delta: P::ScalarField,
    pub(crate) lv_at_delta: P::ScalarField,
    pub(crate) rv_at_delta: P::ScalarField,
    pub(crate) pl_at_delta: P::ScalarField,
    pub(crate) pr_at_delta: P::ScalarField,
    pub(crate) pm_at_delta: P::ScalarField,
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
            let eval = beta + witness.left_values[i] + gamma * pp.position_mappings[&i];
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
    let g1_affine_l = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_l).into_affine();

    // Construct the quotient polynomial of the left half.
    let coset_eval_list_l = pp.domain_coset_l.fft(&poly_l);
    let coset_eval_list_left_values = pp.domain_coset_l.fft(&witness.poly_left_values);
    let mut coset_eval_list_ql: Vec<P::ScalarField> = coset_eval_list_l
        .par_iter()
        .zip(coset_eval_list_left_values.par_iter())
        .zip(pp.coset_eval_list_positions_left.par_iter())
        .zip(pp.coset_eval_list_position_mappings.par_iter())
        .map(|(((&l, &v), &p), &m)| l * (beta + v + gamma * m) - p)
        .collect();
    pp.domain_coset_l.ifft_in_place(&mut coset_eval_list_ql);
    let mut poly_coset_coeff_list_ql = coset_eval_list_ql;
    divide_by_vanishing_poly_on_coset_in_place::<P::G1>(&pp.domain_l, &mut
        poly_coset_coeff_list_ql)?;
    let coeff_ql = poly_coset_coeff_list_ql;
    let poly_ql = DensePolynomial::from_coefficients_vec(coeff_ql);
    let g1_affine_ql = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_ql).into_affine();

    // Construct the polynomial representing the right half.
    let mut poly_eval_r = vec![P::ScalarField::zero(); pp.size_right_values];
    let roots_of_unity_r = roots_of_unity::<P>(&pp.domain_r);
    let non_zero_eval_list: Result<Vec<(usize, P::ScalarField)>, Error> = pp.positions_right.par_iter().map(|&i| {
        let eval = beta + witness.right_values[i] + gamma * roots_of_unity_r[i];
        let inv = eval.inverse().ok_or(Error::FailedToInverseFieldElement)?;

        Ok((i, inv))
    }).collect();
    let non_zero_eval_list = non_zero_eval_list?;
    non_zero_eval_list.iter().for_each(|(i, eval)| {
        poly_eval_r[*i] = *eval;
    });
    let poly_coeff_r = pp.domain_r.ifft(&poly_eval_r);
    let poly_r = DensePolynomial::from_coefficients_vec(poly_coeff_r);
    let g1_affine_r = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_r).into_affine();

    // Construct the quotient polynomial of the right half.
    let coset_eval_list_r = pp.domain_coset_r.fft(&poly_r);
    let coset_eval_list_right_values = pp.domain_coset_r.fft(&witness.poly_right_values);
    let mut coset_eval_list_qr: Vec<P::ScalarField> = coset_eval_list_r
        .par_iter()
        .zip(coset_eval_list_right_values.par_iter())
        .zip(pp.coset_eval_list_positions_right.par_iter())
        .zip(pp.roots_of_unity_coset_r.par_iter())
        .map(|(((&r, &e), &p), &c)| r * (beta + e + gamma * c) - p)
        .collect();
    pp.domain_coset_r.ifft_in_place(&mut coset_eval_list_qr);
    let mut poly_coset_coeff_list_qr = coset_eval_list_qr;
    divide_by_vanishing_poly_on_coset_in_place::<P::G1>(&pp.domain_r, &mut poly_coset_coeff_list_qr)?;
    let coeff_qr = poly_coset_coeff_list_qr;
    let poly_qr = DensePolynomial::from_coefficients_vec(coeff_qr);
    let g1_affine_qr = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_qr).into_affine();

    transcript.append_elements(
        &[
            (Label::G1L, g1_affine_l),
            (Label::G1R, g1_affine_r),
            (Label::G1Ql, g1_affine_ql),
            (Label::G1Qr, g1_affine_qr),
        ]
    )?;

    // Sample random delta, phi.
    let delta = transcript.squeeze_challenge(Label::ChallengeDelta)?;
    let epsilon = transcript.squeeze_challenge(Label::ChallengeEpsilon)?;

    let batch_proof_at_rand_point = Kzg::<P::G1>::batch_open(
        &pp.g1_affine_srs,
        &[
            &poly_l,
            &poly_r,
            &poly_ql,
            &poly_qr,
            &witness.poly_left_values,
            &witness.poly_right_values,
            &pp.poly_positions_left,
            &pp.poly_positions_right,
            &pp.poly_position_mappings,
        ],
        delta,
        epsilon,
    );

    transcript.append_element(Label::G1BatchProofAtRandPoint, &batch_proof_at_rand_point)?;

    let l_at_delta = poly_l.evaluate(&delta);
    let r_at_delta = poly_r.evaluate(&delta);
    let lv_at_delta = witness.poly_left_values.evaluate(&delta);
    let rv_at_delta = witness.poly_right_values.evaluate(&delta);
    let pl_at_delta = pp.poly_positions_left.evaluate(&delta);
    let pr_at_delta = pp.poly_positions_right.evaluate(&delta);
    let pm_at_delta = pp.poly_position_mappings.evaluate(&delta);
    
    let fr_zero = P::ScalarField::zero();
    let l_at_zero = poly_l.evaluate(&fr_zero);
    let r_at_zero = poly_r.evaluate(&fr_zero);

    transcript.append_elements(
        &[
            (Label::FrLAtDelta, l_at_delta),
            (Label::FrRAtDelta, r_at_delta),
            (Label::FrLvAtDelta, lv_at_delta),
            (Label::FrRvAtDelta, rv_at_delta),
            (Label::FrPlAtDelta, pl_at_delta),
            (Label::FrPrAtDelta, pr_at_delta),
            (Label::FrPmAtDelta, pm_at_delta),
            (Label::FrLAtZero, l_at_zero),
            (Label::FrRAtZero, r_at_zero),
        ]
    )?;

    let zeta = transcript.squeeze_challenge(Label::ChallengeZeta)?;

    let batch_proof_at_zero = Kzg::<P::G1>::batch_open(
        &pp.g1_affine_srs,
        &[&poly_l, &poly_r],
        fr_zero,
        zeta,
    );


    Ok(Proof {
        g1_affine_l,
        g1_affine_r,
        g1_affine_ql,
        g1_affine_qr,
        batch_proof_at_rand_point,
        batch_proof_at_zero,
        l_at_zero,
        r_at_zero,
        l_at_delta,
        r_at_delta,
        lv_at_delta,
        rv_at_delta,
        pl_at_delta,
        pr_at_delta,
        pm_at_delta,
    })
}