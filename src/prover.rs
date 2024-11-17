use crate::domain::divide_by_vanishing_poly_checked;
use crate::error::Error;
use crate::kzg::Kzg;
use crate::public_parameters::PublicParameters;
use crate::statement::Statement;
use crate::transcript::{Label, Transcript};
use crate::witness::Witness;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial};

pub struct Proof<P: Pairing> {
    pub(crate) g1_affine_sv: P::G1Affine,
    pub(crate) g1_affine_qlv: P::G1Affine,
    pub(crate) g1_affine_qrv: P::G1Affine,
    pub(crate) batch_proof: P::G1Affine,
    pub(crate) sv_at_beta: P::ScalarField,
    pub(crate) qlv_at_beta: P::ScalarField,
    pub(crate) qrv_at_beta: P::ScalarField,
    pub(crate) lv_at_beta: P::ScalarField,
    pub(crate) rv_at_beta: P::ScalarField,
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

    let size_positions = pp.size_positions;
    let left_value_gap = pp.size_left_values / size_positions;
    let subset_left_values = (0..size_positions).map(|i| witness.left_values[i * left_value_gap]).collect::<Vec<_>>();

    let subset_values = subset_left_values;
    let coeff_subset_values = pp.sub_domain.ifft(&subset_values);
    let poly_sv = DensePolynomial::from_coefficients_vec(coeff_subset_values);
    let g1_affine_sv = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_sv).into_affine();


    let poly_qlv = witness.poly_left_values.clone() - &poly_sv;
    let poly_qlv = divide_by_vanishing_poly_checked::<P>(&pp.sub_domain, &poly_qlv)?;
    let g1_affine_qlv = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_qlv).into_affine();

    let poly_qrv = witness.poly_right_values.clone() - &poly_sv;
    let poly_qrv = divide_by_vanishing_poly_checked::<P>(&pp.sub_domain, &poly_qrv)?;
    let g1_affine_qrv = Kzg::<P::G1>::commit(&pp.g1_affine_srs, &poly_qrv).into_affine();

    transcript.append_elements(&[
        (Label::G1Sv, g1_affine_sv),
        (Label::G1Qlv, g1_affine_qlv),
        (Label::G1Qrv, g1_affine_qrv),
    ])?;

    // Sample random beta, gamma.
    let beta = transcript.squeeze_challenge(Label::ChallengeBeta)?;
    let gamma = transcript.squeeze_challenge(Label::ChallengeGamma)?;

    let sv_at_beta = poly_sv.evaluate(&beta);
    let qlv_at_beta = poly_qlv.evaluate(&beta);
    let qrv_at_beta = poly_qrv.evaluate(&beta);
    let lv_at_beta = witness.poly_left_values.evaluate(&beta);
    let rv_at_beta = witness.poly_right_values.evaluate(&beta);

    let batch_proof = Kzg::<P::G1>::batch_open(
        &pp.g1_affine_srs,
        &[
            &witness.poly_left_values,
            &witness.poly_right_values,
            &poly_sv,
            &poly_qlv,
            &poly_qrv,
        ],
        beta,
        gamma,
    );


    Ok(Proof {
        g1_affine_sv,
        g1_affine_qlv,
        g1_affine_qrv,
        batch_proof,
        sv_at_beta,
        qlv_at_beta,
        qrv_at_beta,
        lv_at_beta,
        rv_at_beta,
    })
}