use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{FftField, One};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_std::rand::Rng;
use ark_std::{UniformRand, Zero};
use rayon::prelude::*;
use std::marker::PhantomData;
use std::ops::Mul;

/// Minimal KZG functionalities needed for the lookup argument.
///
/// This module provides functionalities to commit to polynomials,
/// open evaluations, and verify proofs using KZG commitments.
///
/// Adapted from:
/// - [geometryxyz/cq](https://github.com/geometryxyz/cq/blob/main/src/kzg.rs)
/// - [caulk-crypto/caulk](https://github.com/caulk-crypto/caulk/blob/main/src/kzg.rs)
pub struct Kzg<C: CurveGroup> {
    _marker: PhantomData<C>,
}

impl<C: CurveGroup> Kzg<C> {
    pub fn commit(affine_srs: &[C::Affine], poly: &DensePolynomial<C::ScalarField>) -> C {
        if affine_srs.len() - 1 < poly.degree() {
            panic!(
                "SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                affine_srs.len()
            );
        }

        VariableBaseMSM::msm_unchecked(affine_srs, &poly.coeffs)
    }

    pub fn open(
        affine_srs: &[C::Affine],
        poly: &DensePolynomial<C::ScalarField>,
        challenge: C::ScalarField,
    ) -> (C::ScalarField, C::Affine) {
        let q =
            poly / &DensePolynomial::from_coefficients_slice(&[-challenge, C::ScalarField::one()]);
        if affine_srs.len() - 1 < q.degree() {
            panic!(
                "Open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                q.degree(),
                affine_srs.len()
            );
        }
        let proof = Self::commit(affine_srs, &q);

        (poly.evaluate(&challenge), proof.into())
    }

    pub fn batch_open(
        affine_srs: &[C::Affine],
        poly_list: &[&DensePolynomial<C::ScalarField>],
        fr_opening: C::ScalarField,
        fr_separation: C::ScalarField,
    ) -> C::Affine {
        let num_polys = poly_list.len();
        let powers_of_sep = powers_of_scalars::<C::ScalarField>(fr_separation, num_polys);

        let mut batched = poly_list[0].clone();
        let rest_batched: DensePolynomial<C::ScalarField> = poly_list[1..]
            .par_iter()
            .zip(powers_of_sep.par_iter().skip(1))
            .map(|(&p_i, &fr_sep_pow_i)| {
                p_i * fr_sep_pow_i // Multiply each polynomial by its
                // corresponding power of `fr_separation`
            })
            .reduce(
                || DensePolynomial::from_coefficients_slice(&[C::ScalarField::zero()]),
                |a, b| a + b,
            );
        batched += &rest_batched;

        let q = &batched
            / &DensePolynomial::from_coefficients_slice(&[-fr_opening, C::ScalarField::one()]);

        if affine_srs.len() - 1 < q.degree() {
            panic!(
                "Batch open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                q.degree(),
                affine_srs.len()
            );
        }

        Self::commit(affine_srs, &q).into()
    }
}

/// Create srs from rng
pub fn unsafe_setup_from_rng<P: Pairing, R: Rng + ?Sized>(
    max_power_g1: usize,
    rng: &mut R,
) -> (
    Vec<P::G1Affine>,
    Vec<P::G2Affine>,
) {
    let tau = P::ScalarField::rand(rng);

    unsafe_setup_from_tau::<P, R>(max_power_g1, tau)
}

/// Create srs from specific tau
pub fn unsafe_setup_from_tau<P: Pairing, R: Rng + ?Sized>(
    max_power_g1: usize,
    tau: P::ScalarField,
) -> (
    Vec<P::G1Affine>,
    Vec<P::G2Affine>,
) {
    let max_power_g2 = max_power_g1 + 1;
    let powers_of_tau_size = max_power_g2 + 1;
    let powers_of_tau = powers_of_scalars::<P::ScalarField>(tau, powers_of_tau_size);
    let g1_srs = srs::<P::G1>(&powers_of_tau, max_power_g1);
    let g2_srs = srs::<P::G2>(&powers_of_tau, max_power_g2);

    (g1_srs, g2_srs)
}

const CHUNK_SIZE: usize = 1024;

fn powers_of_scalars<F: FftField>(s: F, size: usize) -> Vec<F> {
    let num_chunks = (size + CHUNK_SIZE - 1) / CHUNK_SIZE;

    let mut result: Vec<F> = (0..num_chunks)
        .into_par_iter()
        .flat_map(|chunk_index| {
            let start_power = chunk_index * CHUNK_SIZE;
            let mut chunk = Vec::with_capacity(CHUNK_SIZE.min(size - start_power));
            let mut power = s.pow(&[start_power as u64]);

            for _ in 0..CHUNK_SIZE.min(size - start_power) {
                chunk.push(power);
                power *= s;
            }
            chunk
        })
        .collect();

    result.truncate(size);

    result
}

fn srs<C: CurveGroup>(powers_of_tau: &[C::ScalarField], max_power: usize) -> Vec<C::Affine> {
    let generator = C::Affine::generator();

    powers_of_tau
        .par_iter()
        .take(max_power + 1)
        .map(|tp| generator.mul(tp).into())
        .collect()
}
