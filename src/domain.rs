use crate::error::Error;
use crate::kzg::Kzg;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use rayon::iter::IntoParallelRefMutIterator;
use rayon::prelude::*;

pub(crate) fn create_domain_with_generator<F: FftField>(
    group_gen: F,
    size: usize,
) -> Result<Radix2EvaluationDomain<F>, Error> {
    if !size.is_power_of_two() {
        return Err(Error::InvalidEvaluationDomainSize(size));
    }
    let log_size_of_group = size.trailing_zeros();

    // libfqfft uses > https://github.com/scipr-lab/libfqfft/blob/e0183b2cef7d4c5deb21a6eaf3fe3b586d738fe0/libfqfft/evaluation_domain/domains/basic_radix2_domain.tcc#L33
    if log_size_of_group > F::TWO_ADICITY {
        return Err(Error::InvalidEvaluationDomainSize(size));
    }

    let size = size as u64;
    // Check that it is indeed the 2^(log_size_of_group) root of unity.
    debug_assert_eq!(group_gen.pow([size]), F::one());
    let size_as_field_element = F::from(size);
    let size_inv = size_as_field_element
        .inverse()
        .ok_or(Error::FailedToInverseFieldElement)?;
    let group_gen_inv = group_gen
        .inverse()
        .ok_or(Error::FailedToInverseFieldElement)?;

    Ok(Radix2EvaluationDomain {
        size,
        log_size_of_group,
        size_as_field_element,
        size_inv,
        group_gen,
        group_gen_inv,
        offset: F::one(),
        offset_inv: F::one(),
        offset_pow_size: F::one(),
    })
}

pub(crate) fn vanishing_poly_commitment_affine<C: CurveGroup>(
    affine_srs: &[C::Affine],
    domain: &Radix2EvaluationDomain<C::ScalarField>,
) -> C::Affine {
    let vanishing_poly: DensePolynomial<C::ScalarField> = domain.vanishing_polynomial().into();

    Kzg::<C>::commit(&affine_srs, &vanishing_poly).into_affine()
}

pub(crate) fn roots_of_unity<P: Pairing>(
    domain: &Radix2EvaluationDomain<P::ScalarField>,
) -> Vec<P::ScalarField> {
    domain.elements().collect()
}

pub fn divide_by_vanishing_poly_on_coset_in_place<C: CurveGroup>(
    domain: &Radix2EvaluationDomain<C::ScalarField>,
    evaluations: &mut [C::ScalarField],
) -> Result<(), Error> {
    let vanishing_poly_eval = domain.evaluate_vanishing_polynomial(C::ScalarField::GENERATOR);
    let inv_vanishing_poly_eval = vanishing_poly_eval
        .inverse()
        .ok_or(Error::FailedToInverseFieldElement)?;
    evaluations
        .par_iter_mut()
        .for_each(|eval| *eval *= &inv_vanishing_poly_eval);

    Ok(())
}