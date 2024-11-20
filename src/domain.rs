use crate::error::Error;
use ark_ec::pairing::Pairing;
use ark_ff::{FftField, Field};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_std::One;

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

pub(crate) fn roots_of_unity<P: Pairing>(
    domain: &Radix2EvaluationDomain<P::ScalarField>,
) -> Vec<P::ScalarField> {
    domain.elements().collect()
}



pub(crate) fn create_sub_domain<P: Pairing>(
    original_domain: &Radix2EvaluationDomain<P::ScalarField>,
    order: usize,
    segment_size: usize,
) -> Result<Radix2EvaluationDomain<P::ScalarField>, Error> {
    if segment_size == 0 {
        return Err(Error::InvalidSegmentSize(segment_size));
    }
    if segment_size == 1 {
        return Ok(original_domain.clone());
    }
    if !order.is_power_of_two() {
        return Err(Error::InvalidEvaluationDomainSize(order));
    }
    let original_order = original_domain.size();
    if segment_size > original_order {
        return Err(Error::InvalidSegmentSize(segment_size));
    }
    if segment_size == original_order {
        let domain_1 = Radix2EvaluationDomain::<P::ScalarField>::new(1)
            .ok_or(Error::FailedToCreateEvaluationDomain)?;
        return Ok(domain_1);
    }

    let size: u64 = order as u64;
    let log_size_of_group = order.trailing_zeros();
    if log_size_of_group > <P::ScalarField as FftField>::TWO_ADICITY {
        return Err(Error::InvalidEvaluationDomainSize(order));
    }

    let roots_of_unity_larger_domain = roots_of_unity::<P>(&original_domain);
    let group_gen = roots_of_unity_larger_domain[segment_size];
    let size_as_field_element = P::ScalarField::from(size);
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
        offset: P::ScalarField::one(),
        offset_inv: P::ScalarField::one(),
        offset_pow_size: P::ScalarField::one(),
    })
}