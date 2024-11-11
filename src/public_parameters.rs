use crate::domain::{create_domain_with_generator, roots_of_unity, vanishing_poly_commitment_affine};
use crate::error::Error;
use crate::kzg::unsafe_setup_from_rng;
use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use ark_std::rand::Rng;
use ark_std::{UniformRand, Zero};
use rayon::prelude::*;
use std::cmp::max;
use std::collections::BTreeMap;

#[derive(Debug)]
pub struct PublicParameters<P: Pairing> {
    left_element_size: usize,
    right_element_size: usize,
    index_size: usize,

    pub g1_affine_srs: Vec<P::G1Affine>,
    pub g2_affine_srs: Vec<P::G2Affine>,

    pub g2_affine_zl: P::G2Affine,
    pub g2_affine_zr: P::G2Affine,
    pub g2_affine_zi: P::G2Affine,

    pub domain_l: Radix2EvaluationDomain<P::ScalarField>,
    pub domain_r: Radix2EvaluationDomain<P::ScalarField>,
    pub domain_i: Radix2EvaluationDomain<P::ScalarField>,

    pub indices_left: Vec<usize>,
    pub indices_right: Vec<usize>,
    pub poly_indices_left: DensePolynomial<P::ScalarField>,
    pub poly_indices_right: DensePolynomial<P::ScalarField>,

    pub index_mapping: BTreeMap<usize, usize>,
    pub poly_index_mapping: DensePolynomial<P::ScalarField>,
}

impl<P: Pairing> PublicParameters<P> {
    pub fn builder() -> PublicParameters<P> {
        PublicParametersBuilder::<P>::defa
    }
}

pub struct PublicParametersBuilder<P: Pairing> {
    left_element_size: Option<usize>,
    right_element_size: Option<usize>,
    index_size: Option<usize>,
    tau: Option<P::ScalarField>,
    domain_generator_l: Option<P::ScalarField>,
    domain_generator_r: Option<P::ScalarField>,
    domain_generator_i: Option<P::ScalarField>,
    indices_left: Option<Vec<usize>>,
    indices_right: Option<Vec<usize>>,
    index_mapping: Option<BTreeMap<usize, usize>>,
}

impl<P: Pairing> PublicParametersBuilder<P> {
    fn default() -> Self {
        Self {
            left_element_size: None,
            right_element_size: None,
            index_size: None,
            tau: None,
            domain_generator_l: None,
            domain_generator_r: None,
            domain_generator_i: None,
            indices_left: None,
            indices_right: None,
            index_mapping: None,
        }
    }

    pub fn left_element_size(mut self, size: usize) -> Self {
        self.left_element_size = Some(size);
        self
    }

    pub fn right_element_size(mut self, size: usize) -> Self {
        self.right_element_size = Some(size);
        self
    }

    pub fn index_size(mut self, size: usize) -> Self {
        self.index_size = Some(size);
        self
    }

    pub fn tau(mut self, tau: P::ScalarField) -> Self {
        self.tau = Some(tau);
        self
    }

    pub fn domain_generator_l(mut self, gen: P::ScalarField) -> Self {
        self.domain_generator_l = Some(gen);
        self
    }

    pub fn domain_generator_r(mut self, gen: P::ScalarField) -> Self {
        self.domain_generator_r = Some(gen);
        self
    }

    pub fn domain_generator_i(mut self, gen: P::ScalarField) -> Self {
        self.domain_generator_i = Some(gen);
        self
    }

    pub fn indices_left(mut self, indices: &[usize]) -> Self {
        self.indices_left = Some(indices.to_vec());
        self
    }

    pub fn indices_right(mut self, indices: &[usize]) -> Self {
        self.indices_right = Some(indices.to_vec());
        self
    }

    pub fn index_mapping(mut self, mapping: &BTreeMap<usize, usize>) -> Self {
        self.index_mapping = Some(mapping.clone());
        self
    }

    pub fn build<R: Rng + ?Sized>(self, rng: &mut R) -> Result<PublicParameters<P>, Error> {
        let left_element_size = self.left_element_size.ok_or(Error::MissingParameter("Left \
        Element Size"))?;
        validate_input(left_element_size, None)?;
        let right_element_size = self.right_element_size.ok_or(Error::MissingParameter("Right \
        Element Size"))?;
        validate_input(right_element_size, None)?;
        let pow_of_tau_g1 = max(left_element_size, right_element_size);
        let index_size = self.index_size.ok_or(Error::MissingParameter("Index Size"))?;
        validate_input(index_size, Some(pow_of_tau_g1))?;

        let tau = self.tau.unwrap_or(P::ScalarField::rand(rng));
        let (g1_affine_srs, g2_affine_srs) = unsafe_setup_from_rng(pow_of_tau_g1, tau);

        let (domain_l, g2_affine_zl) = create_domain_and_vanishing_poly_commitment(self
                                                                                       .domain_generator_l, left_element_size, &g2_affine_srs)?;
        let (domain_r, g2_affine_zr) = create_domain_and_vanishing_poly_commitment(self
                                                                                       .domain_generator_r, right_element_size, &g2_affine_srs)?;
        let (domain_i, g2_affine_zi) = create_domain_and_vanishing_poly_commitment(self
                                                                                       .domain_generator_i, index_size, &g2_affine_srs)?;

        let indices_left = self.indices_left.ok_or(Error::LeftIndicesCannotBeNone)?;
        let indices_right = self.indices_right.ok_or(Error::RightIndicesCannotBeNone)?;
        let index_mapping = self.index_mapping.ok_or(Error::IndexMappingCannotBeNone)?;

        let coeff_indices_left = domain_i.ifft(&indices_left);
        let poly_indices_left = DensePolynomial::from_coefficients_vec(coeff_indices_left);
        let coeff_indices_right = domain_i.ifft(&indices_right);
        let poly_indices_right = DensePolynomial::from_coefficients_vec(coeff_indices_right);

        let mut poly_eval_index_mapping: Vec<P::ScalarField> = vec![P::ScalarField::zero();
                                                                    left_element_size];
        let roots_of_unity_r = roots_of_unity(&domain_r);
        index_mapping.par_iter().for_each(|key, value| {
            poly_eval_index_mapping[key] = roots_of_unity_r[value];
        });
        domain_l.ifft_in_place(&mut poly_eval_index_mapping);
        let coeff_index_mapping = poly_eval_index_mapping;
        let poly_index_mapping = DensePolynomial::from_coefficients_vec(coeff_index_mapping);

        Ok(PublicParameters {
            left_element_size,
            right_element_size,
            index_size,
            g1_affine_srs,
            g2_affine_srs,
            g2_affine_zl,
            g2_affine_zr,
            g2_affine_zi,
            domain_l,
            domain_r,
            domain_i,
            indices_left,
            indices_right,
            poly_indices_left,
            poly_indices_right,
            index_mapping,
            poly_index_mapping,
        })
    }
}

fn validate_input(input: usize, max_limit: Option<usize>) -> Result<(), Error> {
    if !input.is_power_of_two() {
        return Err(Error::InputShouldBePowerOfTwo(input));
    }

    if max_limit.map_or(false, |max| input > max) {
        return Err(Error::InputIsTooLarge(input));
    }

    Ok(())
}

fn create_domain_and_vanishing_poly_commitment<P: Pairing>(
    domain_generator: Option<P::ScalarField>,
    domain_size: usize,
    g2_affine_srs: &[P::G2Affine],
) -> Result<(Radix2EvaluationDomain<P::ScalarField>, P::G2Affine), Error> {
    let domain = domain_generator
        .map_or_else(
            || Radix2EvaluationDomain::<P::ScalarField>::new(domain_size).ok_or(Error::FailedToCreateEvaluationDomain),
            |generator| create_domain_with_generator::<P::ScalarField>(generator, domain_size),
        )?;

    let g2_affine_vanishing_poly = vanishing_poly_commitment_affine(g2_affine_srs, &domain);

    Ok((domain, g2_affine_vanishing_poly))
}