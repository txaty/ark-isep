use crate::domain::{create_domain_with_generator, roots_of_unity, vanishing_poly_commitment_affine};
use crate::error::Error;
use crate::kzg::{unsafe_setup_from_tau, Kzg};
use crate::COMPRESS_MOD;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use ark_std::{One, UniformRand, Zero};
use blake2::{Blake2b512, Digest};
use rayon::prelude::*;
use std::cmp::max;
use std::collections::BTreeMap;

#[derive(Debug)]
pub struct PublicParameters<P: Pairing> {
    pub(crate) left_element_size: usize,
    pub(crate) right_element_size: usize,
    position_vec_size: usize,

    pub g1_affine_srs: Vec<P::G1Affine>,
    pub g2_affine_srs: Vec<P::G2Affine>,

    pub g2_affine_zl: P::G2Affine,
    pub g2_affine_zr: P::G2Affine,
    pub g2_affine_zi: P::G2Affine,

    pub domain_l: Radix2EvaluationDomain<P::ScalarField>,
    pub domain_r: Radix2EvaluationDomain<P::ScalarField>,
    pub domain_i: Radix2EvaluationDomain<P::ScalarField>,

    pub positions_left: Vec<usize>,
    pub positions_right: Vec<usize>,
    pub poly_positions_left: DensePolynomial<P::ScalarField>,
    pub poly_positions_right: DensePolynomial<P::ScalarField>,
    pub g2_affine_positions_left: P::G2Affine,
    pub g2_affine_positions_right: P::G2Affine,

    pub position_mappings: BTreeMap<usize, usize>,
    pub poly_position_mappings: DensePolynomial<P::ScalarField>,
    pub g1_affine_position_mappings: P::G1Affine,

    pub hash_representation: Vec<u8>,
}

impl<P: Pairing> PublicParameters<P> {
    pub fn builder() -> PublicParametersBuilder<P> {
        PublicParametersBuilder::<P>::default()
    }
}

pub struct PublicParametersBuilder<P: Pairing> {
    left_element_size: Option<usize>,
    right_element_size: Option<usize>,
    position_vec_size: Option<usize>,
    tau: Option<P::ScalarField>,
    domain_generator_l: Option<P::ScalarField>,
    domain_generator_r: Option<P::ScalarField>,
    domain_generator_i: Option<P::ScalarField>,
    positions_left: Option<Vec<usize>>,
    positions_right: Option<Vec<usize>>,
    position_mappings: Option<BTreeMap<usize, usize>>,
}

impl<P: Pairing> PublicParametersBuilder<P> {
    fn default() -> Self {
        Self {
            left_element_size: None,
            right_element_size: None,
            position_vec_size: None,
            tau: None,
            domain_generator_l: None,
            domain_generator_r: None,
            domain_generator_i: None,
            positions_left: None,
            positions_right: None,
            position_mappings: None,
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

    pub fn position_vec_size(mut self, size: usize) -> Self {
        self.position_vec_size = Some(size);
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

    pub fn positions_left(mut self, indices: &[usize]) -> Self {
        self.positions_left = Some(indices.to_vec());
        self
    }

    pub fn positions_right(mut self, indices: &[usize]) -> Self {
        self.positions_right = Some(indices.to_vec());
        self
    }

    pub fn position_mappings(mut self, mappings: &BTreeMap<usize, usize>) -> Self {
        self.position_mappings = Some(mappings.clone());
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
        let position_vec_size = self.position_vec_size.ok_or(Error::MissingParameter("Index \
        Size"))?;
        validate_input(position_vec_size, Some(pow_of_tau_g1))?;

        let tau = self.tau.unwrap_or(P::ScalarField::rand(rng));
        let (g1_affine_srs, g2_affine_srs) = unsafe_setup_from_tau::<P, R>(pow_of_tau_g1, tau);

        let (domain_l, g2_affine_zl) = create_domain_and_vanishing_poly_commitment::<P>(self
                                                                                            .domain_generator_l, left_element_size, &g2_affine_srs)?;
        let (domain_r, g2_affine_zr) = create_domain_and_vanishing_poly_commitment::<P>(self
                                                                                            .domain_generator_r, right_element_size, &g2_affine_srs)?;
        let (domain_i, g2_affine_zi) = create_domain_and_vanishing_poly_commitment::<P>(self
                                                                                            .domain_generator_i, position_vec_size, &g2_affine_srs)?;

        let positions_left = self.positions_left.ok_or(Error::LeftIndicesCannotBeNone)?;
        let positions_right = self.positions_right.ok_or(Error::RightIndicesCannotBeNone)?;
        let position_mappings = self.position_mappings.ok_or(Error::IndexMappingCannotBeNone)?;

        let mut poly_eval_positions_left = vec![P::ScalarField::zero(); left_element_size];
        // TODO: Optimize this.
        positions_left.iter().for_each(|&i| {
            poly_eval_positions_left[i] = P::ScalarField::one();
        });
        let coeff_positions_left = domain_i.ifft(&poly_eval_positions_left);
        let poly_positions_left = DensePolynomial::from_coefficients_vec(coeff_positions_left);
        let g2_affine_positions_left = Kzg::<P::G2>::commit(&g2_affine_srs, &poly_positions_left)
            .into_affine();

        let mut poly_eval_positions_right = vec![P::ScalarField::zero(); right_element_size];
        // TODO: Optimize this.
        positions_right.iter().for_each(|&i| {
            poly_eval_positions_right[i] = P::ScalarField::one();
        });
        let coeff_positions_right = domain_i.ifft(&poly_eval_positions_right);
        let poly_positions_right = DensePolynomial::from_coefficients_vec(coeff_positions_right);
        let g2_affine_positions_right = Kzg::<P::G2>::commit(&g2_affine_srs,
                                                             &poly_positions_right).into_affine();

        let mut poly_eval_position_mappings: Vec<P::ScalarField> = vec![P::ScalarField::zero();
                                                                        left_element_size];
        let roots_of_unity_r = roots_of_unity::<P>(&domain_r);
        // TODO: Optimize this.
        position_mappings.iter().for_each(|(&key, &value)| {
            poly_eval_position_mappings[key] = roots_of_unity_r[value];
        });
        domain_l.ifft_in_place(&mut poly_eval_position_mappings);
        let coeff_position_mappings = poly_eval_position_mappings;
        let poly_position_mappings = DensePolynomial::from_coefficients_vec(coeff_position_mappings);
        let g1_affine_position_mappings = Kzg::<P::G1>::commit(&g1_affine_srs,
                                                               &poly_position_mappings).into_affine();

        // Construct Hash Representation.
        let mut buf = Vec::new();
        serialize_usize(left_element_size, &mut buf);
        serialize_usize(right_element_size, &mut buf);
        serialize_usize(position_vec_size, &mut buf);

        g1_affine_srs.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        g2_affine_srs.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        g2_affine_zl.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        g2_affine_zr.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        g2_affine_zi.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;

        domain_l.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        domain_r.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        domain_i.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;

        positions_left.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        positions_right.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        poly_positions_left.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        poly_positions_right.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        g2_affine_positions_left.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        g2_affine_positions_right.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        position_mappings.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        poly_position_mappings.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        g1_affine_position_mappings.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;

        let mut blake2b_hasher = Blake2b512::new();
        blake2b_hasher.update(&buf);
        let hash_representation = blake2b_hasher.finalize().to_vec();


        Ok(PublicParameters {
            left_element_size,
            right_element_size,
            position_vec_size,
            g1_affine_srs,
            g2_affine_srs,
            g2_affine_zl,
            g2_affine_zr,
            g2_affine_zi,
            domain_l,
            domain_r,
            domain_i,
            positions_left,
            positions_right,
            poly_positions_left,
            poly_positions_right,
            g2_affine_positions_left,
            g2_affine_positions_right,
            position_mappings,
            poly_position_mappings,
            g1_affine_position_mappings,
            hash_representation,
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

    let g2_affine_vanishing_poly = vanishing_poly_commitment_affine::<P::G2>(g2_affine_srs,
                                                                             &domain);

    Ok((domain, g2_affine_vanishing_poly))
}

fn serialize_usize(input: usize, buf: &mut Vec<u8>) {
    buf.extend_from_slice(&input.to_le_bytes());
}