use crate::domain::{create_domain, roots_of_unity};
use crate::error::Error;
use crate::kzg::{unsafe_setup_from_tau, Kzg};
use crate::COMPRESS_MOD;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use ark_std::{One, UniformRand, Zero};
use blake2::{Blake2b512, Digest};
use std::cmp::max;
use std::collections::BTreeMap;

#[derive(Debug)]
pub struct PublicParameters<P: Pairing> {
    pub size_left_values: usize,
    pub size_right_values: usize,

    pub g1_affine_srs: Vec<P::G1Affine>,
    pub g2_affine_srs: Vec<P::G2Affine>,
    
    pub domain_l: Radix2EvaluationDomain<P::ScalarField>,
    pub domain_r: Radix2EvaluationDomain<P::ScalarField>,

    pub positions_left: Vec<usize>,
    pub positions_right: Vec<usize>,
    pub poly_positions_left: DensePolynomial<P::ScalarField>,
    pub poly_positions_right: DensePolynomial<P::ScalarField>,
    pub g1_affine_positions_left: P::G1Affine,
    pub g1_affine_positions_right: P::G1Affine,

    pub position_mappings: BTreeMap<usize, P::ScalarField>,
    pub poly_position_mappings: DensePolynomial<P::ScalarField>,
    pub g1_affine_position_mappings: P::G1Affine,

    pub domain_coset_l: Radix2EvaluationDomain<P::ScalarField>,
    pub domain_coset_r: Radix2EvaluationDomain<P::ScalarField>,

    pub coset_eval_list_positions_left: Vec<P::ScalarField>,
    pub coset_eval_list_positions_right: Vec<P::ScalarField>,
    pub coset_eval_list_position_mappings: Vec<P::ScalarField>,
    pub roots_of_unity_coset_r: Vec<P::ScalarField>,

    pub(crate) hash_representation: Vec<u8>,
}

impl<P: Pairing> PublicParameters<P> {
    pub fn builder() -> PublicParametersBuilder<P> {
        PublicParametersBuilder::<P>::default()
    }
}

pub struct PublicParametersBuilder<P: Pairing> {
    size_left_values: Option<usize>,
    size_right_values: Option<usize>,
    tau: Option<P::ScalarField>,
    domain_generator_l: Option<P::ScalarField>,
    domain_generator_r: Option<P::ScalarField>,
    position_mappings: Option<BTreeMap<usize, usize>>,
}

impl<P: Pairing> PublicParametersBuilder<P> {
    fn default() -> Self {
        Self {
            size_left_values: None,
            size_right_values: None,
            tau: None,
            domain_generator_l: None,
            domain_generator_r: None,
            position_mappings: None,
        }
    }

    pub fn size_left_values(mut self, size: usize) -> Self {
        self.size_left_values = Some(size);
        self
    }

    pub fn size_right_values(mut self, size: usize) -> Self {
        self.size_right_values = Some(size);
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
    
    pub fn position_mappings(mut self, mappings: &BTreeMap<usize, usize>) -> Self {
        self.position_mappings = Some(mappings.clone());
        self
    }

    pub fn build<R: Rng + ?Sized>(self, rng: &mut R) -> Result<PublicParameters<P>, Error> {
        let size_left_values = self.size_left_values.ok_or(Error::MissingParameter("Left \
        Element Size"))?;
        validate_input(size_left_values, None)?;
        let size_right_values = self.size_right_values.ok_or(Error::MissingParameter("Right \
        Element Size"))?;
        validate_input(size_right_values, None)?;
        let pow_of_tau_g1 = max(size_left_values, size_right_values);

        let tau = self.tau.unwrap_or(P::ScalarField::rand(rng));
        let (g1_affine_srs, g2_affine_srs) = unsafe_setup_from_tau::<P, R>(pow_of_tau_g1, tau);

        let domain_l = create_domain::<P>(self.domain_generator_l, size_left_values)?;
        let domain_r = create_domain::<P>(self.domain_generator_r, size_right_values)?;

        let position_mappings = self.position_mappings.ok_or(Error::IndexMappingCannotBeNone)?;
        let (positions_left, positions_right): (Vec<_>, Vec<_>) = position_mappings.iter()
            .map(|(&key, &value)| (key, value))
            .unzip();

        let fr_zero = P::ScalarField::zero();
        let fr_one = P::ScalarField::one();
        let mut poly_eval_positions_left = vec![fr_zero; size_left_values];
        positions_left.iter().for_each(|&i| {
            poly_eval_positions_left[i] = fr_one;
        });
        let coeff_positions_left = domain_l.ifft(&poly_eval_positions_left);
        let poly_positions_left = DensePolynomial::from_coefficients_vec(coeff_positions_left);
        let g2_affine_positions_left = Kzg::<P::G1>::commit(&g1_affine_srs, &poly_positions_left)
            .into_affine();

        let mut poly_eval_positions_right = vec![fr_zero; size_right_values];
        positions_right.iter().for_each(|&i| {
            poly_eval_positions_right[i] = fr_one;
        });
        let coeff_positions_right = domain_r.ifft(&poly_eval_positions_right);
        let poly_positions_right = DensePolynomial::from_coefficients_vec(coeff_positions_right);
        let g2_affine_positions_right = Kzg::<P::G1>::commit(&g1_affine_srs,
                                                             &poly_positions_right).into_affine();

        let mut poly_eval_position_mappings: Vec<P::ScalarField> = vec![fr_zero;
                                                                        size_left_values];
        let roots_of_unity_r = roots_of_unity::<P>(&domain_r);
        let mut fr_position_mappings = BTreeMap::new();
        position_mappings.iter().for_each(|(&key, &value)| {
            let eval = roots_of_unity_r[value];
            poly_eval_position_mappings[key] = eval;
            fr_position_mappings.insert(key, eval);
        });
        domain_l.ifft_in_place(&mut poly_eval_position_mappings);
        let coeff_position_mappings = poly_eval_position_mappings;
        let poly_position_mappings = DensePolynomial::from_coefficients_vec(coeff_position_mappings);
        let g1_affine_position_mappings = Kzg::<P::G1>::commit(&g1_affine_srs,
                                                               &poly_position_mappings).into_affine();

        let domain_coset_l = domain_l.get_coset(P::ScalarField::GENERATOR)
            .ok_or(Error::FailedToCreateCosetOfEvaluationDomain)?;
        let domain_coset_r = domain_r.get_coset(P::ScalarField::GENERATOR)
            .ok_or(Error::FailedToCreateCosetOfEvaluationDomain)?;
        let coset_eval_list_positions_left = domain_coset_l.fft(&poly_positions_left);
        let coset_eval_list_positions_right = domain_coset_r.fft(&poly_positions_right);
        let coset_eval_list_position_mappings = domain_coset_l.fft(&poly_position_mappings);
        let roots_of_unity_coset_r = roots_of_unity::<P>(&domain_coset_r);

        // Construct Hash Representation.
        let mut blake2b_hasher = Blake2b512::new();
        let mut buf = Vec::new();
        serialize_usize(size_left_values, &mut buf);
        serialize_usize(size_right_values, &mut buf);
        g2_affine_positions_left.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        g2_affine_positions_right.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        g1_affine_position_mappings.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        domain_l.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        domain_r.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        domain_coset_l.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        domain_coset_r.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        g1_affine_srs.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        g2_affine_srs.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        positions_left.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        positions_right.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);

        poly_positions_left.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        poly_positions_right.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        position_mappings.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        poly_position_mappings.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        coset_eval_list_positions_left.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        coset_eval_list_positions_right.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        coset_eval_list_position_mappings.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        roots_of_unity_coset_r.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);

        let hash_representation = blake2b_hasher.finalize().to_vec();

        Ok(PublicParameters {
            size_left_values,
            size_right_values,
            g1_affine_srs,
            g2_affine_srs,
            domain_l,
            domain_r,
            positions_left,
            positions_right,
            poly_positions_left,
            poly_positions_right,
            g1_affine_positions_left: g2_affine_positions_left,
            g1_affine_positions_right: g2_affine_positions_right,
            position_mappings: fr_position_mappings,
            poly_position_mappings,
            g1_affine_position_mappings,
            hash_representation,
            domain_coset_l,
            domain_coset_r,
            coset_eval_list_positions_left,
            coset_eval_list_positions_right,
            coset_eval_list_position_mappings,
            roots_of_unity_coset_r,
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


fn serialize_usize(input: usize, buf: &mut Vec<u8>) {
    buf.extend_from_slice(&input.to_le_bytes());
}