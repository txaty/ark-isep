use crate::domain::{create_domain_with_generator, create_sub_domain};
use crate::error::Error;
use crate::kzg::unsafe_setup_from_tau;
use crate::COMPRESS_MOD;
use ark_ec::pairing::Pairing;
use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use blake2::{Blake2b512, Digest};
use std::cmp::{max, min};

#[derive(Debug)]
pub struct PublicParameters<P: Pairing> {
    pub size_left_values: usize,
    pub size_right_values: usize,
    pub size_positions: usize,

    pub g1_affine_srs: Vec<P::G1Affine>,
    pub g2_affine_srs: Vec<P::G2Affine>,

    pub domain_l: Radix2EvaluationDomain<P::ScalarField>,
    pub domain_r: Radix2EvaluationDomain<P::ScalarField>,
    pub sub_domain: Radix2EvaluationDomain<P::ScalarField>,

    pub sub_domain_coset: Radix2EvaluationDomain<P::ScalarField>,

    pub vanishing_poly_sub: DensePolynomial<P::ScalarField>,

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
    size_positions: Option<usize>,
    tau: Option<P::ScalarField>,
    domain_generator_l: Option<P::ScalarField>,
    domain_generator_r: Option<P::ScalarField>,
}

impl<P: Pairing> PublicParametersBuilder<P> {
    fn default() -> Self {
        Self {
            size_left_values: None,
            size_right_values: None,
            size_positions: None,
            tau: None,
            domain_generator_l: None,
            domain_generator_r: None,
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

    pub fn size_positions(mut self, size: usize) -> Self {
        self.size_positions = Some(size);
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

    pub fn build<R: Rng + ?Sized>(self, rng: &mut R) -> Result<PublicParameters<P>, Error> {
        let size_left_values = self.size_left_values.ok_or(Error::MissingParameter("Left \
        Element Size"))?;
        validate_input(size_left_values, None)?;

        let size_right_values = self.size_right_values.ok_or(Error::MissingParameter("Right \
        Element Size"))?;
        validate_input(size_right_values, None)?;

        let size_positions = self.size_positions.ok_or(Error::MissingParameter("Position Size"))?;
        validate_input(size_positions, Some(min(size_left_values, size_right_values)))?;

        let domain_l = create_domain::<P>(self.domain_generator_l, size_left_values)?;
        let domain_r = create_domain::<P>(self.domain_generator_r, size_right_values)?;

        let sub_domain_l = create_sub_domain::<P>(
            &domain_l,
            size_positions,
            size_left_values / size_positions,
        )?;

        let sub_domain_r = create_sub_domain::<P>(
            &domain_r,
            size_positions,
            size_right_values / size_positions,
        )?;

        if sub_domain_l != sub_domain_r {
            return Err(Error::FailedToGenerateSubdomain);
        }

        let sub_domain = sub_domain_l;

        let pow_of_tau_g1 = max(size_left_values, size_right_values);
        let tau = self.tau.unwrap_or(P::ScalarField::rand(rng));
        let (g1_affine_srs, g2_affine_srs) = unsafe_setup_from_tau::<P, R>(pow_of_tau_g1, tau);
        
        let sub_domain_coset = sub_domain.get_coset(P::ScalarField::GENERATOR)
            .ok_or(Error::FailedToCreateCosetOfEvaluationDomain)?;

        let vanishing_poly_sub: DensePolynomial<P::ScalarField> = sub_domain.vanishing_polynomial()
            .into();


        // Construct Hash Representation.
        let mut blake2b_hasher = Blake2b512::new();
        let mut buf = Vec::new();
        serialize_usize(size_left_values, &mut buf);
        serialize_usize(size_right_values, &mut buf);
        serialize_usize(size_positions, &mut buf);
        domain_l.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        domain_r.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        sub_domain.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        sub_domain_coset.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_|
            Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        g1_affine_srs.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);
        buf.clear();

        g2_affine_srs.serialize_with_mode(&mut buf, COMPRESS_MOD).map_err(|_| Error::FailedToSerializeElement)?;
        blake2b_hasher.update(&buf);

        let hash_representation = blake2b_hasher.finalize().to_vec();

        Ok(PublicParameters {
            size_left_values,
            size_right_values,
            size_positions,
            g1_affine_srs,
            g2_affine_srs,
            domain_l,
            domain_r,
            sub_domain,
            sub_domain_coset,
            vanishing_poly_sub,
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

fn create_domain<P: Pairing>(
    domain_generator: Option<P::ScalarField>,
    domain_size: usize,
) -> Result<Radix2EvaluationDomain<P::ScalarField>, Error> {
    let domain = domain_generator
        .map_or_else(
            || Radix2EvaluationDomain::<P::ScalarField>::new(domain_size).ok_or(Error::FailedToCreateEvaluationDomain),
            |generator| create_domain_with_generator::<P::ScalarField>(generator, domain_size),
        )?;

    Ok(domain)
}

fn serialize_usize(input: usize, buf: &mut Vec<u8>) {
    buf.extend_from_slice(&input.to_le_bytes());
}