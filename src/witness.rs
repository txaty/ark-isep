use crate::error::Error;
use crate::public_parameters::PublicParameters;
use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain};

pub struct Witness<P: Pairing> {
    pub(crate) left_values: Vec<P::ScalarField>,
    pub(crate) right_values: Vec<P::ScalarField>,
    pub(crate) poly_left_values: DensePolynomial<P::ScalarField>,
    pub(crate) poly_right_values: DensePolynomial<P::ScalarField>,
}

impl<P: Pairing> Witness<P> {
    pub fn new(
        pp: &PublicParameters<P>,
        left_values: &[P::ScalarField],
        right_values: &[P::ScalarField],
    ) -> Result<Self, Error> {
        if left_values.len() != pp.size_left_values {
            return Err(Error::WrongNumberOfLeftValues(left_values.len()));
        }

        if right_values.len() != pp.size_right_values {
            return Err(Error::WrongNumberOfRightValues(right_values.len()));
        }

        let coeff_left_values = pp.domain_l.ifft(left_values);
        let poly_left_values = DensePolynomial::from_coefficients_vec(coeff_left_values);
        let coeff_right_values = pp.domain_r.ifft(right_values);
        let poly_right_values = DensePolynomial::from_coefficients_vec(coeff_right_values);

        Ok(Self {
            left_values: left_values.to_vec(),
            right_values: right_values.to_vec(),
            poly_left_values,
            poly_right_values,
        })
    }
}



