use crate::error::Error;
use crate::public_parameters::PublicParameters;
use ark_ec::pairing::Pairing;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain};

pub struct Witness<P: Pairing> {
    pub(crate) left_elements: Vec<P::ScalarField>,
    pub(crate) right_elements: Vec<P::ScalarField>,
    pub(crate) poly_left_elements: DensePolynomial<P::ScalarField>,
    pub(crate) poly_right_elements: DensePolynomial<P::ScalarField>,
}

impl<P: Pairing> Witness<P> {
    pub fn new(
        pp: &PublicParameters<P>,
        left_elements: &[P::ScalarField],
        right_elements: &[P::ScalarField],
    ) -> Result<Self, Error> {
        if left_elements.len() != pp.size_left_values {
            return Err(Error::WrongNumberOfLeftElements(left_elements.len()));
        }

        if right_elements.len() != pp.size_right_values {
            return Err(Error::WrongNumberOfRightElements(right_elements.len()));
        }

        let coeff_left_elements = pp.domain_l.ifft(left_elements);
        let poly_left_elements = DensePolynomial::from_coefficients_vec(coeff_left_elements);
        let coeff_right_elements = pp.domain_r.ifft(right_elements);
        let poly_right_elements = DensePolynomial::from_coefficients_vec(coeff_right_elements);

        Ok(Self {
            left_elements: left_elements.to_vec(),
            right_elements: right_elements.to_vec(),
            poly_left_elements,
            poly_right_elements,
        })
    }
}



