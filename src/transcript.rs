use crate::error::Error;
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use merlin::Transcript as MerlinTranscript;
use std::marker::PhantomData;

/// Modified from https://github.com/caulk-crypto/caulk/blob/main/src/transcript.rs

#[derive(Copy, Clone)]
pub(crate) enum Label {
    ChallengeBeta,
    // ChallengeDelta,
    ChallengeGamma,
    // ChallengeEta,
    // ChallengeCaulkAlpha,
    // ChallengeCaulkBeta,

    PublicParameters,
    Statement,

    // G1M,
    // G1MDivW,
    // G1Qm,
    // G1L,
    // G1LDivV,
    // G1Ql,
    // G1D,
    // G1Qd,
    // G1A,
    // G1Qa,
    // G1Qb,
    // G1A0,
    // G1B0,
    // G1Px,
    // 
    // FrB0AtGamma,
    // FrFAtGamma,
    // FrLAtGamma,
    // FrAAtZero,
    // FrLAtGammaDivV,
    // FrQlAtGamma,
    // FrDAtGamma,
    // FrQdAtGamma,
    // 
    // CaulkG1D,
    // CaulkG1UBar,
    // CaulkG1H2,
    // CaulkG1H1,
}

impl Label {
    pub fn as_bytes(&self) -> &'static [u8] {
        match self {
            Label::ChallengeBeta => b"beta",
            Label::ChallengeGamma => b"gamma",
            Label::PublicParameters => b"common_inputs",
            Label::Statement => b"statement",
        }
    }
}

pub(crate) struct Transcript<F: PrimeField> {
    transcript: MerlinTranscript,
    _marker: PhantomData<F>,
}

impl<F: PrimeField> Default for Transcript<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: PrimeField> Transcript<F> {
    pub(crate) fn new() -> Self {
        Self {
            transcript: MerlinTranscript::new(b"Init SegLookup Transcript"),
            _marker: PhantomData::default(),
        }
    }

    /// Get a uniform random field element for field size < 384
    pub(crate) fn get_and_append_challenge(&mut self, label: Label) -> Result<F, Error> {
        let mut bytes = [0u8; 64];
        self.transcript
            .challenge_bytes(label.as_bytes(), &mut bytes);
        let challenge = F::from_le_bytes_mod_order(bytes.as_ref());
        self.append_element(label, &challenge)?;

        Ok(challenge)
    }

    /// Append a field/group element to the transcript
    pub(crate) fn append_element<T: CanonicalSerialize>(
        &mut self,
        label: Label,
        element: &T,
    ) -> Result<(), Error> {
        let mut buf = vec![];
        element
            .serialize_uncompressed(&mut buf)
            .map_err(|_| Error::FailedToSerializeElement)?;
        self.transcript
            .append_message(label.as_bytes(), buf.as_ref());

        Ok(())
    }

    pub(crate) fn append_elements<T: CanonicalSerialize>(
        &mut self,
        labels_and_elements: &[(Label, T)],
    ) -> Result<(), Error> {
        for (label, element) in labels_and_elements {
            self.append_element(*label, element)?;
        }

        Ok(())
    }
}
