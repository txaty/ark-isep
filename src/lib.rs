use ark_serialize::Compress;

mod kzg;
mod error;
mod public_parameters;
mod domain;
mod prover;
mod verifier;
mod witness;
mod statement;
mod transcript;

const COMPRESS_MOD: Compress = Compress::Yes;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
