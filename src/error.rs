#[derive(Debug)]
pub enum Error {
    MissingParameter(&'static str),
    InvalidEvaluationDomainSize(usize),
    InputShouldBePowerOfTwo(usize),
    InputIsTooLarge(usize),
    FailedToCreateEvaluationDomain,
    FailedToInverseFieldElement,
    LeftIndicesCannotBeNone,
    RightIndicesCannotBeNone,
    IndexMappingCannotBeNone,
    WrongNumberOfLeftElements(usize),
    WrongNumberOfRightElements(usize),
    FailedToSerializeElement,
    RemainderAfterDivisionIsNonZero,
    FailedToCreateCosetOfEvaluationDomain,

    Pairing1Failed,
    Pairing2Failed,
    Pairing3Failed,
    EqualityCheckFailed,
}