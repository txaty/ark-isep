#[derive(Debug)]
pub enum Error {
    MissingParameter(&'static str),
    InvalidEvaluationDomainSize(usize),
    InputShouldBePowerOfTwo(usize),
    InvalidSegmentSize(usize),
    InvalidSizesOfPositions(usize, usize),
    InputIsTooLarge(usize),
    FailedToCreateEvaluationDomain,
    FailedToInverseFieldElement,
    LeftIndicesCannotBeNone,
    RightIndicesCannotBeNone,
    IndexMappingCannotBeNone,
    WrongNumberOfLeftValues(usize),
    WrongNumberOfRightValues(usize),
    FailedToSerializeElement,
    RemainderAfterDivisionIsNonZero,
    FailedToCreateCosetOfEvaluationDomain,
    
    FailedToGenerateSubdomain,

    ValuesDoNotMatch,
    
    Check1Failed,
    Check2Failed,
    PairingFailed,
}