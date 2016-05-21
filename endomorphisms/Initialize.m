// Global constants; should perhaps just copy these
prec0 := 200;
epscomp0 := 10^(-prec0 + 30);
epsLLL0 := 5^(-prec0 + 7);
epsinv0 := 2^(-4);

// FIXME: Temporary bug fix
forward AutomorphismGroupFixed;
forward FixedGroupFixed;

// Linear algebra routines
forward NumericalSolve;
forward InvertibleSubmatrix;
forward SplitPeriodMatrix;
forward CombinePeriodMatrix;
forward IntegralKernel;
//forward OptimizedPeriodMatrix;
forward ConjugateMatrix;
forward MatrixInBasis;
forward MatrixRatInBasisOverNF;

// Functionality for recognizing complex numbers as algebraic numbers and
// related optimizations
forward PolynomializeElement;
forward PolynomializeMatrix;
forward FractionalApproximation;
forward AlgebraizeElementInField;
forward NearbyRoot;
forward AlgebraizeMatricesInField;
forward IntegralRepresentationNF;
//forward PartialLMFDBLabel;

// Finding an approximate basis for the geometric endomorphism ring through LLL
forward ComplexStructure;
forward RationalEndomorphismEquations;
forward AnalyticRepresentation;
forward GeometricEndomorphismBasisFromPeriods;
forward GeometricEndomorphismBasis;

// Algebraizing the basis
forward EndomorphismBasisOverSubfield;
forward EndomorphismLatticeG2;
forward EndomorphismData;
forward EndomorphismAlgebraFactorsGeneric;
forward EndomorphismRingGeneric;
forward EndomorphismAlgebraFactorsQQG2;
forward EndomorphismAlgebraFactorsRRG2;
forward EndomorphismRingG2;
forward SatoTateGroupG2;
forward SplittingIdempotentsG2;
forward GeoEndRRShorthand;

// Determining the Sato-Tate group
forward SatoTateGroupG2;

// Elliptic curves from splittings
forward SmallIdempotent;
forward InducedEmbedding;
forward LatticesFromIdempotents;
//forward EllipticCurveFromEisensteinValues;

// The actual loading
load "BugFixes.m";
load "Linear.m";
load "Recognition.m";
load "Analytic.m";
load "Algebraic.m";
load "SatoTate.m";
load "Splitting.m";
