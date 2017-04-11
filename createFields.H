#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating twoPhaseSystem\n" << endl;

twoPhaseSystem fluid(mesh, g);

phaseModel& phase1 = fluid.phase1();
phaseModel& phase2 = fluid.phase2();

volScalarField& alpha1 = phase1;
volVectorField& U1 = phase1.U();
volVectorField& U2 = phase2.U();

volScalarField& p = phase1.thermo().p();

word eos1(word(phase1.thermo().subDict("thermoType").lookup("equationOfState")));
word eos2(word(phase2.thermo().subDict("thermoType").lookup("equationOfState")));

if  (!( ( (equal(eos1, word("rhoConst"))) && (equal(eos2, word("perfectGas"))) ) ||
    ( (equal(eos2, word("rhoConst"))) && (equal(eos1, word("perfectGas"))) ) ) )
{
	FatalError
        << "Gas and solid phases not found" << nl
        << exit(FatalError);
}

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"
#include "readpRef.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField mu1
(
    IOobject
    (
        "mu1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phase1.thermo().mu() 
);


volScalarField mu2
(
    IOobject
    (
        "mu2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phase2.thermo().mu() 
);

Info<< "min mu1 " << min(mu1).value() <<
                   " max mu1 " << max(mu1).value() << endl;

Info<< "min mu2 " << min(mu2).value() <<
                   " max mu2 " << max(mu2).value() << endl;

volScalarField mu_temp(0.0*phase1.thermo().mu());



if ( equal(eos1, word("rhoConst")) )
{
    // phase 2 is the gas phase
	// Ishii & Mishima 1982	 
	mu_temp = phase2.thermo().mu() * pow( 1.0 - phase1 / 0.62 , -1.55);
}
else
{
    // phase 1 is the gas phase
	// Ishii & Mishima 1982	 
	mu_temp = phase1.thermo().mu() * pow( 1.0 - phase2 / 0.62 , -1.55);
}

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mu_temp 
);

mesh.setFluxRequired(p_rgh.name());


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ( phase1 * phase1.thermo().rho() + phase2 * phase2.thermo().rho() )
);



volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.U()
);

volVectorField Umix
(
    IOobject
    (
        "Umix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    ( phase1 * phase1.thermo().rho() * U1 + phase2 * phase2.thermo().rho() * U2 ) / rho
);

word kinematicCloudName("kinematicCloud");
args.optionReadIfPresent("cloudName",kinematicCloudName);

Info<< "contructing kinematicCloud " << kinematicCloudName << endl;

basicKinematicTypeCloud kinematicCloud
(
    kinematicCloudName,
    rho,
    Umix,
    mu,
    g
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);


Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K1(IOobject::groupName("K", phase1.name()), 0.5*magSqr(U1));
volScalarField K2(IOobject::groupName("K", phase2.name()), 0.5*magSqr(U2));

#include "createMRF.H"
