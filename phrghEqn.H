if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialization", false))
{
    volScalarField& ph_rgh = regIOobject::store
    (
        new volScalarField
        (
            IOobject
            (
                "ph_rgh",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    if (equal(runTime.value(), 0))
    {
        p = pMin;

        phase1.thermo().correct();
        phase2.thermo().correct();
        volScalarField rho("rho", fluid.rho());
        volScalarField rho_0("rho_0", rho);
        volScalarField delta_rho("delta_rho", rho);

        volVectorField U("U", fluid.U());
        // Limit static pressure
        p = max(ph_rgh + rho*gh + pRef, pMin);

        Info<< "min p " << min(p).value() <<
               " max p " << max(p).value() << endl;


        Info<< "min rho " << min(rho).value() <<
               " max rho " << max(rho).value() << endl;

        for (label i=0; i<30; i++)
        {
 
            phase1.thermo().correct();
            phase2.thermo().correct();
            rho = fluid.rho();
            delta_rho = rho - rho_0;
            rho = rho_0 + 0.5 * delta_rho;
            rho_0 = rho;
            // Limit static pressure
            p = max(ph_rgh + rho*gh + pRef, pMin);

            // Limit ph_rgh
            ph_rgh = p - rho*gh - pRef;

            Info<< "min p " << min(p).value() <<
                   " max p " << max(p).value() << endl;

            Info<< "min rho " << min(rho).value() <<
                   " max rho " << max(rho).value() << endl;

            Info<< "min delta_rho " << min(delta_rho).value() <<
                   " max delta_rho " << max(delta_rho).value() << endl;

        }


        // Limit ph_rgh
        ph_rgh = p - rho*gh - pRef;


        // Take the cell centres
         volVectorField centres
        (
            IOobject
            (
                "centres",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p.mesh().C()
        );
        centres.write();

        label nCorr
        (
            pimple.dict().lookupOrDefault<label>("nHydrostaticCorrectors", 5)
        );

            Info<< "min p " << min(p).value() <<
                   " max p " << max(p).value() << endl;


            Info<< "min rho " << min(rho).value() <<
                   " max rho " << max(rho).value() << endl;

        volScalarField ph_rgh_0("ph_rgh_0", ph_rgh);


        for (label i=0; i<nCorr; i++)
        {
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));

            surfaceScalarField phig
            (
                "phig",
               -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(ph_rgh, rho, U, phig, rhof);

            fvScalarMatrix ph_rghEqn
            (
                fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
            );
  
            //ph_rghEqn.relax();

            ph_rghEqn.solve();

            // Update and limit the static pressure
            p = max(ph_rgh + rho*gh + pRef, pMin);

            // Limit ph_rgh
            ph_rgh = p - rho*gh - pRef;

            phase1.thermo().correct();
            phase2.thermo().correct();
            
            // Correct p_rgh for consistency with p and the updated densities
            rho = fluid.rho();


            ph_rgh = p - rho*gh - pRef;

            ph_rgh_0 = ph_rgh;

            // phase1.thermo().correct();
            // phase2.thermo().correct();
            // rho = fluid.rho();

            Info<< "Hydrostatic pressure variation "
                << (max(ph_rgh) - min(ph_rgh)).value() << endl;

            Info<< "min p " << min(p).value() <<
                   " max p " << max(p).value() << endl;


            Info<< "min rho " << min(rho).value() <<
                   " max rho " << max(rho).value() << endl;


        }

        ph_rgh.write();

        p_rgh = ph_rgh;
    }
}
