
rho = fluid.rho(); 

Umix = ( phase1 * phase1.thermo().rho() * U1 + phase2 * phase2.thermo().rho() * U2 ) / rho;

if ( equal(eos1, word("rhoConst")) )
{
    // phase 2 is the gas phase
	// Ishii & Mishima 1982	 
	mu = phase2.thermo().mu() * pow( 1.0 - phase1 / 0.62 , -1.55);
}
else
{
    // phase 1 is the gas phase
	// Ishii & Mishima 1982	 
	mu = phase1.thermo().mu() * pow( 1.0 - phase2 / 0.62 , -1.55);
}

Info<< "min Umix " << min(Umix).value() <<
       " max Umix " << max(Umix).value() << endl;


Info<< "min mu " << min(mu).value() <<
       " max mu " << max(mu).value() << endl;


kinematicCloud.evolve();

