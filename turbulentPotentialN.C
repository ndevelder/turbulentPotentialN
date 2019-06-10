/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "turbulentPotentialN.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "components.H"
#include "fvCFD.H"
#include "volFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Hello

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentPotentialN, 0);
addToRunTimeSelectionTable(RASModel, turbulentPotentialN, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> turbulentPotentialN::Ts() const
{ 
	if(tslimiter_ == "true")
	{
        return max(k_/(epsilon_ + epsilonSmall_), 6.0*sqrt(nu()/(epsilon_ + epsilonSmall_)));
	}
	
    return ((k_+k0_)/(epsilon_ + epsilonSmall_));
}


tmp<volScalarField> turbulentPotentialN::Ls() const
{
	
	volScalarField trueL = pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
	
	if(lslimiter_ == "true")
	{
		return cL1_*max(pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_),cL2_*pow(pow3(nu())/(epsilon_ + epsilonSmall_),0.25));
	}
	
	//Info << "Max trueL: " << gMax(trueL) << " Min trueL: " << gMin(trueL) << endl;
	//Info << "Dims: " << trueL.dimensions() << endl;
	
	return pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPotentialN::turbulentPotentialN
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),


    cEp1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.45
        )
    ),
    cEp2con_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp2con",
            coeffDict_,
            1.83
        )
    ),
    cEp3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp3",
            coeffDict_,
            0.15
        )
    ),
    cEp4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp4",
            coeffDict_,
            0.042
        )
    ),
    cD1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.5
        )
    ),
    cD2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            0.33
        )
    ),
    cD3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD3",
       	    coeffDict_,
            3.0
        )
    ),
    cD4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD4",
       	    coeffDict_,
            2.4
        )
    ),
    cP1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            2.0
        )
    ),
    cP1beta_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1beta",
            coeffDict_,
            0.0
        )
    ),
    cP1chi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1chi",
            coeffDict_,
            0.67
        )
    ),
    cP1pow_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1pow",
            coeffDict_,
            1.0
        )
    ),
    cP1type_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1type",
            coeffDict_,
            1.0
        )
    ),
    cP2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.6
        )
    ),
    cP3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.42
        )
    ),
    cP4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.85714
        )
    ),
    cP5_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP5",
            coeffDict_,
            0.02
        )
    ),
    cL1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL1",
            coeffDict_,
            0.36
        )
    ),
    cL2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL2",
            coeffDict_,
            85.0
        )
    ),
    cMu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.21
        )
    ),
    betaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaK",
            coeffDict_,
            0.09
        )
    ),
    cT_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cT",
            coeffDict_,
            0.02
        )
    ),
    cA_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cA",
            coeffDict_,
            1.0
        )
    ),
    cG_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cG",
            coeffDict_,
            0.07
        )
    ),
    cGw_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cGw",
            coeffDict_,
            0.67
        )
    ),
    cGp_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cGp",
            coeffDict_,
            1.0
        )
    ),
    cEhmM_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmM",
            coeffDict_,
            10.0
        )
    ),
	cNF_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cNF",
            coeffDict_,
            1.0
        )
    ),
	cN1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cN1",
            coeffDict_,
            0.6
        )
    ),
	cN2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cN2",
            coeffDict_,
            2.2
        )
    ),
	cND1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cND1",
            coeffDict_,
            0.5
        )
    ),
	cND2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cND2",
            coeffDict_,
            0.5
        )
    ),
	pMix_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "pMix",
            coeffDict_,
            0.4
        )
    ),
	cPrK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrK",
            coeffDict_,
            0.6
        )
    ),
	cPrP_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrP",
            coeffDict_,
            1.0
        )
    ),
	cNL_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cNL",
            coeffDict_,
            0.00001
        )
    ),
	nutRatMax_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutRatMax",
            coeffDict_,
            1.0e5
        )
    ),
    sigmaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            0.833
        )
    ),
    sigmaPhi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPhi",
            coeffDict_,
            0.33
        )
    ),
    sigmaPhiC_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPhiC",
            coeffDict_,
            0.7
        )
    ),
    sigmaPsi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPsi",
            coeffDict_,
            1.0
        )
    ),
    sigmaPsiC_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPsiC",
            coeffDict_,
            0.5
        )
    ),
	prodType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"prodType",
			coeffDict_,
			1.0
		)
	),
	dampType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"dampType",
			coeffDict_,
			1.0
		)
	),
	cEp1Type_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"cEp1Type",
			coeffDict_,
			1.0
		)
	),
	gType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"gType",
			coeffDict_,
			1.0
		)
	),
	nutType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"nutType",
			coeffDict_,
			1.0
		)
	),
	
   solveK_
   (
       coeffDict_.lookup("solveK")
   ),

   solveEps_
   (
       coeffDict_.lookup("solveEps")
   ),

   solvePsi_
   (
       coeffDict_.lookup("solvePsi")
   ),

   solvePhi_
   (
       coeffDict_.lookup("solvePhi")
   ),

   solveNut_
   (
       coeffDict_.lookup("solveNut")
   ),

   eqncEp1_
   (
       coeffDict_.lookup("eqncEp1")
   ),
   
   eqncEp2_
   (
       coeffDict_.lookup("eqncEp2")
   ),

   eqnEpsHat_
   (
       coeffDict_.lookup("eqnEpsHat")
   ),

   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   tslimiter_
   (
       coeffDict_.lookup("tslimiter")
   ),
   lslimiter_
   (
       coeffDict_.lookup("lslimiter")
   ),
   eqncMu_
   (
       coeffDict_.lookup("eqncMu")
   ),
   nutReal_
   (
       coeffDict_.lookup("nutReal")
   ),
   y_
   (
   mesh_
   ),
    
	k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(k_))
    ),
    
	epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	tpphi_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
	
	tpphiSqrt_
    (
        IOobject
        (
            "tpphiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(tpphi_))
    ),
	
	ndsPhi_
    (
        IOobject
        (
            "ndsPhi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
	
	ndsPsi_
    (
        IOobject
        (
            "ndsPsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
	
	vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::curl(U_))
    ),
    
	tppsi_
    (
        IOobject
        (
            "tppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::grad(U_))
    ),
	
	epsHat_
    (
        IOobject
        (
            "epsHat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (epsilon_/(k_ + k0_))
    ),
    
	kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(k_))
    ),
	
	alpha_
    (
        IOobject
        (
            "alpha",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (1.0/(1.0 + 1.5*tpphi_))
    ),
	
	gamma_
    (
        IOobject
        (
            "gamma",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (1.0/(1.0 + 0.667*nut_/nu()))
    ),
	
	tpchi_
    (
        IOobject
        (
            "tpchi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        2.0*alpha_*pow(sqr(tpphi_) + (tppsi_ & tppsi_) + SMALL,0.5)
    ),
	
	tpupsilon_
    (
        IOobject
        (
            "tpupsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        2.0 - tpphi_ - tpchi_
    ),
	
	phiSqrt_
    (
        IOobject
        (
            "phiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(tpphi_*k_))
    ),
    
	gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(kSqrt_))
    ),
    
	cEp2_
    (
        IOobject
        (
            "cEp2",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (cEp2con_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
    
	tpProd_
    (
        IOobject
        (
            "tpProd",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((2*nut_*magSqr(symm(fvc::grad(U_)))/k_))
    ),
    
	cP1eqn_
    (
        IOobject
        (
            "cP1eqn",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (2.0*(0.5+0.5*((tpProd_*k_)/epsilon_)))
    ),
    
	gradTpphi_
    (
        IOobject
        (
            "gradTpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(tpphi_))
    ),
    
	
	divIk_
    (
        IOobject
        (
            "divIk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad((2.0/3.0)*k_))
    ),
	
	gradTppsi_
    (
        IOobject
        (
            "gradTppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(tppsi_))
    ),
    
	tpProdSqr_
    (
        IOobject
        (
            "tpProdSqr",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqr(tppsi_ & vorticity_))
    ),
    
	tpProd3d_
    (
        IOobject
        (
            "tpProd3d",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (mag(psiReal() ^ vorticity_))
    )
{

    //Info<< "Made it past constructors " << endl;
	
	
    //*************************************//	
    // Eddy viscosity - diffusion only
    //*************************************//
	nut_ = cMu_*k_*tpphi_*Ts();	
	nut_ = min(nut_,nutRatMax_*nu());        
	nut_.correctBoundaryConditions();
	bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), SMALL));       
	
	
    //*************************************//	
    // Epsilon-hat
    //*************************************//   
	epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	
	
    //Info<< "solveK is: " <<solveK_ <<endl;
    //Info<< "solveEps is: " <<solveEps_ <<endl;
    //Info<< "solvePhi is: " <<solvePhi_ <<endl;
    //Info<< "solvePsi is: " <<solvePsi_ <<endl;

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialN::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_))
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialN::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

// Term that is directly added to the momentum equation
tmp<fvVectorMatrix> turbulentPotentialN::divDevReff() const
{
    return
    (
       fvc::grad(phiReal())
     + fvc::curl(psiReal())
     + fvc::laplacian(nut_, U_, "laplacian(nuEff,U_)")
     - fvm::laplacian(nuEff(), U_)
    );
}


bool turbulentPotentialN::read()
{
    if (RASModel::read())
    {
        cEp1_.readIfPresent(coeffDict());
        cEp2con_.readIfPresent(coeffDict());
        cEp3_.readIfPresent(coeffDict());
		cEp4_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
		cP1beta_.readIfPresent(coeffDict());
		cP1chi_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
		cP4_.readIfPresent(coeffDict());
		cL1_.readIfPresent(coeffDict());
		cL2_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		cEhmM_.readIfPresent(coeffDict());
		cPrK_.readIfPresent(coeffDict());
		cPrP_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
		cD3_.readIfPresent(coeffDict());
		cD4_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		cA_.readIfPresent(coeffDict());
		cG_.readIfPresent(coeffDict());
		cGw_.readIfPresent(coeffDict());
		cNL_.readIfPresent(coeffDict());
		cNF_.readIfPresent(coeffDict());
		cN1_.readIfPresent(coeffDict());
		cN2_.readIfPresent(coeffDict());
		cND1_.readIfPresent(coeffDict());
		cND2_.readIfPresent(coeffDict());
		sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaPhi_.readIfPresent(coeffDict());
		sigmaPsi_.readIfPresent(coeffDict());
		prodType_.readIfPresent(coeffDict());
		cEp1Type_.readIfPresent(coeffDict());
		gType_.readIfPresent(coeffDict());
		nutType_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void turbulentPotentialN::correct()
{

    //**********************************************//	
    // Bounding values not already defined by model
    //**********************************************//
	
    const dimensionedScalar eH0("minEpsHat", epsHat_.dimensions(), SMALL);
	const dimensionedScalar nut0("minNut", nut_.dimensions(), SMALL);
	const dimensionedScalar nutSmall("smallNut", nut_.dimensions(), SMALL);
	const dimensionedScalar tph0("minTpphi", tpphi_.dimensions(), SMALL);
	const dimensionedScalar L0("lMin", dimensionSet(0,1,0,0,0,0,0), SMALL);

    if (mesh_.changing())
    {
        y_.correct();
        bound(k_, k0_);
        bound(epsilon_, epsilonSmall_);
		bound(tpphi_,tph0);
		bound(nut_,nut0);
    }
	
	
    RASModel::correct();

	
    if (!turbulence_)
    {
        return;
    }
	
	
    //*************************************//	
    // Timestep - for use in elliptic switch
    //*************************************//
	
    //dimensionedScalar cTime = U_.mesh().time().value();
	//word sMMdebug = runTime_.controlDict().lookup("showMaxMin");

    //*************************************//	
    // Vorticity and Gradient
    //*************************************//
    
	vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_);

    //*************************************//	
    // Length and Time Scales
    //*************************************//	
	
	//const volScalarField L("Length",Ls());
	//const volScalarField L2("Lsqr",sqr(L));
	volScalarField T("Time",Ts());	
		
	//const volScalarField nutPsi("nutPsi", mag(tppsi_)*k_/(mag(vorticity_) + (cNL_/Ts())));

	
	//*************************************//	
    // Misc Terms
    //*************************************//

	const volVectorField gradPhi_("gradPhi", fvc::grad(phiReal()));		
	//const volScalarField gradgradPhi_("gradgradPhi", fvc::laplacian(DphiEff(),phiReal()));

	tpphiSqrt_ = sqrt(tpphi_ + ROOTVSMALL);
	const volVectorField gradTpphiSqrt("gradTpphiSqrt",fvc::grad(tpphiSqrt_));
	  
	
    //gradTpphi_ = fvc::grad(tpphi_);
    //phiSqrt_ = sqrt(tpphi_*k_ + k0_);
	//const volVectorField gradPhiSqrt_("gradPhiSqrt",fvc::grad(phiSqrt_));	
	
	

	

	//*************************************//	
    // Update Alpha
    //*************************************//
    
	alpha_ = 1.0/(1.0 + 1.5*tpphi_);
    volScalarField beta("beta", min(alpha_/alpha_,7.0*sqr(1.0-alpha_)));
	
	//volScalarField IIb("IIb", alpha_*(2.0*alpha_-1.0));	
    volScalarField IIb("IIb", sqr(2.0*alpha_-1.0) + 2.0*(tppsi_ & tppsi_));
	bound(IIb, SMALL);
	
	volScalarField phiActual("phiActual",tpphi_*k_);
	volVectorField psiActual("psiActual",tppsi_*k_);
	
	volScalarField nutExact("nutExact", mag(psiActual)/(mag(vorticity_) + (cNL_/Ts())));

	//volScalarField gammaNut("gammaNut", (alpha_*(psiActual & psiActual) + 0.57*(1.0-alpha_)*phiActual*phiActual)/(epsHat_*k_));
	volScalarField gammaNut("gammaNut", pow(IIb, 0.5)*nutExact + (1.0-pow(IIb, 0.5))*cMu_*phiActual/epsHat_);
	volScalarField gammaWall("gammaWall", 3.0*nu()*(gradTpphiSqrt & gradTpphiSqrt)*k_/epsilon_); 

	//gamma_ = 1.0/(1.0 + cG_*pow(gammaNut/nu(),cGp_) + cGw_*gammaWall);

    gamma_ = 1.0/(1.0 + cG_*(1.0-alpha_)*sqrt(nut_/nu()));
	
	

	//*************************************//	
    // Anisotropy and Invariants
    //*************************************//	
	tpchi_ = 2.0*alpha_*pow(sqr(tpphi_) + (tppsi_ & tppsi_) + SMALL,0.5);
	tpupsilon_ = min(2.0 - tpphi_ - tpchi_, 2.0);
	
	bound(tpchi_, SMALL);
	bound(tpupsilon_, SMALL);
	
	volScalarField bphi("bphi", tpphi_ - (2.0/3.0));
	volScalarField bchi("bchi", tpchi_ - (2.0/3.0));
	volScalarField bups("bups", tpupsilon_ - (2.0/3.0));  
	
	//volScalarField D("D", (27.0/8.0)*(tpupsilon_*tpphi_*tpchi_ - tpchi_*(tppsi_ & tppsi_)));
	//volScalarField II("II", 0.5*(sqr(tpupsilon_) + sqr(tpphi_) + sqr(tpchi_) + 2.0*(tppsi_ & tppsi_)));
	//volScalarField IIb("IIb", 0.5*(sqr(bups) + sqr(bphi) + sqr(bchi) + 2.0*(tppsi_ & tppsi_)));

	
	
    //*************************************//	
    // K Production
    //*************************************//
 
	const volScalarField S2 = 2*magSqr(dev(symm(uGrad_)));
	const volScalarField magS = sqrt(S2);
	volScalarField G("RASModel::G", nut_*S2); 
	volScalarField GdK("GdK", G/(k_ + k0_));
	

	if(prodType_.value() == 1.0){
		tpProd_ = mag(tppsi_ & vorticity_);
		G = tpProd_*k_;	
		GdK = tpProd_;
	}else if(prodType_.value() == 1.1){
        tpProd_ = tppsi_ & vorticity_;
        G = tpProd_*k_;
        GdK = tpProd_;  
    }else if(prodType_.value() == 2.0){
		tpProd_ = GdK;
	}else if(prodType_.value() == 3.0){
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + pMix_*(1.0-alpha_)*cPrK_*magS + (1.0 - pMix_)*(1.0 - alpha_)*cPrP_*tpphi_*magS;
		G = tpProd_*k_;
		GdK = tpProd_;	
    }else if(prodType_.value() == 4.0){
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + 0.94*(1.0-alpha_)*GdK;
		G = tpProd_*k_;
		GdK = tpProd_;	
    }else if(prodType_.value() == 5.0){
		tpProd_ = alpha_*mag(tppsi_ & vorticity_) + (1.0-alpha_)*(0.2 + 0.4*(1.0-sqrt(IIb)))*magS;
		G = tpProd_*k_;
		GdK = tpProd_;	
    }else{
		tpProd_ = mag(tppsi_ & vorticity_);
		G = tpProd_*k_;
		GdK = tpProd_;		
	}
    
	// tpProdSqr_ = sqr(tpProd_);
	// tpProd3d_ = mag(psiReal() ^ vorticity_);
	
	
	
    //*************************************//
    //Dissipation equation
    //*************************************//


    tmp<fvScalarMatrix> epsEqn   
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
       (cEp1_ + cEp4_*(2.0*alpha_-1.0))*G*epsHat_
     - fvm::Sp(cEp2con_*epsHat_,epsilon_)
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_,epsilonSmall_);
	
	
	
	
    //*************************************//
    // Turbulent kinetic energy equation
    //*************************************//
    
    tmp<fvScalarMatrix> kEqn
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_/(k_+k0_),k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_,k0_);
 
	
	
	kSqrt_ = sqrt(mag(k_)+k0_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), sqrt(ROOTVSMALL)));

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);
	
	
	
	//*************************************//	
    // Epsilon-hat
    //*************************************//
    
    epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)) + k0_);
    bound(epsHat_,eH0);

		
	
    //*************************************//
    // Fix for separated flows and transition
    //*************************************// 	

	volScalarField cTexp("cTexp", cT_*(1.0-gamma_)*sqrt((((nu()/100.0)+nut_)/nu())));
    volScalarField transPhi("transPhi", cTexp*cA_*((2.0/3.0) - tpphi_)*tpProd_);	
	volVectorField transPsi("transPsi", cTexp*((1.0 - alpha_)*vorticity_ - cA_*tppsi_*tpProd_));	
	

	
	//*************************************//
    // Phi/K equation 
    //*************************************//

	cP1eqn_ = (cP1_-1.0)*(1.0-gamma_);

    tmp<fvScalarMatrix> tpphiEqn  
    ( 
        fvm::ddt(tpphi_)
      + fvm::div(phi_, tpphi_)
      + fvm::SuSp(-fvc::div(phi_), tpphi_)
      - fvm::laplacian(DphiEff(), tpphi_) 
      ==   
	  // Pressure Strain Slow 
	    cP1eqn_*((2.0/3.0) - tpphi_)*epsHat_
	  + cP2_*tpProd_*tpphi_

	  // Non-linear pressure strain
	  + cD3_*nutFrac()*(sqr(bphi) + (tppsi_ & tppsi_) - (2.0/3.0)*IIb)*epsHat_
	  
	  // From K eqn + pressure strain fast
      - fvm::Sp(tpProd_,tpphi_)
	  
	  // Dissipation (Commented as taken into account in cP1eqn)
	  
	  // Transition
      + transPhi
    );

    tpphiEqn().relax();
    solve(tpphiEqn);  
    bound(tpphi_,tph0);
 	

    //*************************************//     
    // Psi Equation
    //*************************************//
      
	volScalarField uudk("uudk", (4.0/3.0)*(1.75*alpha_-0.375));
	  
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::Sp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_)   

      ==  

	  // Production
	    tpphi_*vorticity_
		
	  // Slow Pressure Strain
      - fvm::Sp(cP1_*(1.0 - beta*gamma_)*epsHat_,tppsi_)
	  + cD4_*nutFrac()*(bphi + (uudk-(2.0/3.0)))*tppsi_*epsHat_

	  // Fast Pressure Strain
	  - cP2_*tpphi_*vorticity_ 
	  + cP2_*tpProd_*tppsi_ 

	  
	  - fvm::Sp((cD2_*alpha_)*tpProd_,tppsi_)
      - cP3_*(1.0 - sqrt(IIb))*vorticity_

	  //+ cP4_*pow(nut_/nu()+SMALL,0.5)*tpProd_*tppsi_  
      //+ cP4_*sqr(1.0-gamma_)*tpProd_*tppsi_  
	  + cP4_*(1.0-gamma_)*pow(nut_/nu()+SMALL,0.5)*(1.0-alpha_)*vorticity_
	  //+ cP4_*(1.0-alpha_)*pow(nut_/nu()+SMALL,0.5)*tpProd_*tppsi_
	  
	  // From K Equation
      - fvm::Sp(tpProd_,tppsi_)  
	  
	  // Dissipation 
	  // Dissipation 
	  + (1.0-beta*gamma_)*tppsi_*epsHat_  
	  //+ cD1_*(2.0*alpha_-1.0)*tpphi_*vorticity_  
	 

	  // Transition Term
      + transPsi  
    );

    tppsiEqn().relax();  
    solve(tppsiEqn);

	

	//*************************************//
    // Calculate eddy viscosity
	// Note: nut is 
    //*************************************//
    T = Ts();

	nut_ = (cN1_ + (1.0 - cN1_)*(psiActual & psiActual)/(cMu_*phiActual*k_))*cMu_*phiActual/epsHat_; 
	nut_ = min(nut_,nutRatMax_*nu()); 
	nut_.correctBoundaryConditions();
    bound(nut_,nut0); 
  

	
    //*************************************//   
    // Output some max values
    //*************************************//
	
	if(debugWrite_ == "true")
	{    
  
	volScalarField meanUz("meanUz",U_.component(2));
	volScalarField uTauSquared((nu() + nut_)*vorticity_.component(2));
	volScalarField cEp1calc(cEp1_ + cEp4_*(2.0*alpha_-1.0));
	volVectorField tpphiVort(tpphi_*vorticity_);
	const volScalarField S2 = 2*magSqr(dev(symm(uGrad_)));
	volScalarField Gnut("Gnut", nut_*S2);
	volScalarField sigPsi("sigPsi", sigmaPsi_*(0.3 + 0.7*min(1.5,tpProd_*k_/epsilon_)));
	volScalarField pd("pd",tpProd_*k_/epsilon_);
	
	Info<< "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl; 
	Info<< "Max SigPsi: " << gMax(sigPsi) << " Min SigPsi: " << gMin(sigPsi) << " Max P/E: " << gMax(pd) << endl;
	Info<< "Max cEp1: " << gMax(cEp1calc) << " Min cEp1: " << gMin(cEp1calc) << endl;  
	Info<< "Max Epsilon: " << gMax(epsilon_) << " Min Epsilon: " << gMin(epsilon_) << " Max G: " << gMax(G) << " Max Gnut: " << gMax(Gnut) <<endl;
	Info<< "Max nut: " << gMax(nut_) << " Max K: " << gMax(k_) << " Max Phi: " << gMax(phiActual) <<endl;
    Info<< "Max Psi: " << gMax(psiActual) << " Min Psi: " << gMin(psiActual)  <<endl;
    Info<< "Max uTauSquared: " << gMax(uTauSquared) <<endl; 
	Info<< "Max vorticity: " << gMax(vorticity_) << " Min vorticity: " << gMin(vorticity_) <<endl;
	Info<< "Max Uz: " << gMax(meanUz) << " Min Uz: " << gMin(meanUz) << endl;
	Info<< "Max gamma: " << gMax(gamma_) << " Min gamma: " << gMin(gamma_) << endl;
	Info<< "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl; 
	}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
