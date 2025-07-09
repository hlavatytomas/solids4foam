/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "viscoPlasticMy.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscoPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, viscoPlastic, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar viscoPlastic::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label viscoPlastic::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar viscoPlastic::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar viscoPlastic::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::viscoPlastic::makeJ()
{
    if (JPtr_.valid())
    {
        FatalErrorIn("void Foam::viscoPlastic::makeJ()")
            << "pointer already set" << abort(FatalError);
    }

    JPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "lawJ",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    // Store the old-time
    JPtr_().oldTime();
}


Foam::volScalarField& Foam::viscoPlastic::J()
{
    if (JPtr_.empty())
    {
        makeJ();
    }

    return JPtr_();
}


void Foam::viscoPlastic::makeJf()
{
    if (JfPtr_.valid())
    {
        FatalErrorIn("void Foam::viscoPlastic::makeJf()")
            << "pointer already set" << abort(FatalError);
    }

    JfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "lawJf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    // Store the old-time
    JfPtr_().oldTime();
}


Foam::surfaceScalarField& Foam::viscoPlastic::Jf()
{
    if (JfPtr_.empty())
    {
        makeJf();
    }

    return JfPtr_();
}


Foam::scalar Foam::viscoPlastic::curYieldStress
(
    const scalar curEpsilonPEq,    // Current equivalent plastic strain
    const scalar J                 // Current Jacobian
) const
{
    // We assume that the stress-strain curve was specifed as Cauchy stress vs
    // true strain, but we want the Kirchhoff (tau) yield stress,
    // so we multiply Cauchy stress by J as tauSigmaY = J*sigmaCauchySigmaY

    return J*stressPlasticStrainSeries_(max(curEpsilonPEq, SMALL));
}


Foam::scalar Foam::viscoPlastic::yieldFunction
(
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar,            // Scaled shear modulus
    const scalar J                 // Current Jacobian
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current Kirchhoff yield stress which is typically a
    // function of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
        magSTrial - 2*muBar*DLambda
      - sqrtTwoOverThree_
           *curYieldStress
            (
                epsilonPEqOld + sqrtTwoOverThree_*DLambda,
                J
            );
}


void Foam::viscoPlastic::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar J,                // Current Jacobian
    const scalar maxMagDEpsilon    // Max strain increment magnitude
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar, J);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluations are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                epsilonPEqOld, magSTrial, DLambda + finiteDiff_, muBar, J
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar,  J);

        if (i == MaxNewtonIter_)
        {
            WarningIn("viscoPlastic::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Note: we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda,  J
        )/J;
}


Foam::tmp<Foam::volScalarField> Foam::viscoPlastic::Ibar
(
    const volSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<volScalarField> tIbar
    (
        new volScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    volScalarField& Ibar = tIbar.ref();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.primitiveFieldRef();
    const symmTensorField devBEbarI = devBEbar.primitiveField();
#else
    volScalarField& Ibar = tIbar();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();
#endif

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
#ifdef OPENFOAM_NOT_EXTEND
            scalarField& IbarP = Ibar.boundaryFieldRef()[patchI];
#else
            scalarField& IbarP = Ibar.boundaryField()[patchI];
#endif
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}


Foam::tmp<Foam::surfaceScalarField> Foam::viscoPlastic::Ibar
(
    const surfaceSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<surfaceScalarField> tIbar
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

#ifdef OPENFOAM_NOT_EXTEND
    surfaceScalarField& Ibar = tIbar.ref();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.primitiveFieldRef();
    const symmTensorField devBEbarI = devBEbar.primitiveField();
#else
    surfaceScalarField& Ibar = tIbar();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();
#endif

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
#ifdef OPENFOAM_NOT_EXTEND
            scalarField& IbarP = Ibar.boundaryFieldRef()[patchI];
#else
            scalarField& IbarP = Ibar.boundaryField()[patchI];
#endif
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

#ifndef OPENFOAM_NOT_EXTEND
    Ibar.correctBoundaryConditions();
#endif

    return tIbar;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::viscoPlastic::viscoPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    JPtr_(),
    JfPtr_(),
    stressPlasticStrainSeries_(dict),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    ),
    sigmaYf_
    (
        IOobject
        (
            "sigmaYf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
           "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    ),
    DSigmaY_
    (
        IOobject
        (
            "DSigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    DSigmaYf_
    (
        IOobject
        (
            "DSigmaYf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonPf_
    (
        IOobject
        (
            "epsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonPf_
    (
        IOobject
        (
            "DEpsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    bEbarTrial_
    (
        IOobject
        (
            "bEbarTrial",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarTrialf_
    (
        IOobject
        (
            "bEbarTrialf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarf_
    (
        IOobject
        (
            "bEbarf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    DEpsilonPEq_
    (
        IOobject
        (
            "DEpsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DEpsilonPEqf_
    (
        IOobject
        (
            "DEpsilonPEqf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambdaf_
    (
        IOobject
        (
            "DLambdaf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEqf_
    (
        IOobject
        (
            "epsilonPEqf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    plasticN_
    (
        IOobject
        (
            "plasticN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    plasticNf_
    (
        IOobject
        (
            "plasticNf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    updateBEbarConsistent_
    (
        mechanicalLaw::dict().lookupOrAddDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    ),
    Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    Info<< "    updateBEbarConsistent: " << updateBEbarConsistent_ << endl;

    // Force storage of old time for adjustable time-step calculations
    plasticN_.storeOldTime();

    // Force the creation of Fs so they are read on restart
    F();
    Ff();

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E/(2.0*(1.0 + nu));

        // Set the bulk modulus
        if (planeStress())
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));
    }
    else
    {
        FatalErrorIn
        (
            "viscoPlastic::viscoPlastic::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Check if plasticity is a nonlinear function of plastic strain
    if (nonLinearPlasticity_)
    {
        Info<< "    Plasticity is nonlinear" << endl;
    }
    else
    {
        if (stressPlasticStrainSeries_.size() == 1)
        {
            Info<< "    Perfect Plasticity" << endl;
        }
        else
        {
            Info<< "    Plasticity is linear" << endl;

            // Define linear plastic modulus
            Hp_ =
                (
                    stressPlasticStrainSeries_[1].second()
                  - stressPlasticStrainSeries_[0].second()
                )
               /(
                    stressPlasticStrainSeries_[1].first()
                  - stressPlasticStrainSeries_[0].first()
                );
        }
    }

    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::viscoPlastic::~viscoPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscoPlastic::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(mu_*dev(bEbarTrial_));

    const volScalarField Ibar(tr(bEbarTrial_)/3.0);
    const volScalarField muBar(Ibar*mu_);

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*muBar*DLambda_/magSTrial));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            scaleFactor*(4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::viscoPlastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    J() = det(F());

    // Calculate the relative Jacobian
    // const volScalarField relJ(J()/J().oldTime());

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    // const volTensorField relFbar(pow(relJ, -1.0/3.0)*relF());

    // Update bE trial
    // bEbarTrial_ = transform(relFbar, bEbar_.oldTime());

    // Calculate trial deviatoric stress
    // const volSymmTensorField sTrial(mu_*dev(bEbarTrial_));

    // const volScalarField Ibar(tr(bEbarTrial_)/3.0);
    // const volScalarField muBar(Ibar*mu_);

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonP_.storePrevIter();

    // prepare DEpsilon
    const Time& time = mesh().time();
    scalar dTimeSc = time.deltaTValue();
    dimensionedScalar dTime("dTime", dimTime, dTimeSc);
    volScalarField T = mesh().lookupObject<volScalarField>("T") / dimensionedScalar("dummyT", dimTemperature, 1) - 276;
    volScalarField alphaS = mesh().lookupObject<volScalarField>("alphaS");
    volScalarField alphaL = mesh().lookupObject<volScalarField>("alphaL");
    // dimensionedScalar T ("temp", dimless, 27);
    volScalarField alphaG = 1 - alphaL - alphaS;

    // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 45) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1); //* (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1) * (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 45) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1); //* (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 35) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1); //* (- Foam::atan(4e4 * alphaG - 4e3) / 1e-3 + 1571.75);
    // volScalarField tau =  dimensionedScalar("dummyTime", dimTime, 1) * (- Foam::atan(1e3 * alphaG - 1e2) / 0.05 + 32.36);
    // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1) ;
    // forAll(tau, cellI)
    // {
    //     if (alphaG[cellI] < 0.1)
    //     {
    //         tau[cellI] = tau[cellI] * (- Foam::atan(5e6 * alphaG[cellI] - 5e5) / 1e-3 + 1571.8);
    //     }
    // }


    // dimensionedScalar tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    // dimensionedScalar tau = dimensionedScalar("dummyTime", dimTime, 2);
    dimensionedScalar E = 9 * mu_ * K_ / (3 * K_ + mu_);
    dimensionedScalar nu = 0.5 * (3 * K_ - 2 * mu_) / (3 * K_ + mu_);
    dimensionedScalar predMatice = 1 / ((1 + nu) * (1 - 2 * nu));

    // volScalarField Epokus = E * (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) ; 
    // volScalarField Epokus = E * (- Foam::atan(1e3 * alphaG - 1e2) / 0.5 + 4.14); 
    volTensorField F = mesh().lookupObject<volTensorField>("F");
    volScalarField J = det(F);
    volTensorField invF = inv(F);
    volTensorField S = J * invF & sigma & invF.T();
    // volSymmTensorField D0 = J() * sigma;
    // volSymmTensorField D0 = sigma;
    volSymmTensorField D0 = 0 * sigma;
    D0.replace(symmTensor::XX, S.component(tensor::XX) - nu * S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    D0.replace(symmTensor::XY, (2 * nu + 2) * S.component(tensor::XY));
    D0.replace(symmTensor::XZ, (2 * nu + 2) * S.component(tensor::XZ));
    D0.replace(symmTensor::YY, - nu * S.component(tensor::XX) + S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    D0.replace(symmTensor::YZ, (2 * nu + 2) * S.component(tensor::YZ));
    D0.replace(symmTensor::ZZ, - nu * S.component(tensor::XX) - nu * S.component(tensor::YY) + S.component(tensor::ZZ));
    // D0.replace(symmTensor::XX, sigma.component(symmTensor::XX) - nu * sigma.component(symmTensor::YY) - nu * sigma.component(symmTensor::ZZ));
    // D0.replace(symmTensor::XY, (2 * nu + 2) * sigma.component(symmTensor::XY));
    // D0.replace(symmTensor::XZ, (2 * nu + 2) * sigma.component(symmTensor::XZ));
    // D0.replace(symmTensor::YY, - nu * sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY) - nu * sigma.component(symmTensor::ZZ));
    // D0.replace(symmTensor::YZ, (2 * nu + 2) * sigma.component(symmTensor::YZ));
    // D0.replace(symmTensor::ZZ, - nu * sigma.component(symmTensor::XX) - nu * sigma.component(symmTensor::YY) + sigma.component(symmTensor::ZZ));

    // Info << "E = " << E << endl;
    // Info << "tau = " << tau << endl;
    // Info << "nu = " << nu << endl;
    // Info << "dTime = " << dTime << endl; 

    // Update DEpsilonP and DEpsilonPEq
    // DEpsilonP_ = 1 / E / tau * D0 * dTime;

    volTensorField pomoc = 1 / J / E / tau * dTime * (F & D0 & F.T());
    DEpsilonP_.replace(symmTensor::XX, pomoc.component(tensor::XX));
    DEpsilonP_.replace(symmTensor::XY, pomoc.component(tensor::XY));
    DEpsilonP_.replace(symmTensor::XZ, pomoc.component(tensor::XZ));
    DEpsilonP_.replace(symmTensor::YY, pomoc.component(tensor::YY));
    DEpsilonP_.replace(symmTensor::YZ, pomoc.component(tensor::YZ));
    DEpsilonP_.replace(symmTensor::ZZ, pomoc.component(tensor::ZZ));


    // DEpsilonP_ = 1 / Epokus / tau * D0 * dTime;
    // DEpsilonP_ = 1 / E / tau * sigma * dTime;
    DEpsilonP_.correctBoundaryConditions();

    const volTensorField& gradDD = mesh().lookupObject<volTensorField>("grad(DD)");

    // volVectorField D = mesh().lookupObject<volVectorField>("D");

    // volVectorField DD = D - D.oldTime();

    // volTensorField gradDD = fvc::grad(DD);

    
    volSymmTensorField dEpsilon = symm(gradDD);
    // volSymmTensorField dEpsilon = F & symm(gradDD) & F.T();
    volTensorField pomoc2 = 1/ J * F & symm(gradDD) & F.T();
    dEpsilon.replace(symmTensor::XX, pomoc2.component(tensor::XX));
    dEpsilon.replace(symmTensor::XY, pomoc2.component(tensor::XY));
    dEpsilon.replace(symmTensor::XZ, pomoc2.component(tensor::XZ));
    dEpsilon.replace(symmTensor::YY, pomoc2.component(tensor::YY));
    dEpsilon.replace(symmTensor::YZ, pomoc2.component(tensor::YZ));
    dEpsilon.replace(symmTensor::ZZ, pomoc2.component(tensor::ZZ));

    // volSymmTensorField dSigma = Epokus * dEpsilon;
    volSymmTensorField dSigma = E * dEpsilon;

    dSigma.replace(symmTensor::XX, E * predMatice * ((1 - nu) * dEpsilon.component(symmTensor::XX) + (nu) * dEpsilon.component(symmTensor::YY) + (nu) * dEpsilon.component(symmTensor::ZZ)));
    dSigma.replace(symmTensor::XY, E * predMatice * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::XY));
    dSigma.replace(symmTensor::XZ, E * predMatice * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::XZ));
    dSigma.replace(symmTensor::YY, E * predMatice * ((nu) * dEpsilon.component(symmTensor::XX) + (1 - nu) * dEpsilon.component(symmTensor::YY) + (nu) * dEpsilon.component(symmTensor::ZZ)));
    dSigma.replace(symmTensor::YZ, E * predMatice * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::YZ));
    dSigma.replace(symmTensor::ZZ, E * predMatice * ((nu) * dEpsilon.component(symmTensor::XX) + (nu) * dEpsilon.component(symmTensor::YY) + (1 - nu) * dEpsilon.component(symmTensor::ZZ)));


    // volSymmTensorField dSigmaP = Epokus * DEpsilonP_;
    volSymmTensorField dSigmaP = E * DEpsilonP_;

    dSigmaP.replace(symmTensor::XX, E * predMatice * ((1 - nu) * DEpsilonP_.component(symmTensor::XX) + (nu) * DEpsilonP_.component(symmTensor::YY) + (nu) * DEpsilonP_.component(symmTensor::ZZ)));
    dSigmaP.replace(symmTensor::XY, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonP_.component(symmTensor::XY));
    dSigmaP.replace(symmTensor::XZ, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonP_.component(symmTensor::XZ));
    dSigmaP.replace(symmTensor::YY, E * predMatice * ((nu) * DEpsilonP_.component(symmTensor::XX) + (1 - nu) * DEpsilonP_.component(symmTensor::YY) + (nu) * DEpsilonP_.component(symmTensor::ZZ)));
    dSigmaP.replace(symmTensor::YZ, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonP_.component(symmTensor::YZ));
    dSigmaP.replace(symmTensor::ZZ, E * predMatice * ((nu) * DEpsilonP_.component(symmTensor::XX) + (nu) * DEpsilonP_.component(symmTensor::YY) + (1 - nu) * DEpsilonP_.component(symmTensor::ZZ)));

    sigma = sigma.oldTime() + (dSigma - dSigmaP);
    

    // prepare DSigmaEq
    // DEpsilonPEq_ = sqrt(2.0/3.0 * (DEpsilonP_ && DEpsilonP_));
    // forAll(DSigmaY_, cellI)
    // {
    //     scalar curSigmaY =
    //         curYieldStress
    //         (
    //             epsilonPEq_[cellI] + DEpsilonPEq_[cellI],  J()[cellI]
    //         )/J()[cellI];
    //     DSigmaY_[cellI] = curSigmaY - sigmaY_[cellI];
    // }
    // DSigmaY_.correctBoundaryConditions();

    // // Calculate deviatoric stress
    // const volSymmTensorField s(sTrial - 2*mu_*DEpsilonP_);

    // // Update bEbar
    // if (updateBEbarConsistent_)
    // {
    //     const volSymmTensorField devBEbar(s/mu_);
    //     bEbar_ = devBEbar + this->Ibar(devBEbar)*I;
    // }
    // else
    // {
    //     bEbar_ = (s/mu_) + Ibar*I;
    // }

    // // Update hydrostatic stress (negative of pressure)
    // updateSigmaHyd
    // (
    //     0.5*K_*(pow(J(), 2.0) - 1.0),
    //     (4.0/3.0)*mu_ + K_
    // );

    // // Update the Cauchy stress
    // sigma = (1.0/J())*(sigmaHyd()*I + s);
}


void Foam::viscoPlastic::correct(surfaceSymmTensorField& sigma)
{
    Info << "Pozor surfaceField" <<endl;
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    Jf() = det(Ff());


    DEpsilonPf_.storePrevIter();

    // prepare DEpsilon
    const Time& time = mesh().time();
    scalar dTimeSc = time.deltaTValue();
    dimensionedScalar dTime("dTime", dimTime, dTimeSc);
    volScalarField T = mesh().lookupObject<volScalarField>("T") / dimensionedScalar("dummyT", dimTemperature, 1) - 276;
    // dimensionedScalar T ("temp", dimless, 27);
    // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    volScalarField tau = (50 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    // dimensionedScalar tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    // dimensionedScalar tau = dimensionedScalar("dummyTime", dimTime, 2);
    dimensionedScalar E = 9 * mu_ * K_ / (3 * K_ + mu_);
    dimensionedScalar nu = 0.5 * (3 * K_ - 2 * mu_) / (3 * K_ + mu_);
    dimensionedScalar predMatice = 1 / ((1 + nu) * (1 - 2 * nu));
    // volTensorField F = mesh().lookupObject<volTensorField>("F");
    // volTensorField invF = inv(F);
    // volTensorField S = J() * invF & sigma & invF.T();
    // volSymmTensorField D0 = J() * sigma;
    // volSymmTensorField D0 = sigma;
    surfaceSymmTensorField D0 = 0 * sigma;
    // D0.replace(symmTensor::XX, S.component(tensor::XX) - nu * S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    // D0.replace(symmTensor::XY, (2 * nu + 2) * S.component(tensor::XY));
    // D0.replace(symmTensor::XZ, (2 * nu + 2) * S.component(tensor::XZ));
    // D0.replace(symmTensor::YY, - nu * S.component(tensor::XX) + S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    // D0.replace(symmTensor::YZ, (2 * nu + 2) * S.component(tensor::YZ));
    // D0.replace(symmTensor::ZZ, - nu * S.component(tensor::XX) - nu * S.component(tensor::YY) + S.component(tensor::ZZ));
    D0.replace(symmTensor::XX, sigma.component(symmTensor::XX) - nu * sigma.component(symmTensor::YY) - nu * sigma.component(symmTensor::ZZ));
    D0.replace(symmTensor::XY, (2 * nu + 2) * sigma.component(symmTensor::XY));
    D0.replace(symmTensor::XZ, (2 * nu + 2) * sigma.component(symmTensor::XZ));
    D0.replace(symmTensor::YY, - nu * sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY) - nu * sigma.component(symmTensor::ZZ));
    D0.replace(symmTensor::YZ, (2 * nu + 2) * sigma.component(symmTensor::YZ));
    D0.replace(symmTensor::ZZ, - nu * sigma.component(symmTensor::XX) - nu * sigma.component(symmTensor::YY) + sigma.component(symmTensor::ZZ));

    // Info << "E = " << E << endl;
    // Info << "tau = " << tau << endl;
    // Info << "nu = " << nu << endl;
    // Info << "dTime = " << dTime << endl; 

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPf_ = 1 / E / fvc::interpolate(tau) * D0 * dTime;
    // DEpsilonP_ = 1 / E / tau * sigma * dTime;
    DEpsilonPf_.correctBoundaryConditions();

    const volTensorField& gradDD = mesh().lookupObject<volTensorField>("grad(DD)");

    
    surfaceSymmTensorField dEpsilonf = symm(fvc::interpolate(gradDD));

    surfaceSymmTensorField dSigmaf = E * dEpsilonf;

    dSigmaf.replace(symmTensor::XX, E * predMatice * ((1 - nu) * dEpsilonf.component(symmTensor::XX) + (nu) * dEpsilonf.component(symmTensor::YY) + (nu) * dEpsilonf.component(symmTensor::ZZ)));
    dSigmaf.replace(symmTensor::XY, E * predMatice * (1 - 2 * nu) / 2 * dEpsilonf.component(symmTensor::XY));
    dSigmaf.replace(symmTensor::XZ, E * predMatice * (1 - 2 * nu) / 2 * dEpsilonf.component(symmTensor::XZ));
    dSigmaf.replace(symmTensor::YY, E * predMatice * ((nu) * dEpsilonf.component(symmTensor::XX) + (1 - nu) * dEpsilonf.component(symmTensor::YY) + (nu) * dEpsilonf.component(symmTensor::ZZ)));
    dSigmaf.replace(symmTensor::YZ, E * predMatice * (1 - 2 * nu) / 2 * dEpsilonf.component(symmTensor::YZ));
    dSigmaf.replace(symmTensor::ZZ, E * predMatice * ((nu) * dEpsilonf.component(symmTensor::XX) + (nu) * dEpsilonf.component(symmTensor::YY) + (1 - nu) * dEpsilonf.component(symmTensor::ZZ)));


    surfaceSymmTensorField dSigmaPf = E * DEpsilonPf_;

    dSigmaPf.replace(symmTensor::XX, E * predMatice * ((1 - nu) * DEpsilonPf_.component(symmTensor::XX) + (nu) * DEpsilonPf_.component(symmTensor::YY) + (nu) * DEpsilonPf_.component(symmTensor::ZZ)));
    dSigmaPf.replace(symmTensor::XY, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonPf_.component(symmTensor::XY));
    dSigmaPf.replace(symmTensor::XZ, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonPf_.component(symmTensor::XZ));
    dSigmaPf.replace(symmTensor::YY, E * predMatice * ((nu) * DEpsilonPf_.component(symmTensor::XX) + (1 - nu) * DEpsilonPf_.component(symmTensor::YY) + (nu) * DEpsilonPf_.component(symmTensor::ZZ)));
    dSigmaPf.replace(symmTensor::YZ, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonPf_.component(symmTensor::YZ));
    dSigmaPf.replace(symmTensor::ZZ, E * predMatice * ((nu) * DEpsilonPf_.component(symmTensor::XX) + (nu) * DEpsilonPf_.component(symmTensor::YY) + (1 - nu) * DEpsilonPf_.component(symmTensor::ZZ)));

    sigma = sigma.oldTime() + (dSigmaf - dSigmaPf);

    // // Calculate the relative Jacobian
    // const surfaceScalarField relJ(Jf()/Jf().oldTime());

    // // Calculate the relative deformation gradient with the volumetric term
    // // removed
    // // const surfaceTensorField relFbar(pow(relJ, -1.0/3.0)*relFf());

    // // // Update bE trial
    // // bEbarTrialf_ = transform(relFbar, bEbarf_.oldTime());

    // // // Calculate trial deviatoric stress
    // // const surfaceSymmTensorField sTrial(mu_*dev(bEbarTrialf_));

    // // const surfaceScalarField Ibar(tr(bEbarTrialf_)/3.0);
    // // const surfaceScalarField muBar(Ibar*mu_);

    // // Check for plastic loading
    // // and calculate increment of plastic equivalent strain
    // // i.e. the plastic multiplier

    // // Store previous iteration for under-relaxation and calculation of plastic
    // // residual in the solver

    // DEpsilonPf_.storePrevIter();

    // // prepare DEpsilon
    // const Time& time = mesh().time();
    // scalar dTimeSc = time.deltaTValue();
    // dimensionedScalar dTime("dTime", dimTime, dTimeSc);
    // volScalarField T = mesh().lookupObject<volScalarField>("T") / dimensionedScalar("dummyT", dimTemperature, 1) - 276;
    // // dimensionedScalar T ("temp", dimless, 27);
    // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    // // dimensionedScalar tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    // // dimensionedScalar tau = dimensionedScalar("dummyTime", dimTime, 2);
    // dimensionedScalar E = 9 * mu_ * K_ / (3 * K_ + mu_);
    // dimensionedScalar nu = 0.5 * (3 * K_ - 2 * mu_) / (3 * K_ + mu_);
    // dimensionedScalar predMatice = 1 / ((1 + nu) * (1 - 2 * nu));
    // // volTensorField F = mesh().lookupObject<volTensorField>("F");
    // // volTensorField invF = inv(F);
    // // volTensorField S = J() * invF & sigma & invF.T();
    // // volSymmTensorField D0 = J() * sigma;
    // // volSymmTensorField D0 = sigma;
    // surfaceSymmTensorField D0 = 0 * sigma;
    // // D0.replace(symmTensor::XX, S.component(tensor::XX) - nu * S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    // // D0.replace(symmTensor::XY, (2 * nu + 2) * S.component(tensor::XY));
    // // D0.replace(symmTensor::XZ, (2 * nu + 2) * S.component(tensor::XZ));
    // // D0.replace(symmTensor::YY, - nu * S.component(tensor::XX) + S.component(tensor::YY) - nu * S.component(tensor::ZZ));
    // // D0.replace(symmTensor::YZ, (2 * nu + 2) * S.component(tensor::YZ));
    // // D0.replace(symmTensor::ZZ, - nu * S.component(tensor::XX) - nu * S.component(tensor::YY) + S.component(tensor::ZZ));
    // D0.replace(symmTensor::XX, sigma.component(symmTensor::XX) - nu * sigma.component(symmTensor::YY) - nu * sigma.component(symmTensor::ZZ));
    // D0.replace(symmTensor::XY, (2 * nu + 2) * sigma.component(symmTensor::XY));
    // D0.replace(symmTensor::XZ, (2 * nu + 2) * sigma.component(symmTensor::XZ));
    // D0.replace(symmTensor::YY, - nu * sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY) - nu * sigma.component(symmTensor::ZZ));
    // D0.replace(symmTensor::YZ, (2 * nu + 2) * sigma.component(symmTensor::YZ));
    // D0.replace(symmTensor::ZZ, - nu * sigma.component(symmTensor::XX) - nu * sigma.component(symmTensor::YY) + sigma.component(symmTensor::ZZ));

    // // Update DEpsilonP and DEpsilonPEq
    // DEpsilonPf_ = 1 / E / fvc::interpolate(tau) * D0 * dTime;
    // // DEpsilonPf_ = 1 / E / tau * D0 * dTime;
    // // DEpsilonP_ = 1 / E / tau * sigma * dTime;
    // DEpsilonPf_.correctBoundaryConditions();

    // // const volTensorField& gradDD = mesh().lookupObject<volTensorField>("grad(DD)");

    
    // // surfaceSymmTensorField dEpsilon = fvc::interpolate(symm(gradDD));

    // // surfaceSymmTensorField dSigma = E * dEpsilon;

    // // dSigma.replace(symmTensor::XX, E * predMatice * ((1 - nu) * dEpsilon.component(symmTensor::XX) + (nu) * dEpsilon.component(symmTensor::YY) + (nu) * dEpsilon.component(symmTensor::ZZ)));
    // // dSigma.replace(symmTensor::XY, E * predMatice * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::XY));
    // // dSigma.replace(symmTensor::XZ, E * predMatice * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::XZ));
    // // dSigma.replace(symmTensor::YY, E * predMatice * ((nu) * dEpsilon.component(symmTensor::XX) + (1 - nu) * dEpsilon.component(symmTensor::YY) + (nu) * dEpsilon.component(symmTensor::ZZ)));
    // // dSigma.replace(symmTensor::YZ, E * predMatice * (1 - 2 * nu) / 2 * dEpsilon.component(symmTensor::YZ));
    // // dSigma.replace(symmTensor::ZZ, E * predMatice * ((nu) * dEpsilon.component(symmTensor::XX) + (nu) * dEpsilon.component(symmTensor::YY) + (1 - nu) * dEpsilon.component(symmTensor::ZZ)));


    // // surfaceSymmTensorField dSigmaP = E * DEpsilonPf_;

    // // dSigmaP.replace(symmTensor::XX, E * predMatice * ((1 - nu) * DEpsilonPf_.component(symmTensor::XX) + (nu) * DEpsilonPf_.component(symmTensor::YY) + (nu) * DEpsilonPf_.component(symmTensor::ZZ)));
    // // dSigmaP.replace(symmTensor::XY, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonPf_.component(symmTensor::XY));
    // // dSigmaP.replace(symmTensor::XZ, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonPf_.component(symmTensor::XZ));
    // // dSigmaP.replace(symmTensor::YY, E * predMatice * ((nu) * DEpsilonPf_.component(symmTensor::XX) + (1 - nu) * DEpsilonPf_.component(symmTensor::YY) + (nu) * DEpsilonPf_.component(symmTensor::ZZ)));
    // // dSigmaP.replace(symmTensor::YZ, E * predMatice * (1 - 2 * nu) / 2 * DEpsilonPf_.component(symmTensor::YZ));
    // // dSigmaP.replace(symmTensor::ZZ, E * predMatice * ((nu) * DEpsilonPf_.component(symmTensor::XX) + (nu) * DEpsilonPf_.component(symmTensor::YY) + (1 - nu) * DEpsilonPf_.component(symmTensor::ZZ)));

    // // sigma = sigma.oldTime() + (dSigma - dSigmaP);


    // // DEpsilonPf_.storePrevIter();

    // // const Time& time = mesh().time();
    // // scalar dTimeSc = time.deltaTValue();
    // // dimensionedScalar dTime("dTime", dimTime, dTimeSc);
    // // volScalarField T = mesh().lookupObject<volScalarField>("T") / dimensionedScalar("dummyT", dimTemperature, 1) - 276;
    // // volScalarField tau = (9 * (2 / 3.14 * Foam::atan((T - 65) / 2) + 1) + 2) * dimensionedScalar("dummyTime", dimTime, 1);
    // // dimensionedScalar E = 9 * mu_ * K_ / (3 * K_ + mu_);
    // // surfaceSymmTensorField D0f = Jf() * sigma;

    // // Update DEpsilonP and DEpsilonPEq
    // // DEpsilonPf_ = 1 / E / fvc::interpolate(tau) * sigma * dTime;
    // // DEpsilonPf_.correctBoundaryConditions();

    // DEpsilonPEqf_ = sqrt(2.0/3.0 * (DEpsilonPf_ && DEpsilonPf_));
    // forAll(DSigmaYf_, faceI)
    // {
    //     scalar curSigmaY =
    //         curYieldStress
    //         (
    //             epsilonPEqf_[faceI] + DEpsilonPEqf_[faceI],  Jf()[faceI]
    //         )/Jf()[faceI];
    //     DSigmaYf_[faceI] = curSigmaY - sigmaYf_[faceI];
    // }
    // DSigmaYf_.correctBoundaryConditions();

    // // Calculate deviatoric stress
    // const surfaceSymmTensorField s(sTrial - 2*mu_*DEpsilonPf_);

    // // Update bEbar
    // if (updateBEbarConsistent_)
    // {
    //     const surfaceSymmTensorField devBEbar(s/mu_);
    //     bEbarf_ = devBEbar + this->Ibar(devBEbar)*I;
    // }
    // else
    // {
    //     bEbarf_ = (s/mu_) + Ibar*I;
    // }

    // // Update the Cauchy stress
    // // Note: updayeSigmaHyd is not implemented for surface fields
    // sigma = (1.0/Jf())*(0.5*K_*(pow(Jf(), 2) - 1)*I + s);
}


Foam::scalar Foam::viscoPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().time().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).foundObject<surfaceTensorField>("Ff")
    )
    {
        return
#ifdef OPENFOAM_NOT_EXTEND
            gMax
            (
                mag
                (
                    DEpsilonPf_.primitiveField()
                  - DEpsilonPf_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().primitiveField()));
#else
            gMax
            (
                mag
                (
                    DEpsilonPf_.internalField()
                  - DEpsilonPf_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().internalField()));
#endif
    }
    else
    {
        return
#ifdef OPENFOAM_NOT_EXTEND
            gMax
            (
                mag
                (
                    DEpsilonP_.primitiveField()
                  - DEpsilonP_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().primitiveField()));
#else
            gMax
            (
                mag
                (
                    DEpsilonP_.internalField()
                  - DEpsilonP_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().internalField()));
#endif
    }
}


void Foam::viscoPlastic::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;
    sigmaYf_ += DSigmaYf_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
    epsilonPEqf_ += DEpsilonPEqf_;
    epsilonP_ += DEpsilonP_;
    epsilonPf_ += DEpsilonPf_;

    // Count cells actively yielding
    int numCellsYielding = 0;

#ifdef OPENFOAM_NOT_EXTEND
    forAll(activeYield_.primitiveField(), celli)
    {
        if (DEpsilonPEq_.primitiveField()[celli] > SMALL)
        {
            activeYield_.primitiveFieldRef()[celli] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.primitiveFieldRef()[celli] = 0.0;
        }
    }
#else
    forAll(activeYield_.internalField(), celli)
    {
        if (DEpsilonPEq_.internalField()[celli] > SMALL)
        {
            activeYield_.internalField()[celli] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.internalField()[celli] = 0.0;
        }
    }
#endif

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchi)
    {
        if (!activeYield_.boundaryField()[patchi].coupled())
        {
            forAll(activeYield_.boundaryField()[patchi], facei)
            {
                if (DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
                {
#ifdef OPENFOAM_NOT_EXTEND
                    activeYield_.boundaryFieldRef()[patchi][facei] = 1.0;
#else
                    activeYield_.boundaryField()[patchi][facei] = 1.0;
#endif
                }
                else
                {
#ifdef OPENFOAM_NOT_EXTEND
                    activeYield_.boundaryFieldRef()[patchi][facei] = 0.0;
#else
                    activeYield_.boundaryField()[patchi][facei] = 0.0;
#endif
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


Foam::scalar Foam::viscoPlastic::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Update the total deformatio gradient: already done by updateF
    // if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    // {
    //     F() = fvc::average(relFf()) & F().oldTime();
    // }
    // else
    // {
    //     F() = relF() & F().oldTime();
    // }

    // Calculate the total true (Hencky) strain
    const volSymmTensorField epsilon(0.5*log(symm(F().T() & F())));

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon))));

    // Take reference to internal fields
#ifdef OPENFOAM_NOT_EXTEND
    const symmTensorField& DEpsilonPI = DEpsilonP_.primitiveField();
    const symmTensorField& plasticNI = plasticN_.primitiveField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().primitiveField();
    const scalarField& epsilonEqI = epsilonEq.primitiveField();
#else
    const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    const symmTensorField& plasticNI = plasticN_.internalField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();
#endif

    // Calculate error field
    const symmTensorField DEpsilonPErrorI
    (
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
       /(epsilonEqI + SMALL)
    );

    // Max error
    const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    if (maxMagDEpsilonPErr > SMALL)
    {
        Info<< "    " << name() << ": max time integration error = "
            << maxMagDEpsilonPErr
            << endl;

        if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
        {
            WarningIn
            (
                "Foam::scalar Foam::viscoPlastic::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is lover 50 times larger "
                << "than the desired value!\n    Consider starting the "
                << "simulation with a smaller initial time-step" << endl;
        }

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


void Foam::viscoPlastic::setRestart()
{
    F().writeOpt() = IOobject::AUTO_WRITE;
    J().writeOpt() = IOobject::AUTO_WRITE;
    bEbar_.writeOpt() = IOobject::AUTO_WRITE;
    epsilonPEq_.writeOpt() = IOobject::AUTO_WRITE;

    Ff().writeOpt() = IOobject::AUTO_WRITE;
    Jf().writeOpt() = IOobject::AUTO_WRITE;
    bEbarf_.writeOpt() = IOobject::AUTO_WRITE;
    epsilonPEqf_.writeOpt() = IOobject::AUTO_WRITE;
}

// ************************************************************************* //
