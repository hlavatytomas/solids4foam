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

#include "nonLinGeomTotalLagTotalDispSolidMy.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(nonLinGeomTotalLagTotalDispSolidMy, 0);
addToRunTimeSelectionTable
(
    solidModel, nonLinGeomTotalLagTotalDispSolidMy, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void nonLinGeomTotalLagTotalDispSolidMy::predict()
{
    Info<< "Linear predictor" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinGeomTotalLagTotalDispSolidMy::nonLinGeomTotalLagTotalDispSolidMy
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        det(F_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    predictor_(solidModelDict().lookupOrDefault<Switch>("predictor", false))
{
    DisRequired();

    // Force all required old-time fields to be created
    fvm::d2dt2(D());

    if (predictor_)
    {
        // Check ddt scheme for D is not steadyState
        const word ddtDScheme
        (
#ifdef OPENFOAM_NOT_EXTEND
            mesh().ddtScheme("ddt(" + D().name() +')')
#else
            mesh().schemesDict().ddtScheme("ddt(" + D().name() +')')
#endif
        );

        if (ddtDScheme == "steadyState")
        {
            FatalErrorIn(type() + "::" + type())
                << "If predictor is turned on, then the ddt(" << D().name()
                << ") scheme should not be 'steadyState'!" << abort(FatalError);
        }
    }

    // For consistent restarts, we will update the relative kinematic fields
    D().correctBoundaryConditions();
    if (restart())
    {
        DD() = D() - D().oldTime();
        mechanical().grad(D(), gradD());
        gradDD() = gradD() - gradD().oldTime();
        F_ = I + gradD().T();
        Finv_ = inv(F_);
        J_ = det(F_);

        gradD().storeOldTime();

        // Let the mechanical law know
        mechanical().setRestart();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool nonLinGeomTotalLagTotalDispSolidMy::evolve()
{
    Info<< "Evolving solid solver" << endl;

    if (predictor_)
    {
        predict();
    }

    int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector> solverPerfD;
    SolverPerformance<vector>::debug = 0;
#else
    lduSolverPerformance solverPerfD;
    blockLduMatrix::debug = 0;
#endif

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    // Momentum equation loop
    do
    {

        volScalarField pG = mesh().lookupObject<volScalarField>("pG");

        // volScalarField deltaP = pG - min(pG);
        // volScalarField deltaP = pG - (min(pG) + max(pG)) / 2;
        volScalarField deltaP = pG - dimensionedScalar("dummyP", dimPressure, 1e5);



        // volScalarField alphaF = mesh().lookupObject<volScalarField>("alphaF");
        // volScalarField phiL = mesh().lookupObject<volScalarField>("phiL");


        // volScalarField rho = 0 * dimensionedScalar("rhoL", dimMass / dimVolume, 1000) * alphaF * phiL + dimensionedScalar("rhoS", dimMass / dimVolume, 700) * (1 - alphaF);
        // volScalarField rho = dimensionedScalar("rhoS", dimMass / dimVolume, 700) * (1 - alphaF);
        // volScalarField deltaP = pG;

        // Info << "min deltaP " << min(deltaP).value() << " max deltaP in DEqn " << max(deltaP).value() <<endl;
        // Info << "max(impKf_) " << max(impKf_) << "min(impKf_)" << min(impKf_) << endl;
        // Info << "max(Finv_) " << max(Finv_) << "min(Finv_)" << min(Finv_) << endl;
        // Info << "max(g()) " << max(g()) << "min(g())" << min(g()) << endl;
        // Info << "max(rho()) " << max(rho()) << "min(rho())" << min(rho()) << endl;
        // Info << "max(J_) " << max(J_) << "min(J_)" << min(J_) << "\n" << endl;
        // Info << "max(D()) " << max(D()) << "min(D())" << min(D()) << endl;


        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        // Momentum equation total displacement total Lagrangian form
        fvVectorMatrix DEqn
        (
            // rho*fvm::d2dt2(D())
            rho()*fvm::d2dt2(D())
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(J_*Finv_ & sigma(), "div(sigma)")
        //   + rho*g()
          + rho()*g()
          - fvc::div(J_*Finv_ & deltaP*symmTensor(I))
          + stabilisation().stabilisation(D(), gradD(), impK_)

        // rho()*fvm::d2dt2(D())
        //  ==
        //   - fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
        //   ==
        //   - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
        //   + fvc::div(J_*Finv_ & sigma(), "div(sigma)")
        //   + rho()*g()
        //   - fvc::div(J_*Finv_ & 1/3*deltaP*symmTensor(I))
        //   + stabilisation().stabilisation(D(), gradD(), impK_)
        );

        // mechanical().updateSigmaHyd();

        // Info << "Tu sigmaEq" <<endl;
        // volScalarField sigmaEq = mesh().lookupObject<volScalarField>("sigmaEq");
        // sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma())));
        // sigmaEq.write();

        // Under-relax the linear system
        DEqn.relax();

        // Enforce any cell displacements
        solidModel::setCellDisps(DEqn);

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Info << "D     : res : " << solverPerfD.initialResidual() << endl;

        // Fixed or adaptive field under-relaxation
        relaxField(D(), iCorr);

        // Increment of displacement
        DD() = D() - D().oldTime();

        U() = fvc::ddt(D());

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        F_ = I + gradD().T();
        // F_ = I + gradD();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

    }
    while
    (
       !converged
        (
            iCorr,
#ifdef OPENFOAM_NOT_EXTEND
            mag(solverPerfD.initialResidual()),
            cmptMax(solverPerfD.nIterations()),
#else
            solverPerfD.initialResidual(),
            solverPerfD.nIterations(),
#endif
            D()
        ) && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), gradD(), pointD());

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector>::debug = 1;
#else
    blockLduMatrix::debug = 1;
#endif

    return true;
}


tmp<vectorField> nonLinGeomTotalLagTotalDispSolidMy::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    vectorField nCurrent(Finv.T() & n);
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
