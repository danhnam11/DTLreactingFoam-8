/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FTMMultiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::PtrList<ThermoType>
Foam::FTMMultiComponentMixture<ThermoType>::readSpeciesData
(
    const dictionary& thermoDict
) const
{
    PtrList<ThermoType> specieThermos(species_.size());

    forAll(species_, i)
    {
        specieThermos.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
	//Info << "muCoeffs for species " << i << " = " << specieThermos[i].coeffs().muCoeffs() << endl;
	//Info << "kappaCoeffs for species " << i << " = " << specieThermos[i].coeffs().kappaCoeffs() << endl;
    }

    return specieThermos;
}


template<class ThermoType>
typename Foam::FTMMultiComponentMixture<ThermoType>::speciesCompositionTable
Foam::FTMMultiComponentMixture<ThermoType>::readSpeciesComposition
(
    const dictionary& thermoDict,
    const speciesTable& species
) const
{
    speciesCompositionTable speciesComposition_;

    Info << "[DEBUG] Number of species passed to readSpeciesComposition: " << species.size() << endl;

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(species, si)
    {
/*
        if (thermoDict.subDict(species[si]).isDict("elements"))
        {
            dictionary currentElements
            (
                thermoDict.subDict(species[si]).subDict("elements")
            );

            wordList currentElementsName(currentElements.toc());
            List<specieElement> currentComposition(currentElementsName.size());

            forAll(currentElementsName, eni)
            {
                currentComposition[eni].name() = currentElementsName[eni];

                currentComposition[eni].nAtoms() =
                    currentElements.lookupOrDefault
                    (
                        currentElementsName[eni],
                        0
                    );
            }

            // Add current specie composition to the hash table
            speciesCompositionTable::iterator specieCompositionIter
            (
                speciesComposition_.find(species[si])
            );

            if (specieCompositionIter != speciesComposition_.end())
            {
                speciesComposition_.erase(specieCompositionIter);
            }

            speciesComposition_.insert(species[si], currentComposition);
        }
*/


    const word& specieName = species[si];

    if (!thermoDict.found(specieName))
    {
        Info << "[MISSING] No subDict found for species: " << specieName << endl;
        continue;
    }

    const dictionary& specieDict = thermoDict.subDict(specieName);

    if (!specieDict.isDict("elements"))
    {
        Info << "[WARNING] No 'elements' dict found for species: " << specieName << endl;
        continue;
    }

    Info << "[OK] Parsing elements for species: " << specieName << endl;

    dictionary currentElements(specieDict.subDict("elements"));
    wordList currentElementsName(currentElements.toc());
    List<specieElement> currentComposition(currentElementsName.size());

    forAll(currentElementsName, eni)
    {
        currentComposition[eni].name() = currentElementsName[eni];

        currentComposition[eni].nAtoms() =
            currentElements.lookupOrDefault(currentElementsName[eni], 0);
    }

    speciesCompositionTable::iterator specieCompositionIter
    (
        speciesComposition_.find(specieName)
    );

    if (specieCompositionIter != speciesComposition_.end())
    {
        speciesComposition_.erase(specieCompositionIter);
    }

    speciesComposition_.insert(specieName, currentComposition);

    }

Info << "[DEBUG] Final speciesComposition_ contents:\n";
forAllConstIter(speciesCompositionTable, speciesComposition_, iter)
{
    const word& specieName = iter.key();
    const List<specieElement>& elements = iter();

    Info << "  " << specieName << ": ";

    forAll(elements, ei)
    {
        Info << elements[ei].name() << ":" << elements[ei].nAtoms();
        if (ei != elements.size() - 1) Info << ", ";
    }

    Info << endl;
}

    return speciesComposition_;
}


template<class ThermoType>
void Foam::FTMMultiComponentMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::FTMMultiComponentMixture<ThermoType>::FTMMultiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    specieThermos_(readSpeciesData(thermoDict)),
    speciesComposition_(readSpeciesComposition(thermoDict, species())),
    mixture_("mixture", specieThermos_[0]),
    mixtureVol_("volMixture", specieThermos_[0]),
    numberOfSpecies_(species_.size()),  //
    // for diffusivity
    ListW_(species_.size()),
    muCoeffsMk_(species_.size()),
    kappaCoeffsMk_(species_.size()),
    DijCoeffsMk_(species_.size()),
    mesh_(mesh)
    //
{
    correctMassFractions();

    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]           = specieThermos_[i].W();
        //muCoeffsMk_[i]      = specieThermos_[i].muCoeffs();
        //kappaCoeffsMk_[i]   = specieThermos_[i].kappaCoeffs();
        //DijCoeffsMk_[i]     = specieThermos_[i].DijCoeffs();
        muCoeffsMk_[i]      = specieThermos_[i].coeffs().muCoeffs();
        kappaCoeffsMk_[i]   = specieThermos_[i].coeffs().kappaCoeffs();
        DijCoeffsMk_[i]     = specieThermos_[i].coeffs().DijCoeffs();
    }

    // End of pre-calculation for kinetic model
    //

    mixture_.updateTRANSFitCoeff
    (
        muCoeffsMk_, kappaCoeffsMk_, DijCoeffsMk_
    );

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*specieThermos_[n];
    }
 
   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(Y_.size()); 
    List<scalar> Y(Y_.size()); 
    scalar sumXb = 0.0;  
    forAll(X, i)
    {
        sumXb = sumXb + Y_[i][celli]/ListW_[i]; 
    }  
    if (sumXb == 0){ sumXb = 1e-30;} 

    forAll(X, i)
    {
        X[i] = (Y_[i][celli]/ListW_[i])/sumXb;
        Y[i] = Y_[i][celli];
        if(X[i] <= 0) { X[i] = 0; }
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }
    
    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }
    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS
    (
        Y, X, ListW_
    );

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Y_[0].boundaryField()[patchi][facei]*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*specieThermos_[n];
    }

   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(Y_.size());
    List<scalar> Y(Y_.size());
    scalar sumXb = 0.0;
    forAll(X, i)
    {
        sumXb = sumXb + Y_[i].boundaryField()[patchi][facei]/ListW_[i];
    }
    if (sumXb == 0){ sumXb = 1e-30;}

    forAll(X, i)
    {
        X[i] = (Y_[i].boundaryField()[patchi][facei]/ListW_[i])/sumXb;
        Y[i] = Y_[i].boundaryField()[patchi][facei];
        if(X[i] <= 0) { X[i] = 0; } 
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }

    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }
   
    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS
    (
        Y, X, ListW_
    );

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(specieThermos_, i)
    {
        rhoInv += Y_[i][celli]/specieThermos_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/specieThermos_[0].rho(p, T)/rhoInv*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/specieThermos_[n].rho(p, T)/rhoInv*specieThermos_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(specieThermos_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/specieThermos_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/specieThermos_[0].rho(p, T)/rhoInv
      * specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/specieThermos_[n].rho(p,T)
          / rhoInv*specieThermos_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
void Foam::FTMMultiComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        specieThermos_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}


// ************************************************************************* //
