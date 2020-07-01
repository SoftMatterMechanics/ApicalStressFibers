/*

Copyright (c) 2005-2020, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "PolarityModifier.hpp"
#include "RandomNumberGenerator.hpp"


template<unsigned DIM>
PolarityModifier<DIM>::PolarityModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mD(0.1),
      mDt(0.1),
      mPolarityMagnitude(0.1),
      mAngleForInitialization(M_PI)
{
}

template<unsigned DIM>
PolarityModifier<DIM>::~PolarityModifier()
{
}

template<unsigned DIM>
void PolarityModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdatePolarityOfCells(rCellPopulation);
}

template<unsigned DIM>
void PolarityModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    InitializePolarityOfCells(rCellPopulation);
}

template<unsigned DIM>
void PolarityModifier<DIM>::InitializePolarityOfCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    srand((unsigned)time(NULL));// if srand() used in main time loop, we may not need it here. Try later!
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
         cell_iter != rCellPopulation.rGetCells().end();
         ++cell_iter)
    {
        double angle = mAngleForInitialization; //angle: 0-PI
        double polarity_angle = (rand()%int(round((M_PI+angle)*1e5)))/1e5-0.5*angle;// Err: 0-PI-->>(0-x)-(PI+x)
        polarity_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf();
        SetPolarityOfCell(*cell_iter, polarity_angle, mPolarityMagnitude);
    }
}

template<unsigned DIM>
void PolarityModifier<DIM>::UpdatePolarityOfCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
         cell_iter != rCellPopulation.rGetCells().end();
         ++cell_iter)
    {
        double x = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
        double polarity_angle = (*cell_iter)->GetCellData()->GetItem("PolarityAngle")+sqrt(2*mD*mDt)*x;
        if (polarity_angle>=2*M_PI)
            polarity_angle -= 2*M_PI;
        else if (polarity_angle < 0.0)
            polarity_angle += 2*M_PI;
        SetPolarityOfCell(*cell_iter, polarity_angle, mPolarityMagnitude);
    }    
}

template<unsigned DIM>
void PolarityModifier<DIM>::SetPolarityOfCell(CellPtr pCell, double polarityAngle, double polarityMagnitude)
{
    pCell->GetCellData()->SetItem("PolarityAngle", polarityAngle);
    pCell->GetCellData()->SetItem("PolarityMagnitude", polarityMagnitude);
}


template<unsigned DIM>
void PolarityModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // *rParamsFile << "\t\t\t<ReferenceTargetArea>" << mReferenceTargetArea << "</ReferenceTargetArea>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityModifier<1>;
template class PolarityModifier<2>;
template class PolarityModifier<3>;
