/*

Copyright (c) 2005-2018, University of Oxford.
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

#ifndef MYSTRESSFIBERTENSIONFORCE_HPP_
#define MYSTRESSFIBERTENSIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <iostream>

/**
 * In this vertex-based simulations, stressfiber tension force refers to the force generated
 * by stress fibers including passive and active components.
 */
template<unsigned DIM>
class MyStressfiberTensionForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    
    std::string mOutputDirectory;      

    /** Results file for reaction forces. */
    out_stream mpPeelingBisectorOrientationFile; 

    out_stream mpAngleBisectorFile;

    out_stream mpStressFiberTensionFile;

    out_stream mpElementInfoFile;

    out_stream mpCellStressFile;
    
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

protected:

    unsigned mAreaSeed;

    bool   mIfEquilibrateForAWhile;

    double mStartTimeForStretching;

    unsigned mFlag;

    double mSfStiffness;

    double mNucleationPerimeterTension;

    double mHalfWidth;

    double mRestLengthOfNucleation;

    double mAdhesionEnergy;

    double mk;

    double mC0;

    double mRatePower;

    double mCellCellAdhesionEnergyParameter;


public:

    /**
     * Constructor.
     */
    MyStressfiberTensionForce();

    /**
     * Destructor.
     */
    virtual ~MyStressfiberTensionForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the morphogenetic force on each node in the vertex-based cell population
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

    // my methods
    void UpdateStressStateOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);

    void SetAreaSeed(unsigned areaSeed);

    void SetIfEquilibrateForAWhile(bool ifEquilibrateForAWhile);

    void SetStartTimeForStretching(double startTimeForStretching);

    void SetFlagForStressfiberCreation(unsigned flag);

    void SetStressfiberStiffness(double sfStiffness);

    void SetNucleationThresholdOfPerimeterTension(double nucleationPerimeterTension);

    void SetHalfWidth(double halfWidth);

    void SetRestLengthOfNucleation(double restLengthOfNucleation);

    void SetPeelingParameters(double adhesionEnergy, double k, double C0, double ratePower);

    void SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter);

    void SetOutputDirectory(std::string outputDirectory);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyStressfiberTensionForce)

#endif /*MYSTRESSFIBERTENSIONFORCE_HPP_*/
