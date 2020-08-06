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

#ifndef FACEVALUEANDSTRESSSTATEMODIFIER_HPP_
#define FACEVALUEANDSTRESSSTATEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"
#include "MutableVertexMesh.hpp"

/**
 * A modifier class in which the target area property of each cell is updated.
 * It is used to implement growth in vertex-based simulations.
 */
template<unsigned DIM>
class FaceValueAndStressStateModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

protected:

    bool mIfConsiderFeedbackOfFaceValues;
    bool mIfConsiderFeedbackOfFaceValuesOnlyForBoundaryCells;
    bool mIfConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells;
    bool mIfConsiderFeedbackOfCellCellAdhesion;

    bool mEMADontDecreaseWhenEdgeShrink;
    bool mCCADontDecreaseWhenEdgeExpand;
    bool mCCAIncreasingHasAThresholdOfEdgeLength;
    double mCCAIncreasingThresholdOfEdgeLengthPercentage;

    double mEdgeLengthAtRest;
    double mFeedbackStrengthForMyosinActivity;
    double mHillCoefficientForMyosinActivity;
    double mFeedbackStrengthForAdhesion;
    double mHillCoefficientForAdhesion;

    bool mIfCalculateStressState;
    unsigned mCaseNumberOfMembraneSurfaceEnergyForm;
    double mFixedTargetArea;
    double mFixedTargetPerimeter;
    double mNagaiHondaDeformationEnergyParameter;
    double mNagaiHondaMembraneSurfaceEnergyParameter;
    double mNagaiHondaCellCellAdhesionEnergyParameter;
    double mNagaiHondaCellBoundaryAdhesionEnergyParameter;

    bool mUseMyDivisionRuleAlongWithModifier;
    double mDivisionTime;

    bool mWriteGroupNumberToCell;

    bool mMarkLeadingCells;
    double mLamellipodiumMaturationRate;
    double mLamellipodiumDestructionRate;

    bool mIfOutputModifierInformation;

public:

    /**
     * Default constructor.
     */
    FaceValueAndStressStateModifier();

    /**
     * Destructor.
     */
    virtual ~FaceValueAndStressStateModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);
    
    // my chagnes
    void UpdateFaceValuesAndStressStates(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void UpdateUnifiedEdgeMyosinActivtyOfFace(MutableVertexMesh<DIM, DIM>* pMesh, unsigned faceIndex);

    void UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(MutableVertexMesh<DIM, DIM>* pMesh, unsigned faceIndex);

    void UpdateStressStateOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);
    
    void UpdateCellAreas(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
    
    void UpdateCellAreaOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);

    void SetupSolveForCellDivision(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void UpdateForCellDivision(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void UpdateGroupNumbers(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
    
    void UpdateGroupNumberOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);

    void SetupSolveForLamellipodiumInfoOfCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void UpdateLamellipodiumInfoOfCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    // feedback form:
    void SetConsiderFeedbackOfFaceValues(bool ifConsiderFeedbackOfFaceValues)
    {
      mIfConsiderFeedbackOfFaceValues = ifConsiderFeedbackOfFaceValues;
    }

    void SetConsiderFeedbackOfFaceValuesOnlyForBoundaryCells(bool ifConsiderFeedbackOfFaceValuesOnlyForBoundaryCells)
    {
      mIfConsiderFeedbackOfFaceValuesOnlyForBoundaryCells = ifConsiderFeedbackOfFaceValuesOnlyForBoundaryCells;
    }

    void SetConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells(bool ifConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells)
    {
      mIfConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells = ifConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells;
    }

    void SetConsiderFeedbackOfCellCellAdhesion(bool ifConsiderFeedbackOfCellCellAdhesion)
    {
      mIfConsiderFeedbackOfCellCellAdhesion = ifConsiderFeedbackOfCellCellAdhesion;
    }

    // if feedback includes only strengthening or not
    void SetEMADontDecrease_CCADontDecrease_HasAThreshold_Threshold(
        bool EMADontDecrease,
        bool CCADontDecrease,
        bool CCAIncreasingHasAThresholdOfEdgeLength,
        double CCAIncreasingThresholdOfEdgeLengthPercentage)
    {
      mEMADontDecreaseWhenEdgeShrink = EMADontDecrease;
      mCCADontDecreaseWhenEdgeExpand = CCADontDecrease;
      mCCAIncreasingHasAThresholdOfEdgeLength = CCAIncreasingHasAThresholdOfEdgeLength;
      mCCAIncreasingThresholdOfEdgeLengthPercentage = CCAIncreasingThresholdOfEdgeLengthPercentage;
    }

    // detailed feedback information:
    void SetEdgeLengthAtRest(double edgeLengthAtRest)
    {
      this->mEdgeLengthAtRest = edgeLengthAtRest;
    }

    void SetFeedbackStrengthForMyosinActivity(double feedbackStrengthForMyosinActivity)
    {
      this->mFeedbackStrengthForMyosinActivity = feedbackStrengthForMyosinActivity;
    }

    void SetHillCoefficientForMyosinActivity(double hillCoefficientForMyosinActivity)
    {
      this->mHillCoefficientForMyosinActivity = hillCoefficientForMyosinActivity;
    }

    void SetFeedbackStrengthForAdhesion(double feedbackStrengthForAdhesion)
    {
      mFeedbackStrengthForAdhesion = feedbackStrengthForAdhesion;
    }
    
    void SetHillCoefficientForAdhesion(double hillCoefficientForAdhesion)
    {
      mHillCoefficientForAdhesion = hillCoefficientForAdhesion;
    }

    // stress state
    void SetCalculateStressStateBoolean(bool ifCalculateStressState)
    {
      mIfCalculateStressState = ifCalculateStressState;
    }

    void SetFixedTargetArea(double fixedTargetArea)
    {
      this->mFixedTargetArea = fixedTargetArea;
    }

    void SetFixedTargetPerimeter(double fixedTargetPerimeter)
    {
      this->mFixedTargetPerimeter = fixedTargetPerimeter;
    }

    void SetCaseNumberOfMembraneSurfaceEnergyForm(double caseNumberOfMembraneSurfaceEnergyForm)
    {
      this->mCaseNumberOfMembraneSurfaceEnergyForm = caseNumberOfMembraneSurfaceEnergyForm;
    }

    void SetNagaiHondaDeformationEnergyParameter(double nagaiHondaDeformationEnergyParameter)
    {
      this->mNagaiHondaDeformationEnergyParameter = nagaiHondaDeformationEnergyParameter;
    }

    void SetNagaiHondaMembraneSurfaceEnergyParameter(double nagaiHondaMembraneSurfaceEnergyParameter)
    {
      this->mNagaiHondaMembraneSurfaceEnergyParameter = nagaiHondaMembraneSurfaceEnergyParameter;
    }

    void SetNagaiHondaCellCellAdhesionEnergyParameter(double nagaiHondaCellCellAdhesionEnergyParameter)
    {
      this->mNagaiHondaCellCellAdhesionEnergyParameter = nagaiHondaCellCellAdhesionEnergyParameter;
    }

    void SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double nagaiHondaCellBoundaryAdhesionEnergyParameter)
    {
      this->mNagaiHondaCellBoundaryAdhesionEnergyParameter = nagaiHondaCellBoundaryAdhesionEnergyParameter;
    }
    
    // cell division
    void SetUseMyDivisionRuleAlongWithModifier(bool useMyDivisionRuleAlongWithModifier)
    {
      mUseMyDivisionRuleAlongWithModifier = useMyDivisionRuleAlongWithModifier;
    }

    void SetDivisionTime(double divisionTime)
    {
      mDivisionTime = divisionTime;
    }

    // group number
    void SetWriteGroupNumberToCell(bool writeGroupNumberToCell)
    {
      mWriteGroupNumberToCell = writeGroupNumberToCell;
    }

    // new SSA distribution rule
    void SetMarkLeadingCells(bool markLeadingCells)
    {
      mMarkLeadingCells = markLeadingCells;
    }

    void SetLamellipodiumMaturationRate(double lamellipodiumMaturationRate)
    {
      mLamellipodiumMaturationRate = lamellipodiumMaturationRate;
    }

    void SetLamellipodiumDestructionRate(double lamellipodiumDestructionRate)
    {
      mLamellipodiumDestructionRate = lamellipodiumDestructionRate;
    }

    // output information
    void SetOutputModifierInformationBoolean(bool ifOutputModifierInformation)
    {
      mIfOutputModifierInformation = ifOutputModifierInformation;
    }

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(FaceValueAndStressStateModifier)

#endif /*FACEVALUEANDSTRESSSTATEMODIFIER_HPP_*/
