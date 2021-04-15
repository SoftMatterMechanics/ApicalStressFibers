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
#include "VertexBasedCellPopulation.hpp"

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

    bool   mIfConsiderFeedbackOfFaceValues;
    bool   mIfConsiderFeedbackOfFaceValuesOnlyForBoundaryCells;
    bool   mIfConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells;
    bool   mApplyFeedbackOfFaceValuesToTopBoundaryCellsAndCellsAboveReservior;
    double mStripWidth;
    double mStripDistance; // change made by Chao
    double mStripStartXLocation;
    double mStripStartYLocation;
    bool   mIfConsiderFeedbackOfCellCellAdhesion;

    bool   mEMADontDecreaseWhenEdgeShrink;
    bool   mCCADontDecreaseWhenEdgeExpand;
    bool   mCCAIncreasingHasAThresholdOfEdgeLength;
    double mCCAIncreasingThresholdOfEdgeLengthPercentage;

    bool   mEMADontDecreaseBelowAThreshold;
    double mEMADontDecreaseBelowThisThreshold;

    double mEdgeLengthAtRest;
    double mKmForMyosinFeedback;
    double mFeedbackRateForMyosinActivity;
    double mHillPowerForMyosinActivity;
    bool   mCellCellAdhesionDontDecrease;
    double mKsForAdhesionFeedback;
    double mFeedbackRateForAdhesion;
    double mHillPowerForAdhesion;
    double mReferenceStress;

    bool   mIfCalculateStressState;
    bool   mIfSetCellDataOfEachForceContributions;
    unsigned mCaseNumberOfMembraneSurfaceEnergyForm;
    bool   mUseFixedTargetArea;
    double mFixedTargetArea;
    double mFixedTargetPerimeter;
    double mNagaiHondaDeformationEnergyParameter;
    double mNagaiHondaMembraneSurfaceEnergyParameter;
    double mNagaiHondaCellCellAdhesionEnergyParameter;
    double mNagaiHondaCellBoundaryAdhesionEnergyParameter;

    bool   mUseMyDivisionRuleAlongWithModifier;
    double mDivisionTime;

    bool   mWriteGroupNumberToCell;
    bool   mWriteVertexVelocityAndForceToCellData;
    bool   mWriteForcesFromNeighboringCellsToCellData;
    
    bool   mMarkLeadingCells;
    bool   mMultipleLeadingCells;
    unsigned mLeadingCellNumber;
    double mLamellipodiumMaturationRate;
    double mLamellipodiumDestructionRate;
    bool   mIfEquilibrateForAWhile;
    double mEndTimeForEquilibrium;

    bool   mIfOutputModifierInformation;

    bool   mHasMyosinActivityDepression;
    double mMyosinActivityDepressedTime;
    double mMyosinActivityDepressingRate;

    double mTimeForChangingFeedback;
    double mChangedKmForMyosinFeedback;
    double mChangedFeedbackRate;
    double mChangedMyosinActivityBaseValue;

    // start of change made by Chao
    double mSmallChangeForAreaCalculation;

    double mStripSubstrateAdhesionParameter;
    double mReservoirSubstrateAdhesionParameter;   

    double mWidth;
    double mCenterOfWidth;
    bool   mConsiderConsistencyForSSA;
    // end of change made by Chao

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
    void UpdateMyosinActivityAndCellCellAdhesionAndStressState(AbstractCellPopulation<DIM,DIM>& rCellPopulation); // change made by Chao

    void UpdateStressStateOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);
    
    void UpdateCellDataOfForcesFromNeighboringCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);
    
    void UpdateMyosinActivtyOfElement(MutableVertexMesh<DIM, DIM>* pMesh, unsigned elem_index);  // change made by Chao

    void UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned faceIndex);

    void SetIfEquilibrateForAWhile(bool ifEquilibrateForAWhile)
    {
      mIfEquilibrateForAWhile = ifEquilibrateForAWhile;
    }

    void SetEndTimeForEquilibrium(double endTimeForEquilibrium)
    {
      mEndTimeForEquilibrium = endTimeForEquilibrium;
    }

  // start of change made by Chao
    c_vector<double, DIM> GetStripSubstrateAdhesionAreaGradientOfElementAtNode(AbstractCellPopulation<DIM>& rCellPopulation,VertexElement<DIM,DIM>* pElement, unsigned localIndex);

    c_vector<double, DIM> GetReservoirSubstrateAdhesionAreaGradientOfElementAtNode(AbstractCellPopulation<DIM>& rCellPopulation, VertexElement<DIM,DIM>* pElement, unsigned localIndex);

    void InitializeShapeTensorOfCell(CellPtr pCell);

    void InitializeStressOfCell(CellPtr pCell);

    void InitializeMyosinActivity(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);

  //end of change made by Chao

    void UpdateCellAreas(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
    
    void UpdateCellAreaOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);

    void SetupSolveForCellDivision(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void UpdateForCellDivision(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void UpdateGroupNumbers(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
    
    void UpdateGroupNumberOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell);

    void SetupLeaderCellAtTheEndOfEquilibrium(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

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

    void SetApplyFeedbackOfFaceValuesToTopBoundaryCellsAndCellsAboveReservior(bool applyFeedbackOfFaceValuesToTopBoundaryCellsAndCellsAboveReservior)
    {
      mApplyFeedbackOfFaceValuesToTopBoundaryCellsAndCellsAboveReservior = applyFeedbackOfFaceValuesToTopBoundaryCellsAndCellsAboveReservior;
    }

    void SetStripWidth(double stripWidth)
    {
      mStripWidth = stripWidth;
    }

    // start of change made by Chao
    void SetStripDistance(double stripDistance)   
    {
      mStripDistance = stripDistance;
    }
    void SetWidth(double width)
    {
      mWidth = width;
    }
    void SetCenterOfWidth(double centerOfWidth)
    {
      mCenterOfWidth = centerOfWidth;
    }
    // end of change made by Chao

    void SetStripStartXLocation(double stripStartXLocation)
    {
      mStripStartXLocation = stripStartXLocation;
    }

    void SetStripStartYLocation(double stripStartYLocation)
    {
      mStripStartYLocation = stripStartYLocation;
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

    void SetEMADontDecreaseBelowAThreshold_ThisThreshold(bool EMADontDecreaseBelowAThreshold, double EMADontDecreaseBelowThisThreshold)
    {
      mEMADontDecreaseBelowAThreshold = EMADontDecreaseBelowAThreshold;
      mEMADontDecreaseBelowThisThreshold = EMADontDecreaseBelowThisThreshold;
    }


    // detailed feedback information:
    void SetEdgeLengthAtRest(double edgeLengthAtRest)
    {
      this->mEdgeLengthAtRest = edgeLengthAtRest;
    }

    void SetKmForMyosinFeedback(double kmForMyosinFeedback)
    {
      mKmForMyosinFeedback = kmForMyosinFeedback;
    }

    void SetFeedbackRateForMyosinActivity(double feedbackRateForMyosinActivity)
    {
      this->mFeedbackRateForMyosinActivity = feedbackRateForMyosinActivity;
    }

    void SetHillPowerForMyosinActivity(double hillPowerForMyosinActivity)
    {
      this->mHillPowerForMyosinActivity = hillPowerForMyosinActivity;
    }

    void SetCellCellAdhesionDontDecrease(bool cellCellAdhesionDontDecrease)
    {
      mCellCellAdhesionDontDecrease = cellCellAdhesionDontDecrease;
    }

    void SetKsForAdhesionFeedback(double ksForAdhesionFeedback)
    {
      mKsForAdhesionFeedback = ksForAdhesionFeedback;
    }

    void SetFeedbackRateForAdhesion(double feedbackRateForAdhesion)
    {
      mFeedbackRateForAdhesion = feedbackRateForAdhesion;
    }
    
    void SetHillPowerForAdhesion(double hillPowerForAdhesion)
    {
      mHillPowerForAdhesion = hillPowerForAdhesion;
    }

    void SetReferenceStress(double referenceStress)
    {
      mReferenceStress = referenceStress;
    }

    // stress state
    void SetCalculateStressStateBoolean(bool ifCalculateStressState)
    {
      mIfCalculateStressState = ifCalculateStressState;
    }

    void SetIfSetCellDataOfEachForceContributions(bool ifSetCellDataOfEachForceContributions)
    {
      mIfSetCellDataOfEachForceContributions = ifSetCellDataOfEachForceContributions;
    }

    void SetFixedTargetArea(double fixedTargetArea)
    {
      this->mFixedTargetArea = fixedTargetArea;
    }

    void SetUseFixedTargetArea(bool useFixedTargetArea)
    {
      mUseFixedTargetArea = useFixedTargetArea;
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

    void SetMultipleLeadingCells(bool multipleLeadingCells)
    {
      mMultipleLeadingCells = multipleLeadingCells;
    }

    void SetLeadingCellNumber(unsigned leadingCellNumber)
    {
      mLeadingCellNumber = leadingCellNumber;
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

    void SetHasMyosinActivityDepression(bool hasMyosinActivityDepression)
    {
      mHasMyosinActivityDepression = hasMyosinActivityDepression;
    }

    void SetMyosinActivityDepressedTime (double myosinActivityDepressedTime)
    {
      mMyosinActivityDepressedTime = myosinActivityDepressedTime;
    }

    void SetMyosinActivityDepressingRate (double myosinActivityDepressingRate)
    {
      mMyosinActivityDepressingRate = myosinActivityDepressingRate;
    }

    // for changing feedback after a particulat time
    void SetTimeForChangingFeedback(double timeForChangingFeedback)
    {
      mTimeForChangingFeedback = timeForChangingFeedback;
    }

    void SetChangedKmForMyosinFeedback(double changedKmForMyosinFeedback)
    {
      mChangedKmForMyosinFeedback = changedKmForMyosinFeedback;
    }

    void SetChangedFeedbackRate(double changedFeedbackRate)
    {
      mChangedFeedbackRate = changedFeedbackRate;
    }

    void SetChangedMyosinActivityBaseValue(double changedMyosinActivityBaseValue)
    {
      mChangedMyosinActivityBaseValue = changedMyosinActivityBaseValue;
    }

    void SetSmallChangeForAreaCalculation(double smallChangeForAreaCalculation)   // change made by Chao
    {
      mSmallChangeForAreaCalculation = smallChangeForAreaCalculation;
    }

    void SetStripSubstrateAdhesionParameter(double stripSubstrateAdhesionParameter)  // change made by Chao
    {
      mStripSubstrateAdhesionParameter = stripSubstrateAdhesionParameter;
    }

    void SetReservoirSubstrateAdhesionParameter(double reservoirSubstrateAdhesionParameter)  // change made by Chao
    {
      mReservoirSubstrateAdhesionParameter = reservoirSubstrateAdhesionParameter;
    }

    void SetConsiderConsistencyForSSA(bool considerConsistencyForSSA)    // change made by Chao
    {
      mConsiderConsistencyForSSA = considerConsistencyForSSA;
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
