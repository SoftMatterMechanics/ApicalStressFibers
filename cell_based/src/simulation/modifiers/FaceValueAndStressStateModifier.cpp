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

#include "FaceValueAndStressStateModifier.hpp"
#include "MyXToroidal2dVertexMesh.hpp"

template<unsigned DIM>
FaceValueAndStressStateModifier<DIM>::FaceValueAndStressStateModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mIfConsiderFeedbackOfElementMyosinActivity(false),
      mIfConsiderFeedbackOfFaceValuesOnlyForBoundaryCells(false),
      mIfConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells(false),
      mApplyFeedbackOfFaceValuesToTopBoundaryCellsAndCellsAboveReservior(false),
      mStripWidth(1.0),
      mStripStartXLocation(0.0),
      mStripStartYLocation(0.0),
      mIfConsiderFeedbackOfCellCellAdhesion(false),

      mEMADontDecreaseWhenEdgeShrink(false),
      mCCADontDecreaseWhenEdgeExpand(false),
      mCCAIncreasingHasAThresholdOfEdgeLength(false),
      mCCAIncreasingThresholdOfEdgeLengthPercentage(0.5),

      mEMADontDecreaseBelowAThreshold(false),
      mEMADontDecreaseBelowThisThreshold(0.0),

      mEdgeLengthAtRest(sqrt(M_PI/(6*sqrt(3)/4))),
      mKmForMyosinFeedback(1.0),
      mFeedbackRateForMyosinActivity(1.0),
      mHillPowerForMyosinActivity(8.0),
      mKsForAdhesionFeedback(1.0),
      mFeedbackRateForAdhesion(1.0),
      mHillPowerForAdhesion(8.0),
      mReferenceStress(1.0),

      mIfCalculateStressState(true),
      mIfSetCellDataOfEachForceContributions(false),
      mCaseNumberOfMembraneSurfaceEnergyForm(0),
      mUseFixedTargetArea(true),
      mFixedTargetArea(M_PI),
      mFixedTargetPerimeter(6*sqrt(M_PI/(6*sqrt(3)/4))),
      mNagaiHondaDeformationEnergyParameter(1.0),
      mNagaiHondaMembraneSurfaceEnergyParameter(0.1),
      mNagaiHondaCellCellAdhesionEnergyParameter(0.0),
      mNagaiHondaCellBoundaryAdhesionEnergyParameter(0.0),

      mUseMyDivisionRuleAlongWithModifier(false),
      mDivisionTime(DOUBLE_UNSET),

      mWriteGroupNumberToCell(false),
      mWriteVertexVelocityAndForceToCellData(false),
      mWriteForcesFromNeighboringCellsToCellData(false),
      
      mMarkLeadingCells(false),
      mMultipleLeadingCells(false),
      mLeadingCellNumber(1),

      mIfOutputModifierInformation(false),

      mTimeForChangingFeedback(DOUBLE_UNSET),
      mChangedKmForMyosinFeedback(1.0),
      mChangedFeedbackRate(0.0),
      mChangedMyosinActivityBaseValue(1.0)

{
}

template<unsigned DIM>
FaceValueAndStressStateModifier<DIM>::~FaceValueAndStressStateModifier()
{
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateMyosinActivityAndCellCellAdhesionAndStressState(rCellPopulation);
    UpdateCellAreas(rCellPopulation);
    if (mUseMyDivisionRuleAlongWithModifier)
        UpdateForCellDivision(rCellPopulation);
    if (mWriteGroupNumberToCell)
        UpdateGroupNumbers(rCellPopulation);
    if (mMarkLeadingCells)
        UpdateLamellipodiumInfoOfCells(rCellPopulation);
    if (mIfEquilibrateForAWhile)
        SetupLeaderCellAtTheEndOfEquilibrium(rCellPopulation);
    // if (true)
    //     UpdateLamellipodiumInfoOfCells
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        InitializeShapeTensorOfCell(*cell_iter);
        InitializeStressOfCell(*cell_iter);
        InitializeMyosinActivity(rCellPopulation, *cell_iter);
    }
    UpdateMyosinActivityAndCellCellAdhesionAndStressState(rCellPopulation);
    UpdateCellAreas(rCellPopulation);
    if (mUseMyDivisionRuleAlongWithModifier)
        SetupSolveForCellDivision(rCellPopulation);
    if (mWriteGroupNumberToCell)
        UpdateGroupNumbers(rCellPopulation);
    if (mMarkLeadingCells)
        SetupSolveForLamellipodiumInfoOfCells(rCellPopulation);
    if (mIfEquilibrateForAWhile)
        SetupLeaderCellAtTheEndOfEquilibrium(rCellPopulation);
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::InitializeShapeTensorOfCell(CellPtr pCell)
{
    double target_area = 0.0;
    if (mUseFixedTargetArea)
        target_area = this->mFixedTargetArea;
    else
        target_area = pCell->GetCellData()->GetItem("target area");
    
    pCell->GetCellData()->SetItem("shape_tensor_xx",target_area/(3*sqrt(3)));
    pCell->GetCellData()->SetItem("shape_tensor_yy",target_area/(3*sqrt(3)));
    pCell->GetCellData()->SetItem("shape_tensor_xy",0);
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::InitializeStressOfCell(CellPtr pCell)
{
    // double feedback_rate = this->mFeedbackRateForAdhesion;
    // assert(feedback_rate!=0);
    // double q = this->mHillPowerForAdhesion;
    // double Ks = this->mKsForAdhesionFeedback;
    // double pho = 2.0*feedback_rate;
    // double ksi = 1.0*feedback_rate;
    // double reference_stress = mReferenceStress;

    pCell->GetCellData()->SetItem("Stress1",mReferenceStress);
    pCell->GetCellData()->SetItem("Stress2",mReferenceStress);
    pCell->GetCellData()->SetItem("PrincipalAxisOfStress",0.0);
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::InitializeMyosinActivity(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell)
{
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }

    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
    VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(elem_index);

    p_element->SetElementMyosinActivity(1.0);
}

/*
template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateFaceValuesAndStressStates(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (mIfConsiderFeedbackOfFaceValues)
    {
        //loop over the face values
        if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
        {
            EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
        }
        MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());
        for (unsigned face_index = 0; face_index < p_mesh->GetNumFaces(); face_index++)
        {
            if (!p_mesh->GetFace(face_index)->IsDeleted())
            {
                if (mIfConsiderFeedbackOfFaceValuesOnlyForBoundaryCells)
                {
                    if (p_mesh->IsFaceContainedByABoundaryElement(face_index))
                    {
                        if (mIfConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells)
                        {
                            if (p_mesh->IsFaceContainedByATopBoundaryElement(face_index))
                            {
                                this->UpdateUnifiedEdgeMyosinActivtyOfFace(p_mesh, face_index);

                                if (mIfConsiderFeedbackOfCellCellAdhesion)
                                    this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                            }

                        }
                        else
                        {
                            this->UpdateUnifiedEdgeMyosinActivtyOfFace(p_mesh, face_index);

                            if (mIfConsiderFeedbackOfCellCellAdhesion)
                                this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                        }
                        
                    }
                }
                else
                {
                    if (mApplyFeedbackOfFaceValuesToTopBoundaryCellsAndCellsAboveReservior)
                    {
                        VertexElement<DIM-1, DIM>* p_face = p_mesh->GetFace(face_index);
                        double face_y_location = 0.5*( p_face->GetNode(0)->rGetLocation()[1] + p_face->GetNode(1)->rGetLocation()[1] );

                        if (p_mesh->IsFaceContainedByATopBoundaryElement(face_index) || face_y_location>mStripStartYLocation)
                        {
                            this->UpdateUnifiedEdgeMyosinActivtyOfFace(p_mesh, face_index);

                            if (mIfConsiderFeedbackOfCellCellAdhesion)
                                this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                        }

                    }
                    else
                    {
                        this->UpdateUnifiedEdgeMyosinActivtyOfFace(p_mesh, face_index);
                        
                        if (mIfConsiderFeedbackOfCellCellAdhesion)
                            this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                    }
                    
                }
            }
        }
    }
*/


template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateMyosinActivityAndCellCellAdhesionAndStressState(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{   
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }

    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    if (mIfConsiderFeedbackOfElementMyosinActivity)
    {
        // loop over the elements
        for (unsigned elem_index = 0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            this->UpdateMyosinActivtyOfElement(p_mesh, elem_index);
        }
    }

    if (mIfConsiderFeedbackOfCellCellAdhesion)
    {
        // loop over the faces
        for (unsigned face_index = 0; face_index < p_mesh->GetNumFaces(); face_index++)
        {
            VertexElement<DIM-1, DIM>* pFace = p_mesh->GetFace(face_index);
            if (!pFace->IsDeleted())
                this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(rCellPopulation, face_index);
        }
    }

    if (mIfCalculateStressState)
    {
        // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
        for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
            cell_iter != rCellPopulation.rGetCells().end();
            ++cell_iter)
        {
            UpdateStressStateOfCell(rCellPopulation, *cell_iter);
        }

        if (mWriteForcesFromNeighboringCellsToCellData)
        {
            for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
                cell_iter != rCellPopulation.rGetCells().end();
                ++cell_iter)
            {
                UpdateCellDataOfForcesFromNeighboringCell(rCellPopulation, *cell_iter);
            }
        }

    }

}

/*
template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateStressStateOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell)
{   
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }

    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    if (mCaseNumberOfMembraneSurfaceEnergyForm >= 0u)
    {
        unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
        VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(elem_index);

        double area = p_mesh->GetVolumeOfElement(elem_index);
        double target_area = 0.0;
        if (mUseFixedTargetArea)
            target_area = mFixedTargetArea;
        else
            target_area = pCell->GetCellData()->GetItem("target area");

        double Ga= this->mNagaiHondaMembraneSurfaceEnergyParameter;

        double sum_XX = 0.0;
        double sum_YY = 0.0;
        double sum_XY = 0.0;
        double sum_adhe_XX = 0.0;
        double sum_adhe_YY = 0.0;
        double sum_adhe_XY = 0.0;
        double sum_shape_tensor_XX = 0.0;
        double sum_shape_tensor_YY = 0.0;
        double sum_shape_tensor_XY = 0.0;

        double myosin_weighted_perimeter = 0.0;
        unsigned local_index_lowest_vertex = 0;
        double y_coor_lowest_vertex = 10000.0;

        c_vector<double,DIM> centroid_relate_to_first_node = zero_vector<double>(DIM);

        double myosin_activity_max = 0.0;
        double myosin_activity_average = 0.0;

        for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
        {
            Node<DIM>* pNode = p_element->GetNode(local_index);
            if (pNode->rGetLocation()[1]<y_coor_lowest_vertex)
            {
                local_index_lowest_vertex = local_index;
                y_coor_lowest_vertex = pNode->rGetLocation()[1];
            }
            Node<DIM>* pNodeA = p_element->GetNode(local_index);
            Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
            double l_ab = p_mesh->GetDistanceBetweenNodes(pNodeA->GetIndex(), pNodeB->GetIndex());

            if (mIfConsiderFeedbackOfFaceValues)
            {
                VertexElement<DIM-1,  DIM>* pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeA->GetIndex(), pNodeB->GetIndex()));
                double m_ab = pFace->GetUnifiedEdgeMyosinActivty();
                myosin_activity_max = std::max(myosin_activity_max, m_ab);
                myosin_activity_average += 1.0/p_element->GetNumNodes()*m_ab;
                myosin_weighted_perimeter += sqrt(m_ab)*l_ab;
            }
            else
                myosin_weighted_perimeter += l_ab;

            centroid_relate_to_first_node += 1.0/p_element->GetNumNodes() * p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(), pNode->rGetLocation());
        }
        pCell->GetCellData()->SetItem("LocalIndexOfLowestVertex", local_index_lowest_vertex);

        for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
        {
            //  --B
            //    |
            // C--A
            Node<DIM>* pNodeC = p_element->GetNode((local_index-1+p_element->GetNumNodes())%p_element->GetNumNodes());
            Node<DIM>* pNodeA = p_element->GetNode(local_index);
            Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
            c_vector<double,DIM> location_c = pNodeC->rGetLocation();
            c_vector<double,DIM> location_a = pNodeA->rGetLocation();
            c_vector<double,DIM> location_b = pNodeB->rGetLocation();
            double l_ca = p_mesh->GetDistanceBetweenNodes(pNodeC->GetIndex(), pNodeA->GetIndex());
            double l_ab = p_mesh->GetDistanceBetweenNodes(pNodeA->GetIndex(), pNodeB->GetIndex());
            c_vector<double,DIM> vec_ca =  p_mesh->GetVectorFromAtoB(location_c, location_a);
            c_vector<double,DIM> vec_ab =  p_mesh->GetVectorFromAtoB(location_a, location_b);

            VertexElement<DIM-1,  DIM>* pFace0 = nullptr;
            VertexElement<DIM-1,  DIM>* pFace = nullptr;
            if (mIfConsiderFeedbackOfFaceValues)
            {
                pFace0 = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeC->GetIndex(), pNodeA->GetIndex()));                    
                pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeA->GetIndex(), pNodeB->GetIndex()));
            }

            std::set<unsigned> elements_containing_nodeC = pNodeC->rGetContainingElementIndices();
            std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
            std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
            // Find common elements
            std::set<unsigned> shared_elements0;
            std::set_intersection(elements_containing_nodeC.begin(),
                                elements_containing_nodeC.end(),
                                elements_containing_nodeA.begin(),
                                elements_containing_nodeA.end(),
                                std::inserter(shared_elements0, shared_elements0.begin()));
            // Check that the nodes have a common edge
            assert(!shared_elements0.empty());
            if (shared_elements0.size() >= 3)
            {
                std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateStressStateOfCell";
                std::cout<< std::endl << "Get shared elements more than 2";
            }
            // Find common elements
            std::set<unsigned> shared_elements;
            std::set_intersection(elements_containing_nodeA.begin(),
                                elements_containing_nodeA.end(),
                                elements_containing_nodeB.begin(),
                                elements_containing_nodeB.end(),
                                std::inserter(shared_elements, shared_elements.begin()));
            // Check that the nodes have a common edge
            assert(!shared_elements.empty());
            if (shared_elements.size() >= 3)
            {
                std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateStressStateOfCell";
                std::cout<< std::endl << "Get shared elements more than 2";
            }

            std::string name_item;
            std::ostringstream oss;

            // output each force contributions on the vertex from the element
            // Pressue from area term
            c_vector<double,DIM> vec_p;
            vec_p[0] = 0.5*(vec_ca+vec_ab)[1];
            vec_p[1] = -0.5*(vec_ca+vec_ab)[0];
            c_vector<double,DIM> F_Press = (-area+target_area)*vec_p;
            if (mIfSetCellDataOfEachForceContributions)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Press_x";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Press[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Press_y";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Press[1]);
            }

            // Tension from perimeter term
            double m_ca = 1.0;
            double m_ab = 1.0;
            if (mIfConsiderFeedbackOfFaceValues)
            {
                m_ca = pFace0->GetUnifiedEdgeMyosinActivty();
                m_ab = pFace->GetUnifiedEdgeMyosinActivty();
            }
            c_vector<double,DIM> myosin_biased_vec_q = sqrt(m_ab)*vec_ab/l_ab - sqrt(m_ca)*vec_ca/l_ca;
            c_vector<double,DIM> F_Tens = Ga*myosin_weighted_perimeter*myosin_biased_vec_q;
            if (mIfSetCellDataOfEachForceContributions)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Tens_x";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Tens[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Tens_y";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Tens[1]);
            }

            // Cell-cell adhesion
            double lambda_ca = 0.0;
            double lambda_ab = 0.0;
            if (shared_elements0.size() == 1)
                lambda_ca = mNagaiHondaCellBoundaryAdhesionEnergyParameter;
            else
                lambda_ca = mNagaiHondaCellCellAdhesionEnergyParameter;
            if (shared_elements.size() == 1)
                lambda_ab = mNagaiHondaCellBoundaryAdhesionEnergyParameter;
            else
                lambda_ab = mNagaiHondaCellCellAdhesionEnergyParameter;
            if (mIfConsiderFeedbackOfFaceValues)
            {
                lambda_ca *= pFace0->GetUnifiedCellCellAdhesionEnergyParameter();
                lambda_ab *= pFace->GetUnifiedCellCellAdhesionEnergyParameter();
            }
            c_vector<double,DIM> F_Adhe = (-lambda_ca)*vec_ca/l_ca +(-lambda_ab)*(-vec_ab)/l_ab;
            if (mIfSetCellDataOfEachForceContributions)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Adhe_x";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Adhe[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Adhe_y";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Adhe[1]);
            }            

            double Fx = F_Press[0] +F_Tens[0] +F_Adhe[0];
            double Fy = F_Press[1] +F_Tens[1] +F_Adhe[1];
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_Fx";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, Fx);
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_Fy";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, Fy);

            // output velocity of vertices
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_vx";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, pNodeA->rGetAppliedForce()[0]);
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_vy";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, pNodeA->rGetAppliedForce()[1]);

            // output face values of edges
            if (mIfConsiderFeedbackOfFaceValues)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_edge_myo";
                name_item = oss.str();
            //    pCell->GetCellData()->SetItem(name_item, pFace->GetUnifiedEdgeMyosinActivity());
                pCell->GetCellData()->SetItem(name_item, p_element->GetElementMyosinActivity());
            }
            if (mIfConsiderFeedbackOfFaceValues && mIfConsiderFeedbackOfCellCellAdhesion)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_unified_cc";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, pFace->GetUnifiedCellCellAdhesionEnergyParameter());
            }

            // next: stress calculation
            double lambda = 0.0;
            if (shared_elements.size() == 1)
                lambda = mNagaiHondaCellBoundaryAdhesionEnergyParameter;
            else
                lambda = mNagaiHondaCellCellAdhesionEnergyParameter;
            if (mIfConsiderFeedbackOfFaceValues)
                lambda *= pFace->GetUnifiedCellCellAdhesionEnergyParameter();
            sum_adhe_XX += lambda*l_ab*(vec_ab[0]/l_ab)*(vec_ab[0]/l_ab);
            sum_adhe_YY += lambda*l_ab*(vec_ab[1]/l_ab)*(vec_ab[1]/l_ab);
            sum_adhe_XY += lambda*l_ab*(vec_ab[0]/l_ab)*(vec_ab[1]/l_ab);

            double myo_ab = 1.0;
            if (mIfConsiderFeedbackOfFaceValues)
                myo_ab = pFace->GetUnifiedEdgeMyosinActivty();
            sum_XX += sqrt(myo_ab)*(vec_ab[0]/l_ab)*vec_ab[0];
            sum_YY += sqrt(myo_ab)*(vec_ab[1]/l_ab)*vec_ab[1];
            sum_XY += sqrt(myo_ab)*(vec_ab[0]/l_ab)*vec_ab[1];
            c_vector<double, DIM> node_location_ralate_to_centroid = p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(), pNodeA->rGetLocation()) - centroid_relate_to_first_node;

            sum_shape_tensor_XX += node_location_ralate_to_centroid[0]*node_location_ralate_to_centroid[0];
            sum_shape_tensor_YY += node_location_ralate_to_centroid[1]*node_location_ralate_to_centroid[1];
            sum_shape_tensor_XY += node_location_ralate_to_centroid[0]*node_location_ralate_to_centroid[1];

        } // end of iteration of vertices of the element, for calculation of force contributions and stress summation.

        double stress_XX = (area-target_area) + Ga/area*myosin_weighted_perimeter*sum_XX + 1/area*sum_adhe_XX;
        double stress_YY = (area-target_area) + Ga/area*myosin_weighted_perimeter*sum_YY + 1/area*sum_adhe_YY;
        double stress_XY = Ga/area*myosin_weighted_perimeter*sum_XY + 1/area*sum_adhe_XY;
        double R = sqrt( pow((stress_XX-stress_YY)/2, 2) + pow(stress_XY, 2) );
        double stress_1 = (stress_XX + stress_YY)/2 + R;
        double stress_2 = (stress_XX + stress_YY)/2 - R;
    
        pCell->GetCellData()->SetItem("StressXX", stress_XX);
        pCell->GetCellData()->SetItem("StressYY", stress_YY);
        pCell->GetCellData()->SetItem("StressXY", stress_XY);
        pCell->GetCellData()->SetItem("Stress1", stress_1);
        pCell->GetCellData()->SetItem("Stress2", stress_2);

        double shape_XX = 1.0/p_element->GetNumNodes()*sum_shape_tensor_XX;
        double shape_YY = 1.0/p_element->GetNumNodes()*sum_shape_tensor_YY;
        double shape_XY = 1.0/p_element->GetNumNodes()*sum_shape_tensor_XY;
        assert((shape_XX + shape_YY + shape_XY) < DOUBLE_UNSET );

        // A=[xx, xy; xy, yy], A-xI = [xx-x, xy; xy, yy-x], |A-xI|=(xx-x)(yy-x)-xy^2=x^2 -(xx+yy)x +(xx*yy-xy*xy)=0
        // x1,x2=1/2*(xx+yy +/-sqrt(delta)); delta=(xx^2+yy^2-2xx*yy+4xy*xy)
        // A-x1*I=[xx-x1, xy;...]; A-x2*I=[];
        // let x1>=x2, vec1=[xy, x1-xx]/sqrt();
        double delta = (shape_XX-shape_YY)*(shape_XX-shape_YY) +4*shape_XY*shape_XY;
        double lambda_1 = 0.5*(shape_XX+shape_YY+sqrt(delta));
        double lambda_2 = 0.5*(shape_XX+shape_YY-sqrt(delta));
        double ratio = sqrt(lambda_1/lambda_2);
        pCell->GetCellData()->SetItem("AspectRatio", ratio);

        c_vector<double, DIM> eigen_vec_1 = zero_vector<double>(DIM);
        eigen_vec_1[0] = shape_XY/sqrt( shape_XY*shape_XY+(lambda_1-shape_XX)*(lambda_1-shape_XX) );
        eigen_vec_1[1] = (lambda_1-shape_XX)/sqrt( shape_XY*shape_XY+(lambda_1-shape_XX)*(lambda_1-shape_XX) );
        double principal_axis_of_shape = abs(atan(eigen_vec_1[1]/eigen_vec_1[0]));
        pCell->GetCellData()->SetItem("PrincipalAxisOfShape", principal_axis_of_shape);

        // double R_shape = sqrt( pow((shape_XX-shape_YY)/2, 2) + pow(shape_XY, 2) );
        // double shape_1 = (shape_XX + shape_YY)/2 + R;
        // double shape_2 = (shape_XX + shape_YY)/2 - R;
        // double aspect_ratio = sqrt(shape_1/shape_2);
        // double principal_axis_angle = shape_XY>0.0 ? -0.5*acos( (shape_XX-0.5*(shape_XX+shape_YY))/R_shape ) : 0.5*acos( (shape_XX-0.5*(shape_XX+shape_YY))/R_shape );
        // pCell->GetCellData()->SetItem("AspectRatio", aspect_ratio);
        // pCell->GetCellData()->SetItem("PrincipalAxis", principal_axis_angle);

        double principal_axis_stress = 0.5*acos( (stress_XX-0.5*(stress_XX+stress_YY))/R );
        pCell->GetCellData()->SetItem("PrincipalAxisOfStress", principal_axis_stress);

        if (mIfConsiderFeedbackOfFaceValues)
        {
            pCell->GetCellData()->SetItem("Myosin Activity Max", myosin_activity_max);
            pCell->GetCellData()->SetItem("Myosin Activity Averaged", myosin_activity_average);
        }
    }

}
*/

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateStressStateOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell)
{   
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }

    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
    VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(elem_index);

    // double dt = SimulationTime::Instance()->GetTimeStep();

    if (mCaseNumberOfMembraneSurfaceEnergyForm >= 0u)
    {
        double static_stress_XX = 0.0;
        double static_stress_YY = 0.0;
        double static_stress_XY = 0.0;
        double static_stress_YX = 0.0;
        double shape_tensor_XX = 0.0;
        double shape_tensor_YY = 0.0;
        double shape_tensor_XY = 0.0;

        double area = p_mesh->GetVolumeOfElement(elem_index);
        double perimeter = p_mesh->GetSurfaceAreaOfElement(elem_index);
        c_vector<double, DIM> centroid = p_cell_population->rGetMesh().GetCentroidOfElement(elem_index);

        double target_area = 0.0;
        if (mUseFixedTargetArea)
            target_area = mFixedTargetArea;
        else
            target_area = pCell->GetCellData()->GetItem("target area");

        double Ga= this->mNagaiHondaMembraneSurfaceEnergyParameter;

        unsigned local_index_lowest_vertex = 0;
        double y_coor_lowest_vertex = 10000.0;

        c_vector<double,DIM> centroid_relate_to_first_node = zero_vector<double>(DIM);

        for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
        {
            Node<DIM>* pNode = p_element->GetNode(local_index);
            if (pNode->rGetLocation()[1]<y_coor_lowest_vertex)
            {
                local_index_lowest_vertex = local_index;
                y_coor_lowest_vertex = pNode->rGetLocation()[1];
            }
//            centroid_relate_to_first_node += 1.0/p_element->GetNumNodes() * p_mesh->GetVectorFromAtoB(p_element->GetNode(0)->rGetLocation(), pNode->rGetLocation());
        }

        pCell->GetCellData()->SetItem("LocalIndexOfLowestVertex", local_index_lowest_vertex);

        for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
        {
            //  --B
            //    |
            // C--A
            Node<DIM>* pNodeC = p_element->GetNode((local_index-1+p_element->GetNumNodes())%p_element->GetNumNodes());
            Node<DIM>* pNodeA = p_element->GetNode(local_index);
            Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
            c_vector<double,DIM> location_c = pNodeC->rGetLocation();
            c_vector<double,DIM> location_a = pNodeA->rGetLocation();
            c_vector<double,DIM> location_b = pNodeB->rGetLocation();
            double l_ca = p_mesh->GetDistanceBetweenNodes(pNodeC->GetIndex(), pNodeA->GetIndex());
            double l_ab = p_mesh->GetDistanceBetweenNodes(pNodeA->GetIndex(), pNodeB->GetIndex());
            c_vector<double,DIM> vec_ca =  p_mesh->GetVectorFromAtoB(location_c, location_a);
            c_vector<double,DIM> vec_ab =  p_mesh->GetVectorFromAtoB(location_a, location_b);

            VertexElement<DIM-1,  DIM>* pFace0 = nullptr;
            VertexElement<DIM-1,  DIM>* pFace = nullptr;
            if (mIfConsiderFeedbackOfCellCellAdhesion)
            {
                pFace0 = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeC->GetIndex(), pNodeA->GetIndex()));                    
                pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeA->GetIndex(), pNodeB->GetIndex()));
            }

            std::set<unsigned> elements_containing_nodeC = pNodeC->rGetContainingElementIndices();
            std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
            std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
            // Find common elements
            std::set<unsigned> shared_elements0;
            std::set_intersection(elements_containing_nodeC.begin(),
                                elements_containing_nodeC.end(),
                                elements_containing_nodeA.begin(),
                                elements_containing_nodeA.end(),
                                std::inserter(shared_elements0, shared_elements0.begin()));
            // Check that the nodes have a common edge
            assert(!shared_elements0.empty());
            if (shared_elements0.size() >= 3)
            {
                std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateStressStateOfCell";
                std::cout<< std::endl << "Get shared elements more than 2";
            }
            // Find common elements
            std::set<unsigned> shared_elements;
            std::set_intersection(elements_containing_nodeA.begin(),
                                elements_containing_nodeA.end(),
                                elements_containing_nodeB.begin(),
                                elements_containing_nodeB.end(),
                                std::inserter(shared_elements, shared_elements.begin()));
            // Check that the nodes have a common edge
            assert(!shared_elements.empty());
            if (shared_elements.size() >= 3)
            {
                std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateStressStateOfCell";
                std::cout<< std::endl << "Get shared elements more than 2";
            }

            std::string name_item;
            std::ostringstream oss;

            // output each force contributions on the vertex from the element
            // Pressue from area term
            c_vector<double,DIM> element_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            c_vector<double,DIM> F_Press = -(area-target_area)*element_area_gradient;

            if (mIfSetCellDataOfEachForceContributions)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Press_x";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Press[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Press_y";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Press[1]);
            }

            // Tension from perimeter term
            double element_myosin_activity = 1.0;
            if (mIfConsiderFeedbackOfElementMyosinActivity)
                element_myosin_activity = p_element->GetElementMyosinActivity();

            c_vector<double,DIM> perimeter_gradient_at_node = p_cell_population->rGetMesh().GetPerimeterGradientOfElementAtNode(p_element, local_index);
        //    std::cout<<"test"<<std::endl;
        //    unsigned previous_node_local_index = (local_index-1+p_element->GetNumNodes())%p_element->GetNumNodes();
        //    c_vector<double, DIM> previous_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
        //    c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);
        //    c_vector<double, DIM> perimeter_gradient_at_node = previous_edge_gradient - next_edge_gradient;
            c_vector<double, DIM> F_Tens = -element_myosin_activity*Ga*perimeter*perimeter_gradient_at_node;

            if (mIfSetCellDataOfEachForceContributions)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Tens_x";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Tens[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Tens_y";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Tens[1]);
            }

            // Cell-cell adhesion
            double lambda_ca = 0.0;
            double lambda_ab = 0.0;
            if (shared_elements0.size() == 1)
                lambda_ca = mNagaiHondaCellBoundaryAdhesionEnergyParameter;
            else
                lambda_ca = mNagaiHondaCellCellAdhesionEnergyParameter;

            if (shared_elements.size() == 1)
                lambda_ab = mNagaiHondaCellBoundaryAdhesionEnergyParameter;
            else
                lambda_ab = mNagaiHondaCellCellAdhesionEnergyParameter;

            if (mIfConsiderFeedbackOfCellCellAdhesion)
            {
                lambda_ca *= pFace0->GetUnifiedCellCellAdhesionEnergyParameter();
                lambda_ab *= pFace->GetUnifiedCellCellAdhesionEnergyParameter();
            }
            c_vector<double,DIM> F_Adhe = -(lambda_ca*vec_ca/l_ca -lambda_ab*vec_ab/l_ab);

            if (mIfSetCellDataOfEachForceContributions)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Adhe_x";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Adhe[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_F_Adhe_y";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, F_Adhe[1]);
            }            

            // // Substrate adhesion
            // // consider periodicity: modify the node location!
            // if (dynamic_cast<MyXToroidal2dVertexMesh*>(&p_cell_population->rGetMesh()) != nullptr)
            // {
            //     bool triangle_straddles_left_right_boundary = false;

            //     if (fabs(vec_ca[0]) > 0.5*mWidth)
            //         triangle_straddles_left_right_boundary = true;
            //     if (fabs(vec_ab[0]) > 0.5*mWidth)
            //         triangle_straddles_left_right_boundary = true;
            //     if (triangle_straddles_left_right_boundary)
            //     {
            //         if (location_c[0] < mCenterOfWidth)
            //             location_c[0] += mWidth;
            //         if (location_a[0] < mCenterOfWidth)
            //             location_a[0] += mWidth;
            //         if (location_b[0] < mCenterOfWidth)
            //             location_b[0] += mWidth;
            //     }
            // }

            // double triangle_box_bottom = std::min(std::min(location_c[1],location_a[1]),location_b[1]);
            // double triangle_box_top = std::max(std::max(location_c[1],location_a[1]),location_b[1]);
            // double triangle_box_left = std::min(std::min(location_c[0],location_a[0]),location_b[0]);
            // double triangle_box_right = std::max(std::max(location_c[0],location_a[0]),location_b[0]);

            // bool full_strip_substrate_adhesion = false;
            // bool partial_strip_substrate_adhesion = false;
            // bool full_reservoir_substrate_adhesion = false;
            // bool partial_reservoir_substrate_adhesion = false;
            // bool no_substrate_adhesion = false;

            // double strip_start_y_location = this->mStripStartYLocation;
            // double strip_right = this->mStripStartXLocation + this->mStripWidth/2.0;
            // double strip_left = this->mStripStartXLocation - this->mStripWidth/2.0;
            // // consider periodicity. Important!!!
            // // if (fabs(triangle_box_left-triangle_box_right) > strip_distance/2)
            // //    partial_strip_substrate_adhesion = false;

            // if (triangle_box_bottom >= strip_start_y_location)
            // {
            //     if ((triangle_box_right <= strip_left) || (triangle_box_left >= strip_right))
            //         no_substrate_adhesion = true;
            //     else if ((triangle_box_left >= strip_left) && (triangle_box_right <= strip_right) )
            //         full_strip_substrate_adhesion = true;
            //     else
            //         partial_strip_substrate_adhesion = true;
            // }
            // else if ((triangle_box_top >= strip_start_y_location) && (triangle_box_bottom <= strip_start_y_location))
            // {
            //     if ((triangle_box_right <= strip_left) || (triangle_box_left >= strip_right))
            //         partial_reservoir_substrate_adhesion = true;
            //     else
            //     {   
            //         partial_strip_substrate_adhesion = true;
            //         partial_reservoir_substrate_adhesion = true;
            //     }
            // }
            // else if ((triangle_box_top <= strip_start_y_location) && (triangle_box_bottom >= 0.0))
            // {
            //     full_reservoir_substrate_adhesion = true;
            // }
            // else if ((triangle_box_top <= strip_start_y_location) && (triangle_box_top >= 0.0) && (triangle_box_bottom <= 0.0))
            // {
            //     partial_reservoir_substrate_adhesion = true;
            // }
            // else
            //     no_substrate_adhesion = true;

            // c_vector<double, DIM> strip_substrate_adhesion_area_gradient = zero_vector<double>(DIM);
            // c_vector<double, DIM> reservoir_substrate_adhesion_area_gradient = zero_vector<double>(DIM);
            // if (!no_substrate_adhesion)
            //     {if (full_strip_substrate_adhesion)
            //         strip_substrate_adhesion_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            //     else if (full_reservoir_substrate_adhesion)
            //         reservoir_substrate_adhesion_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);     
            //     else 
            //     {
            //         if (partial_strip_substrate_adhesion)
            //             strip_substrate_adhesion_area_gradient = this->GetStripSubstrateAdhesionAreaGradientOfElementAtNode(rCellPopulation, p_element, local_index);
            //         if (partial_reservoir_substrate_adhesion)
            //             reservoir_substrate_adhesion_area_gradient = this->GetReservoirSubstrateAdhesionAreaGradientOfElementAtNode(rCellPopulation, p_element, local_index);
            //     }
            // }
            
            // c_vector<double,DIM> F_Strip_Adh = -mStripSubstrateAdhesionParameter*strip_substrate_adhesion_area_gradient;
            // c_vector<double,DIM> F_Res_Adh = -mReservoirSubstrateAdhesionParameter*reservoir_substrate_adhesion_area_gradient;

            // if (mIfSetCellDataOfEachForceContributions)
            // {
            //     oss.str("");
            //     oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
            //             << "_F_Strip_Adh_x";
            //     name_item = oss.str();
            //     pCell->GetCellData()->SetItem(name_item, F_Strip_Adh[0]);
            //     oss.str("");
            //     oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
            //             << "_F_Strip_Adh_y";
            //     name_item = oss.str();
            //     pCell->GetCellData()->SetItem(name_item, F_Strip_Adh[1]);

            //     oss.str("");
            //     oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
            //             << "_F_Res_Adh_x";
            //     name_item = oss.str();
            //     pCell->GetCellData()->SetItem(name_item, F_Res_Adh[0]);
            //     oss.str("");
            //     oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
            //             << "_F_Res_Adh_y";
            //     name_item = oss.str();
            //     pCell->GetCellData()->SetItem(name_item, F_Res_Adh[1]);
            // } 

            // c_vector<double,DIM> F_Strip_Adh = -mStripSubstrateAdhesionParameter*strip_substrate_adhesion_area_gradient;
            // c_vector<double,DIM> F_Res_Adh = -mReservoirSubstrateAdhesionParameter*reservoir_substrate_adhesion_area_gradient;

            c_vector<double,DIM> Force = zero_vector<double>(DIM);
            Force[0] = F_Press[0] +F_Tens[0] +F_Adhe[0];// +F_Strip_Adh[0] +F_Res_Adh[0];
            Force[1] = F_Press[1] +F_Tens[1] +F_Adhe[1];// +F_Strip_Adh[1] +F_Res_Adh[1];

            // next: stress calculation
            c_vector<double,DIM> node_location_ralate_to_centroid =  p_mesh->GetVectorFromAtoB(centroid, location_a);
            shape_tensor_XX += 1.0/p_element->GetNumNodes()*node_location_ralate_to_centroid[0]*node_location_ralate_to_centroid[0];
            shape_tensor_YY += 1.0/p_element->GetNumNodes()*node_location_ralate_to_centroid[1]*node_location_ralate_to_centroid[1];
            shape_tensor_XY += 1.0/p_element->GetNumNodes()*node_location_ralate_to_centroid[0]*node_location_ralate_to_centroid[1];

            static_stress_XX += -1.0/area*node_location_ralate_to_centroid[0]*Force[0];
            static_stress_YY += -1.0/area*node_location_ralate_to_centroid[1]*Force[1];
            static_stress_XY += -1.0/area*node_location_ralate_to_centroid[0]*Force[1];
            static_stress_YX += -1.0/area*node_location_ralate_to_centroid[1]*Force[0];

            if (mWriteVertexVelocityAndForceToCellData)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_Fx";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, Force[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_Fy";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, Force[1]);

                // output velocity of vertices
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_vx";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, pNodeA->rGetAppliedForce()[0]);
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_vy";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, pNodeA->rGetAppliedForce()[1]);
            }

            // output face values of edges
            if (mIfConsiderFeedbackOfCellCellAdhesion)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_unified_cc";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, pFace->GetUnifiedCellCellAdhesionEnergyParameter());
            }

        } // end of iteration of vertices of the element, for calculation of force contributions and stress summation.
        
        // output myosin activity of this cell(element) to CellData.     
        if (true)
        {
            std::ostringstream oss;
            std::string name_item;
            oss.str("");
            oss << "elem_myo";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, p_element->GetElementMyosinActivity());
        }
        
        assert((shape_tensor_XX + shape_tensor_YY + shape_tensor_XY) < DOUBLE_UNSET );

    //     double shape_tensor_XX_old = pCell->GetCellData()->GetItem("shape_tensor_xx");
    //     double shape_tensor_YY_old = pCell->GetCellData()->GetItem("shape_tensor_yy");
    //     double shape_tensor_XY_old = pCell->GetCellData()->GetItem("shape_tensor_xy");

    //     double shape_tensor_XX_rate = (shape_tensor_XX - shape_tensor_XX_old)/dt;
    //     double shape_tensor_YY_rate = (shape_tensor_YY - shape_tensor_YY_old)/dt;
    //     double shape_tensor_XY_rate = (shape_tensor_XY - shape_tensor_XY_old)/dt;
    // //    double shape_tensor_YX_rate = shape_tensor_XY_rate;

    //     double dynamic_stress_XX = -1.0/2*shape_tensor_XX_rate;
    //     double dynamic_stress_YY = -1.0/2*shape_tensor_YY_rate;
    //     double dynamic_stress_XY = -1.0/2*shape_tensor_XY_rate;
    // //   double dynamic_stress_YX = -1.0/2*shape_tensor_YX_rate;

        double stress_XX = static_stress_XX; // + dynamic_stress_XX;
        double stress_YY = static_stress_YY; // + dynamic_stress_YY;
        double stress_XY = static_stress_XY; // + dynamic_stress_XY;

//        double R = sqrt( pow((stress_XX-stress_YY)/2, 2) + pow(stress_XY, 2) );
//       double stress_1 = (stress_XX + stress_YY)/2 + R;
//       double stress_2 = (stress_XX + stress_YY)/2 - R;

        double stress_delta = (stress_XX-stress_YY)*(stress_XX-stress_YY) + 4*stress_XY*stress_XY;
        double stress_1 = 0.5*(stress_XX + stress_YY + sqrt(stress_delta));
        double stress_2 = 0.5*(stress_XX + stress_YY - sqrt(stress_delta));
//        if ((stress_XX + stress_YY)>=0.0)
//            stress_lambda1 = 0.5*(stress_XX + stress_YY + sqrt(stress_delta));
//        else
//            stress_lambda1 = 0.5*(stress_XX + stress_YY - sqrt(stress_delta));

        c_vector<double, DIM> stress_eigen_vec_1 = zero_vector<double>(DIM);
        double principal_axis_of_stress = 0.0;
        if (stress_1==stress_XX)
            principal_axis_of_stress = 0.0;
        else if (stress_1==stress_YY)
            principal_axis_of_stress = M_PI/2;
        else
        {
            stress_eigen_vec_1[0] = 1.0;
            stress_eigen_vec_1[1] = (stress_1-stress_XX)/stress_XY;
            principal_axis_of_stress = atan(stress_eigen_vec_1[1]/stress_eigen_vec_1[0]);
        }

        double shape_delta = (shape_tensor_XX-shape_tensor_YY)*(shape_tensor_XX-shape_tensor_YY) +4*shape_tensor_XY*shape_tensor_XY;
        double shape_lambda1 = 0.5*(shape_tensor_XX+shape_tensor_YY+sqrt(shape_delta));
        double shape_lambda2 = 0.5*(shape_tensor_XX+shape_tensor_YY-sqrt(shape_delta));
        double aspect_ratio = sqrt(shape_lambda1/shape_lambda2);

        c_vector<double, DIM> shape_eigen_vec_1 = zero_vector<double>(DIM);
        double principal_axis_of_shape = 0.0;
        if (shape_lambda1==shape_tensor_XX)
            principal_axis_of_shape = 0.0;
        else if (shape_lambda1==shape_tensor_YY)
            principal_axis_of_shape = M_PI/2;
        else
        {
            shape_eigen_vec_1[0] = 1.0;
            shape_eigen_vec_1[1] = (shape_lambda1-shape_tensor_XX)/shape_tensor_XY;
            principal_axis_of_shape = atan(shape_eigen_vec_1[1]/shape_eigen_vec_1[0]);
        }

        // double R_shape = sqrt( pow((shape_XX-shape_YY)/2, 2) + pow(shape_XY, 2) );
        // double shape_1 = (shape_XX + shape_YY)/2 + R;
        // double shape_2 = (shape_XX + shape_YY)/2 - R;
        // double aspect_ratio = sqrt(shape_1/shape_2);
        // double principal_axis_angle = shape_XY>0.0 ? -0.5*acos( (shape_XX-0.5*(shape_XX+shape_YY))/R_shape ) : 0.5*acos( (shape_XX-0.5*(shape_XX+shape_YY))/R_shape );
        // pCell->GetCellData()->SetItem("AspectRatio", aspect_ratio);
        // pCell->GetCellData()->SetItem("PrincipalAxis", principal_axis_angle);

        // c_vector<double, DIM> eigen_vec_1 = zero_vector<double>(DIM);
        // eigen_vec_1[0] = shape_tensor_XY/sqrt( shape_tensor_XY*shape_tensor_XY+(shape_lambda1-shape_tensor_XX)*(shape_lambda1-shape_tensor_XX) );
        // eigen_vec_1[1] = (shape_lambda1-shape_tensor_XX)/sqrt( shape_tensor_XY*shape_tensor_XY+(shape_lambda1-shape_tensor_XX)*(shape_lambda1-shape_tensor_XX) );
        // double principal_axis_of_shape = abs(atan(eigen_vec_1[1]/eigen_vec_1[0]));
        // pCell->GetCellData()->SetItem("PrincipalAxisOfShape", principal_axis_of_shape);

        // double principal_axis_of_stress = 0.5*acos( (stress_XX-0.5*(stress_XX+stress_YY))/R );
        // pCell->GetCellData()->SetItem("PrincipalAxisOfStress", principal_axis_of_stress);

        pCell->GetCellData()->SetItem("StressXX", stress_XX);
        pCell->GetCellData()->SetItem("StressYY", stress_YY);
        pCell->GetCellData()->SetItem("StressXY", stress_XY);
        pCell->GetCellData()->SetItem("Stress1", stress_1);
        pCell->GetCellData()->SetItem("Stress2", stress_2);
        pCell->GetCellData()->SetItem("PrincipalAxisOfStress", principal_axis_of_stress);
        pCell->GetCellData()->SetItem("shape_tensor_xx", shape_tensor_XX);
        pCell->GetCellData()->SetItem("shape_tensor_yy", shape_tensor_YY);
        pCell->GetCellData()->SetItem("shape_tensor_xy", shape_tensor_XY);
        pCell->GetCellData()->SetItem("AspectRatio", aspect_ratio);
        pCell->GetCellData()->SetItem("PrincipalAxisOfShape", principal_axis_of_shape);

        // method for calculating principal stresses and principal axises:
        // A=[xx, xy; xy, yy], A-xI = [xx-x, xy; xy, yy-x], |A-xI|=(xx-x)(yy-x)-xy^2=x^2 -(xx+yy)x +(xx*yy-xy*xy)=0
        // x1,x2=1/2*(xx+yy +/-sqrt(delta)); delta=(xx^2+yy^2-2xx*yy+4xy*xy)
        // A-x1*I=[xx-x1, xy; xy, yy-x1]; A-x2*I=[xx-x2, xy; xy, yy-x2];
        // let x1>=x2, vec1=[1, x1-xx]/xy;
    }

}


template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateCellDataOfForcesFromNeighboringCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell)
{
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }

    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
    VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(elem_index);

    double pulling_force_y_from_upper_cells =0.0;
    double pulling_force_y_from_lower_cells =0.0;

    for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
    {
        //  --B
        //    |
        //  --A
        Node<DIM>* pNodeA = p_element->GetNode(local_index);
        Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
        c_vector<double,DIM> location_a = pNodeA->rGetLocation();
        c_vector<double,DIM> location_b = pNodeB->rGetLocation();
        c_vector<double,DIM> vec_ab =  p_mesh->GetVectorFromAtoB(location_a, location_b);

        std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
        std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
        // Find common elements
        std::set<unsigned> shared_elements;
        std::set_intersection(elements_containing_nodeA.begin(),
                            elements_containing_nodeA.end(),
                            elements_containing_nodeB.begin(),
                            elements_containing_nodeB.end(),
                            std::inserter(shared_elements, shared_elements.begin()));
        // Check that the nodes have a common edge
        assert(!shared_elements.empty());
        if (shared_elements.size() >= 3 || shared_elements.size()==0)
        {
            std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateCellDataOfForcesFromNeighboringCell";
            std::cout<< std::endl << "Get shared elements more than 2 or less than 1.";
            std::cout<< std::endl;
        }

        if (shared_elements.size()==1)
        {
            assert( *shared_elements.begin()==p_element->GetIndex() );
            // break;// such a tupid mistake!!! break will end the loop!!!
        }
        // pulling force along y direction from the neighboring upper cells and lower cells:
        else if (shared_elements.size()==2)
        {
            bool neighboring_cell_is_upward = false;
            bool neighboring_cell_is_downward = false;
            if (vec_ab[0]<=0.0) 
                neighboring_cell_is_upward = true;
            else
                neighboring_cell_is_downward = true;
            
            double Fy_from_this_neighbor_cell = 0.0;

            VertexElement<DIM, DIM>* p_neighbor_element = ( *shared_elements.begin()!=p_element->GetIndex() )?( p_mesh->GetElement(*shared_elements.begin()) ):( p_mesh->GetElement(*(++shared_elements.begin())) );
            CellPtr p_neighbor_cell = rCellPopulation.GetCellUsingLocationIndex(p_neighbor_element->GetIndex());
            unsigned local_index_of_lowest_vertex_neighbor_cell = p_neighbor_cell->GetCellData()->GetItem("LocalIndexOfLowestVertex");
            
            if (false)
            {
                if (elem_index==38)
                std::cout << std::endl << "timeStep=" << SimulationTime::Instance()->GetTimeStepsElapsed()
                        << ", elem_index=" << elem_index << ", local_index=" << local_index << ", neighbor_elem_index=" << p_neighbor_element
                        << std::endl;
            }

            std::string name_item;
            std::ostringstream oss;

            unsigned local_index_neighbor_elem = p_neighbor_element->GetNodeLocalIndex(pNodeA->GetIndex());
            oss.str("");
            oss << (local_index_neighbor_elem+p_neighbor_element->GetNumNodes()-local_index_of_lowest_vertex_neighbor_cell)%(p_neighbor_element->GetNumNodes()) 
                    << "_Fy";
            name_item = oss.str();
            Fy_from_this_neighbor_cell += p_neighbor_cell->GetCellData()->GetItem(name_item);

            if (false)// if (SimulationTime::Instance()->GetTimeStepsElapsed()==1411 || SimulationTime::Instance()->GetTimeStepsElapsed()==1412 || SimulationTime::Instance()->GetTimeStepsElapsed()==1413)
            {
                if (elem_index==38)
                std::cout << "local_indexA_neighbor_elem=" << local_index_neighbor_elem
                        << ", FyA=" << (p_neighbor_cell->GetCellData()->GetItem(name_item))
                        << std::endl;
            }

            local_index_neighbor_elem = p_neighbor_element->GetNodeLocalIndex(pNodeB->GetIndex());
            oss.str("");
            oss << (local_index_neighbor_elem+p_neighbor_element->GetNumNodes()-local_index_of_lowest_vertex_neighbor_cell)%(p_neighbor_element->GetNumNodes()) 
                    << "_Fy";
            name_item = oss.str();
            Fy_from_this_neighbor_cell += p_neighbor_cell->GetCellData()->GetItem(name_item);
            
            if (false)
            {
                if (elem_index==38)
                std::cout << "local_indexB_neighbor_elem=" << local_index_neighbor_elem
                        << ", FyB=" << (p_neighbor_cell->GetCellData()->GetItem(name_item))
                        << std::endl;
                std::cout << "Fy_from_this_neighbor_cell=" << Fy_from_this_neighbor_cell << std::endl;
            }

            if (neighboring_cell_is_upward)
                pulling_force_y_from_upper_cells += Fy_from_this_neighbor_cell;
            else
            {
                assert(neighboring_cell_is_downward);
                pulling_force_y_from_lower_cells += -Fy_from_this_neighbor_cell;
            }
        }
    }

    pCell->GetCellData()->SetItem("PullingForceFromUpperCellsY", pulling_force_y_from_upper_cells);
    pCell->GetCellData()->SetItem("PullingForceFromLowerCellsY", pulling_force_y_from_lower_cells);

}

/*
template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateUnifiedEdgeMyosinActivtyOfFace(MutableVertexMesh<DIM, DIM>* pMesh, unsigned faceIndex)
{
    double time_now = SimulationTime::Instance()->GetTime();
    double dt = SimulationTime::Instance()->GetTimeStep();
    double edge_length_at_rest = this->mEdgeLengthAtRest;
    double KL = this->mKLForFeedback;
    if ( time_now>mTimeForChangingFeedback && mChangedKLForFeedback!=mKLForFeedback)
        KL = mChangedKLForFeedback;
    double feedback_rate = this->mFeedbackRateForMyosinActivity;
    if ( time_now>mTimeForChangingFeedback && mChangedFeedbackRate!=mFeedbackRateForMyosinActivity)
        feedback_rate = mChangedFeedbackRate;
    double alpha = (1.0+KL)*feedback_rate;
    if ( time_now>mTimeForChangingFeedback && mChangedMyosinActivityBaseValue!=1.0 )
        alpha *= mChangedMyosinActivityBaseValue;
    double beta = 1.0*feedback_rate;
    double n = this->mHillPowerForMyosinActivity;

    VertexElement<DIM-1, DIM>* p_face = pMesh->GetFace(faceIndex);
    double edge_length = pMesh->GetDistanceBetweenNodes(p_face->GetNodeGlobalIndex(0), p_face->GetNodeGlobalIndex(1));
    double lambda = edge_length/ edge_length_at_rest;//(1.0*sqrt(M_PI/(6*sqrt(3)/4)));

    // // tmp
    // std::cout << "lambda=" << lambda << std::endl;
    double unified_edge_myosin_activty = p_face->GetUnifiedEdgeMyosinActivty();
    double changing_rate = alpha*pow(lambda,n)/(KL+pow(lambda,n))-beta*unified_edge_myosin_activty;
    if (mEMADontDecreaseWhenEdgeShrink && lambda<1.0)
    {
        changing_rate = 0.0;
    }
    if (mEMADontDecreaseBelowAThreshold && unified_edge_myosin_activty<mEMADontDecreaseBelowThisThreshold && changing_rate<0.0)
    {
        changing_rate = 0.0;
    }
    unified_edge_myosin_activty += dt*changing_rate;
    p_face->SetUnifiedEdgeMyosinActivty(unified_edge_myosin_activty);
    if (mIfOutputModifierInformation && random()%100000==0)
    {
        std::cout << std::endl<< "In FaceValueAndStressStateModifier::UpdateUnifiedEdgeMyosinActivityOfFace: the new UnifiedEdgeMyosinActivty=";
        std::cout << p_face->GetUnifiedEdgeMyosinActivty();
        std::cout << std::endl<< "feedback_rate=" << feedback_rate << " alpha=" << alpha << " beta=" <<beta << " KL=" << KL;
        std::cout << " n=" << n << " edge_length_at_rest=" << edge_length_at_rest << " the edge length=" << edge_length << " lambda=" << lambda;
        std::cout << " changing_rate=" << changing_rate;
    }
}
*/

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateMyosinActivtyOfElement(MutableVertexMesh<DIM, DIM>* pMesh, unsigned elem_index)
{
    double time_now = SimulationTime::Instance()->GetTime();
    double dt = SimulationTime::Instance()->GetTimeStep();

    double Km = this->mKmForMyosinFeedback;
    double feedback_rate = this->mFeedbackRateForMyosinActivity;
    double alpha = (1.0+Km)*feedback_rate;
    double beta = 1.0*feedback_rate;
    double n = this->mHillPowerForMyosinActivity;
    double FixedTargetPerimeter = this->mFixedTargetPerimeter;

    VertexElement<DIM, DIM>* pElement = pMesh->GetElement(elem_index);
    double element_perimeter = pMesh->GetSurfaceAreaOfElement(elem_index);
    double lambda = element_perimeter/FixedTargetPerimeter;

    double element_myosin_activity = pElement->GetElementMyosinActivity();

    if (time_now>mTimeForChangingFeedback)
    {
        Km = this->mChangedKmForMyosinFeedback;
        beta = this->mChangedFeedbackRate;
        alpha = (1+this->mChangedKmForMyosinFeedback)*this->mChangedFeedbackRate;
        element_myosin_activity = pElement->GetElementMyosinActivity()*(1.0/mChangedMyosinActivityBaseValue);
    }

    double changing_rate = alpha*pow(lambda,n)/(Km+pow(lambda,n))-beta*element_myosin_activity;
    element_myosin_activity += dt*changing_rate;

    element_myosin_activity = element_myosin_activity*mChangedMyosinActivityBaseValue;

    pElement->SetElementMyosinActivity(element_myosin_activity);

    if (mIfOutputModifierInformation && random()%100000==0)
    {
        std::cout << std::endl<< "In FaceValueAndStressStateModifier::UpdateMyosinActivityOfElement: the new ElementMyosinActivty=";
        std::cout << pElement->GetElementMyosinActivity();
        std::cout << std::endl<< "feedback_rate=" << feedback_rate << " alpha=" << alpha << " beta=" <<beta << " Km=" << Km;
        std::cout << " n=" << n << " fixed_target_periment=" << FixedTargetPerimeter << " the perimeter=" << element_perimeter << " lambda=" << lambda;
        std::cout << " changing_rate=" << changing_rate;
    }
}

// /*
// template<unsigned DIM>
// void FaceValueAndStressStateModifier<DIM>::UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(MutableVertexMesh<DIM, DIM>* pMesh, unsigned faceIndex)
// {
//     double dt = SimulationTime::Instance()->GetTimeStep();
//     double feedback_rate = this->mFeedbackRateForAdhesion;
//     double n = this->mHillPowerForAdhesion;
//     double edge_length_at_rest = this->mEdgeLengthAtRest;
//     double alpha = 2.0*feedback_rate;
//     double beta = 1.0*feedback_rate;
//     double KL = 1.0;

//     VertexElement<DIM-1, DIM>* p_face = pMesh->GetFace(faceIndex);
//     double edge_length = pMesh->GetDistanceBetweenNodes(p_face->GetNodeGlobalIndex(0), p_face->GetNodeGlobalIndex(1));
//     double lambda = edge_length/ edge_length_at_rest;//(1.0*sqrt(M_PI/(6*sqrt(3)/4)));
//     // // tmp
//     // std::cout << "lambda=" << lambda << std::endl;
//     double corresponding_lambda = 2.0-lambda;
//     double unified_cell_cell_adhesion_energy_parameter = p_face->GetUnifiedCellCellAdhesionEnergyParameter();
//     double changing_rate = alpha*pow(corresponding_lambda,n)/(KL+pow(corresponding_lambda,n))-beta*unified_cell_cell_adhesion_energy_parameter;
//     if (mCCADontDecreaseWhenEdgeExpand && lambda>1.0)
//     {
//         changing_rate = 0.0;
//     }
//     if (mCCAIncreasingHasAThresholdOfEdgeLength && lambda>mCCAIncreasingThresholdOfEdgeLengthPercentage)
//     {
//         changing_rate = 0.0;
//     }
//     unified_cell_cell_adhesion_energy_parameter += changing_rate*dt;
//     p_face->SetUnifiedCellCellAdhesionEnergyParameter(unified_cell_cell_adhesion_energy_parameter);
// }
// */

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned faceIndex)
{
    MutableVertexMesh<DIM, DIM>* pMesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());
    double dt = SimulationTime::Instance()->GetTimeStep();

    double feedback_rate = this->mFeedbackRateForAdhesion;
    assert(feedback_rate!=0);
    double q = this->mHillPowerForAdhesion;
    double Ks = this->mKsForAdhesionFeedback;
    double pho = (1+Ks)*feedback_rate;
    double ksi = 1.0*feedback_rate;
    double reference_stress = mReferenceStress;

    VertexElement<DIM-1, DIM>* pFace = pMesh->GetFace(faceIndex);
    Node<DIM>* pNodeA = pFace->GetNode(0);
    Node<DIM>* pNodeB = pFace->GetNode(1);

    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));
    // Check that the nodes have a common edge
    // if (shared_elements.empty())
    // {
    //     c_vector<double, DIM> locationA = pNodeA->rGetLocation();
    //     c_vector<double, DIM> locationB = pNodeB->rGetLocation();
    //     std::cout<<"NodeA="<<locationA[0]<<","<<locationB[1]<<std::endl;
    //     std::cout<<"NodeB="<<locationB[0]<<","<<locationB[1]<<std::endl;
    // }

    assert(!shared_elements.empty());
    if (shared_elements.size() == 1)
    {
        double unified_cell_cell_adhesion_energy_parameter = 0.0;
        pFace->SetUnifiedCellCellAdhesionEnergyParameter(unified_cell_cell_adhesion_energy_parameter);
    }
    else if (shared_elements.size() == 2)
    {
        c_vector<double, DIM> vector_of_face = zero_vector<double>(DIM);
        vector_of_face = pMesh->GetVectorFromAtoB(pNodeA->rGetLocation(), pNodeB->rGetLocation());
        c_vector<double, DIM> vector_normal_to_face = zero_vector<double>(DIM);
        vector_normal_to_face[0] = vector_of_face[1];
        vector_normal_to_face[1] = -vector_of_face[0];
        double length_of_vector_normal_to_face = norm_2(vector_normal_to_face);

        VertexElement<DIM, DIM>* pElement = pMesh->GetElement(*shared_elements.begin());
        CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(pElement->GetIndex());
        double stress1 = pCell->GetCellData()->GetItem("Stress1");
        double stress2 = pCell->GetCellData()->GetItem("Stress2");
        double principal_axis_of_stress = pCell->GetCellData()->GetItem("PrincipalAxisOfStress");
        c_vector<double, DIM> vector_of_principal_axis = zero_vector<double>(DIM);
        vector_of_principal_axis[0] = cos(principal_axis_of_stress);
        vector_of_principal_axis[1] = sin(principal_axis_of_stress);
        double length_of_vector_of_principal_axis = norm_2(vector_of_principal_axis);
        double dot_product = vector_of_principal_axis[0]*vector_normal_to_face[0]+vector_of_principal_axis[1]*vector_normal_to_face[1];
        double angle = acos(dot_product/(length_of_vector_normal_to_face*length_of_vector_of_principal_axis));
        double stress_normal_to_face_in_first_element = (stress1+stress2)/2 + (stress1-stress2)/2*cos(2*angle);

//        std::cout<<"stress1="<<stress1<<std::endl;
//        std::cout<<"stress2="<<stress2<<std::endl;
//        std::cout<<"angle="<<angle<<std::endl;
//        std::cout<<"dot_product="<<dot_product<<std::endl;
//        std::cout<<"normal_to_face="<<length_of_vector_normal_to_face<<std::endl;
//        std::cout<<"principal_axis="<<length_of_vector_of_principal_axis<<std::endl;

        // pElement = pMesh->GetElement(*shared_elements.end());
        // // testing .end().
        // std::cout << std::endl << "EleIndex by *set.end(): " << *shared_elements.end() << ". EleIndex by *++set.begin(): " << *(++shared_elements.begin()) << std::endl;
        pElement = pMesh->GetElement(*(++shared_elements.begin()));
        pCell = rCellPopulation.GetCellUsingLocationIndex(pElement->GetIndex());
        stress1 = pCell->GetCellData()->GetItem("Stress1");
        stress2 = pCell->GetCellData()->GetItem("Stress2");
        principal_axis_of_stress = pCell->GetCellData()->GetItem("PrincipalAxisOfStress");
        vector_of_principal_axis[0] = cos(principal_axis_of_stress);
        vector_of_principal_axis[1] = sin(principal_axis_of_stress);
        length_of_vector_of_principal_axis = norm_2(vector_of_principal_axis);
        dot_product = vector_of_principal_axis[0]*vector_normal_to_face[0]+vector_of_principal_axis[1]*vector_normal_to_face[1];
        angle = acos(dot_product/(length_of_vector_normal_to_face*length_of_vector_of_principal_axis));
        double stress_normal_to_face_in_second_element = (stress1+stress2)/2 + (stress1-stress2)/2*cos(2*angle);

        double sigma = 1.0/2*(fabs(stress_normal_to_face_in_first_element)+fabs(stress_normal_to_face_in_second_element))/reference_stress;
        double unified_cell_cell_adhesion_energy_parameter = pFace->GetUnifiedCellCellAdhesionEnergyParameter();
        double changing_rate = pho*pow(sigma,q)/(Ks+pow(sigma,q))-ksi*unified_cell_cell_adhesion_energy_parameter;

        if (mCellCellAdhesionDontDecrease && changing_rate<0.0 && unified_cell_cell_adhesion_energy_parameter<=1.0)
            unified_cell_cell_adhesion_energy_parameter += 0.0;
        else
            unified_cell_cell_adhesion_energy_parameter += changing_rate*dt;
        pFace->SetUnifiedCellCellAdhesionEnergyParameter(unified_cell_cell_adhesion_energy_parameter);

    //    std::cout<<"initial_sigma="<<reference_stress<<std::endl;
    //    std::cout<<"stress_in_element1="<<fabs(stress_normal_to_face_in_first_element)<<std::endl;
    //    std::cout<<"stress_in_element2="<<fabs(stress_normal_to_face_in_second_element)<<std::endl;
    //    std::cout<<"sigma="<<sigma<<std::endl;
    }
    else
    {
        std::cout<< "Get error in FaceValueAndStressStateModifier::UpdateUnifiedCellCellAdhesionEnergyParameterOfFace"<< std::endl ;
        std::cout<< "Get shared elements more than 2"<< std::endl;
    }

}


template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateCellAreas(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        UpdateCellAreaOfCell(rCellPopulation, *cell_iter);
    }
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateCellAreaOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell)
{
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }

    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    double cell_area = p_mesh->GetVolumeOfElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );
    pCell->GetCellData()->SetItem("cell area", cell_area);

    double cell_perimeter = p_mesh->GetSurfaceAreaOfElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );
    pCell->GetCellData()->SetItem("cell perimeter", cell_perimeter);

    double centroid_x = p_mesh->GetCentroidOfElement( rCellPopulation.GetLocationIndexUsingCell(pCell) )(0);
    double centroid_y = p_mesh->GetCentroidOfElement( rCellPopulation.GetLocationIndexUsingCell(pCell) )(1);
    pCell->GetCellData()->SetItem("centroid x", centroid_x);
    pCell->GetCellData()->SetItem("centroid y", centroid_y);

    VertexElement<DIM, DIM>* pLeadingElement = nullptr;
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        VertexElement<DIM, DIM>* pElement = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(*cell_iter) );
        if (pElement->GetIsLeadingCell())
        {
            pLeadingElement = pElement;
            break;
        }
    }
    if (pLeadingElement!=nullptr)
    {
        double distance_to_leading_cell = norm_2( p_mesh->GetCentroidOfElement( rCellPopulation.GetLocationIndexUsingCell(pCell) )-p_mesh->GetCentroidOfElement(pLeadingElement->GetIndex()) );
        pCell->GetCellData()->SetItem("Distance to leading cell", distance_to_leading_cell);
    }

}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::SetupSolveForCellDivision(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        (*cell_iter)->SetUseMyDivisionRuleAlongWithModifier(true);
    }
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateForCellDivision(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
        SimulationTime* p_time = SimulationTime::Instance();
        bool at_division_time = (p_time->GetTime() -mDivisionTime/2.0 - mDivisionTime*floor((p_time->GetTime()-mDivisionTime/2.0)/mDivisionTime)) < p_time->GetTimeStep();
        if (at_division_time)
        {
            unsigned division_element_index = (unsigned)floor(RandomNumberGenerator::Instance()->ranf()*( (double)(rCellPopulation.GetNumAllCells())-0.01 ));
            assert(division_element_index < rCellPopulation.GetNumAllCells());
            std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
            for (unsigned i=0; i< division_element_index; i++)
                cell_iter++;
            (*cell_iter)->SetCanDivide(true);
        }
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateGroupNumbers(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        UpdateGroupNumberOfCell(rCellPopulation, *cell_iter);
    }
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateGroupNumberOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell)
{
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }

    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    unsigned group_number = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(pCell) )->GetGroupNumber();
    pCell->GetCellData()->SetItem("group number", group_number);
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::SetupLeaderCellAtTheEndOfEquilibrium(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    SimulationTime* p_time = SimulationTime::Instance();
    if (! (mIfEquilibrateForAWhile && p_time->GetTime()<mEndTimeForEquilibrium && p_time->GetTime()>(mEndTimeForEquilibrium-2*p_time->GetTimeStep()) ) )
        return;

    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }
    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    assert(!mMultipleLeadingCells);
    unsigned LeaderCellElementIndex = 10000000;
    double LeaderCellDisToStrip = 100000.0;
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        CellPtr pCell = *cell_iter;
        VertexElement<DIM, DIM>* pElement = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );
        unsigned ele_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
        for (unsigned index=0; index< pElement->GetNumNodes(); index++)
        {
            if (pElement->GetNode(index)->IsBoundaryNode())
            {
                // is_boundary_cell = true;
                double dis = fabs(p_mesh->GetCentroidOfElement(ele_index)[0] - mStripStartXLocation);
                if (p_mesh->GetCentroidOfElement(ele_index)[1]>mStripStartYLocation/2.0 && dis<LeaderCellDisToStrip)
                {
                    LeaderCellElementIndex = ele_index;
                    LeaderCellDisToStrip = dis;
                }
                break;
            }
        }
    }
    assert(LeaderCellElementIndex!=10000000);

    // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        CellPtr pCell = *cell_iter;
        unsigned ele_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
        VertexElement<DIM, DIM>* pElement = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );

        if (ele_index == LeaderCellElementIndex)
        {
            pElement->SetIsLeadingCell(true);
            pElement->SetIsLeadingCellTop(true);
            pElement->SetIsLeadingCellBottom(false);
            pElement->SetIsJustReAttached(false);
            pElement->SetLamellipodiumStrength(1.0);
            pCell->GetCellData()->SetItem("is_leading_cell", 1);
            pCell->GetCellData()->SetItem("is_leading_cell_top", 1);
            pCell->GetCellData()->SetItem("is_leading_cell_bottom", 0);
            pCell->GetCellData()->SetItem("is_just_reattached", 0);
            pCell->GetCellData()->SetItem("lamellipodium_strength", 1.0);
        }
        else
        {
            pElement->SetIsLeadingCell(false);
            pElement->SetIsLeadingCellTop(false);
            pElement->SetIsLeadingCellBottom(false);
            pElement->SetIsJustReAttached(false);
            pElement->SetLamellipodiumStrength(0.0);
            pCell->GetCellData()->SetItem("is_leading_cell", 0);
            pCell->GetCellData()->SetItem("is_leading_cell_top", 0);
            pCell->GetCellData()->SetItem("is_leading_cell_bottom", 0);
            pCell->GetCellData()->SetItem("is_just_reattached", 0);
            pCell->GetCellData()->SetItem("lamellipodium_strength", 0.0);
        }
    }
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::SetupSolveForLamellipodiumInfoOfCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }
    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        CellPtr pCell = *cell_iter;
        VertexElement<DIM, DIM>* pElement = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );

        unsigned ele_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
        if (    ( (p_mesh->GetCentroidOfElement(ele_index)[1] > mStripStartYLocation*5.0/6.0)&&(!mMultipleLeadingCells)&&(fabs(p_mesh->GetCentroidOfElement(ele_index)[0] - mStripStartXLocation)< 1e-1) ) 
             || ( (p_mesh->GetCentroidOfElement(ele_index)[1] > mStripStartYLocation*5.0/6.0)&&mMultipleLeadingCells&&(fabs(p_mesh->GetCentroidOfElement(ele_index)[0] - mStripStartXLocation)< mStripWidth/2.0) ) )
        {
            bool is_leading_cell = false;
            for (unsigned index=0; index< pElement->GetNumNodes(); index++)
            {
                if (pElement->GetNode(index)->IsBoundaryNode())
                {
                    pElement->SetIsLeadingCell(true);
                    pElement->SetIsLeadingCellTop(true);
                    pElement->SetIsLeadingCellBottom(false);
                    pElement->SetIsJustReAttached(false);
                    pElement->SetLamellipodiumStrength(1.0);
                    pCell->GetCellData()->SetItem("is_leading_cell", 1);
                    pCell->GetCellData()->SetItem("is_leading_cell_top", 1);
                    pCell->GetCellData()->SetItem("is_leading_cell_bottom", 0);
                    pCell->GetCellData()->SetItem("is_just_reattached", 0);
                    pCell->GetCellData()->SetItem("lamellipodium_strength", 1.0);
                    is_leading_cell = true;
                    break;
                }
            }
            if (!is_leading_cell)
            {
                pElement->SetIsLeadingCell(false);
                pElement->SetIsLeadingCellTop(false);
                pElement->SetIsLeadingCellBottom(false);
                pElement->SetIsJustReAttached(false);
                pElement->SetLamellipodiumStrength(0.0);
                pCell->GetCellData()->SetItem("is_leading_cell", 0);
                pCell->GetCellData()->SetItem("is_leading_cell_top", 0);
                pCell->GetCellData()->SetItem("is_leading_cell_bottom", 0);
                pCell->GetCellData()->SetItem("is_just_reattached", 0);
                pCell->GetCellData()->SetItem("lamellipodium_strength", 0.0);
            }
        }
        else
        {
            pElement->SetIsLeadingCell(false);
            pElement->SetIsLeadingCellTop(false);
            pElement->SetIsLeadingCellBottom(false);
            pElement->SetIsJustReAttached(false);
            pElement->SetLamellipodiumStrength(0.0);
            pCell->GetCellData()->SetItem("is_leading_cell", 0);
            pCell->GetCellData()->SetItem("is_leading_cell_top", 0);
            pCell->GetCellData()->SetItem("is_leading_cell_bottom", 0);
            pCell->GetCellData()->SetItem("is_just_reattached", 0);
            pCell->GetCellData()->SetItem("lamellipodium_strength", 0.0);
        }
    }

}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateLamellipodiumInfoOfCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh()) == nullptr)
    {
        EXCEPTION("MutableVertexMesh should to be used in the FaceValueAndStressStateModifier");
    }
    MutableVertexMesh<DIM, DIM>* p_mesh = static_cast<MutableVertexMesh<DIM, DIM>*>(& rCellPopulation.rGetMesh());

    double dt = SimulationTime::Instance()->GetTimeStep();

    // Get the upper-boundary LeaderElementVector. And set value.
    if (mMultipleLeadingCells)
    {
        std::vector< VertexElement<DIM, DIM>* > ElementVector;
        std::vector< VertexElement<DIM, DIM>* > LeaderElementVector;

        // Get the upper-boundary ElementVector.
        for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
            cell_iter != rCellPopulation.rGetCells().end();
            ++cell_iter)
        {

            CellPtr pCell = *cell_iter;
            VertexElement<DIM, DIM>* pElement = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );
            unsigned ele_index = rCellPopulation.GetLocationIndexUsingCell(pCell);

            if ( (p_mesh->GetCentroidOfElement(ele_index)[1] > 5.0/6.0*mStripStartYLocation) && (fabs(p_mesh->GetCentroidOfElement(ele_index)[0] - mStripStartXLocation)<= mStripWidth/2.0) )
            {
                bool is_leading_cell = false;
                for (unsigned index=0; index< pElement->GetNumNodes(); index++)
                {
                    if (pElement->GetNode(index)->IsBoundaryNode())
                    {
                        pElement->SetIsLeadingCell(true);
                        pElement->SetIsLeadingCellTop(true);
                        is_leading_cell = true;
                        break;
                    }
                }
                if (is_leading_cell)
                    ElementVector.push_back(pElement);
            }
        }

        LeaderElementVector = ElementVector;

        double LeftHigh = 0.0;
        double RightHigh = 0.0;
        for (typename std::vector<VertexElement<DIM, DIM>*>::iterator itt=ElementVector.begin(); itt!= ElementVector.end();itt++)
        {
            VertexElement<DIM, DIM>* current_elementtt = *itt;
            unsigned current_ele_indextt = current_elementtt->GetIndex();
            double x_location = p_mesh->GetCentroidOfElement(current_ele_indextt)[0];
            double y_location = p_mesh->GetCentroidOfElement(current_ele_indextt)[1];
            if ( (x_location- mStripStartXLocation) < (-mStripWidth/2.0+0.95*2) )
            {
                if (y_location>LeftHigh)
                    LeftHigh = y_location;
            }
            if ( (x_location- mStripStartXLocation) > (mStripWidth/2.0-0.95*2) )
            {
                if (y_location>RightHigh)
                    RightHigh = y_location;
            }
        }
        
        for (typename std::vector<VertexElement<DIM, DIM>*>::iterator itt=ElementVector.begin(); itt!= ElementVector.end();itt++)
        {
            bool is_leader =true;
            VertexElement<DIM, DIM>* current_elementtt = *itt;
            unsigned current_ele_indextt = current_elementtt->GetIndex();
            double x_location = p_mesh->GetCentroidOfElement(current_ele_indextt)[0];
            double y_location = p_mesh->GetCentroidOfElement(current_ele_indextt)[1];
            if ( (x_location- mStripStartXLocation) < (-mStripWidth/2.0+0.95*2) )
            {
                if (y_location<(LeftHigh-0.1))
                    is_leader = false;
            }
            if ( (x_location- mStripStartXLocation) > (mStripWidth/2.0-0.95*2) )
            {
                if (y_location<(RightHigh-0.1))
                    is_leader = false;
            }
            if (!is_leader)
            {
                for (typename std::vector<VertexElement<DIM, DIM>*>::iterator it=LeaderElementVector.begin(); it!=LeaderElementVector.end(); it++)
                {
                    if ( current_elementtt == *it )
                    {
                        LeaderElementVector.erase(it);
                        break;
                    }
                }
            }
        }// Get LeaderElementVec finished.



        // // select N cells from them.
        // for (unsigned i=0; i<mLeadingCellNumber; i++)
        // {
        //     for (typename std::vector<VertexElement<DIM, DIM>*>::iterator it=ElementVector.begin(); it!= ElementVector.end();it++)
        //     {
        //         bool is_current_nearest = true;
        //         VertexElement<DIM, DIM>* current_element = *it;
        //         unsigned current_ele_index = current_element->GetIndex();
        //         double dis = fabs(p_mesh->GetCentroidOfElement(current_ele_index)[0]);
        //         for (typename std::vector<VertexElement<DIM, DIM>*>::iterator itt=ElementVector.begin(); itt!= ElementVector.end();itt++)
        //         {
        //             VertexElement<DIM, DIM>* current_elementtt = *itt;
        //             unsigned current_ele_indextt = current_elementtt->GetIndex();
        //             double distt = fabs(p_mesh->GetCentroidOfElement(current_ele_indextt)[0]);
        //             if (distt<dis)
        //             {
        //                 is_current_nearest = false;
        //                 break;
        //             }
                    
        //         }
        //         if (is_current_nearest)
        //         {
        //             LeaderElementVector.push_back(current_element);
        //             ElementVector.erase(it);
        //             break;
        //         }
        //     }
        // } // Get LeaderElementVec finished.

        // write leader info to element.
        for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
            cell_iter != rCellPopulation.rGetCells().end();
            ++cell_iter)
        {
            CellPtr pCell = *cell_iter;
            VertexElement<DIM, DIM>* pElement = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );

            bool is_leading_cell =false;
            for (typename std::vector<VertexElement<DIM, DIM>*>::iterator it=LeaderElementVector.begin(); it!=LeaderElementVector.end(); it++)
            {
                if ( pElement == *it )
                {
                    is_leading_cell = true;
                    break;
                }
            }
            if (is_leading_cell)
            {
                pElement->SetIsLeadingCell(true);
                pElement->SetIsLeadingCellTop(true);
            }
            else
            {
                pElement->SetIsLeadingCell(false);
                pElement->SetIsLeadingCellTop(false);
            }

        }

    }

    // write it to CellData.
    // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
    for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
        cell_iter != rCellPopulation.rGetCells().end();
        ++cell_iter)
    {
        CellPtr pCell = *cell_iter;
        VertexElement<DIM, DIM>* pElement = p_mesh->GetElement( rCellPopulation.GetLocationIndexUsingCell(pCell) );
        
        pCell->GetCellData()->SetItem("is_leading_cell", pElement->GetIsLeadingCell());
        pCell->GetCellData()->SetItem("is_leading_cell_top", pElement->GetIsLeadingCellTop());
        pCell->GetCellData()->SetItem("is_leading_cell_bottom", pElement->GetIsLeadingCellBottom());
        pCell->GetCellData()->SetItem("is_just_reattached", pElement->GetIsJustReAttached());
        pCell->GetCellData()->SetItem("lamellipodium_strength", pElement->GetLamellipodiumStrength());

        if (pElement->GetLamellipodiumStrength()!=0.0 && !pElement->GetIsLeadingCell())
        {
            pElement->SetIsJustReAttached(true);
            pCell->GetCellData()->SetItem("is_just_reattached", true);
        }
        
        if (pElement->GetIsLeadingCell()) // leading cell
        {
            double new_lamellipodium_strength = std::min(1.0, (pElement->GetLamellipodiumStrength()+mLamellipodiumMaturationRate*dt));
            pElement->SetLamellipodiumStrength(new_lamellipodium_strength);
            pCell->GetCellData()->SetItem("lamellipodium_strength", new_lamellipodium_strength);
        }
        else if (pElement->GetIsJustReAttached()) // intermediate cell but still has lamellipodium!
        {
            double new_lamellipodium_strength = std::max(0.0, (pElement->GetLamellipodiumStrength()-mLamellipodiumDestructionRate*dt));
            pElement->SetLamellipodiumStrength(new_lamellipodium_strength);
            pCell->GetCellData()->SetItem("lamellipodium_strength", new_lamellipodium_strength);
            if (new_lamellipodium_strength==0.0)
            {
                pElement->SetIsJustReAttached(false);
                pCell->GetCellData()->SetItem("is_just_reattached", false);
            }
        }
        else
            assert(pElement->GetLamellipodiumStrength()==0.0);
    }
}

template <unsigned DIM>
c_vector<double, DIM> FaceValueAndStressStateModifier<DIM>::GetStripSubstrateAdhesionAreaGradientOfElementAtNode(AbstractCellPopulation<DIM>& rCellPopulation, VertexElement<DIM, DIM>* pElement, unsigned localIndex)
{
    assert(DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

//    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    c_vector<double, DIM> strip_substrate_adhesion_area_gradient = zero_vector<double>(DIM);

    // Parameters
    double strip_width = this->mStripWidth;
    double strip_distance = this->mStripDistance;
    double strip_start_x_location = this->mStripStartXLocation;
    double strip_start_y_location = this->mStripStartYLocation;

    double small_change = this->mSmallChangeForAreaCalculation;
    double preferred_sample_dis = small_change/5.0;

//    c_vector<double, DIM> centroid = p_cell_population->rGetMesh().GetCentroidOfElement(elem_index);
    ///\todo This should probably be returning the nearest node
    c_vector<double, DIM> centroid = zero_vector<double>(DIM);
    unsigned num_nodes_elem = pElement->GetNumNodes();
    for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
    {
        // Find location of current node and add it to the centroid
        centroid += pElement->GetNodeLocation(local_index);
    }
    centroid /= num_nodes_elem;

    int strip_num = round((centroid[0] - strip_start_x_location)/strip_distance);

    double strip_location = strip_start_x_location + strip_distance*strip_num;
    double strip_left = strip_location - strip_width/2;
    double strip_right = strip_location + strip_width/2;

    double sample_area_bottom = 0.0;
    double sample_area_top = 0.0;
    double sample_area_left = 0.0;
    double sample_area_right = 0.0;

    Node<DIM>* p_this_node = pElement->GetNode(localIndex);
    unsigned previous_node_local_index = (num_nodes_elem+localIndex-1)%num_nodes_elem;
    Node<DIM>* p_previous_node = pElement->GetNode(previous_node_local_index);
    unsigned next_node_local_index = (localIndex+1)%num_nodes_elem;
    Node<DIM>* p_next_node = pElement->GetNode(next_node_local_index);

    c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
    c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
    c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
    
    double expanded_triangle_box_bottom = -small_change+std::min(std::min(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_left = -small_change+std::min(std::min(previous_node_location[0],this_node_location[0]),next_node_location[0]);
    double expanded_triangle_box_right = small_change+std::max(std::max(previous_node_location[0],this_node_location[0]),next_node_location[0]);
    
    bool has_substrate_adhesion_area = true;
    // consider periodicity. Important!!!
    if (fabs(expanded_triangle_box_left-expanded_triangle_box_right) > strip_distance/2)
        has_substrate_adhesion_area = false;

    if (expanded_triangle_box_top<strip_start_y_location)
        has_substrate_adhesion_area = false;
    else
    {
        sample_area_top = expanded_triangle_box_top;
        sample_area_bottom = expanded_triangle_box_bottom>strip_start_y_location ? expanded_triangle_box_bottom : strip_start_y_location;
    }

    if ((expanded_triangle_box_left <= strip_left)&&(expanded_triangle_box_right >= strip_left)&&(expanded_triangle_box_right <= strip_right))
    {   
        sample_area_left = strip_left;
        sample_area_right = expanded_triangle_box_right;
    }
    else if ((expanded_triangle_box_left <= strip_left)&&(expanded_triangle_box_right >= strip_right))
    {   
        sample_area_left = strip_left;
        sample_area_right = strip_right;
    }
    else if ((expanded_triangle_box_left >= strip_left)&&(expanded_triangle_box_right <= strip_right))
    {
        sample_area_left = expanded_triangle_box_left;
        sample_area_right = expanded_triangle_box_right;
    }
    else if ((expanded_triangle_box_left >= strip_left)&&(expanded_triangle_box_left <= strip_right)&&(expanded_triangle_box_right >= strip_right))
    {   
        sample_area_left = expanded_triangle_box_left;
        sample_area_right = strip_right;
    }
    else
        has_substrate_adhesion_area = false;

    unsigned num_across = 0;
    unsigned num_up = 0;
    unsigned sample_num = 0;
    double   sample_area = 0.0;

    if (has_substrate_adhesion_area)
    {
        num_across = round((sample_area_right - sample_area_left) / preferred_sample_dis);
        num_up = round((sample_area_top - sample_area_bottom) / preferred_sample_dis);
        sample_num = num_across * num_up;
        sample_area = (sample_area_top - sample_area_bottom)*(sample_area_right - sample_area_left);
    }
    
    // Calculate strip_substrate_adhesion_area_gradient!
    if (has_substrate_adhesion_area && (num_across>0) && (num_up>0))
    {
        c_vector<c_vector<double, DIM>, 3> points;
        points[0] = previous_node_location;
        points[1] = this_node_location;
        points[2] = next_node_location;
        c_vector<double, DIM> vec1 = this_node_location-previous_node_location;
        c_vector<double, DIM> vec2 = next_node_location-this_node_location;

        // Calculate initial adhesive area
        double adhesive_sample_num = 0.0;
        for (unsigned i = 0; i <sample_num; i++)
        {
            double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
            double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

            c_vector<double, DIM> point;
            point[0] = x_coord;
            point[1] = y_coord;
            c_vector<bool, 3> point_at_left_of_vector;
            for (unsigned j = 0; j < 3; j++)
            {
                c_vector<double, DIM> vec1 = point - points[j];
                c_vector<double, DIM> vec2 = points[(j+1)%3] - points[j];

                if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                    point_at_left_of_vector[j] = false;
                else
                    point_at_left_of_vector[j] = true;
            }
            if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                adhesive_sample_num += 1.0;
            else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                adhesive_sample_num += -1.0;
        }
        double substrate_adhesion_area = adhesive_sample_num/double(sample_num) * sample_area;

        // Calculate adhesive area *changed* with small displacement of the node along the x or y axis
        for (unsigned j = 0; j<2; j++)
        {
            c_vector<double, DIM> unit_vector_in_small_change_direction;
            unit_vector_in_small_change_direction[0] = cos(j*M_PI/2);
            unit_vector_in_small_change_direction[1] = sin(j*M_PI/2);
            // my new changes for cosistent movement of nodes at edges of the strip.
            if (this->mConsiderConsistencyForSSA)
            {
                if (points[1][0]-strip_right>-small_change && points[1][0]-strip_right<small_change)
                    unit_vector_in_small_change_direction[0] *= -1.0;
            }

            c_vector<c_vector<double, DIM>, 3> new_points = points;
            new_points[1] += small_change*unit_vector_in_small_change_direction;

            double adhesive_sample_num = 0.0;

            for (unsigned i = 0; i<sample_num; i++)
            {
                double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                c_vector<double, DIM> point;
                point[0] = x_coord;
                point[1] = y_coord;
                c_vector<bool, 3> point_at_left_of_vector;
                for (unsigned j = 0; j < 3; j++)
                {
                    c_vector<double, DIM> vec1 = point - new_points[j];
                    c_vector<double, DIM> vec2 = new_points[(j+1)%3] - new_points[j];

                    if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                        point_at_left_of_vector[j] = false;
                    else
                        point_at_left_of_vector[j] = true;
                }
                if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                    adhesive_sample_num += 1.0;
                else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                    adhesive_sample_num += -1.0;
            }
            double substrate_adhesion_area_new = adhesive_sample_num/double(sample_num) * sample_area;
            
            strip_substrate_adhesion_area_gradient += (substrate_adhesion_area_new - substrate_adhesion_area)/small_change*unit_vector_in_small_change_direction; 
        }
    }
    return strip_substrate_adhesion_area_gradient;
}

template <unsigned DIM>
c_vector<double, DIM> FaceValueAndStressStateModifier<DIM>::GetReservoirSubstrateAdhesionAreaGradientOfElementAtNode(AbstractCellPopulation<DIM>& rCellPopulation, VertexElement<DIM,DIM>* pElement, unsigned localIndex)
{
    assert(DIM == 2); // LCOV_EXCL_LINE - code will be removed at compile time

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    c_vector<double, DIM> reservoir_substrate_adhesion_area_gradient = zero_vector<double>(DIM);

    // Parameters
    double strip_start_y_location = this->mStripStartYLocation;
    double small_change = this->mSmallChangeForAreaCalculation;
    double preferred_sample_dis = small_change/5.0;

    c_vector<double, DIM> centroid = zero_vector<double>(DIM);
    unsigned num_nodes_elem = pElement->GetNumNodes();
    for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
    {
        // Find location of current node and add it to the centroid
        centroid += pElement->GetNodeLocation(local_index);
    }
    centroid /= num_nodes_elem;

    double sample_area_bottom = 0.0;
    double sample_area_top = 0.0;
    double sample_area_left = 0.0;
    double sample_area_right = 0.0;

    Node<DIM>* p_this_node = pElement->GetNode(localIndex);
    unsigned previous_node_local_index = (num_nodes_elem+localIndex-1)%num_nodes_elem;
    Node<DIM>* p_previous_node = pElement->GetNode(previous_node_local_index);
    unsigned next_node_local_index = (localIndex+1)%num_nodes_elem;
    Node<DIM>* p_next_node = pElement->GetNode(next_node_local_index);

    c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
    c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
    c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
    
    // consider periodicity: modify the node location!
    if (dynamic_cast<MyXToroidal2dVertexMesh*>(&p_cell_population->rGetMesh()) != nullptr)
    {
        bool triangle_straddles_left_right_boundary = false;

        c_vector<double, 2> vector1 = this_node_location - previous_node_location;
        if (fabs(vector1[0]) > 0.5*mWidth)
            triangle_straddles_left_right_boundary = true;
        c_vector<double, 2> vector2 = next_node_location - this_node_location;
        if (fabs(vector2[0]) > 0.5*mWidth)
            triangle_straddles_left_right_boundary = true;
        if (triangle_straddles_left_right_boundary)
        {
            if (previous_node_location[0] < mCenterOfWidth)
                previous_node_location[0] += mWidth;
            if (this_node_location[0] < mCenterOfWidth)
                this_node_location[0] += mWidth;
            if (next_node_location[0] < mCenterOfWidth)
                next_node_location[0] += mWidth;
        }
    }
    
    double expanded_triangle_box_bottom = -small_change+std::min(std::min(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
    double expanded_triangle_box_left = -small_change+std::min(std::min(previous_node_location[0],this_node_location[0]),next_node_location[0]);
    double expanded_triangle_box_right = small_change+std::max(std::max(previous_node_location[0],this_node_location[0]),next_node_location[0]);

    bool has_substrate_adhesion_area = true;
    if (expanded_triangle_box_bottom > strip_start_y_location)
        has_substrate_adhesion_area = false;
    else if ((expanded_triangle_box_top >= strip_start_y_location)&&(expanded_triangle_box_bottom >= 0.0)&&(expanded_triangle_box_bottom <= strip_start_y_location))
    {
        sample_area_top = strip_start_y_location;
        sample_area_bottom = expanded_triangle_box_bottom;
    }
    else if ((expanded_triangle_box_top <= strip_start_y_location)&&(expanded_triangle_box_bottom >= 0.0)) // note: we get errors here previously.
    {
        sample_area_top = expanded_triangle_box_top;
        sample_area_bottom = expanded_triangle_box_bottom;
    }
    else if ((expanded_triangle_box_top <= strip_start_y_location)&&(expanded_triangle_box_top >= 0.0)&&(expanded_triangle_box_bottom <= 0.0))
    {
        sample_area_top = expanded_triangle_box_top;
        sample_area_bottom = 0.0;
    }
    else if (expanded_triangle_box_top <= 0.0)
        has_substrate_adhesion_area = false;       
    else if ((expanded_triangle_box_top >= strip_start_y_location)&&(expanded_triangle_box_bottom <= 0.0))                   
    {
        sample_area_top = strip_start_y_location;
        sample_area_bottom = 0.0;
    }
    else
        has_substrate_adhesion_area = false;

    sample_area_left = expanded_triangle_box_left;
    sample_area_right = expanded_triangle_box_right;

    unsigned num_across = 0;
    unsigned num_up = 0;
    unsigned sample_num = 0;
    double sample_area = 0.0;
    if (has_substrate_adhesion_area)
    {
        num_across = (unsigned)round((sample_area_right - sample_area_left) / preferred_sample_dis);
        num_up = (unsigned)round((sample_area_top - sample_area_bottom) / preferred_sample_dis);
        sample_num = num_across * num_up;
        sample_area = (sample_area_top - sample_area_bottom)*(sample_area_right - sample_area_left);
    }

    // Calculate reservoir_substrate_adhesion_area_gradient!
    if (has_substrate_adhesion_area && (num_across>0) && (num_up>0))
    {
        c_vector<c_vector<double, DIM>, 3> points;
        points[0] = previous_node_location;
        points[1] = this_node_location;
        points[2] = next_node_location;
        c_vector<double, DIM> vec1 = this_node_location-previous_node_location;
        c_vector<double, DIM> vec2 = next_node_location-this_node_location;                    

        // Calculate initial adhesive area
        double adhesive_sample_num = 0.0;
        for (unsigned i = 0; i<sample_num; i++)
        {
            double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
            double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

            c_vector<double, DIM> point;
            point[0] = x_coord;
            point[1] = y_coord;
            c_vector<bool, 3> point_at_left_of_vector;
            for (unsigned j = 0; j < 3; j++)
            {
                c_vector<double, DIM> vec1 = point - points[j];
                c_vector<double, DIM> vec2 = points[(j+1)%3] - points[j];

                if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                    point_at_left_of_vector[j] = false;
                else
                    point_at_left_of_vector[j] = true;
            }
            if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                adhesive_sample_num += 1.0;
            else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                adhesive_sample_num += -1.0;
        }
        double substrate_adhesion_area = adhesive_sample_num/double(sample_num) * sample_area;

        // Calculate adhesive area *changed* with small displacement of the node along the x or y axis
        for (unsigned j = 0; j<2; j++)
        {
            c_vector<double, DIM> unit_vector_in_small_change_direction;
            unit_vector_in_small_change_direction[0] = cos(j*M_PI/2);
            unit_vector_in_small_change_direction[1] = sin(j*M_PI/2);

            c_vector<c_vector<double, DIM>, 3> new_points = points;
            new_points[1] += small_change*unit_vector_in_small_change_direction;

            double adhesive_sample_num = 0.0;

            for (unsigned i = 0; i <sample_num; i++)
            {
                double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                c_vector<double, DIM> point;
                point[0] = x_coord;
                point[1] = y_coord;
                c_vector<bool, 3> point_at_left_of_vector;
                for (unsigned j = 0; j < 3; j++)
                {
                    c_vector<double, DIM> vec1 = point - new_points[j];
                    c_vector<double, DIM> vec2 = new_points[(j+1)%3] - new_points[j];

                    if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                        point_at_left_of_vector[j] = false;
                    else
                        point_at_left_of_vector[j] = true;
                }
                if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                    adhesive_sample_num += 1.0;
                else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                    adhesive_sample_num += -1.0;
            }
            double substrate_adhesion_area_new = adhesive_sample_num/double(sample_num) * sample_area;

            reservoir_substrate_adhesion_area_gradient += (substrate_adhesion_area_new - substrate_adhesion_area)/small_change*unit_vector_in_small_change_direction;
        }// end of calculate adhesive area *changed* with small displacement of the node along the x or y axis  
    }// end of statement 'if (has_substrate_adhesion_area)'

    return reservoir_substrate_adhesion_area_gradient;
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // *rParamsFile << "\t\t\t<ReferenceTargetArea>" << mReferenceTargetArea << "</ReferenceTargetArea>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class FaceValueAndStressStateModifier<1>;
template class FaceValueAndStressStateModifier<2>;
template class FaceValueAndStressStateModifier<3>;
