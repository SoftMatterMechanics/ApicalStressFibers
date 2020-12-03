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

template<unsigned DIM>
FaceValueAndStressStateModifier<DIM>::FaceValueAndStressStateModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mIfConsiderFeedbackOfFaceValues(false),
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

      mEdgeLengthAtRest(sqrt(M_PI/(6*sqrt(3)/4))),
      mKLForFeedback(1.0),
      mFeedbackStrengthForMyosinActivity(1.0),
      mHillCoefficientForMyosinActivity(8.0),
      mFeedbackStrengthForAdhesion(1.0),
      mHillCoefficientForAdhesion(8.0),

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

      mMarkLeadingCells(false),
      mMultipleLeadingCells(false),
      mLeadingCellNumber(1),

      mIfOutputModifierInformation(false)
{
}

template<unsigned DIM>
FaceValueAndStressStateModifier<DIM>::~FaceValueAndStressStateModifier()
{
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateFaceValuesAndStressStates(rCellPopulation);
    UpdateCellAreas(rCellPopulation);
    if (mUseMyDivisionRuleAlongWithModifier)
        UpdateForCellDivision(rCellPopulation);
    if (mWriteGroupNumberToCell)
        UpdateGroupNumbers(rCellPopulation);
    if (mMarkLeadingCells)
        UpdateLamellipodiumInfoOfCells(rCellPopulation);
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateFaceValuesAndStressStates(rCellPopulation);
    UpdateCellAreas(rCellPopulation);
    if (mUseMyDivisionRuleAlongWithModifier)
        SetupSolveForCellDivision(rCellPopulation);
    if (mWriteGroupNumberToCell)
        UpdateGroupNumbers(rCellPopulation);
    if (mMarkLeadingCells)
        SetupSolveForLamellipodiumInfoOfCells(rCellPopulation);
}

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

    if (mIfCalculateStressState)
    {
        // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
        for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
            cell_iter != rCellPopulation.rGetCells().end();
            ++cell_iter)
        {
            UpdateStressStateOfCell(rCellPopulation, *cell_iter);
        }

        for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
            cell_iter != rCellPopulation.rGetCells().end();
            ++cell_iter)
        {
            UpdateCellDataOfForcesFromNeighboringCell(rCellPopulation, *cell_iter);
        }
        
    }

}

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
                pCell->GetCellData()->SetItem(name_item, pFace->GetUnifiedEdgeMyosinActivty());
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
            std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateStressStateOfCell.";
            std::cout<< std::endl << "Get shared elements more than 2 or less than 1.";
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

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateUnifiedEdgeMyosinActivtyOfFace(MutableVertexMesh<DIM, DIM>* pMesh, unsigned faceIndex)
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    double edge_length_at_rest = this->mEdgeLengthAtRest;
    double KL = this->mKLForFeedback;
    double n = this->mHillCoefficientForMyosinActivity;
    double feedback_strength = this->mFeedbackStrengthForMyosinActivity;
    double alpha = (1.0+KL)*feedback_strength;
    double beta = 1.0*feedback_strength;

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
    unified_edge_myosin_activty += dt*changing_rate;
    p_face->SetUnifiedEdgeMyosinActivty(unified_edge_myosin_activty);
    if (mIfOutputModifierInformation && random()%100000==0)
    {
        std::cout << std::endl<< "In FaceValueAndStressStateModifier::UpdateUnifiedEdgeMyosinActivityOfFace: the new UnifiedEdgeMyosinActivty=";
        std::cout << p_face->GetUnifiedEdgeMyosinActivty();
        std::cout << std::endl<< "feedback_strength=" << feedback_strength << " alpha=" << alpha << " beta=" <<beta << " KL=" << KL;
        std::cout << " n=" << n << " edge_length_at_rest=" << edge_length_at_rest << " the edge length=" << edge_length << " lambda=" << lambda;
        std::cout << " changing_rate=" << changing_rate;
    }
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(MutableVertexMesh<DIM, DIM>* pMesh, unsigned faceIndex)
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    double feedback_strength = this->mFeedbackStrengthForAdhesion;
    double n = this->mHillCoefficientForAdhesion;
    double edge_length_at_rest = this->mEdgeLengthAtRest;
    double alpha = 2.0*feedback_strength;
    double beta = 1.0*feedback_strength;
    double KL = 1.0;

    VertexElement<DIM-1, DIM>* p_face = pMesh->GetFace(faceIndex);
    double edge_length = pMesh->GetDistanceBetweenNodes(p_face->GetNodeGlobalIndex(0), p_face->GetNodeGlobalIndex(1));
    double lambda = edge_length/ edge_length_at_rest;//(1.0*sqrt(M_PI/(6*sqrt(3)/4)));
    // // tmp
    // std::cout << "lambda=" << lambda << std::endl;
    double corresponding_lambda = 2.0-lambda;
    double unified_cell_cell_adhesion_energy_parameter = p_face->GetUnifiedCellCellAdhesionEnergyParameter();
    double changing_rate = alpha*pow(corresponding_lambda,n)/(KL+pow(corresponding_lambda,n))-beta*unified_cell_cell_adhesion_energy_parameter;
    if (mCCADontDecreaseWhenEdgeExpand && lambda>1.0)
    {
        changing_rate = 0.0;
    }
    if (mCCAIncreasingHasAThresholdOfEdgeLength && lambda>mCCAIncreasingThresholdOfEdgeLengthPercentage)
    {
        changing_rate = 0.0;
    }
    unified_cell_cell_adhesion_energy_parameter += changing_rate*dt;
    p_face->SetUnifiedCellCellAdhesionEnergyParameter(unified_cell_cell_adhesion_energy_parameter);
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
        if (    ( (p_mesh->GetCentroidOfElement(ele_index)[1] > mStripStartYLocation/2.0)&&(!mMultipleLeadingCells)&&(fabs(p_mesh->GetCentroidOfElement(ele_index)[0] - mStripStartXLocation)< 1e-1) ) 
             || ( (p_mesh->GetCentroidOfElement(ele_index)[1] > mStripStartYLocation/2.0)&&mMultipleLeadingCells&&(fabs(p_mesh->GetCentroidOfElement(ele_index)[0] - mStripStartXLocation)< mStripWidth/2.0) ) )
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
