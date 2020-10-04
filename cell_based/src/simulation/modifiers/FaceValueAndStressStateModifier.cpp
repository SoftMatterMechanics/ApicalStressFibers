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
      mIfConsiderFeedbackOfCellCellAdhesion(false),

      mEMADontDecreaseWhenEdgeShrink(false),
      mCCADontDecreaseWhenEdgeExpand(false),
      mCCAIncreasingHasAThresholdOfEdgeLength(false),
      mCCAIncreasingThresholdOfEdgeLengthPercentage(0.5),

      mEdgeLengthAtRest(sqrt(M_PI/(6*sqrt(3)/4))),
      mFeedbackStrengthForMyosinActivity(1.0),
      mHillCoefficientForMyosinActivity(8.0),
      mFeedbackStrengthForAdhesion(1.0),
      mHillCoefficientForAdhesion(8.0),

      mIfCalculateStressState(true),
      mCaseNumberOfMembraneSurfaceEnergyForm(0),
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
    if (mIfCalculateStressState)
    {
        // Loop over the list of cells, rather than using the population iterator, so as to include(exclude??) dead cells
        for (std::list<CellPtr>::iterator cell_iter = rCellPopulation.rGetCells().begin();
            cell_iter != rCellPopulation.rGetCells().end();
            ++cell_iter)
        {
            UpdateStressStateOfCell(rCellPopulation, *cell_iter);
        }
    }

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
                    this->UpdateUnifiedEdgeMyosinActivtyOfFace(p_mesh, face_index);
                    
                    if (mIfConsiderFeedbackOfCellCellAdhesion)
                        this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                }
            }
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

    double stress_XX = 0.0;
    double stress_YY = 0.0;
    double stress_XY = 0.0;
    double stress_1 = 0.0;
    double stress_2 = 0.0;
    if (mCaseNumberOfMembraneSurfaceEnergyForm >= 0u)
    {
        unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
        VertexElement<DIM, DIM>* p_element = p_mesh->GetElement(elem_index);
        double area = p_mesh->GetVolumeOfElement(elem_index);
        // double perimeter = p_mesh->GetSurfaceAreaOfElement(elem_index);
        double Pressure_ = area/(this->mFixedTargetArea)-1;
        double Ga_= 1/(this->mFixedTargetArea)*this->mNagaiHondaMembraneSurfaceEnergyParameter;
        double myosin_weighted_perimeter = 0.0;
        double sum_XX = 0.0;
        double sum_YY = 0.0;
        double sum_XY = 0.0;
        double sum_adhe_XX = 0.0;
        double sum_adhe_YY = 0.0;
        double sum_adhe_XY = 0.0;

        // for (unsigned i=0; i<8; i++)
        // {
        //     std::string name_item;
        //     std::ostringstream oss;
        //     oss.str("");
        //     oss << i << "_vx";
        //     name_item = oss.str();
        //     if (pCell->GetCellData()->HasItem(name_item))
        //         pCell->GetCellData()->DeleteItem(name_item);
        //     oss.str("");
        //     oss << i << "_vy";
        //     name_item = oss.str();
        //     if (pCell->GetCellData()->HasItem(name_item))
        //         pCell->GetCellData()->DeleteItem(name_item);
        // }

        unsigned local_index_lowest_vertex = 0;
        double y_coor_lowest_vertex = 10000.0;
        for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
        {
            Node<DIM>* pNode = p_element->GetNode(local_index);
            if (pNode->rGetLocation()[1]<y_coor_lowest_vertex)
            {
                local_index_lowest_vertex = local_index;
                y_coor_lowest_vertex = pNode->rGetLocation()[1];
            }
        }

        for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
        {
            Node<DIM>* pNodeA = p_element->GetNode(local_index);
            Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
            c_vector<double,DIM> location_a = pNodeA->rGetLocation();
            c_vector<double,DIM> location_b = pNodeB->rGetLocation();
            double l_ab = p_mesh->GetDistanceBetweenNodes(pNodeA->GetIndex(), pNodeB->GetIndex());
            c_vector<double,DIM> vec_ab =  p_mesh->GetVectorFromAtoB(location_a, location_b);

            VertexElement<DIM-1,  DIM>* pFace = nullptr;
            if (mIfConsiderFeedbackOfFaceValues)
                pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeA->GetIndex(), pNodeB->GetIndex()));

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
            if (shared_elements.size() >= 3)
            {
                std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateStressStateOfCell";
                std::cout<< std::endl << "Get shared elements more than 2";
            }
            // output velocity of vertices
            std::string name_item;
            std::ostringstream oss;

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

            if (mIfConsiderFeedbackOfFaceValues)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_myo";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, pFace->GetUnifiedEdgeMyosinActivty());
            }
            if (mIfConsiderFeedbackOfFaceValues && mIfConsiderFeedbackOfCellCellAdhesion)
            {
                oss.str("");
                oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                        << "_cc";
                name_item = oss.str();
                pCell->GetCellData()->SetItem(name_item, pFace->GetUnifiedCellCellAdhesionEnergyParameter());
            }

            double lambda = 0.0;
            if (shared_elements.size() == 1)
                lambda = 0.0;
            else
                lambda = mNagaiHondaCellCellAdhesionEnergyParameter;
            if (mIfConsiderFeedbackOfFaceValues)
                lambda *= pFace->GetUnifiedCellCellAdhesionEnergyParameter();
            sum_adhe_XX += lambda/pow(this->mFixedTargetArea,1.5)*l_ab/sqrt(this->mFixedTargetArea)*(vec_ab[0]/l_ab)*(vec_ab[0]/l_ab);
            sum_adhe_YY += lambda/pow(this->mFixedTargetArea,1.5)*l_ab/sqrt(this->mFixedTargetArea)*(vec_ab[1]/l_ab)*(vec_ab[1]/l_ab);
            sum_adhe_XY += lambda/pow(this->mFixedTargetArea,1.5)*l_ab/sqrt(this->mFixedTargetArea)*(vec_ab[0]/l_ab)*(vec_ab[1]/l_ab);

            if (mIfConsiderFeedbackOfFaceValues)
            {
                double m_ab = pFace->GetUnifiedEdgeMyosinActivty();
                myosin_weighted_perimeter += sqrt(m_ab)*l_ab/sqrt(this->mFixedTargetArea);

                sum_XX += sqrt(m_ab)*(vec_ab[0]/l_ab)*vec_ab[0]/sqrt(this->mFixedTargetArea);
                sum_YY += sqrt(m_ab)*(vec_ab[1]/l_ab)*vec_ab[1]/sqrt(this->mFixedTargetArea);
                sum_XY += sqrt(m_ab)*(vec_ab[0]/l_ab)*vec_ab[1]/sqrt(this->mFixedTargetArea);
            }
            else
            {
                myosin_weighted_perimeter += l_ab/sqrt(this->mFixedTargetArea);
                sum_XX += (vec_ab[0]/l_ab)*vec_ab[0]/sqrt(this->mFixedTargetArea);
                sum_YY += (vec_ab[1]/l_ab)*vec_ab[1]/sqrt(this->mFixedTargetArea);
                sum_XY += (vec_ab[0]/l_ab)*vec_ab[1]/sqrt(this->mFixedTargetArea);
            }
        } // end of iteration of vertices of the element, for calculation of summation.

        stress_XX = Pressure_ + Ga_/(area/(this->mFixedTargetArea))*myosin_weighted_perimeter*sum_XX + 1/(area/(this->mFixedTargetArea))*sum_adhe_XX;
        stress_YY = Pressure_ + Ga_/(area/(this->mFixedTargetArea))*myosin_weighted_perimeter*sum_YY + 1/(area/(this->mFixedTargetArea))*sum_adhe_YY;
        stress_XY = Ga_/(area/(this->mFixedTargetArea))*myosin_weighted_perimeter*sum_XY + 1/(area/(this->mFixedTargetArea))*sum_adhe_XY;
        // up is the stress normalized using the literature method(A0_=1); below is the stress normalized using our method(A0_=PI);
        stress_XX *= (this->mFixedTargetArea);
        stress_YY *= (this->mFixedTargetArea);
        stress_XY *= (this->mFixedTargetArea);

        double R = sqrt( pow((stress_XX-stress_YY)/2, 2) + pow(stress_XY, 2) );
        stress_1 = (stress_XX + stress_YY)/2 + R;
        stress_2 = (stress_XX + stress_YY)/2 - R;

        if (mIfOutputModifierInformation && random()%10000==0)
        {
            double t = SimulationTime::Instance()->GetTime();
            c_vector<double, DIM> centroid = p_mesh->GetCentroidOfElement(elem_index);
            std::cout << std::endl << "1||t=" << t << "|elem_index=" << elem_index << "|centroid=" << centroid[0] << ", " << centroid[1];
            std::cout << std::endl << "2||stressXX=" << stress_XX << "|stressYY=" << stress_YY << "|stressXY=" << stress_XY;
        }

    }
    
    
    pCell->GetCellData()->SetItem("StressXX", stress_XX);
    pCell->GetCellData()->SetItem("StressYY", stress_YY);
    pCell->GetCellData()->SetItem("StressXY", stress_XY);
    pCell->GetCellData()->SetItem("Stress1", stress_1);
    pCell->GetCellData()->SetItem("Stress2", stress_2);

}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::UpdateUnifiedEdgeMyosinActivtyOfFace(MutableVertexMesh<DIM, DIM>* pMesh, unsigned faceIndex)
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    double feedback_strength = this->mFeedbackStrengthForMyosinActivity;
    double n = this->mHillCoefficientForMyosinActivity;
    double edge_length_at_rest = this->mEdgeLengthAtRest;
    double alpha = 2.0*feedback_strength;
    double beta = 1.0*feedback_strength;
    double KL = 1.0;

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
        if (fabs(p_mesh->GetCentroidOfElement(ele_index)[0] - 0.0)< 1e-1)
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
