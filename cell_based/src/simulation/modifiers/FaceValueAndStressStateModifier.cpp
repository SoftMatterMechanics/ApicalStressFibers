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
      mNagaiHondaDeformationEnergyParameter(1.0),
      mNagaiHondaMembraneSurfaceEnergyParameter(0.1),
      mNagaiHondaCellCellAdhesionEnergyParameter(0.0),
      mNagaiHondaCellBoundaryAdhesionEnergyParameter(0.0),

      mCaseNumberOfMembraneSurfaceEnergyForm(0),
      mFixedTargetArea(M_PI),
      mFixedTargetPerimeter(6*sqrt(M_PI/(6*sqrt(3)/4))),

      mFeedbackStrengthForMyosinActivity(1.0),
      mHillCoefficientForMyosinActivity(2.0),
      mEdgeLengthAtRest(sqrt(M_PI/(6*sqrt(3)/4))),
      mFeedbackStrengthForAdhesion(1.0),
      mHillCoefficientForAdhesion(2.0),
      mEMADontDecreaseWhenEdgeShrink(true),
      mCCADontDecreaseWhenEdgeExpand(true),
      mCCADontInreaseWhenEdgeShrink(false),
      mCCADontInreaseUntilShorterThanAThreshold(false),
      mCCADontInreaseUntilShorterThanThisValue(0.2),
      mIfOutputModifierInformation(false),
      mIfCalculateStressState(true),
      mIfConsiderFeedbackOfFaceValues(false),
      mIfConsiderFeedbackOfFaceValuesOnlyForBoundaryCells(false),
      mIfConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells(false)
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
}

template<unsigned DIM>
void FaceValueAndStressStateModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateFaceValuesAndStressStates(rCellPopulation);
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

                                this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                            }

                        }
                        else
                        {
                            this->UpdateUnifiedEdgeMyosinActivtyOfFace(p_mesh, face_index);

                            this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                        }
                        
                    }
                }
                else
                {
                    this->UpdateUnifiedEdgeMyosinActivtyOfFace(p_mesh, face_index);

                    this->UpdateUnifiedCellCellAdhesionEnergyParameterOfFace(p_mesh, face_index);
                }
            }
        }
    }
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
    if (mCCADontInreaseUntilShorterThanAThreshold && lambda>mCCADontInreaseUntilShorterThanThisValue)
    {
        changing_rate = 0.0;
    }
    if (mCCADontInreaseWhenEdgeShrink && lambda<1.0)
    {
        changing_rate = 0.0;
    }

    unified_cell_cell_adhesion_energy_parameter += changing_rate*dt;
    p_face->SetUnifiedCellCellAdhesionEnergyParameter(unified_cell_cell_adhesion_energy_parameter);
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
        double Pai = -this->mNagaiHondaDeformationEnergyParameter*(area-this->mFixedTargetArea);
        double sum_XX = 0.0;
        double sum_YY = 0.0;
        double sum_XY = 0.0;
        for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
        {
            Node<DIM>* pNodeA = p_element->GetNode(local_index);
            Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
            c_vector<double,DIM> location_a = pNodeA->rGetLocation();
            c_vector<double,DIM> location_b = pNodeB->rGetLocation();
            double l_ab = p_mesh->GetDistanceBetweenNodes(pNodeA->GetIndex(), pNodeB->GetIndex());
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
            double elem1_perimeter = p_mesh->GetSurfaceAreaOfElement(*shared_elements.begin());
            double elem2_perimeter = 0.0;
            double T_ab = this->mNagaiHondaMembraneSurfaceEnergyParameter*(elem1_perimeter-this->mFixedTargetPerimeter);
            if (shared_elements.size() == 2)
            {
                elem2_perimeter = p_mesh->GetSurfaceAreaOfElement(*shared_elements.end());
                T_ab += this->mNagaiHondaMembraneSurfaceEnergyParameter*(elem2_perimeter-this->mFixedTargetPerimeter);
            }
            if (shared_elements.size() >= 3)
            {
                std::cout<< std::endl << "Get error in FaceValueAndStressStateModifier::UpdateStressStateOfCell";
                std::cout<< std::endl << "Get shared elements more than 2";
            }
            sum_XX += (T_ab*vec_ab[0]/l_ab)*vec_ab[0];
            sum_YY += (T_ab*vec_ab[1]/l_ab)*vec_ab[1];
            sum_XY += (T_ab*vec_ab[0]/l_ab)*vec_ab[1];

            // if (random()%10000==0)
            // {
            //     double t = SimulationTime::Instance()->GetTime();
            //     std::cout << std::endl << "1||t=" << t << "|elem_index=" << elem_index << "|area=" << area;
            //     std::cout << "|Pai=" << Pai << "|local_index=" << local_index << "|l_ab=" << l_ab << "|vec_ab=" << vec_ab[0] << ", " << vec_ab[1];
            //     std::cout << "|elem1_perimeter=" << elem1_perimeter << "|elem2_perimeter=" << elem2_perimeter << "|T_ab=" << T_ab;
            // }
        }
        stress_XX = -Pai + 1/(2*area)*sum_XX;
        stress_YY = -Pai + 1/(2*area)*sum_YY;
        stress_XY = 1/(2*area)*sum_XY;
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
