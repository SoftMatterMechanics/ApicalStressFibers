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

#include "MyXToroidal2dVertexMesh.hpp"

MyXToroidal2dVertexMesh::MyXToroidal2dVertexMesh(double centerOfWidth,
                                           double width,
                                           std::vector<Node<2>*> nodes,
                                           std::vector<VertexElement<2, 2>*> vertexElements,
                                           double cellRearrangementThreshold,
                                           double t2Threshold)
    : MutableVertexMesh<2,2>(nodes, vertexElements, cellRearrangementThreshold, t2Threshold),
      mWidth(width),
      mCenterOfWidth(centerOfWidth),
      mpMeshForVtk(nullptr)
{
    // Call ReMesh() to remove any deleted nodes and relabel
    ReMesh();
}

MyXToroidal2dVertexMesh::MyXToroidal2dVertexMesh()
{
}

MyXToroidal2dVertexMesh::~MyXToroidal2dVertexMesh()
{
    if (mpMeshForVtk != NULL)
    {
         delete mpMeshForVtk;
    }
}

c_vector<double, 2> MyXToroidal2dVertexMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);

    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);

    // If the points are more than halfway across the domain, measure the other way
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }

    return vector;
}

void MyXToroidal2dVertexMesh::SetNode(unsigned nodeIndex, ChastePoint<2> point)
{
    double x_coord = point.rGetLocation()[0];

    // Perform a periodic movement if necessary
    if (x_coord >= (mCenterOfWidth+mWidth/2))
    {
        // Move point left
        point.SetCoordinate(0, x_coord - mWidth);
    }
    else if (x_coord < (mCenterOfWidth-mWidth/2))
    {
        // Move point right
        point.SetCoordinate(0, x_coord + mWidth);
    }

    // Update the node's location
    MutableVertexMesh<2,2>::SetNode(nodeIndex, point);
}

double MyXToroidal2dVertexMesh::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension==0 || rDimension==1);

    double width = mWidth;

    return width;
}

void MyXToroidal2dVertexMesh::SetWidth(double width)
{
    assert(width > 0);
    mWidth = width;
}

unsigned MyXToroidal2dVertexMesh::AddNode(Node<2>* pNewNode)
{
    unsigned node_index = MutableVertexMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back onto the torus
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;
}

VertexMesh<2, 2>* MyXToroidal2dVertexMesh::GetMeshForVtk()
{
    unsigned num_nodes = GetNumNodes();

    std::vector<Node<2>*> temp_nodes(2*num_nodes);
    std::vector<VertexElement<2, 2>*> elements;

    // Create two copies of each node
    for (unsigned index=0; index<num_nodes; index++)
    {
        c_vector<double, 2> location;
        location = GetNode(index)->rGetLocation();

        // Node copy at original location
        Node<2>* p_node = new Node<2>(index, false, (location[0] + this->mMoveMeshRightForNPeriods*mWidth)*(this->mMultiplyResultsBy), location[1]*(this->mMultiplyResultsBy));
        temp_nodes[index] = p_node;

        // Node copy shifted right
        p_node = new Node<2>(num_nodes + index, false, (location[0] + mWidth + this->mMoveMeshRightForNPeriods*mWidth)*(this->mMultiplyResultsBy), location[1]*(this->mMultiplyResultsBy));
        temp_nodes[num_nodes + index] = p_node;

    }

    // Iterate over elements
    for (VertexMesh<2,2>::VertexElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        unsigned num_nodes_in_elem = elem_iter->GetNumNodes();

        std::vector<Node<2>*> elem_nodes;

        // Compute whether the element straddles either periodic boundary
        bool element_straddles_left_right_boundary = false;

        const c_vector<double, 2>& r_this_node_location = elem_iter->GetNode(0)->rGetLocation();
        for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        {
            const c_vector<double, 2>& r_next_node_location = elem_iter->GetNode((local_index+1)%num_nodes_in_elem)->rGetLocation();
            c_vector<double, 2> vector;
            vector = r_next_node_location - r_this_node_location;

            if (fabs(vector[0]) > 0.5*mWidth)
            {
                element_straddles_left_right_boundary = true;
            }
        }
        // // My changes
        // bool element_at_left_of_width_center = true;
        // for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        // {
        //     const c_vector<double, 2>& r_node_location = elem_iter->GetNode(local_index)->rGetLocation();
        //     if (r_node_location[0] > this->mCenterOfWidth)
        //     {
        //         element_at_left_of_width_center = false;
        //         break;
        //     }
        // }


        // Use the above information when duplicating the element
        for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(local_index);

            // If the element straddles the left/right periodic boundary...
            if (element_straddles_left_right_boundary)
            {
                // ...and this node is located to the left of the centre of the mesh...
                bool node_is_left_of_centre = (elem_iter->GetNode(local_index)->rGetLocation()[0] - mCenterOfWidth < 0);
                if (node_is_left_of_centre)
                {
                    // ...then choose the equivalent node to the right
                    this_node_index += num_nodes;
                }
            }
            // // my changes
            // if (element_at_left_of_width_center)
            // {
            //     this_node_index += num_nodes;
            // }
            
            elem_nodes.push_back(temp_nodes[this_node_index]);
        }

        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index, elem_nodes);
        elements.push_back(p_element);
    }


    // Now delete any nodes from the mesh for VTK that are not contained in any elements
    std::vector<Node<2>*> nodes;
    unsigned count = 0;
    for (unsigned index=0; index<temp_nodes.size(); index++)
    {
        unsigned num_elems_containing_this_node = temp_nodes[index]->rGetContainingElementIndices().size();

        if (num_elems_containing_this_node == 0)
        {
            // Avoid memory leak
            delete temp_nodes[index];
        }
        else
        {
            temp_nodes[index]->SetIndex(count);
            nodes.push_back(temp_nodes[index]);
            count++;
        }
    }

    // //My changes:
    // for (unsigned node_index = 0; node_index != num_nodes; node_index++)
    // {
    //     c_vector<double, 2> location_right;
    //     location_right = nodes[node_index]->rGetLocation();
    //     location_right[0] += mWidth;
    //     Node<2>* p_node = new Node<2>(num_nodes + nodes[node_index]->GetIndex(), false, location_right[0], location_right[1]);
    //     nodes.push_back(p_node);//****************
    // }

    // unsigned num_elements = elements.size();
    // for (unsigned index = 0; index < num_elements; index++)
    // {
    //     VertexElement<2,2>* p_element = elements[index];
    //     unsigned elem_index = p_element->GetIndex();
    //     unsigned num_nodes_in_elem = p_element->GetNumNodes();
    //     std::vector<Node<2>*> elem_nodes;
    //     for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
    //     {
    //         unsigned this_node_index = p_element->GetNodeGlobalIndex(local_index);
    //         elem_nodes.push_back(nodes[this_node_index+num_nodes]);
    //     }
    //     VertexElement<2,2>* p_new_element = new VertexElement<2,2>(num_elements+elem_index, elem_nodes);
    //     elements.push_back(p_new_element);//***************************
    // }


    if (mpMeshForVtk != nullptr)
    {
        delete mpMeshForVtk;
    }

    mpMeshForVtk = new VertexMesh<2,2>(nodes, elements);
    return mpMeshForVtk;
}

void MyXToroidal2dVertexMesh::ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader, double width)
{
    assert(rMeshReader.HasNodePermutation() == false);

    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
    {
        node_data = rMeshReader.GetNextNode();
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(node_idx, node_data, false));
    }

    rMeshReader.Reset();

    // Reserve memory for elements
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_idx = 0 ; elem_idx < num_elements ; elem_idx++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j=0; j<num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_idx, nodes);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = (unsigned) element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }

    /*
     * Set width and height from function arguments, and validate by checking area is correct
     */
    this->mWidth = width;

    double total_surface_area = 0.0;
    for (unsigned elem_idx = 0 ; elem_idx < num_elements ; elem_idx++)
    {
        total_surface_area += this->GetVolumeOfElement(elem_idx);
    }

    // Set default parameter values
    this->mCellRearrangementRatio = 1.5;
    this->mCellRearrangementThreshold = 0.01;
    this->mT2Threshold = 0.001;
    this->mMeshChangesDuringSimulation = true;
    this->mpMeshForVtk = nullptr;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyXToroidal2dVertexMesh)
