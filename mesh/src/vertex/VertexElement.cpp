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
#include "VertexElement.hpp"
#include <cassert>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> // the choice of constructor in my simulation!
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
      mFaces(rFaces),
      mOrientations(rOrientations),
      mGroupNumber(0),
      mIsLeadingCell(false),
      mIsLeadingCellTop(false),
      mIsLeadingCellBottom(false),
      mIsJustReAttached(false),
      mLamellipodiumStrength(0.0),
      mElementMyosinActivity(1.0)
{
    // This constructor should only be used in 3D
    assert(SPACE_DIM == 3 || (ELEMENT_DIM==2 && SPACE_DIM==2));    // LCOV_EXCL_LINE - code will be removed at compile time

    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    if (SPACE_DIM == ELEMENT_DIM)
    {
        // Register element with nodes
        this->RegisterWithNodes();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index),
      mFaces(rFaces),
      mOrientations(rOrientations)
{
    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    // Make a set of nodes with mFaces
    std::set<Node<SPACE_DIM>* > nodes_set;
    for (unsigned face_index=0; face_index<mFaces.size(); face_index++)
    {
        for (unsigned node_index=0; node_index<mFaces[face_index]->GetNumNodes(); node_index++)
        {
            nodes_set.insert(mFaces[face_index]->GetNode(node_index));
        }
    }

    // Populate mNodes
    for (typename std::set< Node<SPACE_DIM>* >::iterator node_iter = nodes_set.begin();
         node_iter != nodes_set.end();
         ++node_iter)
    {
         this->mNodes.push_back(*node_iter);
    }

    // Register element with nodes
    this->RegisterWithNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
      mGroupNumber(0),
      mIsLeadingCell(false),
      mIsLeadingCellTop(false),
      mIsLeadingCellBottom(false),
      mIsJustReAttached(false),
      mLamellipodiumStrength(0.0),
      mElementMyosinActivity(1.0)
{
}

// my addition
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> // the choice of constructor in my simulation!
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rStressfibers)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
      mFaces(rFaces),
      mOrientations(rOrientations),
      mStressfibers(rStressfibers),
      mGroupNumber(0),
      mIsLeadingCell(false),
      mIsLeadingCellTop(false),
      mIsLeadingCellBottom(false),
      mIsJustReAttached(false),
      mLamellipodiumStrength(0.0),
      mElementMyosinActivity(1.0)
{
    // This constructor should only be used in 3D
    assert(SPACE_DIM == 3 || (ELEMENT_DIM==2 && SPACE_DIM==2));    // LCOV_EXCL_LINE - code will be removed at compile time

    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    if (SPACE_DIM == ELEMENT_DIM)
    {
        // Register element with nodes
        this->RegisterWithNodes();
    }
}

// my addition for construction of a stress fiber
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index, 
                                                    const std::vector<Node<SPACE_DIM>*>& rStressfiberNodes, 
                                                    c_vector<double,2> rEndpointsratio,
                                                    double rRestLength,
                                                    bool rIsPeeled,
                                                    double rPeelingTime)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index)
{
    mStressfiberGlobalIndex = index;
    mStressfiberNodes=rStressfiberNodes;
    mEndpointsratio=rEndpointsratio;
    mRestLength=rRestLength;
    mIsPeeled=rIsPeeled;
    mPeelingTime=rPeelingTime;
} 

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::~VertexElement()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddFace(VertexElement<ELEMENT_DIM-1,SPACE_DIM>* pFace)
{
    // Add pFace to the end of mFaces
    this->mFaces.push_back(pFace);

    // Create a set of indices of nodes currently owned by this element
    std::set<unsigned> node_indices;
    for (unsigned local_index=0; local_index<this->GetNumNodes(); local_index++)
    {
        node_indices.insert(this->GetNodeGlobalIndex(local_index));
    }

    // Loop over nodes owned by pFace
    unsigned end_index = this->GetNumNodes()-1;
    for (unsigned local_index=0; local_index<pFace->GetNumNodes(); local_index++)
    {
        // If this node is not already owned by this element...
        unsigned global_index = pFace->GetNodeGlobalIndex(local_index);
        if (node_indices.find(global_index) == node_indices.end())
        {
            // ... then add it to the element (and vice versa)
            this->AddNode(pFace->GetNode(local_index), end_index);
            end_index++;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM-1,  SPACE_DIM>* VertexElement<ELEMENT_DIM, SPACE_DIM>::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
    return mFaces[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceIsOrientatedClockwise(unsigned index) const
{
    assert(index < mOrientations.size());
    return mOrientations[index];
}

// my addition
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddStressfiber(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pStressfiber)
{
    // Add pStressfiber to the end of mStressfibers
    this->mStressfibers.push_back(pStressfiber);
}

// my addition
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteStressfiber(unsigned index)
{
    // delete the stress fiber with specified local index in the element
    this->mStressfibers.erase(this->mStressfibers.begin()+index);
}

// my addition
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM-1, SPACE_DIM>* VertexElement<ELEMENT_DIM, SPACE_DIM>::GetStressfiber(unsigned index) const
{
    // return a pointer to the stress fiber with specified local index
    assert(index < mStressfibers.size());
    return mStressfibers[index];
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetStressfiberGlobalIndex()
{
    return 0;
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNumStressfibers() const
{
    return mStressfibers.size();
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* VertexElement<ELEMENT_DIM, SPACE_DIM>::GetStressfiberNode(unsigned index)
{
    return nullptr;
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,2> VertexElement<ELEMENT_DIM, SPACE_DIM>::GetStressfiberEndpointsratio()
{
    return zero_vector<double>(2);
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetStressfiberRestLength()
{
    return 0.0;
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::GetStressfiberPeelStatus()
{
    return false;
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetStressfiberPeelingTime()
{
    return 0.0;
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateStressfiberEndpointsratio(c_vector<double,2> endpointsratio)
{
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateStressfiberNode(Node<SPACE_DIM>* newStressfiberNode, unsigned index)
{
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateStressfiberRestLength(double newRestLength)
{
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateStressfiberPeelStatus(bool newPeelStatus)
{
}

// my addition
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateStressfiberPeelingTime(double newPeelingTime)
{
}

//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

// my addition
/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
VertexElement<1, SPACE_DIM>::VertexElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<1, SPACE_DIM>(index, rNodes)
{
    if (SPACE_DIM==2)
    {
        mCellCellAdhesionEnergyParameter = 1.0;
        mEdgeMyosinActivty = 1.0;
        mIsDeleted = false;
    }
}// this is the typical constructor for a face object!

// my addition
template<unsigned SPACE_DIM>
VertexElement<1, SPACE_DIM>::VertexElement(unsigned index, 
                                           const std::vector<Node<SPACE_DIM>*>& rStressfiberNodes, 
                                           c_vector<double,2> rEndpointsratio,
                                           double rRestLength,
                                           bool rIsPeeled,
                                           double rPeelingTime)
    : MutableElement<1, SPACE_DIM>(index)
{
    if (SPACE_DIM==2)
    {
        mStressfiberGlobalIndex = index;
        mStressfiberNodes=rStressfiberNodes;
        mEndpointsratio=rEndpointsratio;
        mRestLength=rRestLength;
        mIsPeeled=rIsPeeled;
        mPeelingTime=rPeelingTime;
    }
} // this is the typical constructor for a stress fiber!


template<unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::GetNumFaces() const
{
    return 0;
}

template<unsigned SPACE_DIM>
VertexElement<0, SPACE_DIM>* VertexElement<1, SPACE_DIM>::GetFace(unsigned index) const
{
    return nullptr;
}

template<unsigned SPACE_DIM>
bool VertexElement<1, SPACE_DIM>::FaceIsOrientatedClockwise(unsigned index) const
{
    return false;
}

// my addition
template<unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::AddStressfiber(VertexElement<1-1, SPACE_DIM>* pStressfiber)
{
}

// my addition
template<unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::DeleteStressfiber(unsigned index)
{
}

// my addition
template<unsigned SPACE_DIM>
VertexElement<0, SPACE_DIM>* VertexElement<1, SPACE_DIM>::GetStressfiber(unsigned index) const
{
    return nullptr;
}

// my addition
template <unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::GetNumStressfibers() const
{
    return 0;
}

// my addition
template <unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::GetStressfiberGlobalIndex()
{
    return mStressfiberGlobalIndex;
}

// my addition
template <unsigned SPACE_DIM>
Node<SPACE_DIM>* VertexElement<1, SPACE_DIM>::GetStressfiberNode(unsigned index)
{
    assert(index < mStressfiberNodes.size());
    return mStressfiberNodes[index];
}

// my addition
template <unsigned SPACE_DIM>
c_vector<double,2> VertexElement<1, SPACE_DIM>::GetStressfiberEndpointsratio()
{
    return mEndpointsratio;
}

// my addition
template <unsigned SPACE_DIM>
double VertexElement<1, SPACE_DIM>::GetStressfiberRestLength()
{
    return mRestLength;
}

// my addition
template <unsigned SPACE_DIM>
bool VertexElement<1, SPACE_DIM>::GetStressfiberPeelStatus()
{
    return mIsPeeled;
}

// my addition
template <unsigned SPACE_DIM>
double VertexElement<1, SPACE_DIM>::GetStressfiberPeelingTime()
{
    return mPeelingTime;
}

// my addition
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::UpdateStressfiberNode(Node<SPACE_DIM>* newStressfiberNode, unsigned index)
{
    assert(index < mStressfiberNodes.size());
    mStressfiberNodes[index] = newStressfiberNode;
}

// my addition
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::UpdateStressfiberEndpointsratio(c_vector<double,2> endpointsratio)
{
    mEndpointsratio = endpointsratio;
}

// my addition
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::UpdateStressfiberRestLength(double newRestLength)
{
    mRestLength= newRestLength;
}

// my addition
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::UpdateStressfiberPeelStatus(bool newPeelStatus)
{
    mIsPeeled= newPeelStatus;
}

// my addition
template <unsigned SPACE_DIM>
void VertexElement<1, SPACE_DIM>::UpdateStressfiberPeelingTime(double newPeelingTime)
{
    mPeelingTime= newPeelingTime;
}

// Explicit instantiation
template class VertexElement<1,1>;
template class VertexElement<1,2>;
template class VertexElement<1,3>;
template class VertexElement<2,2>;
template class VertexElement<2,3>;
template class VertexElement<3,3>;
