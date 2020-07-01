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
#ifndef VERTEXELEMENT_HPP_
#define VERTEXELEMENT_HPP_

#include "MutableElement.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * An element class for use in the VertexMesh class. The main
 * difference between this and the Element class is that a
 * VertexElement can have a variable number of nodes associated
 * with it.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexElement : public MutableElement<ELEMENT_DIM, SPACE_DIM>
{
private:

    /**
     * Faces of the VertexElement, which should be distinct.
     */
    std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*> mFaces;

    /**
     * How each face is oriented. From the perspective of the centre
     * of the element, the vertices of each face should be ordered
     * anti clockwise. If and only if this is false, the order of vertices
     * in the corresponding face should be reversed.
     *
     * N.B. Most faces belong to two elements, but with opposite
     * orientations. This allows us to reuse the face data across the
     * two cells.
     */
    std::vector<bool> mOrientations;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This needs to be first so that MeshBasedCellPopulation::Validate() doesn't go mental.
        archive & mFaces;
        archive & mOrientations;
        archive & boost::serialization::base_object<MutableElement<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rFaces vector of faces associated with the element
     * @param rOrientations vector of orientations of the faces associated with the element
     */
    VertexElement(unsigned index,
                  const std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*>& rFaces,
                  const std::vector<bool>& rOrientations);

    /**
     *
     * Alternative constructor.
     *
     * When constructing a VertexMesh as the Voronoi dual to a Delaunay mesh,
     * each VertexElement is initially constructed without nodes.
     *
     * @param index global index of the element
     */
    VertexElement(unsigned index);

    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     */
    VertexElement(unsigned index,
                  const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Constructor used to specify the element completely. This ensures that
     * the nodes and faces are owned by the element *in a specified order*.
     * See #1076 and #1377 for more details.
     *
     * @param index global index of the element
     * @param rFaces vector of faces associated with the element
     * @param rOrientations vector of orientations of the faces associated with the element
     * @param rNodes vector of Nodes associated with the element
     */
    VertexElement(unsigned index,
                  const std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*>& rFaces,
                  const std::vector<bool>& rOrientations,
                  const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Destructor.
     */
    ~VertexElement();

    /**
     * @return the number of faces owned by this element.
     */
    unsigned GetNumFaces() const;

    /**
     * Add a face to the element.
     *
     * @param pFace a pointer to the new face
     */
    void AddFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace);

    /**
     * @param index the local index of a specified face
     *
     * @return a pointer to the face
     */
    VertexElement<ELEMENT_DIM-1, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * @return whether the face with a given index is oriented clockwise.
     *
     * @param index the index of the face
     */
    bool FaceIsOrientatedClockwise(unsigned index) const;

    // my changes
    bool GetOrientation(unsigned faceLocalIndex)
    {
      return this->mOrientations[faceLocalIndex];
    }

    unsigned GetFaceLocalIndex(unsigned globalIndex) const
    {
      for (unsigned face_local_index = 0; face_local_index < this->GetNumFaces(); face_local_index++)
      {
        if (globalIndex == mFaces[face_local_index]->GetIndex())
        {
          return face_local_index;
          break;
        }
        if (face_local_index == this->GetNumFaces()-1)
        {
          std::cout << std::endl << "ERR: Method VertexElement::GetFaceLocalIndexGetFaceLocalIndex";
          std::cout << std::endl << "Not Reachable!";
        }
      }
      return 0;
    }

    unsigned GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(unsigned startNodeGlobalIndex, unsigned endNodeGlobalIndex)
    {
      for (unsigned face_local_index = 0; face_local_index < this->GetNumFaces(); face_local_index++ )
      {
        if (this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==endNodeGlobalIndex) 
        {
          if (this->mOrientations[face_local_index]==false)
          {
            std::cout << std::endl << "Err in VertexElement::GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex";
            std::cout << std::endl << "Orientation value is wrong!" << std::endl;
            std::cout << "Details: startNodeGlobalIndex=" << startNodeGlobalIndex << " endNodeGlobalIndex=" << endNodeGlobalIndex;
            std::cout <<  std::endl << "ElementIndex: " << this->GetIndex();
            std::cout <<  std::endl << "Nodes: ";
            for(unsigned index =0; index< this->GetNumNodes(); index++)
            {
                std::cout << "NodeLocalIndex"<< index << "_" << this->GetNodeGlobalIndex(index) << ' ';
            }
            std::cout <<  std::endl << "Faces: ";
            for(unsigned index =0; index< this->GetNumFaces(); index++)
            {
                std::cout << "FaceLocalIndex" << index << "_" << this->GetFace(index)->GetIndex();
                std::cout << " orien_" << this->GetOrientation(index);
                std::cout << " fir_node_" << this->GetFace(index)->GetNodeGlobalIndex(0) << " sec_node_" << this->GetFace(index)->GetNodeGlobalIndex(1) << " || ";
            }
            std::cout << std::endl;

          }
          assert(this->mOrientations[face_local_index]==true);
          return face_local_index;
          break;
        }
        else if (this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==endNodeGlobalIndex)
        {
          if (this->mOrientations[face_local_index]==true)
          {
            std::cout << std::endl << "Err in VertexElement::GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex";
            std::cout << std::endl << "Orientation value is wrong!" << std::endl;
            std::cout << "Details: startNodeGlobalIndex=" << startNodeGlobalIndex << " endNodeGlobalIndex=" << endNodeGlobalIndex;
            std::cout <<  std::endl << "ElementIndex: " << this->GetIndex();
            std::cout <<  std::endl << "Nodes: ";
            for(unsigned index =0; index< this->GetNumNodes(); index++)
            {
                std::cout << "NodeLocalIndex"<< index << "_" << this->GetNodeGlobalIndex(index) << ' ';
            }
            std::cout <<  std::endl << "Faces: ";
            for(unsigned index =0; index< this->GetNumFaces(); index++)
            {
                std::cout << "FaceLocalIndex" << index << "_" << this->GetFace(index)->GetIndex();
                std::cout << " orien_" << this->GetOrientation(index);
                std::cout << " fir_node_" << this->GetFace(index)->GetNodeGlobalIndex(0) << " sec_node_" << this->GetFace(index)->GetNodeGlobalIndex(1) << " || ";
            }
            std::cout << std::endl;

          }
          assert(this->mOrientations[face_local_index]==false);
          return face_local_index;
          break;
        }
        if (face_local_index == this->GetNumFaces()-1)
        {
          std::cout << std::endl << "ERR: Method VertexElement::GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex";
          std::cout << std::endl << "Not Reachable!" << std::endl;
          std::cout << "Details: startNodeGlobalIndex=" << startNodeGlobalIndex << " endNodeGlobalIndex=" << endNodeGlobalIndex;
          std::cout <<  std::endl << "ElementIndex: " << this->GetIndex();
          std::cout <<  std::endl << "Nodes: ";
          for(unsigned index =0; index< this->GetNumNodes(); index++)
          {
              std::cout << "NodeLocalIndex"<< index << "_" << this->GetNodeGlobalIndex(index) << ' ';
          }
          std::cout <<  std::endl << "Faces: ";
          for(unsigned index =0; index< this->GetNumFaces(); index++)
          {
              std::cout << "FaceLocalIndex" << index << "_" << this->GetFace(index)->GetIndex();
              std::cout << " orien_" << this->GetOrientation(index);
              std::cout << " fir_node_" << this->GetFace(index)->GetNodeGlobalIndex(0) << " sec_node_" << this->GetFace(index)->GetNodeGlobalIndex(1) << " || ";
          }
          std::cout << std::endl;
          assert(!(face_local_index == this->GetNumFaces()-1));

        }
      }
      return 0;
    }

    bool CheckIfHasThisFace(unsigned startNodeGlobalIndex, unsigned endNodeGlobalIndex)
    {
      for (unsigned face_local_index = 0; face_local_index < this->GetNumFaces(); face_local_index++ )
      {
        if (this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==endNodeGlobalIndex) 
        {
          if (this->mOrientations[face_local_index]==false)
          {
            std::cout << std::endl << "Err in VertexElement::GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex";
            std::cout << std::endl << "Orientation value is wrong!" << std::endl;
            std::cout << "Details: startNodeGlobalIndex=" << startNodeGlobalIndex << " endNodeGlobalIndex=" << endNodeGlobalIndex;
            std::cout <<  std::endl << "ElementIndex: " << this->GetIndex();
            std::cout <<  std::endl << "Nodes: ";
            for(unsigned index =0; index< this->GetNumNodes(); index++)
            {
                std::cout << "NodeLocalIndex"<< index << "_" << this->GetNodeGlobalIndex(index) << ' ';
            }
            std::cout <<  std::endl << "Faces: ";
            for(unsigned index =0; index< this->GetNumFaces(); index++)
            {
                std::cout << "FaceLocalIndex" << index << "_" << this->GetFace(index)->GetIndex();
                std::cout << " orien_" << this->GetOrientation(index);
                std::cout << " fir_node_" << this->GetFace(index)->GetNodeGlobalIndex(0) << " sec_node_" << this->GetFace(index)->GetNodeGlobalIndex(1) << " || ";
            }
            std::cout << std::endl;

          }
          assert(this->mOrientations[face_local_index]==true);
          return true;
          break;
        }
        else if (this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==endNodeGlobalIndex)
        {
          if (this->mOrientations[face_local_index]==true)
          {
            std::cout << std::endl << "Err in VertexElement::GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex";
            std::cout << std::endl << "Orientation value is wrong!" << std::endl;
            std::cout << "Details: startNodeGlobalIndex=" << startNodeGlobalIndex << " endNodeGlobalIndex=" << endNodeGlobalIndex;
            std::cout <<  std::endl << "ElementIndex: " << this->GetIndex();
            std::cout <<  std::endl << "Nodes: ";
            for(unsigned index =0; index< this->GetNumNodes(); index++)
            {
                std::cout << "NodeLocalIndex"<< index << "_" << this->GetNodeGlobalIndex(index) << ' ';
            }
            std::cout <<  std::endl << "Faces: ";
            for(unsigned index =0; index< this->GetNumFaces(); index++)
            {
                std::cout << "FaceLocalIndex" << index << "_" << this->GetFace(index)->GetIndex();
                std::cout << " orien_" << this->GetOrientation(index);
                std::cout << " fir_node_" << this->GetFace(index)->GetNodeGlobalIndex(0) << " sec_node_" << this->GetFace(index)->GetNodeGlobalIndex(1) << " || ";
            }
            std::cout << std::endl;

          }
          assert(this->mOrientations[face_local_index]==false);
          return true;
          break;
        }
      }
      return false;

    }

    // Note: previous node global index should be able to be found in this element!
    // Better not use this method!
    unsigned GetFaceLocalIndexUsingEndNodeGlobalIndex(unsigned endNodeGlobalIndex)
    {
      unsigned endNodeLocalIndex = this->GetNodeLocalIndex(endNodeGlobalIndex);
      unsigned startNodeLocalIndex = (endNodeLocalIndex+this->GetNumNodes()-1)%this->GetNumNodes();
      unsigned startNodeGlobalIndex = this->GetNodeGlobalIndex(startNodeLocalIndex);
      for (unsigned face_local_index = 0; face_local_index < this->GetNumFaces(); face_local_index++ )
      {
        if ((this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==endNodeGlobalIndex) || (this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==endNodeGlobalIndex))
        {
          return face_local_index;
          break;
        }
        if (face_local_index == this->GetNumFaces()-1)
        {
          std::cout << std::endl << "ERR: Method VertexElement::GetFaceLocalIndexUsingEndNodeGlobalIndex";
          std::cout << std::endl << "Not Reachable!";
        }
      }
      return 0;
    }

    unsigned GetFaceLocalIndexUsingStartNodeGlobalIndex(unsigned startNodeGlobalIndex)
    {
      unsigned startNodeLocalIndex = this->GetNodeLocalIndex(startNodeGlobalIndex);
      unsigned endNodeLocalIndex = (startNodeLocalIndex+1)%this->GetNumNodes();
      unsigned endNodeGlobalIndex = this->GetNodeGlobalIndex(endNodeLocalIndex);
      for (unsigned face_local_index = 0; face_local_index < this->GetNumFaces(); face_local_index++ )
      {
        if ((this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==endNodeGlobalIndex) || (this->mFaces[face_local_index]->GetNodeGlobalIndex(1)==startNodeGlobalIndex && this->mFaces[face_local_index]->GetNodeGlobalIndex(0)==endNodeGlobalIndex))
        {
          return face_local_index;
          break;
        }
        if (face_local_index == this->GetNumFaces()-1)
        {
          std::cout << std::endl << "ERR: Method VertexElement::GetFaceLocalIndexUsingStartNodeGlobalIndex";
          std::cout << std::endl << "Not Reachable!";
        }
      }
      return 0;
    }

    bool IfFaceIsCounterClockwise(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace)
    {
      for (unsigned index = 0; index < this->GetNumNodes(); index++)
      {
        unsigned global_index_this_node = this->GetNodeGlobalIndex(index);
        unsigned global_index_next_node = this->GetNodeGlobalIndex((index+1)%this->GetNumNodes());
        if (pFace->GetNodeGlobalIndex(0)==global_index_this_node && pFace->GetNodeGlobalIndex(1)==global_index_next_node)
        {
          return true;
        }
        else if (pFace->GetNodeGlobalIndex(1)==global_index_this_node && pFace->GetNodeGlobalIndex(0)==global_index_next_node)
        {
          return false;
        }
        if (index == this->GetNumNodes()-1)
        {
          std::cout << std::endl << "ERR: Method VertexElement::GetFaceLocalIndexUsingStartNodeGlobalIndex";
          std::cout << std::endl << "Not Reachable!";
        }
      }
      return true;
    }

    void AddFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace, unsigned afterThisFaceLocalIndex)
    {
      this->mFaces.insert(mFaces.begin()+afterThisFaceLocalIndex+1,pFace);
      bool counter_clockwise = this->IfFaceIsCounterClockwise(pFace);
      this->mOrientations.insert(mOrientations.begin()+afterThisFaceLocalIndex+1,counter_clockwise);
    }

    void DeleteFace(unsigned faceLocalIndex)
    {
      this->mFaces.erase(this->mFaces.begin()+faceLocalIndex);
      this->mOrientations.erase(this->mOrientations.begin()+faceLocalIndex);
    }

    // For solving instantialization problem here! this should only be used as a face method!
    double GetUnifiedEdgeMyosinActivty()
    {
      return 0.0;
    }

    double GetUnifiedCellCellAdhesionEnergyParameter()
    {
      return 0.0;
    }

    bool IsDeleted()
    {
      return false;
    }

    void MarkAsDeleted()
    {
    }

    void SetUnifiedEdgeMyosinActivty(double unifiedEdgeMyosinActivty)
    {
    }

    void SetUnifiedCellCellAdhesionEnergyParameter(double unifiedCellCellAdhesionEnergyParameter)
    {
    }

    void ResetFaceValues()
    {
    }

    void ReplaceOneNodeBy(Node<SPACE_DIM>* pNodeA,Node<SPACE_DIM>* pNodeB)
    {
    }
};


//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
class VertexElement<1, SPACE_DIM> : public MutableElement<1,SPACE_DIM>
{
private:
    // My changes
    double mUnifiedEdgeMyosinActivty;

    double mUnifiedCellCellAdhesionEnergyParameter;

    double mUnifiedCellBoundaryAdhesionEnergyParameter;

    bool mIsDeleted;

public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes the nodes owned by the element
     */
    VertexElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * @return the number of faces owned by this element.
     */
    unsigned GetNumFaces() const;

    /**
     * @param index the global index of a specified face
     *
     * @return a pointer to the face
     */
    VertexElement<0, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * @return whether the face with a given index is oriented clockwise.
     *
     * @param index the index of the face
     */
    bool FaceIsOrientatedClockwise(unsigned index) const;

    // My changes
    // void RegisterWithFaces();
    // //return a set of indices of elements containing this element as a face.
    // std::set<unsigned>& rGetContainingElementIndices();

    double GetUnifiedEdgeMyosinActivty()
    {
      return mUnifiedEdgeMyosinActivty;
    }

    double GetUnifiedCellCellAdhesionEnergyParameter()
    {
      return mUnifiedCellCellAdhesionEnergyParameter;
    }

    bool IsDeleted()
    {
      return mIsDeleted;
    }

    void MarkAsDeleted()
    {
      mIsDeleted = true;
    }

    void SetUnifiedEdgeMyosinActivty(double unifiedEdgeMyosinActivty)
    {
      this->mUnifiedEdgeMyosinActivty = unifiedEdgeMyosinActivty;
    }

    void SetUnifiedCellCellAdhesionEnergyParameter(double unifiedCellCellAdhesionEnergyParameter)
    {
      this->mUnifiedCellCellAdhesionEnergyParameter = unifiedCellCellAdhesionEnergyParameter;
    }

    void ResetFaceValues()
    {
      this->mUnifiedEdgeMyosinActivty = 1.0;
      this->mUnifiedCellCellAdhesionEnergyParameter = 1.0;
    }

    void ReplaceOneNodeBy(Node<SPACE_DIM>* pNodeA,Node<SPACE_DIM>* pNodeB)
    {
      if (this->mNodes[0]->GetIndex()==pNodeA->GetIndex())
      {
        this->mNodes[0] = pNodeB;
      }
      else if (this->mNodes[1]->GetIndex()==pNodeA->GetIndex())
      {
        this->mNodes[1] = pNodeB;
      }
      else
      {
          std::cout << std::endl << "ERR: Method VertexElement::ReplaceOneNodeBy";
          std::cout << std::endl << "Not Reachable!";
      }
    }

    // For solving instantialization problem here in 1d element! this should only be used as an element method!
    VertexElement(unsigned index,
              const std::vector<VertexElement<1-1, SPACE_DIM>*>& rFaces,
              const std::vector<bool>& rOrientations,
              const std::vector<Node<SPACE_DIM>*>& rNodes)
      : MutableElement<1, SPACE_DIM>(index, rNodes)// this statement is indispensable!
    {      
    }

    bool GetOrientation(unsigned faceLocalIndex)
    {
      return true;
    }

    unsigned GetFaceLocalIndex(unsigned globalIndex) const
    {
      return 0;
    }

    unsigned GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(unsigned startNodeGlobalIndex, unsigned endNodeGlobalIndex)
    {
      return 0;
    }

    bool CheckIfHasThisFace(unsigned startNodeGlobalIndex, unsigned endNodeGlobalIndex)
    {
      return true;
    }

    unsigned GetFaceLocalIndexUsingEndNodeGlobalIndex(unsigned endNodeGlobalIndex)
    {
      return 0;
    }

    unsigned GetFaceLocalIndexUsingStartNodeGlobalIndex(unsigned startNodeGlobalIndex)
    {
      return 0;
    }

    bool IfFaceIsCounterClockwise(VertexElement<1-1, SPACE_DIM>* pFace)
    {
      return true;
    }

    void AddFace(VertexElement<1-1, SPACE_DIM>* pFace, unsigned afterThisFaceLocalIndex)
    {
    }

    void DeleteFace(unsigned faceLocalIndex)
    {
    }

};

#endif /*VERTEXELEMENT_HPP_*/
