/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef POTTSMESHGENERATOR_HPP_
#define POTTSMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "PottsMesh.hpp"

/**
 *  Generator of regular potts meshes, used as starting points for many simulations.
 *
 *  This class takes in options such as width, height,
 *
 *
 *  NOTE: the user should delete the mesh after use to manage memory.
 */
class PottsMeshGenerator
{
protected:

    /** A pointer to the mesh this class creates */
    PottsMesh* mpMesh;

public:

    /**
     * Constructor.
     *
     * @param numNodesAcross  The number of columns of nodes in the mesh
     * @param numNodesUp  The number of rows of nodes in the mesh
     * @param numElementsAcross  The number of columns of elements in the mesh
     * @param numElementsUp  The number of rows of elements in the mesh
     * @param elementWidth  The number of rows of nodes in each element
     * @param elementHeight  The number of rows of nodes in each element
     */
    PottsMeshGenerator(unsigned numNodesAcross, unsigned numNodesUp,
                       unsigned numElementsAcross, unsigned numElementsUp,
                       unsigned elementWidth, unsigned elementHeight);

    /**
     * Null constructor for derived classes to call...
     */
    PottsMeshGenerator()
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer
     */
    virtual ~PottsMeshGenerator();

    /**
     * @return a potts mesh based on a 2D plane.
     */
    virtual PottsMesh* GetMesh();
};

#endif /*POTTSMESHGENERATOR_HPP_*/
