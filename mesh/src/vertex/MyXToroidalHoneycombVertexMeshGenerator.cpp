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

#include "MyXToroidalHoneycombVertexMeshGenerator.hpp"

MyXToroidalHoneycombVertexMeshGenerator::MyXToroidalHoneycombVertexMeshGenerator(unsigned numElementsAcross,
   unsigned numElementsUp,
   double initialArea,
   double cellRearrangementThreshold,
   double t2Threshold)
{
    // numElementsAcross and numElementsUp must be even for toroidal meshes
    assert(numElementsAcross > 1);
    assert(numElementsUp > 1);
    assert(numElementsAcross%2 == 0);///\todo This should be an exception
    assert(numElementsUp%2 == 0);///\todo This should be an exception

    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);

    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;
    std::vector<VertexElement<1,2>*> faces;

    unsigned node_index = 0;
    unsigned node_indices[6];
    unsigned element_index;

    // Create the nodes
    for (unsigned j=0; j<2*numElementsUp+2; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            double x_coord = (((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i) - numElementsAcross/2.0;
            double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3.0);
            x_coord = x_coord*sqrt(initialArea/(sqrt(3)/2));
            y_coord = y_coord*sqrt(initialArea/(sqrt(3)/2));
            bool is_boundary_node = (j==0) || (j==1) || (j == 2*numElementsUp) || (j == (2*numElementsUp+1));

            Node<2>* p_node = new Node<2>(node_index, is_boundary_node , x_coord, y_coord);
            nodes.push_back(p_node);
            node_index++;
        }
    }

    // my changes: create the faces
    unsigned face_index =0;
    for (unsigned j=0; j<2*numElementsUp+1; j++)
    {
        if (j%4==0)
        {
            for (unsigned i=0; i<2*numElementsAcross; i++)
            {
                std::vector<Node<2>*> face_nodes;
                if (i%2==0)
                {
                    face_nodes.push_back(nodes[j*numElementsAcross+i/2+numElementsAcross]);
                    face_nodes.push_back(nodes[j*numElementsAcross+i/2]);
                }
                else
                {
                    if (i == 2*numElementsAcross-1)
                    {
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2]);
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2+numElementsAcross+1-numElementsAcross]);
                    }
                    else
                    {
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2]);
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2+numElementsAcross+1]);
                    }
                }
                VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                faces.push_back(p_face);
                face_index +=1;
            }
        }
        if (j%4==1)
        {
            for (unsigned i=0; i<numElementsAcross; i++)
            {
                std::vector<Node<2>*> face_nodes;
                face_nodes.push_back(nodes[j*numElementsAcross+i]);
                face_nodes.push_back(nodes[j*numElementsAcross+i+numElementsAcross]);                
                VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                faces.push_back(p_face);
                face_index +=1;
            }
        }
        if (j%4==2)
        {
            for (unsigned i=0; i<2*numElementsAcross; i++)
            {
                std::vector<Node<2>*> face_nodes;
                if (i%2==0)
                {
                    face_nodes.push_back(nodes[j*numElementsAcross+i/2]);
                    face_nodes.push_back(nodes[j*numElementsAcross+i/2+numElementsAcross]);
                }
                else
                {
                    if (i == 2*numElementsAcross-1)
                    {
                        face_nodes.push_back(nodes[(j+1)*numElementsAcross+i/2]);
                        face_nodes.push_back(nodes[(j+1)*numElementsAcross+i/2-numElementsAcross+1-numElementsAcross]);
                    }
                    else
                    {
                        face_nodes.push_back(nodes[(j+1)*numElementsAcross+i/2]);
                        face_nodes.push_back(nodes[(j+1)*numElementsAcross+i/2-numElementsAcross+1]);
                    }
                }
                VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                faces.push_back(p_face);
                face_index +=1;
            }
        }
        if (j%4==3)
        {
            for (unsigned i=0; i<numElementsAcross; i++)
            {
                std::vector<Node<2>*> face_nodes;
                face_nodes.push_back(nodes[j*numElementsAcross+i]);
                face_nodes.push_back(nodes[j*numElementsAcross+i+numElementsAcross]);                
                VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                faces.push_back(p_face);
                face_index +=1;
            }
        }

    }

    if (face_index != (3*numElementsUp+2)*numElementsAcross)
        std::cout<< std::endl << "Face building got error!" ;
    else
        std::cout<< std::endl << "Number of faces: " << face_index ;

    assert(face_index == (3*numElementsUp+2)*numElementsAcross);

    // std::cout << std::endl << "Faces information:"
    // for (unsigned index = 0; index < faces.size(); index++)
    // {
    //     std::cout << std::endl << index << ' ' << faces[index]->GetNodeGlobalIndex(0) << ' ' << faces[index]->GetNodeGlobalIndex(1);
    // }

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    for (unsigned j=0; j<numElementsUp; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            element_index = j*numElementsAcross + i;

            node_indices[0] = 2*j*numElementsAcross + i + 1*(j%2==1);
            node_indices[1] = node_indices[0] + numElementsAcross + 1*(j%2==0);
            node_indices[2] = node_indices[0] + 2*numElementsAcross + 1*(j%2==0);
            node_indices[3] = node_indices[0] + 3*numElementsAcross;
            node_indices[4] = node_indices[0] + 2*numElementsAcross - 1*(j%2==1);
            node_indices[5] = node_indices[0] + numElementsAcross - 1*(j%2==1);

            if (i == numElementsAcross-1) // on far right
            {
                node_indices[0] -= numElementsAcross*(j%2==1);
                node_indices[1] -= numElementsAcross;
                node_indices[2] -= numElementsAcross;
                node_indices[3] -= numElementsAcross*(j%2==1);
            }

            std::vector<Node<2>*> element_nodes;
            for (unsigned k=0; k<6; k++)
            {
               element_nodes.push_back(nodes[node_indices[k]]);//node_indices[k] is global index of Node.
            }

            std::vector<VertexElement<1,2>*> element_faces;
            std::vector<bool> element_orientations;
            // element_faces[0] = faces[6*(j/2)*numElementsAcross+2*i+1];//j%2==0
            // element_faces[0] = faces[6*(j/2)*numElementsAcross+2*i+3*numElementsAcross+2];//j%2==1
            element_faces.push_back(faces[6*(j/2)*numElementsAcross+2*i+1+(3*numElementsAcross+1)*(j%2==1)]);
            //element_faces[1] = faces[6*(j/2)*numElementsAcross+i+2*numElementsAcross+1];
            //element_faces[1] = faces[6*(j/2)*numElementsAcross+i+5*numElementsAcross+1];
            element_faces.push_back(faces[6*(j/2)*numElementsAcross+i+2*numElementsAcross+1+3*numElementsAcross*(j%2==1)]);
            // element_faces[2] = faces[6*(j/2)*numElementsAcross+2*i+3*numElementsAcross+1];
            // element_faces[2] = faces[6*(j/2)*numElementsAcross+2*i+6*numElementsAcross+2];
            element_faces.push_back(faces[6*(j/2)*numElementsAcross+2*i+3*numElementsAcross+1+(3*numElementsAcross+1)*(j%2==1)]);
            // element_faces[3] = faces[6*(j/2)*numElementsAcross+2*i+3*numElementsAcross];
            // element_faces[3] = faces[6*(j/2)*numElementsAcross+2*i+6*numElementsAcross+1];
            element_faces.push_back(faces[6*(j/2)*numElementsAcross+2*i+3*numElementsAcross+(3*numElementsAcross+1)*(j%2==1)]);
            // element_faces[4] = faces[6*(j/2)*numElementsAcross+i+2*numElementsAcross];
            // element_faces[4] = faces[6*(j/2)*numElementsAcross+i+5*numElementsAcross];
            element_faces.push_back(faces[6*(j/2)*numElementsAcross+i+2*numElementsAcross+3*numElementsAcross*(j%2==1)]);
            // element_faces[5] = faces[6*(j/2)*numElementsAcross+2*i];
            // element_faces[5] = faces[6*(j/2)*numElementsAcross+2*i+3*numElementsAcross+1];
            element_faces.push_back(faces[6*(j/2)*numElementsAcross+2*i+(3*numElementsAcross+1)*(j%2==1)]);

            // on far right
            if (i== numElementsAcross-1)
            {
                if (j%2==0)
                {
                    element_faces[1]= faces[6*(j/2)*numElementsAcross+2*numElementsAcross];
                }
                else
                {
                    element_faces[0] = faces[6*(j/2)*numElementsAcross+3*numElementsAcross];
                    element_faces[1] = faces[6*(j/2)*numElementsAcross+5*numElementsAcross];
                    element_faces[2] = faces[6*(j/2)*numElementsAcross+6*numElementsAcross];

                }
            }

            element_orientations.push_back(true);
            element_orientations.push_back(true);
            element_orientations.push_back(false);
            element_orientations.push_back(false);
            element_orientations.push_back(false);
            element_orientations.push_back(true);

            VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_faces, element_orientations,element_nodes);
            // VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
            elements.push_back(p_element);
        }
    }

    double center_of_width = 0.0;
    double mesh_width = numElementsAcross*sqrt(initialArea/(sqrt(3)/2));
    mpMesh = new MyXToroidal2dVertexMesh(center_of_width,mesh_width, nodes, elements, cellRearrangementThreshold, t2Threshold);
}

MutableVertexMesh<2,2>* MyXToroidalHoneycombVertexMeshGenerator::GetMesh()
{
    EXCEPTION("A toroidal mesh was created but a normal mesh is being requested.");
    return mpMesh; // Not really
}

MyXToroidal2dVertexMesh* MyXToroidalHoneycombVertexMeshGenerator::GetToroidalMesh()
{
    return (MyXToroidal2dVertexMesh*) mpMesh;
}
