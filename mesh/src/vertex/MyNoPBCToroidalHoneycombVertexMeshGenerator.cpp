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

#include "MyNoPBCToroidalHoneycombVertexMeshGenerator.hpp"

MyNoPBCToroidalHoneycombVertexMeshGenerator::MyNoPBCToroidalHoneycombVertexMeshGenerator(unsigned numElementsAcross,
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
    std::vector<VertexElement<1,2>*>  faces;  // changes made by Chao
    std::vector<VertexElement<2,2>*>  elements;
    
    unsigned node_index = 0;
    unsigned face_index = 0; // changes made by Chao
    unsigned node_indices[6];
    unsigned element_index;

    // Create the nodes
    for (unsigned j=0; j<2*numElementsUp+2; j++)
    {
        unsigned numLineNode = ((j==0||j==2*numElementsUp+1)? numElementsAcross:numElementsAcross+1);
        // std::cout << "line " << j << " node number " << numLineNode << std::endl;

        for (unsigned i=0; i<numLineNode; i++)
        {
            double x_coord = (((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i) - numElementsAcross/2.0 + 1*(j==2*numElementsUp+1);
            double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3.0);
            x_coord = x_coord*sqrt(initialArea/(sqrt(3)/2));
            y_coord = y_coord*sqrt(initialArea/(sqrt(3)/2));
            bool is_boundary_node = (j==0) || (j==1) || (j==2*numElementsUp) || (j==(2*numElementsUp+1)) || (i==0) || (i==(numLineNode-1));

            // if (is_boundary_node)
            // {
            //     std::cout << "boundary node " << node_index << std::endl;
            // }

            Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord, y_coord);
            nodes.push_back(p_node);
            node_index++;
        }
    }

    // Create the faces
    for (unsigned j=0; j<2*numElementsUp+1; j++)
    {
        if (j%2==0)
        {
            if (j==0)
            {
                for (unsigned i=0; i<2*numElementsAcross; i++)
                {
                    std::vector<Node<2>*> face_nodes;
                    if (i%2==0)
                    {
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2+numElementsAcross]);
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2]);

                        // std::cout << "face " << face_index << " nodes " << j*numElementsAcross+i/2+numElementsAcross << " " << j*numElementsAcross+i/2 << std::endl;
                    }
                    else
                    {
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2]);
                        face_nodes.push_back(nodes[j*numElementsAcross+i/2+numElementsAcross+1]);

                        // std::cout << "face " << face_index << " nodes " << j*numElementsAcross+i/2 << " " << j*numElementsAcross+i/2+numElementsAcross+1 << std::endl;
                    }
                    VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                    faces.push_back(p_face);

                    face_index ++;
                }
            }
            else if (j==2*numElementsUp)
            {
                for (unsigned i=0; i<2*numElementsAcross; i++)
                {
                    std::vector<Node<2>*> face_nodes;
                    if (i%2==0)
                    {
                        face_nodes.push_back(nodes[j*(numElementsAcross+1)-1+i/2]);
                        face_nodes.push_back(nodes[(j+1)*(numElementsAcross+1)-1+i/2]);

                        // std::cout << "face " << face_index << " nodes " << j*(numElementsAcross+1)-1+i/2 << " " << (j+1)*(numElementsAcross+1)-1+i/2 << std::endl;
                    }
                    else
                    {
                        face_nodes.push_back(nodes[(j+1)*(numElementsAcross+1)-1+i/2]);
                        face_nodes.push_back(nodes[j*(numElementsAcross+1)-1+(i+1)/2]);

                        // std::cout << "face " << face_index << " nodes " << (j+1)*(numElementsAcross+1)-1+i/2 << " " << j*(numElementsAcross+1)-1+(i+1)/2 << std::endl;
                    }
                    VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                    faces.push_back(p_face);
                    face_index ++;
                }
            }
            else if (j%4==0)
            {
                for (unsigned i=0; i<2*numElementsAcross+1; i++)
                {
                    std::vector<Node<2>*> face_nodes;
                    if (i%2==0)
                    {
                        face_nodes.push_back(nodes[(j+1)*(numElementsAcross+1)-1+i/2]);
                        face_nodes.push_back(nodes[j*(numElementsAcross+1)-1+i/2]);

                        // std::cout << "face " << face_index << " nodes " << (j+1)*(numElementsAcross+1)-1+i/2 << " " << j*(numElementsAcross+1)-1+i/2 << std::endl;
                    }
                    else
                    {
                        face_nodes.push_back(nodes[j*(numElementsAcross+1)-1+i/2]);
                        face_nodes.push_back(nodes[(j+1)*(numElementsAcross+1)-1+(i+1)/2]);

                        // std::cout << "face " << face_index << " nodes " << j*(numElementsAcross+1)-1+i/2 << " " << (j+1)*(numElementsAcross+1)-1+(i+1)/2 << std::endl;
                    }
                    VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                    faces.push_back(p_face);
                    face_index ++;
                }
            }
            else
            {
                for (unsigned i=0; i<2*numElementsAcross+1; i++)
                {
                    std::vector<Node<2>*> face_nodes;
                    if (i%2==0)
                    {
                        face_nodes.push_back(nodes[j*(numElementsAcross+1)-1+i/2]);
                        face_nodes.push_back(nodes[(j+1)*(numElementsAcross+1)-1+i/2]);

                        // std::cout << "face " << face_index << " nodes " << j*(numElementsAcross+1)-1+i/2 << " " << (j+1)*(numElementsAcross+1)-1+i/2 << std::endl;
                    }
                    else
                    {
                        face_nodes.push_back(nodes[(j+1)*(numElementsAcross+1)-1+i/2]);
                        face_nodes.push_back(nodes[j*(numElementsAcross+1)-1+(i+1)/2]);

                        // std::cout << "face " << face_index << " nodes " << (j+1)*(numElementsAcross+1)-1+i/2 << " " << j*(numElementsAcross+1)-1+(i+1)/2 << std::endl;
                    }
                    VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                    faces.push_back(p_face);
                    face_index ++;
                }
            }
        }

        if (j%2==1)
        {
            for (unsigned i=0; i<numElementsAcross+1; i++)
            {
                std::vector<Node<2>*> face_nodes;
                face_nodes.push_back(nodes[j*(numElementsAcross+1)-1+i]);
                face_nodes.push_back(nodes[(j+1)*(numElementsAcross+1)-1+i]); 

                // std::cout << "face " << face_index << "nodes " << j*(numElementsAcross+1)-1+i << " " << (j+1)*(numElementsAcross+1)-1+i << std::endl;

                VertexElement<1,2>* p_face = new VertexElement<1,2>(face_index, face_nodes);
                faces.push_back(p_face);
                face_index ++;
            }
        }

    }

    if (face_index != (3*numElementsAcross+2)*numElementsUp+2*numElementsAcross-1)
        std::cout<< std::endl << "Face building got error in MyToroidalHoneycombVertexMeshGenerator!" << std::endl;
    // else
    //     std::cout<< std::endl << "Successfully built " << face_index << " faces in MyToroidalHoneycombVertexMeshGenerator." << std::endl;

    assert(face_index == (3*numElementsAcross+2)*numElementsUp+2*numElementsAcross-1);

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    for (unsigned j=0; j<numElementsUp; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            element_index = j*numElementsAcross + i;

            if(j%2==0)
            {
                node_indices[0] = j*2*(numElementsAcross+1)-1*(j==0? 0:1) + i ;
                node_indices[1] = (2*j+1)*(numElementsAcross+1) + i ;
                node_indices[2] = 2*(j+1)*(numElementsAcross+1) + i ;
                node_indices[3] = (2*(j+1)+1)*(numElementsAcross+1) -1 + i ;
                node_indices[4] = 2*(j+1)*(numElementsAcross+1) -1 + i ;
                node_indices[5] = (2*(j+1)-1)*(numElementsAcross+1) -1 + i ;
            }
            if (j%2==1)
            {
                node_indices[0] = 2*j*(numElementsAcross+1) + i ;
                node_indices[1] = (2*j+1)*(numElementsAcross+1) + i ;
                node_indices[2] = 2*(j+1)*(numElementsAcross+1) + i ;
                node_indices[3] = (2*(j+1)+1)*(numElementsAcross+1) -1*(j==(numElementsUp-1)?1:0) + i ;
                node_indices[4] = 2*(j+1)*(numElementsAcross+1) -1 + i ;
                node_indices[5] = (2*(j+1)-1)*(numElementsAcross+1) -1 + i ;
            }

            std::vector<Node<2>*> element_nodes;
            for (unsigned k=0; k<6; k++)
            {
                element_nodes.push_back(nodes[node_indices[k]]);//node_indices[k] is global index of Node.
            }

            std::vector<VertexElement<1,2>*> element_faces;
            element_faces.push_back(faces[j*(3*numElementsAcross+2)-1 + ((j==0||j%2==1)? 2:1) + 2*i]);
            element_faces.push_back(faces[j*(3*numElementsAcross+2)-1 + 2*numElementsAcross+1 +1 + i]);
            element_faces.push_back(faces[(j+1)*(3*numElementsAcross+2)-1 + ((j==(numElementsUp-1)||j%2==0)? 1:2) + 2*i]);
            element_faces.push_back(faces[(j+1)*(3*numElementsAcross+2)-1 + ((j==(numElementsUp-1)||j%2==0)? 0:1) + 2*i]);
            element_faces.push_back(faces[j*(3*numElementsAcross+2) + 2*numElementsAcross+1 -1 + i]);
            element_faces.push_back(faces[j*(3*numElementsAcross+2)-1 + ((j==0||j%2==1)? 1:0) + 2*i]);

            std::vector<bool> element_orientations;
            element_orientations.push_back(true);
            element_orientations.push_back(true);
            element_orientations.push_back(false);
            element_orientations.push_back(false);
            element_orientations.push_back(false);
            element_orientations.push_back(true);

            std::vector<VertexElement<1,2>*> element_stress_fibers;
            VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_faces, element_orientations, element_nodes, element_stress_fibers);

            // VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_faces, element_orientations, element_nodes);
            // VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
            elements.push_back(p_element);

            // std::cout << "element " << element_index << " " 
            // << "faces " << j*(3*numElementsAcross+2)-1 + ((j==0||j%2==1)? 2:1) + 2*i << " " 
            // << j*(3*numElementsAcross+2)-1 + 2*numElementsAcross+1 +1 + i << " " 
            // << (j+1)*(3*numElementsAcross+2)-1 + ((j==(numElementsUp-1)||j%2==0)? 1:2) + 2*i << " "
            // << (j+1)*(3*numElementsAcross+2)-1 + ((j==(numElementsUp-1)||j%2==0)? 0:1) + 2*i << " " 
            // << j*(3*numElementsAcross+2) + 2*numElementsAcross+1 -1 + i << " " 
            // << j*(3*numElementsAcross+2)-1 + ((j==0||j%2==1)? 1:0) + 2*i << " "
            // << "nodes " << node_indices[0] << " " << node_indices[1] << " " << node_indices[2] << " "
            // << node_indices[3] << " " << node_indices[4] << " " << node_indices[5] << std::endl;

        }
    }

    double center_of_width = 0.0;
    double mesh_width = numElementsAcross*sqrt(initialArea/(sqrt(3)/2));
    mpMesh = new MyNoPBCToroidal2dVertexMesh(center_of_width, mesh_width, nodes, elements, cellRearrangementThreshold, t2Threshold);
}

MutableVertexMesh<2,2>* MyNoPBCToroidalHoneycombVertexMeshGenerator::GetMesh()
{
    EXCEPTION("A toroidal mesh was created but a normal mesh is being requested.");
    return mpMesh; // Not really
}

MyNoPBCToroidal2dVertexMesh* MyNoPBCToroidalHoneycombVertexMeshGenerator::GetToroidalMesh()
{
    return (MyNoPBCToroidal2dVertexMesh*) mpMesh;
}
