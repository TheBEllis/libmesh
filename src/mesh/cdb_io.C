// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/cdb_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/int_range.h"
#include "libmesh/utility.h" // map_find
#include "libmesh/enum_to_string.h"

// C++ includes
#include <fstream>
#include <set>
#include <cstring> // std::memcpy
#include <numeric>
#include <unordered_map>
#include <cstddef>

namespace libMesh
{

// Initialize the static data member
CDBIO::ElementMaps CDBIO::_element_maps = CDBIO::build_element_maps();

// Initialize the static data member
CDBIO::Ansys2LibmeshElementMaps CDBIO::_ansys2libmesh_element_maps = CDBIO::build_ansys_element_maps();

// Definition of the static function which constructs the ElementMaps object.
CDBIO::ElementMaps CDBIO::build_ansys_element_maps()
{
    std::map<> element_map;

    element_map[226] = HEX20; 
    element_map[227] = TET10;

    return element_map;
}

// Definition of the static function which constructs the ElementMaps object.
CDBIO::ElementMaps CDBIO::build_element_maps()
{
  // Object to be filled up
  ElementMaps em;

  // POINT (import only)
  em.in.emplace(15, ElementDefinition(NODEELEM, 15, 0, 1));

  // Add elements with trivial node mappings
  em.add_def(ElementDefinition(EDGE2, 1, 1, 2));
  em.add_def(ElementDefinition(EDGE3, 8, 1, 3));
  em.add_def(ElementDefinition(TRI3, 2, 2, 3));
  em.add_def(ElementDefinition(TRI6, 9, 2, 6));
  em.add_def(ElementDefinition(QUAD4, 3, 2, 4));
  em.add_def(ElementDefinition(QUAD8, 16, 2, 8));
  em.add_def(ElementDefinition(QUAD9, 10, 2, 9));
  em.add_def(ElementDefinition(HEX8, 5, 3, 8));
  em.add_def(ElementDefinition(TET4, 4, 3, 4));
  em.add_def(ElementDefinition(PRISM6, 6, 3, 6));
  em.add_def(ElementDefinition(PYRAMID5, 7, 3, 5));

  // Add elements with non-trivial node mappings

  // HEX20
  {
    ElementDefinition eledef(HEX20, 17, 3, 20);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,15,16,19,17,18};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // HEX27
  {
    ElementDefinition eledef(HEX27, 12, 3, 27);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,
                                  15,16,19,17,18,20,21,24,22,23,25,26};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // TET10
  {
    ElementDefinition eledef(TET10, 11, 3, 10);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,7,9,8};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // PRISM15
  {
    ElementDefinition eledef(PRISM15, 18, 3, 15);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,12,14,13};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // PRISM18
  {
    ElementDefinition eledef(PRISM18, 13, 3, 18);
    const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,12,14,13,15,17,16};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  return em;
}


std::vector<std::string> tokenize(const std::string to_split) 
{
  std::stringstream ss(to_split);
  std::vector<std::string> result;

  while( ss.good() )
  {
    std::string substr;
    std::getline( ss, substr, ',' );
    result.push_back( substr );
  }
  return result;
}

void CDBIO::read (const std::string & name)
{
  std::ifstream in (name.c_str());
  this->read_mesh (in);
}

void CDBIO::read_mesh(std::istream & in)
{
  // This is a serial-only process for now;
  // the Mesh should be read on processor 0 and
  // broadcast later
  libmesh_assert_equal_to (MeshOutput<MeshBase>::mesh().processor_id(), 0);

  libmesh_assert(in.good());

  // clear any data in the mesh
  MeshBase & mesh = MeshInput<MeshBase>::mesh();
  mesh.clear();

  // For reading the file line by line
  std::string s;

  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(in, s);

      if (in)
        {
          // Process s...
          if (s.find("NBLOCK,6,SOLID") == static_cast<std::string::size_type>(0)) 
            {
            // we have found the NBLOCK keyword, following this is some garbage 
            // and following that are the node ids and positions

            // next line isn't really needed but could be used to set the format.
            std::getline(in, s);

            // regex to match the nodes
            std::regex regexp("^ \d+\s+\d+\s+\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+$"); 
            // s should point to a string
            std::getline(in,s);
            // now we can start reading the nodes
            while (std::regex_match(s,regexp) ) 
            {
                int id;
                int dump1,dump2;
                Real x, y, z;
                // read the data 
                in >> id >> dump1 >> dump2 >> x >> y >> z;

                // make the node
                mesh.add_point (Point(x, y, z), idx);

              
                // these two should likely be class member varuiables (maps) too
                // maybe only need on of them and search for the key given the 
                // value
                // add cdb node id to map
                nodeidx2id[idx] = id;
                nodeid2idx[id] = idx;

                // increment the node index
                idx++;

                // update string
                std::getline(in, s);
            }           
          } else if (s.find("ET,") == static_cast<std::string::size_type>(0))  {
            // we have found the Element Type keyword - there may be more than one
            std::vector<std::string> tokens = tokenize(s);
            // the first token is ET
            int element_id = std::to_int(tokens[1]);
            int element_type = std::to_int(tokens[2]);

            _element_types_present[element_id] = element_type;
            // update string
            std::getline(in, s);
           } else if (s.find("TYPE,") == static_cast<std::string::size_type>(0)) {
            // we have found the TYPE keyword this is followed by the 
            // EBLOCK keyword
            std::getline(in,s);
            // the next line is a misleading fortran style formatting statement
            std::getline(in,s);
            // the next line is the real start of the data
            while (true) 
            {
              // if we read a -1, then that indicates the end of the element data
              if(s.find("-1") == static_cast<std::string::size_type(0))
                break;

              int garbage;

              // trash the first 8 items
              for (int i = 0 ; i < 8 ; i++) 
              {
                s >> garbage;
              }  
              // the next item is the number of nodes the element has
              int nel;
              s >> nel;

              // the next item is garbage
              s >> garbage;

              // the next item is the CDB element ID
              int elid;
              s >> elid;

              std::array<int,nel> nodes; 
              // we can now read 8 nodes
              for ( int i = 0 ; i < 8 ; i++)
              {
                s >> nodes[i];
              } 
              // read another line
              std::getline(in,s);

              // now we can read the remaining nel-8 nodes
              for ( int i = 0 ; i < nel-8 ; i++) 
              {
                s >> nodes[i+8];
              }

              // now we can make a new element
              Elem* elem = mesh.add_elem(Elem::build_with_id(eletype.type, iel++));
              
              // iel should likely be a class member variable which can be incremented

              // loop over the nodes adding them to the element
              for ( int i = 0 ; i < nel ; i++)
              {
                // setup the element with the correct node ids;
                // i have no idea if the node ordering is the same for 
                // ansys elements and libmesh elements? 
                // im sure life will be cruel
                elem->set_node(i) = mesh.node_ptr(nodeid2idx[nodes[i]]);
              }

              // read a new line
              std::getline(in,s);
            }    

           // we should now have a collection of elements in the map
           // and these can now belong to a block 
           std::getline(in,s);
           // tokenize the string
           std::vector<std::string> tokens = tokenize(s);
           // the new block name is given by the 2nd token
           std::string block_name = tokens[1];
           // now add all the elements that we just created into that block
           

           }

        }


      // If !in, check to see if EOF was set.  If so, break out
      // of while loop.
      if (in.eof())
        break;

      // If !in and !in.eof(), stream is in a bad state!
      libmesh_error_msg("Stream is bad! Perhaps the file does not exist?");

    } // while true
}

} // namespace libMesh
