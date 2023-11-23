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
#include <regex>

// Testing includes 
#include"unistd.h"

namespace libMesh
{

// Initialize the static data member
CDBIO::CDBMaps CDBIO::build_element_maps()
{
  // Object to be filled up
  CDBMaps em;

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
    ElementDefinition eledef(HEX20, 226, 3, 20);
    
    //  const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
    // const unsigned int nodes[] = {0,3,2,1,7,4,5,6,11,8,9,10,19,16,17,18,15,12,13,14};
    // const unsigned int nodes[] = {3,0,1,2,7,4,5,6,11,8,9,10,19,16,17,18,15,12,13,14};
    // const unsigned int nodes[] = {3,2,6,7,0,1,5,4,10,18,14,19,11,9,13,15,8,17,12,16};
    const unsigned int nodes[] = {7,4,0,3,6,5,1,2,15,16,11,19,14,12,8,10,13,17,9,18};
    std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
    em.add_def(eledef);
  }

  // // HEX27
  // {
  //   ElementDefinition eledef(HEX27, 12, 3, 27);
  //   const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,11,12,9,13,10,14,
  //                                 15,16,19,17,1 8,20,21,24,22,23,25,26};
  //   std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
  //   em.add_def(eledef);
  // }

  // // TET10
  // {
  //   ElementDefinition eledef(TET10, 11, 3, 10);
  //   const unsigned int nodes[] = {0,1,2,3,4,5,6,7,9,8};
  //   std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
  //   em.add_def(eledef);
  // }

  // // PRISM15
  // {
  //   ElementDefinition eledef(PRISM15, 18, 3, 15);
  //   const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,12,14,13};
  //   std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
  //   em.add_def(eledef);
  // }

  // // PRISM18
  // {
  //   ElementDefinition eledef(PRISM18, 13, 3, 18);
  //   const unsigned int nodes[] = {0,1,2,3,4,5,6,8,9,7,10,11,12,14,13,15,17,16};
  //   std::vector<unsigned int>(nodes, nodes+eledef.nnodes).swap(eledef.nodes); // swap trick
  //   em.add_def(eledef);
  // }

  return em;
}

CDBIO::CDBMaps CDBIO::_cdb_maps = CDBIO::build_element_maps();


CDBIO::CDBIO (const MeshBase & mesh) :
  MeshOutput<MeshBase>(mesh)
{;}

CDBIO::CDBIO (MeshBase & mesh) :
  MeshInput<MeshBase>  (mesh),
  MeshOutput<MeshBase> (mesh)
{;}


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
  // Need node id initialised outside of loop
  int id = 0;
  int iel = 0;
  int block_num = 1;
  int nodeset_id = 1;

  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(in, s);

      // If the stream is working...
      if (in)
      {
        // Process s...
        if (s.find("NBLOCK,6,SOLID") == static_cast<std::string::size_type>(0)) 
        {
          // we have found the NBLOCK keyword, following this is some garbage 
          // and following that are the node ids and positions
          std::getline(in, s);  

          // regex to match the nodes
          std::regex regexp("^\\s+\\d+\\s+\\d+\\s+\\d+\\s+[-+]?\\d+\\.\\d+(E{0,1}[-\\+]\\d+)?\\s+[-+]?\\d+\\.\\d+(E{0,1}[-\\+]\\d+)?\\s+[-+]?\\d+\\.\\d+(E{0,1}[-\\+]\\d+)?\\s*?\\r*?$");
          // s should point to a string
          // now we can start reading the nodes

                

          // Set s to node line 1, so regex matches correctly
          std::getline(in, s);  

          // Set a position at start of line 2
          std::streampos sol = in.tellg();

          // Set istream to start of line 1, so we can correctly read data on first iteration
          in.seekg(sol);

          while(std::regex_match(s, regexp)) 
          {
            // Create stream from current string
            std::stringstream line_stream(s);
            int ansys_id;
            double x, y, z;
            int dump1, dump2;

            // Input data from stringstream
            line_stream >> ansys_id >> dump1 >> dump2 >> x >> y >> z;

            // make the node
            mesh.add_point(Point(x, y, z), id);

            // these two should likely be class member varuiables (maps) too
            // maybe only need on of them and search for the key given the 
            // value
            // add cdb node id to map
            ansys2LibmeshIdMap[ansys_id] = id;
            // nodeid2idx[id] = idx;

            // increment the node index
            id++;

            // Get next line
            std::getline(in, s);
          }           
        } 
        
        else if (s.find("ET,") == static_cast<std::string::size_type>(0))  
        {
          // we have found the Element Type keyword - there may be more than one
          std::vector<std::string> tokens = tokenize(s);
          // the first token is ET
          int element_id = std::stoi(tokens[1]);
          int element_type = std::stoi(tokens[2]);

          // _element_types_present[element_id] = element_type;
          // update string
          std::getline(in, s);
        } 
        
        else if (s.find("TYPE,") == static_cast<std::string::size_type>(0)) 
        {
          // we have found the TYPE keyword this is followed by the 
          // EBLOCK keyword
          std::getline(in,s);
          // the next line is a misleading fortran style formatting statement
          std::getline(in,s);
          // the next line is the real start of the data
          std::getline(in,s);

          while(true)
          {
            // if we read a -1, then that indicates the end of the element data
            if(s.find("-1") != std::string::npos)
            {
              break;
            }
            
            std::stringstream strm(s);
            int garbage;
            // trash the first 8 items
            for (int i = 0 ; i < 8 ; i++) 
            {
              strm >> garbage;
            }  
            // the next item is the number of nodes the element has
            int n_cdb_nodes;
            
            strm >> n_cdb_nodes;

            // the next item is garbage
            strm >> garbage;

            // the next item is the CDB element ID
            int elid;

            strm >> elid;

            std::vector<int> nodes;
            nodes.reserve(n_cdb_nodes); 
            // we can now read the first 8 nodes
            
            for ( int i = 0 ; i < 8 ; i++)
            {
              strm >> nodes[i];
            } 
            
            // read another line to get the next line of node ID's
            std::getline(in,s);

            std::stringstream strm_node_line_two(s);         
            // now we can read the remaining n_cdb_nodes-8 nodes
            for ( int i = 0 ; i < n_cdb_nodes-8 ; i++) 
            {
              strm_node_line_two >> nodes[i+8];
            }

            // Get the element_map corresponding to the ansys element type we have identified
            std::vector<unsigned int> elem_map = _cdb_maps.in[226].nodes;

            // now we can make a new element
            Elem* elem = mesh.add_elem(Elem::build_with_id(_cdb_maps.in[226].type, iel++));
            
            // iel should likely be a class member variable which can be incremented

            // loop over the nodes adding them to the element
            for ( int i = 0 ; i < n_cdb_nodes ; i++)
            {
              // setup the element with the correct node ids;
              // i have no idea if the node ordering is the same for 
              // ansys elements and libmesh elements? 
              // im sure life will be cruel
              elem->set_node(i) = mesh.node_ptr(ansys2LibmeshIdMap[nodes[elem_map[i]]]);
              elem->subdomain_id() = block_num;
            }
            mesh.add_elem(elem);
            // read a new line
            std::getline(in, s);
          }    
          // we should now have a collection of elements in the map
          // and these can now belong to a block 
          std::getline(in, s);
          // tokenize the string
          std::vector<std::string> tokens = tokenize(s);
          // the new block name is given by the 2nd token
          std::string block_name = tokens[1];
          boost::algorithm::trim(block_name);
          // now correctly name the block
          mesh.subdomain_name(block_num) = block_name;
          // Increment block id
          block_num++;
        } 
        
        else if (s.find("CMBLOCK,") != std::string::npos) 
        {
          // Semi global ints to track how many nodes have been added to the "nodes" vector
          int nodeset_counter = 0;

          // This conditional block gets the nodesets from the ansys model
          std::cout << "CMBLOCK" << std::endl;
          
          // Tokenize the string, as it contains the number of nodes in the nodeset
          std::vector<std::string> tokens = tokenize(s);

          // Make it long just in case there's a massive nodeset
          long int n_nodeset_nodes = std::stoi(tokens[3]);

          // Second token is the name of the nodeset
          std::string nodeset_name = tokens[1];
          boost::algorithm::trim(nodeset_name);
          

          // Create a vector to store nodes
          std::vector<int> nodes;
          
          // Next line is more ansys garbage
          std::getline(in, s);
          // The next line is also garbage
          std::getline(in, s);

          // The last nodeset ends the file, there is no nice ending character, just eof
          std::regex nodeset_regex("^(\\s+[-]?\\d+)+(\\s+)?\\r?$");
          // The CMBLOCK doesn't seem to have a nice "-1" ending like the element blocks.
          // Henceforth, we need to store a streampos object that stores the position of the previous line,
          //  in case we go too far
          std::streampos prev_line;
          
          while(std::regex_match(s, nodeset_regex))
          {
            // We are now not getting garbage
            // Create a stream of the current line from a string
            std::stringstream current_line(s);

            //  Nodesets can specify ranges of nodes.
            //  They do this using a minus symbol, for example, the line
            //  "1  2  4  -6",
            //  indicates that nodes 1 2 4 5 and 6 are all part of the nodeset.
            // This horrible for loop helps deal with that.
            for(long int i = 0; i < 8; i++)
            {
              long int node = 0;
              current_line >> node;
              // If what we just input was not an int, then cry
              if(current_line.fail())
              {
                break;
              }

              nodes.push_back(node);

              // If the value we just added was negative, then we need to add a range
              if(nodes.back() < 0)
              {
                // Set node to it's correct positive value
                nodes.back() = std::abs(nodes.back());
                
                // Get the element we just added, which is our upper bound
                int upper_bound = nodes.back();

                // Get the second to last element added, which is our lower bound
                int lower_bound = *(&(nodes.back()) - 1);

                // Loop over all the values in the specified range
                for(int range = lower_bound; range < upper_bound; range++)
                {
                  nodes.push_back(range);
                }
              }
            }
            // Store position of the line we just got data from
            prev_line = in.tellg();

            // Get the next line
            std::getline(in, s);
          }

          
          for(auto& node: nodes)
          {
            mesh.get_boundary_info().add_node(ansys2LibmeshIdMap[node], nodeset_id);
          }
          
          // Set the correct nodeset name and increment the nodeset_id for next nodeset
          mesh.get_boundary_info().nodeset_name(nodeset_id++) = nodeset_name;

          // If we have exited the while loop, then we have found the next CMBLOCK line
          //  Because we do an std::nextline at the start of the loop, we need to go back
          //   a line so it is not skipped
          in.seekg(prev_line);
        }
      }
      else if (in.eof())
      {
        break;
      }
      else
      {
        // If !in and !in.eof(), stream is in a bad state!
        libmesh_error_msg("Stream is bad! Perhaps the file does not exist?");
      }
      

    } // while true
}

} // namespace libMesh
