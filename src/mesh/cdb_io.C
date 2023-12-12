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

  // Add elements with non-trivial node mappings

  // ANSYS SOLID226
  {
    AnsysElementDefinition SOLID226(226, 3);
    
    // std::vector<unsigned int> hex_ordering = {7,4,0,3,6,5,1,2,15,16,11,19,14,12,8,10,13,17,9,18};
    std::vector<unsigned int> hex_ordering = {3,0,1,2,7,4,5,6,11,8,9,10,19,16,17,18,15,12,13,14};
    std::vector<unsigned int> tet_ordering = {2,0,1,3,6,4,5,9,7,8};
    std::vector<unsigned int> prism_ordering = {2,0,1,5,3,4,8,6,7,14,12,13,11,9,10};
    std::vector<unsigned int> pyramid_ordering = {3,0,1,2,4,8,5,6,7,12,9,10,11};

    SOLID226.add_elem_sub_mapping(hex_ordering, libMesh::ElemType::HEX20, "HEX20");
    SOLID226.add_elem_sub_mapping(tet_ordering, libMesh::ElemType::TET10, "TET10");
    SOLID226.add_elem_sub_mapping(prism_ordering, libMesh::ElemType::PRISM15, "PRISM15");
    SOLID226.add_elem_sub_mapping(pyramid_ordering, libMesh::ElemType::PYRAMID13, "PYR13");
    em.add_def(SOLID226);
  }

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
  int ansys_element_type;

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
          
          // Now we can start reading the nodes
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

            // Add cdb node id to map
            ansys_to_libmesh_node_id_map[ansys_id] = id;

            // increment the node index
            id++;

            // Get next line
            std::getline(in, s);
          }           
        } 
        
        else if (s.find("ET,") == static_cast<std::string::size_type>(0))  
        {
          // We have found the Element Type keyword - there may be more than one
          std::vector<std::string> tokens = tokenize(s);

          // The second token is the ansys element type i.e. SOLID226
          ansys_element_type = std::stoi(tokens[2]);

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

          // Variable outside of while loop scope, to store last iterations element node count
          int prev_elem_nodes = -1;

          // Vector to store element types that occur in the ANSYS block
          std::vector<std::string> block_elem_types;

          // Vector to store block id's that will be used for this ANSYS block
          std::vector<unsigned int> block_nums;
          block_nums.push_back(block_num);
          
          while(true)
          {
            // String to keep track of elem type
            std::string block_elem_type;

            // If we read a -1, then that indicates the end of the element data
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
            int ansys_elem_id;
            strm >> ansys_elem_id;

            // Nodes vector to store the ansys to libmesh node mapping for this element
            std::vector<int> nodes;
                
            // We can now read the first 8 nodes
            for ( int i = 0 ; i < 8 ; i++)
            {
              int temp;
              strm >> temp;
              nodes.push_back(temp);
            } 
            
            // If there are more than 8 nodes, we need to read another line
            if(n_cdb_nodes > 8)
            {
              // read another line to get the next line of node ID's
              std::getline(in,s);
              std::stringstream strm_node_line_two(s);         
              // now we can read the remaining n_cdb_nodes - 8 nodes
              for ( int i = 0 ; i < n_cdb_nodes - 8 ; i++) 
              {
                int temp;
                strm_node_line_two >> temp;
                nodes.push_back(temp);
              }
            }

            // Remove potential duplicate entries from our nodes vector
            {
              std::set<unsigned int> set;            
              for (auto iter = nodes.begin(); iter != nodes.end();) {
                if (set.find(*iter) == set.end()) {
                  set.insert(*iter);
                  iter++;
                }
                else {
                  iter = nodes.erase(iter);
                }
              }
            }

            /**
             * Because of ANSYS duplicate node numbering, n_cdb_nodes can be a lie,
             * hence use the size of the nodes vector with duplicates removed
             * to determine number of nodes in an element.
            */
            int n_elem_nodes = nodes.size();
    
            // Get the element_map corresponding to the ansys element type we have identified
            std::vector<unsigned int> elem_map;
            libMesh::ElemType elem_type;
            std::string elem_type_string;

            elem_map = _cdb_maps[ansys_element_type].ansys_node_ordering_map[n_elem_nodes];
            elem_type = _cdb_maps[ansys_element_type].ansys_to_libmesh_elem_type_map[n_elem_nodes];
            block_elem_type = _cdb_maps[ansys_element_type].ansys_to_libmesh_elem_type_string_map[n_elem_nodes];
          
            // Ansys can change element type in the middle of a block. If this happens we need to increment the block number,
            // because exodus can only deal with one element type per block
            if(prev_elem_nodes == -1)
            {
              block_elem_types.push_back(block_elem_type);
            }

            // Ansys can change element type in the middle of a block. If this happens we need to increment the block number,
            // because exodus can only deal with one element type per block
            if((prev_elem_nodes != -1) && (prev_elem_nodes != n_elem_nodes))
            {
              block_nums.push_back(++block_num);
              block_elem_types.push_back(block_elem_type);
              std::cout << "ELEMENT TYPE CHANGE MID BLOCK" << std::endl;
              std::cout << "THE GUILTY BLOCK IS: Block " << block_num << ", " << std::endl;
            }

            // Now we can make a new element
            Elem* elem = mesh.add_elem(Elem::build_with_id(elem_type, iel++));
            elem->subdomain_id() = block_num;
            
            // loop over the nodes adding them to the element
            for ( int i = 0 ; i < n_elem_nodes ; i++)
            {
              // setup the element with the correct node ids;
              // i have no idea if the node ordering is the same for 
              // ansys elements and libmesh elements? 
              // im sure life will be cruel

              // Note: life was cruel
              elem->set_node(i) = mesh.node_ptr(ansys_to_libmesh_node_id_map[nodes[elem_map[i]]]);
            }
            // Add element to mesh
            mesh.add_elem(elem);
            // Read a new line
            std::getline(in, s);
            // Update previous element node count
            prev_elem_nodes = n_elem_nodes;
          }    
          // we should now have a collection of elements in the map
          // and these can now belong to a block 
          std::getline(in, s);
          // tokenize the string
          std::vector<std::string> tokens = tokenize(s);
          // the new block name is given by the 2nd token
          std::string block_name = tokens[1];
          // Since libmesh is using boost anyway, use boost trim to remove whitespace
          boost::algorithm::trim(block_name);
          // Block names with element type identifier
          if(block_elem_types.size() > 0)
          {
            int it = 0;
            for(auto& elem_type_str : block_elem_types)
            {
              mesh.subdomain_name(block_nums[it]) = block_name + "_" + elem_type_str;
              it++;
            }
          }
          // Increment block id
          block_num++;
        } 
        
        else if (s.find("CMBLOCK,") != std::string::npos) 
        {
          // Semi global ints to track how many nodes have been added to the "nodes" vector
          int nodeset_counter = 0;

          // Tokenize the string, as it contains the number of nodes in the nodeset
          std::vector<std::string> tokens = tokenize(s);

          // Make it long just in case there's a massive nodeset
          long int n_nodeset_nodes = std::stoi(tokens[3]);

          // Second token is the name of the nodeset
          std::string nodeset_name = tokens[1];

          // Given libmesh uses boost, use boost whitespace trim
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
            mesh.get_boundary_info().add_node(ansys_to_libmesh_node_id_map[node], nodeset_id);
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
