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



#ifndef LIBMESH_CDB_IO_H
#define LIBMESH_CDB_IO_H

#include "boost/algorithm/string/trim.hpp"

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"

// C++ includes
#include <cstddef>

namespace libMesh
{

// Forward declarations
class MeshBase;

/**
 * \brief Reading and writing meshes in the Ansys CDB format.
 *
 *
 * \author Andrew Davis
 * \date 2023
 */
class CDBIO : public MeshInput<MeshBase>,
               public MeshOutput<MeshBase>
{
public:

  CDBIO(){};
  /**
   * Constructor.  Takes a non-const Mesh reference which it
   * will fill up with elements via the read() command.
  */
  explicit
  CDBIO (MeshBase & mesh);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
  */
  explicit
  CDBIO (const MeshBase & mesh);

  /**
   * This method should implement writing a mesh to a specified file
   * in the *.cdb format. But it is currently not implemented.
  */
  virtual void write (const std::string & name) {};

  /**
   * Reads in a mesh in the Ansys *.cdb format from the ASCII file
   * given by name.
  */
  virtual void read (const std::string & name) override;

private:
  /**
   * Implementation of the read() function.  This function
   * is called by the public interface function and implements
   * reading the file.
  */
  void read_mesh (std::istream & in);

  /**
   * Defines mapping from libMesh element types to Gmsh element types or vice-versa.
  */
  struct AnsysElementDefinition {
    AnsysElementDefinition(unsigned int ansys_type_in,
                           unsigned int dim_in) :
      ansys_type(ansys_type_in),
      dim(dim_in)
    {}
    
    // Needs default constructor otherwise the compiler throws a hissy fit for unknown reasons??
    AnsysElementDefinition(){}

    // Struct data
    unsigned int ansys_type;
    unsigned int dim;

    /**
     *  Each ansys element type can actually refer to multiple element types,
     *  for example SOLID226 can be a HEX20, TET10, PYRAMID13 or PRISM15. 
     *  Therefore we need a map to store all of potential element mappings 
     *  and libmesh element types.
     * */
    std::map<unsigned int, std::vector<unsigned int>> ansys_node_ordering_map;

    std::map<unsigned int, libMesh::ElemType> ansys_to_libmesh_elem_type_map;

    std::map<unsigned int, std::string> ansys_to_libmesh_elem_type_string_map;

    /**
     *  Each ansys element type can actually refer to multiple element types,
     *  for example SOLID226 can be a HEX20, TET10, PYRAMID13 or PRISM15. 
     *  Therefore we need a map to store all of potential element "sub" mappings.
     */
    void add_elem_sub_mapping(const std::vector<unsigned int>& ansys_node_ordering,
                              const libMesh::ElemType& elem_type, const std::string& elem_type_string)
    {
      ansys_node_ordering_map.emplace(ansys_node_ordering.size(), ansys_node_ordering);
      ansys_to_libmesh_elem_type_map.emplace(ansys_node_ordering.size(), elem_type); 
      ansys_to_libmesh_elem_type_string_map.emplace(ansys_node_ordering.size(), elem_type_string); 
    }
  };

  /**
   * struct which holds a map from CDB to libMesh element numberings
   * and vice-versa.
   */
  struct CDBMaps : public std::map<unsigned int, AnsysElementDefinition>
  {
    // Helper function to add a (key, value) pair to both maps
    void add_def(const AnsysElementDefinition & eledef)
    {
      this->emplace(eledef.ansys_type, eledef);
    }
  };

  /**
   * A static ElementMaps object that is built statically and used by
   * all instances of this class.
   */
  static CDBMaps _cdb_maps;

  /**
   * A static function used to construct the _element_maps struct,
   * statically.
   */
  static CDBMaps build_element_maps();

  /**
   * Map to store mapping of global libmesh node ID's to global CDB node ID's 
  */
  std::map<unsigned int, unsigned int> ansys_to_libmesh_node_id_map;
};

} // namespace libMesh

#endif // LIBMESH_CDB_IO_H