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
   * Reads in a mesh in the Ansys *.cdb format from the ASCII file
   * given by name.
   *
   * \note The user is responsible for calling Mesh::prepare_for_use()
   * after reading the mesh and before using it.
   *
   * The physical group names defined in the Gmsh-file are stored, depending
   * on the element dimension, as either subdomain name or as side name. The
   * IDs of the former can be retrieved by using the MeshBase::get_id_by_name
   * method; the IDs of the latter can be retrieved by using the
   * MeshBase::get_boundary_info and the BoundaryInfo::get_id_by_name methods.
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
   * A static ElementMaps object that is built statically and used by
   * all instances of this class.
   */
  static ElementMaps _element_maps;

  /**
   * A static function used to construct the _element_maps struct,
   * statically.
   */
  static ElementMaps build_element_maps();
};

} // namespace libMesh

#endif // LIBMESH_CDB_IO_H