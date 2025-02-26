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



#ifndef LIBMESH_ELEM_QUALITY_H
#define LIBMESH_ELEM_QUALITY_H

// C++ includes
#include <vector>
#include <string>

namespace libMesh
{

// Forward declarations
enum ElemType : int;
enum ElemQuality : int;

/**
 * A namespace for quality utility functions.
 *
 * \author John W. Peterson
 * \date 2002
 * \brief Utility functions for computing element quality indicators.
 */
namespace Quality
{
/**
 * The number of element quality types we have
 * defined.  This needs to be updated if you
 * add one.
 */
const unsigned int num_quals = 16;

/**
 * \returns A descriptive name for a \p ElemQuality
 * \p enum
 */
std::string name (const ElemQuality q);

/**
 * \returns A description for a \p ElemQuality
 * \p enum
 */
std::string describe (const ElemQuality q);

/**
 * \returns The valid \p ElemQuality metrics for a given
 * \p ElemType element type.
 */
std::vector<ElemQuality> valid (const ElemType t);
}

} // namespace libMesh

#endif // LIBMESH_ELEM_QUALITY_H
