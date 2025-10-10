/** @file gsLibMeshb.h

    @brief Provides declaration of functions writing LibMeshb files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Gobe
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsExport.h>

#include <string>

namespace gismo {

/// \brief Export a gsMultiPatch geometry to a .meshb file as an unstructured mesh
///
/// This function samples each patch of a multipatch geometry at a regular grid
/// and writes the resulting unstructured mesh to a LibMeshb format file.
/// The mesh consists of vertices and cells (edges for 1D, quads for 2D, hexahedra for 3D).
///
/// \param mPatch The multipatch geometry to export
/// \param fn Filename where the .meshb file will be written (without extension)
/// \param npts Number of sampling points per dimension for each patch
/// \param version LibMeshb version (2 or 3, default: 2)
///
/// \ingroup IO
template<class T>
void gsWriteLibMeshb(const gsMultiPatch<T> & mPatch,
                     std::string const & fn,
                     unsigned npts = 1000,
                     int version = 2);

/// \brief Export a gsField (solution on geometry) to a .meshb file as an unstructured mesh with solution data
///
/// This function samples each patch of a field at a regular grid and writes
/// the resulting unstructured mesh with solution data to a LibMeshb format file.
///
/// \param field The field object containing geometry and solution
/// \param fn Filename where the .meshb file will be written (without extension)
/// \param npts Number of sampling points per dimension for each patch
/// \param version LibMeshb version (2 or 3, default: 2)
///
/// \ingroup IO
template<class T>
void gsWriteLibMeshb(const gsField<T> & field,
                     std::string const & fn,
                     unsigned npts = 1000,
                     int version = 2);

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsLibMeshb.hpp)
#endif
