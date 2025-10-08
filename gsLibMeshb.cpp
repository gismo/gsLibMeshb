/** @file gsLibMeshb.cpp

    @brief Provides implementation of functions writing LibMeshb files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Gobe
*/

#include <gsCore/gsTemplateTools.h>
#include <gsLibMeshb/gsLibMeshb.h>
#include <gsLibMeshb/gsLibMeshb.hpp>

namespace gismo
{
    // Explicit instantiations for common types
    TEMPLATE_INST void gsWriteLibMeshb(const gsMultiPatch<real_t>&, std::string const&, unsigned, int);
    TEMPLATE_INST void gsWriteLibMeshb(const gsField<real_t>&, std::string const&, unsigned, int);

} // namespace gismo
