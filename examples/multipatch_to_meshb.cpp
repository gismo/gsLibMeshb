/** @file helloSubmodule.cpp

    @brief First example of submodule

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    std::string fn("");
    gsCmdLine cmd("Export meshb file.");
    cmd.addPlainString("filename", "File containing geometry data.", fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    if (!fn.empty())
    {
        gsReadFile<>(fn, mp);
        
    }
  
  return EXIT_SUCCESS;
}
