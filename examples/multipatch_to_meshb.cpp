/** @file multipatch_to_meshb.cpp

    @brief Export a gsMultiPatch geometry to a LibMeshb .meshb file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Gobe
*/

#include <gismo.h>
#include <gsLibMeshb/gsLibMeshb.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    std::string fn("");
    index_t npts = 10;

    gsCmdLine cmd("Export a multipatch geometry to LibMeshb .meshb format.");
    cmd.addPlainString("filename", "File containing geometry data (XML, G2, etc.).", fn);
    cmd.addInt("n", "npts", "Number of sampling points per dimension", npts);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read geometry]
    if (fn.empty())
    {
        gsInfo << "No input file specified. Please provide a geometry file.\n";
        gsInfo << "Usage: " << argv[0] << " <geometry_file> [-n npts] [-v version]\n";
        return EXIT_FAILURE;
    }

    gsMultiPatch<> mp;
    gsReadFile<>(fn, mp);

    gsInfo << "Loaded multipatch with " << mp.nPatches() << " patches\n";
    gsInfo << "Domain dimension: " << mp.domainDim() << "\n";
    gsInfo << "Target dimension: " << mp.targetDim() << "\n";
    //! [Read geometry]

    //! [Export to meshb]
    std::string outFilename = fn;
    // Remove extension from input filename
    size_t lastDot = outFilename.find_last_of(".");
    if (lastDot != std::string::npos)
        outFilename = outFilename.substr(0, lastDot);

    gsInfo << "\nExporting to " << outFilename << ".meshb with " << npts
           << " points per dimension...\n";

    // Use LibMeshb version 2 (default)
    gsWriteLibMeshb(mp, outFilename, static_cast<unsigned>(npts), 2);

    gsInfo << "Export complete!\n";
    //! [Export to meshb]

    return EXIT_SUCCESS;
}
