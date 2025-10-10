/** @file gsLibMeshb.hpp

    @brief Provides implementation of functions writing LibMeshb files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Gobe
*/

#pragma once

#include <libmeshb7.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsField.h>
#include <gsCore/gsGeometry.h>
#include <gsIO/gsOptionList.h>
#include <vector>
#include <set>

namespace gismo
{

namespace internal
{
    /// Helper function to get uniform sample counts
    template<typename T>
    inline gsVector<unsigned> uniformSampleCount(const gsVector<T>& a,
                                                  const gsVector<T>& b,
                                                  unsigned npts)
    {
        const int d = a.size();
        gsVector<unsigned> np(d);

        if (d == 1)
        {
            np[0] = npts;
        }
        else
        {
            gsVector<T> range = b - a;
            T total = range.sum();
            for (int i = 0; i < d; ++i)
            {
                np[i] = std::max(2u, static_cast<unsigned>(std::ceil(npts * range[i] / total)));
            }
        }
        return np;
    }

    /// Helper function to create a uniformly spaced grid of points
    template<typename T>
    inline gsMatrix<T> gsPointGrid(const gsVector<T>& a,
                                   const gsVector<T>& b,
                                   const gsVector<unsigned>& np)
    {
        const int d = a.size();
        std::vector<gsVector<T>> coords(d);
        index_t totalPts = 1;

        for (int i = 0; i < d; ++i)
        {
            coords[i].resize(np[i]);
            for (unsigned j = 0; j < np[i]; ++j)
            {
                coords[i][j] = a[i] + (b[i] - a[i]) * static_cast<T>(j) / static_cast<T>(np[i] - 1);
            }
            totalPts *= np[i];
        }

        gsMatrix<T> pts(d, totalPts);
        index_t idx = 0;

        if (d == 1)
        {
            for (unsigned i = 0; i < np[0]; ++i)
            {
                pts(0, idx++) = coords[0][i];
            }
        }
        else if (d == 2)
        {
            for (unsigned j = 0; j < np[1]; ++j)
            {
                for (unsigned i = 0; i < np[0]; ++i)
                {
                    pts(0, idx) = coords[0][i];
                    pts(1, idx) = coords[1][j];
                    ++idx;
                }
            }
        }
        else if (d == 3)
        {
            for (unsigned k = 0; k < np[2]; ++k)
            {
                for (unsigned j = 0; j < np[1]; ++j)
                {
                    for (unsigned i = 0; i < np[0]; ++i)
                    {
                        pts(0, idx) = coords[0][i];
                        pts(1, idx) = coords[1][j];
                        pts(2, idx) = coords[2][k];
                        ++idx;
                    }
                }
            }
        }

        return pts;
    }
} // namespace internal

/// Export a gsMultiPatch to a .meshb file as an unstructured mesh
template<class T>
void gsWriteLibMeshb(const gsMultiPatch<T> & mPatch,
                     std::string const & fn,
                     unsigned npts,
                     int version)
{
    std::string filename = fn;
    if (filename.substr(filename.length() - 6) != ".meshb")
        filename += ".meshb";

    // Collect data from all patches
    std::vector<gsMatrix<T>> patchPoints;
    std::vector<gsVector<unsigned>> patchNp;
    std::vector<index_t> patchPointCounts;
    std::vector<index_t> patchCells;
    index_t totalPoints = 0;
    index_t totalCells = 0;
    int maxDim = 0;

    // First pass: evaluate all patches and count points/cells
    for (index_t i = 0; i < mPatch.nPatches(); ++i)
    {
        const gsGeometry<T>& geo = mPatch.patch(i);
        const int d = geo.domainDim();
        const int n = geo.targetDim();
        maxDim = std::max(maxDim, n);

        gsMatrix<T> ab = geo.support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        gsVector<unsigned> np = internal::uniformSampleCount(a, b, npts);
        patchNp.push_back(np);

        gsMatrix<T> pts = internal::gsPointGrid(a, b, np);
        gsMatrix<T> eval_geo = geo.eval(pts);

        // Pad to 3D if needed (meshb requires 3D coordinates)
        if (n < 3)
        {
            eval_geo.conservativeResize(3, eval_geo.cols());
            for (index_t row = n; row < 3; ++row)
                eval_geo.row(row).setZero();
        }

        patchPoints.push_back(eval_geo);
        totalPoints += eval_geo.cols();
        patchPointCounts.push_back(eval_geo.cols());

        // Calculate number of cells for this patch
        index_t cellsInPatch = 1;
        for (int dim = 0; dim < d; ++dim)
            cellsInPatch *= (np[dim] - 1);
        totalCells += cellsInPatch;
        patchCells.push_back(cellsInPatch);
    }

    // Determine dimension for the mesh file
    int dim = std::max(maxDim, 2);
    int domainDim = mPatch.patch(0).domainDim();

    // Open mesh file for writing
    int64_t meshIdx = GmfOpenMesh(filename.c_str(), GmfWrite, version, dim);
    if (!meshIdx)
    {
        gsWarn << "gsWriteLibMeshb: Cannot open file \"" << filename << "\" for writing\n";
        return;
    }

    // Write vertices with boundary tags
    GmfSetKwd(meshIdx, GmfVertices, totalPoints);
    index_t pointOffset = 0;
    for (size_t patchIdx = 0; patchIdx < patchPoints.size(); ++patchIdx)
    {
        const gsMatrix<T>& points = patchPoints[patchIdx];
        const gsVector<unsigned>& np = patchNp[patchIdx];

        for (index_t j = 0; j < points.cols(); ++j)
        {
            // Determine if vertex is on boundary
            int ref = static_cast<int>(patchIdx) + 1; // Patch ID as reference

            // Mark boundary vertices
            if (domainDim == 2)
            {
                index_t i = j % np[0];
                index_t jj = j / np[0];
                if (i == 0 || i == np[0]-1 || jj == 0 || jj == np[1]-1)
                    ref = -ref; // Negative for boundary
            }
            else if (domainDim == 3)
            {
                index_t i = j % np[0];
                index_t jj = (j / np[0]) % np[1];
                index_t k = j / (np[0] * np[1]);
                if (i == 0 || i == np[0]-1 || jj == 0 || jj == np[1]-1 || k == 0 || k == np[2]-1)
                    ref = -ref; // Negative for boundary
            }

            GmfSetLin(meshIdx, GmfVertices,
                     static_cast<double>(points(0, j)),
                     static_cast<double>(points(1, j)),
                     static_cast<double>(points(2, j)),
                     ref);
        }
    }

    // Write cells
    pointOffset = 0;

    for (size_t patchIdx = 0; patchIdx < patchPoints.size(); ++patchIdx)
    {
        const gsVector<unsigned>& np = patchNp[patchIdx];
        const int d = mPatch.patch(patchIdx).domainDim();
        const index_t cellsInPatch = patchCells[patchIdx];

        if (d == 1)
        {
            // 1D: Write edges
            if (patchIdx == 0)
                GmfSetKwd(meshIdx, GmfEdges, totalCells);

            for (index_t i = 0; i < np[0] - 1; ++i)
            {
                index_t i0 = pointOffset + i + 1;      // 1-based indexing
                index_t i1 = pointOffset + i + 2;
                GmfSetLin(meshIdx, GmfEdges, i0, i1, static_cast<int>(patchIdx) + 1);
            }
        }
        else if (d == 2)
        {
            // 2D: Write quadrilaterals
            if (patchIdx == 0)
                GmfSetKwd(meshIdx, GmfQuadrilaterals, totalCells);

            for (index_t j = 0; j < np[1] - 1; ++j)
            {
                for (index_t i = 0; i < np[0] - 1; ++i)
                {
                    index_t i0 = pointOffset + i + j * np[0] + 1;     // 1-based indexing
                    index_t i1 = i0 + 1;
                    index_t i2 = i0 + np[0] + 1;
                    index_t i3 = i0 + np[0];

                    GmfSetLin(meshIdx, GmfQuadrilaterals, i0, i1, i2, i3,
                             static_cast<int>(patchIdx) + 1);
                }
            }
        }
        else if (d == 3)
        {
            // 3D: Write hexahedra
            if (patchIdx == 0)
                GmfSetKwd(meshIdx, GmfHexahedra, totalCells);

            for (index_t k = 0; k < np[2] - 1; ++k)
            {
                for (index_t j = 0; j < np[1] - 1; ++j)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + j * np[0] + k * np[0] * np[1] + 1; // 1-based
                        index_t i1 = i0 + 1;
                        index_t i2 = i0 + np[0] + 1;
                        index_t i3 = i0 + np[0];
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i5 = i4 + 1;
                        index_t i6 = i4 + np[0] + 1;
                        index_t i7 = i4 + np[0];

                        GmfSetLin(meshIdx, GmfHexahedra, i0, i1, i2, i3, i4, i5, i6, i7,
                                 static_cast<int>(patchIdx) + 1);
                    }
                }
            }
        }

        pointOffset += patchPointCounts[patchIdx];
    }

    // Extract and write boundary quads for 3D meshes
    if (maxDim == 3)
    {
        std::vector<std::array<index_t, 5>> boundaryQuads; // 4 vertices + 1 ref
        pointOffset = 0;

        for (size_t patchIdx = 0; patchIdx < patchPoints.size(); ++patchIdx)
        {
            const gsVector<unsigned>& np = patchNp[patchIdx];
            const int d = mPatch.patch(patchIdx).domainDim();

            if (d == 3)
            {
                // Extract 6 boundary faces of the structured grid
                // Face at k=0 (bottom)
                for (index_t j = 0; j < np[1] - 1; ++j)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + j * np[0] + 1;
                        index_t i1 = i0 + 1;
                        index_t i2 = i0 + np[0] + 1;
                        index_t i3 = i0 + np[0];
                        boundaryQuads.push_back({i0, i1, i2, i3, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at k=np[2]-1 (top)
                index_t kOffset = (np[2] - 1) * np[0] * np[1];
                for (index_t j = 0; j < np[1] - 1; ++j)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + j * np[0] + kOffset + 1;
                        index_t i1 = i0 + 1;
                        index_t i2 = i0 + np[0] + 1;
                        index_t i3 = i0 + np[0];
                        boundaryQuads.push_back({i0, i1, i2, i3, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at j=0 (front)
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + k * np[0] * np[1] + 1;
                        index_t i1 = i0 + 1;
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i5 = i4 + 1;
                        boundaryQuads.push_back({i0, i1, i5, i4, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at j=np[1]-1 (back)
                index_t jOffset = (np[1] - 1) * np[0];
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + jOffset + k * np[0] * np[1] + 1;
                        index_t i1 = i0 + 1;
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i5 = i4 + 1;
                        boundaryQuads.push_back({i0, i1, i5, i4, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at i=0 (left)
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t j = 0; j < np[1] - 1; ++j)
                    {
                        index_t i0 = pointOffset + j * np[0] + k * np[0] * np[1] + 1;
                        index_t i3 = i0 + np[0];
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i7 = i4 + np[0];
                        boundaryQuads.push_back({i0, i3, i7, i4, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at i=np[0]-1 (right)
                index_t iOffset = np[0] - 1;
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t j = 0; j < np[1] - 1; ++j)
                    {
                        index_t i1 = pointOffset + iOffset + j * np[0] + k * np[0] * np[1] + 1;
                        index_t i2 = i1 + np[0];
                        index_t i5 = i1 + np[0] * np[1];
                        index_t i6 = i5 + np[0];
                        boundaryQuads.push_back({i1, i2, i6, i5, static_cast<index_t>(patchIdx) + 1});
                    }
                }
            }

            pointOffset += patchPointCounts[patchIdx];
        }

        // Write all boundary quads to a separate section
        if (!boundaryQuads.empty())
        {
            GmfSetKwd(meshIdx, GmfQuadrilaterals, boundaryQuads.size());
            for (const auto& quad : boundaryQuads)
            {
                GmfSetLin(meshIdx, GmfQuadrilaterals, quad[0], quad[1], quad[2], quad[3], (int)quad[4]);
            }
            gsInfo << "Written " << boundaryQuads.size() << " boundary quads\n";
        }
    }

    // Close the mesh file
    GmfCloseMesh(meshIdx);

    gsInfo << "Written " << totalPoints << " vertices and " << totalCells
           << " cells to " << filename << "\n";
}

/// Export a gsField to a .meshb file with solution data
template<class T>
void gsWriteLibMeshb(const gsField<T> & field,
                     std::string const & fn,
                     unsigned npts,
                     int version)
{
    std::string meshFilename = fn;
    if (meshFilename.substr(meshFilename.length() - 6) != ".meshb")
        meshFilename += ".meshb";

    std::string solFilename = fn;
    if (solFilename.substr(solFilename.length() - 5) != ".solb")
    {
        // Remove .meshb if present
        if (solFilename.substr(solFilename.length() - 6) == ".meshb")
            solFilename = solFilename.substr(0, solFilename.length() - 6);
        solFilename += ".solb";
    }

    // Collect data from all patches
    std::vector<gsMatrix<T>> patchPoints;
    std::vector<gsMatrix<T>> patchFields;
    std::vector<gsVector<unsigned>> patchNp;
    std::vector<index_t> patchPointCounts;
    std::vector<index_t> patchCells;
    index_t totalPoints = 0;
    index_t totalCells = 0;
    int maxDim = 0;
    int fieldDim = 0;

    // First pass: evaluate all patches and count points/cells
    for (index_t i = 0; i < field.nPieces(); ++i)
    {
        const gsFunction<T>& geo = field.patch(i);
        const gsFunction<T>& func = field.function(i);
        const int d = geo.domainDim();
        const int n = geo.targetDim();
        maxDim = std::max(maxDim, n);
        fieldDim = func.targetDim();

        gsMatrix<T> ab = geo.support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        gsVector<unsigned> np = internal::uniformSampleCount(a, b, npts);
        patchNp.push_back(np);

        gsMatrix<T> pts = internal::gsPointGrid(a, b, np);
        gsMatrix<T> eval_geo = geo.eval(pts);
        gsMatrix<T> eval_field;
        if (field.isParametric()) {
            eval_field = func.eval(pts);
        } else {
            eval_field = func.eval(eval_geo);
        }
        
        // Pad geometry to 3D if needed
        if (n < 3)
        {
            eval_geo.conservativeResize(3, eval_geo.cols());
            for (index_t row = n; row < 3; ++row)
                eval_geo.row(row).setZero();
        }

        patchPoints.push_back(eval_geo);
        patchFields.push_back(eval_field);
        totalPoints += eval_geo.cols();
        patchPointCounts.push_back(eval_geo.cols());

        // Calculate number of cells
        index_t cellsInPatch = 1;
        for (int dim = 0; dim < d; ++dim)
            cellsInPatch *= (np[dim] - 1);
        totalCells += cellsInPatch;
        patchCells.push_back(cellsInPatch);
    }

    int dim = std::max(maxDim, 2);
    int domainDim = field.patch(0).domainDim();

    // Write mesh file
    int64_t meshIdx = GmfOpenMesh(meshFilename.c_str(), GmfWrite, version, dim);
    if (!meshIdx)
    {
        gsWarn << "gsWriteLibMeshb: Cannot open file \"" << meshFilename << "\" for writing\n";
        return;
    }

    // Write vertices
    GmfSetKwd(meshIdx, GmfVertices, totalPoints);
    for (size_t patchIdx = 0; patchIdx < patchPoints.size(); ++patchIdx)
    {
        const gsMatrix<T>& points = patchPoints[patchIdx];
        for (index_t j = 0; j < points.cols(); ++j)
        {
            GmfSetLin(meshIdx, GmfVertices,
                     static_cast<double>(points(0, j)),
                     static_cast<double>(points(1, j)),
                     static_cast<double>(points(2, j)),
                     static_cast<int>(patchIdx) + 1);
        }
    }

    // Write cells
    index_t pointOffset = 0;
    for (size_t patchIdx = 0; patchIdx < patchPoints.size(); ++patchIdx)
    {
        const gsVector<unsigned>& np = patchNp[patchIdx];
        const int d = field.patch(patchIdx).domainDim();
        const index_t cellsInPatch = patchCells[patchIdx];

        if (d == 1)
        {
            if (patchIdx == 0)
                GmfSetKwd(meshIdx, GmfEdges, totalCells);

            for (index_t i = 0; i < np[0] - 1; ++i)
            {
                index_t i0 = pointOffset + i + 1;
                index_t i1 = pointOffset + i + 2;
                GmfSetLin(meshIdx, GmfEdges, i0, i1, static_cast<int>(patchIdx) + 1);
            }
        }
        else if (d == 2)
        {
            if (patchIdx == 0)
                GmfSetKwd(meshIdx, GmfQuadrilaterals, totalCells);

            for (index_t j = 0; j < np[1] - 1; ++j)
            {
                for (index_t i = 0; i < np[0] - 1; ++i)
                {
                    index_t i0 = pointOffset + i + j * np[0] + 1;
                    index_t i1 = i0 + 1;
                    index_t i2 = i0 + np[0] + 1;
                    index_t i3 = i0 + np[0];
                    GmfSetLin(meshIdx, GmfQuadrilaterals, i0, i1, i2, i3,
                             static_cast<int>(patchIdx) + 1);
                }
            }
        }
        else if (d == 3)
        {
            if (patchIdx == 0)
                GmfSetKwd(meshIdx, GmfHexahedra, totalCells);

            for (index_t k = 0; k < np[2] - 1; ++k)
            {
                for (index_t j = 0; j < np[1] - 1; ++j)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + j * np[0] + k * np[0] * np[1] + 1;
                        index_t i1 = i0 + 1;
                        index_t i2 = i0 + np[0] + 1;
                        index_t i3 = i0 + np[0];
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i5 = i4 + 1;
                        index_t i6 = i4 + np[0] + 1;
                        index_t i7 = i4 + np[0];
                        GmfSetLin(meshIdx, GmfHexahedra, i0, i1, i2, i3, i4, i5, i6, i7,
                                 static_cast<int>(patchIdx) + 1);
                    }
                }
            }
        }

        pointOffset += patchPointCounts[patchIdx];
    }

    // Extract and write boundary quads for 3D meshes
    if (maxDim == 3)
    {
        std::vector<std::array<index_t, 5>> boundaryQuads; // 4 vertices + 1 ref
        pointOffset = 0;

        for (size_t patchIdx = 0; patchIdx < patchPoints.size(); ++patchIdx)
        {
            const gsVector<unsigned>& np = patchNp[patchIdx];
            const int d = field.patch(patchIdx).domainDim();

            if (d == 3)
            {
                // Extract 6 boundary faces of the structured grid
                // Face at k=0 (bottom)
                for (index_t j = 0; j < np[1] - 1; ++j)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + j * np[0] + 1;
                        index_t i1 = i0 + 1;
                        index_t i2 = i0 + np[0] + 1;
                        index_t i3 = i0 + np[0];
                        boundaryQuads.push_back({i0, i1, i2, i3, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at k=np[2]-1 (top)
                index_t kOffset = (np[2] - 1) * np[0] * np[1];
                for (index_t j = 0; j < np[1] - 1; ++j)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + j * np[0] + kOffset + 1;
                        index_t i1 = i0 + 1;
                        index_t i2 = i0 + np[0] + 1;
                        index_t i3 = i0 + np[0];
                        boundaryQuads.push_back({i0, i1, i2, i3, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at j=0 (front)
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + k * np[0] * np[1] + 1;
                        index_t i1 = i0 + 1;
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i5 = i4 + 1;
                        boundaryQuads.push_back({i0, i1, i5, i4, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at j=np[1]-1 (back)
                index_t jOffset = (np[1] - 1) * np[0];
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t i = 0; i < np[0] - 1; ++i)
                    {
                        index_t i0 = pointOffset + i + jOffset + k * np[0] * np[1] + 1;
                        index_t i1 = i0 + 1;
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i5 = i4 + 1;
                        boundaryQuads.push_back({i0, i1, i5, i4, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at i=0 (left)
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t j = 0; j < np[1] - 1; ++j)
                    {
                        index_t i0 = pointOffset + j * np[0] + k * np[0] * np[1] + 1;
                        index_t i3 = i0 + np[0];
                        index_t i4 = i0 + np[0] * np[1];
                        index_t i7 = i4 + np[0];
                        boundaryQuads.push_back({i0, i3, i7, i4, static_cast<index_t>(patchIdx) + 1});
                    }
                }

                // Face at i=np[0]-1 (right)
                index_t iOffset = np[0] - 1;
                for (index_t k = 0; k < np[2] - 1; ++k)
                {
                    for (index_t j = 0; j < np[1] - 1; ++j)
                    {
                        index_t i1 = pointOffset + iOffset + j * np[0] + k * np[0] * np[1] + 1;
                        index_t i2 = i1 + np[0];
                        index_t i5 = i1 + np[0] * np[1];
                        index_t i6 = i5 + np[0];
                        boundaryQuads.push_back({i1, i2, i6, i5, static_cast<index_t>(patchIdx) + 1});
                    }
                }
            }

            pointOffset += patchPointCounts[patchIdx];
        }

        // Write all boundary quads to a separate section
        if (!boundaryQuads.empty())
        {
            GmfSetKwd(meshIdx, GmfQuadrilaterals, boundaryQuads.size());
            for (const auto& quad : boundaryQuads)
            {
                GmfSetLin(meshIdx, GmfQuadrilaterals, quad[0], quad[1], quad[2], quad[3], (int)quad[4]);
            }
            gsInfo << "Written " << boundaryQuads.size() << " boundary quads\n";
        }
    }

    GmfCloseMesh(meshIdx);

    // Write solution file
    int64_t solIdx = GmfOpenMesh(solFilename.c_str(), GmfWrite, version, dim);
    if (!solIdx)
    {
        gsWarn << "gsWriteLibMeshb: Cannot open file \"" << solFilename << "\" for writing\n";
        return;
    }

    // Determine solution type based on field dimension
    int solType;

    if (fieldDim == 1)
    {
        solType = GmfSca;
    }
    else if (fieldDim >= 2 && fieldDim <= 3)
    {
        solType = GmfVec;
    }
    else
    {
        // For higher dimensions, write as scalar (first component)
        solType = GmfSca;
        gsWarn << "Field dimension " << fieldDim << " > 3, writing only first component\n";
    }

    // Write solution at vertices - LibMeshb expects array of types
    int types[1] = {solType};
    GmfSetKwd(solIdx, GmfSolAtVertices, (int64_t)totalPoints, (int64_t)1, types);

    index_t vertexCount = 0;
    for (size_t patchIdx = 0; patchIdx < patchFields.size(); ++patchIdx)
    {
        const gsMatrix<T>& fieldData = patchFields[patchIdx];
        for (index_t j = 0; j < fieldData.cols(); ++j)
        {
            if (fieldDim == 1)
            {
                // For scalar: pass pointer to double value
                double val = static_cast<double>(fieldData(0, j));
                GmfSetLin(solIdx, GmfSolAtVertices, &val);
            }
            else if (fieldDim == 2)
            {
                // For vector: pass array pointer (pad 2D to 3D)
                double vecData[3] = {
                    static_cast<double>(fieldData(0, j)),
                    static_cast<double>(fieldData(1, j)),
                    0.0
                };
                GmfSetLin(solIdx, GmfSolAtVertices, vecData);
            }
            else if (fieldDim == 3)
            {
                // For vector: pass array pointer
                double vecData[3] = {
                    static_cast<double>(fieldData(0, j)),
                    static_cast<double>(fieldData(1, j)),
                    static_cast<double>(fieldData(2, j))
                };
                GmfSetLin(solIdx, GmfSolAtVertices, vecData);
            }
            else
            {
                // Write first component only
                double val = static_cast<double>(fieldData(0, j));
                GmfSetLin(solIdx, GmfSolAtVertices, &val);
            }
            vertexCount++;
        }
    }

    GmfCloseMesh(solIdx);

    gsInfo << "Written " << totalPoints << " vertices and " << totalCells
           << " cells to " << meshFilename << "\n";
    gsInfo << "Written " << vertexCount << " solution values (" << fieldDim << "D) to " << solFilename << "\n";
}

} // namespace gismo
