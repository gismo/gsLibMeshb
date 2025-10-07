/** @file multipatch_to_surface.cpp

    @brief
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.


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
        gsReadFile<>(fn, mp);
    else
        return EXIT_SUCCESS;

    

    //-------------------- Standard libmeshb example


       int i, NmbVer, NmbQad, ver, dim, *RefTab, (*QadTab)[5];
   long long InpMsh, OutMsh;
   double (*VerTab)[3];


   /*-----------------------------------*/
   /* Open mesh file "quad.meshb"       */
   /*-----------------------------------*/

   if(!(InpMsh = GmfOpenMesh("../sample_meshes/quad.mesh", GmfRead, &ver, &dim)))
      return(1);

   printf("InpMsh: idx = %lld, version = %d, dimension = %d\n", InpMsh, ver, dim);

   if( (ver != 2) || (dim != 3) )
      exit(1);

   // Read the number of vertices and allocate memory
   NmbVer = (int)GmfStatKwd(InpMsh, GmfVertices);
   printf("InpMsh: nmb vertices = %d\n", NmbVer);
   VerTab = malloc((NmbVer+1) * 3 * sizeof(double));
   RefTab = malloc((NmbVer+1) * sizeof(int));

   // Read the number of quads and allocate memory
   NmbQad = (int)GmfStatKwd(InpMsh, GmfQuadrilaterals);
   printf("InpMsh: nmb quads = %d\n", NmbQad);
   QadTab = malloc((NmbQad+1) * 5 * sizeof(int));

   // Read the vertices
   GmfGotoKwd(InpMsh, GmfVertices);

   for(i=1;i<=NmbVer;i++)
      GmfGetLin(  InpMsh, GmfVertices, &VerTab[i][0], &VerTab[i][1],
                  &VerTab[i][2], &RefTab[i] );

   // Read the quads
   GmfGotoKwd(InpMsh, GmfQuadrilaterals);

   for(i=1;i<=NmbQad;i++)
      GmfGetLin(  InpMsh, GmfQuadrilaterals, &QadTab[i][0], &QadTab[i][1],
                  &QadTab[i][2], &QadTab[i][3], &QadTab[i][4] );

   // Close the quad mesh
   GmfCloseMesh(InpMsh);


   /*-----------------------------------*/
   /* Write the triangle mesh           */
   /*-----------------------------------*/

   if(!(OutMsh = GmfOpenMesh("tri.mesh", GmfWrite, ver, dim)))
      return(1);

   // Write the vertices
   GmfSetKwd(OutMsh, GmfVertices, NmbVer);

   for(i=1;i<=NmbVer;i++)
      GmfSetLin(  OutMsh, GmfVertices, VerTab[i][0],
                  VerTab[i][1], VerTab[i][2], RefTab[i] );

   // Write the triangles
   GmfSetKwd(OutMsh, GmfTriangles, 2*NmbQad);

   // Split each quad into two triangles on the fly
   for(i=1;i<=NmbQad;i++)
   {
      GmfSetLin(  OutMsh, GmfTriangles, QadTab[i][0],
                  QadTab[i][1], QadTab[i][2], QadTab[i][4] );
      GmfSetLin(  OutMsh, GmfTriangles, QadTab[i][0],
                  QadTab[i][2], QadTab[i][3], QadTab[i][4] );
   }

   // Do not forget to close the mesh file
   GmfCloseMesh(OutMsh);
   printf("OutMsh: nmb triangles = %d\n", 2 * NmbQad);

   free(QadTab);
   free(RefTab);
   free(VerTab);
   
    
  return EXIT_SUCCESS;
}
