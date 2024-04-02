// GetDP - Copyright (C) 1997-2022 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.
//
// Contributor(s):
//   Christophe Trophime
//

#include <math.h>
#include "GetDPConfig.h"
#include "Pos_Element.h"
#include "GeoData.h"
#include "Get_Geometry.h"
#include "Get_DofOfElement.h"
#include "Cal_Value.h"
#include "MallocUtils.h"
#include "Message.h"
#include "ProDefine.h"

#if defined(HAVE_GMSH)
#include <gmsh/BasisFactory.h>
#include <gmsh/nodalBasis.h>
#endif

extern struct CurrentData Current;

/* ------------------------------------------------------------------------ */
/*  Create/Destroy/Compare                                                  */
/* ------------------------------------------------------------------------ */

void Alloc_PostElement(struct PostElement *PostElement)
{
  PostElement->NumNodes = (int *)Malloc(PostElement->NbrNodes * sizeof(int));
  /* allocate as much as possible in one step */
  PostElement->u = (double *)Malloc(6 * PostElement->NbrNodes * sizeof(double));
  PostElement->v = &PostElement->u[PostElement->NbrNodes];
  PostElement->w = &PostElement->u[2 * PostElement->NbrNodes];
  PostElement->x = &PostElement->u[3 * PostElement->NbrNodes];
  PostElement->y = &PostElement->u[4 * PostElement->NbrNodes];
  PostElement->z = &PostElement->u[5 * PostElement->NbrNodes];
  PostElement->Value =
    (struct Value *)Malloc(PostElement->NbrNodes * sizeof(struct Value));
}

struct PostElement *Create_PostElement(int Index, int Type, int NbrNodes,
                                       int Depth)
{
  struct PostElement *PostElement;

  PostElement = (struct PostElement *)Malloc(sizeof(struct PostElement));
  PostElement->Index = Index;
  PostElement->Type = Type;
  PostElement->Depth = Depth;
  PostElement->NbrNodes = NbrNodes;
  if(NbrNodes > 0) Alloc_PostElement(PostElement);

  return PostElement;
}

void Destroy_PostElement(struct PostElement *PostElement)
{
  if(PostElement->NbrNodes > 0) {
    Free(PostElement->NumNodes);
    if(PostElement->u)
      Free(PostElement->u); /* normal case */
    else if(PostElement->x)
      Free(PostElement->x); /* partial copy */
    Free(PostElement->Value);
  }
  Free(PostElement);
}

struct PostElement *NodeCopy_PostElement(struct PostElement *PostElement)
{
  struct PostElement *Copy;
  int i;

  Copy = (struct PostElement *)Malloc(sizeof(struct PostElement));
  Copy->Index = PostElement->Index;
  Copy->Type = PostElement->Type;
  Copy->Depth = PostElement->Depth;
  Copy->NbrNodes = PostElement->NbrNodes;
  if(Copy->NbrNodes > 0) {
    Alloc_PostElement(Copy);
    for(i = 0; i < Copy->NbrNodes; i++) {
      Copy->NumNodes[i] = PostElement->NumNodes[i];
      Copy->u[i] = PostElement->u[i];
      Copy->v[i] = PostElement->v[i];
      Copy->w[i] = PostElement->w[i];
      Copy->x[i] = PostElement->x[i];
      Copy->y[i] = PostElement->y[i];
      Copy->z[i] = PostElement->z[i];
    }
  }

  return Copy;
}

struct PostElement *PartialCopy_PostElement(struct PostElement *PostElement)
{
  struct PostElement *Copy;
  int i;

  Copy = (struct PostElement *)Malloc(sizeof(struct PostElement));
  Copy->Index = PostElement->Index;
  Copy->Type = PostElement->Type;
  Copy->Depth = PostElement->Depth;
  Copy->NbrNodes = PostElement->NbrNodes;
  if(Copy->NbrNodes > 0) {
    Copy->NumNodes = NULL;
    Copy->u = Copy->v = Copy->w = NULL;
    /* allocate as much as possible in one step */
    Copy->x = (double *)Malloc(3 * Copy->NbrNodes * sizeof(double));
    Copy->y = &Copy->x[Copy->NbrNodes];
    Copy->z = &Copy->x[2 * Copy->NbrNodes];
    Copy->Value = (struct Value *)Malloc(Copy->NbrNodes * sizeof(struct Value));
    for(i = 0; i < Copy->NbrNodes; i++) {
      Copy->x[i] = PostElement->x[i];
      Copy->y[i] = PostElement->y[i];
      Copy->z[i] = PostElement->z[i];
      Cal_CopyValue(&PostElement->Value[i], &Copy->Value[i]);
    }
  }

  return Copy;
}

/* 2 PostElements never have the same barycenter unless they are identical */

int fcmp_PostElement(const void *a, const void *b)
{
  struct PostElement *PE1, *PE2;
  double s1, s2, TOL = Current.GeoData->CharacteristicLength * 1.e-12;
  int i;

  PE1 = *(struct PostElement **)a;
  PE2 = *(struct PostElement **)b;

  if(PE1->NbrNodes != PE2->NbrNodes) return PE1->NbrNodes - PE2->NbrNodes;

  s1 = s2 = 0.0;
  for(i = 0; i < PE1->NbrNodes; i++) {
    s1 += PE1->x[i];
    s2 += PE2->x[i];
  }
  if(s1 - s2 > TOL)
    return 1;
  else if(s1 - s2 < -TOL)
    return -1;

  s1 = s2 = 0.0;
  for(i = 0; i < PE1->NbrNodes; i++) {
    s1 += PE1->y[i];
    s2 += PE2->y[i];
  }
  if(s1 - s2 > TOL)
    return 1;
  else if(s1 - s2 < -TOL)
    return -1;

  s1 = s2 = 0.0;
  for(i = 0; i < PE1->NbrNodes; i++) {
    s1 += PE1->z[i];
    s2 += PE2->z[i];
  }
  if(s1 - s2 > TOL)
    return 1;
  else if(s1 - s2 < -TOL)
    return -1;

  return 0;
}

int fcmp_PostElement_v0(const void *a, const void *b)
{
  return (int)((*(struct PostElement **)a)->v[0] -
               (*(struct PostElement **)b)->v[0]);
}

int fcmp_PostElement_absu0(const void *a, const void *b)
{
  return (int)(fabs((*(struct PostElement **)b)->u[0]) -
               fabs((*(struct PostElement **)a)->u[0]));
}

/* ------------------------------------------------------------------------ */
/*  C u t _ P o s t E l e m e n t                                           */
/* ------------------------------------------------------------------------ */

void Cut_PostElement(struct PostElement *PE, struct Geo_Element *GE,
                     List_T *PE_L, int Index, int Depth, int Skin,
                     int DecomposeInSimplex)
{
  struct Element E;
  struct PostElement *C[8];

  double u01, u02, u03, u12, u13, u23;
  double v01, v02, v03, v12, v13, v23;
  double w01, w02, w03, w12, w13, w23;
  int i, j, NbCut = 0;

  /* Recursive division */

  if(PE->Depth < Depth) {
    switch(PE->Type) {
    case POINT_ELEMENT:
      Message::Error("Impossible to divide a Point recursively");
      break;

    case LINE:
    case LINE_2:
    case LINE_3:
    case LINE_4:
      u01 = .5 * (PE->u[0] + PE->u[1]);
      v01 = .5 * (PE->v[0] + PE->v[1]);
      w01 = .5 * (PE->w[0] + PE->w[1]);

      C[0] = Create_PostElement(Index, LINE, 2, PE->Depth);
      C[0]->u[0] = PE->u[0];
      C[0]->v[0] = PE->v[0];
      C[0]->w[0] = PE->w[0];
      C[0]->u[1] = u01;
      C[0]->v[1] = v01;
      C[0]->w[1] = w01;

      C[1] = PE;
      C[1]->u[0] = u01;
      C[1]->v[0] = v01;
      C[1]->w[0] = w01;

      NbCut = 2;
      break;

    case TRIANGLE:
    case TRIANGLE_2:
    case TRIANGLE_3:
    case TRIANGLE_4:
      u01 = .5 * (PE->u[0] + PE->u[1]);
      u02 = .5 * (PE->u[0] + PE->u[2]);
      v01 = .5 * (PE->v[0] + PE->v[1]);
      v02 = .5 * (PE->v[0] + PE->v[2]);
      w01 = .5 * (PE->w[0] + PE->w[1]);
      w02 = .5 * (PE->w[0] + PE->w[2]);

      u12 = .5 * (PE->u[1] + PE->u[2]);
      v12 = .5 * (PE->v[1] + PE->v[2]);
      w12 = .5 * (PE->w[1] + PE->w[2]);

      C[0] = Create_PostElement(Index, TRIANGLE, 3, PE->Depth);
      C[0]->u[0] = PE->u[0];
      C[0]->v[0] = PE->v[0];
      C[0]->w[0] = PE->w[0];
      C[0]->u[1] = u01;
      C[0]->v[1] = v01;
      C[0]->w[1] = w01;
      C[0]->u[2] = u02;
      C[0]->v[2] = v02;
      C[0]->w[2] = w02;

      C[1] = Create_PostElement(Index, TRIANGLE, 3, PE->Depth);
      C[1]->u[0] = u01;
      C[1]->v[0] = v01;
      C[1]->w[0] = w01;
      C[1]->u[1] = PE->u[1];
      C[1]->v[1] = PE->v[1];
      C[1]->w[1] = PE->w[1];
      C[1]->u[2] = u12;
      C[1]->v[2] = v12;
      C[1]->w[2] = w12;

      C[2] = Create_PostElement(Index, TRIANGLE, 3, PE->Depth);
      C[2]->u[0] = u02;
      C[2]->v[0] = v02;
      C[2]->w[0] = w02;
      C[2]->u[1] = u12;
      C[2]->v[1] = v12;
      C[2]->w[1] = w12;
      C[2]->u[2] = PE->u[2];
      C[2]->v[2] = PE->v[2];
      C[2]->w[2] = PE->w[2];

      C[3] = PE;
      C[3]->u[0] = u01;
      C[3]->v[0] = v01;
      C[3]->w[0] = w01;
      C[3]->u[1] = u12;
      C[3]->v[1] = v12;
      C[3]->w[1] = w12;
      C[3]->u[2] = u02;
      C[3]->v[2] = v02;
      C[3]->w[2] = w02;

      NbCut = 4;
      break;

    case QUADRANGLE:
    case QUADRANGLE_2:
    case QUADRANGLE_3:
    case QUADRANGLE_4:
      u01 = .5 * (PE->u[0] + PE->u[1]);
      u12 = .5 * (PE->u[1] + PE->u[2]);
      v01 = .5 * (PE->v[0] + PE->v[1]);
      v12 = .5 * (PE->v[1] + PE->v[2]);
      w01 = .5 * (PE->w[0] + PE->w[1]);
      w12 = .5 * (PE->w[1] + PE->w[2]);

      u23 = .5 * (PE->u[2] + PE->u[3]);
      u03 = .5 * (PE->u[0] + PE->u[3]);
      v23 = .5 * (PE->v[2] + PE->v[3]);
      v03 = .5 * (PE->v[0] + PE->v[3]);
      w23 = .5 * (PE->w[2] + PE->w[3]);
      w03 = .5 * (PE->w[0] + PE->w[3]);

      u02 = .5 * (PE->u[0] + PE->u[2]);
      v02 = .5 * (PE->v[0] + PE->v[2]);
      w02 = .5 * (PE->w[0] + PE->w[2]);

      C[0] = Create_PostElement(Index, QUADRANGLE, 4, PE->Depth);
      C[0]->u[0] = PE->u[0];
      C[0]->v[0] = PE->v[0];
      C[0]->w[0] = PE->w[0];
      C[0]->u[1] = u01;
      C[0]->v[1] = v01;
      C[0]->w[1] = w01;
      C[0]->u[2] = u02;
      C[0]->v[2] = v02;
      C[0]->w[2] = w02;
      C[0]->u[3] = u03;
      C[0]->v[3] = v03;
      C[0]->w[3] = w03;

      C[1] = Create_PostElement(Index, QUADRANGLE, 4, PE->Depth);
      C[1]->u[0] = u01;
      C[1]->v[0] = v01;
      C[1]->w[0] = w01;
      C[1]->u[1] = PE->u[1];
      C[1]->v[1] = PE->v[1];
      C[1]->w[1] = PE->w[1];
      C[1]->u[2] = u12;
      C[1]->v[2] = v12;
      C[1]->w[2] = w12;
      C[1]->u[3] = u02;
      C[1]->v[3] = v02;
      C[1]->w[3] = w02;

      C[2] = Create_PostElement(Index, QUADRANGLE, 4, PE->Depth);
      C[2]->u[0] = u02;
      C[2]->v[0] = v02;
      C[2]->w[0] = w02;
      C[2]->u[1] = u12;
      C[2]->v[1] = v12;
      C[2]->w[1] = w12;
      C[2]->u[2] = PE->u[2];
      C[2]->v[2] = PE->v[2];
      C[2]->w[2] = PE->w[2];
      C[2]->u[3] = u23;
      C[2]->v[3] = v23;
      C[2]->w[3] = w23;

      C[3] = Create_PostElement(Index, QUADRANGLE, 4, PE->Depth);
      C[3]->u[0] = u03;
      C[3]->v[0] = v03;
      C[3]->w[0] = w03;
      C[3]->u[1] = u02;
      C[3]->v[1] = v02;
      C[3]->w[1] = w02;
      C[3]->u[2] = u23;
      C[3]->v[2] = v23;
      C[3]->w[2] = w23;
      C[3]->u[3] = PE->u[3];
      C[3]->v[3] = PE->v[3];
      C[3]->w[3] = PE->w[3];

      NbCut = 4;
      break;

    case TETRAHEDRON:
    case TETRAHEDRON_2:
    case TETRAHEDRON_3:
    case TETRAHEDRON_4:
      u01 = .5 * (PE->u[0] + PE->u[1]);
      u02 = .5 * (PE->u[0] + PE->u[2]);
      v01 = .5 * (PE->v[0] + PE->v[1]);
      v02 = .5 * (PE->v[0] + PE->v[2]);
      w01 = .5 * (PE->w[0] + PE->w[1]);
      w02 = .5 * (PE->w[0] + PE->w[2]);

      u03 = .5 * (PE->u[0] + PE->u[3]);
      u12 = .5 * (PE->u[1] + PE->u[2]);
      v03 = .5 * (PE->v[0] + PE->v[3]);
      v12 = .5 * (PE->v[1] + PE->v[2]);
      w03 = .5 * (PE->w[0] + PE->w[3]);
      w12 = .5 * (PE->w[1] + PE->w[2]);

      u13 = .5 * (PE->u[1] + PE->u[3]);
      u23 = .5 * (PE->u[2] + PE->u[3]);
      v13 = .5 * (PE->v[1] + PE->v[3]);
      v23 = .5 * (PE->v[2] + PE->v[3]);
      w13 = .5 * (PE->w[1] + PE->w[3]);
      w23 = .5 * (PE->w[2] + PE->w[3]);

      C[0] = Create_PostElement(Index, TETRAHEDRON, 4, PE->Depth);
      C[0]->u[0] = PE->u[0];
      C[0]->v[0] = PE->v[0];
      C[0]->w[0] = PE->w[0];
      C[0]->u[1] = u01;
      C[0]->v[1] = v01;
      C[0]->w[1] = w01;
      C[0]->u[2] = u02;
      C[0]->v[2] = v02;
      C[0]->w[2] = w02;
      C[0]->u[3] = u03;
      C[0]->v[3] = v03;
      C[0]->w[3] = w03;

      C[1] = Create_PostElement(Index, TETRAHEDRON, 4, PE->Depth);
      C[1]->u[0] = PE->u[1];
      C[1]->v[0] = PE->v[1];
      C[1]->w[0] = PE->w[1];
      C[1]->u[1] = u01;
      C[1]->v[1] = v01;
      C[1]->w[1] = w01;
      C[1]->u[2] = u12;
      C[1]->v[2] = v12;
      C[1]->w[2] = w12;
      C[1]->u[3] = u13;
      C[1]->v[3] = v13;
      C[1]->w[3] = w13;

      C[2] = Create_PostElement(Index, TETRAHEDRON, 4, PE->Depth);
      C[2]->u[0] = PE->u[2];
      C[2]->v[0] = PE->v[2];
      C[2]->w[0] = PE->w[2];
      C[2]->u[1] = u02;
      C[2]->v[1] = v02;
      C[2]->w[1] = w02;
      C[2]->u[2] = u12;
      C[2]->v[2] = v12;
      C[2]->w[2] = w12;
      C[2]->u[3] = u23;
      C[2]->v[3] = v23;
      C[2]->w[3] = w23;

      C[3] = Create_PostElement(Index, TETRAHEDRON, 4, PE->Depth);
      C[3]->u[0] = PE->u[3];
      C[3]->v[0] = PE->v[3];
      C[3]->w[0] = PE->w[3];
      C[3]->u[1] = u03;
      C[3]->v[1] = v03;
      C[3]->w[1] = w03;
      C[3]->u[2] = u13;
      C[3]->v[2] = v13;
      C[3]->w[2] = w13;
      C[3]->u[3] = u23;
      C[3]->v[3] = v23;
      C[3]->w[3] = w23;

      C[4] = Create_PostElement(Index, TETRAHEDRON, 4, PE->Depth);
      C[4]->u[0] = u01;
      C[4]->v[0] = v01;
      C[4]->w[0] = w01;
      C[4]->u[1] = u02;
      C[4]->v[1] = v02;
      C[4]->w[1] = w02;
      C[4]->u[2] = u03;
      C[4]->v[2] = v03;
      C[4]->w[2] = w03;
      C[4]->u[3] = u23;
      C[4]->v[3] = v23;
      C[4]->w[3] = w23;

      C[5] = Create_PostElement(Index, TETRAHEDRON, 4, PE->Depth);
      C[5]->u[0] = u01;
      C[5]->v[0] = v01;
      C[5]->w[0] = w01;
      C[5]->u[1] = u02;
      C[5]->v[1] = v02;
      C[5]->w[1] = w02;
      C[5]->u[2] = u12;
      C[5]->v[2] = v12;
      C[5]->w[2] = w12;
      C[5]->u[3] = u23;
      C[5]->v[3] = v23;
      C[5]->w[3] = w23;

      C[6] = Create_PostElement(Index, TETRAHEDRON, 4, PE->Depth);
      C[6]->u[0] = u01;
      C[6]->v[0] = v01;
      C[6]->w[0] = w01;
      C[6]->u[1] = u12;
      C[6]->v[1] = v12;
      C[6]->w[1] = w12;
      C[6]->u[2] = u13;
      C[6]->v[2] = v13;
      C[6]->w[2] = w13;
      C[6]->u[3] = u23;
      C[6]->v[3] = v23;
      C[6]->w[3] = w23;

      C[7] = PE;
      C[7]->u[0] = u01;
      C[7]->v[0] = v01;
      C[7]->w[0] = w01;
      C[7]->u[1] = u03;
      C[7]->v[1] = v03;
      C[7]->w[1] = w03;
      C[7]->u[2] = u13;
      C[7]->v[2] = v13;
      C[7]->w[2] = w13;
      C[7]->u[3] = u23;
      C[7]->v[3] = v23;
      C[7]->w[3] = w23;

      NbCut = 8;
      break;

    default:
      Message::Error(
        "Recursive division not implemented for Quadrangles, Hexahedra, "
        "Prisms and Pyramids");
      break;
    }

    for(i = 0; i < NbCut; i++) {
      C[i]->Depth++;
      for(j = 0; j < C[i]->NbrNodes; j++) C[i]->NumNodes[j] = -1;
      Cut_PostElement(C[i], GE, PE_L, Index, Depth, Skin, DecomposeInSimplex);
    }
  }
  else {
    Get_InitDofOfElement(&E);
    E.GeoElement = GE;
    E.Num = E.GeoElement->Num;
    E.Type = E.GeoElement->Type;
    E.Region = E.GeoElement->Region;
    Get_NodesCoordinatesOfElement(&E);

    for(i = 0; i < PE->NbrNodes; i++) {
      if(Skin == 0 && PE->Depth == 1 &&
         (DecomposeInSimplex == 0 || E.GeoElement->Type == LINE ||
          E.GeoElement->Type == TRIANGLE ||
          E.GeoElement->Type == TETRAHEDRON)) {
        PE->x[i] = E.x[i];
        PE->y[i] = E.y[i];
        PE->z[i] = E.z[i];
      }
      else {
        Get_BFGeoElement(&E, PE->u[i], PE->v[i], PE->w[i]);

        PE->x[i] = PE->y[i] = PE->z[i] = 0.;
        for(j = 0; j < E.GeoElement->NbrNodes; j++) {
          PE->x[i] += E.x[j] * E.n[j];
          PE->y[i] += E.y[j] * E.n[j];
          PE->z[i] += E.z[j] * E.n[j];
        }
      }
    }

    List_Add(PE_L, &PE);
  }
}

/* ------------------------------------------------------------------------ */
/*  F i l l _ P o s t E l e m e n t                                         */
/* ------------------------------------------------------------------------ */

static struct PostElement *Create_HighOrderPostElement(struct Geo_Element *GE,
                                                       int Index)
{
  struct PostElement *PE = NULL;
#if defined(HAVE_GMSH)
  const nodalBasis *basis = BasisFactory::getNodalBasis(GetDP2Gmsh(GE->Type));
  if(basis) {
    int NbrNodes = basis->points.size1();
    PE = Create_PostElement(Index, GE->Type, NbrNodes, 1);
    for(int i = 0; i < NbrNodes; i++) {
      PE->NumNodes[i] = GE->NumNodes[i];
      int dim = basis->points.size2();
      PE->u[i] = (dim >= 1) ? basis->points(i, 0) : 0;
      PE->v[i] = (dim >= 2) ? basis->points(i, 1) : 0;
      PE->w[i] = (dim >= 3) ? basis->points(i, 2) : 0;
    }
  }
#endif
  return PE;
}

#define POS_CUT_FILL                                                           \
  Cut_PostElement(PE, GE, PE_L, Index, Depth, 0, DecomposeInSimplex)
#define POS_CUT_SKIN                                                           \
  Cut_PostElement(PE, GE, PE_L, Index, Depth, 1, DecomposeInSimplex)

void Fill_PostElement(struct Geo_Element *GE, List_T *PE_L, int Index,
                      int Depth, int Skin, int DecomposeInSimplex,
                      int HighOrder, int Gauss)
{
  struct PostElement *PE;

  if(Gauss > 0) {
    Depth = 0;
    int error;
    void (*f)(int, int, double *, double *, double *, double *);
    Get_FunctionForDefine(FunctionForGauss, GE->Type, &error,
                          (void (**)()) & f);
    if(!error) {
      double dummy;
      for(int i = 0; i < Gauss; i++) {
        PE = Create_PostElement(Index, POINT_ELEMENT, 1, 0);
        f(Gauss, i, &PE->u[0], &PE->v[0], &PE->w[0], &dummy);
        POS_CUT_FILL;
      }
    }
  }
  else if(!Depth) {
    PE = Create_PostElement(Index, POINT_ELEMENT, 1, 0);
    switch(GE->Type) {
    case POINT_ELEMENT:
      PE->u[0] = 0.;
      PE->v[0] = 0.;
      PE->w[0] = 0.;
      break;
    case LINE:
    case LINE_2:
    case LINE_3:
    case LINE_4:
      PE->u[0] = 0.;
      PE->v[0] = 0.;
      PE->w[0] = 0.;
      break;
    case TRIANGLE:
    case TRIANGLE_2:
    case TRIANGLE_3:
    case TRIANGLE_4:
      PE->u[0] = 1. / 3.;
      PE->v[0] = 1. / 3.;
      PE->w[0] = 0.;
      break;
    case QUADRANGLE:
    case QUADRANGLE_2:
    case QUADRANGLE_2_8N:
    case QUADRANGLE_3:
    case QUADRANGLE_4:
      PE->u[0] = 0.;
      PE->v[0] = 0.;
      PE->w[0] = 0.;
      break;
    case TETRAHEDRON:
    case TETRAHEDRON_2:
    case TETRAHEDRON_3:
    case TETRAHEDRON_4:
      PE->u[0] = 0.25;
      PE->v[0] = 0.25;
      PE->w[0] = 0.25;
      break;
    case HEXAHEDRON:
    case HEXAHEDRON_2:
    case HEXAHEDRON_2_20N:
    case HEXAHEDRON_3:
    case HEXAHEDRON_4:
      PE->u[0] = 0.;
      PE->v[0] = 0.;
      PE->w[0] = 0.;
      break;
    case PRISM:
    case PRISM_2:
    case PRISM_2_15N:
    case PRISM_3:
    case PRISM_4:
      PE->u[0] = 1. / 3.;
      PE->v[0] = 1. / 3.;
      PE->w[0] = 0.;
      break;
    case PYRAMID:
    case PYRAMID_2:
    case PYRAMID_2_13N:
    case PYRAMID_3:
      // case PYRAMID_4 :
      PE->u[0] = 0.;
      PE->v[0] = 0.;
      PE->w[0] = 1. / 3.;
      break;
    }
    POS_CUT_FILL;
  }
  else {
    if(!Skin) {
      PE = NULL;

      switch(GE->Type) {
      case POINT_ELEMENT:
        PE = Create_PostElement(Index, POINT_ELEMENT, 1, 1); /* node 1 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->u[0] = 0.;
        PE->v[0] = 0.;
        PE->w[0] = 0.;
        POS_CUT_FILL;
        break;

      case LINE:
      case LINE_2:
      case LINE_3:
      case LINE_4:
        if(HighOrder && GE->Type != LINE)
          PE = Create_HighOrderPostElement(GE, Index);
        if(!PE) {
          PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 1 2 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[1];
          PE->u[0] = -1.;
          PE->v[0] = 0.;
          PE->w[0] = 0.;
          PE->u[1] = 1.;
          PE->v[1] = 0.;
          PE->w[1] = 0.;
        }
        POS_CUT_FILL;
        break;

      case TRIANGLE:
      case TRIANGLE_2:
      case TRIANGLE_3:
      case TRIANGLE_4:
        if(HighOrder && GE->Type != TRIANGLE)
          PE = Create_HighOrderPostElement(GE, Index);
        if(!PE) {
          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 2 3 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[1];
          PE->NumNodes[2] = GE->NumNodes[2];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = 0.;
          PE->u[1] = 1.;
          PE->v[1] = 0.;
          PE->w[1] = 0.;
          PE->u[2] = 0.;
          PE->v[2] = 1.;
          PE->w[2] = 0.;
        }
        POS_CUT_FILL;
        break;

      case QUADRANGLE:
      case QUADRANGLE_2:
      case QUADRANGLE_2_8N:
      case QUADRANGLE_3:
      case QUADRANGLE_4:
        if(HighOrder && GE->Type != QUADRANGLE)
          PE = Create_HighOrderPostElement(GE, Index);
        if(!PE) {
          if(DecomposeInSimplex) {
            PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 2 4 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[3];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = 0.;
            PE->u[1] = 1.;
            PE->v[1] = -1.;
            PE->w[1] = 0.;
            PE->u[2] = -1.;
            PE->v[2] = 1.;
            PE->w[2] = 0.;
            POS_CUT_FILL;
            PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 2 3 4 */
            PE->NumNodes[0] = GE->NumNodes[1];
            PE->NumNodes[1] = GE->NumNodes[2];
            PE->NumNodes[2] = GE->NumNodes[3];
            PE->u[0] = 1.;
            PE->v[0] = -1.;
            PE->w[0] = 0.;
            PE->u[1] = 1.;
            PE->v[1] = 1.;
            PE->w[1] = 0.;
            PE->u[2] = -1.;
            PE->v[2] = 1.;
            PE->w[2] = 0.;
          }
          else {
            PE =
              Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 1 2 3 4 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[2];
            PE->NumNodes[3] = GE->NumNodes[3];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = 0.;
            PE->u[1] = 1.;
            PE->v[1] = -1.;
            PE->w[1] = 0.;
            PE->u[2] = 1.;
            PE->v[2] = 1.;
            PE->w[2] = 0.;
            PE->u[3] = -1.;
            PE->v[3] = 1.;
            PE->w[3] = 0.;
          }
        }
        POS_CUT_FILL;
        break;

      case TETRAHEDRON:
      case TETRAHEDRON_2:
      case TETRAHEDRON_3:
      case TETRAHEDRON_4:
        if(HighOrder && GE->Type != TETRAHEDRON)
          PE = Create_HighOrderPostElement(GE, Index);
        if(!PE) {
          PE = Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 2 3 4 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[1];
          PE->NumNodes[2] = GE->NumNodes[2];
          PE->NumNodes[3] = GE->NumNodes[3];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = 0.;
          PE->u[1] = 1.;
          PE->v[1] = 0.;
          PE->w[1] = 0.;
          PE->u[2] = 0.;
          PE->v[2] = 1.;
          PE->w[2] = 0.;
          PE->u[3] = 0.;
          PE->v[3] = 0.;
          PE->w[3] = 1.;
        }
        POS_CUT_FILL;
        break;

      case HEXAHEDRON:
      case HEXAHEDRON_2:
      case HEXAHEDRON_2_20N:
      case HEXAHEDRON_3:
      case HEXAHEDRON_4:
        if(HighOrder && GE->Type != HEXAHEDRON)
          PE = Create_HighOrderPostElement(GE, Index);
        if(!PE) {
          if(DecomposeInSimplex) {
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 2 3 6 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[2];
            PE->NumNodes[3] = GE->NumNodes[5];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = -1.;
            PE->u[1] = 1.;
            PE->v[1] = -1.;
            PE->w[1] = -1.;
            PE->u[2] = 1.;
            PE->v[2] = 1.;
            PE->w[2] = -1.;
            PE->u[3] = 1.;
            PE->v[3] = -1.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 3 6 7 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[2];
            PE->NumNodes[2] = GE->NumNodes[5];
            PE->NumNodes[3] = GE->NumNodes[6];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = -1.;
            PE->u[1] = 1.;
            PE->v[1] = 1.;
            PE->w[1] = -1.;
            PE->u[2] = 1.;
            PE->v[2] = -1.;
            PE->w[2] = 1.;
            PE->u[3] = 1.;
            PE->v[3] = 1.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 5 6 7 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[4];
            PE->NumNodes[2] = GE->NumNodes[5];
            PE->NumNodes[3] = GE->NumNodes[6];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = -1.;
            PE->u[1] = -1.;
            PE->v[1] = -1.;
            PE->w[1] = 1.;
            PE->u[2] = 1.;
            PE->v[2] = -1.;
            PE->w[2] = 1.;
            PE->u[3] = 1.;
            PE->v[3] = 1.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 3 4 7 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[2];
            PE->NumNodes[2] = GE->NumNodes[3];
            PE->NumNodes[3] = GE->NumNodes[6];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = -1.;
            PE->u[1] = 1.;
            PE->v[1] = 1.;
            PE->w[1] = -1.;
            PE->u[2] = -1.;
            PE->v[2] = 1.;
            PE->w[2] = -1.;
            PE->u[3] = 1.;
            PE->v[3] = 1.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 5 7 8 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[4];
            PE->NumNodes[2] = GE->NumNodes[6];
            PE->NumNodes[3] = GE->NumNodes[7];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = -1.;
            PE->u[1] = -1.;
            PE->v[1] = -1.;
            PE->w[1] = 1.;
            PE->u[2] = 1.;
            PE->v[2] = 1.;
            PE->w[2] = 1.;
            PE->u[3] = -1.;
            PE->v[3] = 1.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 4 7 8 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[3];
            PE->NumNodes[2] = GE->NumNodes[6];
            PE->NumNodes[3] = GE->NumNodes[7];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = -1.;
            PE->u[1] = -1.;
            PE->v[1] = 1.;
            PE->w[1] = -1.;
            PE->u[2] = 1.;
            PE->v[2] = 1.;
            PE->w[2] = 1.;
            PE->u[3] = -1.;
            PE->v[3] = 1.;
            PE->w[3] = 1.;
          }
          else {
            PE = Create_PostElement(Index, HEXAHEDRON, 8,
                                    1); /* nodes 1 2 3 4 5 6 7 8 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[2];
            PE->NumNodes[3] = GE->NumNodes[3];
            PE->NumNodes[4] = GE->NumNodes[4];
            PE->NumNodes[5] = GE->NumNodes[5];
            PE->NumNodes[6] = GE->NumNodes[6];
            PE->NumNodes[7] = GE->NumNodes[7];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = -1.;
            PE->u[1] = 1.;
            PE->v[1] = -1.;
            PE->w[1] = -1.;
            PE->u[2] = 1.;
            PE->v[2] = 1.;
            PE->w[2] = -1.;
            PE->u[3] = -1.;
            PE->v[3] = 1.;
            PE->w[3] = -1.;
            PE->u[4] = -1.;
            PE->v[4] = -1.;
            PE->w[4] = 1.;
            PE->u[5] = 1.;
            PE->v[5] = -1.;
            PE->w[5] = 1.;
            PE->u[6] = 1.;
            PE->v[6] = 1.;
            PE->w[6] = 1.;
            PE->u[7] = -1.;
            PE->v[7] = 1.;
            PE->w[7] = 1.;
          }
        }
        POS_CUT_FILL;
        break;

      case PRISM:
      case PRISM_2:
      case PRISM_2_15N:
      case PRISM_3:
      case PRISM_4:
        if(HighOrder && GE->Type != PRISM)
          PE = Create_HighOrderPostElement(GE, Index);
        if(!PE) {
          if(DecomposeInSimplex) {
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 2 3 5 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[2];
            PE->NumNodes[3] = GE->NumNodes[4];
            PE->u[0] = 0.;
            PE->v[0] = 0.;
            PE->w[0] = -1.;
            PE->u[1] = 1.;
            PE->v[1] = 0.;
            PE->w[1] = -1.;
            PE->u[2] = 0.;
            PE->v[2] = 1.;
            PE->w[2] = -1.;
            PE->u[3] = 1.;
            PE->v[3] = 0.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 3 5 6 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[2];
            PE->NumNodes[2] = GE->NumNodes[4];
            PE->NumNodes[3] = GE->NumNodes[5];
            PE->u[0] = 0.;
            PE->v[0] = 0.;
            PE->w[0] = -1.;
            PE->u[1] = 0.;
            PE->v[1] = 1.;
            PE->w[1] = -1.;
            PE->u[2] = 1.;
            PE->v[2] = 0.;
            PE->w[2] = 1.;
            PE->u[3] = 0.;
            PE->v[3] = 1.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 4 5 6 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[3];
            PE->NumNodes[2] = GE->NumNodes[4];
            PE->NumNodes[3] = GE->NumNodes[5];
            PE->u[0] = 0.;
            PE->v[0] = 0.;
            PE->w[0] = -1.;
            PE->u[1] = 0.;
            PE->v[1] = 0.;
            PE->w[1] = 1.;
            PE->u[2] = 1.;
            PE->v[2] = 0.;
            PE->w[2] = 1.;
            PE->u[3] = 0.;
            PE->v[3] = 1.;
            PE->w[3] = 1.;
          }
          else {
            PE = Create_PostElement(Index, PRISM, 6, 1); /* nodes 1 2 3 4 5 6 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[2];
            PE->NumNodes[3] = GE->NumNodes[3];
            PE->NumNodes[4] = GE->NumNodes[4];
            PE->NumNodes[5] = GE->NumNodes[5];
            PE->u[0] = 0.;
            PE->v[0] = 0.;
            PE->w[0] = -1.;
            PE->u[1] = 1.;
            PE->v[1] = 0.;
            PE->w[1] = -1.;
            PE->u[2] = 0.;
            PE->v[2] = 1.;
            PE->w[2] = -1.;
            PE->u[3] = 0.;
            PE->v[3] = 0.;
            PE->w[3] = 1.;
            PE->u[4] = 1.;
            PE->v[4] = 0.;
            PE->w[4] = 1.;
            PE->u[5] = 0.;
            PE->v[5] = 1.;
            PE->w[5] = 1.;
          }
        }
        POS_CUT_FILL;
        break;

      case PYRAMID:
      case PYRAMID_2:
      case PYRAMID_2_13N:
      case PYRAMID_3: // case PYRAMID_4 :
        if(HighOrder && GE->Type != PYRAMID)
          PE = Create_HighOrderPostElement(GE, Index);
        if(!PE) {
          if(DecomposeInSimplex) {
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 1 2 4 5 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[3];
            PE->NumNodes[3] = GE->NumNodes[4];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = 0.;
            PE->u[1] = 1.;
            PE->v[1] = -1.;
            PE->w[1] = 0.;
            PE->u[2] = -1.;
            PE->v[2] = 1.;
            PE->w[2] = 0.;
            PE->u[3] = 0.;
            PE->v[3] = 0.;
            PE->w[3] = 1.;
            POS_CUT_FILL;
            PE =
              Create_PostElement(Index, TETRAHEDRON, 4, 1); /* nodes 2 3 4 5 */
            PE->NumNodes[0] = GE->NumNodes[1];
            PE->NumNodes[1] = GE->NumNodes[2];
            PE->NumNodes[2] = GE->NumNodes[3];
            PE->NumNodes[3] = GE->NumNodes[4];
            PE->u[0] = 1.;
            PE->v[0] = -1.;
            PE->w[0] = 0.;
            PE->u[1] = 1.;
            PE->v[1] = 1.;
            PE->w[1] = 0.;
            PE->u[2] = -1.;
            PE->v[2] = 1.;
            PE->w[2] = 0.;
            PE->u[3] = 0.;
            PE->v[3] = 0.;
            PE->w[3] = 1.;
          }
          else {
            PE = Create_PostElement(Index, PYRAMID, 5, 1); /* nodes 1 2 3 4 5 */
            PE->NumNodes[0] = GE->NumNodes[0];
            PE->NumNodes[1] = GE->NumNodes[1];
            PE->NumNodes[2] = GE->NumNodes[2];
            PE->NumNodes[3] = GE->NumNodes[3];
            PE->NumNodes[4] = GE->NumNodes[4];
            PE->u[0] = -1.;
            PE->v[0] = -1.;
            PE->w[0] = 0.;
            PE->u[1] = 1.;
            PE->v[1] = -1.;
            PE->w[1] = 0.;
            PE->u[2] = 1.;
            PE->v[2] = 1.;
            PE->w[2] = 0.;
            PE->u[3] = -1.;
            PE->v[3] = 1.;
            PE->w[3] = 0.;
            PE->u[4] = 0.;
            PE->v[4] = 0.;
            PE->w[4] = 1.;
          }
        }
        POS_CUT_FILL;
        break;
      }
    }
    else { /* Skin: facets oriented with normals pointing outwards */

      switch(GE->Type) {
      case TRIANGLE:
        PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 1 2 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->NumNodes[1] = GE->NumNodes[1];
        PE->u[0] = 0.;
        PE->v[0] = 0.;
        PE->w[0] = 0.;
        PE->u[1] = 1.;
        PE->v[1] = 0.;
        PE->w[1] = 0.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 2 3 */
        PE->NumNodes[0] = GE->NumNodes[1];
        PE->NumNodes[1] = GE->NumNodes[2];
        PE->u[0] = 1.;
        PE->v[0] = 0.;
        PE->w[0] = 0.;
        PE->u[1] = 0.;
        PE->v[1] = 1.;
        PE->w[1] = 0.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 3 1 */
        PE->NumNodes[0] = GE->NumNodes[2];
        PE->NumNodes[1] = GE->NumNodes[0];
        PE->u[0] = 0.;
        PE->v[0] = 1.;
        PE->w[0] = 0.;
        PE->u[1] = 0.;
        PE->v[1] = 0.;
        PE->w[1] = 0.;
        POS_CUT_SKIN;
        break;

      case QUADRANGLE:
        PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 1 2 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->NumNodes[1] = GE->NumNodes[1];
        PE->u[0] = -1.;
        PE->v[0] = -1.;
        PE->w[0] = 0.;
        PE->u[1] = 1.;
        PE->v[1] = -1.;
        PE->w[1] = 0.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 2 3 */
        PE->NumNodes[0] = GE->NumNodes[1];
        PE->NumNodes[1] = GE->NumNodes[2];
        PE->u[0] = 1.;
        PE->v[0] = -1.;
        PE->w[0] = 0.;
        PE->u[1] = 1.;
        PE->v[1] = 1.;
        PE->w[1] = 0.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 3 4 */
        PE->NumNodes[0] = GE->NumNodes[2];
        PE->NumNodes[1] = GE->NumNodes[3];
        PE->u[0] = 1.;
        PE->v[0] = 1.;
        PE->w[0] = 0.;
        PE->u[1] = -1.;
        PE->v[1] = 1.;
        PE->w[1] = 0.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, LINE, 2, 1); /* nodes 4 1 */
        PE->NumNodes[0] = GE->NumNodes[3];
        PE->NumNodes[1] = GE->NumNodes[0];
        PE->u[0] = -1.;
        PE->v[0] = 1.;
        PE->w[0] = 0.;
        PE->u[1] = -1.;
        PE->v[1] = -1.;
        PE->w[1] = 0.;
        POS_CUT_SKIN;
        break;

      case TETRAHEDRON:
        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 2 4 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->NumNodes[1] = GE->NumNodes[1];
        PE->NumNodes[2] = GE->NumNodes[3];
        PE->u[0] = 0.;
        PE->v[0] = 0.;
        PE->w[0] = 0.;
        PE->u[1] = 1.;
        PE->v[1] = 0.;
        PE->w[1] = 0.;
        PE->u[2] = 0.;
        PE->v[2] = 0.;
        PE->w[2] = 1.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 3 2 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->NumNodes[1] = GE->NumNodes[2];
        PE->NumNodes[2] = GE->NumNodes[1];
        PE->u[0] = 0.;
        PE->v[0] = 0.;
        PE->w[0] = 0.;
        PE->u[1] = 0.;
        PE->v[1] = 1.;
        PE->w[1] = 0.;
        PE->u[2] = 1.;
        PE->v[2] = 0.;
        PE->w[2] = 0.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 4 3 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->NumNodes[1] = GE->NumNodes[3];
        PE->NumNodes[2] = GE->NumNodes[2];
        PE->u[0] = 0.;
        PE->v[0] = 0.;
        PE->w[0] = 0.;
        PE->u[1] = 0.;
        PE->v[1] = 0.;
        PE->w[1] = 1.;
        PE->u[2] = 0.;
        PE->v[2] = 1.;
        PE->w[2] = 0.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 2 3 4 */
        PE->NumNodes[0] = GE->NumNodes[1];
        PE->NumNodes[1] = GE->NumNodes[2];
        PE->NumNodes[2] = GE->NumNodes[3];
        PE->u[0] = 1.;
        PE->v[0] = 0.;
        PE->w[0] = 0.;
        PE->u[1] = 0.;
        PE->v[1] = 1.;
        PE->w[1] = 0.;
        PE->u[2] = 0.;
        PE->v[2] = 0.;
        PE->w[2] = 1.;
        POS_CUT_SKIN;
        break;

      case HEXAHEDRON:
        if(DecomposeInSimplex) {
          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 4 2 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[1];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = -1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = -1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 2 4 3 */
          PE->NumNodes[0] = GE->NumNodes[1];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[2];
          PE->u[0] = 1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = -1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = 1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 5 6 8 */
          PE->NumNodes[0] = GE->NumNodes[4];
          PE->NumNodes[1] = GE->NumNodes[5];
          PE->NumNodes[2] = GE->NumNodes[7];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = 1.;
          PE->u[1] = 1.;
          PE->v[1] = -1.;
          PE->w[1] = 1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 6 7 8 */
          PE->NumNodes[0] = GE->NumNodes[5];
          PE->NumNodes[1] = GE->NumNodes[6];
          PE->NumNodes[2] = GE->NumNodes[7];
          PE->u[0] = 1.;
          PE->v[0] = -1.;
          PE->w[0] = 1.;
          PE->u[1] = 1.;
          PE->v[1] = 1.;
          PE->w[1] = 1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 5 4 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[4];
          PE->NumNodes[2] = GE->NumNodes[3];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = -1.;
          PE->v[1] = -1.;
          PE->w[1] = 1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 5 8 4 */
          PE->NumNodes[0] = GE->NumNodes[4];
          PE->NumNodes[1] = GE->NumNodes[7];
          PE->NumNodes[2] = GE->NumNodes[3];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = 1.;
          PE->u[1] = -1.;
          PE->v[1] = 1.;
          PE->w[1] = 1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 2 3 6 */
          PE->NumNodes[0] = GE->NumNodes[1];
          PE->NumNodes[1] = GE->NumNodes[2];
          PE->NumNodes[2] = GE->NumNodes[5];
          PE->u[0] = 1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = 1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = -1.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 6 3 7 */
          PE->NumNodes[0] = GE->NumNodes[5];
          PE->NumNodes[1] = GE->NumNodes[2];
          PE->NumNodes[2] = GE->NumNodes[6];
          PE->u[0] = 1.;
          PE->v[0] = -1.;
          PE->w[0] = 1.;
          PE->u[1] = 1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 5 1 6 */
          PE->NumNodes[0] = GE->NumNodes[4];
          PE->NumNodes[1] = GE->NumNodes[0];
          PE->NumNodes[2] = GE->NumNodes[5];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = 1.;
          PE->u[1] = -1.;
          PE->v[1] = -1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = -1.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 6 1 2 */
          PE->NumNodes[0] = GE->NumNodes[5];
          PE->NumNodes[1] = GE->NumNodes[0];
          PE->NumNodes[2] = GE->NumNodes[1];
          PE->u[0] = 1.;
          PE->v[0] = -1.;
          PE->w[0] = 1.;
          PE->u[1] = -1.;
          PE->v[1] = -1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = -1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 8 7 4 */
          PE->NumNodes[0] = GE->NumNodes[7];
          PE->NumNodes[1] = GE->NumNodes[6];
          PE->NumNodes[2] = GE->NumNodes[3];
          PE->u[0] = -1.;
          PE->v[0] = 1.;
          PE->w[0] = 1.;
          PE->u[1] = 1.;
          PE->v[1] = 1.;
          PE->w[1] = 1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 7 3 4 */
          PE->NumNodes[0] = GE->NumNodes[6];
          PE->NumNodes[1] = GE->NumNodes[2];
          PE->NumNodes[2] = GE->NumNodes[3];
          PE->u[0] = 1.;
          PE->v[0] = 1.;
          PE->w[0] = 1.;
          PE->u[1] = 1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;
        }
        else {
          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 1 2 6 5 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[1];
          PE->NumNodes[2] = GE->NumNodes[5];
          PE->NumNodes[3] = GE->NumNodes[4];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = 1.;
          PE->v[1] = -1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = -1.;
          PE->w[2] = 1.;
          PE->u[3] = -1.;
          PE->v[3] = -1.;
          PE->w[3] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes  1 4 3 2 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[2];
          PE->NumNodes[3] = GE->NumNodes[1];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = -1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = 1.;
          PE->w[2] = -1.;
          PE->u[3] = 1.;
          PE->v[3] = -1.;
          PE->w[3] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 1 5 8 4 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[4];
          PE->NumNodes[2] = GE->NumNodes[7];
          PE->NumNodes[3] = GE->NumNodes[3];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = -1.;
          PE->v[1] = -1.;
          PE->w[1] = 1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          PE->u[3] = -1.;
          PE->v[3] = 1.;
          PE->w[3] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 2 3 7 6 */
          PE->NumNodes[0] = GE->NumNodes[1];
          PE->NumNodes[1] = GE->NumNodes[2];
          PE->NumNodes[2] = GE->NumNodes[6];
          PE->NumNodes[3] = GE->NumNodes[5];
          PE->u[0] = 1.;
          PE->v[0] = -1.;
          PE->w[0] = -1.;
          PE->u[1] = 1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          PE->u[3] = 1.;
          PE->v[3] = -1.;
          PE->w[3] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 3 4 8 7 */
          PE->NumNodes[0] = GE->NumNodes[2];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[7];
          PE->NumNodes[3] = GE->NumNodes[6];
          PE->u[0] = 1.;
          PE->v[0] = 1.;
          PE->w[0] = -1.;
          PE->u[1] = -1.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = -1.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          PE->u[3] = 1.;
          PE->v[3] = 1.;
          PE->w[3] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 5 6 7 8 */
          PE->NumNodes[0] = GE->NumNodes[4];
          PE->NumNodes[1] = GE->NumNodes[5];
          PE->NumNodes[2] = GE->NumNodes[6];
          PE->NumNodes[3] = GE->NumNodes[7];
          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = 1.;
          PE->u[1] = 1.;
          PE->v[1] = -1.;
          PE->w[1] = 1.;
          PE->u[2] = 1.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          PE->u[3] = -1.;
          PE->v[3] = 1.;
          PE->w[3] = 1.;
          POS_CUT_SKIN;
        }
        break;

      case PRISM:
        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 3 2 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->NumNodes[1] = GE->NumNodes[2];
        PE->NumNodes[2] = GE->NumNodes[1];
        PE->u[0] = 0.;
        PE->v[0] = 0.;
        PE->w[0] = -1.;
        PE->u[1] = 0.;
        PE->v[1] = 1.;
        PE->w[1] = -1.;
        PE->u[2] = 1.;
        PE->v[2] = 0.;
        PE->w[2] = -1.;
        POS_CUT_SKIN;

        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 4 5 6 */
        PE->NumNodes[0] = GE->NumNodes[3];
        PE->NumNodes[1] = GE->NumNodes[4];
        PE->NumNodes[2] = GE->NumNodes[5];
        PE->u[0] = 0.;
        PE->v[0] = 0.;
        PE->w[0] = 1.;
        PE->u[1] = 1.;
        PE->v[1] = 0.;
        PE->w[1] = 1.;
        PE->u[2] = 0.;
        PE->v[2] = 1.;
        PE->w[2] = 1.;
        POS_CUT_SKIN;

        if(DecomposeInSimplex) {
          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 2 5 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[1];
          PE->NumNodes[2] = GE->NumNodes[4];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 1.;
          PE->v[1] = 0.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = 0.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 5 4 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[4];
          PE->NumNodes[2] = GE->NumNodes[3];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 1.;
          PE->v[1] = 0.;
          PE->w[1] = 1.;
          PE->u[2] = 0.;
          PE->v[2] = 0.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 6 3 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[5];
          PE->NumNodes[2] = GE->NumNodes[2];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 0.;
          PE->v[1] = 1.;
          PE->w[1] = 1.;
          PE->u[2] = 0.;
          PE->v[2] = 1.;
          PE->w[2] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 4 6 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[5];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 0.;
          PE->v[1] = 0.;
          PE->w[1] = 1.;
          PE->u[2] = 0.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 2 3 5 */
          PE->NumNodes[0] = GE->NumNodes[1];
          PE->NumNodes[1] = GE->NumNodes[2];
          PE->NumNodes[2] = GE->NumNodes[4];
          PE->u[0] = 1.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 0.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = 0.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 3 6 5 */
          PE->NumNodes[0] = GE->NumNodes[2];
          PE->NumNodes[1] = GE->NumNodes[5];
          PE->NumNodes[2] = GE->NumNodes[4];
          PE->u[0] = 0.;
          PE->v[0] = 1.;
          PE->w[0] = -1.;
          PE->u[1] = 0.;
          PE->v[1] = 1.;
          PE->w[1] = 1.;
          PE->u[2] = 1.;
          PE->v[2] = 0.;
          PE->w[2] = 1.;
          POS_CUT_SKIN;
        }
        else {
          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 1 2 5 4 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[1];
          PE->NumNodes[2] = GE->NumNodes[4];
          PE->NumNodes[3] = GE->NumNodes[3];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 1.;
          PE->v[1] = 0.;
          PE->w[1] = -1.;
          PE->u[2] = 1.;
          PE->v[2] = 0.;
          PE->w[2] = 1.;
          PE->u[3] = 0.;
          PE->v[3] = 0.;
          PE->w[3] = 1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 1 4 6 3 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[5];
          PE->NumNodes[3] = GE->NumNodes[2];
          PE->u[0] = 0.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 0.;
          PE->v[1] = 0.;
          PE->w[1] = 1.;
          PE->u[2] = 0.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          PE->u[3] = 0.;
          PE->v[3] = 1.;
          PE->w[3] = -1.;
          POS_CUT_SKIN;

          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 2 3 6 5 */
          PE->NumNodes[0] = GE->NumNodes[1];
          PE->NumNodes[1] = GE->NumNodes[2];
          PE->NumNodes[2] = GE->NumNodes[5];
          PE->NumNodes[3] = GE->NumNodes[4];
          PE->u[0] = 1.;
          PE->v[0] = 0.;
          PE->w[0] = -1.;
          PE->u[1] = 0.;
          PE->v[1] = 1.;
          PE->w[1] = -1.;
          PE->u[2] = 0.;
          PE->v[2] = 1.;
          PE->w[2] = 1.;
          PE->u[3] = 1.;
          PE->v[3] = 0.;
          PE->w[3] = 1.;
          POS_CUT_SKIN;
        }
        break;

      case PYRAMID:
        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 2 5 */
        PE->NumNodes[0] = GE->NumNodes[0];
        PE->NumNodes[1] = GE->NumNodes[1];
        PE->NumNodes[2] = GE->NumNodes[4];

        PE->u[0] = -1.;
        PE->v[0] = -1.;
        PE->w[0] = 0.;
        PE->u[1] = 1.;
        PE->v[1] = -1.;
        PE->w[1] = 0.;
        PE->u[2] = 0.;
        PE->v[2] = 0.;
        PE->w[2] = 1.;

        POS_CUT_SKIN;

        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 2 3 5 */
        PE->NumNodes[0] = GE->NumNodes[1];
        PE->NumNodes[1] = GE->NumNodes[2];
        PE->NumNodes[2] = GE->NumNodes[4];

        PE->u[0] = 1.;
        PE->v[0] = -1.;
        PE->w[0] = 0.;
        PE->u[1] = 1.;
        PE->v[1] = 1.;
        PE->w[1] = 0.;
        PE->u[2] = 0.;
        PE->v[2] = 0.;
        PE->w[2] = 1.;

        POS_CUT_SKIN;

        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 3 4 5 */
        PE->NumNodes[0] = GE->NumNodes[2];
        PE->NumNodes[1] = GE->NumNodes[3];
        PE->NumNodes[2] = GE->NumNodes[4];

        PE->u[0] = 1.;
        PE->v[0] = 1.;
        PE->w[0] = 0.;
        PE->u[1] = -1.;
        PE->v[1] = 1.;
        PE->w[1] = 0.;
        PE->u[2] = 0.;
        PE->v[2] = 0.;
        PE->w[2] = 1.;

        POS_CUT_SKIN;

        PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 4 1 5 */
        PE->NumNodes[0] = GE->NumNodes[3];
        PE->NumNodes[1] = GE->NumNodes[0];
        PE->NumNodes[2] = GE->NumNodes[4];

        PE->u[0] = -1.;
        PE->v[0] = 1.;
        PE->w[0] = 0.;
        PE->u[1] = -1.;
        PE->v[1] = -1.;
        PE->w[1] = 0.;
        PE->u[2] = 0.;
        PE->v[2] = 0.;
        PE->w[2] = 1.;

        POS_CUT_SKIN;

        if(DecomposeInSimplex) {
          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 3 2 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[2];
          PE->NumNodes[2] = GE->NumNodes[1];

          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = 0.;
          PE->u[1] = 1.;
          PE->v[1] = 1.;
          PE->w[1] = 0.;
          PE->u[2] = 1.;
          PE->v[2] = -1.;
          PE->w[2] = 0.;

          POS_CUT_SKIN;

          PE = Create_PostElement(Index, TRIANGLE, 3, 1); /* nodes 1 4 3 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[2];

          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = 0.;
          PE->u[1] = -1.;
          PE->v[1] = 1.;
          PE->w[1] = 0.;
          PE->u[2] = 1.;
          PE->v[2] = 1.;
          PE->w[2] = 0.;

          POS_CUT_SKIN;
        }
        else {
          PE = Create_PostElement(Index, QUADRANGLE, 4, 1); /* nodes 1 4 3 2 */
          PE->NumNodes[0] = GE->NumNodes[0];
          PE->NumNodes[1] = GE->NumNodes[3];
          PE->NumNodes[2] = GE->NumNodes[2];
          PE->NumNodes[3] = GE->NumNodes[1];

          PE->u[0] = -1.;
          PE->v[0] = -1.;
          PE->w[0] = 0.;
          PE->u[1] = -1.;
          PE->v[1] = 1.;
          PE->w[1] = 0.;
          PE->u[2] = 1.;
          PE->v[2] = 1.;
          PE->w[2] = 0.;
          PE->u[3] = 1.;
          PE->v[3] = -1.;
          PE->w[3] = 0.;

          POS_CUT_SKIN;
        }
        break;
      }
    }
  }
}

#undef POS_CUT_FILL
#undef POS_CUT_SKIN

/* ------------------------------------------------------------------------ */
/*  S o r t B y C o n n e c t i v i t y                                     */
/* ------------------------------------------------------------------------ */

int Compare_PostElement_Node(struct PostElement *PE1, int n1,
                             struct PostElement *PE2, int n2)
{
  double TOL = Current.GeoData->CharacteristicLength * 1.e-8;
  if((fabs(PE1->x[n1] - PE2->x[n2]) < TOL) &&
     (fabs(PE1->y[n1] - PE2->y[n2]) < TOL) &&
     (fabs(PE1->z[n1] - PE2->z[n2]) < TOL))
    return 1;
  return 0;
}

void Sort_PostElement_Connectivity(List_T *PostElement_L)
{
  int ii, jj, start, end, iPost, NbrPost;
  struct PostElement *PE, *PE2;

  NbrPost = List_Nbr(PostElement_L);

  /*
     u[0] = 1 if the element is in the ordered list, with natural orientation
           -1 if the element is in the ordered list, but with opposite
     orientation 0 if the element is not in the list v[0] = relative index (to
     the first element) in the ordered list
  */

  for(ii = 0; ii < NbrPost; ii++) {
    PE = *(struct PostElement **)List_Pointer(PostElement_L, ii);
    if(PE->NbrNodes != 2) {
      Message::Error("Connectivity sorting impossible for %d-noded elements",
                     PE->NbrNodes);
      return;
    }
    PE->u[0] = 0.;
  }

  PE = *(struct PostElement **)List_Pointer(PostElement_L, 0);
  PE->u[0] = 1.;
  PE->v[0] = 0.;

  iPost = 1;
  while(iPost < NbrPost) {
    for(ii = 0; ii < NbrPost; ii++) {
      PE = *(struct PostElement **)List_Pointer(PostElement_L, ii);
      if(PE->u[0]) {
        if(PE->u[0] > 0) {
          start = 0;
          end = 1;
        }
        else {
          start = 1;
          end = 0;
        }
        for(jj = 0; jj < NbrPost; jj++) {
          if(jj != ii) {
            PE2 = *(struct PostElement **)List_Pointer(PostElement_L, jj);
            if(!PE2->u[0]) {
              if(Compare_PostElement_Node(PE, end, PE2, 0)) {
                PE2->u[0] = 1.;
                PE2->v[0] = PE->v[0] + 1.;
                iPost++;
              }
              else if(Compare_PostElement_Node(PE, start, PE2, 0)) {
                PE2->u[0] = -1.;
                PE2->v[0] = PE->v[0] - 1.;
                iPost++;
              }
              else if(Compare_PostElement_Node(PE, start, PE2, 1)) {
                PE2->u[0] = 1.;
                PE2->v[0] = PE->v[0] - 1.;
                iPost++;
              }
              else if(Compare_PostElement_Node(PE, end, PE2, 1)) {
                PE2->u[0] = -1.;
                PE2->v[0] = PE->v[0] + 1.;
                iPost++;
              }
            }
          }
        }
      }
    }
    List_Sort(PostElement_L, fcmp_PostElement_absu0);
  }
  List_Sort(PostElement_L, fcmp_PostElement_v0);
}
