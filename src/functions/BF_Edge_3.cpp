// GetDP - Copyright (C) 1997-2022 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.

#include <stdlib.h>
#include "GetDPConfig.h"
#include "ProData.h"
#include "Message.h"

#if defined(HAVE_KERNEL)
#include "GeoData.h"
#endif

/* The non-symmetric facet functions are selected according to the
   NumIndex^th smallest global node number */

int fcmp_Int2(const void *a, const void *b)
{
  return ((struct TwoInt *)a)->Int2 - ((struct TwoInt *)b)->Int2;
}

int Get_FacetFunctionIndex(struct Element *Element, int NumEntity, int NumIndex)
{
#if !defined(HAVE_KERNEL)
  Message::Error("Get_FacetFunctionIndex requires Kernel");
  return 0;
#else

  int i, j, *NumNodes;

  if(Element->NumLastElementForSortedNodesByFacet != Element->Num) {
    for(i = 0; i < Element->GeoElement->NbrFacets; i++) {
      NumNodes = Geo_GetNodesOfFacetInElement(Element->GeoElement, i);
      j = 0;
      while(NumNodes[j]) {
        Element->SortedNodesByFacet[i][j].Int1 = NumNodes[j];
        Element->SortedNodesByFacet[i][j].Int2 =
          Element->GeoElement->NumNodes[NumNodes[j] - 1];
        j++;
      }
      qsort(Element->SortedNodesByFacet[i], j, sizeof(struct TwoInt),
            fcmp_Int2);
    }

    Element->NumLastElementForSortedNodesByFacet = Element->Num;
  }

  return Element->SortedNodesByFacet[NumEntity - 1][NumIndex - 1].Int1;
#endif
}

/* ------------------------------------------------------------------------ */
/*  B F _ E d g e _ 3                                                       */
/* ------------------------------------------------------------------------ */

/* ------- */
/*  Edges  */
/* ------- */

#define WrongNumEntity Message::Error("Wrong Edge number in 'BF_Edge_3E'")

void BF_Edge_3E(struct Element *Element, int NumEntity, double u, double v,
                double w, double s[])
{
  Message::Error("You should never end up here!");
}

#undef WrongNumEntity

/* -------- */
/*  Facets  */
/* -------- */

#define WrongNumEntity Message::Error("Wrong Face number in 'BF_Edge_3F'")

void BF_Edge_3F(struct Element *Element, int NumEntity, int Index, double u,
                double v, double w, double s[])
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4: Message::Error("You should never end up here!"); break;

  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
    switch(NumEntity) {
    case 1:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 3:
        s[0] = 0.;
        s[1] = (1 - u - v) * u;
        s[2] = 0.;
        break;
      case 1:
        s[0] = -u * v;
        s[1] = -u * v;
        s[2] = 0.;
        break;
      case 2:
        s[0] = (1 - u - v) * v;
        s[1] = 0.;
        s[2] = 0.;
        break;
      }
      break;
    default: WrongNumEntity;
    }
    break;

  case QUADRANGLE:
  case QUADRANGLE_2:
  case QUADRANGLE_2_8N:
  case QUADRANGLE_3:
  case QUADRANGLE_4:
    switch(NumEntity) {
    case 1:
      // WARNING: this should be generalized for 3D elements (hexahedrons), by
      // using the Get_FacetFunctionIndex function the 2D basis functions should
      // be the restriction of the 3D basis functions on each facet for now, the
      // BF_Edge_3F_c function is not implemented
      // switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      switch(Index) {
      case 1: // BF_Edge_3F_a
        s[0] = 45. / 16. * (1 - u) * (1 - v * v);
        s[1] = 45. / 16. * (1 - v) * (1 - u * u);
        s[2] = 0.;
        break;
      case 2: // BF_Edge_3F_b
        s[0] = 45. / 16. * (1 + u) * (v * v - 1);
        s[1] = 45. / 16. * (1 - v) * (1 - u * u);
        s[2] = 0.;
        break;
      case 3: // BF_Edge_3F_c
        Message::Error("You should never end up here!");
        break;
      }
      break;
    default: WrongNumEntity;
    }
    break;

  case TETRAHEDRON:
  case TETRAHEDRON_2:
  case TETRAHEDRON_3:
  case TETRAHEDRON_4:
    switch(NumEntity) {
    case 1:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 4:
        s[0] = 0.;
        s[1] = 0.;
        s[2] = u * (1. - u - v - w);
        break;
      case 1:
        s[0] = -u * w;
        s[1] = -u * w;
        s[2] = -u * w;
        break;
      case 2:
        s[0] = (1. - u - v - w) * w;
        s[1] = 0.;
        s[2] = 0.;
        break;
      }
      break;
    case 2:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 2:
        s[0] = v * (1 - u - v - w);
        s[1] = 0.;
        s[2] = 0.;
        break;
      case 1:
        s[0] = -u * v;
        s[1] = -u * v;
        s[2] = -u * v;
        break;
      case 3:
        s[0] = 0.;
        s[1] = u * (1. - u - v - w);
        s[2] = 0.;
        break;
      }
      break;
    case 3:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 3:
        s[0] = 0.;
        s[1] = (1. - u - v - w) * w;
        s[2] = 0.;
        break;
      case 1:
        s[0] = -v * w;
        s[1] = -v * w;
        s[2] = -v * w;
        break;
      case 4:
        s[0] = 0.;
        s[1] = 0.;
        s[2] = v * (1. - u - v - w);
        break;
      }
      break;
    case 4:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 4:
        s[0] = 0.;
        s[1] = 0.;
        s[2] = u * v;
        break;
      case 2:
        s[0] = v * w;
        s[1] = 0.;
        s[2] = 0.;
        break;
      case 3:
        s[0] = 0.;
        s[1] = u * w;
        s[2] = 0.;
        break;
      }
      break;
    default: WrongNumEntity;
    }
    break;

  case HEXAHEDRON:
  case HEXAHEDRON_2:
  case HEXAHEDRON_2_20N:
  case HEXAHEDRON_3:
  case HEXAHEDRON_4:
    switch(NumEntity) {
    default: Message::Error("BF_Edge_3F not ready for HEXAHEDRON");
    }
    break;

  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
    switch(NumEntity) {
    default: Message::Error("BF_Edge_3F not ready for PRISM");
    }
    break;

  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    switch(NumEntity) {
    default: Message::Error("BF_Edge_3F not ready for PYRAMID");
    }
    break;

  default: Message::Error("Unknown type of Element in BF_Edge_3F"); break;
  }
}

#undef WrongNumEntity

void BF_Edge_3F_a(struct Element *Element, int NumEntity, double u, double v,
                  double w, double s[])
{
  BF_Edge_3F(Element, NumEntity, 1, u, v, w, s);
}

void BF_Edge_3F_b(struct Element *Element, int NumEntity, double u, double v,
                  double w, double s[])
{
  BF_Edge_3F(Element, NumEntity, 2, u, v, w, s);
}

void BF_Edge_3F_c(struct Element *Element, int NumEntity, double u, double v,
                  double w, double s[])
{
  BF_Edge_3F(Element, NumEntity, 3, u, v, w, s);
}

/* -------- */
/*  Volume  */
/* -------- */

void BF_Edge_3V(struct Element *Element, int NumEntity, double u, double v,
                double w, double s[])
{
  Message::Error("You should never end up here!");
}

/* ------------------------------------------------------------------------ */
/*  B F _ C u r l E d g e _ 3                                               */
/* ------------------------------------------------------------------------ */

/* ------- */
/*  Edges  */
/* ------- */

#define WrongNumEntity Message::Error("Wrong Edge number in 'BF_CurlEdge_3E'")

void BF_CurlEdge_3E(struct Element *Element, int NumEntity, double u, double v,
                    double w, double s[])
{
  Message::Error("You should never end up here!");
}

#undef WrongNumEntity

/* -------- */
/*  Facets  */
/* -------- */

#define WrongNumEntity Message::Error("Wrong Face number in 'BF_CurlEdge_3F'")

void BF_CurlEdge_3F(struct Element *Element, int NumEntity, int Index, double u,
                    double v, double w, double s[])
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4: Message::Error("You should never end up here!"); break;

  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
    switch(NumEntity) {
    case 1:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 3:
        s[0] = 0.;
        s[1] = 0.;
        s[2] = -2.0 * u + 1.0 - v;
        break;
      case 1:
        s[0] = 0.;
        s[1] = 0.;
        s[2] = -v + u;
        break;
      case 2:
        s[0] = 0.;
        s[1] = 0.;
        s[2] = 2.0 * v - 1.0 + u;
        break;
      }
      break;
    default: WrongNumEntity;
    }
    break;

  case QUADRANGLE:
  case QUADRANGLE_2:
  case QUADRANGLE_2_8N:
  case QUADRANGLE_3:
  case QUADRANGLE_4:
    switch(NumEntity) {
    case 1:
      // WARNING: this should be generalized for 3D elements (hexahedrons), by
      // using the Get_FacetFunctionIndex function the 2D basis functions should
      // be the restriction of the 3D basis functions on each facet for now, the
      // BF_CurlEdge_3F_c function is not implemented
      // switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      switch(Index) {
      case 1: // BF_CurlEdge_3F_a
        s[0] = 0.;
        s[1] = 0.;
        s[2] = 45 / 8 * (v - u);
        break;
      case 2: // BF_CurlEdge_3F_b
        s[0] = 0.;
        s[1] = 0.;
        s[2] = -45 / 8 * (v + u);
        break;
      case 3: // BF_CurlEdge_3F_c
        Message::Error("You should never end up here!");
        break;
      }
      break;
    default: WrongNumEntity;
    }
    break;

  case TETRAHEDRON:
  case TETRAHEDRON_2:
  case TETRAHEDRON_3:
  case TETRAHEDRON_4:
    switch(NumEntity) {
    case 1:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 4:
        s[0] = -u;
        s[1] = -1. + 2. * u + v + w;
        s[2] = 0.;
        break;
      case 1:
        s[0] = u;
        s[1] = -u + w;
        s[2] = -w;
        break;
      case 2:
        s[0] = 0.;
        s[1] = 1. - u - v - 2. * w;
        s[2] = w;
        break;
      }
      break;
    case 2:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 2:
        s[0] = 0.;
        s[1] = -v;
        s[2] = -1. + u + 2. * v + w;
        break;
      case 1:
        s[0] = -u;
        s[1] = v;
        s[2] = u - v;
        break;
      case 3:
        s[0] = u;
        s[1] = 0.;
        s[2] = 1. - 2. * u - v - w;
        break;
      }
      break;
    case 3:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 3:
        s[0] = -1. + u + v + 2. * w;
        s[1] = 0.;
        s[2] = -w;
        break;
      case 1:
        s[0] = v - w;
        s[1] = -v;
        s[2] = w;
        break;
      case 4:
        s[0] = 1. - u - 2. * v - w;
        s[1] = v;
        s[2] = 0.;
        break;
      }
      break;
    case 4:
      switch(Get_FacetFunctionIndex(Element, NumEntity, Index)) {
      case 4:
        s[0] = u;
        s[1] = -v;
        s[2] = 0.;
        break;
      case 2:
        s[0] = 0.;
        s[1] = v;
        s[2] = -w;
        break;
      case 3:
        s[0] = -u;
        s[1] = 0.;
        s[2] = w;
        break;
      }
      break;
    default: WrongNumEntity;
    }
    break;

  case HEXAHEDRON:
  case HEXAHEDRON_2:
  case HEXAHEDRON_2_20N:
  case HEXAHEDRON_3:
  case HEXAHEDRON_4:
    switch(NumEntity) {
    default: Message::Error("BF_CurlEdge_3F not ready for HEXAHEDRON");
    }
    break;

  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
    switch(NumEntity) {
    default: Message::Error("BF_CurlEdge_3F not ready for PRISM");
    }
    break;

  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    switch(NumEntity) {
    default: Message::Error("BF_CurlEdge_3F not ready for PYRAMID");
    }
    break;

  default: Message::Error("Unknown type of Element in BF_CurlEdge_3F"); break;
  }
}

#undef WrongNumEntity

void BF_CurlEdge_3F_a(struct Element *Element, int NumEntity, double u,
                      double v, double w, double s[])
{
  BF_CurlEdge_3F(Element, NumEntity, 1, u, v, w, s);
}

void BF_CurlEdge_3F_b(struct Element *Element, int NumEntity, double u,
                      double v, double w, double s[])
{
  BF_CurlEdge_3F(Element, NumEntity, 2, u, v, w, s);
}

void BF_CurlEdge_3F_c(struct Element *Element, int NumEntity, double u,
                      double v, double w, double s[])
{
  BF_CurlEdge_3F(Element, NumEntity, 3, u, v, w, s);
}

/* -------- */
/*  Volume  */
/* -------- */

void BF_CurlEdge_3V(struct Element *Element, int NumEntity, double u, double v,
                    double w, double s[])
{
  Message::Error("You should never end up here!");
}
