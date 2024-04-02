// GetDP - Copyright (C) 1997-2022 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.

#include "ProData.h"
#include "Message.h"

/* ------------------------------------------------------------------------ */
/*  B F _ N o d e _ 2                                                       */
/* ------------------------------------------------------------------------ */

/* ------- */
/*  Edges  */
/* ------- */

#define WrongNumEntity Message::Error("Wrong Edge number in 'BF_Node_2E'")

void BF_Node_2E(struct Element *Element, int NumEntity, double u, double v,
                double w, double *s)
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4:
    switch(NumEntity) {
    case 1: *s = 0.25 * (1. - u) * (1. + u); break;
    default: WrongNumEntity;
    }
    break;

  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
    switch(NumEntity) {
    case 1: *s = (1. - u - v) * u; break;
    case 2: *s = (1. - u - v) * v; break;
    case 3: *s = u * v; break;
    default: WrongNumEntity;
    }
    break;

  case QUADRANGLE:
  case QUADRANGLE_2:
  case QUADRANGLE_2_8N:
  case QUADRANGLE_3:
  case QUADRANGLE_4:
    switch(NumEntity) {
    case 1: *s = 0.0625 * (1. - u * u) * (1. - v) * (1. - v); break;
    case 2: *s = 0.0625 * (1. - u) * (1. - u) * (1. - v * v); break;
    case 3: *s = 0.0625 * (1. + u) * (1. + u) * (1. - v * v); break;
    case 4: *s = 0.0625 * (1. - u * u) * (1. + v) * (1. + v); break;
    default: WrongNumEntity;
    }
    break;

  case TETRAHEDRON:
  case TETRAHEDRON_2:
  case TETRAHEDRON_3:
  case TETRAHEDRON_4:
    switch(NumEntity) {
    case 1: *s = (1. - u - v - w) * u; break;
    case 2: *s = (1. - u - v - w) * v; break;
    case 3: *s = (1. - u - v - w) * w; break;
    case 4: *s = u * v; break;
    case 5: *s = u * w; break;
    case 6: *s = v * w; break;
    default: WrongNumEntity;
    }
    break;

  case HEXAHEDRON:
  case HEXAHEDRON_2:
  case HEXAHEDRON_2_20N:
  case HEXAHEDRON_3:
  case HEXAHEDRON_4:
    switch(NumEntity) {
    case 1:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. - w) * (1. + u) * (1. - v) *
           (1. - w);
      break;
    case 2:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. - w) * (1. - u) * (1. + v) *
           (1. - w);
      break;
    case 3:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. - w) * (1. - u) * (1. - v) *
           (1. + w);
      break;
    case 4:
      *s = 0.015625 * (1. + u) * (1. - v) * (1. - w) * (1. + u) * (1. + v) *
           (1. - w);
      break;
    case 5:
      *s = 0.015625 * (1. + u) * (1. - v) * (1. - w) * (1. + u) * (1. - v) *
           (1. + w);
      break;
    case 6:
      *s = 0.015625 * (1. + u) * (1. + v) * (1. - w) * (1. - u) * (1. + v) *
           (1. - w);
      break;
    case 7:
      *s = 0.015625 * (1. + u) * (1. + v) * (1. - w) * (1. + u) * (1. + v) *
           (1. + w);
      break;
    case 8:
      *s = 0.015625 * (1. - u) * (1. + v) * (1. - w) * (1. - u) * (1. + v) *
           (1. + w);
      break;
    case 9:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. + w) * (1. + u) * (1. - v) *
           (1. + w);
      break;
    case 10:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. + w) * (1. - u) * (1. + v) *
           (1. + w);
      break;
    case 11:
      *s = 0.015625 * (1. + u) * (1. - v) * (1. + w) * (1. + u) * (1. + v) *
           (1. + w);
      break;
    case 12:
      *s = 0.015625 * (1. + u) * (1. + v) * (1. + w) * (1. - u) * (1. + v) *
           (1. + w);
      break;
    default: WrongNumEntity;
    }
    break;

  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
    // FIXME: not tested!
    switch(NumEntity) {
    case 1: *s = 0.25 * (1. - u - v) * (1. - w) * u * (1. - w); break;
    case 2: *s = 0.25 * (1. - u - v) * (1. - w) * v * (1. - w); break;
    case 3:
      *s = 0.25 * (1. - u - v) * (1. - w) * (1. - u - v) * (1. + w);
      break;
    case 4: *s = 0.25 * u * (1. - w) * v * (1. - w); break;
    case 5: *s = 0.25 * u * (1. - w) * u * (1. + w); break;
    case 6: *s = 0.25 * v * (1. - w) * v * (1. + w); break;
    case 7: *s = 0.25 * (1. - u - v) * (1. + w) * u * (1. + w); break;
    case 8: *s = 0.25 * (1. - u - v) * (1. + w) * v * (1. + w); break;
    case 9: *s = 0.25 * u * (1. + w) * v * (1. + w); break;
    default: WrongNumEntity;
    }
    break;

  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    switch(NumEntity) {
    default: Message::Error("BF_Node_2E not ready for PYRAMID");
    }
    break;

  default: Message::Error("Unknown type of Element in BF_Node_2E"); break;
  }
}

#undef WrongNumEntity

/* -------- */
/*  Facets  */
/* -------- */

#define WrongNumEntity Message::Error("Wrong Face number in 'BF_Node_2F'")

void BF_Node_2F(struct Element *Element, int NumEntity, double u, double v,
                double w, double *s)
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4:
  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
  case TETRAHEDRON:
  case TETRAHEDRON_2:
  case TETRAHEDRON_3:
  case TETRAHEDRON_4:
    Message::Error("BF_Node_2F cannot be associated with this type of element");
    break;

  case QUADRANGLE:
  case QUADRANGLE_2:
  case QUADRANGLE_2_8N:
  case QUADRANGLE_3:
  case QUADRANGLE_4:
    *s = 0.0625 * (1. - u) * (1. - v) * (1. + u) * (1. + v);
    break;

  case HEXAHEDRON:
  case HEXAHEDRON_2:
  case HEXAHEDRON_2_20N:
  case HEXAHEDRON_3:
  case HEXAHEDRON_4:
    switch(NumEntity) {
    case 1:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. - w) * (1. + u) * (1. - v) *
           (1. + w);
      break;
    case 2:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. - w) * (1. + u) * (1. + v) *
           (1. - w);
      break;
    case 3:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. - w) * (1. - u) * (1. + v) *
           (1. + w);
      break;
    case 4:
      *s = 0.015625 * (1. + u) * (1. - v) * (1. - w) * (1. + u) * (1. + v) *
           (1. + w);
      break;
    case 5:
      *s = 0.015625 * (1. + u) * (1. + v) * (1. - w) * (1. - u) * (1. + v) *
           (1. + w);
      break;
    case 6:
      *s = 0.015625 * (1. - u) * (1. - v) * (1. + w) * (1. + u) * (1. + v) *
           (1. + w);
      break;
    default: WrongNumEntity;
    }
    break;

  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
    switch(NumEntity) {
    default:
      Message::Error("BF_Node_2F not ready for PRISM");
      // cannot do this, as dofs should only be associated with quad-faces; use
      // 15 and 18-node prisms instead, i.e. through BF_Node (which calls Gmsh).
      /*
    case 1  : *s = 0.25 * (1.-u-v) * (1.-w) *     u    * (1.+w) ; break ;
    case 3  : *s = 0.25 * (1.-u-v) * (1.-w) *       v  * (1.+w) ; break ;
    case 4  : *s = 0.25 *     u    * (1.-w) *       v  * (1.+w) ; break ;
    default : WrongNumEntity ;
      */
    }
    break;

  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    switch(NumEntity) {
    default: Message::Error("BF_Node_2F not ready for PYRAMID");
    }
    break;

  default: Message::Error("Unknown Element Type in BF_Node_2F"); break;
  }
}

#undef WrongNumEntity

/* -------- */
/*  Volume  */
/* -------- */

void BF_Node_2V(struct Element *Element, int NumEntity, double u, double v,
                double w, double *s)
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4:
  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
  case QUADRANGLE:
  case QUADRANGLE_2:
  case QUADRANGLE_2_8N:
  case QUADRANGLE_3:
  case QUADRANGLE_4:
  case TETRAHEDRON:
  case TETRAHEDRON_2:
  case TETRAHEDRON_3:
  case TETRAHEDRON_4:
  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    Message::Error("BF_Node_2V cannot be associated with this type of element");
    break;

  case HEXAHEDRON:
  case HEXAHEDRON_2:
  case HEXAHEDRON_2_20N:
  case HEXAHEDRON_3:
  case HEXAHEDRON_4:
    *s = 0.015625 * (1. - u) * (1. - v) * (1. - w) * (1. + u) * (1. + v) *
         (1. + w);
    break;

  default: Message::Error("Unknown type of Element in BF_Node_2V"); break;
  }
}

/* ------------------------------------------------------------------------ */
/*  B F _ G r a d N o d e _ 2                                               */
/* ------------------------------------------------------------------------ */

/* ------- */
/*  Edges  */
/* ------- */

#define WrongNumEntity Message::Error("Wrong Edge number in 'BF_GradNode_2E'")

void BF_GradNode_2E(struct Element *Element, int NumEntity, double u, double v,
                    double w, double s[])
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4:
    switch(NumEntity) {
    case 1:
      s[0] = -0.5 * u;
      s[1] = 0.;
      s[2] = 0.;
      break;
    default: WrongNumEntity;
    }
    break;

  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
    switch(NumEntity) {
    case 1:
      s[0] = 1. - 2. * u - v;
      s[1] = -u;
      s[2] = 0.;
      break;
    case 2:
      s[0] = -v;
      s[1] = 1. - u - 2. * v;
      s[2] = 0.;
      break;
    case 3:
      s[0] = v;
      s[1] = u;
      s[2] = 0.;
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
      s[0] = 0.0625 * (-2. * u) * (1. - v) * (1. - v);
      s[1] = 0.0625 * (1. - u * u) * (-2.) * (1. - v);
      s[2] = 0.;
      break;
    case 2:
      s[0] = 0.0625 * (-2.) * (1. - u) * (1. - v * v);
      s[1] = 0.0625 * (1. - u) * (1. - u) * (-2. * v);
      s[2] = 0.;
      break;
    case 3:
      s[0] = 0.0625 * (2.) * (1. + u) * (1. - v * v);
      s[1] = 0.0625 * (1. + u) * (1. + u) * (-2. * v);
      s[2] = 0.;
      break;
    case 4:
      s[0] = 0.0625 * (-2. * u) * (1. + v) * (1. + v);
      s[1] = 0.0625 * (1. - u * u) * (2.) * (1. + v);
      s[2] = 0.;
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
      s[0] = 1. - 2. * u - v - w;
      s[1] = -u;
      s[2] = -u;
      break;
    case 2:
      s[0] = -v;
      s[1] = 1. - u - 2. * v - w;
      s[2] = -v;
      break;
    case 3:
      s[0] = -w;
      s[1] = -w;
      s[2] = 1. - u - v - 2. * w;
      break;
    case 4:
      s[0] = v;
      s[1] = u;
      s[2] = 0.;
      break;
    case 5:
      s[0] = w;
      s[1] = 0.;
      s[2] = u;
      break;
    case 6:
      s[0] = 0.;
      s[1] = w;
      s[2] = v;
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
    case 1:
      s[0] = 0.015625 * (-2. * u) * (1. - v) * (1. - w) * (1. - v) * (1. - w);
      s[1] =
        0.015625 * (1. - u) * (-2.) * (1. - v) * (1. - w) * (1. + u) * (1. - w);
      s[2] =
        0.015625 * (1. - u) * (1. - v) * (-2.) * (1. - w) * (1. + u) * (1. - v);
      break;

    case 2:
      s[0] =
        0.015625 * (-2.) * (1. - u) * (1. - v) * (1. - w) * (1. + v) * (1. - w);
      s[1] = 0.015625 * (1. - u) * (-2. * v) * (1. - w) * (1. - u) * (1. - w);
      s[2] =
        0.015625 * (1. - u) * (1. - v) * (-2.) * (1. - w) * (1. - u) * (1. + v);
      break;

    case 3:
      s[0] =
        0.015625 * (-2.) * (1. - u) * (1. - v) * (1. - w) * (1. - v) * (1. + w);
      s[1] =
        0.015625 * (1. - u) * (-2.) * (1. - v) * (1. - w) * (1. - u) * (1. + w);
      s[2] = 0.015625 * (1. - u) * (1. - v) * (-2. * w) * (1. - u) * (1. - v);
      break;

    case 4:
      s[0] =
        0.015625 * 2. * (1. + u) * (1. - v) * (1. - w) * (1. + v) * (1. - w);
      s[1] = 0.015625 * (1. + u) * (-2. * v) * (1. - w) * (1. + u) * (1. - w);
      s[2] =
        0.015625 * (1. + u) * (1. - v) * (-2.) * (1. - w) * (1. + u) * (1. + v);
      break;

    case 5:
      s[0] =
        0.015625 * 2. * (1. + u) * (1. - v) * (1. - w) * (1. - v) * (1. + w);
      s[1] =
        0.015625 * (1. + u) * (-2.) * (1. - v) * (1. - w) * (1. + u) * (1. + w);
      s[2] = 0.015625 * (1. + u) * (1. - v) * (-2. * w) * (1. + u) * (1. - v);
      break;

    case 6:
      s[0] = 0.015625 * (-2. * u) * (1. + v) * (1. - w) * (1. + v) * (1. - w);
      s[1] =
        0.015625 * (1. + u) * 2. * (1. + v) * (1. - w) * (1. - u) * (1. - w);
      s[2] =
        0.015625 * (1. + u) * (1. + v) * (-2.) * (1. - w) * (1. - u) * (1. + v);
      break;

    case 7:
      s[0] =
        0.015625 * 2. * (1. + u) * (1. + v) * (1. - w) * (1. + v) * (1. + w);
      s[1] =
        0.015625 * (1. + u) * 2. * (1. + v) * (1. - w) * (1. + u) * (1. + w);
      s[2] = 0.015625 * (1. + u) * (1. + v) * (-2. * w) * (1. + u) * (1. + v);
      break;

    case 8:
      s[0] =
        0.015625 * (-2.) * (1. - u) * (1. + v) * (1. - w) * (1. + v) * (1. + w);
      s[1] =
        0.015625 * (1. - u) * 2. * (1. + v) * (1. - w) * (1. - u) * (1. + w);
      s[2] = 0.015625 * (1. - u) * (1. + v) * (-2. * w) * (1. - u) * (1. + v);
      break;

    case 9:
      s[0] = 0.015625 * (-2. * u) * (1. - v) * (1. + w) * (1. - v) * (1. + w);
      s[1] =
        0.015625 * (1. - u) * (-2.) * (1. - v) * (1. + w) * (1. + u) * (1. + w);
      s[2] =
        0.015625 * (1. - u) * (1. - v) * 2. * (1. + w) * (1. + u) * (1. - v);
      break;

    case 10:
      s[0] =
        0.015625 * (-2.) * (1. - u) * (1. - v) * (1. + w) * (1. + v) * (1. + w);
      s[1] = 0.015625 * (1. - u) * (-2. * v) * (1. + w) * (1. - u) * (1. + w);
      s[2] =
        0.015625 * (1. - u) * (1. - v) * 2. * (1. + w) * (1. - u) * (1. + v);
      break;

    case 11:
      s[0] =
        0.015625 * 2. * (1. + u) * (1. - v) * (1. + w) * (1. + v) * (1. + w);
      s[1] = 0.015625 * (1. + u) * (-2. * v) * (1. + w) * (1. + u) * (1. + w);
      s[2] =
        0.015625 * (1. + u) * (1. - v) * 2. * (1. + w) * (1. + u) * (1. + v);
      break;

    case 12:
      s[0] = 0.015625 * (-2. * u) * (1. + v) * (1. + w) * (1. + v) * (1. + w);
      s[1] =
        0.015625 * (1. + u) * 2. * (1. + v) * (1. + w) * (1. - u) * (1. + w);
      s[2] =
        0.015625 * (1. + u) * (1. + v) * 2. * (1. + w) * (1. - u) * (1. + v);
      break;

    default: WrongNumEntity;
    }
    break;

  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
    // FIXME: not tested!
    switch(NumEntity) {
    case 1:
      s[0] = 0.25 * (1. - 2. * u - v) * (1. - w) * (1. - w);
      s[1] = 0.25 * (-u) * (1. - w) * (1. - w);
      s[2] = 0.25 * (u - u * u - v * u) * (-2.) * (1. - w);
      break;

    case 2:
      s[0] = 0.25 * (-v) * (1. - w) * (1. - w);
      s[1] = 0.25 * (1 - u - 2 * v) * (1. - w) * (1. - w);
      s[2] = 0.25 * (v - u * v - v * v) * (-2.) * (1. - w);
      break;

    case 3:
      s[0] = 0.25 * (-2. + 2. * u + 2. * v) * (1. - w) * (1. + w);
      s[1] = 0.25 * (-2. + 2. * u + 2. * v) * (1. - w) * (1. + w);
      s[2] =
        0.25 * (1. - 2. * u - 2. * v + u * u + 2. * u * v + v * v) * (-2. * w);
      break;

    case 4:
      s[0] = 0.25 * (1. - w) * v * (1. - w);
      s[1] = 0.25 * u * (1. - w) * (1. - w);
      s[2] = 0.25 * u * (-2.) * (1. - w) * v;
      break;

    case 5:
      s[0] = 0.25 * 2. * u * (1. - w) * (1. + w);
      s[1] = 0.;
      s[2] = 0.25 * u * u * (-2. * w);
      break;

    case 6:
      s[0] = 0.;
      s[1] = 0.25 * 2. * v * (1. - w) * (1. + w);
      s[2] = 0.25 * v * v * (-2. * w);
      break;

    case 7:
      s[0] = 0.25 * (1. - 2. * u - v) * (1. + w) * (1. + w);
      s[1] = 0.25 * (-u) * (1. + w) * (1. + w);
      s[2] = 0.25 * (u - u * u - u * v) * 2. * (1. + w);
      break;

    case 8:
      s[0] = 0.25 * (-v) * (1. + w) * (1. + w);
      s[1] = 0.25 * (1. - u - 2. * v) * (1. + w) * (1. + w);
      s[2] = 0.25 * (v - u * v - v * v) * 2. * (1. + w);
      break;

    case 9:
      s[0] = 0.25 * (1. + w) * v * (1. + w);
      s[1] = 0.25 * u * (1. + w) * (1. + w);
      s[2] = 0.25 * u * 2. * (1. + w) * v;
      break;

    default: WrongNumEntity;
    }
    break;

  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    switch(NumEntity) {
    default: Message::Error("BF_GradNode_2E not ready for PYRAMID");
    }
    break;

  default: Message::Error("Unknown type of Element in BF_GradNode_2E"); break;
  }
}

#undef WrongNumEntity

/* -------- */
/*  Facets  */
/* -------- */

#define WrongNumEntity Message::Error("Wrong Face number in 'BF_GradNode_2F'")

void BF_GradNode_2F(struct Element *Element, int NumEntity, double u, double v,
                    double w, double s[])
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4:
  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
  case TETRAHEDRON:
  case TETRAHEDRON_2:
  case TETRAHEDRON_3:
  case TETRAHEDRON_4:
    Message::Error(
      "BF_GradNode_2F cannot be associated with this type of element");
    break;

  case QUADRANGLE:
  case QUADRANGLE_2:
  case QUADRANGLE_2_8N:
  case QUADRANGLE_3:
  case QUADRANGLE_4:
    s[0] = 0.0625 * (-2. * u) * (1. - v) * (1. + v);
    s[1] = 0.0625 * (1. - u) * (-2. * v) * (1. + u);
    s[2] = 0.;
    break;

  case HEXAHEDRON:
  case HEXAHEDRON_2:
  case HEXAHEDRON_2_20N:
  case HEXAHEDRON_3:
  case HEXAHEDRON_4:
    switch(NumEntity) {
    case 1:
      s[0] = 0.015625 * (-2. * u) * (1. - v) * (1. - w) * (1. - v) * (1. + w);
      s[1] =
        0.015625 * (1. - u) * (-2.) * (1. - v) * (1. - w) * (1. + u) * (1. + w);
      s[2] = 0.015625 * (1. - u) * (1. - v) * (-2. * w) * (1. + u) * (1. - v);
      break;

    case 2:
      s[0] = 0.015625 * (-2. * u) * (1. - v) * (1. - w) * (1. + v) * (1. - w);
      s[1] = 0.015625 * (1. - u) * (-2. * v) * (1. - w) * (1. + u) * (1. - w);
      s[2] =
        0.015625 * (1. - u) * (1. - v) * (-2.) * (1. - w) * (1. + u) * (1. + v);
      break;

    case 3:
      s[0] =
        0.015625 * (-2.) * (1. - u) * (1. - v) * (1. - w) * (1. + v) * (1. + w);
      s[1] = 0.015625 * (1. - u) * (-2. * v) * (1. - w) * (1. - u) * (1. + w);
      s[2] = 0.015625 * (1. - u) * (1. - v) * (-2. * w) * (1. - u) * (1. + v);
      break;

    case 4:
      s[0] =
        0.015625 * 2. * (1. + u) * (1. - v) * (1. - w) * (1. + v) * (1. + w);
      s[1] = 0.015625 * (1. + u) * (-2. * v) * (1. - w) * (1. + u) * (1. + w);
      s[2] = 0.015625 * (1. + u) * (1. - v) * (-2. * w) * (1. + u) * (1. + v);
      break;

    case 5:
      s[0] = 0.015625 * (-2. * u) * (1. + v) * (1. - w) * (1. + v) * (1. + w);
      s[1] =
        0.015625 * (1. + u) * 2. * (1. + v) * (1. - w) * (1. - u) * (1. + w);
      s[2] = 0.015625 * (1. + u) * (1. + v) * (-2. * w) * (1. - u) * (1. + v);
      break;

    case 6:
      s[0] = 0.015625 * (-2. * u) * (1. - v) * (1. + w) * (1. + v) * (1. + w);
      s[1] = 0.015625 * (1. - u) * (-2. * v) * (1. + w) * (1. + u) * (1. + w);
      s[2] =
        0.015625 * (1. - u) * (1. - v) * 2. * (1. + w) * (1. + u) * (1. + v);
      break;
    default: WrongNumEntity;
    }
    break;

  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
    switch(NumEntity) {
    default:
      Message::Error("BF_GradNode_2F not ready for PRISM");
      // cannot do this yet in getdp, as dofs should only be
      // associated with quad-faces: if really necessary we could
      // implement actual 15 and 18-node prisms
      /*
    case 1   : s[0] = 0.25 * (1.-2.*u-v) * (1.-w) * (1.+w) ;
               s[1] = 0.25 * (-u) * (1.-w) * (1.+w);
               s[2] = 0.25 * (u-u*u-u*v) * (-2.*w); break ;

    case 3   : s[0] = 0.25 * (-v) * (1.-w) * (1.+w) ;
               s[1] = 0.25 * (1.-u-2.*v) * (1.-w) * (1.+w) ;
               s[2] = 0.25 * (v-u*v-v*v) * (-2.*w) ; break;

    case 4   : s[0] = 0.25 * (1.-w) * v * (1.+w) ;
               s[1] = 0.25 * u * (1.-w) * (1.+w) ;
               s[2] = 0.25 * u * (-2.*w) * v ; break;
      */
    }
    break;

  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    switch(NumEntity) {
    default: Message::Error("BF_GradNode_2F not ready for PYRAMID");
    }
    break;

  default: Message::Error("Unknown type of Element in BF_GradNode_2F"); break;
  }
}

#undef WrongNumEntity

/* -------- */
/*  Volume  */
/* -------- */

void BF_GradNode_2V(struct Element *Element, int NumEntity, double u, double v,
                    double w, double s[])
{
  switch(Element->Type) {
  case LINE:
  case LINE_2:
  case LINE_3:
  case LINE_4:
  case TRIANGLE:
  case TRIANGLE_2:
  case TRIANGLE_3:
  case TRIANGLE_4:
  case QUADRANGLE:
  case QUADRANGLE_2:
  case QUADRANGLE_2_8N:
  case QUADRANGLE_3:
  case QUADRANGLE_4:
  case TETRAHEDRON:
  case TETRAHEDRON_2:
  case TETRAHEDRON_3:
  case TETRAHEDRON_4:
  case PRISM:
  case PRISM_2:
  case PRISM_2_15N:
  case PRISM_3:
  case PRISM_4:
  case PYRAMID:
  case PYRAMID_2:
  case PYRAMID_2_13N:
  case PYRAMID_3: // case PYRAMID_4
    Message::Error(
      "BF_GradNode_2V cannot be associated with this type of element");
    break;

  case HEXAHEDRON:
  case HEXAHEDRON_2:
  case HEXAHEDRON_2_20N:
  case HEXAHEDRON_3:
  case HEXAHEDRON_4:
    s[0] = 0.015625 * (-2. * u) * (1. - v) * (1. - w) * (1. + v) * (1. + w);
    s[1] = 0.015625 * (1. - u) * (-2. * v) * (1. - w) * (1. + u) * (1. + w);
    s[2] = 0.015625 * (1. - u) * (1. - v) * (-2. * w) * (1. + u) * (1. + v);
    break;

  default: Message::Error("Unknown type of Element in BF_GradNode_2V"); break;
  }
}
