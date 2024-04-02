// GetDP - Copyright (C) 1997-2022 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.

/* 8 integration points from Coulomb et al., IEEE tr.mag. 32(3) May
   1996, p.1395

   2 plans // a la base quadrangulaire, 4 points par plan suffisant
   pour integrer exactement nodal degre 2

   cf. ../utils/pyram.c */

static double upyr8[8] = {0.2631840555694285, -0.2631840555694285,
                          0.2631840555694285, -0.2631840555694285,
                          0.5066163033492386, -0.5066163033492386,
                          0.5066163033492386, -0.5066163033492386};
static double vpyr8[8] = {0.2631840555694285,  0.2631840555694285,
                          -0.2631840555694285, -0.2631840555694285,
                          0.5066163033492386,  0.5066163033492386,
                          -0.5066163033492386, -0.5066163033492386};
static double wpyr8[8] = {
  0.544151844011225, 0.544151844011225, 0.544151844011225, 0.544151844011225,
  0.122514822655441, 0.122514822655441, 0.122514822655441, 0.122514822655441};
static double ppyr8[8] = {
  0.100785882079825, 0.100785882079825, 0.100785882079825, 0.100785882079825,
  0.232547451253508, 0.232547451253508, 0.232547451253508, 0.232547451253508};
