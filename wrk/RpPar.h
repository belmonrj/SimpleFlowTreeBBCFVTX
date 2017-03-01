#ifndef RPPAR_H_
#define RPPAR_H_
static const int NDET = 34; // number of detectors (see below)
//  0 -- BBCS  Psi2 from node tree
//  1 -- FVTXS Psi2 from node tree
//  2 -- BBCS Psi2 calculated
//  3 -- FVTXS Psi2 calculated (with weight)
//  4 -- FVTXS layer 0 Psi2 calculated (with weight)
//  5 -- FVTXS layer 1 Psi2 calculated (with weight)
//  6 -- FVTXS layer 2 Psi2 calculated (with weight)
//  7 -- FVTXS layer 3 Psi2 calculated (with weight)
//  8 -- FVTXN Psi2 calculated (with weight)
//  9 -- FVTXN layer 0 Psi2 calculated (with weight)
// 10 -- FVTXN layer 1 Psi2 calculated (with weight)
// 11 -- FVTXN layer 2 Psi2 calculated (with weight)
// 12 -- FVTXN layer 3 Psi2 calculated (with weight)
// 13 -- FVTXS Psi2 calculated (without weight)
// 14 -- FVTXS layer 0 Psi2 calculated (without weight)
// 15 -- FVTXS layer 1 Psi2 calculated (without weight)
// 16 -- FVTXS layer 2 Psi2 calculated (without weight)
// 17 -- FVTXS layer 3 Psi2 calculated (without weight)
// 18 -- FVTXN Psi2 calculated (without weight)
// 19 -- FVTXN layer 0 Psi2 calculated (without weight)
// 20 -- FVTXN layer 1 Psi2 calculated (without weight)
// 21 -- FVTXN layer 2 Psi2 calculated (without weight)
// 22 -- FVTXN layer 3 Psi2 calculated (without weight)
// 23 -- BBCS Psi2 calculated (without gain correction)
// 24 -- FVTXS tracks all
// 25 -- FVTXS tracks inner
// 26 -- FVTXS tracks outer
// 27 -- FVTXN tracks all
// 28 -- FVTXN tracks inner
// 29 -- FVTXN tracks outer
// 30 -- FVTXS clusters 0 1 2
// 31 -- FVTXS clusters 0 1 3
// 32 -- FVTXS clusters 0 2 3
// 33 -- FVTXS clusters 1 2 3
static const int NDETSHORT = 3; // number of detectors (see below)
//  0 -- BBCS Psi2 calculated
//  1 -- FVTXS Psi2 calculated
//  2 -- FVTXN Psi2 calculated
static const int NHAR = 4; // number of harmonics (1, 2, 3)
static const int NMUL = 6; // centrality bins 0-5, 5-10, 10-20, 20-40, 40-60, 60-88
static const int NZPS = 10; // number of z-vertex bins for flattening (very important for FVTX, maybe not so important for BBC)
static const int NORD = 12; // number of orders for fourier fit of Psi2 distribution (sin and cos?)
// --- all below here are irrelevant for now?
// --- all for central arm track stuff
static const int NZED = 10; // matching ???
static const int NMAT = 8; // all of these are for cnt matching ???
static const int NPAR = 20; // for cnt matching ???
static const int NETA = 8; // ???
static const int NTRK = 50; // ???

const float phbeta = 29.9792458;
const float mpion  = 0.139570;
const float mkaon  = 0.493677;
const float mproton= 0.938270;
const float pi = 3.1415927;

#endif /*RPPAR_H_*/
