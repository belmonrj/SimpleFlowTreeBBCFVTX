#ifndef RPPAR_H_
#define RPPAR_H_
static const int NDET = 8; // number of detectors (BBCs Psi2, FVTX Psi2, BBCs corrected Psi2, FVTX corrected Psi2 for all, then four layers)
static const int NHAR = 3; // number of harmonics (1, 2, 3)
static const int NMUL = 1; // Multiplicity bins (centrality selections)
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
