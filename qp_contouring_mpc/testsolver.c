/* Produced by CVXGEN, 2021-03-07 13:52:47 -0500.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: testsolver.c. */
/* Description: Basic test harness for solver.c. */
#include "solver.h"

Vars vars;
Params params;
Workspace work;
Settings settings;

double* perform_solve(int T, 
  double _Qk[100 * (T+2)], 
  double _Rk[3 * (T+1)], 
  double _fk[4 * (T+2)],
  double _Ak[40 * (T+1)],
  double _Bk[16 * (T+1)],
  double _gk[6 * (T+1)],
  double _Ck[2 * (T+2)],
  double _ug[1 * (T+2)],
  double _lg[1 * (T+2)],
  double _x[10]) {

  set_defaults();
  setup_indexing();
  load_default_data();
  settings.verbose = 0;

  params.x_init[0] = _x[0];
  params.x_init[1] = _x[1];
  params.x_init[2] = _x[2];
  params.x_init[3] = _x[3];
  params.x_init[4] = _x[4];
  params.x_init[5] = _x[5];
  params.x_init[6] = _x[6];
  params.x_init[7] = _x[7];
  params.x_init[8] = _x[8];
  params.x_init[9] = _x[9];

  int data_len;

  data_len = 100;
  for(int i = 0; i < data_len * (T+2); i++) {
    params.Qk[i/data_len][i%data_len] = _Qk[i];
  }

  data_len = 3;
  for(int i = 0; i < data_len * (T+1); i++) {
    params.Rk[i/data_len][i%data_len] = _Rk[i];
  }

  data_len = 4;
  for(int i = 0; i < data_len * (T+2); i++) {
    params.fk[i/data_len][i%data_len] = _fk[i];
  }

  data_len = 40;
  for(int i = 0; i < data_len * (T+1); i++) {
    params.Ak[i/data_len][i%data_len] = _Ak[i];
  }

  data_len = 16;
  for(int i = 0; i < data_len * (T+1); i++) {
    params.Bk[i/data_len][i%data_len] = _Bk[i];
  }

  data_len = 6;
  for(int i = 0; i < data_len * (T+1); i++) {
    params.gk[i/data_len][i%data_len] = _gk[i];
  }

  data_len = 2;
  for(int i = 0; i < data_len * (T+2); i++) {
    params.Ck[i/data_len][i%data_len] = _Ck[i];
  }

  data_len = 1;
  for(int i = 0; i < data_len * (T+2); i++) {
    params.ug[i/data_len][i%data_len] = _ug[i];
  }

  data_len = 1;
  for(int i = 0; i < data_len * (T+2); i++) {
    params.lg[i/data_len][i%data_len] = _lg[i];
  }

  solve();

  /*
  for(int i=0; i<10; i++) {
    printf("[%f %f %f %f %f %f %f %f %f %f]\n", vars.x[i][0], vars.x[i][1], vars.x[i][2], vars.x[i][3], vars.x[i][4], vars.x[i][5], vars.x[i][6], vars.x[i][7], vars.x[i][8], vars.x[i][9]);
  }
  printf("[%f %f %f]\n", vars.u[0][0], vars.u[0][1], vars.u[0][2]);
  */
  
  return vars.x;
}

void load_default_data(void) {
  /*
  params.Qk_0[0] = 1.5507979025745755;
  params.Qk_0[10] = 0;
  params.Qk_0[20] = 0;
  params.Qk_0[30] = 0;
  params.Qk_0[40] = 0;
  params.Qk_0[50] = 0;
  params.Qk_0[60] = 0;
  params.Qk_0[70] = 0;
  params.Qk_0[80] = 0;
  params.Qk_0[90] = 0;
  params.Qk_0[1] = 0;
  params.Qk_0[11] = 1.7081478226181048;
  params.Qk_0[21] = 0;
  params.Qk_0[31] = 0;
  params.Qk_0[41] = 0;
  params.Qk_0[51] = 0;
  params.Qk_0[61] = 0;
  params.Qk_0[71] = 0;
  params.Qk_0[81] = 0;
  params.Qk_0[91] = 0;
  params.Qk_0[2] = 0;
  params.Qk_0[12] = 0;
  params.Qk_0[22] = 1.2909047389129444;
  params.Qk_0[32] = 0;
  params.Qk_0[42] = 0;
  params.Qk_0[52] = 0;
  params.Qk_0[62] = 0;
  params.Qk_0[72] = 0;
  params.Qk_0[82] = 0;
  params.Qk_0[92] = 0;
  params.Qk_0[3] = 0;
  params.Qk_0[13] = 0;
  params.Qk_0[23] = 0;
  params.Qk_0[33] = 1.510827605197663;
  params.Qk_0[43] = 0;
  params.Qk_0[53] = 0;
  params.Qk_0[63] = 0;
  params.Qk_0[73] = 0;
  params.Qk_0[83] = 0;
  params.Qk_0[93] = 0;
  params.Qk_0[4] = 0;
  params.Qk_0[14] = 0;
  params.Qk_0[24] = 0;
  params.Qk_0[34] = 0;
  params.Qk_0[44] = 1.8929469543476547;
  params.Qk_0[54] = 0;
  params.Qk_0[64] = 0;
  params.Qk_0[74] = 0;
  params.Qk_0[84] = 0;
  params.Qk_0[94] = 0;
  params.Qk_0[5] = 0;
  params.Qk_0[15] = 0;
  params.Qk_0[25] = 0;
  params.Qk_0[35] = 0;
  params.Qk_0[45] = 0;
  params.Qk_0[55] = 1.896293088933438;
  params.Qk_0[65] = 0;
  params.Qk_0[75] = 0;
  params.Qk_0[85] = 0;
  params.Qk_0[95] = 0;
  params.Qk_0[6] = 0;
  params.Qk_0[16] = 0;
  params.Qk_0[26] = 0;
  params.Qk_0[36] = 0;
  params.Qk_0[46] = 0;
  params.Qk_0[56] = 0;
  params.Qk_0[66] = 1.1255853104638363;
  params.Qk_0[76] = 0;
  params.Qk_0[86] = 0;
  params.Qk_0[96] = 0;
  params.Qk_0[7] = 0;
  params.Qk_0[17] = 0;
  params.Qk_0[27] = 0;
  params.Qk_0[37] = 0;
  params.Qk_0[47] = 0;
  params.Qk_0[57] = 0;
  params.Qk_0[67] = 0;
  params.Qk_0[77] = 1.2072428781381868;
  params.Qk_0[87] = 0;
  params.Qk_0[97] = 0;
  params.Qk_0[8] = 0;
  params.Qk_0[18] = 0;
  params.Qk_0[28] = 0;
  params.Qk_0[38] = 0;
  params.Qk_0[48] = 0;
  params.Qk_0[58] = 0;
  params.Qk_0[68] = 0;
  params.Qk_0[78] = 0;
  params.Qk_0[88] = 1.0514672033008299;
  params.Qk_0[98] = 0;
  params.Qk_0[9] = 0;
  params.Qk_0[19] = 0;
  params.Qk_0[29] = 0;
  params.Qk_0[39] = 0;
  params.Qk_0[49] = 0;
  params.Qk_0[59] = 0;
  params.Qk_0[69] = 0;
  params.Qk_0[79] = 0;
  params.Qk_0[89] = 0;
  params.Qk_0[99] = 1.4408098436506365;
  params.Rk_0[0] = 1.0298762108785668;
  params.Rk_0[1] = 1.456833224394711;
  params.Rk_0[2] = 1.6491440476147607;
  params.fk_0[0] = -0.8860508694080989;
  params.fk_0[1] = 0.7050196079205251;
  params.fk_0[2] = 0.3634512696654033;
  params.fk_0[3] = -1.9040724704913385;
  params.Qk_1[0] = 1.5588540879908819;
  params.Qk_1[10] = 0;
  params.Qk_1[20] = 0;
  params.Qk_1[30] = 0;
  params.Qk_1[40] = 0;
  params.Qk_1[50] = 0;
  params.Qk_1[60] = 0;
  params.Qk_1[70] = 0;
  params.Qk_1[80] = 0;
  params.Qk_1[90] = 0;
  params.Qk_1[1] = 0;
  params.Qk_1[11] = 1.2592524469074653;
  params.Qk_1[21] = 0;
  params.Qk_1[31] = 0;
  params.Qk_1[41] = 0;
  params.Qk_1[51] = 0;
  params.Qk_1[61] = 0;
  params.Qk_1[71] = 0;
  params.Qk_1[81] = 0;
  params.Qk_1[91] = 0;
  params.Qk_1[2] = 0;
  params.Qk_1[12] = 0;
  params.Qk_1[22] = 1.4151011970100695;
  params.Qk_1[32] = 0;
  params.Qk_1[42] = 0;
  params.Qk_1[52] = 0;
  params.Qk_1[62] = 0;
  params.Qk_1[72] = 0;
  params.Qk_1[82] = 0;
  params.Qk_1[92] = 0;
  params.Qk_1[3] = 0;
  params.Qk_1[13] = 0;
  params.Qk_1[23] = 0;
  params.Qk_1[33] = 1.2835250817713186;
  params.Qk_1[43] = 0;
  params.Qk_1[53] = 0;
  params.Qk_1[63] = 0;
  params.Qk_1[73] = 0;
  params.Qk_1[83] = 0;
  params.Qk_1[93] = 0;
  params.Qk_1[4] = 0;
  params.Qk_1[14] = 0;
  params.Qk_1[24] = 0;
  params.Qk_1[34] = 0;
  params.Qk_1[44] = 1.6931379183129964;
  params.Qk_1[54] = 0;
  params.Qk_1[64] = 0;
  params.Qk_1[74] = 0;
  params.Qk_1[84] = 0;
  params.Qk_1[94] = 0;
  params.Qk_1[5] = 0;
  params.Qk_1[15] = 0;
  params.Qk_1[25] = 0;
  params.Qk_1[35] = 0;
  params.Qk_1[45] = 0;
  params.Qk_1[55] = 1.4404537176707395;
  params.Qk_1[65] = 0;
  params.Qk_1[75] = 0;
  params.Qk_1[85] = 0;
  params.Qk_1[95] = 0;
  params.Qk_1[6] = 0;
  params.Qk_1[16] = 0;
  params.Qk_1[26] = 0;
  params.Qk_1[36] = 0;
  params.Qk_1[46] = 0;
  params.Qk_1[56] = 0;
  params.Qk_1[66] = 1.1568677384749633;
  params.Qk_1[76] = 0;
  params.Qk_1[86] = 0;
  params.Qk_1[96] = 0;
  params.Qk_1[7] = 0;
  params.Qk_1[17] = 0;
  params.Qk_1[27] = 0;
  params.Qk_1[37] = 0;
  params.Qk_1[47] = 0;
  params.Qk_1[57] = 0;
  params.Qk_1[67] = 0;
  params.Qk_1[77] = 1.5446490180318446;
  params.Qk_1[87] = 0;
  params.Qk_1[97] = 0;
  params.Qk_1[8] = 0;
  params.Qk_1[18] = 0;
  params.Qk_1[28] = 0;
  params.Qk_1[38] = 0;
  params.Qk_1[48] = 0;
  params.Qk_1[58] = 0;
  params.Qk_1[68] = 0;
  params.Qk_1[78] = 0;
  params.Qk_1[88] = 1.780314764511367;
  params.Qk_1[98] = 0;
  params.Qk_1[9] = 0;
  params.Qk_1[19] = 0;
  params.Qk_1[29] = 0;
  params.Qk_1[39] = 0;
  params.Qk_1[49] = 0;
  params.Qk_1[59] = 0;
  params.Qk_1[69] = 0;
  params.Qk_1[79] = 0;
  params.Qk_1[89] = 0;
  params.Qk_1[99] = 1.3063635323761797;
  params.fk_1[0] = -1.1121684642712744;
  params.fk_1[1] = -0.44811496977740495;
  params.fk_1[2] = 1.7455345994417217;
  params.fk_1[3] = 1.9039816898917352;
  params.Ak_0[0] = 0.6895347036512547;
  params.Ak_0[1] = 1.6113364341535923;
  params.Ak_0[2] = 1.383003485172717;
  params.Ak_0[3] = -0.48802383468444344;
  params.Ak_0[4] = -1.631131964513103;
  params.Ak_0[5] = 0.6136436100941447;
  params.Ak_0[6] = 0.2313630495538037;
  params.Ak_0[7] = -0.5537409477496875;
  params.Ak_0[8] = -1.0997819806406723;
  params.Ak_0[9] = -0.3739203344950055;
  params.Ak_0[10] = -0.12423900520332376;
  params.Ak_0[11] = -0.923057686995755;
  params.Ak_0[12] = -0.8328289030982696;
  params.Ak_0[13] = -0.16925440270808823;
  params.Ak_0[14] = 1.442135651787706;
  params.Ak_0[15] = 0.34501161787128565;
  params.Ak_0[16] = -0.8660485502711608;
  params.Ak_0[17] = -0.8880899735055947;
  params.Ak_0[18] = -0.1815116979122129;
  params.Ak_0[19] = -1.17835862158005;
  params.Ak_0[20] = -1.1944851558277074;
  params.Ak_0[21] = 0.05614023926976763;
  params.Ak_0[22] = -1.6510825248767813;
  params.Ak_0[23] = -0.06565787059365391;
  params.Ak_0[24] = -0.5512951504486665;
  params.Ak_0[25] = 0.8307464872626844;
  params.Ak_0[26] = 0.9869848924080182;
  params.Ak_0[27] = 0.7643716874230573;
  params.Ak_0[28] = 0.7567216550196565;
  params.Ak_0[29] = -0.5055995034042868;
  params.Ak_0[30] = 0.6725392189410702;
  params.Ak_0[31] = -0.6406053441727284;
  params.Ak_0[32] = 0.29117547947550015;
  params.Ak_0[33] = -0.6967713677405021;
  params.Ak_0[34] = -0.21941980294587182;
  params.Ak_0[35] = -1.753884276680243;
  params.Ak_0[36] = -1.0292983112626475;
  params.Ak_0[37] = 1.8864104246942706;
  params.Ak_0[38] = -1.077663182579704;
  params.Ak_0[39] = 0.7659100437893209;
  params.Bk_0[0] = 0.6019074328549583;
  params.Bk_0[1] = 0.8957565577499285;
  params.Bk_0[2] = -0.09964555746227477;
  params.Bk_0[3] = 0.38665509840745127;
  params.Bk_0[4] = -1.7321223042686946;
  params.Bk_0[5] = -1.7097514487110663;
  params.Bk_0[6] = -1.2040958948116867;
  params.Bk_0[7] = -1.3925560119658358;
  params.Bk_0[8] = -1.5995826216742213;
  params.Bk_0[9] = -1.4828245415645833;
  params.Bk_0[10] = 0.21311092723061398;
  params.Bk_0[11] = -1.248740700304487;
  params.Bk_0[12] = 1.808404972124833;
  params.Bk_0[13] = 0.7264471152297065;
  params.Bk_0[14] = 0.16407869343908477;
  params.Bk_0[15] = 0.8287224032315907;
  params.gk_0[0] = -0.9444533161899464;
  params.gk_0[1] = 1.7069027370149112;
  params.gk_0[2] = 1.3567722311998827;
  params.gk_0[3] = 0.9052779937121489;
  params.gk_0[4] = -0.07904017565835986;
  params.gk_0[5] = 1.3684127435065871;
  params.Ck_0[0] = 0.979009293697437;
  params.Ck_0[1] = 0.6413036255984501;
  params.ug_0[0] = 1.6559010680237511;
  params.ug_0[1] = 0.5346622551502991;
  params.ug_0[2] = -0.5362376605895625;
  params.ug_0[3] = 0.2113782926017822;
  params.ug_0[4] = -1.2144776931994525;
  params.ug_0[5] = -1.2317108144255875;
  params.ug_0[6] = 0.9026784957312834;
  params.ug_0[7] = 1.1397468137245244;
  params.ug_0[8] = 1.8883934547350631;
  params.ug_0[9] = 1.4038856681660068;
  params.Ck_1[0] = 0.17437730638329096;
  params.Ck_1[1] = -1.6408365219077408;
  params.ug_1[0] = -0.04450702153554875;
  params.ug_1[1] = 1.7117453902485025;
  params.ug_1[2] = 1.1504727980139053;
  params.ug_1[3] = -0.05962309578364744;
  params.ug_1[4] = -0.1788825540764547;
  params.ug_1[5] = -1.1280569263625857;
  params.ug_1[6] = -1.2911464767927057;
  params.ug_1[7] = -1.7055053231225696;
  params.ug_1[8] = 1.56957275034837;
  params.ug_1[9] = 0.5607064675962357;
  params.lg_0[0] = -1.4266707301147146;
  params.lg_0[1] = -0.3434923211351708;
  params.lg_0[2] = -1.8035643024085055;
  params.lg_0[3] = -1.1625066019105454;
  params.lg_0[4] = 0.9228324965161532;
  params.lg_0[5] = 0.6044910817663975;
  params.lg_0[6] = -0.0840868104920891;
  params.lg_0[7] = -0.900877978017443;
  params.lg_0[8] = 0.608892500264739;
  params.lg_0[9] = 1.8257980452695217;
  params.lg_1[0] = -0.25791777529922877;
  params.lg_1[1] = -1.7194699796493191;
  params.lg_1[2] = -1.7690740487081298;
  params.lg_1[3] = -1.6685159248097703;
  params.lg_1[4] = 1.8388287490128845;
  params.lg_1[5] = 0.16304334474597537;
  params.lg_1[6] = 1.3498497306788897;
  params.lg_1[7] = -1.3198658230514613;
  params.lg_1[8] = -0.9586197090843394;
  params.lg_1[9] = 0.7679100474913709;
  */

  params.lx_bounds[0] = -1;
  params.lx_bounds[1] = -1;
  params.lx_bounds[2] = -3;
  params.lx_bounds[3] = 0;
  params.lx_bounds[4] = -1;
  params.lx_bounds[5] = -1;
  params.lx_bounds[6] = 0;
  params.lx_bounds[7] = -0.1;
  params.lx_bounds[8] = -1;
  params.lx_bounds[9] = 0;
  params.lu_bounds[0] = -3;
  params.lu_bounds[1] = -3;
  params.lu_bounds[2] = -5;

  params.ux_bounds[0] = 1;
  params.ux_bounds[1] = 1;
  params.ux_bounds[2] = 3;
  params.ux_bounds[3] = 10;
  params.ux_bounds[4] = 1;
  params.ux_bounds[5] = 1;
  params.ux_bounds[6] = 1;
  params.ux_bounds[7] = 1;
  params.ux_bounds[8] = 1;
  params.ux_bounds[9] = 1;
  params.uu_bounds[0] = 3;
  params.uu_bounds[1] = 3;
  params.uu_bounds[2] = 5;
}
