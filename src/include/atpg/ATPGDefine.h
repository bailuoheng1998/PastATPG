// SPDX-FileCopyrightText: 2023 CASTest Corporation Limited
// SPDX-License-Identifier: LGPL-v3

#ifndef ICTEST_ATPGDEFINE_H
#define ICTEST_ATPGDEFINE_H

namespace ictest {
///////////////////////////////
// define supported logic value for D atpg
const int LOGIC_0 = 0;   // (good/bad) (0/0)
const int LOGIC_1 = 1;   // (good/bad) (1/1)
const int LOGIC_x = 2;   // (good/bad) (X/X)
const int LOGIC_S1 = 3;  // (good/bad) (1/0)
const int LOGIC_S0 = 4;  // (good/bad) (0/1)
const int LOGIC_G0 = 5;  // (good/bad) (0/X)
const int LOGIC_F0 = 6;  // (good/bad) (X/0)
const int LOGIC_G1 = 7;  // (good/bad) (1/X)
const int LOGIC_F1 = 8;  // (good/bad) (X/1)
const int LOGIC_CF = 9;  // conflict logic
const int LOGIC_UJ = 10; // unjustified logic

///////////////////////////////
// define supported logic value for podem atpg
// faulty/good
const int XX = 0xa;   // 10 (x)
const int ZERO = 0x0;   // 0  (0)
const int ONE =  0xf;   // 15 (1)
const int S1  =  0xc;   // 12 (DB)
const int S0  =  0x3;   // 3  (D)
const int S0X =  0x2;   // 2
const int S1X =  0xe;   // 14
const int X0  =  0x8;   // 8
const int X1  =  0xb;   // 11
const int ZZ  =  0x5;   // 5
const int Z0   = 0x4;   // 4
const int Z1  =  0x7;   // 7
const int ZX  =  0x6;   // 6
const int S1Z  = 0xd;   // 13
const int S0Z =  0x1;   // 1
const int XZ  =  0x9;   // 9

// Test generation status.
const int TEST_FOUND = 0;
const int OVER_BACKTRACK = 1;
const int NO_TEST = 2;

inline bool deter_val(int val) {
  if (val == LOGIC_0 || val == LOGIC_1) {
    return true;
  }

  return false;
}

inline bool unknown_val(int val) {
  if (val == LOGIC_x) {
    return true;
  }

  return false;
}

inline bool faulty_val(int val) {
  if (val == LOGIC_S0 || val == LOGIC_S1) {
    return true;
  }

  return false;
}

const int GOOD_VALUE[9] = {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_1, LOGIC_0,
                           LOGIC_0, LOGIC_x, LOGIC_1, LOGIC_x};
const int BAD_VALUE[9] = {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_0, LOGIC_1,
                          LOGIC_x, LOGIC_0, LOGIC_x, LOGIC_1};
// [dim1-good][dim2-bad]
const int NINE_VALUE[3][3] = {
    // 0        1         x   --dim2
    {LOGIC_0, LOGIC_S0, LOGIC_G0}, // 0 -- dim1
    {LOGIC_S1, LOGIC_1, LOGIC_G1}, // 1
    {LOGIC_F0, LOGIC_F1, LOGIC_x}  // x
};

// five-logic system
const int and_table_five[5][5] = {
    //  0        1        x        S1       S0
    {LOGIC_0, LOGIC_0, LOGIC_0, LOGIC_0, LOGIC_0},   // 0
    {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_S0}, // 1
    {LOGIC_0, LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x},   // x
    {LOGIC_0, LOGIC_S1, LOGIC_x, LOGIC_S1, LOGIC_0}, // S1
    {LOGIC_0, LOGIC_S0, LOGIC_x, LOGIC_0, LOGIC_S0}, // S0
};
const int or_table_five[5][5] = {
    //  0        1         x        S1        S0
    {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_S0}, // 0
    {LOGIC_1, LOGIC_1, LOGIC_1, LOGIC_1, LOGIC_1},   // 1
    {LOGIC_x, LOGIC_1, LOGIC_x, LOGIC_x, LOGIC_x},   // x
    {LOGIC_S1, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_1}, // S1
    {LOGIC_S0, LOGIC_1, LOGIC_x, LOGIC_1, LOGIC_S0}, // S0
};
const int xor_table_five[5][5] = {
    //  0        1        x        S1       S0
    {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_S0}, // 0
    {LOGIC_1, LOGIC_0, LOGIC_x, LOGIC_S0, LOGIC_S1}, // 1
    {LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x},   // x
    {LOGIC_S1, LOGIC_S0, LOGIC_x, LOGIC_0, LOGIC_1}, // S1
    {LOGIC_S0, LOGIC_S1, LOGIC_x, LOGIC_1, LOGIC_0}, // S0
};

const int not_table_five[5] = {LOGIC_1, LOGIC_0, LOGIC_x, LOGIC_S0, LOGIC_S1};
const int buf_table_five[5] = {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_S0};

// 0-LOGIC_0, 1-LOGIC_1, 2-LOGIC_x, 3-LOGIC_S1, 4-LOGIC_S0,
// 5-LOGIC_G0, 6-LOGIC_F0, 7-LOGIC_G1, 8-LOGIC_F1
const int and_table_nine[9][9] = {
    //       0        1        x        S1       S0       G0       F0 G1 F1
    {LOGIC_0, LOGIC_0, LOGIC_0, LOGIC_0, LOGIC_0, LOGIC_0, LOGIC_0, LOGIC_0,
     LOGIC_0}, // 0
    {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_S0, LOGIC_G0, LOGIC_F0,
     LOGIC_G1, LOGIC_F1}, // 1
    {LOGIC_0, LOGIC_x, LOGIC_x, LOGIC_F0, LOGIC_G0, LOGIC_G0, LOGIC_F0, LOGIC_x,
     LOGIC_x}, // x
    {LOGIC_0, LOGIC_S1, LOGIC_F0, LOGIC_S1, LOGIC_0, LOGIC_0, LOGIC_F0,
     LOGIC_S1, LOGIC_F0}, // S1
    {LOGIC_0, LOGIC_S0, LOGIC_G0, LOGIC_0, LOGIC_S0, LOGIC_G0, LOGIC_0,
     LOGIC_G0, LOGIC_S0}, // S0
    {LOGIC_0, LOGIC_G0, LOGIC_G0, LOGIC_0, LOGIC_G0, LOGIC_G0, LOGIC_0,
     LOGIC_G0, LOGIC_G0}, // G0
    {LOGIC_0, LOGIC_F0, LOGIC_F0, LOGIC_F0, LOGIC_0, LOGIC_0, LOGIC_F0,
     LOGIC_F0, LOGIC_F0}, // F0
    {LOGIC_0, LOGIC_G1, LOGIC_x, LOGIC_S1, LOGIC_G0, LOGIC_G0, LOGIC_F0,
     LOGIC_G1, LOGIC_x}, // G1
    {LOGIC_0, LOGIC_F1, LOGIC_x, LOGIC_F0, LOGIC_S0, LOGIC_G0, LOGIC_F0,
     LOGIC_x, LOGIC_F1} // F1
};
const int or_table_nine[9][9] = {
    //       0        1         x        S1        S0       G0       F0       G1
    //       F1
    {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_S0, LOGIC_G0, LOGIC_F0,
     LOGIC_G1, LOGIC_F1}, // 0
    {LOGIC_1, LOGIC_1, LOGIC_1, LOGIC_1, LOGIC_1, LOGIC_1, LOGIC_1, LOGIC_1,
     LOGIC_1}, // 1
    {LOGIC_x, LOGIC_1, LOGIC_x, LOGIC_G1, LOGIC_F1, LOGIC_x, LOGIC_x, LOGIC_G1,
     LOGIC_F1}, // x
    {LOGIC_S1, LOGIC_1, LOGIC_G1, LOGIC_S1, LOGIC_1, LOGIC_G1, LOGIC_S1,
     LOGIC_G1, LOGIC_1}, // S1
    {LOGIC_S0, LOGIC_1, LOGIC_F1, LOGIC_1, LOGIC_S0, LOGIC_S0, LOGIC_F1,
     LOGIC_1, LOGIC_F1}, // S0
    {LOGIC_G0, LOGIC_1, LOGIC_x, LOGIC_G1, LOGIC_S0, LOGIC_G0, LOGIC_x,
     LOGIC_G1, LOGIC_F1}, // G0
    {LOGIC_F0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_F1, LOGIC_x, LOGIC_F0,
     LOGIC_G1, LOGIC_F1}, // F0
    {LOGIC_G1, LOGIC_1, LOGIC_G1, LOGIC_G1, LOGIC_1, LOGIC_G1, LOGIC_G1,
     LOGIC_G1, LOGIC_1}, // G1
    {LOGIC_F1, LOGIC_1, LOGIC_F1, LOGIC_1, LOGIC_F1, LOGIC_F1, LOGIC_F1,
     LOGIC_1, LOGIC_F1} // F1
};
const int xor_table_nine[9][9] = {
    //       0        1        x        S1       S0       G0       F0 G1 F1
    {LOGIC_0, LOGIC_1, LOGIC_x, LOGIC_S1, LOGIC_S0, LOGIC_G0, LOGIC_F0,
     LOGIC_G1, LOGIC_F1}, // 0
    {LOGIC_1, LOGIC_0, LOGIC_x, LOGIC_S0, LOGIC_S1, LOGIC_G1, LOGIC_F1,
     LOGIC_G0, LOGIC_F0}, // 1
    {LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x, LOGIC_x,
     LOGIC_x}, // x
    {LOGIC_S1, LOGIC_S0, LOGIC_x, LOGIC_0, LOGIC_1, LOGIC_G1, LOGIC_F0,
     LOGIC_G0, LOGIC_F1}, // S1
    {LOGIC_S0, LOGIC_S1, LOGIC_x, LOGIC_1, LOGIC_0, LOGIC_G0, LOGIC_F1,
     LOGIC_G1, LOGIC_F0}, // S0
    {LOGIC_G0, LOGIC_G1, LOGIC_x, LOGIC_G1, LOGIC_G0, LOGIC_G0, LOGIC_x,
     LOGIC_G1, LOGIC_x}, // G0
    {LOGIC_F0, LOGIC_F1, LOGIC_x, LOGIC_F0, LOGIC_F1, LOGIC_x, LOGIC_F0,
     LOGIC_x, LOGIC_F1}, // F0
    {LOGIC_G1, LOGIC_G0, LOGIC_x, LOGIC_G0, LOGIC_G1, LOGIC_G1, LOGIC_x,
     LOGIC_G0, LOGIC_x}, // G1
    {LOGIC_F1, LOGIC_F0, LOGIC_x, LOGIC_F1, LOGIC_F0, LOGIC_x, LOGIC_F1,
     LOGIC_x, LOGIC_F0} // F1
};

// [SE][In0][In1]
const int three_value_mux_table[3][3][3] = {
    // 000      001      00x       010      011      01x       0x0      0x1 0xx
    {{LOGIC_0, LOGIC_0, LOGIC_0},
     {LOGIC_1, LOGIC_1, LOGIC_1},
     {LOGIC_x, LOGIC_x, LOGIC_x}},
    // 100      101      10x       110      111      11x       1x0      1x1 1xx
    {{LOGIC_0, LOGIC_1, LOGIC_x},
     {LOGIC_0, LOGIC_1, LOGIC_x},
     {LOGIC_0, LOGIC_1, LOGIC_x}},
    // x00      x01      x0x       x10      x11      x1x       xx0      xx1 xxx
    {{LOGIC_0, LOGIC_x, LOGIC_x},
     {LOGIC_x, LOGIC_1, LOGIC_x},
     {LOGIC_x, LOGIC_x, LOGIC_x}}};

const int not_table_nine[9] = {LOGIC_1,  LOGIC_0,  LOGIC_x,  LOGIC_S0, LOGIC_S1,
                               LOGIC_G1, LOGIC_F1, LOGIC_G0, LOGIC_F0};
const int buf_table_nine[9] = {LOGIC_0,  LOGIC_1,  LOGIC_x,  LOGIC_S1, LOGIC_S0,
                               LOGIC_G0, LOGIC_F0, LOGIC_G1, LOGIC_F1};

} // namespace ictest

#endif // ICTEST_ATPGDEFINE_H
