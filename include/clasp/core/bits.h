#ifndef core_bits_H
#define core_bits_H

namespace core {

// Line-end comments are the respective bit-array operator DEFUN's as named by the CL standard.
// We will reuse these boole values when implementing the bit-array operations, which is more convenient than correct.
typedef enum { boole_and = 0, // bit-and
               boole_andc1,   // bit-andc1
               boole_andc2,   // bit-andc2
               boole_nand,    // bit-nand
               boole_xor,     // bit-xor
               boole_ior,     // bit-ior
               boole_nor,     // bit-nor
               boole_orc1,    // bit-orc1
               boole_orc2,    // bit-orc2
               boole_eqv,     // bit-eqv
               boole_c1,      // bit-not
               boole_c2,
               boole_1,
               boole_2,
               boole_clr,
               boole_set,
               booleOperatorCount,
               bitArrayOperatorCount = boole_c2
} boole_op;

void initialize_bits();

};

#endif
