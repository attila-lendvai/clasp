#include <clasp/core/foundation.h>
#include <clasp/core/corePackage.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/bformat.h>
#include <clasp/core/array.h>
#include <clasp/core/bits.h>
#include <clasp/core/wrappers.h>

/*
    num_log.c  -- Logical operations on numbers.
*/
/*
    Copyright (c) 1984, Taiichi Yuasa and Masami Hagiya.
    Copyright (c) 1990, Giuseppe Attardi.
    Copyright (c) 2001, Juan Jose Garcia Ripoll.

    ECL is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    See file '../Copyright' for full details.
*/

namespace core {

/*
 * BIT OPERATIONS FOR FIXNUMS
 */

#define FIXNUM_BITWISE_OP(name) static gctools::Fixnum fixnum_bitwise_##name(gctools::Fixnum i, gctools::Fixnum j)

FIXNUM_BITWISE_OP(and) {
  return i & j;
}

FIXNUM_BITWISE_OP(andc1) {
  return (~i) & j;
}

FIXNUM_BITWISE_OP(andc2) {
  return i & (~j);
}

FIXNUM_BITWISE_OP(nand) {
  return ~(i & j);
}

FIXNUM_BITWISE_OP(xor) {
  return i ^ j;
}

FIXNUM_BITWISE_OP(ior) {
  return i | j;
}

FIXNUM_BITWISE_OP(nor) {
  return ~(i | j);
}

FIXNUM_BITWISE_OP(orc1) {
  return (~i) | j;
}

FIXNUM_BITWISE_OP(orc2) {
  return i | (~j);
}

FIXNUM_BITWISE_OP(eqv) {
  return ~(i ^ j);
}

FIXNUM_BITWISE_OP(c1) {
  return ~i;
}

FIXNUM_BITWISE_OP(c2) {
  return ~j;
}

FIXNUM_BITWISE_OP(1) {
  return i;
}

FIXNUM_BITWISE_OP(2) {
  return j;
}

FIXNUM_BITWISE_OP(clr) {
  return 0;
}

FIXNUM_BITWISE_OP(set) {
  return -1;
}

typedef gctools::Fixnum (*fixnum_bitwise_operator_t)(gctools::Fixnum, gctools::Fixnum);

// NOTE: the order here must match the boole_op enum in bits.h
static fixnum_bitwise_operator_t fixnum_bitwise_operators[booleOperatorCount] = {
    fixnum_bitwise_and,
    fixnum_bitwise_andc1,
    fixnum_bitwise_andc2,
    fixnum_bitwise_nand,
    fixnum_bitwise_xor,
    fixnum_bitwise_ior,
    fixnum_bitwise_nor,
    fixnum_bitwise_orc1,
    fixnum_bitwise_orc2,
    fixnum_bitwise_eqv,
    fixnum_bitwise_c1,
    fixnum_bitwise_c2,
    fixnum_bitwise_1,
    fixnum_bitwise_2,
    fixnum_bitwise_clr,
    fixnum_bitwise_set
};

// ----------------------------------------------------------------------

#define BIGNUM_BITWISE_OP(name) static void bignum_bitwise_##name(Bignum_sp out, Bignum_sp i, Bignum_sp j)

BIGNUM_BITWISE_OP(orc1); // need to declare this signature to be able to retain the order of the enum below

BIGNUM_BITWISE_OP(and) {
  mpz_and(out->get().get_mpz_t(), i->get().get_mpz_t(), j->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(andc1) {
  mpz_com(out->get().get_mpz_t(), i->get().get_mpz_t());
  mpz_and(out->get().get_mpz_t(), out->get().get_mpz_t(), j->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(andc2) {
  /* (i & ~j) = ~((~i) | j) */
  bignum_bitwise_orc1(out, i, j);
  mpz_com(out->get().get_mpz_t(), out->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(nand) {
  mpz_and(out->get().get_mpz_t(), i->get().get_mpz_t(), j->get().get_mpz_t());
  mpz_com(out->get().get_mpz_t(), out->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(xor) {
  mpz_xor(out->get().get_mpz_t(), i->get().get_mpz_t(), j->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(ior) {
  mpz_ior(out->get().get_mpz_t(), i->get().get_mpz_t(), j->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(nor) {
  mpz_ior(out->get().get_mpz_t(), i->get().get_mpz_t(), j->get().get_mpz_t());
  mpz_com(out->get().get_mpz_t(), out->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(orc1) {
  mpz_com(out->get().get_mpz_t(), i->get().get_mpz_t());
  mpz_ior(out->get().get_mpz_t(), out->get().get_mpz_t(), j->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(orc2) {
  // (i | ~j) = ~((~i) & j)
  bignum_bitwise_andc1(out, i, j);
  mpz_com(out->get().get_mpz_t(), out->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(eqv) {
  mpz_xor(out->get().get_mpz_t(), i->get().get_mpz_t(), j->get().get_mpz_t());
  mpz_com(out->get().get_mpz_t(), out->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(c1) {
  mpz_com(out->get().get_mpz_t(), i->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(c2) {
  mpz_com(out->get().get_mpz_t(), j->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(1) {
  if (i != out)
    mpz_set(out->get().get_mpz_t(), i->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(2) {
  if (j != out)
    mpz_set(out->get().get_mpz_t(), j->get().get_mpz_t());
}

BIGNUM_BITWISE_OP(clr) {
  mpz_set_si(out->get().get_mpz_t(), 0);
}

BIGNUM_BITWISE_OP(set) {
  mpz_set_si(out->get().get_mpz_t(), -1);
}

typedef void (*_bignum_bitwise_operator_t)(Bignum_sp out, Bignum_sp o1, Bignum_sp o2);

// NOTE: the order here must match the boole_op enum in bits.h
static _bignum_bitwise_operator_t bignum_bitwise_operators[booleOperatorCount] = {
    bignum_bitwise_and,
    bignum_bitwise_andc1,
    bignum_bitwise_andc2,
    bignum_bitwise_nand,
    bignum_bitwise_xor,
    bignum_bitwise_ior,
    bignum_bitwise_nor,
    bignum_bitwise_orc1,
    bignum_bitwise_orc2,
    bignum_bitwise_eqv,
    bignum_bitwise_c1,
    bignum_bitwise_c2,
    bignum_bitwise_1,
    bignum_bitwise_2,
    bignum_bitwise_clr,
    bignum_bitwise_set
};

// ----------------------------------------------------------------------

T_sp clasp_boole(int op, T_sp x, T_sp y) {
  if (x.nilp() || y.nilp()) {
    SIMPLE_ERROR(BF("boole cannot accept nil"));
  }
  if (op < 0 || op >= booleOperatorCount)
     ERROR_WRONG_TYPE_NTH_ARG(cl::_sym_boole, 1, x, Cons_O::createList(cl::_sym_Integer_O, make_fixnum(0), make_fixnum(booleOperatorCount - 1)));
  if (x.fixnump()) {
    Fixnum_sp fnx = gc::As<Fixnum_sp>(x);
    if (y.fixnump()) { //Fixnum_sp fny = y.asOrNull<Fixnum_O>() ) {
      Fixnum_sp fny = gc::As<Fixnum_sp>(y);
      gctools::Fixnum z = fixnum_bitwise_operators[op](unbox_fixnum(fnx), unbox_fixnum(fny));
      return make_fixnum(z);
    } else if (Bignum_sp bny = y.asOrNull<Bignum_O>()) {
      Bignum_sp x_copy = my_thread->bigRegister0();
      x_copy->setFixnum(unbox_fixnum(fnx));
      bignum_bitwise_operators[op](x_copy, x_copy, bny);
      return _clasp_big_register_normalize(x_copy);
    } else {
      ERROR_WRONG_TYPE_NTH_ARG(cl::_sym_boole, 2, y, cl::_sym_integer);
    }
  } else if (Bignum_sp bnx = x.asOrNull<Bignum_O>()) {
    Bignum_sp x_copy = my_thread->bigRegister0();
    if (y.fixnump()) { // Fixnum_sp fny = y.asOrNull<Fixnum_O>() ) {
      Fixnum_sp fny(gc::As<Fixnum_sp>(y));
      Bignum_sp bny = my_thread->bigRegister1();
      bny->setFixnum(unbox_fixnum(fny));
      bignum_bitwise_operators[op](x_copy, bnx, bny);
      clasp_big_register_free(bny);
    } else if (Bignum_sp bny = y.asOrNull<Bignum_O>()) {
      bignum_bitwise_operators[op](x_copy, x, bny);
    } else {
      ERROR_WRONG_TYPE_NTH_ARG(cl::_sym_boole, 2, y, cl::_sym_integer);
    }
    return _clasp_big_register_normalize(x_copy);
  } else {
    ERROR_WRONG_TYPE_NTH_ARG(cl::_sym_boole, 1, x, cl::_sym_integer);
  }
  return x;
}

#if BIT_ARRAY_BYTE_SIZE==8

CL_LAMBDA(op x y &optional r);
CL_DECLARE();
CL_DOCSTRING("bitArrayOp");
CL_DEFUN T_sp core__bit_array_op(int opid, T_sp tx, T_sp ty, T_sp tr) {
  gctools::Fixnum i, j, n, d;
  SimpleBitVector_sp r0;
  size_t startr0 = 0;
  fixnum_bitwise_operator_t op;
  bool replace = false;
  int xi, yi, ri;
  byte8_t *xp, *yp, *rp;
  int xo, yo, ro;
  AbstractSimpleVector_sp ax;
  size_t startx, endx;
  AbstractSimpleVector_sp ay;
  size_t starty, endy;
  Array_sp array_x = gc::As<Array_sp>(tx);
  array_x->asAbstractSimpleVectorRange(ax, startx, endx);
  SimpleBitVector_sp x = gc::As_unsafe<SimpleBitVector_sp>(ax);
  Array_sp array_y = gc::As<Array_sp>(ty);
  array_y->asAbstractSimpleVectorRange(ay, starty, endy);
  SimpleBitVector_sp y = gc::As_unsafe<SimpleBitVector_sp>(ay);
  SimpleBitVector_sp r;
  size_t startr, endr;

  if (opid < 0 || opid >= bitArrayOperatorCount)
    goto ERROR;

  d = (endx - startx); // x->arrayTotalSize();
  xp = x->bytes();
  xo = startx; // x->offset();
  if (d != array_y->arrayTotalSize())
    goto ERROR;
  yp = y->bytes();
  yo = starty; // y->offset();
  if (tr == _lisp->_true())
    tr = x;
  if (tr.notnilp()) {
    AbstractSimpleVector_sp ar;
    Array_sp array_r = gc::As<Array_sp>(tr);
    array_r->asAbstractSimpleVectorRange(ar, startr, endr);
    r = gc::As_unsafe<SimpleBitVector_sp>(ar);
    if (!r) {
      ERROR_WRONG_TYPE_NTH_ARG(core::_sym_bitArrayOp, 4, tr, cl::_sym_SimpleBitVector_O);
    }
    if (endr-startr != d) //(r->arrayTotalSize() != d)
      goto ERROR;
//    i = (r->bytes() - xp) * 8 + (r->offset() - xo);
    i = (r->bytes()-xp)*8+startr-xo;
    if ((i > 0 && i < d) || (i < 0 && -i < d)) {
      r0 = r;
      startr0 = startr;
      tr = _Nil<T_O>();
      replace = true;
      goto L1;
    }
    // i = (r->bytes() - yp) * 8 + (r->offset() - yo);
    i = (r->bytes() - yp) * 8 + (startr - yo);
    if ((i > 0 && i < d) || (i < 0 && -i < d)) {
      r0 = r;
      startr0 = startr;
      tr = _Nil<T_O>();
      replace = true;
    }
  }
L1:
  if (tr.nilp()) {
    startr = 0;
    endr = d;
    r = SimpleBitVector_O::make(d);
  }
  rp = r->bytes();
  ro = startr; // r->offset();
  op = fixnum_bitwise_operators[opid];

#define set_high8(place, nbits, value) \
  (place) = ((place) & ~(-0400 >> (nbits))) | ((value) & (-0400 >> (nbits)));

#define set_low8(place, nbits, value) \
  (place) = ((place) & (-0400 >> (8 - (nbits)))) | ((value) & ~(-0400 >> (8 - (nbits))));

#define extract_byte8(integer, pointer, index, offset) \
  (integer) = (pointer)[(index)+1] & 0377;            \
  (integer) = ((pointer)[index] << (offset)) | ((integer) >> (8 - (offset)));

#define store_byte8(pointer, index, offset, value)               \
  set_low8((pointer)[index], 8 - (offset), (value) >> (offset)); \
  set_high8((pointer)[(index)+1], offset, (value) << (8 - (offset)));

  //
  if (xo == 0 && yo == 0 && ro == 0) {
    for (n = d / 8, i = 0; i < n; i++) {
      rp[i] = (*op)(xp[i], yp[i]);
    }
    if ((j = d % 8) > 0) {
      byte rpt = (*op)(xp[n], yp[n]);
      set_high8(rp[n], j, rpt);
    }
    if (!replace)
      return r;
  } else {
    for (n = d / 8, i = 0; i <= n; i++) {
      extract_byte8(xi, xp, i, xo);
      extract_byte8(yi, yp, i, yo);
      if (i == n) {
        if ((j = d % 8) == 0)
          break;
        extract_byte8(ri, rp, n, ro);
        set_high8(ri, j, (*op)(xi, yi));
      } else {
        ri = (*op)(xi, yi);
      }
      store_byte8(rp, i, ro, ri);
    }
    if (!replace)
      return r;
  }
  rp = r0->bytes();
  ro = startr0; // r0->offset();
  for (n = d / 8, i = 0; i <= n; i++) {
    if (i == n) {
      if ((j = d % 8) == 0)
        break;
      extract_byte8(ri, rp, n, ro);
      set_high8(ri, j, r->bytes()[n]);
    } else
      ri = r->bytes()[i];
    store_byte8(rp, i, ro, ri);
  }
  return r0;
ERROR:
  SIMPLE_ERROR(BF("Illegal arguments for bit-array operation."));
}

#endif // BIT_ARRAY_BYTE_SIZE==8


#if BIT_ARRAY_BYTE_SIZE==32

#define mask32 0xFFFFFFFF00000000
template <typename Place, typename Nbits, typename Value>
inline void set_high32(Place& place, Nbits nbits, Value value)
{
  (place) = ((place) & ~(mask32 >> (nbits))) | ((value) & (mask32>> (nbits)));
}

template <typename Place, typename Nbits, typename Value>
inline void set_low32(Place& place, Nbits nbits, Value value) {
  (place) = ((place) & (mask32 >> (32 - (nbits)))) | ((value) & ~(mask32 >> (32 - (nbits))));
}

template <typename Integer, typename Pointer, typename Index, typename Offset>
inline void extract_byte32(Integer& integer, Pointer pointer, Index index, Offset offset) {
  (integer) = (pointer)[(index)+1] & (~mask32);
  (integer) = ((pointer)[index] << (offset)) | ((integer) >> (32 - (offset)));
}

template <typename Pointer, typename Index, typename Offset, typename Value>
inline void store_byte32(Pointer pointer, Index index, Offset offset, Value value) {
  set_low32((pointer)[index], 32 - (offset), (value) >> (offset));
  set_high32((pointer)[(index)+1], offset, (value) << (32 - (offset)));
}

//#define TEMPLATE_BIT_ARRAY_OP 1
#ifndef TEMPLATE_BIT_ARRAY_OP

class Dispatcher
{
public:
    virtual Array_sp dispatcher (Array_O *target, SimpleBitVector_sp result) {
      //std::cout << "listclass is  " << lisp_classNameAsString(core::instance_class(target->asSmartPtr())) << '\n';
      // Und bist du nicht willig, so brauche ich Gewalt
      // Don't understand why target->create_result_bitarray(result) does not work
      SimpleBitVector_O * sbvptr = dynamic_cast<SimpleBitVector_O *>(target);
      MDArrayBit_O * mdabptr = dynamic_cast<MDArrayBit_O *>(target);
      SimpleMDArrayBit_O * smbabptr = dynamic_cast<SimpleMDArrayBit_O *>(target);
      BitVectorNs_O * bvnptr =  dynamic_cast<BitVectorNs_O *>(target);
      if (sbvptr)
        return result;
      else if (mdabptr)
        return mdabptr->create_result_bitarray(result);
      else if (smbabptr)
        return smbabptr->create_result_bitarray(result);
      else if (bvnptr)
        return bvnptr->create_result_bitarray(result);
      else
        // most likely return a SimpleBitVector_sp but perhaps better than crashing
        return result;
    }
};

// without these c++ seems always to call create_result_bitarray of Array_0
//Array_sp return_correct_bitvector (SimpleBitVector_sp target, SimpleBitVector_sp result) {
//  return (* target).create_result_bitarray(result);
//}


CL_LAMBDA(op x y &optional r);
CL_DECLARE();
CL_DOCSTRING("bitArrayOp");
CL_DEFUN T_sp core__bit_array_op(int opid, Array_sp tx, Array_sp ty, T_sp tr) {
  gctools::Fixnum i, j, n, d;
  SimpleBitVector_sp r0;
  size_t startr0 = 0;
  fixnum_bitwise_operator_t op;
  bool replace = false;
  byte64_t xi, yi, ri;
  byte32_t *xp, *yp, *rp;
  byte64_t xo, yo, ro;
  AbstractSimpleVector_sp ax;
  size_t startx, endx;
  AbstractSimpleVector_sp ay;
  size_t starty, endy;
  Dispatcher dispatcher;
  Array_sp array_x = gc::As<Array_sp>(tx);
  array_x->asAbstractSimpleVectorRange(ax, startx, endx);
  SimpleBitVector_sp x = gc::As_unsafe<SimpleBitVector_sp>(ax);
  Array_sp array_y = gc::As<Array_sp>(ty);
  array_y->asAbstractSimpleVectorRange(ay, starty, endy);
  SimpleBitVector_sp y = gc::As_unsafe<SimpleBitVector_sp>(ay);
  SimpleBitVector_sp r;
  size_t startr, endr;

  if (opid < 0 || opid >= bitArrayOperatorCount)
    goto ERROR;

  d = (endx - startx); // x->arrayTotalSize();
  xp = x->bytes();
  xo = startx; // x->offset();
  if (d != (endy - starty)) //(d != array_y->arrayTotalSize()) gives incorrect result for multi-dimensional arrays
    goto ERROR;
  yp = y->bytes();
  yo = starty; // y->offset();
  if (tr == _lisp->_true())
    tr = tx;
  if (tr.notnilp()) {
    AbstractSimpleVector_sp ar;
    Array_sp array_r = gc::As<Array_sp>(tr);
    array_r->asAbstractSimpleVectorRange(ar, startr, endr);
    r = gc::As_unsafe<SimpleBitVector_sp>(ar);
    if (!r) {
      ERROR_WRONG_TYPE_NTH_ARG(core::_sym_bitArrayOp, 4, tr, cl::_sym_SimpleBitVector_O);
    }
    if (endr-startr != d) //(r->arrayTotalSize() != d)
      goto ERROR;
    i = (r->bytes()-xp)*32+startr-xo;
    if ((i > 0 && i < d) || (i < 0 && -i < d)) {
      r0 = r;
      startr0 = startr;
      tr = _Nil<T_O>();
      replace = true;
      goto L1;
    }
    i = (r->bytes() - yp) * 32 + (startr - yo);
    if ((i > 0 && i < d) || (i < 0 && -i < d)) {
      r0 = r;
      startr0 = startr;
      tr = _Nil<T_O>();
      replace = true;
    }
  }
L1:
  if (tr.nilp()) {
    startr = 0;
    endr = d;
    r = SimpleBitVector_O::make(d);
  }
  rp = r->bytes();
  ro = startr; // r->offset();
  op = fixnum_bitwise_operators[opid];

  if (xo == 0 && yo == 0 && ro == 0) {
    for (n = d / 32, i = 0; i < n; i++) {
      rp[i] = (*op)(xp[i], yp[i]);
    }
    if ((j = d % 32) > 0) {
      byte64_t rpt = (*op)(xp[n], yp[n]);
      set_high32(rp[n], j, rpt);
    }
    if (!replace)
      return dispatcher.dispatcher(&(* tx),r);
  } else {
    for (n = d / 32, i = 0; i <= n; i++) {
      extract_byte32(xi, xp, i, xo);
      extract_byte32(yi, yp, i, yo);
      if (i == n) {
        if ((j = d % 32) == 0)
          break;
        extract_byte32(ri, rp, n, ro);
        set_high32(ri, j, (*op)(xi, yi));
      } else {
        ri = (*op)(xi, yi);
      }
      store_byte32(rp, i, ro, ri);
    }
    if (!replace)
      return dispatcher.dispatcher(&(* tx),r);
  }
  rp = r0->bytes();
  ro = startr0; // r0->offset();
  for (n = d / 32, i = 0; i <= n; i++) {
    if (i == n) {
      if ((j = d % 32) == 0)
        break;
      extract_byte32(ri, rp, n, ro);
      set_high32(ri, j, r->bytes()[n]);
    } else
      ri = r->bytes()[i];
    store_byte32(rp, i, ro, ri);
  }
  return dispatcher.dispatcher(&(* tx),r0);
ERROR:
  SIMPLE_ERROR(BF("Illegal arguments for bit-array operation."));
}

// As of now the lisp side calls bit_array_op directly, i.e. for now TEMPLATE_BIT_ARRAY_OP is not a drop-in replacement.

#define DO(name) CL_DEFUN T_sp core__bit_array_op_##name(T_sp tx, T_sp ty, T_sp tr) { \
    return core__bit_array_op(boole_##name, tx, ty, tr); \
  };

DO(and)
DO(andc1)
DO(andc2)
DO(nand)
DO(xor)
DO(ior)
DO(nor)
DO(orc1)
DO(orc2)
DO(eqv)
DO(c1)
DO(c2)
DO(1)
DO(2)
DO(clr)
DO(set)

#else // ifndef TEMPLATE_BIT_ARRAY_OP

template <boole_op op> struct do_bitwise_op {};

#define DO(name) \
  template <> struct do_bitwise_op<boole_##name> {static gc::Fixnum do_it(gc::Fixnum i, gc::Fixnum j) { return fixnum_bitwise_##name(i,j);};};

DO(and)
DO(andc1)
DO(andc2)
DO(nand)
DO(xor)
DO(ior)
DO(nor)
DO(orc1)
DO(orc2)
DO(eqv)
DO(c1)
DO(c2)
DO(1)
DO(2)
DO(clr)
DO(set)

template <boole_op OP>
T_sp template_bit_array_op(T_sp tx, T_sp ty, T_sp tr) {
  gctools::Fixnum i, j, n, d;
  SimpleBitVector_sp r0;
  size_t startr0 = 0;
  bool replace = false;
  byte64_t xi, yi, ri;
  byte32_t *xp, *yp, *rp;
  byte64_t xo, yo, ro;
  AbstractSimpleVector_sp ax;
  size_t startx, endx;
  AbstractSimpleVector_sp ay;
  size_t starty, endy;
  Array_sp array_x = gc::As<Array_sp>(tx);
  array_x->asAbstractSimpleVectorRange(ax, startx, endx);
  SimpleBitVector_sp x = gc::As_unsafe<SimpleBitVector_sp>(ax);
  Array_sp array_y = gc::As<Array_sp>(ty);
  array_y->asAbstractSimpleVectorRange(ay, starty, endy);
  SimpleBitVector_sp y = gc::As_unsafe<SimpleBitVector_sp>(ay);
  SimpleBitVector_sp r;
  size_t startr, endr;
  d = (endx - startx); // x->arrayTotalSize();
  xp = x->bytes();
  xo = startx; // x->offset();
  if (d != array_y->arrayTotalSize())
    goto ERROR;
  yp = y->bytes();
  yo = starty; // y->offset();
  if (tr == _lisp->_true())
    tr = x;
  if (tr.notnilp()) {
    AbstractSimpleVector_sp ar;
    Array_sp array_r = gc::As<Array_sp>(tr);
    array_r->asAbstractSimpleVectorRange(ar, startr, endr);
    r = gc::As_unsafe<SimpleBitVector_sp>(ar);
    if (!r) {
      ERROR_WRONG_TYPE_NTH_ARG(core::_sym_bitArrayOp, 4, tr, cl::_sym_SimpleBitVector_O);
    }
    if (endr-startr != d) //(r->arrayTotalSize() != d)
      goto ERROR;
    i = (r->bytes()-xp)*32+startr-xo;
    if ((i > 0 && i < d) || (i < 0 && -i < d)) {
      r0 = r;
      startr0 = startr;
      tr = _Nil<T_O>();
      replace = true;
      goto L1;
    }
    i = (r->bytes() - yp) * 32 + (startr - yo);
    if ((i > 0 && i < d) || (i < 0 && -i < d)) {
      r0 = r;
      startr0 = startr;
      tr = _Nil<T_O>();
      replace = true;
    }
  }
L1:
  if (tr.nilp()) {
    startr = 0;
    endr = d;
    r = SimpleBitVector_O::make(d);
  }
  rp = r->bytes();
  ro = startr; // r->offset();
  if (xo == 0 && yo == 0 && ro == 0) {
    for (n = d / 32, i = 0; i < n; i++) {
      rp[i] = do_bitwise_op<OP>::do_it(xp[i], yp[i]);
    }
    if ((j = d % 32) > 0) {
      byte64_t rpt = do_bitwise_op<OP>::do_it(xp[n], yp[n]);
      set_high32(rp[n], j, rpt);
    }
    if (!replace)
      return r;
  } else {
    for (n = d / 32, i = 0; i <= n; i++) {
      extract_byte32(xi, xp, i, xo);
      extract_byte32(yi, yp, i, yo);
      if (i == n) {
        if ((j = d % 32) == 0)
          break;
        extract_byte32(ri, rp, n, ro);
        set_high32(ri, j, do_bitwise_op<OP>::do_it(xi, yi));
      } else {
        ri = do_bitwise_op<OP>::do_it(xi, yi);
      }
      store_byte32(rp, i, ro, ri);
    }
    if (!replace)
      return r;
  }
  rp = r0->bytes();
  ro = startr0; // r0->offset();
  for (n = d / 32, i = 0; i <= n; i++) {
    if (i == n) {
      if ((j = d % 32) == 0)
        break;
      extract_byte32(ri, rp, n, ro);
      set_high32(ri, j, r->bytes()[n]);
    } else
      ri = r->bytes()[i];
    store_byte32(rp, i, ro, ri);
  }
  return r0;
ERROR:
  SIMPLE_ERROR(BF("Illegal arguments for bit-array operation."));
}

// Let's expand the template into functions that the lisp side can call

#define DO(name) \
  CL_DEFUN T_sp core__bit_array_op_##name(T_sp tx, T_sp ty, T_sp tr) { \
    return template_bit_array_op<boole_##name>(tx,ty,tr); \
  }

DO(and)
DO(andc1)
DO(andc2)
DO(nand)
DO(xor)
DO(ior)
DO(nor)
DO(orc1)
DO(orc2)
DO(eqv)
DO(c1)
DO(c2)
DO(1)
DO(2)
DO(clr)
DO(set)

#endif // TEMPLATE_BIT_ARRAY_OP
#endif // BIT_ARRAY_BYTE_SIZE==32

/*! Copied from ECL */
CL_DEFUN T_sp cl__logbitp(Integer_sp p, Integer_sp x) {
  bool i;
  if (p.fixnump()) {
    cl_index n = clasp_to_size(p);
    if (x.fixnump()) {
      gctools::Fixnum y = clasp_fixnum(x);
      if (n >= FIXNUM_BITS) {
        i = (y < 0);
      } else {
        i = ((y >> n) & 1);
      }
    } else {
      i = mpz_tstbit(gc::As<Bignum_sp>(x)->as_mpz_().get_mpz_t(), n);
    }
  } else {
    IMPLEMENT_MEF("Convert the code below to something Clasp can use");
#if 0
    assert_type_non_negative_integer(p);
    if (CLASP_FIXNUMP(x))
      i = (clasp_fixnum(x) < 0);
    else
      i = (_clasp_big_sign(x) < 0);
#endif
  }
  return i ? _lisp->_true() : _Nil<T_O>();
}

CL_LAMBDA(op arg1 arg2);
CL_DECLARE();
CL_DOCSTRING("boole");
CL_DEFUN T_sp cl__boole(T_sp op, T_sp arg1, T_sp arg2) {
  if (op.nilp()) {
    ERROR_WRONG_TYPE_NTH_ARG(cl::_sym_boole, 1, op, cl::_sym_integer);
  }
  Fixnum_sp fnop = gc::As<Fixnum_sp>(op);
  return clasp_boole(unbox_fixnum(fnop), arg1, arg2);
};

void initialize_bits() {

#define DO(name) \
  cl::_sym_boole_##name -> defconstant(make_fixnum(boole_##name)); \
  SYMBOL_EXPORT_SC_(ClPkg, boole_##name)

DO(and);
DO(andc1);
DO(andc2);
DO(nand);
DO(xor);
DO(ior);
DO(nor);
DO(orc1);
DO(orc2);
DO(eqv);
DO(c1);
DO(c2);
DO(1);
DO(2);
DO(clr);
DO(set);

//  af_def(ClPkg, "logbitp", &cl_logbitp);
};
};
