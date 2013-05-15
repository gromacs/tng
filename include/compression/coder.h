/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */

#ifndef CODER_H
#define CODER_H

struct coder
{
    unsigned int pack_temporary;
    int pack_temporary_bits;
    int stat_overflow;
    int stat_numval;
};

struct coder *Ptngc_coder_init(void);
void Ptngc_coder_deinit(struct coder *coder);
unsigned char *Ptngc_pack_array(struct coder *coder,int *input, int *length, int coding, int coding_parameter, int natoms, int speed);
int Ptngc_unpack_array(struct coder *coder,unsigned char *packed,int *output, int length, int coding, int coding_parameter, int natoms);
unsigned char *Ptngc_pack_array_xtc2(struct coder *coder,int *input, int *length);
int Ptngc_unpack_array_xtc2(struct coder *coder,unsigned char *packed,int *output, int length);
unsigned char *Ptngc_pack_array_xtc3(int *input, int *length, int natoms, int speed);
int Ptngc_unpack_array_xtc3(unsigned char *packed,int *output, int length, int natoms);

void Ptngc_out8bits(struct coder *coder, unsigned char **output);
void Ptngc_pack_flush(struct coder *coder,unsigned char **output);
void Ptngc_write_pattern(struct coder *coder,unsigned int pattern, int nbits, unsigned char **output);

void Ptngc_writebits(struct coder *coder,unsigned int value,int nbits, unsigned char **output_ptr);
void Ptngc_write32bits(struct coder *coder,unsigned int value,int nbits, unsigned char **output_ptr);
void Ptngc_writemanybits(struct coder *coder,unsigned char *value,int nbits, unsigned char **output_ptr);


#endif
