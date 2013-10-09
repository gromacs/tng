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


#include <stdio.h>
#include <stdlib.h>
#include "tng_compress.h"
#include "bwlzh.h"
#include "coder.h"
#include "warnmalloc.h"

#ifndef USE_WINDOWS
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define USE_WINDOWS
#endif /* win32... */
#endif /* not defined USE_WINDOWS */

#ifdef USE_WINDOWS
#define TNG_INLINE __inline
#define TNG_SNPRINTF _snprintf
#else
#define TNG_INLINE inline
#define TNG_SNPRINTF snprintf
#endif

struct coder *Ptngc_coder_init(void)
{
    struct coder *coder=warnmalloc(sizeof *coder);
    coder->pack_temporary_bits=0;
    return coder;
}

void Ptngc_coder_deinit(struct coder *coder)
{
    free(coder);
}

TNG_INLINE void Ptngc_out8bits(struct coder *coder, unsigned char **output)
{
  int pack_temporary_bits=coder->pack_temporary_bits;
  unsigned int pack_temporary=coder->pack_temporary;
  while (pack_temporary_bits>=8)
    {
      unsigned int mask=~(0xFFU<<(pack_temporary_bits-8));
      unsigned char out=(unsigned char)(pack_temporary>>(pack_temporary_bits-8));
      **output=out;
      (*output)++;
      pack_temporary_bits-=8;
      pack_temporary&=mask;
    }
  coder->pack_temporary_bits=pack_temporary_bits;
  coder->pack_temporary=pack_temporary;
}

void Ptngc_write_pattern(struct coder *coder,unsigned int pattern, int nbits, unsigned char **output)
{
    unsigned int mask1,mask2;
    mask1=1;
    mask2=1<<(nbits-1);
    coder->pack_temporary<<=nbits; /* Make room for new data. */
    coder->pack_temporary_bits+=nbits;
    while (nbits)
    {
	if (pattern & mask1)
	    coder->pack_temporary|=mask2;
	nbits--;
	mask1<<=1;
	mask2>>=1;
    }
    Ptngc_out8bits(coder,output);
}

/* Write up to 24 bits */
TNG_INLINE void Ptngc_writebits(struct coder *coder,unsigned int value,int nbits, unsigned char **output_ptr)
{
  /* Make room for the bits. */
  coder->pack_temporary<<=nbits;
  coder->pack_temporary_bits+=nbits;
  coder->pack_temporary|=value;
  Ptngc_out8bits(coder,output_ptr);
}

/* Write up to 32 bits */
void Ptngc_write32bits(struct coder *coder,unsigned int value,int nbits, unsigned char **output_ptr)
{
  unsigned int mask;
  if (nbits>=8)
    mask=0xFFU<<(nbits-8);
  else
    mask=0xFFU>>(8-nbits);
  while (nbits>8)
    {
      /* Make room for the bits. */
      coder->pack_temporary<<=8;
      coder->pack_temporary_bits+=8;
      coder->pack_temporary|=(value&mask)>>(nbits-8);
      Ptngc_out8bits(coder,output_ptr);
      nbits-=8;
      mask>>=8;
    }
  if (nbits)
    Ptngc_writebits(coder,value&mask,nbits,output_ptr);
}

/* Write "arbitrary" number of bits */
void Ptngc_writemanybits(struct coder *coder,unsigned char *value,int nbits, unsigned char **output_ptr)
{
  int vptr=0;
  while (nbits>=24)
    {
      unsigned int v=((((unsigned int)value[vptr])<<16)|
		      (((unsigned int)value[vptr+1])<<8)|
		      (((unsigned int)value[vptr+2])));
      Ptngc_writebits(coder,v,24,output_ptr);
      vptr+=3;
      nbits-=24;
    }
  while (nbits>=8)
    {
      Ptngc_writebits(coder,(unsigned int)value[vptr],8,output_ptr);
      vptr++;
      nbits-=8;
    }
  if (nbits)
    {
      Ptngc_writebits(coder,(unsigned int)value[vptr],nbits,output_ptr);
    }
}

static int write_stop_bit_code(struct coder *coder, unsigned int s,unsigned int coding_parameter, unsigned char **output)
{
  do {
    unsigned int extract=~(0xffffffffU<<coding_parameter);
    unsigned int this=(s&extract)<<1;
    s>>=coding_parameter;
    if (s)
      {
	this|=1U;
	coder->stat_overflow++;
      }
    coder->pack_temporary<<=(coding_parameter+1);
    coder->pack_temporary_bits+=coding_parameter+1;
    coder->pack_temporary|=this;
    Ptngc_out8bits(coder,output);
    if (s)
      {
	coding_parameter>>=1;
	if (coding_parameter<1)
	  coding_parameter=1;
      }
  } while (s);
  coder->stat_numval++;
  return 0;
}

static int pack_stopbits_item(struct coder *coder,int item, unsigned char **output, int coding_parameter)
{
    /* Find this symbol in table. */
    int s=0;
    if (item>0)
	s=1+(item-1)*2;
    else if (item<0)
	s=2+(-item-1)*2;
    return write_stop_bit_code(coder,s,coding_parameter,output);
    return 0;
}

static int pack_triplet(struct coder *coder,unsigned int *s, unsigned char **output, int coding_parameter,
			unsigned int max_base, int maxbits)
{
  /* Determine base for this triplet. */
  unsigned int min_base=1U<<coding_parameter;
  unsigned int this_base=min_base;
  int i;
  unsigned int jbase=0;
  unsigned int bits_per_value;
  for (i=0; i<3; i++)
    while (s[i]>=this_base)
      {
	this_base*=2;
	jbase++;
      }
  bits_per_value=coding_parameter+jbase;
  if (jbase>=3)
    {
      if (this_base>max_base)
	return 1;
      this_base=max_base;
      bits_per_value=maxbits;
      jbase=3;
    }
  /* 2 bits selects the base */
  coder->pack_temporary<<=2;
  coder->pack_temporary_bits+=2;
  coder->pack_temporary|=jbase;
  Ptngc_out8bits(coder,output);
  for (i=0; i<3; i++)
    Ptngc_write32bits(coder,s[i],bits_per_value,output);
  return 0;
}

void Ptngc_pack_flush(struct coder *coder,unsigned char **output)
{
  /* Zero-fill just enough. */
  if (coder->pack_temporary_bits>0)
    Ptngc_write_pattern(coder,0,8-coder->pack_temporary_bits,output);
}

unsigned char *Ptngc_pack_array(struct coder *coder,int *input, int *length, int coding, int coding_parameter,int natoms, int speed)
{
  if ((coding==TNG_COMPRESS_ALGO_BWLZH1) || (coding==TNG_COMPRESS_ALGO_BWLZH2))
    {
      unsigned char *output=warnmalloc(4+bwlzh_get_buflen(*length));
      int i,j,k,n=*length;
      unsigned int *pval=warnmalloc(n*sizeof *pval);
      int nframes=n/natoms/3;
      int cnt=0;
      int most_negative=2147483647;
      for (i=0; i<n; i++)
	if (input[i]<most_negative)
	  most_negative=input[i];
      most_negative=-most_negative;
      output[0]=((unsigned int)most_negative)&0xFFU;
      output[1]=(((unsigned int)most_negative)>>8)&0xFFU;
      output[2]=(((unsigned int)most_negative)>>16)&0xFFU;
      output[3]=(((unsigned int)most_negative)>>24)&0xFFU;
      for (i=0; i<natoms; i++)
	for (j=0; j<3; j++)
	  for (k=0; k<nframes; k++)
	    {
	      int item=input[k*3*natoms+i*3+j];
	      pval[cnt++]=(unsigned int)(item+most_negative);

	    }
      if (speed>=5)
	bwlzh_compress(pval,n,output+4,length);
      else
	bwlzh_compress_no_lz77(pval,n,output+4,length);
      (*length)+=4;
      free(pval);
      return output;
    }
  else if (coding==TNG_COMPRESS_ALGO_POS_XTC3)
    return Ptngc_pack_array_xtc3(input,length,natoms,speed);
  else if (coding==TNG_COMPRESS_ALGO_POS_XTC2)
    return Ptngc_pack_array_xtc2(coder,input,length);
  else
    {
      unsigned char *output=NULL;
      unsigned char *output_ptr=NULL;
      int i;
      int output_length=0;

      coder->stat_numval=0;
      coder->stat_overflow=0;
      /* Allocate enough memory for output */
      output=warnmalloc(8* *length*sizeof *output);
      output_ptr=output;
      if ((coding==TNG_COMPRESS_ALGO_TRIPLET) ||
	  (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA) ||
	  (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE))
	{
	  /* Pack triplets. */
	  int ntriplets=*length/3;
	  /* Determine max base and maxbits */
	  unsigned int max_base=1U<<coding_parameter;
	  unsigned int maxbits=coding_parameter;
	  unsigned int intmax=0;
	  for (i=0; i<*length; i++)
	    {
	      int item=input[i];
	      unsigned int s=0;
	      if (item>0)
		s=1+(item-1)*2;
	      else if (item<0)
		s=2+(-item-1)*2;
	      if (s>intmax)
		intmax=s;
	    }
	  /* Store intmax */
	  coder->pack_temporary_bits=32;
	  coder->pack_temporary=intmax;
	  Ptngc_out8bits(coder,&output_ptr);
	  while (intmax>=max_base)
	    {
	      max_base*=2;
	      maxbits++;
	    }
	  for (i=0; i<ntriplets; i++)
	    {
	      int j;
	      unsigned int s[3];
	      for (j=0; j<3; j++)
		{
		  int item=input[i*3+j];
		  /* Find this symbol in table. */
		  s[j]=0;
		  if (item>0)
		    s[j]=1+(item-1)*2;
		  else if (item<0)
		    s[j]=2+(-item-1)*2;
		}
	      if (pack_triplet(coder,s,&output_ptr,coding_parameter,max_base,maxbits))
		{
		  free(output);
		  return NULL;
		}
	    }
	}
      else
	for (i=0; i<*length; i++)
	  if (pack_stopbits_item(coder,input[i],&output_ptr,coding_parameter))
	    {
	      free(output);
	      return NULL;
	    }
      Ptngc_pack_flush(coder,&output_ptr);
      output_length=(int)(output_ptr-output);
      *length=output_length;
      return output;
    }
}

static int unpack_array_stop_bits(struct coder *coder,unsigned char *packed,int *output, int length, int coding_parameter)
{
  (void)coder;
  int i,j;
  unsigned int extract_mask=0x80;
  unsigned char *ptr=packed;
  for (i=0; i<length; i++)
    {
      unsigned int pattern=0;
      int numbits=coding_parameter;
      unsigned int bit;
      int s;
      unsigned int insert_mask=1U<<(numbits-1);
      int inserted_bits=numbits;
      do {
	for (j=0; j<numbits; j++)
	  {
	    bit=*ptr & extract_mask;
	    if (bit)
	      pattern|=insert_mask;
	    insert_mask>>=1;
	    extract_mask>>=1;
	    if (!extract_mask)
	      {
		extract_mask=0x80;
		ptr++;
	      }
	  }
	/* Check stop bit */
	bit=*ptr & extract_mask;
	extract_mask>>=1;
	if (!extract_mask)
	  {
	    extract_mask=0x80;
	    ptr++;
	  }
	if (bit)
	  {
	    numbits>>=1;
	    if (numbits<1)
	      numbits=1;
	    inserted_bits+=numbits;
	    insert_mask=1U<<(inserted_bits-1);
	  }
      } while (bit);
      s=(pattern+1)/2;
      if ((pattern%2)==0)
	s=-s;
      output[i]=s;
    }
  return 0;
}

static int unpack_array_triplet(struct coder *coder,unsigned char *packed,int *output, int length, int coding_parameter)
{
  (void)coder;
  int i,j;
  unsigned int extract_mask=0x80;
  unsigned char *ptr=packed;
  /* Determine max base and maxbits */
  unsigned int max_base=1U<<coding_parameter;
  unsigned int maxbits=coding_parameter;
  unsigned int intmax;
  /* Get intmax */
  intmax=((unsigned int)ptr[0])<<24|
    ((unsigned int)ptr[1])<<16|
    ((unsigned int)ptr[2])<<8|
    ((unsigned int)ptr[3]);
  ptr+=4;
  while (intmax>=max_base)
    {
      max_base*=2;
      maxbits++;
    }
  length/=3;
  for (i=0; i<length; i++)
    {
      /* Find base */
      unsigned int jbase=0;
      unsigned int numbits;
      unsigned int bit;
      for (j=0; j<2; j++)
	{
	  bit=*ptr & extract_mask;
	  jbase<<=1;
	  if (bit)
	    jbase|=1U;
	  extract_mask>>=1;
	  if (!extract_mask)
	    {
	      extract_mask=0x80;
	      ptr++;
	    }
	}
      if (jbase==3)
	numbits=maxbits;
      else
	numbits=coding_parameter+jbase;
      for (j=0; j<3; j++)
	{
	  int s;
	  unsigned int jbit;
	  unsigned int pattern=0;
	  for (jbit=0; jbit<numbits; jbit++)
	    {
	      bit=*ptr & extract_mask;
	      pattern<<=1;
	      if (bit)
		pattern|=1U;
	      extract_mask>>=1;
	      if (!extract_mask)
		{
		  extract_mask=0x80;
		  ptr++;
		}
	    }
	  s=(pattern+1)/2;
	  if ((pattern%2)==0)
	    s=-s;
	  output[i*3+j]=s;
	}
    }
  return 0;
}

static int unpack_array_bwlzh(struct coder *coder,unsigned char *packed,int *output, int length, int natoms)
{
  (void)coder;
  int i,j,k,n=length;
  unsigned int *pval=warnmalloc(n*sizeof *pval);
  int nframes=n/natoms/3;
  int cnt=0;
  int most_negative=(int)(((unsigned int)packed[0]) |
			  (((unsigned int)packed[1])<<8) |
			  (((unsigned int)packed[2])<<16) |
			  (((unsigned int)packed[3])<<24));
  bwlzh_decompress(packed+4,length,pval);
  for (i=0; i<natoms; i++)
    for (j=0; j<3; j++)
      for (k=0; k<nframes; k++)
	{
	  unsigned int s=pval[cnt++];
	  output[k*3*natoms+i*3+j]=(int)s-most_negative;
	}
  free(pval);
  return 0;
}

int Ptngc_unpack_array(struct coder *coder,unsigned char *packed,int *output, int length, int coding, int coding_parameter, int natoms)
{
  if ((coding==TNG_COMPRESS_ALGO_STOPBIT) ||
      (coding==TNG_COMPRESS_ALGO_VEL_STOPBIT_INTER))
    return unpack_array_stop_bits(coder, packed, output, length, coding_parameter);
  else if ((coding==TNG_COMPRESS_ALGO_TRIPLET) ||
	   (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_INTRA) ||
	   (coding==TNG_COMPRESS_ALGO_POS_TRIPLET_ONETOONE))
    return unpack_array_triplet(coder, packed, output, length, coding_parameter);
  else if (coding==TNG_COMPRESS_ALGO_POS_XTC2)
    return Ptngc_unpack_array_xtc2(coder, packed, output, length);
  else if ((coding==TNG_COMPRESS_ALGO_BWLZH1) || (coding==TNG_COMPRESS_ALGO_BWLZH2))
    return unpack_array_bwlzh(coder, packed, output, length,natoms);
  else if (coding==TNG_COMPRESS_ALGO_POS_XTC3)
    return Ptngc_unpack_array_xtc3(packed, output, length,natoms);
  return 1;
}

