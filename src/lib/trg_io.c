/* This code is part of the trg binary trajectory format.
 *
 *                      VERSION 1.0
 *
 * Written by Magnus Lundborg
 * Copyright (c) 2012, The GROMACS development team.
 * check out http://www.gromacs.org for more information.
 *
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */


/* TODO: Make sure UTF-8 works properly */

#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "trg_io.h"
#include "md5.h"

#ifndef LOGIN_NAME_MAX
#define LOGIN_NAME_MAX 33
#endif
#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX 257
#endif


/* This function swaps the byte order of a 32 bit numerical variable.
   It does not only work with integer, but e.g. floats need casting */
static inline trg_function_status trg_swap_byte_order_32
                (const struct trg_trajectory *trg_data, uint32_t *v)
{
    switch(trg_data->endianness_32)
    {
    case TRG_LITTLE_ENDIAN_32: // Byte order is reversed.
        *v = (*v << 24) |              // Move first byte to end
            ((*v & 0x00FF0000) >> 8) | // Move 2nd byte to pos 3
            ((*v & 0x0000FF00) << 8) | // Move 3rd byte to pos 2
            (*v >> 24);                // Move last byte to first

        return(TRG_SUCCESS);

    case TRG_BYTE_PAIR_SWAP_32: // byte pair swap
        *v = ((*v & 0xFFFF0000) >> 16) |
            ((*v & 0x0000FFFF) << 16);

        return(TRG_SUCCESS);

    case TRG_BIG_ENDIAN_32: // Already correct
        return(TRG_SUCCESS);

    default:
        return(TRG_FAILURE);
    }
}

/* This function swaps the byte order of a 64 bit numerical variable.
   It does not only work with integer, but e.g. floats need casting
   The byte order swapping routine can convert five different byte
   orders to big endian. */
static inline trg_function_status trg_swap_byte_order_64
                (const struct trg_trajectory *trg_data, uint64_t *v)
{
    switch(trg_data->endianness_64)
    {
    case TRG_LITTLE_ENDIAN_64: // Byte order is reversed.
        *v = (*v >> 56) |                       // Move first byte to end
            ((*v & 0x00FF000000000000) >> 40) | // Move 2nd byte to pos 7
            ((*v & 0x0000FF0000000000) >> 24) | // Move 3rd byte to pos 6
            ((*v & 0x000000FF00000000) >> 8 ) | // Move 4th byte to pos 5
            ((*v & 0x00000000FF000000) << 8 ) | // Move 5th byte to pos 4
            ((*v & 0x0000000000FF0000) << 24) | // Move 6th byte to pos 3
            ((*v & 0x000000000000FF00) << 40) | // Move 7th byte to pos 2
            (*v << 56);                         // Move last byte to first

        return(TRG_SUCCESS);

    case TRG_QUAD_SWAP_64: // Byte quad swap
        *v = ((*v & 0xFFFFFFFF00000000) >> 32) |
            ((*v & 0x00000000FFFFFFFF) << 32);

        return(TRG_SUCCESS);

    case TRG_BYTE_PAIR_SWAP_64: // Byte pair swap
        *v = ((*v & 0xFFFF0000FFFF0000) >> 16) |
            ((*v & 0x0000FFFF0000FFFF) << 16);

        return(TRG_SUCCESS);

    case TRG_BYTE_SWAP_64: // Byte swap
        *v = ((*v & 0xFF00FF00FF00FF00) >> 8) |
            ((*v & 0x00FF00FF00FF00FF) << 8);

        return(TRG_SUCCESS);

    case TRG_BIG_ENDIAN_64: // Already correct
        return(TRG_SUCCESS);

    default:
        return(TRG_FAILURE);
    }
}

/* Generate the md5 hash of a block.
   The hash is created based on the actual block contents. */
static trg_function_status trg_generate_block_hash(struct trg_gen_block *block)
{
    md5_state_t md5_state;

    md5_init(&md5_state);
    md5_append(&md5_state, (md5_byte_t *)block->block_contents,
               block->block_contents_size);
    md5_finish(&md5_state, (md5_byte_t *)block->hash);

    return(TRG_SUCCESS);
}

/* Compare the current block hash (e.g. read from file) with the hash
   calculated from the current contents.
   If the current hash is all zeros skip the comparison.
   If the hashes match results is set to TRUE, otherwise it is set to
   FALSE. */
static trg_function_status trg_verify_hash_match(struct trg_gen_block *block,
                                                 trg_bool *results)
{
    md5_state_t md5_state;
    char hash[TRG_HASH_LEN];

    if(strcmp(block->hash,"\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0") == 0)
    {
        *results = TRUE;
        return(TRG_SUCCESS);
    }
    md5_init(&md5_state);
    md5_append(&md5_state, (md5_byte_t *)block->block_contents,
               block->block_contents_size);
    md5_finish(&md5_state, (md5_byte_t *)hash);

    if(strncmp(block->hash, hash, 16) != 0)
    {
        *results = FALSE;
    }
    else
    {
        *results = TRUE;
    }
    
    return(TRG_SUCCESS);
}


static trg_function_status trg_init_input_file(struct trg_trajectory *trg_data,
                                               const trg_bool update_read_pos)
{
    if(!trg_data->input_file)
    {
        if(!trg_data->input_file_path)
        {
            printf("No file specified for reading. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        trg_data->input_file = fopen(trg_data->input_file_path, "r");
        if(!trg_data->input_file)
        {
            printf("Cannot open file %s. %s: %d\n",
                   trg_data->input_file_path, __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        if(fseek(trg_data->input_file, trg_data->input_file_pos,
                 SEEK_SET) != 0)
        {
            printf("Cannot specify position in file %s. %s: %d\n",
                   trg_data->input_file_path, __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }
    else if(update_read_pos && fseek(trg_data->input_file,
                                     trg_data->input_file_pos, SEEK_SET) != 0)
    {
        printf("Cannot specify position in file %s. %s: %d\n",
               trg_data->input_file_path, __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    return(TRG_SUCCESS);
}

static trg_function_status trg_init_output_file
                (struct trg_trajectory *trg_data,
                 const trg_bool update_write_pos)
{
    if(!trg_data->output_file)
    {
        if(!trg_data->output_file_path)
        {
            printf("No file specified for writing. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }

        if(trg_data->output_file_pos <= 0)
        {
            trg_data->output_file = fopen(trg_data->output_file_path, "w+");
        }
        else
        {
            trg_data->output_file = fopen(trg_data->output_file_path, "a+");
        }
        
        if(!trg_data->output_file)
        {
            printf("Cannot open file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        if(fseek(trg_data->output_file, 0, SEEK_SET) != 0)
        {
            printf("Cannot specify position in file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }
    else if(update_write_pos && fseek(trg_data->output_file, 0,
                                      SEEK_SET) != 0)
    {
        printf("Cannot specify position in file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    return(TRG_SUCCESS);
}

static trg_function_status trg_read_block_header
                (struct trg_trajectory *trg_data, struct trg_gen_block *block)
{
    int len, offset = 0;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        return(TRG_CRITICAL);
    }
    
    /* First read the header size to be able to read the whole header. */
    if(fread(&block->header_contents_size, sizeof(block->header_contents_size),
        1, trg_data->input_file) == 0)
    {
        printf("Cannot read header size. %s: %d\n",
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  &block->header_contents_size) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }

    /* Move the reading position to the beginning of the header. */
    fseek(trg_data->input_file, -sizeof(block->header_contents_size),
          SEEK_CUR);

    /* If there is already memory allocated for the contents free it (we do not
     * know if the size is correct). */
    if(block->header_contents)
    {
        free(block->header_contents);
    }

    block->header_contents = (char *) malloc(block->header_contents_size);
    if(!block->header_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->header_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* Read the whole header into header_contents. This way it can be saved
     * even if it cannot be interpreted
     * for one reason or another. */
    if(fread(block->header_contents, block->header_contents_size, 1,
        trg_data->input_file) == 0)
    {
        printf("Cannot read header. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* The header contents size has already been read. Skip ahead. */
    offset = sizeof(block->header_contents_size);


    /* Copy the respective parameters from the header contents block */
    memcpy(&block->block_contents_size, block->header_contents+offset,
           sizeof(block->block_contents_size));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &block->block_contents_size)
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->block_contents_size);

    memcpy(&block->id, block->header_contents+offset, sizeof(block->id));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &block->id) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->id);
    
    memcpy(block->hash, block->header_contents+offset, TRG_HASH_LEN);
    offset += TRG_HASH_LEN;

    if(block->name && strcmp(block->name, block->header_contents+offset) != 0)
    {
        free(block->name);
        block->name = 0;
    }
    len = min(strlen(block->header_contents+offset) + 1, TRG_MAX_STR_LEN);
    if(!block->name)
    {
        block->name = (char *) malloc(len);
        if(!block->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                    __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        strncpy(block->name, block->header_contents+offset, len);
    }
    offset += len;

    memcpy(&block->block_version, block->header_contents+offset,
           sizeof(block->block_version));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &block->block_version)
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->block_version);

    return(TRG_SUCCESS);
}

static trg_function_status trg_write_block_header
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block, write_mode mode)
{
    int name_len, offset = 0;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    

    /* Just dump the full header contents to file */
    if(mode == TRG_COPY_EXISTING)
    {
        if(!block->header_contents)
        {
            printf("No contents to write. %s: %d\n", __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->header_contents, block->header_contents_size, 1,
                  trg_data->output_file) != 1)
        {
            printf("Could not write all header data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    if(!block->name)
    {
        block->name = (char *) malloc(1);
        if(!block->name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        block->name[0] = 0;
    }

    name_len = min(strlen(block->name) + 1, TRG_MAX_STR_LEN);

    trg_generate_block_hash(block);

    /* Calculate the size of the header to write */
    block->header_contents_size = sizeof(block->header_contents_size) +
                                  sizeof(block->block_contents_size) +
                                  sizeof(block->id) +
                                  sizeof(block->block_version) +
                                  TRG_HASH_LEN +
                                  name_len;

    if(block->header_contents)
    {
        free(block->header_contents);
    }

    block->header_contents = (char *) malloc(block->header_contents_size);
    if(!block->header_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->header_contents_size, __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    /* First copy all data into the header_contents block and finally write
     * the whole block at once. */
    memcpy(block->header_contents, &block->header_contents_size,
           sizeof(block->header_contents_size));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  (uint64_t *)(block->header_contents))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->header_contents_size);
    
    memcpy(block->header_contents+offset, &block->block_contents_size,
           sizeof(block->block_contents_size));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  (uint64_t *)(block->header_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->block_contents_size);
    
    memcpy(block->header_contents+offset, &block->id, sizeof(block->id));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  (uint64_t *)(block->header_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->id);

    memcpy(block->header_contents+offset, block->hash, TRG_HASH_LEN);
    offset += TRG_HASH_LEN;

    strncpy(block->header_contents+offset, block->name, name_len);
    offset += name_len;

    memcpy(block->header_contents+offset, &block->block_version,
           sizeof(block->block_version));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  (uint64_t *)(block->header_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(block->block_version);

    if(fwrite(block->header_contents, block->header_contents_size,
       1, trg_data->output_file) != 1)
    {
        printf("Could not write all header data. %s: %d\n", __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    return(TRG_SUCCESS);
}


static trg_function_status trg_read_general_info_block
                (struct trg_trajectory *trg_data, struct trg_gen_block *block)
{
    int len, offset = 0;
    trg_bool same_hash;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        return(TRG_CRITICAL);
    }
    
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     * even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             trg_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    
    if(trg_verify_hash_match(block, &same_hash) != TRG_SUCCESS)
    {
        printf("Error comparing hashes. %s: %d\n", __FILE__, __LINE__);
        return(TRG_FAILURE);
    }
    if(same_hash != TRUE)
    {
        printf("General info block contents corrupt. Hashes do not match. "
               "%s: %d\n",
               __FILE__, __LINE__);
//         return(TRG_FAILURE);
    }

    len = min(strlen(block->block_contents) + 1, TRG_MAX_STR_LEN);
    trg_data->program_name = (char *) malloc(len);
    if(!trg_data->program_name)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }
    strncpy(trg_data->program_name, block->block_contents, len);
    offset += len;
    
    len = min(strlen(block->block_contents+offset) + 1, TRG_MAX_STR_LEN);
    trg_data->forcefield_name = (char *) malloc(len);
    if(!trg_data->forcefield_name)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }
    strncpy(trg_data->forcefield_name, block->block_contents+offset, len);
    offset += len;
    
    len = min(strlen(block->block_contents+offset) + 1, TRG_MAX_STR_LEN);
    trg_data->user_name = (char *) malloc(len);
    if(!trg_data->user_name)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }
    strncpy(trg_data->user_name, block->block_contents+offset, len);
    offset += len;
    
    memcpy(&trg_data->time, block->block_contents+offset,
           sizeof(trg_data->time));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &trg_data->time) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->time);

    len = min(strlen(block->block_contents+offset) + 1, TRG_MAX_STR_LEN);
    trg_data->computer_name = (char *) malloc(len);
    if(!trg_data->computer_name)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    strncpy(trg_data->computer_name, block->block_contents+offset, len);
    offset += len;

    len = min(strlen(block->block_contents+offset) + 1, TRG_MAX_STR_LEN);
    trg_data->pgp_signature = (char *) malloc(len);
    if(!trg_data->pgp_signature)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    strncpy(trg_data->pgp_signature, block->block_contents+offset, len);
    offset += len;

    memcpy(&trg_data->var_num_atoms_flag, block->block_contents+offset,
           sizeof(trg_data->var_num_atoms_flag));
    offset += sizeof(trg_data->var_num_atoms_flag);

    memcpy(&trg_data->frame_set_n_frames, block->block_contents+offset,
           sizeof(trg_data->frame_set_n_frames));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &trg_data->frame_set_n_frames)
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->frame_set_n_frames);
    
    memcpy(&trg_data->first_trajectory_frame_set_input_file_pos,
           block->block_contents+offset,
           sizeof(trg_data->first_trajectory_frame_set_input_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                          &trg_data->first_trajectory_frame_set_input_file_pos)
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->first_trajectory_frame_set_input_file_pos);
    trg_data->current_trajectory_frame_set.next_frame_set_file_pos =
    trg_data->first_trajectory_frame_set_input_file_pos;
    
    
    memcpy(&trg_data->last_trajectory_frame_set_input_file_pos,
           block->block_contents+offset,
           sizeof(trg_data->last_trajectory_frame_set_input_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                           &trg_data->last_trajectory_frame_set_input_file_pos)
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->last_trajectory_frame_set_input_file_pos);

    memcpy(&trg_data->stride_length, block->block_contents+offset,
           sizeof(trg_data->stride_length));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &trg_data->stride_length)
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }


    return(TRG_SUCCESS);
}

static trg_function_status trg_write_general_info_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block,
                 write_mode mode)
{
    int program_name_len, forcefield_name_len, user_name_len;
    int computer_name_len, pgp_signature_len;
    int offset = 0;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        return(TRG_CRITICAL);
    }
    
    if(!trg_data->program_name)
    {
        trg_data->program_name = (char *) malloc(1);
        if(!trg_data->program_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        trg_data->program_name[0] = 0;
    }
    if(!trg_data->forcefield_name)
    {
        trg_data->forcefield_name = (char *) malloc(1);
        if(!trg_data->forcefield_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        trg_data->forcefield_name[0] = 0;
    }
    if(!trg_data->user_name)
    {
        trg_data->user_name = (char *) malloc(1);
        if(!trg_data->user_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        trg_data->user_name[0] = 0;
    }
    if(!trg_data->computer_name)
    {
        trg_data->computer_name = (char *) malloc(1);
        if(!trg_data->computer_name)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        trg_data->computer_name[0] = 0;
    }
    if(!trg_data->pgp_signature)
    {
        trg_data->pgp_signature = (char *) malloc(1);
        if(!trg_data->pgp_signature)
        {
            printf("Cannot allocate memory (1 byte). %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        trg_data->pgp_signature[0] = 0;
    }

    program_name_len = min(strlen(trg_data->program_name) + 1,
                           TRG_MAX_STR_LEN);
    forcefield_name_len = min(strlen(trg_data->forcefield_name) + 1,
                              TRG_MAX_STR_LEN);
    user_name_len = min(strlen(trg_data->user_name) + 1,
                        TRG_MAX_STR_LEN);
    computer_name_len = min(strlen(trg_data->computer_name) + 1,
                            TRG_MAX_STR_LEN);
    pgp_signature_len = min(strlen(trg_data->pgp_signature) + 1,
                            TRG_MAX_STR_LEN);

    /* If just dumping the whole block_contents it is not certain that the
     * contents are known beforehand (e.g. due to different file versions) */
    if(mode == TRG_COPY_EXISTING)
    {
        if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
        {
            printf("Cannot write header of file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        if(!block->block_contents)
        {
            printf("No block data to write. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->block_contents, block->block_contents_size, 1,
                  trg_data->output_file) != 1)
        {
            printf("Could not write all block data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    block->block_contents_size = sizeof(trg_data->time) +
                sizeof(trg_data->var_num_atoms_flag) +
                sizeof(trg_data->frame_set_n_frames) +
                sizeof(trg_data->first_trajectory_frame_set_input_file_pos) +
                sizeof(trg_data->last_trajectory_frame_set_input_file_pos) +
                sizeof(trg_data->stride_length) +
                program_name_len +
                forcefield_name_len +
                user_name_len +
                computer_name_len +
                pgp_signature_len;

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    strncpy(block->block_contents, trg_data->program_name, program_name_len);
    offset += program_name_len;

    strncpy(block->block_contents+offset, trg_data->forcefield_name,
            forcefield_name_len);
    offset += forcefield_name_len;

    strncpy(block->block_contents+offset, trg_data->user_name, user_name_len);
    offset += user_name_len;

    memcpy(block->block_contents+offset, &trg_data->time,
           sizeof(trg_data->time));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->time);

    strncpy(block->block_contents+offset, trg_data->computer_name,
            computer_name_len);
    offset += computer_name_len;

    strncpy(block->block_contents+offset, trg_data->pgp_signature,
            pgp_signature_len);
    offset += pgp_signature_len;

    memcpy(block->block_contents+offset, &trg_data->var_num_atoms_flag,
           sizeof(trg_data->var_num_atoms_flag));
    offset += sizeof(trg_data->var_num_atoms_flag);
    
    memcpy(block->block_contents+offset, &trg_data->frame_set_n_frames,
           sizeof(trg_data->frame_set_n_frames));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->frame_set_n_frames);

    memcpy(block->block_contents+offset,
           &trg_data->first_trajectory_frame_set_input_file_pos,
           sizeof(trg_data->first_trajectory_frame_set_input_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->first_trajectory_frame_set_input_file_pos);

    memcpy(block->block_contents+offset,
           &trg_data->last_trajectory_frame_set_input_file_pos,
           sizeof(trg_data->last_trajectory_frame_set_input_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->last_trajectory_frame_set_input_file_pos);

    
    memcpy(block->block_contents+offset, &trg_data->stride_length,
           sizeof(trg_data->stride_length));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }


    if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
        trg_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    return(TRG_SUCCESS);
}

/* FIXME: Update this according to the new specs */
static trg_function_status trg_read_molecules_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block)
{
    int i, j, k, l, len, offset = 0;
    struct trg_molecule *molecule;
    struct trg_chain *chain;
    struct trg_residue *residue;
    struct trg_atom *atom;
    struct trg_bond *bond;
    trg_bool same_hash;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     * even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             trg_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(trg_verify_hash_match(block, &same_hash) != TRG_SUCCESS)
    {
        printf("Error comparing hashes. %s: %d\n", __FILE__, __LINE__);
        return(TRG_FAILURE);
    }
    if(same_hash != TRUE)
    {
        printf("Molecules block contents corrupt. Hashes do not match. "
               "%s: %d\n",
               __FILE__, __LINE__);
//         return(TRG_FAILURE);
    }

    memcpy(&trg_data->n_molecules, block->block_contents,
           sizeof(trg_data->n_molecules));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &trg_data->n_molecules)
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->n_molecules);

    if(trg_data->molecules)
    {
        free(trg_data->molecules);
    }

    trg_data->n_particles = 0;

    trg_data->molecules = (struct trg_molecule *)
                          malloc(trg_data->n_molecules *
                          sizeof(struct trg_molecule));
    if(!trg_data->molecules)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               trg_data->n_molecules * sizeof(struct trg_molecule),
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(!trg_data->var_num_atoms_flag)
    {
        if(trg_data->molecule_cnt_list)
        {
            free(trg_data->molecule_cnt_list);
        }
        trg_data->molecule_cnt_list = (uint64_t *) malloc(sizeof(uint64_t) *
                                      trg_data->n_molecules);
        if(!trg_data->molecule_cnt_list)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   trg_data->n_molecules * sizeof(struct trg_molecule),
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
    }
    

    for(i=0; i < trg_data->n_molecules; i++)
    {
        molecule = &trg_data->molecules[i];

        memcpy(&molecule->id, block->block_contents+offset,
               sizeof(molecule->id));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &molecule->id) != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->id);
        
//         printf("Read id: %"PRId64" offset: %d\n", molecule->id, offset);
        len = min(strlen(block->block_contents+offset) + 1, TRG_MAX_STR_LEN);
        molecule->name = (char *) malloc(len);
        strncpy(molecule->name, block->block_contents+offset, len);
        offset += len;

        memcpy(&molecule->quaternary_str, block->block_contents+offset,
               sizeof(molecule->quaternary_str));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &molecule->quaternary_str) != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->quaternary_str);

        if(!trg_data->var_num_atoms_flag)
        {
            memcpy(&trg_data->molecule_cnt_list[i],
                   block->block_contents+offset,
                   sizeof(uint64_t));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data,
                   &trg_data->molecule_cnt_list[i]) != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(uint64_t);
        }
        

        memcpy(&molecule->n_chains, block->block_contents+offset,
               sizeof(molecule->n_chains));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &molecule->n_chains)
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_chains);

        memcpy(&molecule->n_residues, block->block_contents+offset,
               sizeof(molecule->n_residues));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &molecule->n_residues)
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_residues);

        memcpy(&molecule->n_atoms, block->block_contents+offset,
               sizeof(molecule->n_atoms));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &molecule->n_atoms)
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_atoms);

        trg_data->n_particles += molecule->n_atoms *
                                 trg_data->molecule_cnt_list[i];

        molecule->chains = (struct trg_chain *) malloc(molecule->n_chains *
                                                      sizeof(struct trg_chain));
        if(!molecule->chains)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   molecule->n_chains * sizeof(struct trg_chain),
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        chain = molecule->chains;

        molecule->residues = (struct trg_residue *)
                             malloc(molecule->n_residues *
                             sizeof(struct trg_residue));
        if(!molecule->residues)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   molecule->n_residues * sizeof(struct trg_residue),
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        residue = molecule->residues;

        molecule->atoms = (struct trg_atom *) malloc(molecule->n_atoms *
                                                     sizeof(struct trg_atom));
        if(!molecule->atoms)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   molecule->n_atoms * sizeof(struct trg_atom),
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        atom = molecule->atoms;

        for(j=molecule->n_chains; j--;)
        {
            chain->molecule = molecule;
            
            memcpy(&chain->id, block->block_contents+offset,
                   sizeof(chain->id));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, &chain->id) != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                        __FILE__, __LINE__);
                }
            }
            offset += sizeof(chain->id);

            len = min(strlen(block->block_contents+offset) + 1,
                    TRG_MAX_STR_LEN);
            chain->name = (char *) malloc(len);
            strncpy(chain->name,
                    block->block_contents+offset, len);
            offset += len;

            memcpy(&chain->n_residues, block->block_contents+offset,
                sizeof(chain->n_residues));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, &chain->n_residues)
                    != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                        __FILE__, __LINE__);
                }
            }
            offset += sizeof(chain->n_residues);

            chain->residues = residue;
            for(k=chain->n_residues; k--;)
            {
                residue->chain = chain;
                memcpy(&residue->id, block->block_contents+offset,
                    sizeof(residue->id));
                if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                {
                    if(trg_swap_byte_order_64(trg_data, &residue->id) != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. %s: %d\n",
                            __FILE__, __LINE__);
                    }
                }
                offset += sizeof(residue->id);

                len = min(strlen(block->block_contents+offset) + 1,
                        TRG_MAX_STR_LEN);
                residue->name = (char *) malloc(len);
                strncpy(residue->name,
                        block->block_contents+offset, len);
                offset += len;

                memcpy(&residue->n_atoms, block->block_contents+offset,
                    sizeof(residue->n_atoms));
                if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                {
                    if(trg_swap_byte_order_64(trg_data, &residue->n_atoms)
                        != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. %s: %d\n",
                            __FILE__, __LINE__);
                    }
                }
                offset += sizeof(residue->n_atoms);

                residue->atoms = atom;
                for(l=residue->n_atoms; l--;)
                {
                    atom->residue = residue;
                    
                    memcpy(&atom->id, block->block_contents+offset,
                        sizeof(atom->id));
                    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                    {
                        if(trg_swap_byte_order_64(trg_data, &atom->id) != TRG_SUCCESS)
                        {
                            printf("Cannot swap byte order to get big endian. %s: %d\n",
                                __FILE__, __LINE__);
                        }
                    }
                    offset += sizeof(atom->id);

                    len = min(strlen(block->block_contents+offset) + 1,
                            TRG_MAX_STR_LEN);
                    atom->name = (char *) malloc(len);
                    strncpy(atom->name,
                            block->block_contents+offset, len);
                    offset += len;

                    len = min(strlen(block->block_contents+offset) + 1,
                            TRG_MAX_STR_LEN);
                    atom->atom_type = (char *) malloc(len);
                    strncpy(atom->atom_type,
                            block->block_contents+offset, len);
                    offset += len;

                    atom++;
                }
                residue++;
            }
            chain++;
        }

        memcpy(&molecule->n_bonds, block->block_contents+offset,
               sizeof(molecule->n_bonds));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &molecule->n_bonds)
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_bonds);
        
        molecule->bonds = (struct trg_bond *) malloc(molecule->n_bonds *
                                                     sizeof(struct trg_bond));
        if(!molecule->bonds)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   molecule->n_bonds * sizeof(struct trg_bond),
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        bond = molecule->bonds;

        for(j=molecule->n_bonds; j--;)
        {
            memcpy(&bond->from_atom_id, block->block_contents+offset,
                sizeof(bond->from_atom_id));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, &bond->from_atom_id)
                    != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->from_atom_id);

            memcpy(&bond->to_atom_id, block->block_contents+offset,
                sizeof(bond->to_atom_id));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, &bond->to_atom_id)
                    != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->to_atom_id);
            
            bond++;
        }
    }
    
    return(TRG_SUCCESS);
}

/* FIXME: Update this according to the new specs */
static trg_function_status trg_write_molecules_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block,
                 write_mode mode)
{
    int len = 0;
    int i, j, k, l, offset = 0;
    struct trg_molecule *molecule;
    struct trg_chain *chain;
    struct trg_residue *residue;
    struct trg_atom *atom;
    struct trg_bond *bond;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(mode != TRG_COPY_EXISTING && !trg_data->molecules)
    {
        return(TRG_SUCCESS);
    }

//     printf("Predicting block size\n");
    /* First predict the size of the block */
    for(i = 0; i < trg_data->n_molecules; i++)
    {
        molecule = &trg_data->molecules[i];
//         printf("mol %s\n", molecule->name);
        if(!molecule->name)
        {
            molecule->name = (char *) malloc(1);
            if(!molecule->name)
            {
                printf("Cannot allocate memory (1 byte). %s: %d\n",
                       __FILE__, __LINE__);
                trg_destroy_block(block);
                return(TRG_CRITICAL);
            }
            molecule->name[0] = 0;
        }
        len += min(strlen(molecule->name) + 1, TRG_MAX_STR_LEN);

        chain = molecule->chains;
        for(j = molecule->n_chains; j--;)
        {
            len += sizeof(chain->id);
            
            if(!chain->name)
            {
                chain->name = (char *) malloc(1);
                if(!chain->name)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    trg_destroy_block(block);
                    return(TRG_CRITICAL);
                }
                chain->name[0] = 0;
            }
            len += min(strlen(chain->name) + 1, TRG_MAX_STR_LEN);

            len += sizeof(chain->n_residues);

            chain++;
        }

        residue = molecule->residues;
        for(j = molecule->n_residues; j--;)
        {
            len += sizeof(residue->id);

            if(!residue->name)
            {
                residue->name = (char *) malloc(1);
                if(!residue->name)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    trg_destroy_block(block);
                    return(TRG_CRITICAL);
                }
                residue->name[0] = 0;
            }
            len += min(strlen(residue->name) + 1, TRG_MAX_STR_LEN);

            len += sizeof(residue->n_atoms);

            residue++;
        }

        atom = molecule->atoms;
        for(j = molecule->n_atoms; j--;)
        {
            len += sizeof(atom->id);
            if(!atom->name)
            {
                atom->name = (char *) malloc(1);
                if(!atom->name)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    trg_destroy_block(block);
                    return(TRG_CRITICAL);
                }
                atom->name[0] = 0;
            }
            len += min(strlen(atom->name) + 1, TRG_MAX_STR_LEN);

            if(!atom->atom_type)
            {
                atom->atom_type = (char *) malloc(1);
                if(!atom->atom_type)
                {
                    printf("Cannot allocate memory (1 byte). %s: %d\n",
                           __FILE__, __LINE__);
                    trg_destroy_block(block);
                    return(TRG_CRITICAL);
                }
                atom->atom_type[0] = 0;
            }
            len += min(strlen(atom->atom_type) + 1, TRG_MAX_STR_LEN);

            atom++;
        }

        bond = molecule->bonds;
        for(j = molecule->n_bonds; j--;)
        {
            len += sizeof(bond->from_atom_id) + sizeof(bond->to_atom_id);
        }
    }

    /* If just dumping the whole block_contents it is not certain that the
     * contents are known beforehand (e.g. due to different file versions) */
    if(mode == TRG_COPY_EXISTING)
    {
        if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
        {
            printf("Cannot write header of file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        if(!block->block_contents)
        {
            printf("No block data to write. %s: %d\n", __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->block_contents, block->block_contents_size, 1,
                  trg_data->output_file) != 1)
        {
            printf("Could not write all block data. %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    block->block_contents_size = sizeof(trg_data->n_molecules) +
                                 (sizeof(molecule->id) +
                                 sizeof(molecule->quaternary_str) +
                                 sizeof(molecule->n_chains) +
                                 sizeof(molecule->n_residues) +
                                 sizeof(molecule->n_atoms) +
                                 sizeof(molecule->n_bonds)) *
                                 trg_data->n_molecules +
                                 len;

    if(!trg_data->var_num_atoms_flag)
    {
        block->block_contents_size += trg_data->n_molecules * sizeof(uint64_t);
    }
    
    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    memcpy(block->block_contents+offset, &trg_data->n_molecules,
           sizeof(trg_data->n_molecules));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  (uint64_t *)(block->block_contents+offset))
            != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(trg_data->n_molecules);

    for(i = 0; i < trg_data->n_molecules; i++)
    {
        molecule = &trg_data->molecules[i];
//         printf("i=%d\n", i);
        memcpy(block->block_contents+offset, &molecule->id,
               sizeof(molecule->id));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data,
                                     (uint64_t *)(block->block_contents+offset))
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->id);

//         printf("Wrote id: %"PRId64" offset: %d\n", molecule->id, offset);
        len = min(strlen(molecule->name) + 1, TRG_MAX_STR_LEN);
        strncpy(block->block_contents + offset, molecule->name, len);
        offset += len;

        memcpy(block->block_contents+offset, &molecule->quaternary_str,
               sizeof(molecule->quaternary_str));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data,
                                     (uint64_t *)(block->block_contents+offset))
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->quaternary_str);
        
        if(!trg_data->var_num_atoms_flag)
        {
            memcpy(block->block_contents+offset,
                   &trg_data->molecule_cnt_list[i], sizeof(uint64_t));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data,
                                     (uint64_t *)(block->block_contents+offset))
                    != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(uint64_t);
        }

        memcpy(block->block_contents+offset, &molecule->n_chains,
               sizeof(molecule->n_chains));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data,
                                     (uint64_t *)(block->block_contents+offset))
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_chains);

        memcpy(block->block_contents+offset, &molecule->n_residues,
               sizeof(molecule->n_residues));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data,
                                     (uint64_t *)(block->block_contents+offset))
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_residues);

        memcpy(block->block_contents+offset, &molecule->n_atoms,
               sizeof(molecule->n_atoms));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data,
                                     (uint64_t *)(block->block_contents+offset))
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_atoms);

        chain = molecule->chains;
        for(j = molecule->n_chains; j--;)
        {
            memcpy(block->block_contents+offset, &chain->id, sizeof(chain->id));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data,
                                     (uint64_t *)(block->block_contents+offset))
                    != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(chain->id);

            len = min(strlen(chain->name) + 1, TRG_MAX_STR_LEN);
            strncpy(block->block_contents + offset, chain->name, len);
            offset += len;
            
            memcpy(block->block_contents+offset, &chain->n_residues,
                sizeof(chain->n_residues));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data,
                                        (uint64_t *)(block->block_contents+offset))
                    != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                        __FILE__, __LINE__);
                }
            }
            offset += sizeof(chain->n_residues);

            residue = chain->residues;
            for(k = chain->n_residues; k--;)
            {
                memcpy(block->block_contents+offset, &residue->id, sizeof(residue->id));
                if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                {
                    if(trg_swap_byte_order_64(trg_data,
                                        (uint64_t *)(block->block_contents+offset))
                        != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. %s: %d\n",
                            __FILE__, __LINE__);
                    }
                }
                offset += sizeof(residue->id);

                len = min(strlen(residue->name) + 1, TRG_MAX_STR_LEN);
                strncpy(block->block_contents + offset, residue->name, len);
                offset += len;
                
                memcpy(block->block_contents+offset, &residue->n_atoms,
                    sizeof(residue->n_atoms));
                if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                {
                    if(trg_swap_byte_order_64(trg_data,
                                            (uint64_t *)(block->block_contents+offset))
                        != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. %s: %d\n",
                            __FILE__, __LINE__);
                    }
                }
                offset += sizeof(residue->n_atoms);

                atom = residue->atoms;
                for(l = residue->n_atoms; l--;)
                {
        //             printf("j=%d\n", j);
                    memcpy(block->block_contents+offset, &atom->id, sizeof(atom->id));
                    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                    {
                        if(trg_swap_byte_order_64(trg_data,
                                            (uint64_t *)(block->block_contents+offset))
                            != TRG_SUCCESS)
                        {
                            printf("Cannot swap byte order to get big endian. %s: %d\n",
                                __FILE__, __LINE__);
                        }
                    }
                    offset += sizeof(atom->id);

                    len = min(strlen(atom->name) + 1, TRG_MAX_STR_LEN);
                    strncpy(block->block_contents + offset, atom->name, len);
                    offset += len;

                    len = min(strlen(atom->atom_type) + 1, TRG_MAX_STR_LEN);
                    strncpy(block->block_contents + offset, atom->atom_type, len);
                    offset += len;

                    atom++;
                }
                residue++;
            }
            chain++;
        }
                
        memcpy(block->block_contents+offset, &molecule->n_bonds,
               sizeof(molecule->n_bonds));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents+offset))
                != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(molecule->n_bonds);
        
        bond = molecule->bonds;
        for(j = molecule->n_bonds; j--;)
        {
            memcpy(block->block_contents+offset, &bond->from_atom_id,
                   sizeof(bond->from_atom_id));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                          (block->block_contents+offset))
                    != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->from_atom_id);
            
            memcpy(block->block_contents+offset, &bond->to_atom_id,
                   sizeof(bond->to_atom_id));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                    (block->block_contents+offset)) != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(bond->to_atom_id);

            bond++;
        }
    }

    if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              trg_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n",
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    return(TRG_SUCCESS);
}


static trg_function_status trg_read_frame_set_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block)
{
    int i, file_pos, offset = 0;
    uint64_t prev_n_particles;
    trg_bool same_hash;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    struct trg_particle_mapping *mapping;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to
     * disk even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             trg_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */

    if(trg_verify_hash_match(block, &same_hash) != TRG_SUCCESS)
    {
        printf("Error comparing hashes. %s: %d\n", __FILE__, __LINE__);
        return(TRG_FAILURE);
    }
    if(same_hash != TRUE)
    {
        printf("Frame set block contents corrupt. Hashes do not match. "
               "%s: %d\n",
               __FILE__, __LINE__);
//         return(TRG_FAILURE);
    }

    file_pos = ftell(trg_data->input_file);

    if(frame_set->n_mapping_blocks && frame_set->mappings)
    {
        for(i = frame_set->n_mapping_blocks; i--;)
        {
            mapping = &frame_set->mappings[i];
            if(mapping->real_particle_numbers)
            {
                free(mapping->real_particle_numbers);
                mapping->real_particle_numbers = 0;
            }
        }
        free(frame_set->mappings);
        frame_set->n_mapping_blocks = 0;
    }

    if(trg_data->first_trajectory_frame_set_input_file_pos <= 0)
    {
        trg_data->first_trajectory_frame_set_input_file_pos = file_pos;
    }
    /* FIXME: Should check the frame number instead of the file_pos, in case
     * frame sets are not in order */
    if(trg_data->last_trajectory_frame_set_input_file_pos < file_pos)
    {
        trg_data->last_trajectory_frame_set_input_file_pos = file_pos;
    }

    memcpy(&frame_set->first_frame, block->block_contents,
           sizeof(frame_set->first_frame));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &frame_set->first_frame) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->first_frame);
    
    memcpy(&frame_set->n_frames, block->block_contents + offset,
           sizeof(frame_set->n_frames));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &frame_set->n_frames) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->n_frames);

    if(trg_data->var_num_atoms_flag)
    {
        prev_n_particles = frame_set->n_particles;
        frame_set->n_particles = 0;
        /* If the list of molecule counts has already been created assume that
         * it is of correct size. */
        if(!frame_set->molecule_cnt_list)
        {
                frame_set->molecule_cnt_list =
                (uint64_t *) malloc(sizeof(uint64_t) * trg_data->n_molecules);
                
                if(!frame_set->molecule_cnt_list)
                {
                    printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                           sizeof(uint64_t) * trg_data->n_molecules,
                           __FILE__, __LINE__);
                    trg_destroy_block(block);
                    return(TRG_CRITICAL);
                }
        }
        for(i = 0; i < trg_data->n_molecules; i++)
        {
            memcpy(&frame_set->molecule_cnt_list[i],
                   block->block_contents + offset,
                   sizeof(uint64_t));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data,
                                          &frame_set->molecule_cnt_list[i]) !=
                   TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(uint64_t);
            frame_set->n_particles += trg_data->molecules[i].n_atoms *
                                      frame_set->molecule_cnt_list[i];
        }
        if(prev_n_particles && frame_set->n_particles != prev_n_particles)
        {
            /* FIXME: Particle dependent data memory management */
        }
    }

    memcpy(&frame_set->next_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->next_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  &frame_set->next_frame_set_file_pos) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->next_frame_set_file_pos);
    
    memcpy(&frame_set->prev_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->prev_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  &frame_set->prev_frame_set_file_pos) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->prev_frame_set_file_pos);

    memcpy(&frame_set->long_stride_next_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->long_stride_next_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  &frame_set->
                                  long_stride_next_frame_set_file_pos) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_next_frame_set_file_pos);

    memcpy(&frame_set->long_stride_prev_frame_set_file_pos,
           block->block_contents + offset,
           sizeof(frame_set->long_stride_prev_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data,
                                  &frame_set->
                                  long_stride_prev_frame_set_file_pos) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_prev_frame_set_file_pos);

//     frame_set->n_double_particle_data_blocks = 0;
//     frame_set->n_float_particle_data_blocks = 0;
//     frame_set->n_int_particle_data_blocks = 0;
//     frame_set->n_gen_particle_data_blocks = 0;
//     frame_set->n_double_data_blocks = 0;
//     frame_set->n_float_data_blocks = 0;
//     frame_set->n_int_data_blocks = 0;
//     frame_set->n_gen_data_blocks = 0;
    
    return(TRG_SUCCESS);
}

static trg_function_status trg_write_frame_set_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block,
                 write_mode mode)
{
    char *temp_name;
    uint64_t i;
    int offset = 0, name_len;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* If just dumping the whole block_contents it is not certain that the
     * contents are known beforehand (e.g. due to different file versions) */
    if(mode == TRG_COPY_EXISTING)
    {
        if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
        {
            printf("Cannot write header of file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        if(!block->block_contents)
        {
            printf("No block data to write. %s: %d\n", __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->block_contents, block->block_contents_size, 1,
                  trg_data->output_file) != 1)
        {
            printf("Could not write all block data. %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    name_len = strlen("TRAJECTORY FRAME SET");

    if(!block->name || strlen(block->name) < name_len)
    {
        temp_name = realloc(block->name, name_len + 1);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   name_len+1, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        block->name = temp_name;
    }
    strcpy(block->name, "TRAJECTORY FRAME SET");

    block->block_contents_size = sizeof(uint64_t) * 6;
    if(trg_data->var_num_atoms_flag)
    {
        block->block_contents_size += sizeof(uint64_t) * trg_data->n_molecules;
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    memcpy(block->block_contents, &frame_set->first_frame,
           sizeof(frame_set->first_frame));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->first_frame);

    memcpy(block->block_contents+offset, &frame_set->n_frames,
           sizeof(frame_set->n_frames));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->n_frames);

    if(trg_data->var_num_atoms_flag)
    {
        for(i = 0; i < trg_data->n_molecules; i++)
        {
            memcpy(block->block_contents+offset,
                   &frame_set->molecule_cnt_list[i],
                   sizeof(uint64_t));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                          (block->block_contents+offset)) !=
                   TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(uint64_t);
        }
    }

    
    memcpy(block->block_contents+offset, &frame_set->next_frame_set_file_pos,
           sizeof(frame_set->next_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->next_frame_set_file_pos);

    memcpy(block->block_contents+offset, &frame_set->prev_frame_set_file_pos,
           sizeof(frame_set->prev_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->prev_frame_set_file_pos);

    memcpy(block->block_contents+offset,
           &frame_set->long_stride_next_frame_set_file_pos,
           sizeof(frame_set->long_stride_next_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_next_frame_set_file_pos);

    memcpy(block->block_contents+offset,
           &frame_set->long_stride_prev_frame_set_file_pos,
           sizeof(frame_set->long_stride_prev_frame_set_file_pos));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(frame_set->long_stride_prev_frame_set_file_pos);

    if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              trg_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    return(TRG_SUCCESS);
}

static trg_function_status trg_read_trajectory_toc_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block)
{
    uint64_t i, old_n_blocks;
    int offset = 0, len;
    trg_bool same_hash;
    struct trg_frame_set_toc *toc =
    &trg_data->current_trajectory_frame_set.contents;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     *  even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             trg_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */


    if(trg_verify_hash_match(block, &same_hash) != TRG_SUCCESS)
    {
        printf("Error comparing hashes. %s: %d\n", __FILE__, __LINE__);
        return(TRG_FAILURE);
    }
    if(same_hash != TRUE)
    {
        printf("Table of contents block contents corrupt. Hashes do not match. "
               "%s: %d\n",
               __FILE__, __LINE__);
//         return(TRG_FAILURE);
    }

    old_n_blocks = toc->n_blocks;
    
    memcpy(&toc->n_blocks, block->block_contents,
           sizeof(toc->n_blocks));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &toc->n_blocks) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(toc->n_blocks);

    if(old_n_blocks != toc->n_blocks)
    {
        if(toc->block_names)
        {
            for(i = old_n_blocks; i--;)
            {
                if(toc->block_names[i])
                {
                    free(toc->block_names[i]);
                }
            }
            free(toc->block_names);
        }
        toc->block_names = (char **) malloc(toc->n_blocks * sizeof(char *));
        if(!toc->block_names)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   toc->n_blocks * sizeof(uint64_t), __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
    }

    for(i = 0; i < toc->n_blocks; i++)
    {
        len = min(strlen(block->block_contents+offset) + 1, TRG_MAX_STR_LEN);

        toc->block_names[i] = (char *) malloc(len);
        if(!toc->block_names[i])
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                    len, __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }

        strncpy(toc->block_names[i], block->block_contents+offset, len);
        
        offset += len;
    }
    
    return(TRG_SUCCESS);
}

static trg_function_status trg_write_trajectory_toc_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block,
                 write_mode mode)
{
    char *temp_name;
    uint64_t i;
    int offset = 0, name_len;
    struct trg_frame_set_toc *toc =
    &trg_data->current_trajectory_frame_set.contents;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* If just dumping the whole block_contents it is not certain that the
     * contents are known beforehand (e.g. due to different file versions) */
    if(mode == TRG_COPY_EXISTING)
    {
        if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
        {
            printf("Cannot write header of file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        if(!block->block_contents)
        {
            printf("No block data to write. %s: %d\n", __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->block_contents, block->block_contents_size, 1,
                  trg_data->output_file) != 1)
        {
            printf("Could not write all block data. %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    name_len = strlen("BLOCK TABLE OF CONTENTS");

    if(!block->name || strlen(block->name) < name_len)
    {
        temp_name = realloc(block->name, name_len + 1);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   name_len+1, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        block->name = temp_name;
    }
    strcpy(block->name, "BLOCK TABLE OF CONTENTS");

    block->block_contents_size = sizeof(uint64_t);

    for(i = toc->n_blocks; i--;)
    {
        block->block_contents_size += strlen(toc->block_names[i]) + 1;
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    memcpy(block->block_contents, &toc->n_blocks, sizeof(toc->n_blocks));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(toc->n_blocks);

    for(i = 0; i < toc->n_blocks; i++)
    {
        name_len = strlen(toc->block_names[i]) + 1;
        strncpy(block->block_contents+offset, toc->block_names[i], name_len);
        offset += name_len;
    }


    if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              trg_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n",
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    return(TRG_SUCCESS);
}

static trg_function_status trg_read_trajectory_mapping_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block)
{
    uint64_t i, old_n_particles;
    int offset = 0;
    trg_bool same_hash;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    
    struct trg_particle_mapping *mapping, *mappings;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to disk
     *  even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
        trg_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */


    if(trg_verify_hash_match(block, &same_hash) != TRG_SUCCESS)
    {
        printf("Error comparing hashes. %s: %d\n", __FILE__, __LINE__);
        return(TRG_FAILURE);
    }
    if(same_hash != TRUE)
    {
        printf("Particle mapping block contents corrupt. Hashes do not match. "
               "%s: %d\n",
               __FILE__, __LINE__);
//         return(TRG_FAILURE);
    }

    frame_set->n_mapping_blocks++;
    mappings = realloc(frame_set->mappings,
                       sizeof(struct trg_particle_mapping) *
                       frame_set->n_mapping_blocks);
    if(!mappings)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }
    frame_set->mappings = mappings;
    mapping = &mappings[frame_set->n_mapping_blocks - 1];
    

    memcpy(&mapping->num_first_particle, block->block_contents+offset,
           sizeof(mapping->num_first_particle));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &mapping->num_first_particle) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->num_first_particle);

    old_n_particles = mapping->n_particles;

    memcpy(&mapping->n_particles, block->block_contents+offset,
           sizeof(mapping->n_particles));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &mapping->n_particles) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->n_particles);

    if(old_n_particles != mapping->n_particles)
    {
        if(mapping->real_particle_numbers)
        {
            free(mapping->real_particle_numbers);
        }
        mapping->real_particle_numbers = malloc(mapping->n_particles *
                                                sizeof(uint64_t));
        if(!mapping->real_particle_numbers)
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                   mapping->n_particles * sizeof(uint64_t), __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
    }

    for(i = 0; i < mapping->n_particles; i++)
    {
        memcpy(&mapping->real_particle_numbers[i],
                block->block_contents + offset,
                sizeof(uint64_t));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data,
                                      &mapping->real_particle_numbers[i]) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(uint64_t);
    }

    return(TRG_SUCCESS);
}

static trg_function_status trg_write_trajectory_mapping_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block,
                 int mapping_block_nr,
                 write_mode mode)
{
    char *temp_name;
    int i, offset = 0, name_len;
    struct trg_particle_mapping *mapping =
    &trg_data->current_trajectory_frame_set.mappings[mapping_block_nr];

    if(mapping_block_nr >=
       trg_data->current_trajectory_frame_set.n_mapping_blocks)
    {
        printf("Mapping block index out of bounds. %s: %d\n",
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* If just dumping the whole block_contents it is not certain that the
     * contents are known beforehand (e.g. due to different file versions) */
    if(mode == TRG_COPY_EXISTING)
    {
        if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
        {
            printf("Cannot write header of file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        if(!block->block_contents)
        {
            printf("No block data to write. %s: %d\n", __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->block_contents, block->block_contents_size, 1,
                  trg_data->output_file) != 1)
        {
            printf("Could not write all block data. %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    name_len = strlen("PARTICLE MAPPING");

    if(!block->name || strlen(block->name) < name_len)
    {
        temp_name = realloc(block->name, name_len + 1);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   name_len+1, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        block->name = temp_name;
    }
    strcpy(block->name, "PARTICLE MAPPING");

    block->block_contents_size = sizeof(uint64_t) * (2 + mapping->n_particles);

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    memcpy(block->block_contents, &mapping->num_first_particle,
           sizeof(mapping->num_first_particle));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->num_first_particle);

    memcpy(block->block_contents+offset, &mapping->n_particles,
           sizeof(mapping->n_particles));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(mapping->n_particles);

    for(i = 0; i < mapping->n_particles; i++)
    {
        memcpy(block->block_contents+offset, &mapping->real_particle_numbers[i],
               sizeof(uint64_t));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents+offset)) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(uint64_t);
    }


    if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              trg_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    return(TRG_SUCCESS);
}

static trg_function_status trg_create_particle_data_block
                (struct trg_trajectory *trg_data,
                 const trg_block_type block_type_flag)
{
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    
    struct trg_particle_data *data;
    
    if(block_type_flag == TRG_TRAJECTORY_BLOCK)
    {
        frame_set->n_particle_data_blocks++;
        data = realloc(frame_set->tr_particle_data,
                    sizeof(struct trg_particle_data) *
                    frame_set->n_particle_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct trg_particle_data) *
                frame_set->n_particle_data_blocks,
                __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        frame_set->tr_particle_data = data;
        data = &frame_set->tr_particle_data[frame_set->
                                            n_particle_data_blocks - 1];
    }
    else
    {
        trg_data->n_particle_data_blocks++;
        data = realloc(trg_data->non_tr_particle_data,
                        sizeof(struct trg_particle_data) *
                        trg_data->n_particle_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                    sizeof(struct trg_particle_data) *
                    trg_data->n_particle_data_blocks,
                    __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        trg_data->non_tr_particle_data = data;
        data = &trg_data->non_tr_particle_data[trg_data->
                                               n_particle_data_blocks - 1];
    }

    return(TRG_SUCCESS);
}

trg_function_status trg_allocate_particle_data_mem
                (struct trg_trajectory  *trg_data,
                 struct trg_particle_data *data,
                 uint64_t n_frames,
                 const uint64_t n_particles,
                 const uint64_t n_values_per_frame)
{
    union data_values ***values;
    uint64_t i, j, k;

    for(i = data->n_frames; i--;)
    {
        for(j = n_particles; j--;)
        {
            free(data->values[i][j]);
        }
        free(data->values[i]);
    }
    data->n_frames = n_frames;
    n_frames = max(1, n_frames);
    data->n_values_per_frame = n_values_per_frame;
    values = (union data_values ***) realloc(data->values,
                                             sizeof(union data_values **) *
                                             n_frames);
    if(!values)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
            sizeof(union data_values **) * n_frames,
            __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    data->values = values;

    for(i = n_frames; i-- ;)
    {
        data->values[i] = (union data_values **)
                            malloc(sizeof(union data_values *) *
                            n_particles);
        if(!data->values[i])
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                sizeof(union data_values *) * n_particles,
                __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        for(j = n_particles; j--;)
        {
            data->values[i][j] = (union data_values *)
                                    malloc(sizeof(union data_values) *
                                    n_values_per_frame);
            if(!data->values[i][j])
            {
                printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                    sizeof(union data_values) * n_values_per_frame,
                    __FILE__, __LINE__);
                return(TRG_CRITICAL);
            }
            if(data->datatype == TRG_CHAR_DATA)
            {
                for(k = n_values_per_frame; k--;)
                {
                    data->values[i][j][k].c = 0;
                }
            }
        }
    }
    return(TRG_SUCCESS);
}

static trg_function_status trg_read_particle_data
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block,
                 int *offset,
                 const char datatype,
                 const uint64_t first_particle_number,
                 const uint64_t n_particles,
                 const uint64_t first_frame_with_data,
                 const uint64_t stride_length,
                 uint64_t n_frames,
                 const uint64_t n_values,
                 const uint64_t codec_id,
                 const uint64_t multiplier)
{
    uint64_t block_index, i, j, k, tot_n_particles;
    int size, len;
    struct trg_particle_data *data;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    trg_block_type block_type_flag;

    switch(datatype)
    {
    case TRG_CHAR_DATA:
        size = 1;
        break;
    case TRG_INT_DATA:
        size = sizeof(uint64_t);
        break;
    case TRG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TRG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }

    if(trg_data->current_trajectory_frame_set_input_file_pos > 0)
    {
        block_type_flag = TRG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TRG_NON_TRAJECTORY_BLOCK;
    }
    
    block_index = -1;
    /* See if there is already a data block of this ID */
    if(block_type_flag == TRG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_particle_data_blocks; i-- ;)
        {
            data = &frame_set->tr_particle_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        for(i = trg_data->n_particle_data_blocks; i-- ;)
        {
            data = &trg_data->non_tr_particle_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }

    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(trg_create_particle_data_block(trg_data, block_type_flag) !=
            TRG_SUCCESS)
        {
            printf("Cannot create particle data block. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        if(block_type_flag == TRG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_particle_data[frame_set->
                                                n_particle_data_blocks - 1];
        }
        else
        {
            data = &trg_data->non_tr_particle_data[trg_data->
                                                   n_particle_data_blocks - 1];
        }
        data->block_id = block->id;
        
        data->block_name = malloc(strlen(block->name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block->name)+1, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        strcpy(data->block_name, block->name);

        data->datatype = datatype;

        data->values = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
        data->compression_multiplier = multiplier;
    }

    if(block_type_flag == TRG_TRAJECTORY_BLOCK &&
       trg_data->var_num_atoms_flag)
    {
        tot_n_particles = frame_set->n_particles;
    }
    else
    {
        tot_n_particles = trg_data->n_particles;
    }
    
    /* Allocate memory */
    if(!data->values || data->n_frames != n_frames ||
       data->n_values_per_frame != n_values)
    {
        if(trg_allocate_particle_data_mem(trg_data, data, n_frames,
                                          tot_n_particles, n_values) !=
           TRG_SUCCESS)
        {
            printf("Cannot allocate memory for particle data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    data->first_frame_with_data = first_frame_with_data;
    data->stride_length = stride_length;
    
    for(i = 0; i < n_frames; i++)
    {
        for(j = first_particle_number; j < n_particles; j++)
        {
            for(k = 0; k < n_values; k++)
            {
                /* FIXME: Not optimal to have this switch statement inside
                 * the for loops. Check if it should be moved outside. */
                switch(datatype)
                {
                case TRG_FLOAT_DATA:
                    memcpy(&data->values[i][j][k].f,
                           block->block_contents+*offset,
                        size);
                    if(trg_data->endianness_32 != TRG_BIG_ENDIAN_32)
                    {
                        if(trg_swap_byte_order_32(trg_data,
                        (uint32_t *) &data->values[i][j][k]) != TRG_SUCCESS)
                        {
                            printf("Cannot swap byte order to get big endian. "
                                   "%s: %d\n",
                                   __FILE__, __LINE__);
                        }
                    }
                    *offset += size;
                    break;
                case TRG_INT_DATA:
                    memcpy(&data->values[i][j][k].i,
                           block->block_contents+*offset,
                        size);
                    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                    {
                        if(trg_swap_byte_order_64(trg_data,
                        (uint64_t *) &data->values[i][j][k].i) != TRG_SUCCESS)
                        {
                            printf("Cannot swap byte order to get big endian. "
                                   "%s: %d\n",
                                   __FILE__, __LINE__);
                        }
                    }
                    *offset += size;
                    break;
                case TRG_CHAR_DATA:
                    len = min(strlen(block->block_contents+*offset) + 1,
                              TRG_MAX_STR_LEN);
                    if(data->values[i][j][k].c)
                    {
                        free(data->values[i][j][k].c);
                    }
                    data->values[i][j][k].c = malloc(len);
                    if(!data->values[i][j][k].c)
                    {
                        printf("Cannot allocate memory (%d bytes). %s: %d\n",
                            len, __FILE__, __LINE__);
                        return(TRG_CRITICAL);
                    }
                    strncpy(data->values[i][j][k].c,
                            block->block_contents+*offset, len);
                    *offset += len;
                    break;
                case TRG_DOUBLE_DATA:
                default:
                    memcpy(&data->values[i][j][k].d,
                           block->block_contents+*offset,
                        size);
                    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                    {
                        if(trg_swap_byte_order_64(trg_data,
                        (uint64_t *) &data->values[i][j][k].d) != TRG_SUCCESS)
                        {
                            printf("Cannot swap byte order to get big endian. "
                                   "%s: %d\n",
                                   __FILE__, __LINE__);
                        }
                    }
                    *offset += size;
                }
            }
        }
    }
    return(TRG_SUCCESS);
}


static trg_function_status trg_write_particle_data_block
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block,
                 const int block_index,
                 const struct trg_particle_mapping *mapping,
                 const write_mode mode)
{
    uint64_t n_particles, num_first_particle, n_frames;
    int i, j, k, offset = 0, size, len;
    char temp, *temp_name;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    
    struct trg_particle_data *data;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(trg_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        data = &frame_set->tr_particle_data[block_index];
    }
    else
    {
        data = &trg_data->non_tr_particle_data[block_index];
    }

    switch(data->datatype)
    {
    case TRG_CHAR_DATA:
        size = 1;
        break;
    case TRG_INT_DATA:
        size = sizeof(uint64_t);
        break;
    case TRG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TRG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }

    len = strlen(data->block_name) + 1;

    /* If just dumping the whole block_contents it is not certain that the
     * contents are known beforehand (e.g. due to different file versions) */
    if(mode == TRG_COPY_EXISTING)
    {
        if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
        {
            printf("Cannot write header of file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        if(!block->block_contents)
        {
            printf("No block data to write. %s: %d\n", __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->block_contents, block->block_contents_size, 1,
            trg_data->output_file) != 1)
        {
            printf("Could not write all block data. %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    if(!block->name || strlen(block->name) < len)
    {
        temp_name = realloc(block->name, len);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        block->name = temp_name;
    }
    strncpy(block->name, data->block_name, len);

    /* If writing frame independent data data->n_frames is be 0, but n_frames
       is used for the loop writing the data (and reserving memory) and needs
       to be at least 1 */
    n_frames = max(1, data->n_frames);

    
    if(mapping && mapping->n_particles != 0)
    {
        n_particles = mapping->n_particles;
        num_first_particle = mapping->num_first_particle;
    }
    else
    {
        num_first_particle = 0;
        if(trg_data->var_num_atoms_flag)
        {
            n_particles = frame_set->n_particles;
        }
        else
        {
            n_particles = trg_data->n_particles;
        }
    }
    
    block->block_contents_size = sizeof(char) * 3 +
                                 sizeof(data->n_values_per_frame) +
                                 sizeof(data->codec_id) +
                                 sizeof(num_first_particle) +
                                 sizeof(n_particles);

    if(data->codec_id != TRG_UNCOMPRESSED)
    {
        block->block_contents_size += sizeof(data->compression_multiplier);
    }

    if(data->datatype == TRG_CHAR_DATA)
    {
        for(i = n_frames; i--;)
        {
            for(j = num_first_particle; j < n_particles; j++)
            {
                for(k = data->n_values_per_frame; k--;)
                {
                    block->block_contents_size +=
                    strlen(data->values[i][j][k].c) + 1;
                }
            }
        }
    }
    else
    {
        block->block_contents_size += size * n_frames * n_particles *
                                      data->n_values_per_frame;
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }


    memcpy(block->block_contents, &data->datatype, sizeof(char));
    offset += sizeof(char);

    if(data->n_frames > 0)
    {
        temp = TRG_FRAME_DEPENDENT + TRG_PARTICLE_DEPENDENT;
    }
    else
    {
        temp = TRG_PARTICLE_DEPENDENT;
    }
    memcpy(block->block_contents+offset, &temp, sizeof(char));
    offset += sizeof(char);

    if(data->n_frames > 0 && data->stride_length > 1)
    {
        temp = 1;
    }
    else
    {
        temp = 0;
    }
    memcpy(block->block_contents+offset, &temp, sizeof(char));
    offset += sizeof(char);

    memcpy(block->block_contents+offset, &data->n_values_per_frame,
           sizeof(data->n_values_per_frame));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents + offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->n_values_per_frame);

    memcpy(block->block_contents+offset, &data->codec_id,
           sizeof(data->codec_id));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents + offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->codec_id);

    if(data->codec_id != TRG_UNCOMPRESSED)
    {
        memcpy(block->block_contents+offset, &data->compression_multiplier,
               sizeof(data->compression_multiplier));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents + offset)) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->compression_multiplier);
    }

    if(data->n_frames > 0 && data->stride_length > 1)
    {
        memcpy(block->block_contents+offset, &data->first_frame_with_data,
               sizeof(data->first_frame_with_data));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents + offset)) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->first_frame_with_data);

        memcpy(block->block_contents+offset, &data->stride_length,
               sizeof(data->stride_length));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents + offset)) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->stride_length);
    }

    
    memcpy(block->block_contents+offset, &num_first_particle,
           sizeof(num_first_particle));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(num_first_particle);

    memcpy(block->block_contents+offset, &n_particles, sizeof(n_particles));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents+offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(n_particles);

    
    for(i = 0; i < data->n_frames; i++)
    {
        for(j = num_first_particle; j < n_particles; j++)
        {
            for(k = 0; k < data->n_values_per_frame; k++)
            {
                switch(data->datatype)
                {
                case TRG_FLOAT_DATA:
                    memcpy(block->block_contents+offset,
                           &data->values[i][j][k].f,
                            size);
                    break;
                case TRG_INT_DATA:
                    memcpy(block->block_contents+offset,
                           &data->values[i][j][k].i,
                            size);
                    break;
                case TRG_CHAR_DATA:
                    len = strlen(data->values[i][j][k].c) + 1;
                    strncpy(block->block_contents+offset,
                            data->values[i][j][k].c, len);
                    offset += len;
                    break;
                case TRG_DOUBLE_DATA:
                default:
                    memcpy(block->block_contents+offset,
                           &data->values[i][j][k].d,
                            size);
                }
                if(data->datatype == TRG_FLOAT_DATA)
                {
                    if(trg_data->endianness_32 != TRG_BIG_ENDIAN_32)
                    {
                        if(trg_swap_byte_order_32(trg_data,
                            (uint32_t *)(block->block_contents+offset))
                            != TRG_SUCCESS)
                        {
                            printf("Cannot swap byte order to get big endian. "
                                   "%s: %d\n",
                                   __FILE__, __LINE__);
                        }
                    }
                    offset += sizeof(size);
                }
                else if(data->datatype == TRG_INT_DATA ||
                        data->datatype == TRG_DOUBLE_DATA)
                {
                    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                    {
                        if(trg_swap_byte_order_64(trg_data,
                            (uint64_t *)(block->block_contents+offset))
                            != TRG_SUCCESS)
                        {
                            printf("Cannot swap byte order to get big endian. "
                                   "%s: %d\n",
                                   __FILE__, __LINE__);
                        }
                    }
                    offset += sizeof(size);
                }
            }
        }
    }



    if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
        trg_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    return(TRG_SUCCESS);
}

/* UNTESTED */
static trg_function_status trg_create_data_block
                (struct trg_trajectory *trg_data,
                 const trg_block_type block_type_flag)
{
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    
    struct trg_data *data;
    
    if(block_type_flag == TRG_TRAJECTORY_BLOCK)
    {
        frame_set->n_data_blocks++;
        data = realloc(frame_set->tr_data, sizeof(struct trg_data) *
                        frame_set->n_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct trg_data) * frame_set->n_data_blocks,
                __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        frame_set->tr_data = data;
        data = &frame_set->tr_data[frame_set->n_data_blocks - 1];
    }
    else
    {
        trg_data->n_data_blocks++;
        data = realloc(trg_data->non_tr_data, sizeof(struct trg_data) *
                        trg_data->n_data_blocks);
        if(!data)
        {
            printf("Cannot allocate memory (%lu bytes). %s: %d\n",
                sizeof(struct trg_data) * trg_data->n_data_blocks,
                __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        trg_data->non_tr_data = data;
        data = &trg_data->non_tr_data[trg_data->n_data_blocks - 1];
    }

    return(TRG_SUCCESS);
}

/* UNTESTED */
trg_function_status trg_allocate_data_mem
                (struct trg_trajectory  *trg_data,
                 struct trg_data *data,
                 uint64_t n_frames,
                 const uint64_t n_values_per_frame)
{
    union data_values **values;
    uint64_t i, j;

    for(i = data->n_frames; i--;)
    {
        if(data->datatype == TRG_CHAR_DATA)
        {
            for(j = data->n_values_per_frame; j--;)
            {
                if(data->values[i][j].c)
                {
                    free(data->values[i][j].c);
                    data->values[i][j].c = 0;
                }
            }
        }
        free(data->values[i]);
    }
    data->n_frames = n_frames;
    n_frames = max(1, n_frames);
    data->n_values_per_frame = n_values_per_frame;
    values = (union data_values **) realloc(data->values,
                                             sizeof(union data_values *) *
                                             n_frames);
    if(!values)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
            sizeof(union data_values **) * n_frames,
            __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    data->values = values;

    for(i = n_frames; i-- ;)
    {
        data->values[i] = (union data_values *)
                           malloc(sizeof(union data_values) *
                           n_values_per_frame);
        if(!data->values[i])
        {
            printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
                sizeof(union data_values) * n_values_per_frame,
                __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        if(data->datatype == TRG_CHAR_DATA)
        {
            for(j = n_values_per_frame; j--;)
            {
                data->values[i][j].c = 0;
            }
        }
    }
    return(TRG_SUCCESS);
}

static trg_function_status trg_read_data(struct trg_trajectory *trg_data,
                                         struct trg_gen_block *block,
                                         int *offset,
                                         const char datatype,
                                         const uint64_t first_frame_with_data,
                                         const uint64_t stride_length,
                                         uint64_t n_frames,
                                         const uint64_t n_values,
                                         const uint64_t codec_id,
                                         const uint64_t multiplier)
{
    uint64_t block_index, i, j;
    int size, len;
    struct trg_data *data;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    trg_block_type block_type_flag;

//     printf("%s\n", block->name);

    if(trg_data->current_trajectory_frame_set_input_file_pos > 0)
    {
        block_type_flag = TRG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TRG_NON_TRAJECTORY_BLOCK;
    }

    switch(datatype)
    {
    case TRG_CHAR_DATA:
        size = 1;
        break;
    case TRG_INT_DATA:
        size = sizeof(uint64_t);
        break;
    case TRG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TRG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }
    
    block_index = -1;
    /* See if there is already a data block of this ID */
    if(block_type_flag == TRG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_data_blocks; i-- ;)
        {
            data = &frame_set->tr_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        for(i = trg_data->n_data_blocks; i-- ;)
        {
            data = &trg_data->non_tr_data[i];
            if(data->block_id == block->id)
            {
                block_index = i;
                break;
            }
        }
    }
    
    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(trg_create_data_block(trg_data, block_type_flag) !=
            TRG_SUCCESS)
        {
            printf("Cannot create particle data block. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        if(block_type_flag == TRG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_data[frame_set->n_data_blocks - 1];
        }
        else
        {
            data = &trg_data->non_tr_data[trg_data->n_data_blocks - 1];
        }
        data->block_id = block->id;

        data->block_name = malloc(strlen(block->name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block->name)+1, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        strcpy(data->block_name, block->name);

        data->datatype = datatype;

        data->values = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
        data->compression_multiplier = multiplier;
    }

    /* Allocate memory */
    if(!data->values || data->n_frames != n_frames ||
       data->n_values_per_frame != n_values)
    {
        if(trg_allocate_data_mem(trg_data, data, n_frames, n_values) !=
           TRG_SUCCESS)
        {
            printf("Cannot allocate memory for data. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    data->first_frame_with_data = first_frame_with_data;
    data->stride_length = stride_length;
    
    n_frames = max(1, n_frames);

    for(i = 0; i < n_frames; i++)
    {
        for(j = 0; j < n_values; j++)
        {
            switch(datatype)
            {
            case TRG_CHAR_DATA:
                len = min(strlen(block->block_contents+*offset) + 1,
                          TRG_MAX_STR_LEN);
                if(data->values[i][j].c)
                {
                    free(data->values[i][j].c);
                }
                data->values[i][j].c = malloc(len);
                if(!data->values[i][j].c)
                {
                    printf("Cannot allocate memory (%d bytes). %s: %d\n",
                           len, __FILE__, __LINE__);
                    return(TRG_CRITICAL);
                }
                strncpy(data->values[i][j].c, block->block_contents+*offset,
                        len);
                *offset += len;
                break;
            case TRG_INT_DATA:
                memcpy(&data->values[i][j].i, block->block_contents+*offset,
                       size);
                if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                {
                    if(trg_swap_byte_order_64(trg_data,
                       (uint64_t *) &data->values[i][j].i) != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. "
                               "%s: %d\n",
                               __FILE__, __LINE__);
                    }
                }
                *offset += size;
                break;
            case TRG_FLOAT_DATA:
                memcpy(&data->values[i][j].f, block->block_contents+*offset,
                       size);
                if(trg_data->endianness_32 != TRG_BIG_ENDIAN_32)
                {
                    if(trg_swap_byte_order_32(trg_data,
                    (uint32_t *) &data->values[i][j]) != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. "
                               "%s: %d\n",
                               __FILE__, __LINE__);
                    }
                }
                *offset += size;
                break;
            case TRG_DOUBLE_DATA:
            default:
                memcpy(&data->values[i][j].d, block->block_contents+*offset,
                       size);
                if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                {
                    if(trg_swap_byte_order_64(trg_data,
                       (uint64_t *) &data->values[i][j].d) != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. "
                               "%s: %d\n",
                               __FILE__, __LINE__);
                    }
                }
                *offset += size;
            }
        }
    }
    return(TRG_SUCCESS);
}

static trg_function_status trg_write_data_block(struct trg_trajectory *trg_data,
                                                struct trg_gen_block *block,
                                                const int block_index,
                                                const write_mode mode)
{
    uint64_t n_frames;
    int i, j, offset = 0, size, len;
    char temp, *temp_name;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    
    struct trg_data *data;
    trg_block_type block_type_flag;

    if(trg_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        block_type_flag = TRG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TRG_NON_TRAJECTORY_BLOCK;
    }

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(block_type_flag == TRG_TRAJECTORY_BLOCK)
    {
        data = &frame_set->tr_data[block_index];
    }
    else
    {
        data = &trg_data->non_tr_data[block_index];
    }

    switch(data->datatype)
    {
    case TRG_CHAR_DATA:
        size = 1;
        break;
    case TRG_INT_DATA:
        size = sizeof(uint64_t);
        break;
    case TRG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TRG_DOUBLE_DATA:
    default:
        size = sizeof(double);
    }

    /* If just dumping the whole block_contents it is not certain that the 
     * contents are known beforehand (e.g. due to different file versions) */
    if(mode == TRG_COPY_EXISTING)
    {
        if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
        {
            printf("Cannot write header of file %s. %s: %d\n",
                   trg_data->output_file_path, __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }

        if(!block->block_contents)
        {
            printf("No block data to write. %s: %d\n", __FILE__, __LINE__);
            return(TRG_FAILURE);
        }
        if(fwrite(block->block_contents, block->block_contents_size, 1,
                  trg_data->output_file) != 1)
        {
            printf("Could not write all block data. %s: %d\n",
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        return(TRG_SUCCESS);
    }

    len = strlen(data->block_name) + 1;

    if(!block->name || strlen(block->name) < len)
    {
        temp_name = realloc(block->name, len);
        if(!temp_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len+1,
                   __FILE__, __LINE__);
            trg_destroy_block(block);
            return(TRG_CRITICAL);
        }
        block->name = temp_name;
    }
    strncpy(block->name, data->block_name, len);

    /* If writing frame independent data data->n_frames is be 0, but n_frames
       is used for the loop writing the data (and reserving memory) and needs
       to be at least 1 */
    n_frames = max(1, data->n_frames);

    
    block->block_contents_size = sizeof(char) * 2 +
                                 sizeof(data->n_values_per_frame) +
                                 sizeof(data->codec_id);

    if(data->codec_id != TRG_UNCOMPRESSED)
    {
        block->block_contents_size += sizeof(data->compression_multiplier);
    }
    
    if(data->datatype == TRG_CHAR_DATA)
    {
        for(i = n_frames; i--;)
        {
            for(j = data->n_values_per_frame; j--;)
            {
                block->block_contents_size += strlen(data->values[i][j].c) + 1;
            }
        }
    }
    else
    {
        block->block_contents_size += size * n_frames *
        data->n_values_per_frame;
    }

    if(data->n_frames > 0 || data->stride_length != 0)
    {
        block->block_contents_size += sizeof(char);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }
    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }


    memcpy(block->block_contents, &data->datatype, sizeof(char));
    offset += sizeof(char);

    if(data->n_frames > 0 || data->stride_length != 0)
    {
        temp = TRG_FRAME_DEPENDENT;
    }
    else
    {
        temp = 0;
    }
    memcpy(block->block_contents+offset, &temp, sizeof(char));
    offset += sizeof(char);

    if(data->n_frames > 0)
    {
        if(data->stride_length > 1)
        {
            temp = 1;
        }
        else
        {
            temp = 0;
        }
        memcpy(block->block_contents+offset, &temp, sizeof(char));
        offset += sizeof(char);
    }

    memcpy(block->block_contents+offset, &data->n_values_per_frame,
           sizeof(data->n_values_per_frame));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents + offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. "
                   "%s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->n_values_per_frame);

    memcpy(block->block_contents+offset, &data->codec_id,
           sizeof(data->codec_id));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                  (block->block_contents + offset)) !=
           TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. "
                   "%s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(data->codec_id);

    if(data->codec_id != TRG_UNCOMPRESSED)
    {
        memcpy(block->block_contents+offset, &data->compression_multiplier,
               sizeof(data->compression_multiplier));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents + offset)) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. "
                       "%s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->compression_multiplier);
    }

    if(data->n_frames > 0 && data->stride_length > 1)
    {
        memcpy(block->block_contents+offset, &data->first_frame_with_data,
               sizeof(data->first_frame_with_data));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents + offset)) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. "
                       "%s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->first_frame_with_data);

        memcpy(block->block_contents+offset, &data->stride_length,
               sizeof(data->stride_length));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)
                                      (block->block_contents + offset)) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. "
                       "%s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(data->stride_length);
    }

    for(i = 0; i < n_frames; i++)
    {
        for(j = 0; j < data->n_values_per_frame; j++)
        {
            switch(data->datatype)
            {
            case TRG_CHAR_DATA:
                len = strlen(data->values[i][j].c) + 1;
                strncpy(block->block_contents+offset, data->values[i][j].c,
                        len);
                offset += len;
                break;
            case TRG_INT_DATA:
                memcpy(block->block_contents+offset, &data->values[i][j].i,
                       size);
                break;
            case TRG_FLOAT_DATA:
                memcpy(block->block_contents+offset, &data->values[i][j].f,
                       size);
                break;
            case TRG_DOUBLE_DATA:
            default:
                memcpy(block->block_contents+offset, &data->values[i][j].d,
                       size);
            }
            if(data->datatype != TRG_CHAR_DATA)
            {
                if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
                {
                    if(trg_swap_byte_order_64(trg_data,
                        (uint64_t *)(block->block_contents+offset))
                        != TRG_SUCCESS)
                    {
                        printf("Cannot swap byte order to get big endian. "
                               "%s: %d\n",
                                __FILE__, __LINE__);
                    }
                }
                offset += size;
            }
        }
    }

    if(trg_write_block_header(trg_data, block, mode) != TRG_SUCCESS)
    {
        printf("Cannot write header of file %s. %s: %d\n",
               trg_data->output_file_path, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(fwrite(block->block_contents, block->block_contents_size, 1,
              trg_data->output_file) != 1)
    {
        printf("Could not write all block data. %s: %d\n",
               __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    return(TRG_SUCCESS);
}

static trg_function_status trg_read_data_block_contents
                (struct trg_trajectory *trg_data,
                 struct trg_gen_block *block)
{
    uint64_t n_values, codec_id, n_frames, first_frame_with_data;
    uint64_t steps_between_data, block_n_particles, first_particle_number;
    double multiplier;
    char datatype, dependency, sparse_data;
    int offset = 0;
    trg_bool same_hash;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    if(block->block_contents)
    {
        free(block->block_contents);
    }

    block->block_contents = (char *) malloc(block->block_contents_size);
    if(!block->block_contents)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               block->block_contents_size, __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* Read the whole block into block_contents to be able to write it to
     * disk even if it cannot be interpreted. */
    if(fread(block->block_contents, block->block_contents_size, 1,
             trg_data->input_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        trg_destroy_block(block);
        return(TRG_CRITICAL);
    }

    /* FIXME: Does not check if the size of the contents matches the expected
     * size or if the contents can be read. */


    if(trg_verify_hash_match(block, &same_hash) != TRG_SUCCESS)
    {
        printf("Error comparing hashes. %s: %d\n", __FILE__, __LINE__);
        return(TRG_FAILURE);
    }
    if(same_hash != TRUE)
    {
        printf("%s\n", block->hash);
        printf("Data block contents corrupt. Hashes do not match. %s: %d\n",
               __FILE__, __LINE__);
//         return(TRG_FAILURE);
    }

    memcpy(&datatype, block->block_contents+offset,
           sizeof(datatype));
    offset += sizeof(datatype);

    memcpy(&dependency, block->block_contents+offset,
           sizeof(dependency));
    offset += sizeof(dependency);

//     memcpy(&var_n_particles, block->block_contents+offset,
//            sizeof(var_n_particles));
//     offset += sizeof(var_n_particles);

    if(dependency & TRG_FRAME_DEPENDENT)
    {
        memcpy(&sparse_data, block->block_contents+offset,
               sizeof(sparse_data));
        offset += sizeof(sparse_data);
    }

//     memcpy(&var_n_values, block->block_contents+offset,
//            sizeof(var_n_values));
//     offset += sizeof(var_n_values);
//     
//     if(!var_n_values)
//     {
    memcpy(&n_values, block->block_contents+offset,
        sizeof(n_values));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &n_values) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(n_values);
//     }

    memcpy(&codec_id, block->block_contents+offset,
        sizeof(codec_id));
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &codec_id) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                   __FILE__, __LINE__);
        }
    }
    offset += sizeof(codec_id);

    if(codec_id != TRG_UNCOMPRESSED)
    {
        memcpy(&multiplier, block->block_contents+offset,
            sizeof(multiplier));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, (uint64_t *)&multiplier) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(multiplier);
    }
    else
    {
        multiplier = 1;
    }

    if(dependency & TRG_FRAME_DEPENDENT)
    {
        if(sparse_data)
        {
            memcpy(&first_frame_with_data, block->block_contents+offset,
                sizeof(first_frame_with_data));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, &first_frame_with_data) !=
                   TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(first_frame_with_data);
            
            memcpy(&steps_between_data, block->block_contents+offset,
                sizeof(steps_between_data));
            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data, &steps_between_data) !=
                   TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                           __FILE__, __LINE__);
                }
            }
            offset += sizeof(steps_between_data);
        }
        else
        {
            first_frame_with_data = 0;
            steps_between_data = 0;
        }
        n_frames = trg_data->current_trajectory_frame_set.n_frames;
    }
    else
    {
        first_frame_with_data = 0;
        steps_between_data = 0;
        n_frames = 0;
    }
    
    if (dependency & TRG_PARTICLE_DEPENDENT)
    {
        memcpy(&first_particle_number, block->block_contents+offset,
            sizeof(first_particle_number));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &first_particle_number) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(first_particle_number);

        memcpy(&block_n_particles, block->block_contents+offset,
            sizeof(block_n_particles));
        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &block_n_particles) !=
               TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                       __FILE__, __LINE__);
            }
        }
        offset += sizeof(block_n_particles);
    }

    if (dependency & TRG_PARTICLE_DEPENDENT)
    {
        return(trg_read_particle_data(trg_data, block,
                                      &offset, datatype,
                                      first_particle_number,
                                      block_n_particles,
                                      first_frame_with_data,
                                      steps_between_data,
                                      n_frames, n_values,
                                      codec_id, multiplier));
    }
    else
    {
        return(trg_read_data(trg_data, block,
                             &offset, datatype,
                             first_frame_with_data,
                             steps_between_data,
                             n_frames, n_values,
                             codec_id, multiplier));
    }
}

static trg_function_status trg_update_md5_hash(struct trg_trajectory *trg_data,
                                               struct trg_gen_block *block,
                                               uint64_t header_start_pos,
                                               uint64_t contents_start_pos)
{
    if(block->block_contents)
    {
        free(block->block_contents);
    }
    
    block->block_contents = malloc(block->block_contents_size);
    fseek(trg_data->output_file, contents_start_pos, SEEK_SET);
    if(fread(block->block_contents, block->block_contents_size, 1,
            trg_data->output_file) == 0)
    {
        printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    trg_generate_block_hash(block);

    fseek(trg_data->output_file, header_start_pos + 3 * sizeof(uint64_t),
          SEEK_SET);
    fwrite(block->hash, TRG_HASH_LEN, 1, trg_data->output_file);

    return(TRG_SUCCESS);
}


static trg_function_status trg_update_header_pointers
                (struct trg_trajectory *trg_data)
{
    struct trg_gen_block block;
    FILE *temp = trg_data->input_file;
    uint64_t pos, contents_start_pos;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    trg_data->input_file = trg_data->output_file;
    
    trg_init_block(&block);

    trg_data->output_file_pos = ftell(trg_data->output_file);
    fseek(trg_data->output_file, 0, SEEK_SET);

    if(trg_read_block_header(trg_data, &block) != TRG_SUCCESS)
    {
        printf("Cannot read general info header. %s: %d\n",
               __FILE__, __LINE__);
        trg_destroy_block(&block);
        trg_data->input_file = temp;
        return(TRG_CRITICAL);
    }

    contents_start_pos = ftell(trg_data->output_file);

    fseek(trg_data->output_file, block.block_contents_size - 3 *
          sizeof(uint64_t), SEEK_CUR);

    trg_data->input_file = temp;

//     printf("Updating header\n");
//     printf("%ld: First frame set %ld\n", ftell(trg_data->output_file),
//            trg_data->first_trajectory_frame_set_output_file_pos);

    pos = trg_data->first_trajectory_frame_set_output_file_pos;
    
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &pos) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }
        
    if(fwrite(&pos, sizeof(uint64_t), 1, trg_data->output_file) != 1)
    {
        trg_destroy_block(&block);
        return(TRG_CRITICAL);
    }
    
//     printf("%ld: Last frame set %ld\n", ftell(trg_data->output_file),
//            trg_data->last_trajectory_frame_set_output_file_pos);

    pos = trg_data->last_trajectory_frame_set_output_file_pos;
    
    if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
    {
        if(trg_swap_byte_order_64(trg_data, &pos) != TRG_SUCCESS)
        {
            printf("Cannot swap byte order to get big endian. %s: %d\n",
                    __FILE__, __LINE__);
        }
    }

    if(fwrite(&pos,
        sizeof(uint64_t), 1, trg_data->output_file) != 1)
    {
        trg_destroy_block(&block);
        return(TRG_CRITICAL);
    }

    trg_update_md5_hash(trg_data, &block, 0, contents_start_pos);
    
    fseek(trg_data->output_file, trg_data->output_file_pos, SEEK_SET);

    trg_destroy_block(&block);

    return(TRG_SUCCESS);
}

static trg_function_status trg_update_frame_set_pointers
                (struct trg_trajectory *trg_data)
{
    struct trg_gen_block block;
    struct trg_trajectory_frame_set *frame_set;
    FILE *temp = trg_data->input_file;
    uint64_t pos, header_start_pos, contents_start_pos;

    if(trg_init_output_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        printf("Cannot initialise destination file. %s: %d\n",
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    trg_init_block(&block);
    trg_data->output_file_pos = ftell(trg_data->output_file);
    
    trg_data->input_file = trg_data->output_file;

    frame_set = &trg_data->current_trajectory_frame_set;
    if(frame_set->prev_frame_set_file_pos != -1ULL &&
       frame_set->prev_frame_set_file_pos != 0)
    {
        fseek(trg_data->output_file, frame_set->prev_frame_set_file_pos,
              SEEK_SET);

        header_start_pos = frame_set->prev_frame_set_file_pos;
        
        if(trg_read_block_header(trg_data, &block) != TRG_SUCCESS)
        {
            printf("Cannot read frame header. %s: %d\n",
                __FILE__, __LINE__);
            trg_destroy_block(&block);
            trg_data->input_file = temp;
            return(TRG_CRITICAL);
        }

        contents_start_pos = ftell(trg_data->output_file);

        fseek(trg_data->output_file, block.block_contents_size - 4 *
            sizeof(uint64_t), SEEK_CUR);

//         printf("Updating frame set\n");
//         printf("%ld: Next frame set %ld\n", ftell(trg_data->output_file),
//             trg_data->current_trajectory_frame_set_output_file_pos);

        pos = trg_data->current_trajectory_frame_set_output_file_pos;

        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &pos) != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fwrite(&pos, sizeof(uint64_t), 1, trg_data->output_file) != 1)
        {
            trg_destroy_block(&block);
            trg_data->input_file = temp;
            return(TRG_CRITICAL);
        }
        
        trg_update_md5_hash(trg_data, &block, header_start_pos,
                            contents_start_pos);
        
        fseek(trg_data->output_file, trg_data->output_file_pos, SEEK_SET);
    }

    if(frame_set->long_stride_prev_frame_set_file_pos != -1ULL &&
       frame_set->long_stride_prev_frame_set_file_pos != 0)
    {
        fseek(trg_data->output_file,
              frame_set->long_stride_prev_frame_set_file_pos,
              SEEK_SET);

        if(trg_read_block_header(trg_data, &block) != TRG_SUCCESS)
        {
            printf("Cannot read frame header. %s: %d\n",
                __FILE__, __LINE__);
            trg_destroy_block(&block);
            trg_data->input_file = temp;
            return(TRG_CRITICAL);
        }

        contents_start_pos = ftell(trg_data->output_file);
        
        fseek(trg_data->output_file, block.block_contents_size - 2 *
            sizeof(uint64_t), SEEK_CUR);

        pos = trg_data->current_trajectory_frame_set_output_file_pos;

        if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
        {
            if(trg_swap_byte_order_64(trg_data, &pos) != TRG_SUCCESS)
            {
                printf("Cannot swap byte order to get big endian. %s: %d\n",
                        __FILE__, __LINE__);
            }
        }

        if(fwrite(&pos, sizeof(uint64_t), 1, trg_data->output_file) != 1)
        {
            trg_destroy_block(&block);
            trg_data->input_file = temp;
            return(TRG_CRITICAL);
        }
        
        trg_update_md5_hash(trg_data, &block,
                            frame_set->long_stride_prev_frame_set_file_pos,
                            contents_start_pos);
    }

    fseek(trg_data->output_file, trg_data->output_file_pos, SEEK_SET);
    
    trg_data->input_file = temp;
    trg_destroy_block(&block);

    return(TRG_SUCCESS);
}


trg_function_status trg_set_block_name(struct trg_trajectory *trg_data,
                                       struct trg_gen_block *block,
                                       const char *new_name)
{
    int len;
    
    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);
    
    if(block->name && strlen(block->name) < len)
    {
        free(block->name);
        block->name = 0;
    }
    if(!block->name)
    {
        block->name = (char *) malloc(len);
        if(!block->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }    
    
    strncpy(block->name, new_name, len);
    
    return(TRG_SUCCESS);
}

trg_function_status trg_init_block(struct trg_gen_block *block)
{
//     printf("In trg_init_block\n");

    block->id = -1;
/*    block->hash_type = TRG_NO_HASH;
    block->hash_name = 0;*/
    memcpy(block->hash, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0", TRG_HASH_LEN);
    block->name = 0;
    block->block_version = TRG_VERSION;
    block->header_contents = 0;
    block->header_contents_size = 0;
    block->block_contents = 0;
    block->block_contents_size = 0;

    return(TRG_SUCCESS);
}


trg_function_status trg_destroy_block(struct trg_gen_block *block)
{
//     printf("Destroying block\n");
/*    if(block->hash_name)
    {
        free(block->hash_name);
        block->hash_name = 0;
    }*/
/*    if(block->hash)
    {
        free(block->hash);
        block->hash = 0;
    }*/
    if(block->name)
    {
        free(block->name);
        block->name = 0;
    }
    if(block->header_contents)
    {
        free(block->header_contents);
        block->header_contents = 0;
    }
    if(block->block_contents)
    {
        free(block->block_contents);
        block->block_contents = 0;
    }

    return(TRG_SUCCESS);
}

trg_function_status trg_set_atom_name(struct trg_trajectory *trg_data,
                                      struct trg_atom *atom,
                                      const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(atom->name && strlen(atom->name) < len)
    {
        free(atom->name);
        atom->name = 0;
    }
    if(!atom->name)
    {
        atom->name = (char *) malloc(len);
        if(!atom->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(atom->name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_set_atom_type(struct trg_trajectory *trg_data,
                                      struct trg_atom *atom,
                                      const char *new_type)
{
    int len;

    len = min(strlen(new_type) + 1, TRG_MAX_STR_LEN);

    if(atom->atom_type && strlen(atom->atom_type) < len)
    {
        free(atom->atom_type);
        atom->atom_type = 0;
    }
    if(!atom->atom_type)
    {
        atom->atom_type = (char *) malloc(len);
        if(!atom->atom_type)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(atom->atom_type, new_type, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_init_atom(struct trg_atom *atom)
{
    atom->name = 0;
    atom->atom_type = 0;
    
    return(TRG_SUCCESS);
}

trg_function_status trg_destroy_atom(struct trg_atom *atom)
{
    if(atom->name)
    {
        free(atom->name);
        atom->name = 0;
    }
    if(atom->atom_type)
    {
        free(atom->atom_type);
        atom->atom_type = 0;
    }
    
    return(TRG_SUCCESS);
}

trg_function_status trg_add_molecule(struct trg_trajectory *trg_data,
                                     const char *name,
                                     struct trg_molecule **molecule)
{
    struct trg_molecule *new_molecules;
    uint64_t *new_molecule_cnt_list;
    int id, i;
    trg_bool found_id = TRUE;
    
    new_molecules = (struct trg_molecule *)realloc(trg_data->molecules,
                                                   sizeof(struct trg_molecule) *
                                                   (trg_data->n_molecules + 1));

    if(!new_molecules)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct trg_molecule) * (trg_data->n_molecules + 1),
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    new_molecule_cnt_list = (uint64_t *) realloc(trg_data->molecule_cnt_list,
                                                 sizeof(uint64_t) *
                                                 (trg_data->n_molecules + 1));
    
    if(!new_molecules)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(uint64_t) * (trg_data->n_molecules + 1),
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    trg_data->molecules = new_molecules;
    trg_data->molecule_cnt_list = new_molecule_cnt_list;

    *molecule = &new_molecules[trg_data->n_molecules];

    trg_init_molecule(*molecule);
    trg_set_molecule_name(trg_data, *molecule, name);

    /* FIXME: Should this be a function argument instead? */
    trg_data->molecule_cnt_list[trg_data->n_molecules] = 0;

    /* Find an unused ID */
    id = 0;
    while(found_id)
    {
        found_id = FALSE;
        for(i = trg_data->n_molecules; i--;)
        {
            if(trg_data->molecules[i].id == id)
            {
                found_id = TRUE;
                i = 0;
            }
        }
        if(found_id)
        {
            id++;
        }
    }
    
    (*molecule)->id = id;

    trg_data->n_molecules++;

    return(TRG_SUCCESS);
}

trg_function_status trg_set_molecule_name(struct trg_trajectory *trg_data,
                                          struct trg_molecule *molecule,
                                          const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(molecule->name && strlen(molecule->name) < len)
    {
        free(molecule->name);
        molecule->name = 0;
    }
    if(!molecule->name)
    {
        molecule->name = (char *) malloc(len);
        if(!molecule->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(molecule->name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_get_molecule_cnt(struct trg_trajectory *trg_data,
                                         struct trg_molecule *molecule,
                                         uint64_t *cnt)
{
    int i, index = -1;

    for(i = trg_data->n_molecules; i--;)
    {
        if(&trg_data->molecules[i] == molecule)
        {
            index = i;
            i = 0;
        }
    }
    if(index == -1)
    {
        return(TRG_FAILURE);
    }
    *cnt = trg_data->molecule_cnt_list[index];

    return(TRG_SUCCESS);
}

trg_function_status trg_set_molecule_cnt(struct trg_trajectory *trg_data,
                                         struct trg_molecule *molecule,
                                         const uint64_t cnt)
{
    int i, index = -1, old_cnt;

    for(i = trg_data->n_molecules; i--;)
    {
        if(&trg_data->molecules[i] == molecule)
        {
            index = i;
            i = 0;
        }
    }
    if(index == -1)
    {
        return(TRG_FAILURE);
    }
    old_cnt = trg_data->molecule_cnt_list[index];
    trg_data->molecule_cnt_list[index] = cnt;

    trg_data->n_particles += (cnt-old_cnt) *
                             trg_data->molecules[index].n_atoms;

    return(TRG_SUCCESS);
}

trg_function_status trg_add_chain_to_molecule(struct trg_trajectory *trg_data,
                                              struct trg_molecule *molecule,
                                              const char *name,
                                              struct trg_chain **chain)
{
    struct trg_chain *new_chains;

    new_chains = (struct trg_chain *) realloc(molecule->chains,
                                              sizeof(struct trg_chain) *
                                              (molecule->n_chains + 1));

    if(!new_chains)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct trg_chain) * (molecule->n_chains + 1),
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    molecule->chains = new_chains;

    *chain = &new_chains[molecule->n_chains];
    (*chain)->name = 0;

    trg_set_chain_name(trg_data, *chain, name);

    (*chain)->molecule = molecule;
    (*chain)->id = molecule->n_chains;
    (*chain)->n_residues = 0;

    molecule->n_chains++;

    return(TRG_SUCCESS);
}

trg_function_status trg_set_chain_name(struct trg_trajectory *trg_data,
                                       struct trg_chain *chain,
                                       const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(chain->name && strlen(chain->name) < len)
    {
        free(chain->name);
        chain->name = 0;
    }
    if(!chain->name)
    {
        chain->name = (char *) malloc(len);
        if(!chain->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(chain->name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_add_residue_to_chain(struct trg_trajectory *trg_data,
                                             struct trg_chain *chain,
                                             const char *name,
                                             struct trg_residue **residue)
{
    int curr_index;
    struct trg_residue *new_residues, *temp_residue, *last_residue;
    struct trg_molecule *molecule = chain->molecule;

    if(chain->n_residues)
    {
        curr_index = (chain->residues - molecule->residues) /
                     sizeof(struct trg_residue);
    }
    else
    {
        curr_index = -1;
    }

    new_residues = (struct trg_residue *) realloc(molecule->residues,
                                                  sizeof(struct trg_residue) *
                                                  (molecule->n_residues + 1));

    if(!new_residues)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct trg_residue) * (molecule->n_residues + 1),
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    molecule->residues = new_residues;

    if(curr_index != -1)
    {
        chain->residues = new_residues + curr_index * sizeof(struct trg_residue);
        if(molecule->n_residues)
        {
            last_residue = &new_residues[molecule->n_atoms-1];

            temp_residue = chain->residues + (chain->n_residues - 1);
            /* Make space in list of residues to add the new residues together with the other
            * residues of this chain */
            if(temp_residue != last_residue)
            {
                temp_residue++;
                memmove(temp_residue + 1, temp_residue,
                        last_residue - temp_residue);
            }
        }
    }
    else
    {
        curr_index = molecule->n_residues;
    }

    *residue = &new_residues[curr_index + chain->n_residues];

    if(!chain->n_residues)
    {
        chain->residues = *residue;
    }
    
    (*residue)->name = 0;
    trg_set_residue_name(trg_data, *residue, name);

    (*residue)->chain = chain;
    (*residue)->id = chain->n_residues;
    (*residue)->n_atoms = 0;
    (*residue)->atoms = 0;

    chain->n_residues++;
    molecule->n_residues++;

    return(TRG_SUCCESS);
}

trg_function_status trg_set_residue_name(struct trg_trajectory *trg_data,
                                         struct trg_residue *residue,
                                         const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(residue->name && strlen(residue->name) < len)
    {
        free(residue->name);
        residue->name = 0;
    }
    if(!residue->name)
    {
        residue->name = (char *) malloc(len);
        if(!residue->name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(residue->name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_add_atom_to_residue(struct trg_trajectory *trg_data,
                                            struct trg_residue *residue,
                                            const char *atom_name,
                                            const char *atom_type,
                                            struct trg_atom **atom)
{
    int curr_index;
    struct trg_atom *new_atoms, *temp_atom, *last_atom;
    struct trg_molecule *molecule = residue->chain->molecule;

    if(residue->n_atoms)
    {
        curr_index = (residue->atoms - molecule->atoms) /
                     sizeof(struct trg_atom);
    }
    else
    {
        curr_index = -1;
    }

    new_atoms = (struct trg_atom *) realloc(molecule->atoms,
                                            sizeof(struct trg_atom) *
                                            (molecule->n_atoms + 1));

    if(!new_atoms)
    {
        printf("Cannot allocate memory (%"PRId64" bytes). %s: %d\n",
               sizeof(struct trg_atom) * (molecule->n_atoms + 1),
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }

    molecule->atoms = new_atoms;

    if(curr_index != -1)
    {
        residue->atoms = new_atoms + curr_index * sizeof(struct trg_atom);
        if(molecule->n_atoms)
        {
            last_atom = &new_atoms[molecule->n_atoms-1];

            temp_atom = residue->atoms + (residue->n_atoms - 1);
            /* Make space in list of atoms to add the new atoms together with the other
            * atoms of this residue */
            if(temp_atom != last_atom)
            {
                temp_atom++;
                memmove(temp_atom + 1, temp_atom,
                        last_atom - temp_atom);
            }
        }
    }
    else
    {
        curr_index = molecule->n_atoms;
    }

    *atom = &new_atoms[curr_index + residue->n_atoms];

    if(!residue->n_atoms)
    {
        residue->atoms = *atom;
    }

    trg_init_atom(*atom);
    trg_set_atom_name(trg_data, *atom, atom_name);
    trg_set_atom_type(trg_data, *atom, atom_type);

    (*atom)->residue = residue;
    (*atom)->id = molecule->n_atoms;

    residue->n_atoms++;
    molecule->n_atoms++;
    
    return(TRG_SUCCESS);
}


trg_function_status trg_init_molecule(struct trg_molecule *molecule)
{
    molecule->name = 0;
    molecule->n_chains = 0;
    molecule->chains = 0;
    molecule->n_residues = 0;
    molecule->residues = 0;
    molecule->n_atoms = 0;
    molecule->atoms = 0;
    molecule->n_bonds = 0;
    molecule->bonds = 0;
    
    return(TRG_SUCCESS);
}

trg_function_status trg_destroy_molecule(struct trg_molecule *molecule)
{
    int i;

    if(molecule->name)
    {
        free(molecule->name);
        molecule->name = 0;
    }
    
    if(molecule->chains)
    {
        for(i = molecule->n_chains; i--;)
        {
            if(molecule->chains[i].name)
            {
                free(molecule->chains[i].name);
                molecule->chains[i].name = 0;
            }
        }
        free(molecule->chains);
        molecule->chains = 0;
    }
    molecule->n_chains = 0;

    if(molecule->residues)
    {
        for(i = molecule->n_residues; i--;)
        {
            if(molecule->residues[i].name)
            {
                free(molecule->residues[i].name);
                molecule->residues[i].name = 0;
            }
        }
        free(molecule->residues);
        molecule->residues = 0;
    }
    molecule->n_residues = 0;

    if(molecule->atoms)
    {
        for(i = molecule->n_atoms; i--;)
        {
            trg_destroy_atom(&molecule->atoms[i]);
        }
        free(molecule->atoms);
        molecule->atoms = 0;
    }
    molecule->n_atoms = 0;
    
    if(molecule->bonds)
    {
        free(molecule->bonds);
        molecule->bonds = 0;
    }
    molecule->n_bonds = 0;

    return(TRG_SUCCESS);
}

trg_function_status trg_init_trajectory(struct trg_trajectory *trg_data)
{
    time_t seconds;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    
    trg_data->input_file_path = 0;
    trg_data->input_file = 0;
    trg_data->input_file_pos = 0;
    trg_data->input_file_len = 0;
    trg_data->output_file_path = 0;
    trg_data->output_file = 0;
    trg_data->output_file_pos = 0;

    trg_data->program_name = 0;
    trg_data->forcefield_name = 0;
    
    /* FIXME: No unistd.h on Windows!! */
//     trg_data->user_name = (char *) malloc(LOGIN_NAME_MAX);
//     if(getlogin_r(trg_data->user_name, LOGIN_NAME_MAX) != 0)
//     {
//         printf("Cannot get user name. %s: %d\n", __FILE__, __LINE__);
//         free(trg_data->user_name);
//         trg_data->user_name = 0;
//     }
    trg_data->user_name = 0;
    
    seconds = time(0);
    if ( seconds == -1)
    {
        printf("Cannot get time. %s: %d\n", __FILE__, __LINE__);
    }
    else
    {
        trg_data->time = seconds;
    }
    
    /* FIXME: No unistd.h on Windows!! */
    /* FIXME: Append operating system to computer_name */
//     trg_data->computer_name = (char *) malloc(HOST_NAME_MAX);
//     if(gethostname(trg_data->computer_name, HOST_NAME_MAX) != 0)
//     {
//         printf("Cannot get computer name. %s: %d\n", __FILE__, __LINE__);
//         free(trg_data->computer_name);
//         trg_data->computer_name = 0;
//     }
    trg_data->computer_name = 0;

    trg_data->pgp_signature = 0;
    trg_data->var_num_atoms_flag = TRG_CONSTANT_N_ATOMS;
    trg_data->first_trajectory_frame_set_input_file_pos = -1;
    trg_data->last_trajectory_frame_set_input_file_pos = -1;
    trg_data->current_trajectory_frame_set_input_file_pos = -1;
    trg_data->first_trajectory_frame_set_output_file_pos = -1;
    trg_data->last_trajectory_frame_set_output_file_pos = -1;
    trg_data->current_trajectory_frame_set_output_file_pos = -1;
    trg_data->frame_set_n_frames = 100;
    trg_data->n_trajectory_frame_sets = 0;
    trg_data->n_trajectory_blocks = 0;
    trg_data->stride_length = 100;

    trg_data->n_particle_data_blocks = 0;
    trg_data->n_data_blocks = 0;

    trg_data->non_tr_particle_data = 0;
    trg_data->non_tr_data = 0;
    
    frame_set->contents.n_blocks = 0;
    frame_set->contents.block_names = 0;
    frame_set->n_mapping_blocks = 0;
    frame_set->mappings = 0;
    frame_set->molecule_cnt_list = 0;
    
    frame_set->n_particle_data_blocks = 0;
    frame_set->n_data_blocks = 0;

    frame_set->tr_particle_data = 0;
    frame_set->tr_data = 0;
    
    trg_data->n_molecules = 0;
    trg_data->molecules = 0;
    trg_data->molecule_cnt_list = 0;
    trg_data->n_particles = 0;

    trg_data->n_id_name_pairs = 0;
    trg_data->id_name_pairs = 0;

    /* Check the endianness of the computer */
    static uint32_t endianness_32 = 0x01234567;
    /* 0x01234567 */
    if ( *(const uint8_t*)&endianness_32 == 0x01 )
    {
        trg_data->endianness_32 = TRG_BIG_ENDIAN_32;
    }

    /* 0x67452301 */
    else if( *(const uint8_t*)&endianness_32 == 0x67 )
    {
        trg_data->endianness_32 = TRG_LITTLE_ENDIAN_32;

    }
    
    /* 0x45670123 */
    else if ( *(const uint8_t*)&endianness_32 == 0x45 )
    {
        trg_data->endianness_32 = TRG_BYTE_PAIR_SWAP_32;
    }

    static uint64_t endianness_64 = 0x0123456789ABCDEF;
    /* 0x0123456789ABCDEF */
    if ( *(const uint8_t*)&endianness_64 == 0x01 )
    {
        trg_data->endianness_64 = TRG_BIG_ENDIAN_64;
    }
    
    /* 0xEFCDAB8967452301 */
    else if ( *(const uint8_t*)&endianness_64 == 0xEF )
    {
        trg_data->endianness_64 = TRG_LITTLE_ENDIAN_64;
    }

    /* 0x89ABCDEF01234567 */
    else if ( *(const uint8_t*)&endianness_64 == 0x89 )
    {
        trg_data->endianness_64 = TRG_QUAD_SWAP_64;
    }

    /* 0x45670123CDEF89AB */
    else if ( *(const uint8_t*)&endianness_64 == 0x45 )
    {
        trg_data->endianness_64 = TRG_BYTE_PAIR_SWAP_64;
    }
    
    /* 0x23016745AB89EFCD */
    else if ( *(const uint8_t*)&endianness_64 == 0x23 )
    {
        trg_data->endianness_64 = TRG_BYTE_SWAP_64;
    }

    
    
    trg_init_block(&trg_data->non_trajectory_blocks[0]);
    trg_data->non_trajectory_blocks[0].id = TRG_GENERAL_INFO;
    trg_set_block_name(trg_data, &trg_data->non_trajectory_blocks[0],
                       "GENERAL INFO");

    trg_data->current_trajectory_frame_set.next_frame_set_file_pos = -1;
    trg_data->current_trajectory_frame_set.prev_frame_set_file_pos = -1;
    
    /* The Endianness and String Length block and the Trajectory Info block
     * are present. */
    trg_data->n_non_trajectory_blocks = 1;
    
    return(TRG_SUCCESS);
}

trg_function_status trg_destroy_trajectory(struct trg_trajectory *trg_data)
{
    uint64_t n_frames, n_particles;
    int i, j, k, l;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;

    struct trg_particle_mapping *mapping;
    
    if(trg_data->input_file_path)
    {
        free(trg_data->input_file_path);
        trg_data->input_file_path = 0;
    }
    
    if(trg_data->input_file)
    {
        fclose(trg_data->input_file);
        trg_data->input_file = 0;
    }
    
    if(trg_data->output_file_path)
    {
        free(trg_data->output_file_path);
        trg_data->output_file_path = 0;
    }
    
    if(trg_data->output_file)
    {
        fclose(trg_data->output_file);
        trg_data->output_file = 0;
    }
    
    if(trg_data->program_name)
    {
        free(trg_data->program_name);
        trg_data->program_name = 0;
    }

    if(trg_data->forcefield_name)
    {
        free(trg_data->forcefield_name);
        trg_data->forcefield_name = 0;
    }
    
    if(trg_data->user_name)
    {
        free(trg_data->user_name);
        trg_data->user_name = 0;
    }
    
    if(trg_data->computer_name)
    {
        free(trg_data->computer_name);
        trg_data->computer_name = 0;
    }

    if(trg_data->pgp_signature)
    {
        free(trg_data->pgp_signature);
        trg_data->pgp_signature = 0;
    }

    if(frame_set->contents.block_names)
    {
        for(i = frame_set->contents.n_blocks; i--;)
        {
            if(frame_set->contents.block_names[i])
            {
                free(frame_set->contents.block_names[i]);
                frame_set->contents.block_names[i] = 0;
            }
        }
        free(frame_set->contents.block_names);
        frame_set->contents.block_names = 0;
    }

    for(i = frame_set->n_mapping_blocks; i--;)
    {
        mapping = &frame_set->mappings[i];
        if(mapping->real_particle_numbers)
        {
            free(mapping->real_particle_numbers);
            mapping->real_particle_numbers = 0;
        }
    }
    free(frame_set->mappings);
    frame_set->n_mapping_blocks = 0;

    if(frame_set->molecule_cnt_list)
    {
        free(frame_set->molecule_cnt_list);
        frame_set->molecule_cnt_list = 0;
    }
    
    for(i=trg_data->n_non_trajectory_blocks; i--;)
    {
        trg_destroy_block(&trg_data->non_trajectory_blocks[i]);
    }
    trg_data->n_trajectory_blocks = 0;

    if(trg_data->var_num_atoms_flag)
    {
        n_particles = trg_data->current_trajectory_frame_set.n_particles;
    }
    else
    {
        n_particles = trg_data->n_particles;
    }

    if(trg_data->non_tr_particle_data)
    {
        for(i = trg_data->n_particle_data_blocks; i--; )
        {
            if(trg_data->non_tr_particle_data[i].values)
            {
                /* Only one frame for non-trajectory data */
                j = 0;
                if(trg_data->non_tr_particle_data[i].values[j])
                {
                    for(k = n_particles; k--;)
                    {
                        if(trg_data->non_tr_particle_data[i].values[j][k])
                        {
                            if(trg_data->non_tr_particle_data[i].datatype ==
                                TRG_CHAR_DATA)
                            {
                                for(l = trg_data->non_tr_particle_data[i].
                                    n_values_per_frame;
                                    l--;)
                                {
                                    if(trg_data->non_tr_particle_data[i].
                                        values[j][k][l].c)
                                    {
                                        free(trg_data->non_tr_particle_data[i].
                                        values[j][k][l].c);
                                        trg_data->non_tr_particle_data[i].
                                        values[j][k][l].c = 0;
                                    }
                                }
                            }
                            free(trg_data->non_tr_particle_data[i].
                            values[j][k]);
                            trg_data->non_tr_particle_data[i].
                            values[j][k] = 0;
                        }
                    }
                    free(trg_data->non_tr_particle_data[i].values[j]);
                    trg_data->non_tr_particle_data[i].values[j] = 0;
                }
                free(trg_data->non_tr_particle_data[i].values);
                trg_data->non_tr_particle_data[i].values = 0;
            }
            if(trg_data->non_tr_particle_data[i].block_name)
            {
                free(trg_data->non_tr_particle_data[i].block_name);
                trg_data->non_tr_particle_data[i].block_name = 0;
            }
        }
        free(trg_data->non_tr_particle_data);
        trg_data->non_tr_particle_data = 0;
    }
    
    if(trg_data->non_tr_data)
    {
        for(i = trg_data->n_data_blocks; i--;)
        {
            if(trg_data->non_tr_data[i].values)
            {
                /* Only one frame for non-trajectory data */
                if(trg_data->non_tr_data[i].values[0])
                {
                    if(trg_data->non_tr_data[i].datatype ==
                        TRG_CHAR_DATA)
                    {
                        for(k = trg_data->non_tr_data[i].n_values_per_frame;
                            k--;)
                        {
                            if(trg_data->non_tr_data[i].values[0][k].c)
                            {
                                free(trg_data->non_tr_data[i].values[0][k].c);
                                trg_data->non_tr_data[i].values[0][k].c = 0;
                            }
                        }
                    }
                    free(trg_data->non_tr_data[i].values[0]);
                    trg_data->non_tr_data[i].values[0] = 0;
                }
                free(trg_data->non_tr_data[i].values);
                trg_data->non_tr_data[i].values = 0;
            }
            if(trg_data->non_tr_data[i].block_name)
            {
                free(trg_data->non_tr_data[i].block_name);
                trg_data->non_tr_data[i].block_name = 0;
            }
        }
        free(trg_data->non_tr_data);
        trg_data->non_tr_data = 0;
    }
    
    trg_data->n_particle_data_blocks = 0;
    trg_data->n_data_blocks = 0;
    
    if(frame_set->tr_particle_data)
    {
        for(i = frame_set->n_particle_data_blocks; i--; )
        {
            if(frame_set->tr_particle_data[i].values)
            {
                n_frames = max(1, frame_set->tr_particle_data[i].n_frames);
                for(j = n_frames; j--;)
                {
                    if(frame_set->tr_particle_data[i].values[j])
                    {
                        for(k = n_particles; k--;)
                        {
                            if(frame_set->tr_particle_data[i].values[j][k])
                            {
                                if(frame_set->tr_particle_data[i].datatype ==
                                    TRG_CHAR_DATA)
                                {
                                    for(l = frame_set->tr_particle_data[i].
                                        n_values_per_frame;
                                        l--;)
                                    {
                                        if(frame_set->tr_particle_data[i].
                                            values[j][k][l].c)
                                        {
                                            free(frame_set->tr_particle_data[i].
                                            values[j][k][l].c);
                                            
                                            frame_set->tr_particle_data[i].
                                            values[j][k][l].c = 0;
                                        }
                                    }
                                }
                                free(frame_set->tr_particle_data[i].
                                     values[j][k]);
                                frame_set->tr_particle_data[i].values[j][k] = 0;
                            }
                        }
                        free(frame_set->tr_particle_data[i].values[j]);
                        frame_set->tr_particle_data[i].values[j] = 0;
                    }
                }
                free(frame_set->tr_particle_data[i].values);
                frame_set->tr_particle_data[i].values = 0;
            }
            if(frame_set->tr_particle_data[i].block_name)
            {
                free(frame_set->tr_particle_data[i].block_name);
                frame_set->tr_particle_data[i].block_name = 0;
            }
        }
        free(frame_set->tr_particle_data);
        frame_set->tr_particle_data = 0;
    }
    
    if(frame_set->tr_data)
    {
        for(i = frame_set->n_data_blocks; i--;)
        {
            if(frame_set->tr_data[i].values)
            {
                n_frames = max(1, frame_set->tr_data[i].n_frames);
                for(j = n_frames; j--;)
                {
                    if(frame_set->tr_data[i].values[j])
                    {
                        if(frame_set->tr_data[i].datatype ==
                            TRG_CHAR_DATA)
                        {
                            for(k = frame_set->tr_data[i].n_values_per_frame;
                                k--;)
                            {
                                if(frame_set->tr_data[i].values[j][k].c)
                                {
                                    free(frame_set->tr_data[i].
                                    values[j][k].c);
                                    
                                    frame_set->tr_data[i].
                                    values[j][k].c = 0;
                                }
                            }
                        }
                        free(frame_set->tr_data[i].values[j]);
                        frame_set->tr_data[i].values[j] = 0;
                    }
                }
                free(frame_set->tr_data[i].values);
                frame_set->tr_data[i].values = 0;
            }
            if(frame_set->tr_data[i].block_name)
            {
                free(frame_set->tr_data[i].block_name);
                frame_set->tr_data[i].block_name = 0;
            }
        }
        free(frame_set->tr_data);
        frame_set->tr_data = 0;
    }
    
    frame_set->n_particle_data_blocks = 0;
    frame_set->n_data_blocks = 0;
    
    if(trg_data->molecules)
    {
        for(i=trg_data->n_molecules; i--;)
        {
            trg_destroy_molecule(&trg_data->molecules[i]);
        }
        free(trg_data->molecules);
        trg_data->molecules = 0;
        trg_data->n_molecules = 0;
    }
    if(trg_data->molecule_cnt_list)
    {
        free(trg_data->molecule_cnt_list);
        trg_data->molecule_cnt_list = 0;
    }
    return(TRG_SUCCESS);
}

trg_function_status trg_set_input_file(struct trg_trajectory *trg_data,
                                       const char *file_name)
{
    int len;
    char *temp;
    
    if(trg_data->input_file_path && strcmp(trg_data->input_file_path,
                                           file_name) == 0)
    {
        return(TRG_SUCCESS);
    }
    
    if(trg_data->input_file)
    {
        fclose(trg_data->input_file);
    }

    len = min(strlen(file_name) + 1, TRG_MAX_STR_LEN);
    temp = realloc(trg_data->input_file_path, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    trg_data->input_file_path = temp;
        
    strncpy(trg_data->input_file_path, file_name, len);

    return(trg_init_input_file(trg_data, FALSE));
}

trg_function_status trg_set_output_file(struct trg_trajectory *trg_data,
                                        const char *file_name)
{
    int len;
    char *temp;

    if(trg_data->output_file_path &&
       strcmp(trg_data->output_file_path, file_name) == 0)
    {
        return(TRG_SUCCESS);
    }

    if(trg_data->output_file)
    {
        fclose(trg_data->output_file);
    }

    len = min(strlen(file_name) + 1, TRG_MAX_STR_LEN);
    temp = realloc(trg_data->output_file_path, len);
    if(!temp)
    {
        printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
               __FILE__, __LINE__);
        return(TRG_CRITICAL);
    }
    trg_data->output_file_path = temp;
        
    strncpy(trg_data->output_file_path, file_name, len);

    return(trg_init_output_file(trg_data, FALSE));
}


trg_function_status trg_set_program_name(struct trg_trajectory *trg_data,
                                         const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(trg_data->program_name && strlen(trg_data->program_name) < len)
    {
        free(trg_data->program_name);
        trg_data->program_name = 0;
    }
    if(!trg_data->program_name)
    {
        trg_data->program_name = (char *) malloc(len);
        if(!trg_data->program_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(trg_data->program_name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_set_forcefield_name(struct trg_trajectory *trg_data,
                                            const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(trg_data->forcefield_name && strlen(trg_data->forcefield_name) < len)
    {
        free(trg_data->forcefield_name);
        trg_data->forcefield_name = 0;
    }
    if(!trg_data->forcefield_name)
    {
        trg_data->forcefield_name = (char *) malloc(len);
        if(!trg_data->forcefield_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(trg_data->forcefield_name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_set_user_name(struct trg_trajectory *trg_data,
                                      const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(trg_data->user_name && strlen(trg_data->user_name) < len)
    {
        free(trg_data->user_name);
        trg_data->user_name = 0;
    }
    if(!trg_data->user_name)
    {
        trg_data->user_name = (char *) malloc(len);
        if(!trg_data->user_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(trg_data->user_name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_set_computer_name(struct trg_trajectory *trg_data,
                                          const char *new_name)
{
    int len;

    len = min(strlen(new_name) + 1, TRG_MAX_STR_LEN);

    if(trg_data->computer_name && strlen(trg_data->computer_name) < len)
    {
        free(trg_data->computer_name);
        trg_data->computer_name = 0;
    }
    if(!trg_data->computer_name)
    {
        trg_data->computer_name = (char *) malloc(len);
        if(!trg_data->computer_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(trg_data->computer_name, new_name, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_set_signature(struct trg_trajectory *trg_data,
                                      const char *signature)
{
    int len;

    len = min(strlen(signature) + 1, TRG_MAX_STR_LEN);

    if(trg_data->pgp_signature && strlen(trg_data->pgp_signature) < len)
    {
        free(trg_data->pgp_signature);
        trg_data->pgp_signature = 0;
    }
    if(!trg_data->pgp_signature)
    {
        trg_data->pgp_signature = (char *) malloc(len);
        if(!trg_data->pgp_signature)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n", len,
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    strncpy(trg_data->pgp_signature, signature, len);

    return(TRG_SUCCESS);
}

trg_function_status trg_read_file_headers(struct trg_trajectory *trg_data,
                                          trg_close_file_flag close_file)
{
    int i, cnt = 0, prev_pos = 0;
    struct trg_gen_block *block = trg_data->non_trajectory_blocks;

    trg_data->input_file_pos = 0;
    
    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        return(TRG_CRITICAL);
    }

    if(!trg_data->input_file_len)
    {
        fseek(trg_data->input_file, 0, SEEK_END);
        trg_data->input_file_len = ftell(trg_data->input_file);
        fseek(trg_data->input_file, 0, SEEK_SET);
    }
    
    for(i = trg_data->n_non_trajectory_blocks; i--;)
    {
        trg_destroy_block(block++);
    }
    trg_data->n_non_trajectory_blocks = 0;

    block = trg_data->non_trajectory_blocks;
    
    trg_init_block(block);
    /* Non trajectory blocks (they come before the trajectory
     * blocks in the file) */
    while (prev_pos < trg_data->input_file_len &&
           trg_read_block_header(trg_data, block) != TRG_CRITICAL &&
           block->id != -1 &&
           block->id != TRG_TRAJECTORY_FRAME_SET &&
           trg_data->n_non_trajectory_blocks < 32)
    {
//         printf("Reading block header %d: %s\n", (int)block->id, block->name);
        if(trg_read_next_block(trg_data, block,
                               TRG_KEEP_FILE_OPEN) == TRG_SUCCESS)
        {
//             printf("Read block %s\n", block->name);
            block++;
            cnt++;
            trg_data->n_non_trajectory_blocks++;
            if(trg_init_block(block) != TRG_SUCCESS)
            {
                return(TRG_CRITICAL);
            }
        }
        else
        {
            trg_destroy_block(block);
        }
        prev_pos = ftell(trg_data->input_file);
    }

    /* Go back if a trajectory block was encountered */
    if(block->id == TRG_TRAJECTORY_FRAME_SET)
    {
        trg_destroy_block(block);
        fseek(trg_data->input_file, prev_pos, SEEK_SET);
    }
    
    if(close_file)
    {
        trg_data->input_file_pos=ftell(trg_data->input_file);
        fclose(trg_data->input_file);
        trg_data->input_file = 0;
    }
    
    return(TRG_SUCCESS);
}

trg_function_status trg_write_file_headers(struct trg_trajectory *trg_data,
                                           trg_close_file_flag close_file)
{
    int i;
    struct trg_gen_block *block, data_block;

    trg_data->output_file_pos = 0;

    if(trg_init_output_file(trg_data, TRUE) != TRG_SUCCESS)
    {
        return(TRG_CRITICAL);
    }


    for(i=0; i<trg_data->n_non_trajectory_blocks; ++i)
    {
        block = &trg_data->non_trajectory_blocks[i];
        if(block->id == TRG_GENERAL_INFO)
        {
            if(trg_write_general_info_block(trg_data, block, TRG_NORMAL_WRITE)
                != TRG_SUCCESS)
            {
                printf("Error writing general info block of file %s. %s: %d\n",
                       trg_data->input_file_path, __FILE__, __LINE__);
                return(TRG_CRITICAL);
            }
            break;
        }
    }
    
    for(i=0; i<trg_data->n_non_trajectory_blocks; ++i)
    {
        block = &trg_data->non_trajectory_blocks[i];
        if(block->id == TRG_MOLECULES)
        {
            if(trg_write_molecules_block(trg_data, block, TRG_NORMAL_WRITE)
                != TRG_SUCCESS)
            {
                printf("Error writing atom names block of file %s. %s: %d\n",
                       trg_data->input_file_path, __FILE__, __LINE__);
                return(TRG_CRITICAL);
            }
            break;
        }
    }

    /* FIXME: Currently writing non-trajectory data blocks here.
     * Should perhaps be moved. */
    trg_init_block(&data_block);
    for(i = 0; i < trg_data->n_data_blocks; i++)
    {
        data_block.id = trg_data->non_tr_data[i].block_id;
        trg_write_data_block(trg_data, &data_block,
                             i, TRG_NORMAL_WRITE);
    }

    for(i = 0; i < trg_data->n_particle_data_blocks; i++)
    {
        data_block.id = trg_data->non_tr_particle_data[i].block_id;
        trg_write_particle_data_block(trg_data, &data_block,
                                      i, 0, TRG_NORMAL_WRITE);
    }

    trg_destroy_block(&data_block);

    if(close_file)
    {
        trg_data->output_file_pos=ftell(trg_data->output_file);
        fclose(trg_data->output_file);
        trg_data->output_file = 0;
    }

    return(TRG_SUCCESS);
}

trg_function_status trg_read_next_block(struct trg_trajectory *trg_data,
                                        struct trg_gen_block *block,
                                        trg_close_file_flag close_file)
{
    switch(block->id)
    {
    case TRG_TRAJECTORY_FRAME_SET:
        return(trg_read_frame_set_block(trg_data, block));
    case TRG_BLOCK_TABLE_OF_CONTENTS:
        return(trg_read_trajectory_toc_block(trg_data, block));
    case TRG_PARTICLE_MAPPING:
        return(trg_read_trajectory_mapping_block(trg_data, block));
    case TRG_GENERAL_INFO:
        return(trg_read_general_info_block(trg_data, block));
    case TRG_MOLECULES:
        return(trg_read_molecules_block(trg_data, block));
    default:
        if(block->id >= TRG_TRAJ_BOX_SHAPE)
        {
            return(trg_read_data_block_contents(trg_data, block));
        }
        else
        {
            /* Skip to the next block */
            fseek(trg_data->input_file, block->block_contents_size, SEEK_CUR);
            return(TRG_FAILURE);
        }
    }

    /* FIXME: Never reached. */
    if(close_file)
    {
        trg_data->input_file_pos=ftell(trg_data->input_file);
        fclose(trg_data->input_file);
        trg_data->input_file = 0;
    }
}

// trg_function_status trg_write_block(struct trg_trajectory *trg_data,
//                                     struct trg_gen_block *block,
//                                     trg_close_file_flag close_file)
// {
//     if(trg_data->output_file)
//     {
//         trg_data->output_file_pos = ftell(trg_data->output_file);
//     }
// 
//     switch(block->id)
//     {
//     case TRG_TRAJECTORY_FRAME_SET:
//         return(trg_write_frame_set_block(trg_data, block, TRG_NORMAL_WRITE));
//         break;
//     case TRG_BLOCK_TABLE_OF_CONTENTS:
//         return(trg_write_trajectory_toc_block(trg_data, block,
//                                               TRG_NORMAL_WRITE));
//         break;
//     case TRG_PARTICLE_MAPPING:
//         return(trg_write_trajectory_mapping_block(trg_data, block,
//                                                   TRG_NORMAL_WRITE));
//         break;
//     case TRG_TRAJ_POSITIONS:
//         break;
//     case TRG_TRAJ_VELOCITIES:
//         break;
//     case TRG_TRAJ_FORCES:
//         break;
//     case TRG_TRAJ_BOX_SHAPE:
//         break;
//     case TRG_GENERAL_INFO:
//         return(trg_write_general_info_block(trg_data, block,
//                                             TRG_NORMAL_WRITE));
//     case TRG_MOLECULES:
//         return(trg_write_molecules_block(trg_data, block,
//                                          TRG_NORMAL_WRITE));
//     default:
//         if(block->id > TRG_TRAJ_FORCES)
//         {
//             /* Implement writing data blocks. */
//             return(TRG_FAILURE);
//         }
//         else
//         {
//             return(TRG_FAILURE);
//         }
//     }
// 
//     /* FIXME: Never reached. */
//     if(close_file)
//     {
//         fclose(trg_data->output_file);
//         trg_data->output_file = 0;
//     }
// 
//     return(TRG_SUCCESS);
// }

trg_function_status trg_read_next_frame_set(struct trg_trajectory *trg_data,
                                            trg_close_file_flag close_file)
{
    long int file_pos;
    struct trg_gen_block block;
    trg_function_status stat = TRG_SUCCESS;

    if(trg_init_input_file(trg_data, FALSE) != TRG_SUCCESS)
    {
        return(TRG_CRITICAL);
    }

    trg_init_block(&block);

    file_pos = trg_data->current_trajectory_frame_set.next_frame_set_file_pos;
    
    if(file_pos > 0)
    {
        fseek(trg_data->input_file,
              file_pos,
              SEEK_SET);
    }
    else
    {
        return(TRG_CRITICAL);
    }
    
    if(!trg_data->input_file_len)
    {
        fseek(trg_data->input_file, 0, SEEK_END);
        trg_data->input_file_len = ftell(trg_data->input_file);
        fseek(trg_data->input_file, file_pos, SEEK_SET);
    }

    stat = trg_read_block_header(trg_data, &block);
    if(stat == TRG_CRITICAL || block.id != TRG_TRAJECTORY_FRAME_SET)
    {
        return(TRG_CRITICAL);
    }

    trg_data->current_trajectory_frame_set_input_file_pos = file_pos;
    
    if(trg_read_next_block(trg_data, &block,
                           TRG_KEEP_FILE_OPEN) == TRG_SUCCESS)
    {
        file_pos = ftell(trg_data->input_file);
        /* Read all blocks until next frame set block */
        stat = trg_read_block_header(trg_data, &block);
        while(file_pos < trg_data->input_file_len &&
              stat != TRG_CRITICAL &&
              block.id != TRG_TRAJECTORY_FRAME_SET)
        {
            stat = trg_read_next_block(trg_data, &block,
                                       TRG_KEEP_FILE_OPEN) == TRG_SUCCESS;

            if(stat != TRG_CRITICAL)
            {
                file_pos = ftell(trg_data->input_file);
                if(file_pos < trg_data->input_file_len)
                {
                    stat = trg_read_block_header(trg_data, &block);
                }
            }
        }
        if(stat == TRG_CRITICAL)
        {
            return(stat);
        }

        if(block.id == TRG_TRAJECTORY_FRAME_SET)
        {
            fseek(trg_data->input_file, file_pos, SEEK_SET);
        }
    }

    trg_data->input_file_pos=ftell(trg_data->input_file);

    trg_destroy_block(&block);
    
    if(close_file)
    {
        fclose(trg_data->input_file);
        trg_data->input_file = 0;
    }
    
    return(TRG_SUCCESS);
}

trg_function_status trg_write_frame_set(struct trg_trajectory *trg_data,
                                        trg_close_file_flag close_file)
{
    int i, j;
    struct trg_gen_block block;
    struct trg_trajectory_frame_set *frame_set =
    &trg_data->current_trajectory_frame_set;
    
    trg_function_status stat;
    

    if(trg_data->output_file)
    {
        trg_data->current_trajectory_frame_set_output_file_pos =
        ftell(trg_data->output_file);
    }
    else
    {
        trg_data->current_trajectory_frame_set_output_file_pos =
        trg_data->output_file_pos;
    }
    
    trg_init_block(&block);
    block.id = TRG_TRAJECTORY_FRAME_SET;

    trg_write_frame_set_block(trg_data, &block, TRG_NORMAL_WRITE);

    if(frame_set->contents.n_blocks > 0)
    {
        block.id = TRG_BLOCK_TABLE_OF_CONTENTS;
        trg_write_trajectory_toc_block(trg_data, &block, TRG_NORMAL_WRITE);
    }
    for(i = 0; i<frame_set->n_data_blocks; i++)
    {
        block.id = frame_set->tr_data[i].block_id;
        trg_write_data_block(trg_data, &block, i, TRG_NORMAL_WRITE);
    }
    if(frame_set->n_mapping_blocks)
    {
        for(i = 0; i < frame_set->n_mapping_blocks; i++)
        {
            block.id = TRG_PARTICLE_MAPPING;
            if(frame_set->mappings[i].n_particles > 0)
            {
                trg_write_trajectory_mapping_block(trg_data, &block, i,
                                                   TRG_NORMAL_WRITE);
                for(j = 0; j<frame_set->n_particle_data_blocks; j++)
                {
                    block.id = frame_set->tr_particle_data[i].block_id;
                    trg_write_particle_data_block(trg_data, &block,
                                                  j, &frame_set->mappings[i],
                                                  TRG_NORMAL_WRITE);
                }
            }
        }
    }
    else
    {
        for(i = 0; i<frame_set->n_particle_data_blocks; i++)
        {
            block.id = frame_set->tr_particle_data[i].block_id;
            trg_write_particle_data_block(trg_data, &block,
                                          i, 0, TRG_NORMAL_WRITE);
        }
    }

    trg_data->output_file_pos = ftell(trg_data->output_file);

    stat = trg_update_header_pointers(trg_data);
    
    if(stat == TRG_SUCCESS)
    {
        stat = trg_update_frame_set_pointers(trg_data);
    }
    
    trg_destroy_block(&block);
    
    if(close_file)
    {
        fclose(trg_data->input_file);
        trg_data->input_file = 0;
    }

    return(stat);
}

trg_function_status trg_new_frame_set(struct trg_trajectory *trg_data,
                                      const uint64_t first_frame,
                                      const uint64_t n_frames)
{
    int i;
    struct trg_gen_block block;
    struct trg_trajectory_frame_set *frame_set;
    struct trg_particle_mapping *mapping;
    FILE *temp = trg_data->input_file;
    uint64_t curr_pos;

    frame_set = &trg_data->current_trajectory_frame_set;

    if(trg_data->n_trajectory_frame_sets)
    {
        frame_set->prev_frame_set_file_pos =
        trg_data->current_trajectory_frame_set_output_file_pos;
    }

    if(trg_data->output_file)
    {
        trg_data->current_trajectory_frame_set_output_file_pos =
        ftell(trg_data->output_file);
    }
    else
    {
        trg_data->current_trajectory_frame_set_output_file_pos =
        trg_data->output_file_pos;
    }

    if(frame_set->n_mapping_blocks && frame_set->mappings)
    {
        for(i = frame_set->n_mapping_blocks; i--;)
        {
            mapping = &frame_set->mappings[i];
            if(mapping->real_particle_numbers)
            {
                free(mapping->real_particle_numbers);
                mapping->real_particle_numbers = 0;
            }
        }
        free(frame_set->mappings);
        frame_set->n_mapping_blocks = 0;
    }

    trg_data->n_trajectory_frame_sets++;

    /* Set the long range pointers */
    if(trg_data->n_trajectory_frame_sets == trg_data->stride_length + 1)
    {
        frame_set->long_stride_prev_frame_set_file_pos =
        trg_data->first_trajectory_frame_set_output_file_pos;
    }
    else if(trg_data->n_trajectory_frame_sets > trg_data->stride_length + 1)
    {
        /* FIXME: Currently only working if the previous frame set has its
         * long stride pointer already set. This might need some fixing. */
        if(frame_set->long_stride_prev_frame_set_file_pos != -1ULL &&
           frame_set->long_stride_prev_frame_set_file_pos != 0)
        {
            trg_init_block(&block);
            trg_data->input_file = trg_data->output_file;
            
            curr_pos = ftell(trg_data->output_file);
            fseek(trg_data->output_file,
                  frame_set->long_stride_prev_frame_set_file_pos,
                  SEEK_SET);
            
            if(trg_read_block_header(trg_data, &block) != TRG_SUCCESS)
            {
                printf("Cannot read frame header. %s: %d\n",
                    __FILE__, __LINE__);
                trg_destroy_block(&block);
                trg_data->input_file = temp;
                return(TRG_CRITICAL);
            }

            /* Read the next frame set from the previous frame set and one
             * long stride step back */
            fseek(trg_data->output_file, block.block_contents_size - 4 *
                sizeof(uint64_t), SEEK_CUR);
            if(fread(&frame_set->long_stride_prev_frame_set_file_pos,
                  sizeof(frame_set->long_stride_prev_frame_set_file_pos),
                  1, trg_data->output_file) == 0)
            {
                printf("Cannot read block. %s: %d\n", __FILE__, __LINE__);
                trg_destroy_block(&block);
                trg_data->input_file = temp;
                return(TRG_CRITICAL);
            }

            if(trg_data->endianness_64 != TRG_BIG_ENDIAN_64)
            {
                if(trg_swap_byte_order_64(trg_data,
                                          &frame_set->
                                          long_stride_prev_frame_set_file_pos)
                   != TRG_SUCCESS)
                {
                    printf("Cannot swap byte order to get big endian. %s: %d\n",
                            __FILE__, __LINE__);
                }
            }

            trg_destroy_block(&block);
            trg_data->input_file = temp;
            fseek(trg_data->output_file, curr_pos, SEEK_SET);
        }
    }

    frame_set->first_frame = first_frame;
    frame_set->n_frames = n_frames;
//     frame_set->n_particle_data_blocks = 0;
//     frame_set->n_data_blocks = 0;

    if(trg_data->first_trajectory_frame_set_output_file_pos == -1ULL ||
       trg_data->first_trajectory_frame_set_output_file_pos == 0)
    {
        trg_data->first_trajectory_frame_set_output_file_pos =
        trg_data->current_trajectory_frame_set_output_file_pos;
    }
    /* FIXME: Should check the frame number instead of the file_pos,
     * in case frame sets are not in order */
    if(trg_data->last_trajectory_frame_set_output_file_pos == -1ULL ||
       trg_data->last_trajectory_frame_set_output_file_pos == 0 ||
       trg_data->last_trajectory_frame_set_output_file_pos <
       trg_data->current_trajectory_frame_set_output_file_pos)
    {
        trg_data->last_trajectory_frame_set_output_file_pos =
        trg_data->current_trajectory_frame_set_output_file_pos;
    }

    return(TRG_SUCCESS);
}

/* UNTESTED */
trg_function_status trg_add_data_block(struct trg_trajectory *trg_data,
                                       const uint64_t id,
                                       const char *block_name,
                                       const char datatype,
                                       const uint64_t n_frames,
                                       const uint64_t n_values_per_frame,
                                       const uint64_t stride_length,
                                       uint64_t codec_id,
                                       void *new_data)
{
    int i, j, block_index, size, len;
    struct trg_trajectory_frame_set *frame_set;
    struct trg_data *data;
    void *orig;
    trg_block_type block_type_flag;

    frame_set = &trg_data->current_trajectory_frame_set;

    if(trg_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        block_type_flag = TRG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TRG_NON_TRAJECTORY_BLOCK;
    }

    block_index = -1;
    /* See if there is already a data block of this ID */
    if(block_type_flag == TRG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_data_blocks; i-- ;)
        {
            data = &frame_set->tr_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        for(i = trg_data->n_data_blocks; i-- ;)
        {
            data = &trg_data->non_tr_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }

    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(trg_create_data_block(trg_data, block_type_flag) !=
            TRG_SUCCESS)
        {
            printf("Cannot create data block. %s: %d\n", __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        if(block_type_flag == TRG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_data[frame_set->n_data_blocks - 1];
        }
        else
        {
            data = &trg_data->non_tr_data[trg_data->n_data_blocks - 1];
        }
        data->block_id = id;

        data->block_name = malloc(strlen(block_name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block_name)+1, __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        strncpy(data->block_name, block_name, strlen(block_name) + 1);

        data->datatype = datatype;

        data->stride_length = 1;
        data->values = 0;
        data->n_values_per_frame = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
    }

    /* Allocate memory */
    if(!data->values || data->n_frames != n_frames ||
       data->n_values_per_frame != n_values_per_frame)
    {
        if(trg_allocate_data_mem(trg_data, data, n_frames,
                                 n_values_per_frame) !=
           TRG_SUCCESS)
        {
            printf("Cannot allocate data memory. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    switch(datatype)
    {
    case TRG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TRG_INT_DATA:
        size = sizeof(uint64_t);
        break;
    case TRG_CHAR_DATA:
        size = sizeof(char);
    case TRG_DOUBLE_DATA:
    default:
        size = sizeof(double);
        break;
    }

    orig = new_data;

    for(i = 0; i < n_frames; i++)
    {
        for(j = 0; j < n_values_per_frame; j++)
        {
            /* FIXME: Not optimal to have this switch statement inside
                * the for loops. Check if it should be moved outside. */
            switch(datatype)
            {
            case TRG_CHAR_DATA:
                len = min(strlen(new_data) + 1,
                            TRG_MAX_STR_LEN);
                if(data->values[i][j].c)
                {
                    free(data->values[i][j].c);
                }
                data->values[i][j].c = malloc(len);
                if(!data->values[i][j].c)
                {
                    printf("Cannot allocate memory (%d bytes). %s: %d\n",
                        len, __FILE__, __LINE__);
                    return(TRG_CRITICAL);
                }
                strncpy(data->values[i][j].c,
                        new_data, len);
                new_data += len;
                break;
            case TRG_INT_DATA:
                memcpy(&data->values[i][j].i,
                        new_data, size);
                new_data += size;
                break;
            case TRG_FLOAT_DATA:
                memcpy(&data->values[i][j].f,
                        new_data, size);
                new_data += size;
                break;
            case TRG_DOUBLE_DATA:
            default:
                memcpy(&data->values[i][j].d,
                        new_data, size);
                new_data += size;
            }
        }
    }

    new_data = orig;

    return(TRG_SUCCESS);
}


trg_function_status trg_add_particle_data_block(struct trg_trajectory *trg_data,
                                        const uint64_t id,
                                        const char *block_name,
                                        const char datatype,
                                        const uint64_t n_frames,
                                        const uint64_t n_values_per_frame,
                                        const uint64_t stride_length,
                                        const uint64_t first_particle_number,
                                        const uint64_t n_particles,
                                        uint64_t codec_id,
                                        void *new_data)
{
    int i, j, k, block_index, size, len;
    uint64_t tot_n_particles;
    struct trg_trajectory_frame_set *frame_set;
    struct trg_particle_data *data;
    void *orig;
    trg_block_type block_type_flag;

    frame_set = &trg_data->current_trajectory_frame_set;
    
    if(trg_data->current_trajectory_frame_set_output_file_pos > 0)
    {
        block_type_flag = TRG_TRAJECTORY_BLOCK;
    }
    else
    {
        block_type_flag = TRG_NON_TRAJECTORY_BLOCK;
    }

    block_index = -1;
    /* See if there is already a data block of this ID */
    if(block_type_flag == TRG_TRAJECTORY_BLOCK)
    {
        for(i = frame_set->n_particle_data_blocks; i-- ;)
        {
            data = &frame_set->tr_particle_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }
    else
    {
        for(i = trg_data->n_particle_data_blocks; i-- ;)
        {
            data = &trg_data->non_tr_particle_data[i];
            if(data->block_id == id)
            {
                block_index = i;
                break;
            }
        }
    }

    /* Otherwise create a data block */
    if(block_index == -1)
    {
        if(trg_create_particle_data_block(trg_data, block_type_flag) !=
            TRG_SUCCESS)
        {
            printf("Cannot create particle data block. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        if(block_type_flag == TRG_TRAJECTORY_BLOCK)
        {
            data = &frame_set->tr_particle_data[frame_set->
                                                n_particle_data_blocks - 1];
        }
        else
        {
            data = &trg_data->non_tr_particle_data[trg_data->
                                                   n_particle_data_blocks - 1];
        }
        data->block_id = id;

        data->block_name = malloc(strlen(block_name) + 1);
        if(!data->block_name)
        {
            printf("Cannot allocate memory (%d bytes). %s: %d\n",
                   (int)strlen(block_name)+1, __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
        strncpy(data->block_name, block_name, strlen(block_name) + 1);

        data->datatype = datatype;

        data->stride_length = 1;
        data->values = 0;
        data->n_values_per_frame = 0;
        data->n_frames = 0;
        data->codec_id = codec_id;
    }

    if(block_type_flag == TRG_TRAJECTORY_BLOCK && trg_data->var_num_atoms_flag)
    {
        tot_n_particles = frame_set->n_particles;
    }
    else
    {
        tot_n_particles = trg_data->n_particles;
    }

    /* Allocate memory */
    if(!data->values || data->n_frames != n_frames ||
       data->n_values_per_frame != n_values_per_frame)
    {
        if(trg_allocate_particle_data_mem(trg_data, data, n_frames,
                                          tot_n_particles,
                                          n_values_per_frame) !=
           TRG_SUCCESS)
        {
            printf("Cannot allocate particle data memory. %s: %d\n",
                   __FILE__, __LINE__);
            return(TRG_CRITICAL);
        }
    }

    switch(datatype)
    {
    case TRG_FLOAT_DATA:
        size = sizeof(float);
        break;
    case TRG_INT_DATA:
        size = sizeof(uint64_t);
        break;
    case TRG_CHAR_DATA:
        size = sizeof(char);
    case TRG_DOUBLE_DATA:
    default:
        size = sizeof(double);
        break;
    }

    orig = new_data;

    for(i = 0; i < n_frames; i++)
    {
        for(j = first_particle_number; j < n_particles; j++)
        {
            for(k = 0; k < n_values_per_frame; k++)
            {
                /* FIXME: Not optimal to have this switch statement inside
                 * the for loops. Check if it should be moved outside. */
                switch(datatype)
                {
                case TRG_CHAR_DATA:
                    len = min(strlen(new_data) + 1,
                              TRG_MAX_STR_LEN);
                    if(data->values[i][j][k].c)
                    {
                        free(data->values[i][j][k].c);
                    }
                    data->values[i][j][k].c = malloc(len);
                    if(!data->values[i][j][k].c)
                    {
                        printf("Cannot allocate memory (%d bytes). %s: %d\n",
                            len, __FILE__, __LINE__);
                        return(TRG_CRITICAL);
                    }
                    strncpy(data->values[i][j][k].c,
                            new_data, len);
                    new_data += len;
                    break;
                case TRG_INT_DATA:
                    memcpy(&data->values[i][j][k].i,
                           new_data, size);
                    new_data += size;
                    break;
                case TRG_FLOAT_DATA:
                    memcpy(&data->values[i][j][k].f,
                           new_data, size);
                    new_data += size;
                    break;
                case TRG_DOUBLE_DATA:
                default:
                    memcpy(&data->values[i][j][k].d,
                           new_data, size);
                    new_data += size;
                }
            }
        }
    }

    new_data = orig;
    
    return(TRG_SUCCESS);
}
                                                

trg_function_status trg_read_next_traj_block(struct trg_trajectory *trg_data,
                                             trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}

trg_function_status trg_write_next_traj_block(struct trg_trajectory *trg_data,
                                              trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}

trg_function_status trg_read_traj_block(struct trg_trajectory *trg_data,
                                        uint64_t block_id,
                                        trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}
        
trg_function_status trg_write_traj_block(struct trg_trajectory *trg_data,
                                         uint64_t block_id,
                                         trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}

trg_function_status trg_read_frame_nr(struct trg_trajectory *trg_data,
                                      uint64_t frame_nr,
                                      trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}
        
trg_function_status trg_write_frame_nr(struct trg_trajectory *trg_data,
                                       uint64_t frame_nr,
                                       trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}

trg_function_status trg_read_frame_nrs(struct trg_trajectory *trg_data,
                                       uint64_t start_frame_nr,
                                       uint64_t end_frame_nr,
                                       trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}

trg_function_status trg_write_frame_nrs(struct trg_trajectory *trg_data,
                                        uint64_t start_frame_nr,
                                        uint64_t end_frame_nr,
                                        trg_close_file_flag close_file)
{
    /* STUB */
    return(TRG_SUCCESS);
}

trg_function_status trg_read_frame_set_nr(struct trg_trajectory *trg_data,
                                          uint64_t frame_set_nr)
{
    return(TRG_SUCCESS);
}


trg_function_status trg_get_time_str(struct trg_trajectory *trg_data,
                                     char *time)
{
    struct tm *time_data;
    time_t secs;

    secs = trg_data->time;

    time_data = localtime(&secs); // Returns a statically allocated variable.
    snprintf(time, TRG_MAX_DATE_STR_LEN,
             "%4d-%02d-%02d %02d:%02d:%02d",
             time_data->tm_year+1900, time_data->tm_mon+1, time_data->tm_mday,
             time_data->tm_hour, time_data->tm_min, time_data->tm_sec);

    return(TRG_SUCCESS);
}

