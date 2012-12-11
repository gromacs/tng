/** @file tng_io.h
 *  @brief API for input and output of tng trajectory files
 */

#ifndef _TNGIO_H
#define _TNGIO_H     1

#include <stdio.h>
#include <inttypes.h>

/** The version of this TNG build */
#define TNG_VERSION 1

/** Flag to indicate particle dependent data. */
#define TNG_PARTICLE_DEPENDENT 1
/** Flag to indicate frame dependent data. */
#define TNG_FRAME_DEPENDENT 2

/** The maximum length of a date string */
#define TNG_MAX_DATE_STR_LEN 24
/** The length of an MD5 hash */
#define TNG_HASH_LEN 16
/** The maximum allowed length of a string */
#define TNG_MAX_STR_LEN 1024


/** Inline function for finding the lowest of two values */
#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a < _b ? _a : _b; })
     
/** Inline function for finding the highest of two values */
#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a > _b ? _a : _b; })


typedef enum {TNG_BIG_ENDIAN_32,
              TNG_LITTLE_ENDIAN_32,
              TNG_BYTE_PAIR_SWAP_32} tng_endianness_32;

typedef enum {TNG_BIG_ENDIAN_64,
              TNG_LITTLE_ENDIAN_64,
              TNG_QUAD_SWAP_64,
              TNG_BYTE_PAIR_SWAP_64,
              TNG_BYTE_SWAP_64} tng_endianness_64;


typedef enum {TNG_UNCOMPRESSED,
              TNG_XTC_COMPRESSION,
              TNG_TNG_COMPRESSION} tng_compression;

typedef enum {TNG_NON_TRAJECTORY_BLOCK, TNG_TRAJECTORY_BLOCK} tng_block_type;

typedef enum {TNG_ENDIANNESS_AND_STRING_LENGTH,
              TNG_GENERAL_INFO,
              TNG_MOLECULES,
              TNG_TRAJECTORY_IDS_AND_NAMES,
              TNG_TRAJECTORY_FRAME_SET,
              TNG_BLOCK_TABLE_OF_CONTENTS,
              TNG_PARTICLE_MAPPING} tng_non_trajectory_block_ids;

typedef enum {TNG_TRAJ_BOX_SHAPE = 10000,
              TNG_TRAJ_POSITIONS,
              TNG_TRAJ_VELOCITIES,
              TNG_TRAJ_FORCES} tng_trajectory_block_ids;
              

typedef enum {TNG_NON_PARTICLE_BLOCK_DATA,
              TNG_PARTICLE_BLOCK_DATA} tng_particle_block_data;
              
/*typedef enum {TNG_NO_HASH, TNG_OTHER_HASH,
              TNG_MD5_HASH, TNG_MD6_HASH,
              TNG_SHA0_HASH, TNG_SHA1_HASH,
              TNG_SHA256_HASH, TNG_SHA512_HASH} tng_hash_type;*/

typedef enum {FALSE, TRUE} tng_bool;

typedef enum {TNG_CONSTANT_N_ATOMS, TNG_VARIABLE_N_ATOMS}
             tng_variable_n_atoms_flag;

typedef enum {TNG_SUCCESS, TNG_FAILURE, TNG_CRITICAL} tng_function_status;

typedef enum {TNG_NORMAL_WRITE, TNG_COPY_EXISTING} tng_write_mode;

typedef enum {TNG_SKIP_HASH, TNG_USE_HASH} tng_hash_mode;

typedef enum {TNG_CHAR_DATA,
              TNG_INT_DATA,
              TNG_FLOAT_DATA,
              TNG_DOUBLE_DATA} tng_data_type;




typedef struct tng_trajectory *tng_trajectory_t;
typedef struct tng_molecule *tng_molecule_t;
typedef struct tng_chain *tng_chain_t;
typedef struct tng_residue *tng_residue_t;
typedef struct tng_atom *tng_atom_t;
typedef struct tng_bond *tng_bond_t;
              
#ifdef __cplusplus
extern "C"
{
#endif

struct tng_bond {
    /** One of the atoms of the bond */
    int64_t from_atom_id;
    /** The other atom of the bond */
    int64_t to_atom_id;
};

struct tng_atom {
    /** The residue containing this atom */
    tng_residue_t residue;
    /** A unique (per molecule) ID number of the atom */
    int64_t id;
    /** The atom_type (depending on the forcefield) */
    char *atom_type;
    /** The name of the atom */
    char *name;
};

struct tng_residue {
    /** The chain containing this residue */
    tng_chain_t chain;
    /** A unique (per chain) ID number of the residue */
    int64_t id;
    /** The name of the residue */
    char *name;
    /** The number of atoms in the residue */
    int64_t n_atoms;
    /** A list of atoms in the residue */
    tng_atom_t atoms;
};

struct tng_chain {
    /** The molecule containing this chain */
    tng_molecule_t molecule;
    /** A unique (per molecule) ID number of the chain */
    int64_t id;
    /** The name of the chain */
    char *name;
    /** The number of residues in the chain */
    int64_t n_residues;
    /** A list of residues in the chain */
    tng_residue_t residues;
};

struct tng_molecule {
    /** A unique ID number of the molecule */
    int64_t id;
    /** Quaternary structure of the molecule.
     *  1 => monomeric
     *  2 => dimeric
     *  3 => trimeric
     *  etc */
    int64_t quaternary_str;
    /** The number of chains in the molecule */
    int64_t n_chains;
    /** The number of residues in the molecule */
    int64_t n_residues;
    /** The number of atoms in the molecule */
    int64_t n_atoms;
    /** The number of bonds in the molecule. If the bonds are not specified this
     * value can be 0. */
    int64_t n_bonds;
    /** The name of the molecule */
    char *name;
    /** A list of chains in the molecule */
    tng_chain_t chains;
    /** A list of residues in the molecule */
    tng_residue_t residues;
    /** A list of the atoms in the molecule */
    tng_atom_t atoms;
    /** A list of the bonds in the molecule */
    tng_bond_t bonds;
};

struct tng_gen_block {
    /** The size of the block header in bytes */
    int64_t header_contents_size;
    /** The size of the block contents in bytes */
    int64_t block_contents_size;
    /** The ID of the block to determine its type */
    int64_t id;
    /** The MD5 hash of the block to verify integrity */
    char hash[TNG_HASH_LEN];
    /** The name of the block */
    char *name;
    /** The library version used to write the block */
    int64_t block_version;
    /** The full block header contents */
    char *header_contents;
    /** The full block contents */
    char *block_contents;
};

struct tng_frame_set_toc {
    /** The number of blocks listed in this table of contents */
    int64_t n_blocks;
    /** A list of block names */
    char **block_names;
};

struct tng_particle_mapping {
    /** The index number of the first particle in this mapping block */
    int64_t num_first_particle;
    /** The number of particles list in this mapping block */
    int64_t n_particles;
    /** the mapping of index numbers to the real particle numbers in the
     * trajectory */
    int64_t *real_particle_numbers;
};

struct tng_trajectory_frame_set {
    /** The table of contents of this frame set */
    struct tng_frame_set_toc contents;
    /** The number of different particle mapping blocks present. */
    int64_t n_mapping_blocks;
    /** The atom mappings of this frame set */
    struct tng_particle_mapping *mappings;
    /** The first frame of this frame set */
    int64_t first_frame;
    /** The number of frames in this frame set */
    int64_t n_frames;
    /** A list of the number of each molecule type - only used when using
     * variable number of atoms */
    int64_t *molecule_cnt_list;
    /** The number of particles/atoms - only used when using variable number
     * of atoms */
    int64_t n_particles;
    /** The file position of the next frame set */
    int64_t next_frame_set_file_pos;
    /** The file position of the previous frame set */
    int64_t prev_frame_set_file_pos;
    /** The file position of the frame set one long stride step ahead */
    int64_t medium_stride_next_frame_set_file_pos;
    /** The file position of the frame set one long stride step behind */
    int64_t medium_stride_prev_frame_set_file_pos;
    /** The file position of the frame set one long stride step ahead */
    int64_t long_stride_next_frame_set_file_pos;
    /** The file position of the frame set one long stride step behind */
    int64_t long_stride_prev_frame_set_file_pos;

    /* The data blocks in a frame set are trajectory data blocks */
    /** The number of trajectory data blocks of particle dependent data */
    int n_particle_data_blocks;
    /** A list of data blocks containing particle dependent data */
    struct tng_particle_data *tr_particle_data;
    /** The number of trajectory data blocks independent of particles */
    int n_data_blocks;
    /** A list of data blocks containing particle indepdendent data */
    struct tng_data *tr_data;
};

/** Data can be either double, float, int or a string */
union data_values {
    double d;
    float f;
    int i;
    char *c;
};

/* FIXME: Should there be a pointer to a tng_gen_block from each data block? */
struct tng_particle_data {
    /** The block ID of the data block containing this particle data.
     *  This is used to determine the kind of data that is stored */
    int64_t block_id;
    /** The name of the data block. This is used to determine the kind of
     *  data that is stored */
    char *block_name;
    /** The type of data stored. */
    tng_data_type datatype;
    /** The first frame number of the first data value */
    int64_t first_frame_with_data;
    /** The number of frames in this frame set */
    int64_t n_frames;
    /** The number of values stored per frame */
    int64_t n_values_per_frame;
    /** The number of frames between each data point - e.g. when
     *  storing sparse data. */
    int64_t stride_length;
    /** ID of the CODEC used for compression 0 == no compression. */
    int64_t codec_id;
    /** The multiplier used for getting integer values for compression */
    double compression_multiplier;
    /** A 3-dimensional array of values, sized
     *  n_frames * n_particles * n_values_per_frame */
    union data_values ***values;
};

struct tng_data {
    /** The ID of the data block */
    int64_t block_id;
    /** The name of the data block. This is used to determine the kind of
     *  data that is stored */
    char *block_name;
    /** The type of data stored. */
    tng_data_type datatype;
    /** The first frame number of the first data value */
    int64_t first_frame_with_data;
    /** The number of frames in this data block */
    int64_t n_frames;
    /** The number of values stored per frame */
    int64_t n_values_per_frame;
    /** The number of frames between each data value, e.g. if storing data
     *  that is not saved every frame. */
    int64_t stride_length;
    /** ID of the CODEC used for compression. 0 == no compression. */
    int64_t codec_id;
    /** Compressed data is stored as integers. This compression multiplier is
     *  the multiplication factor to convert from integer to float/double */
    double compression_multiplier;
    /** A 2-dimensional array of values, sized n_frames * n_values_per_frame */
    union data_values **values;
};



struct tng_trajectory {
    /** The path of the input trajectory file */
    char *input_file_path;
    /** A handle to the input file */
    FILE *input_file;
    /** The reading position of the file */
    long int input_file_pos;
    /** The length of the input file */
    long int input_file_len;
    /** The path of the output trajectory file */
    char *output_file_path;
    /** A handle to the output file */
    FILE *output_file;
    /** The writing position of the file */
    long int output_file_pos;        
    /** The endianness of 32 bit values of the current computer */
    tng_endianness_32 endianness_32;
    /** The endianness of 64 bit values of the current computer */
    tng_endianness_64 endianness_64;

    /** The name of the program producing this trajectory */
    char *program_name;
    /** The forcefield used in the simulations */
    char *forcefield_name;
    /** The name of the user running the simulations */
    char *user_name;
    /** The time (n seconds since 1970) when the file was created */
    int64_t time;
    /** The name of the computer on which the simulations were performed */
    char *computer_name;
    /** The PGP signature of the user creating the file. */
    char *pgp_signature;
    /** A flag indicating if the number of atoms can vary throughout the
     *  simulation, e.g. using a grand canonical ensemble */
    char var_num_atoms_flag;
    /** The number of frames in a frame set. It is allowed to have frame sets
     *  with fewer frames, but this will help searching for specific frames */
    int64_t frame_set_n_frames;
    /** The number of frame sets in a medium stride step */
    int64_t medium_stride_length;
    /** The number of frame sets in a long stride step */
    int64_t long_stride_length;
    
    /** The number of different kinds of molecules in the trajectory */
    int64_t n_molecules;
    /** A list of molecules in the trajectory */
    tng_molecule_t molecules;
    /** A list of the count of each molecule - if using variable number of
     *  particles this will be specified in each frame set */
    int64_t *molecule_cnt_list;
    /** The total number of particles/atoms. If using variable number of
     *  particles this will be specified in each frame set */
    int64_t n_particles;

     /** The pos in the src file of the first frame set */
    int64_t first_trajectory_frame_set_input_file_pos;
    /** The pos in the dest file of the first frame set */
    int64_t first_trajectory_frame_set_output_file_pos;
    /** The pos in the src file of the last frame set */
    int64_t last_trajectory_frame_set_input_file_pos;
    /** The pos in the dest file of the last frame set */
    int64_t last_trajectory_frame_set_output_file_pos;
    /** The currently active frame set */
    struct tng_trajectory_frame_set current_trajectory_frame_set;
    /** The pos in the src file of the current frame set */
    long int current_trajectory_frame_set_input_file_pos;
    /** The pos in the dest file of the current frame set */
    long int current_trajectory_frame_set_output_file_pos;
    /** The number of frame sets in the trajectory */
    int64_t n_trajectory_frame_sets;

    /** The number of trajectory blocks in the file */
    int64_t n_trajectory_blocks;
    /** The number of non-trajectory blocks in the file */
    int64_t n_non_trajectory_blocks;
    /** A list of non-trajectory blocks */
    struct tng_gen_block non_trajectory_blocks[32];
                                                       
    /* These data blocks are non-trajectory data blocks */
    /** The number of non-frame dependent particle dependent data blocks */
    int n_particle_data_blocks;
    /** A list of data blocks containing particle dependent data */
    struct tng_particle_data *non_tr_particle_data;

    /** The number of frame and particle independent data blocks */
    int n_data_blocks;
    /** A list of frame and particle indepdendent data blocks */
    struct tng_data *non_tr_data;
};



/**
 * @brief Setup a trajectory data container.
 * @param tng_data pre-allocated memory to initialise as a trajectory.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_trajectory_init(tng_trajectory_t tng_data);

/**
 * @brief Clean up a trajectory data container.
 * @param tng_data the trajectory data to destroy.
 * @details All allocated memory in the data structure is freed, but not the memory
 * of tng_data itself.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_trajectory_destroy(tng_trajectory_t tng_data);

/**
 * @brief Set the name of the input file.
 * @param tng_data the trajectory of which to set the input file name.
 * @param file_name the name of the input file.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_input_file_set(tng_trajectory_t tng_data,
                                       const char *file_name);

/**
 * @brief Set the name of the output file.
 * @param tng_data the trajectory of which to set the output file name.
 * @param file_name the name of the output file.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_output_file_set(tng_trajectory_t tng_data,
                                        const char *file_name);

/**
 * @brief Set the name of the program used when creating the trajectory.
 * @param tng_data the trajectory of which to set the program name.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_program_name_set(tng_trajectory_t tng_data,
                                         const char *new_name);

/**
 * @brief Set the name of the forcefield used in the trajectory.
 * @param tng_data the trajectory of which to set the forcefield name.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_forcefield_name_set(tng_trajectory_t tng_data,
                                            const char *new_name);

/**
 * @brief Set the name of the user who created the trajectory.
 * @param tng_data the trajectory of which to set the user name.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_user_name_set(tng_trajectory_t tng_data,
                                      const char *new_name);

/**
 * @brief Set the name of the computer used when creating the trajectory.
 * @param tng_data the trajectory of which to set the computer name.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_computer_name_set(tng_trajectory_t tng_data,
                                          const char *new_name);

/**
 * @brief Set the pgp_signature of the trajectory.
 * @param tng_data the trajectory of which to set the computer name.
 * @param signature is a string containing the pgp_signature.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_signature_set(tng_trajectory_t tng_data,
                                      const char *signature);

/**
 * @brief Setup a molecule container.
 * @param molecule is the molecule to initialise. Memory must be preallocated.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_molecule_init(tng_molecule_t molecule);

/**
 * @brief Clean up a molecule container.
 * @param molecule is the molecule to destroy.
 * @details All allocated memory in the data structure is freed, but not the
 * memory of molecule itself.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_molecule_destroy(tng_molecule_t molecule);

/**
 * @brief Setup a data block.
 * @param block is a pointer to pre-allocated memory.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_block_init(struct tng_gen_block *block);

/**
 * @brief Clean up a data block.
 * @param block is a pointer to pre-allocated memory.
 * @details All allocated memory in the data structure is freed, but not the
 * memory of block itself.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_block_destroy(struct tng_gen_block *block);

/**
 * @brief Set the name of a data block.
 * @param tng_data is the trajectory data container containing the block..
 * @param block is a pointer to the block to rename.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_block_name_set(tng_trajectory_t tng_data,
                                       struct tng_gen_block *block,
                                       const char *new_name);

/**
 * @brief Add a molecule to the trajectory.
 * @param tng_data is the trajectory data container containing the block..
 * @param name is a pointer to the string containing the name of the new molecule.
 * @param molecule is a pointer to the newly created molecule.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_molecule_add(tng_trajectory_t tng_data,
                                     const char *name,
                                     tng_molecule_t *molecule);

/**
 * @brief Set the name of a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule to rename.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_molecule_name_set(tng_trajectory_t tng_data,
                                          tng_molecule_t molecule,
                                          const char *new_name);

/**
 * @brief Get the count of a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule of which to get the count.
 * @param cnt is a pointer to the variable to be populated with the count.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_molecule_cnt_get(tng_trajectory_t tng_data,
                                         tng_molecule_t molecule,
                                         int64_t *cnt);

/**
 * @brief Set the count of a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule of which to set the count.
 * @param cnt is the number of instances of this molecule.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_molecule_cnt_set(tng_trajectory_t tng_data,
                                         tng_molecule_t molecule,
                                         int64_t cnt);

/**
 * @brief Add a chain to a molecule.
 * @param tng_data is the trajectory data container containing the molecule..
 * @param molecule is the molecule to add a chain to.
 * @param name is a string containing the name of the chain.
 * @param chain is a pointer to the newly created chain.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_molecule_chain_add(tng_trajectory_t tng_data,
                                              tng_molecule_t molecule,
                                              const char *name,
                                              tng_chain_t *chain);

/**
 * @brief Set the name of a chain.
 * @param tng_data is the trajectory data container containing the atom..
 * @param chain is the chain to rename.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_chain_name_set(tng_trajectory_t tng_data,
                                       tng_chain_t chain,
                                       const char *new_name);

/**
 * @brief Add a residue to a chain.
 * @param tng_data is the trajectory data container containing the chain..
 * @param chain is the chain to add a residue to.
 * @param name is a string containing the name of the residue.
 * @param residue is a pointer to the newly created residue.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_chain_residue_add(tng_trajectory_t tng_data,
                                             tng_chain_t chain,
                                             const char *name,
                                             tng_residue_t *residue);

/**
 * @brief Set the name of a residue.
 * @param tng_data is the trajectory data container containing the residue.
 * @param residue is the residue to rename.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_residue_name_set(tng_trajectory_t tng_data,
                                         tng_residue_t residue,
                                         const char *new_name);

/**
 * @brief Add an atom to a residue.
 * @param tng_data is the trajectory containing the residue.
 * @param residue is the residue to add an atom to.
 * @param atom_name is a string containing the name of the atom.
 * @param atom_type is a string containing the atom type of the atom.
 * @param atom is a pointer to the newly created atom.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_residue_atom_add(tng_trajectory_t tng_data,
                                            tng_residue_t residue,
                                            const char *atom_name,
                                            const char *atom_type,
                                            tng_atom_t *atom);

/**
 * @brief Set the name of an atom.
 * @param tng_data is the trajectory data container containing the atom..
 * @param atom is the atom to rename.
 * @param new_name is a string containing the wanted name.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_atom_name_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_name);

/**
 * @brief Set the atom type of an atom.
 * @param tng_data is the trajectory data container containing the atom..
 * @param atom is the atom to change.
 * @param new_name is a string containing the atom type.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_atom_type_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_type);

/**
 * @brief Read the header blocks from the input_file of tng_data.
 * @details The trajectory blocks must be read separately and iteratively in chunks
 * to fit in memory.
 * @param tng_data is a trajectory data container.
 * @details tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @details If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_file_headers_read(tng_trajectory_t tng_data,
                                          const tng_hash_mode hash_mode);

/**
 * @brief Write the header blocks to the output_file of tng_data.
 * @details The trajectory blocks must be written separately and iteratively in chunks
 * to fit in memory.
 * @param tng_data is a trajectory data container.
 * @details tng_data->output_file_path
 * specifies which file to write to. If the file (output_file) is not open it
 * will be opened.
 * @details If hash_mode == TNG_USE_HASH an md5 hash for each header block will be generated.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_file_headers_write(tng_trajectory_t tng_data,
                                           const tng_hash_mode hash_mode);


/**
 * @brief Read one (the next) block (of any kind) from the input_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_data is a pointer to the struct which will be populated with the
 * data.
 * @details If block_data->input_file_pos > 0 it is the position from where the
 * reading starts otherwise it starts from the current position.
 * @details If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_block_read_next(tng_trajectory_t tng_data,
                                        struct tng_gen_block *block_data,
                                        const tng_hash_mode hash_mode);


/**
 * @brief Read one (the next) frame set, including toc, mapping and related data blocks
 * from the input_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details  tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @details If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
 * compared to the md5 hash of the read contents to ensure valid data.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_frame_set_read_next(tng_trajectory_t tng_data,
                                            const tng_hash_mode hash_mode);

/**
 * @brief Write one frame set, including toc, mapping and related data blocks
 * to the output_file of tng_data.
 * @param tng_data is a trajectory data container.
 * @details  tng_data->output_file_path specifies
 * which file to write to. If the file (output_file) is not open it will be
 * opened.
 * @details If hash_mode == TNG_USE_HASH an md5 hash for each block will be generated.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_frame_set_write(tng_trajectory_t tng_data,
                                        const tng_hash_mode hash_mode);

/**
 * @brief Create and initialise a frame set.
 * @param tng_data is the trajectory data container in which to add the frame
 * set
 * @param first_frame is the first frame of the frame set.
 * @param n_frames is the number of frames in the frame set.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_frame_set_new(tng_trajectory_t tng_data,
                                      const int64_t first_frame,
                                      const int64_t n_frames);

/**
 * @brief Add a non-particle dependent data block.
 * @param tng_data is the trajectory data container in which to add the data
 * block
 * @param id is the block ID of the block to add.
 * @param block_name is a descriptive name of the block to add
 * @param datatype is the datatype of the data in the block (e.g. int/float)
 * @param block_type_flag indicates if this is a non-trajectory block (added
 * directly to tng_data) or if it is a trajectory block (added to the
 * frame set)
 * @param n_frames is the number of frames of the data block (automatically
 * set to 1 if adding a non-trajectory data block)
 * @param n_values_per_frame is how many values a stored each frame (e.g. 9
 * for a box shape block)
 * @param stride_length is how many frames are between each entry in the
 * data block
 * @param codec_id is the ID of the codec to compress the data.
 * @param new_data is an array of data values to add.
 */
tng_function_status tng_data_block_add(tng_trajectory_t tng_data,
                                        const int64_t id,
                                        const char *block_name,
                                        const char datatype,
                                        const tng_block_type block_type_flag,
                                        int64_t n_frames,
                                        const int64_t n_values_per_frame,
                                        const int64_t stride_length,
                                        int64_t codec_id,
                                        void *new_data);
                                       

/**
 * @brief Add a particle dependent data block.
 * @param tng_data is the trajectory data container in which to add the data
 * block
 * @param id is the block ID of the block to add.
 * @param block_name is a descriptive name of the block to add
 * @param datatype is the datatype of the data in the block (e.g. int/float)
 * @param block_type_flag indicates if this is a non-trajectory block (added
 * directly to tng_data) or if it is a trajectory block (added to the
 * frame set)
 * @param n_frames is the number of frames of the data block (automatically
 * set to 1 if adding a non-trajectory data block)
 * @param n_values_per_frame is how many values a stored each frame (e.g. 9
 * for a box shape block)
 * @param stride_length is how many frames are between each entry in the
 * data block
 * @param first_particle_number is the number of the first particle stored
 * in this data block
 * @param n_particles is the number of particles stored in this data block
 * @param codec_id is the ID of the codec to compress the data.
 * @param new_data is an array of data values to add.
 */
tng_function_status tng_particle_data_block_add(tng_trajectory_t tng_data,
                                        const int64_t id,
                                        const char *block_name,
                                        const char datatype,
                                        const tng_block_type block_type_flag,
                                        int64_t n_frames,
                                        const int64_t n_values_per_frame,
                                        const int64_t stride_length,
                                        const int64_t first_particle_number,
                                        const int64_t n_particles,
                                        int64_t codec_id,
                                        void *new_data);


/**
 * @brief Read one trajectory block from the input_file of tng_data.
 * @param tng_data is a trajectory data container. tng_data->input_file_path
 * specifies which file to read from. If the file (input_file) is not open it
 * will be opened.
 * @param block_id is the ID of the block to read.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_traj_block_read(tng_trajectory_t tng_data,
                                        int64_t block_id);
        
/**
 * @brief Write one trajectory block to the output_file of tng_data.
 * @param tng_data is a trajectory data container. tng_data->output_file_path
 * specifies which file to write to. If the file (output_file) is not open it
 * will be opened.
 * @param block_id is the ID of the block to write.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_traj_block_write(tng_trajectory_t tng_data,
                                         int64_t block_id);

/**
 * @brief Read a requested frame set.
 * @param tng_data is a trajectory data container.
 * tng_data->current_trajectory_frame_set will be the read frame set.
 * @param frame_set_nr is the number of the frame set to return (starting from 0).
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_frame_set_read_nr(tng_trajectory_t tng_data,
                                          int64_t frame_set_nr);

/**
 * @brief Read a number of consecutive trajectory frames from the input_file of
 * tng_data.
 * @param tng_data is a trajectory data container.
 * tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param start_frame_nr is the index number of the first frame to read.
 * @param end_frame_nr is the index number of the last frame to read.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_frame_read_interval(tng_trajectory_t tng_data,
                                       int64_t start_frame_nr,
                                       int64_t end_frame_nr);

/**
 * @brief Write a number of consecutive trajectory frames to the output_file of
 * tng_data.
 * @param tng_data is a trajectory data container. tng_data->output_file_path specifies
 * which file to write to. If the file (output_file) is not open it will be
 * opened.
 * @param start_frame_nr is the index number of the first frame to write.
 * @param end_frame_nr is the index number of the last frame to write.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_frame_write_interval(tng_trajectory_t tng_data,
                                        int64_t start_frame_nr,
                                        int64_t end_frame_nr);

/**
 * @brief Free data is an array of values (2D).
 * @param values is the 2D array to free and will be set to 0 afterwards.
 * @param n_frames is the number of frames in the data array.
 * @param n_values_per_frame is the number of values per frame in the data array.
 * @param type is the data type of the data in the array (e.g. int/float/char).
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_data_values_free(union data_values **values,
                                         int64_t n_frames,
                                         int64_t n_values_per_frame,
                                         tng_data_type type);

/**
 * @brief Free data is an array of values (3D).
 * @param values is the array to free and will be set to 0 afterwards.
 * @param n_frames is the number of frames in the data array.
 * @param n_particles is the number of particles in the data array.
 * @param n_values_per_frame is the number of values per frame in the data array.
 * @param type is the data type of the data in the array (e.g. int/float/char).
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_particle_data_values_free(union data_values ***values,
                                                  int64_t n_frames,
                                                  int64_t n_particles,
                                                  int64_t n_values_per_frame,
                                                  tng_data_type type);

/**
 * @brief Retrieve non-particle data, from the last read frame set.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param values is a pointer to a 2-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_values_per_frame).
 * Since ***values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_frames is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_data_get(tng_trajectory_t tng_data,
                                 int64_t block_id,
                                 union data_values ***values,
                                 int64_t *n_frames,
                                 int64_t *n_values_per_frame,
                                 tng_data_type *type);

/**
 * @brief Read and retrieve non-particle data, in a specific interval.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param start_frame_nr is the index number of the first frame to read.
 * @param end_frame_nr is the index number of the last frame to read.
 * @param values is a pointer to a 2-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_values_per_frame).
 * Since ***values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_data_interval_get(tng_trajectory_t tng_data,
                                          int64_t block_id,
                                          int64_t start_frame_nr,
                                          int64_t end_frame_nr,
                                          union data_values ***values,
                                          int64_t *n_values_per_frame,
                                          tng_data_type *type);

/**
 * @brief Retrieve particle data, from the last read frame set.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param values is a pointer to a 3-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_particles * n_values_per_frame).
 * Since ****values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_frames is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_particles is set to the number of particles in the returned data. This is
 * needed to properly reach and/or free the data afterwards.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_particle_data_get(tng_trajectory_t tng_data,
                                          int64_t block_id,
                                          union data_values ****values,
                                          int64_t *n_frames,
                                          int64_t *n_particles,
                                          int64_t *n_values_per_frame,
                                          tng_data_type *type);

/**
 * @brief Read and retrieve particle data, in a specific interval.
 * @param tng_data is a trajectory data container. tng_data->input_file_path specifies
 * which file to read from. If the file (input_file) is not open it will be
 * opened.
 * @param block_id is the id number of the particle data block to read.
 * @param start_frame_nr is the index number of the first frame to read.
 * @param end_frame_nr is the index number of the last frame to read.
 * @param first_particle_number is the number of the first particle, for which to
 * retrieve data.
 * @param last_particle_number is the number of the last particle, for which to
 * retrieve data.
 * @param values is a pointer to a 3-dimensional array (memory unallocated), which
 * will be filled with data. The array will be sized
 * (n_frames * n_particles * n_values_per_frame).
 * Since ****values is allocated in this function it is the callers
 * responsibility to free the memory.
 * @param n_values_per_frame is set to the number of values per frame in the data.
 * This is needed to properly reach and/or free the data afterwards.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_particle_data_interval_get(tng_trajectory_t tng_data,
                                                  int64_t block_id,
                                                  int64_t start_frame_nr,
                                                  int64_t end_frame_nr,
                                                  int64_t first_particle_number,
                                                  int64_t last_particle_number,
                                                  union data_values ****values,
                                                  int64_t *n_values_per_frame,
                                                  tng_data_type *type);

/** @brief Get the date and time of initial file creation in ISO format (string).
 *  @param tng_data is a trajectory data container.
 *  @param time is a pointer to the string in which the date will be stored. Memory
   must be reserved beforehand.
 * @return TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
 * has occurred or TNG_CRITICAL (2) if a major error has occured.
 */
tng_function_status tng_time_get_str(tng_trajectory_t tng_data, char *time);

#ifdef __cplusplus
}  /* end extern "C" */
#endif

#endif /* _TNGIO_H */