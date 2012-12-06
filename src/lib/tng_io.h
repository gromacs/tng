#ifndef _TNGIO_H
#define _TNGIO_H     1

#include <stdio.h>
#include <inttypes.h>

#define TNG_VERSION 1

#define TNG_PARTICLE_DEPENDENT 1
#define TNG_FRAME_DEPENDENT 2

// #define TNG_MAX_BLOCK_PARTICLES 1000
#define TNG_MAX_DATE_STR_LEN 24
#define TNG_HASH_LEN 16
#define TNG_MAX_STR_LEN 1024

typedef enum {TNG_BIG_ENDIAN_32,
              TNG_LITTLE_ENDIAN_32,
              TNG_BYTE_PAIR_SWAP_32} tng_endianness_32;

typedef enum {TNG_BIG_ENDIAN_64,
              TNG_LITTLE_ENDIAN_64,
              TNG_QUAD_SWAP_64,
              TNG_BYTE_PAIR_SWAP_64,
              TNG_BYTE_SWAP_64} tng_endianness_64;


#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a < _b ? _a : _b; })
     
#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a > _b ? _a : _b; })

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

typedef enum {TNG_CONSTANT_N_ATOMS, TNG_VARIABLE_N_ATOMS} tng_variable_n_atoms_flag;

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
    int64_t from_atom_id;      /* One of the atoms of the bond */
    int64_t to_atom_id;        /* The other atom of the bond */
};

struct tng_atom {
    tng_residue_t residue;  /* The molecule containing this atom */
    int64_t id;                /* A unique (per molecule) ID number of the atom */
    char *atom_type;            /* The atom_type (depending on the forcefield) */
    char *name;                 /* The name of the atom */
};

struct tng_residue {
    tng_chain_t chain; /* The chain containing this residue */
    int64_t id;                /* A unique (per chain) ID number of the residue */
    char *name;                 /* The name of the residue */
    int64_t n_atoms;           /* The number of atoms in the residue */
    tng_atom_t atoms;     /* A list of atoms in the residue */
};

struct tng_chain {
    tng_molecule_t molecule;  /* The molecule containing this chain */
    int64_t id;                    /* A unique (per molecule) ID number of the chain */
    char *name;                     /* The name of the chain */
    int64_t n_residues;            /* The number of residues in the chain */
    tng_residue_t residues;   /* A list of residues in the chain */
};

struct tng_molecule {
    int64_t id;                    /* A unique ID number of the molecule */
    int64_t quaternary_str;        /* Quaternary structure of the molecule.
                                       1 => monomeric
                                       2 => dimeric
                                       3 => trimeric
                                       etc */
    int64_t n_chains;              /* The number of chains in the molecule */
    int64_t n_residues;            /* The number of residues in the molecule */
    int64_t n_atoms;               /* The number of atoms in the molecule */
    int64_t n_bonds;               /* The number of bonds in the molecule.
                                       If the bonds are not specified this
                                       value can be 0. */
    char *name;                     /* The name of the molecule */
    tng_chain_t chains;       /* A list of chains in the molecule */
    tng_residue_t residues;   /* A list of residues in the molecule */
    tng_atom_t atoms;         /* A list of the atoms in the molecule */
    tng_bond_t bonds;         /* A list of the bonds in the molecule */
};

struct tng_gen_block {
    int64_t header_contents_size;      /* The size of the block header in bytes */
    int64_t block_contents_size;       /* The size of the block contents in bytes */
    int64_t id;                        /* The ID of the block to determine its type */
    char hash[TNG_HASH_LEN];            /* The MD5 hash of the block to verify integrity */
    char *name;                         /* The name of the block */
    int64_t block_version;             /* The library version used to write the block */
    char *header_contents;              /* The full block header contents */
    char *block_contents;               /* The full block contents */
};

struct tng_frame_set_toc {
    int64_t n_blocks;                  /* The number of blocks listed in this
                                           table of contents */
    char **block_names;             /* A list of block names */
};

struct tng_particle_mapping {
    int64_t num_first_particle;        /* The index number of the first particle
                                           in this mapping block */
    int64_t n_particles;               /* The number of particles list in this
                                           mapping block */
    int64_t *real_particle_numbers;    /* the mapping of index numbers to the
                                           real particle numbers in the
                                           trajectory */
};

struct tng_trajectory_frame_set {
    struct tng_frame_set_toc contents;  /* The table of contents of this frame set */
    int64_t n_mapping_blocks;          /* The number of different particle
                                           mapping blocks present. */
    struct tng_particle_mapping *mappings; /* The atom mappings of this frame set */
    int64_t first_frame;               /* The first frame of this frame set */
    int64_t n_frames;                  /* The number of frames in this frame set */
    int64_t *molecule_cnt_list;        /* A list of the number of each molecule
                                           type - only used when using variable
                                           number of atoms */
    int64_t n_particles;               /* The number of particles/atoms - only
                                           used when using variable number of
                                           atoms */
    int64_t next_frame_set_file_pos;   /* The file position of the next frame set */
    int64_t prev_frame_set_file_pos;   /* The file position of the previous
                                           frame set */
    int64_t long_stride_next_frame_set_file_pos; /* The file position of the frame
                                                     set one stride step ahead */
    int64_t long_stride_prev_frame_set_file_pos; /* The file position of the frame
                                                     set one stride step behind */

    /* The data blocks in a frame set are trajectory data blocks */
    int n_particle_data_blocks;            /* The number of data blocks
                                              of particle dependent data */
    struct tng_particle_data *tr_particle_data; /* A list of data blocks
                                                containing particle
                                                dependent data */
    int n_data_blocks;
    struct tng_data *tr_data;
};

/* Data can be either double, float, int or a string */
union data_values {
    double d;
    float f;
    int i;
    char *c;
};

/* FIXME: Should there be a pointer to a tng_gen_block from each data block? */
struct tng_particle_data {
    int64_t block_id;                  /* The block ID of the data block
                                           containing this particle data.
                                           This is used to determine the
                                           kind of data that is stored */
    char *block_name;                   /* The name of the data block.
                                           This is used to determine the
                                           kind of data that is stored */
    tng_data_type datatype;             /* The type of data stored. */
    int64_t first_frame_with_data;     /* The first frame number of the
                                           first data point */
    int64_t n_frames;                  /* The number of frames in this
                                           frame set */
    int64_t n_values_per_frame;        /* The number of values stored per
                                           frame */
    int64_t stride_length;             /* The number of frames between
                                           each data point - e.g. when
                                           storing sparse data. */
    int64_t codec_id;                  /* ID of the CODEC used for compression
                                           0 == no compression. */
    double compression_multiplier;      /* The multiplier used for getting
                                           integer values for compression */
    union data_values ***values;        /* A 3-dimensional array of values,
                                           sized n_frames x n_particles x
                                           n_values_per_frame */
};

struct tng_data {
    int64_t block_id;
    char *block_name;                   /* The name of the data block.
                                           This is used to determine the
                                           kind of data that is stored */
    tng_data_type datatype;             /* The type of data stored. */
    int64_t first_frame_with_data;
    int64_t n_frames;
    int64_t n_values_per_frame;
    int64_t stride_length;
    int64_t codec_id;                  /* ID of the CODEC used for compression
                                           0 == no compression. */
    double compression_multiplier;
    union data_values **values;         /* A 2-dimensional array of values,
                                           sized n_frames x n_values_per_frame */
};



struct tng_trajectory {
    char *input_file_path;                /* The path of the input trajectory file */
    FILE *input_file;                     /* A handle to the input file */
    long int input_file_pos;              /* The reading position of the file */
    long int input_file_len;              /* The length of the input file */
    char *output_file_path;               /* The path of the output trajectory file */
    FILE *output_file;                    /* A handle to the output file */
    long int output_file_pos;             /* The writing position of the file */
    
    tng_endianness_32 endianness_32;    /* The endianness of 32 bit values of
                                           the current computer */
    tng_endianness_64 endianness_64;    /* The endianness of 64 bit values of
                                           the current computer */
    
    char *program_name;                 /* The name of the program producing
                                           this trajectory */
    char *forcefield_name;              /* The forcefield used in the simulations */
    char *user_name;                    /* The name of the user running the
                                           simulations */
    int64_t time;                      /* The time (n seconds since 1970) when
                                           the file was created */
    char *computer_name;                /* The name of the computer on which the
                                           simulations were performed */
    char *pgp_signature;                /* The PGP signature of the user creating
                                           the file. */
    char var_num_atoms_flag;            /* A flag indicating if the number of atoms
                                           can vary throughout the simulation, e.g.
                                           using a grand canonical ensemble */
    int64_t frame_set_n_frames;        /* The number of frames in a frame set.
                                           It is allowed to have frame sets with
                                           fewer frames, but this will help searching
                                           for specific frames */
    int64_t stride_length;             /* The number of frame sets in a long stride
                                           step */

    int64_t n_molecules;               /* The number of different kinds of
                                           molecules in the trajectory */
    tng_molecule_t molecules;          /* A list of molecules in the trajectory */
    int64_t *molecule_cnt_list;        /* A list of the count of each molecule -
                                           if using variable number of particles
                                           this will be specified in each frame set */
    int64_t n_particles;               /* The total number of particles/atoms. If
                                           using variable number of particles this
                                           will be specified in each frame set */
    
    int64_t n_id_name_pairs;                           /* The number of ID-name
                                                           pairs */
    struct tng_block_id_name_pair *id_name_pairs;       /* A list of ID-name
                                                           pairs */

    int64_t first_trajectory_frame_set_input_file_pos;   /* The pos in the src
                                                           file of the first
                                                           frame set */
    int64_t first_trajectory_frame_set_output_file_pos;  /* The pos in the dest
                                                           file of the first
                                                           frame set */
    int64_t last_trajectory_frame_set_input_file_pos;    /* The pos in the src
                                                           file of the last
                                                           frame set */
    int64_t last_trajectory_frame_set_output_file_pos;   /* The pos in the dest
                                                           file of the last
                                                           frame set */
    struct tng_trajectory_frame_set current_trajectory_frame_set; /* The currently
                                                                     active frame
                                                                     set */
    long int current_trajectory_frame_set_input_file_pos; /* The pos in the src
                                                           file of the current
                                                           frame set */
    long int current_trajectory_frame_set_output_file_pos;/* The pos in the dest
                                                           file of the current
                                                           frame set */
    int64_t n_trajectory_frame_sets;                   /* The number of frame sets
                                                           in the trajectory */
    
    int64_t n_trajectory_blocks;       /* The number of trajectory blocks
                                           in the file */
    int64_t n_non_trajectory_blocks;   /* The number of non-trajectory blocks
                                           in the file */
    struct tng_gen_block non_trajectory_blocks[32]; /* A list of non-trajectory
                                                       blocks */
                                                       
    /* These data blocks are non-trajectory data blocks */
    int n_particle_data_blocks;                     /* The number of non-frame dependent
                                                       particle dependent data blocks */
    struct tng_particle_data *non_tr_particle_data; /* A list of data blocks containing particle
                                                       dependent data */
    
    int n_data_blocks;                  /* The number of frame and particle independent
                                           data blocks */
    struct tng_data *non_tr_data;       /* A list of frame and particle indepdendent data blocks */
};



/* Setup a trajectory data container.
   *tng_data is a pointer to pre-allocated memory.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_trajectory_init(tng_trajectory_t tng_data);

/* Clean up a trajectory data container.
   *tng_data is a pointer to pre-allocated memory containing trajectory data.
   All allocated memory in the data structure is freed.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_trajectory_destroy(tng_trajectory_t tng_data);

/* Set the name of the input file.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_input_file_set(tng_trajectory_t tng_data,
                                       const char *file_name);

/* Set the name of the output file.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_output_file_set(tng_trajectory_t tng_data,
                                        const char *file_name);

/* Set the program name used when creating the trajectory.
   *tng_data is a pointer to pre-allocated memory containing trajectory data.
   *new_name is a pointer to the string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_program_name_set(tng_trajectory_t tng_data,
                                         const char *new_name);

/* Set the name of the forcefield used in the trajectory.
   *tng_data is a pointer to pre-allocated memory containing trajectory data.
   *new_name is a pointer to the string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_forcefield_name_set(tng_trajectory_t tng_data,
                                            const char *new_name);

/* Set the name of the user creating the trajectory.
   *tng_data is a pointer to pre-allocated memory containing trajectory data.
   *new_name is a pointer to the string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_user_name_set(tng_trajectory_t tng_data,
                                      const char *new_name);

/* Set the name of the computer used when creating the trajectory.
   *tng_data is a pointer to pre-allocated memory containing trajectory data.
   *new_name is a pointer to the string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_computer_name_set(tng_trajectory_t tng_data,
                                          const char *new_name);

/* Set the pgp_signature of the trajectory.
   *tng_data is a pointer to pre-allocated memory containing trajectory data.
   *signature is a pointer to the string containing the pgp_signature.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_signature_set(tng_trajectory_t tng_data,
                                      const char *signature);

/* Setup a molecule container.
   *molecule is a pointer to pre-allocated memory.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_molecule_init(tng_molecule_t molecule);

/* Clean up a molecule container.
   *molecule is a pointer to pre-allocated memory containing a molecule.
   All allocated memory in the data structure is freed.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_molecule_destroy(tng_molecule_t molecule);

/* Setup a data block.
   *block is a pointer to pre-allocated memory.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_block_init(struct tng_gen_block *block);

/* Clean up a data block.
   *block is a pointer to pre-allocated memory.
   All allocated memory in the data structure is freed.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_block_destroy(struct tng_gen_block *block);

/* Set the name of a data block.
   tng_data is the trajectory data container containing the block..
   *block is a pointer to the block to rename.
   *new_name is a string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_block_name_set(tng_trajectory_t tng_data,
                                       struct tng_gen_block *block,
                                       const char *new_name);

/* Add a molecule to the trajectory.
   tng_data is the trajectory data container containing the block..
   *name is a pointer to the string containing the name of the new molecule.
   **molecule is a pointer to a pointer to the newly created molecule.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_molecule_add(tng_trajectory_t tng_data,
                                     const char *name,
                                     tng_molecule_t *molecule);

/* Set the name of a molecule.
   tng_data is the trajectory data container containing the molecule..
   *molecule is a pointer to the molecule to rename.
   *new_name is a string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_molecule_name_set(tng_trajectory_t tng_data,
                                          tng_molecule_t molecule,
                                          const char *new_name);

/* Set the count of a molecule.
   tng_data is the trajectory data container containing the molecule..
   *molecule is a pointer to the molecule to rename.
   *cnt is a pointer to the variable to be populated with the count.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_molecule_cnt_get(tng_trajectory_t tng_data,
                                         tng_molecule_t molecule,
                                         int64_t *cnt);

/* Set the count of a molecule.
   tng_data is the trajectory data container containing the molecule..
   *molecule is a pointer to the molecule to rename.
   cnt is the number of instances of this molecule.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_molecule_cnt_set(tng_trajectory_t tng_data,
                                         tng_molecule_t molecule,
                                         int64_t cnt);

/* Add a chain to a molecule.
   tng_data is the trajectory data container containing the molecule..
   *molecule is a pointer to the molecule to add a chain to.
   *name is a pointer to a string containing the name of the chain.
   **chain is a pointer to a pointer to the newly created chain.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_molecule_chain_add(tng_trajectory_t tng_data,
                                              tng_molecule_t molecule,
                                              const char *name,
                                              tng_chain_t *chain);

/* Set the name of a chain.
   tng_data is the trajectory data container containing the atom..
   *chain is a pointer to the chain to rename.
   *new_name is a string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_chain_name_set(tng_trajectory_t tng_data,
                                       tng_chain_t chain,
                                       const char *new_name);

/* Add a residue to a chain.
   tng_data is the trajectory data container containing the chain..
   *chain is a pointer to the chain to add a residue to.
   *name is a pointer to a string containing the name of the residue.
   **residue is a pointer to a pointer to the newly created residue.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_chain_residue_add(tng_trajectory_t tng_data,
                                             tng_chain_t chain,
                                             const char *name,
                                             tng_residue_t *residue);

/* Set the name of a residue.
   tng_data is the trajectory data container containing the residue.due.
   *residue is a pointer to the residue to rename.
   *new_name is a string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_residue_name_set(tng_trajectory_t tng_data,
                                         tng_residue_t residue,
                                         const char *new_name);

/* Add an atom to a residue.
   *tng_data is a pointer to the trajectory containing the residue.
   *residue is a pointer to the residue to add an atom to.
   *atom_name is a pointer to a string containing the name of the atom.
   *atom_type is a pointer to a string containing the atom type of the atom.
   **atom is a pointer to a pointer to the newly created atom.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_residue_atom_add(tng_trajectory_t tng_data,
                                            tng_residue_t residue,
                                            const char *atom_name,
                                            const char *atom_type,
                                            tng_atom_t *atom);

/* Set the name of an atom.
   tng_data is the trajectory data container containing the atom..
   *atom is a pointer to the atom to rename.
   *new_name is a string containing the wanted name.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_atom_name_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_name);

/* Set the atom type of an atom.
   tng_data is the trajectory data container containing the atom..
   *atom is a pointer to the atom to change.
   *new_name is a string containing the atom type.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_atom_type_set(tng_trajectory_t tng_data,
                                      tng_atom_t atom,
                                      const char *new_type);

/* Read the header blocks from the input_file of tng_data.
   The trajectory blocks must be read separately and iteratively in chunks
   to fit in memory.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
   compared to the md5 hash of the read contents to ensure valid data.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_file_headers_read(tng_trajectory_t tng_data,
                                          const tng_hash_mode hash_mode);

/* Write the header blocks to the output_file of tng_data.
   The trajectory blocks must be written separately and iteratively in chunks
   to fit in memory.
   tng_data is a trajectory data container. tng_data->output_file_path specifies which
   file to write to. If the file (output_file) is not open it will be opened.
   If hash_mode == TNG_USE_HASH an md5 hash for each header block will be generated.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_file_headers_write(tng_trajectory_t tng_data,
                                           const tng_hash_mode hash_mode);


/* Read one (the next) block (of any kind) from the input_file of tng_data.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   *block_data is a pointer to the struct which will be populated with the
   data. If block_data->input_file_pos > 0 it is the position from where the reading
   starts otherwise it starts from the current position.
   If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
   compared to the md5 hash of the read contents to ensure valid data.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor
   error has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_block_read_next(tng_trajectory_t tng_data,
                                        struct tng_gen_block *block_data,
                                        const tng_hash_mode hash_mode);


/* Read one (the next) frame set, including toc, mapping and related data blocks
   from the input_file of tng_data.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   If hash_mode == TNG_USE_HASH the written md5 hash in the file will be
   compared to the md5 hash of the read contents to ensure valid data.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_frame_set_read_next(tng_trajectory_t tng_data,
                                            const tng_hash_mode hash_mode);

/* Write one frame set, including toc, mapping and related data blocks
   to the output_file of tng_data.
   tng_data is a trajectory data container. tng_data->output_file_path specifies which
   file to write to. If the file (output_file) is not open it will be opened.
   If hash_mode == TNG_USE_HASH an md5 hash for each block will be generated.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_frame_set_write(tng_trajectory_t tng_data,
                                        const tng_hash_mode hash_mode);

tng_function_status tng_frame_set_new(tng_trajectory_t tng_data,
                                      const int64_t first_frame,
                                      const int64_t n_frames);

tng_function_status tng_particle_data_block_add(tng_trajectory_t tng_data,
                                        const int64_t id,
                                        const char *block_name,
                                        const char datatype,
                                        const int64_t n_frames,
                                        const int64_t n_values_per_frame,
                                        const int64_t stride_length,
                                        const int64_t first_particle_number,
                                        const int64_t n_particles,
                                        int64_t codec_id,
                                        void *new_data);


/* Read one trajectory block from the input_file of tng_data.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   block_id is the ID of the block to read.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_traj_block_read(tng_trajectory_t tng_data,
                                        int64_t block_id);
        
/* Write one trajectory block to the output_file of tng_data.
   tng_data is a trajectory data container. tng_data->output_file_path specifies which
   file to write to. If the file (output_file) is not open it will be opened.
   block_id is the ID of the block to write.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_traj_block_write(tng_trajectory_t tng_data,
                                         int64_t block_id);

/* Read a requested frame set.
   tng_data is a trajectory data container. tng_data->current_trajectory_frame_set
   will be the read frame set.
   frame_set_nr is the number of the frame set to return (starting from 0).
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_frame_set_read_nr(tng_trajectory_t tng_data,
                                          int64_t frame_set_nr);

/* Read a number of consecutive trajectory frames from the input_file of tng_data.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   start_frame_nr is the index number of the first frame to read.
   end_frame_nr is the index number of the last frame to read.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_frame_read_interval(tng_trajectory_t tng_data,
                                       int64_t start_frame_nr,
                                       int64_t end_frame_nr);

/* Write a number of consecutive trajectory frames to the output_file of tng_data.
   tng_data is a trajectory data container. tng_data->output_file_path specifies which
   file to write to. If the file (output_file) is not open it will be opened.
   start_frame_nr is the index number of the first frame to write.
   end_frame_nr is the index number of the last frame to write.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_frame_write_interval(tng_trajectory_t tng_data,
                                        int64_t start_frame_nr,
                                        int64_t end_frame_nr);

/* Read and retrieve non-particle data from the last read frame set.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   block_id is the id number of the data block to read.
   ***values is a pointer to a 2-dimensional array (memory unallocated), which will
   be filled with data. The array will be sized (n_frames * n_values_per_frame).
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_data_get(tng_trajectory_t tng_data,
                                 int64_t block_id,
                                 union data_values ***values);

/* Read and retrieve non-particle data, in a specific interval, from the trajectory.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   block_id is the id number of the data block to read.
   start_frame_nr is the index number of the first frame to read.
   end_frame_nr is the index number of the last frame to read.
   ***values is a pointer to a 2-dimensional array (memory unallocated), which will
   be filled with data. The array will be sized (n_frames * n_values_per_frame).
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_data_interval_get(tng_trajectory_t tng_data,
                                          int64_t block_id,
                                          int64_t start_frame_nr,
                                          int64_t end_frame_nr,
                                          union data_values ***values);

/* Read and retrieve particle data, in a specific interval, from the last read frame set.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   block_id is the id number of the particle data block to read.
   ****values is a pointer to a 3-dimensional array (memory unallocated), which will
   be filled with data. The array will be sized (n_frames * n_particles * n_values_per_frame).
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_particle_data_get(tng_trajectory_t tng_data,
                                          int64_t block_id,
                                          union data_values ****values);

/* Read and retrieve particle data, in a specific interval, from the trajectory.
   tng_data is a trajectory data container. tng_data->input_file_path specifies which
   file to read from. If the file (input_file) is not open it will be opened.
   block_id is the id number of the particle data block to read.
   start_frame_nr is the index number of the first frame to read.
   end_frame_nr is the index number of the last frame to read.
   first_particle_number is the number of the first particle, for which to retrive data.
   last_particle_number is the number of the last particle, for which to retrieve data.
   ****values is a pointer to a 3-dimensional array (memory unallocated), which will
   be filled with data. The array will be sized (n_frames * n_particles * n_values_per_frame).
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_particle_data_interval_get(tng_trajectory_t tng_data,
                                                  int64_t block_id,
                                                  int64_t start_frame_nr,
                                                  int64_t end_frame_nr,
                                                  int64_t first_particle_number,
                                                  int64_t last_particle_number,
                                                  union data_values ****values);

/* Get the date and time of initial file creation in ISO format (string).
   tng_data is a trajectory data container.
   *time is a pointer to the string in which the date will be stored. Memory
   must be reserved beforehand.
   Returns TNG_SUCCESS (0) if successful, TNG_FAILURE (1) if a minor error
   has occurred or TNG_CRITICAL (2) if a major error has occured. */
tng_function_status tng_time_get_str(tng_trajectory_t tng_data, char *time);

#ifdef __cplusplus
}  /* end extern "C" */
#endif

#endif /* _TNGIO_H */