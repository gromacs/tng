Trajectory format specification

1.  Release notes
    =============

General notes:
--------------

\* Currently the API is lacking many data retrieval and getter
functions. This will be fixed soon.

\* It might be problematic searching for a specific time if frames are
not consecutive. This problem is somewhat hypothetical, but might be
good to point out.

\* Currently we cannot track a specific molecule in a grand canonical
ensemble.

\* We should request comments on the specs from other groups - after v1.

Erik's comments:

\* Keep the number of calls small

\* Use a prefix for the calls, e.g. tng\_open, tng\_close, etc.; use
e.g. adv\_ prefix for advanced API

\* Help routine to return the full list of atom types

\* In general info:

- number of stride pointers

- for each pointer, the number of frame sets that it skips

=\> have 3 stride pointers with defaults 1, 100 and 10,000 frame sets

\* include zlib as a supported compressor

To do:
------

\* Make a drawing

\* Include API function for signing of the trajectory

\* Sanity checklist:

1.  Block header size
2.  Block contents size
3.  Compare hashes (currently doesn’t abort only prints warning)

1.  maybe add a flag to choose between abort/warning
2.  empty hash is always accepted

4.  

Version 0.98
------------

\* Use signed int instead of unsigned.

\* Changed names from ‘trg’ to ‘tng’.

Version 0.97
------------

\* Added chains and residues to the molecule description.

Version 0.96
------------

\*Renamed “Trajectory Box Shape”, “Trajectory Positions” etc.

\*Renamed “Trajectory Index Block” to “Particle Mapping Block”

Version 0.95
------------

\* Changed name of ‘Atom name block’ to ‘Molecule block’ and included
connectivity in that block. The number of each molecule is also
specified there, unless using variable number of particles (grand
canonical ensemble).

\* Removed ‘variable number of particles’ and ‘variable number of
values’ flags from data blocks.

\* Updated the API and headers

Version 0.9
-----------

\* Removed endianness/string block

\* Added PGP signature to General Info block

\* Removed reserved ‘user data blocks’. Those are already handled as
normal ‘data blocks’. Reserved are only box, positions, velocities and
forces.

Version 0.8

\* Use only MD5 hashes

\* Changed time (in general info block) to 64 bit int.

Version 0.7
-----------

\* New hierarchy rules

\* Changed "Trajectory Info" to "General Info" block

\* Moved "BOX SHAPE" before the index and trajectory frame blocks

\* general, int., bin32, bin64 user data is explicitly not frame
dependent

\* Variable number of atoms will be supported in version 1

\* Made more clear the nesting in the "trajectory group blocks"

\* Allow header only files, i.e. no trajectory blocks

\* Removed number of frames from the "index block" description

\* Removed ASCII recommendation from the string description, i.e. [1]

\* Added molecule ID to “atom name block”

Version 0.6
-----------

Changed endianness block to endianness and string length block.

Changed order of hash specifications in the header block.

Added initial API specifications

Version 0.5
-----------

Included an additional header block that specifies the endianness. It
comes before all other blocks.

2.  Specifications
    ==============

General specifications are given below. Others are described in the
relevant block sections.

1.  The file contains a number of blocks.
2.  The order of the blocks follows the order specified in section 4
    “Description of blocks”
3.  All integers and floating point (floats and doubles) values are
    stored using big endian byte ordering. Conversions to and from the
    native format of the computer is performed during reading and
    writing of numerical fields.
4.  MD5 hashes are used to verify the integrity of the data.
5.  Strings are limited to a max length of 1023 characters and are
    terminated by a null character (‘\\000’). If longer text data must
    be saved a data block containing multiple entries of general
    (character/string) data can be used.
6.  If a trajectory converter program encounters errors during reading a
    block which format is not recognized, the block is to be written out
    as binary object without modification.
7.  Each group of blocks to contain a Table of Contents which lists the
    included blocks within the group.^[[a]](#cmnt1)^
8.  No compression (in ver.1) for general data streams such as integers,
    floating point numbers, particle indices etc.

3.  Each block contains the following fields (header):
    ==================================================

1.  64 bit size of the header
2.  64 bit size of the block contents (except header)
3.  64 bit block type identifier
4.  16 characters MD5 Hash (or 16 “\\0” characters)
5.  name[1]
6.  64 bit version of the block ^[[b]](#cmnt2)^(allows addition of more
    fields in the future to existing blocks, although old fields should
    never be removed, to allow older readers read new files)

4.  Description of blocks (each with a unique 64 bit identifier and a matching "name"):
    ===================================================================================

1.  info block (1) "GENERAL INFO" (required)
2.  molecules block (2) "MOLECULES" (optional)
3.  trajectory ids and names (4) “TRAJECTORY IDS AND NAMES” (optional)
4.  trajectory frames, box shape block (10000) "BOX SHAPE" (optional,
    can be present before the frame sets if it does not change or inside
    the frame sets if it varies)
5.  trajectory frame set block (5) "TRAJECTORY FRAME SET" (required)
    (multiple “trajectory frame sets” are allowed)

1.  trajectory table of contents block (6) “BLOCK TABLE OF CONTENTS”
    (required)
2.  trajectory frames, box shape block (10000) "BOX SHAPE" (optional,
    can be present before the frame sets if it does not change or inside
    the frame sets if it varies)
3.  trajectory particle mapping block (7) "PARTICLE MAPPING" (required
    if there are trajectory frames blocks) (multiple particle mapping
    blocks with corresponding trajectory frames blocks are allowed, e.g.
    to allow parallel writes of different atom sets)

1.  trajectory frames, positions, block (10001) "POSITIONS" (optional)
2.  trajectory frames, velocities, block (10002) "VELOCITIES" (optional)
3.  trajectory frames, forces, block (10003) "FORCES" (optional)

6.  ...other specified blocks, both non-trajectory and trajectory
    blocks, each with unique id & name

Data blocks can be used to store whatever data is needed. Data blocks
with IDs in the range 10000 to 10999 are reserved for standard data
(such as box shape, positions, velocities etc.), whereas IDs from 11000
and above can be used for any kind of user data.

5.  NOTES ABOUT BLOCK KINDS:
    ========================

There can be only one block at the beginning of the file in the standard
case with fixed charges throughout the simulation (that’s the case for
version 1 of the format). For simulations where charges vary each frame
will include a block with the values.

The trajectory frame blocks must be collected into frame sets, each

such frame set has as its first block the "trajectory frame set block".
Each frame set will contain (multiple) particle mapping blocks,
positions, velocities, forces etc.^[[c]](#cmnt3)^

6.  Requirements on block order:
    ============================

1.  The order follows section 4. “Description of blocks” with the
    corresponding nesting of multiple “trajectory frame sets” and
    “particle mapping blocks”.
2.  All non-trajectory frame blocks (e.g. user ones) must appear before
    the trajectory blocks
3.  Trajectory particle mapping blocks are optional. If they are
    present, they must appear before the corresponding trajectory data
    blocks. If there are multiple trajectory data blocks, the
    corresponding particle mapping blocks come right before them. I.e.
    ParticleMappingBlock1-\>DataBlock1-\>P.MappingBlock2-\>DataBlock2.
4.  Blocks within a group of blocks are ordered by their ID.

7.  Other requirements:
    ===================

1.  Most blocks are optional except for the “general info” blocks
2.  No limit on the number of times that trajectory related blocks are
    allowed to appear

8.  Specification of the block contents (all blocks have the same header as described above) for version 1 of each block type.
    ==========================================================================================================================

BLOCK: general info block
-------------------------

1.  name and version of the program used to perform the simulation (upon
    file creation)[1]
2.  name and version of the program used when finishing the file[1]
3.  name of the force field used to perform the simulation
    [1]^[[d]](#cmnt4)^
4.  name of the person who created the file [1]
5.  name of the person who last modified the file [1]
6.  64 bit time of initial file creation, seconds since 1970
7.  64 bit time of completing the simulation, seconds since 1970
8.  name of computer/other info where the file was created [1]
9.  name of computer/other info where the file was completed [1]
10. PGP signature (optional and 0 terminated string)
11. 8 bit flag Use variable number of atoms.
12. 64 bit number of frames in each frame set (this is the expected
    number of frames in each set, but it does not have to be constant,
    it is OK to have frame sets with fewer or more frames, e.g. after
    concatenating multiple trajectory files. This avoids the need to
    recompress all data after a concatenation, but it means that
    searching for a specific frame might need a few more steps between
    frame sets.). For simulations using a grand canonical ensemble it is
    best to set this to 1 so that the number of atoms in the frame sets
    can be updated regularly.
13. 64 bit pointer from the beginning of the info block to the beginning
    of the first trajectory frame set [2]
14. 64 bit pointer from the beginning of the info block to the beginning
    of the last trajectory frame set [2] (updated when finishing writing
    the trajectory file - otherwise set to -1)
15. 64 bit length of steps (number of “trajectory frame set blocks”) for
    long stride pointers (default 100 “trajectory frame set blocks”).

BLOCK: molecules block (optional)
---------------------------------

1.  64 bit number of molecules
2.  For each molecule:

1.  64 bit Molecule ID
2.  Molecule name [1]
3.  64 bit quaternary structure, e.g. 1 means monomeric, 4 means
    tetrameric etc.
4.  64 bit number of molecules of this kind - only if not using
    “variable number of atoms” in the “general info block”.
5.  64 bit number of chains in the molecule
6.  64 bit number of residues in the molecule
7.  64 bit number of atoms in the molecule^[[e]](#cmnt5)^
8.  For each chain:

1.  64 bit Chain ID (unique in molecule)
2.  Chain name [1]
3.  64 bit number of residues in the chain
4.  For each residue:

1.  64 bit Residue ID (unique in the chain)
2.  Residue name [1]
3.  64 bit number of atoms in the residue
4.  For each atom:

1.  64 bit Atom ID (unique in the molecule)
2.  Atom name [1]
3.  Atom type [1]

9.  64 bit number of bonds in the molecule
10. For each bond:

5.  64 bit integer From Atom ID.
6.  64 bit integer To Atom ID.

BLOCK: trajectory frame set block
---------------------------------

1.  64 bit number of first frame (zero based numbering)
2.  64 bit number of frames (NF)
3.  Array of 64 bit integers specifying the count of each molecule type.
    The molecule types are listed in the “Atom names block” and should
    be listed in the same order here. This should only be present when
    the variable number of atoms flag in the “General info block” is set
    to TRUE. This is used for e.g. simulations using a grand canonical
    ensemble (in which case the number of frames in each frame set
    should be 1).
4.  64 bit pointer to the next “trajectory frame set block”.
5.  64 bit pointer the previous “trajectory frame set block”.
6.  64 bit long stride pointer to the next e.g. 100th “trajectory frame
    set block”. (Stride length specified in “general info” block.)
7.  64 bit long stride pointer to the previous e.g. 100th “trajectory
    frame set block”.

BLOCK: trajectory table of contents
-----------------------------------

1.  64 bit number of blocks

Contains a listing of all data blocks \_present\_ in the frame set. It
is possible to have multiple blocks with the same ID, but the ID is only
listed once in the “trajectory table of contents” block.

It includes for each block type:

1.  Block name [1]

BLOCK: data blocks
------------------

Frame dependent data blocks should come after the frame set block to
which it belongs. Frame and particle dependent data blocks should come
after the relevant particle mapping block (if using any particle mapping
block).

1.  Char data type flag. 0 = character/string data, 1 = 64 bit integer
    data, 2 = float data (32 bit), 3 = double data (64 bit)
2.  Char dependency flag. 1 = frame dependent, 2 = particle dependent.
    Can be combined, i.e. 3 = frame and particle dependent.
3.  Char sparse data flag to signify if not all frames in the frame sets
    have data entries in this data block, e.g. energies and positions
    might be saved at different intervals meaning that at least one of
    them would be saved as sparse data. Only present if the data is
    frame dependent.
4.  64 bit number of values.
5.  64 bit id of the CODEC used to store the positions
6.  Double (64 bit) multiplier for integers to obtain the appropriate
    floating point number, for compressed frames [3] [\*\*] (only
    present if the above CODEC id is \> 0 and if the data type is double
    or float)

If using sparse data the following fields are required:

1.  64 bit number of first frame containing data.
2.  64 bit number of frames between data points

Particle dependent data blocks contain the following fields:

1.  64 bit number of first particle as stored in the trajectory, zero
    based numbering) (J), this must be the same as in the preceding
    trajectory particle mapping block, if present.

1.  64 bit number of particles in block, this must be the same as in the
    preceding trajectory particle mapping block, if present.

Example 1:

Box shape block (10000) in a frame set with frames 0-99:

1.  Data type: 3 (double)
2.  Dependency: 1 (frame dependent)
3.  Sparse data: 1
4.  Number of values: 9
5.  Codec ID: 0
6.  First frame containing data: 0
7.  Number of frames between data points: 50
8.  For each frame (2 frames with data in this block):

1.  9 double (64 bit) values describing the shape of the block

Example 2:

Positions block (ID 10001) in a frame set with frames 1000-1099:

1.  Data type: 2 (float)
2.  Dependency: 3 (frame and particle dependent)
3.  Sparse data: 1
4.  Number of values: 3 (x, y and z)
5.  Coded ID: 0
6.  First frame containing data: 100
7.  Number of frames between data points: 10
8.  Number of first particle: 0
9.  Number of particles in block: 1000
10. For each frame (10 frames with data in this block):

1.  For each particle (1000 particles):

1.  32 bit float x coordinate
2.  32 bit float y coordinate
3.  32 bit float z coordinate

Example 3:

Forces block (ID 10003) in a frame set with frames 0-99:

11. Data type: 2 (float)
12. Dependency: 3 (frame and particle dependent)
13. Sparse data: 0
14. Number of values: 3 (x, y and z)
15. Coded ID: 0
16. Number of first particle: 0
17. Number of particles in block: 100
18. For each frame (100 frames with data in this block):

2.  For each particle (100 particles):

1.  32 bit float x coordinate
2.  32 bit float y coordinate
3.  32 bit float z coordinate

BLOCK: particle mapping block
-----------------------------

1.  64 bit number of first particle (particle number as stored in the
    trajectory, zero based numbering) (J)
2.  64 bit number of particles in this particle mapping block (M)
3.  64 bit array of particle numbers^[[f]](#cmnt6)^ (M values):

1.  Each value is the number of the real particle corresponding to the
    particle number as stored in the trajectory.

Should no particle mapping block be present, the mapping is the number

of the real particle == the particle number as stored in the trajectory.

It is possible to have several trajectory/velocities etc. frame blocks
within a frame set, e.g. when faster parallel writes or memory
considerations are needed. In that case a separate particle mapping
block is needed for each of the trajectory/velocities etc. blocks.

Relation between trajectory blocks:
===================================

Particle mapping blocks contain the remapping of actual particle index
and the particle index as appearing in the trajectory file. They are
optional. If they are not given, there is no remapping. All trajectory
blocks for the same set of particles must follow each other, i.e.
positions for particle 0-99, then velocities for particle 0-99, then
positions for particle 100-199, then velocities for particle 100-199.
All non-particle trajectory blocks must appear before any particle
containing trajectory blocks.

Limitations on the number of particles in trajectory frame blocks: In
order to be able to read and uncompress data there must be a limit on
the number of particles in each trajectory frame block, therefore most
trajectory frame sets will contain multiple particle mapping / positions
/ velocities / ... blocks. The limit on the number of particles per
trajectory frame blocks should be XXXX. This should be a good value, and
not allowed to be set by the user, since this may prevent reading of the
files on smaller memory machines.

CODEC specifications (id) "name"
================================

1.  uncompressed (0) "UNCOMPRESSED"
2.  XTC positions (1) "XTC"
3.  TNG (2) "TNG"
4.  …

[](#)

[\*\*] Storage of compressed positions / velocities / ...: These are now
all converted to integers before stored. In order to facilitate
recompression without loss of precision it is essential that these are
visible as integers. Therefore the compression blocks all must contain
somewhere a conversion factor from integer to float.

Notes
=====

[1] UTF-8 text string. Make all text strings zero terminated.

[2] 64 bit pointer format: -1UL (all ones), means "not set", which is
what should be written whenever a pointer needs to be written when the
appropriate value is not yet known, while 0 (all zeros), typically means
the end of the list.

[3] Floating point format is big-endian IEEE-754, float (32 bit) or
double (64 bit).

API
===

(The API should be separated into one high- and one low-level API, using
e.g. a tng\_low tag for the low-level functions.)

/\* Setup a trajectory data container.

\*tng\_data is a pointer to pre-allocated memory.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_trajectory\_init(tng\_trajectory\_t
tng\_data);

/\* Clean up a trajectory data container.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

All allocated memory in the data structure is freed.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_trajectory\_destroy(tng\_trajectory\_t
tng\_data);

/\* Set the name of the input file.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_input\_file\_set(tng\_trajectory\_t
tng\_data,

const char \*file\_name);

/\* Set the name of the output file.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_output\_file\_set(tng\_trajectory\_t
tng\_data,

const char \*file\_name);

/\* Set the program name used when creating the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_program\_name\_set(tng\_trajectory\_t
tng\_data,

const char \*new\_name);

/\* Set the name of the forcefield used in the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_forcefield\_name\_set(tng\_trajectory\_t
tng\_data,

const char \*new\_name);

/\* Set the name of the user creating the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_user\_name\_set(tng\_trajectory\_t tng\_data,

const char \*new\_name);

/\* Set the name of the computer used when creating the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_computer\_name\_set(tng\_trajectory\_t
tng\_data,

const char \*new\_name);

/\* Set the pgp\_signature of the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*signature is a pointer to the string containing the pgp\_signature.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_signature\_set(tng\_trajectory\_t tng\_data,

const char \*signature);

/\* Setup a molecule container.

\*molecule is a pointer to pre-allocated memory.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_init(tng\_molecule\_t molecule);

/\* Clean up a molecule container.

\*molecule is a pointer to pre-allocated memory containing a molecule.

All allocated memory in the data structure is freed.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_destroy(tng\_molecule\_t molecule);

/\* Setup a data block.

\*block is a pointer to pre-allocated memory.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_block\_init(struct tng\_gen\_block \*block);

/\* Clean up a data block.

\*block is a pointer to pre-allocated memory.

All allocated memory in the data structure is freed.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_block\_destroy(struct tng\_gen\_block
\*block);

/\* Set the name of a data block.

tng\_data is the trajectory data container containing the block..

\*block is a pointer to the block to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_block\_name\_set(tng\_trajectory\_t
tng\_data,

struct tng\_gen\_block \*block,

const char \*new\_name);

/\* Add a molecule to the trajectory.

tng\_data is the trajectory data container containing the block..

\*name is a pointer to the string containing the name of the new
molecule.

\*\*molecule is a pointer to a pointer to the newly created molecule.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_add(tng\_trajectory\_t tng\_data,

const char \*name,

tng\_molecule\_t \*molecule);

/\* Set the name of a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_name\_set(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

const char \*new\_name);

/\* Set the count of a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to rename.

\*cnt is a pointer to the variable to be populated with the count.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_cnt\_get(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

int64\_t \*cnt);

/\* Set the count of a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to rename.

cnt is the number of instances of this molecule.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_cnt\_set(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

int64\_t cnt);

/\* Add a chain to a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to add a chain to.

\*name is a pointer to a string containing the name of the chain.

\*\*chain is a pointer to a pointer to the newly created chain.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_chain\_add(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

const char \*name,

tng\_chain\_t \*chain);

/\* Set the name of a chain.

tng\_data is the trajectory data container containing the atom..

\*chain is a pointer to the chain to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_chain\_name\_set(tng\_trajectory\_t
tng\_data,

tng\_chain\_t chain,

const char \*new\_name);

/\* Add a residue to a chain.

tng\_data is the trajectory data container containing the chain..

\*chain is a pointer to the chain to add a residue to.

\*name is a pointer to a string containing the name of the residue.

\*\*residue is a pointer to a pointer to the newly created residue.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_chain\_residue\_add(tng\_trajectory\_t
tng\_data,

tng\_chain\_t chain,

const char \*name,

tng\_residue\_t \*residue);

/\* Set the name of a residue.

tng\_data is the trajectory data container containing the residue.due.

\*residue is a pointer to the residue to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_residue\_name\_set(tng\_trajectory\_t
tng\_data,

tng\_residue\_t residue,

const char \*new\_name);

/\* Add an atom to a residue.

\*tng\_data is a pointer to the trajectory containing the residue.

\*residue is a pointer to the residue to add an atom to.

\*atom\_name is a pointer to a string containing the name of the atom.

\*atom\_type is a pointer to a string containing the atom type of the
atom.

\*\*atom is a pointer to a pointer to the newly created atom.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_residue\_atom\_add(tng\_trajectory\_t
tng\_data,

tng\_residue\_t residue,

const char \*atom\_name,

const char \*atom\_type,

tng\_atom\_t \*atom);

/\* Set the name of an atom.

tng\_data is the trajectory data container containing the atom..

\*atom is a pointer to the atom to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_atom\_name\_set(tng\_trajectory\_t tng\_data,

tng\_atom\_t atom,

const char \*new\_name);

/\* Set the atom type of an atom.

tng\_data is the trajectory data container containing the atom..

\*atom is a pointer to the atom to change.

\*new\_name is a string containing the atom type.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_atom\_type\_set(tng\_trajectory\_t tng\_data,

tng\_atom\_t atom,

const char \*new\_type);

/\* Read the header blocks from the input\_file of tng\_data.

The trajectory blocks must be read separately and iteratively in chunks

to fit in memory.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

If hash\_mode == TNG\_USE\_HASH the written md5 hash in the file will be

compared to the md5 hash of the read contents to ensure valid data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_file\_headers\_read(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

/\* Write the header blocks to the output\_file of tng\_data.

The trajectory blocks must be written separately and iteratively in
chunks

to fit in memory.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path

specifies which file to write to. If the file (output\_file) is not open
it

will be opened.

If hash\_mode == TNG\_USE\_HASH an md5 hash for each header block will
be generated.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_file\_headers\_write(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

/\* Read one (the next) block (of any kind) from the input\_file of
tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

\*block\_data is a pointer to the struct which will be populated with
the

data. If block\_data-\>input\_file\_pos \> 0 it is the position from
where the

reading starts otherwise it starts from the current position.

If hash\_mode == TNG\_USE\_HASH the written md5 hash in the file will be

compared to the md5 hash of the read contents to ensure valid data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor

error has occurred or TNG\_CRITICAL (2) if a major error has occured.
\*/

tng\_function\_status tng\_block\_read\_next(tng\_trajectory\_t
tng\_data,

struct tng\_gen\_block \*block\_data,

const tng\_hash\_mode hash\_mode);

/\* Read one (the next) frame set, including toc, mapping and related
data blocks

from the input\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

If hash\_mode == TNG\_USE\_HASH the written md5 hash in the file will be

compared to the md5 hash of the read contents to ensure valid data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_set\_read\_next(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

/\* Write one frame set, including toc, mapping and related data blocks

to the output\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path
specifies

which file to write to. If the file (output\_file) is not open it will
be

opened.

If hash\_mode == TNG\_USE\_HASH an md5 hash for each block will be
generated.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_set\_write(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

tng\_function\_status tng\_frame\_set\_new(tng\_trajectory\_t tng\_data,

const int64\_t first\_frame,

const int64\_t n\_frames);

tng\_function\_status tng\_particle\_data\_block\_add(tng\_trajectory\_t
tng\_data,

const int64\_t id,

const char \*block\_name,

const char datatype,

const int64\_t n\_frames,

const int64\_t n\_values\_per\_frame,

const int64\_t stride\_length,

const int64\_t first\_particle\_number,

const int64\_t n\_particles,

int64\_t codec\_id,

void \*new\_data);

/\* Read one trajectory block from the input\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the ID of the block to read.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_traj\_block\_read(tng\_trajectory\_t
tng\_data,

int64\_t block\_id);

/\* Write one trajectory block to the output\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path
specifies

which file to write to. If the file (output\_file) is not open it will
be

opened.

block\_id is the ID of the block to write.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_traj\_block\_write(tng\_trajectory\_t
tng\_data,

int64\_t block\_id);

/\* Read a requested frame set.

tng\_data is a trajectory data container.

tng\_data-\>current\_trajectory\_frame\_set will be the read frame set.

frame\_set\_nr is the number of the frame set to return (starting from
0).

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_set\_read\_nr(tng\_trajectory\_t
tng\_data,

int64\_t frame\_set\_nr);

/\* Read a number of consecutive trajectory frames from the input\_file
of

tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

start\_frame\_nr is the index number of the first frame to read.

end\_frame\_nr is the index number of the last frame to read.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_read\_interval(tng\_trajectory\_t
tng\_data,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr);

/\* Write a number of consecutive trajectory frames to the output\_file
of

tng\_data.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path
specifies

which file to write to. If the file (output\_file) is not open it will
be

opened.

start\_frame\_nr is the index number of the first frame to write.

end\_frame\_nr is the index number of the last frame to write.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_write\_interval(tng\_trajectory\_t
tng\_data,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr);

/\* Free data is an array of values (2D).

\*\*values is the array to free and will be set to 0 afterwards.

n\_frames is the number of frames in the data array.

n\_values\_per\_frame is the number of values per frame in the data
array.

type is the data type of the data in the array.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_data\_values\_free(union data\_values
\*\*values,

int64\_t n\_frames,

int64\_t n\_values\_per\_frame,

tng\_data\_type type);

/\* Free data is an array of values (3D).

\*\*\*values is the array to free and will be set to 0 afterwards.

n\_frames is the number of frames in the data array.

n\_particles is the number of particles in the data array.

n\_values\_per\_frame is the number of values per frame in the data
array.

type is the data type of the data in the array.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_particle\_data\_values\_free(union
data\_values \*\*\*values,

int64\_t n\_frames,

int64\_t n\_particles,

int64\_t n\_values\_per\_frame,

tng\_data\_type type);

/\* Retrieve non-particle data from the last read frame set.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the data block to read.

\*\*\*values is a pointer to a 2-dimensional array (memory unallocated),
which

will point to the data of the requested data block. The array will be
sized

(n\_frames \* n\_values\_per\_frame).

Since \*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_frames is set to the number of frames in the data. This is needed
to

properly reach and/or free the data.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_data\_get(tng\_trajectory\_t tng\_data,

int64\_t block\_id,

union data\_values \*\*\*values,

int64\_t \*n\_frames,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Read and retrieve non-particle data, in a specific interval, from
the

trajectory.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the data block to read.

start\_frame\_nr is the index number of the first frame to read.

end\_frame\_nr is the index number of the last frame to read.

\*\*\*values is a pointer to a 2-dimensional array (memory unallocated),
which

will be filled with data. The array will be sized

(n\_frames \* n\_values\_per\_frame).

Since \*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_data\_interval\_get(tng\_trajectory\_t
tng\_data,

int64\_t block\_id,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr,

union data\_values \*\*\*values,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Retrieve particle data, from the last read frame set.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the particle data block to read.

\*\*\*\*values is a pointer to a 3-dimensional array (memory
unallocated), which

will point to the data of the requested data block. The array will be
sized

(n\_frames \* n\_particles \* n\_values\_per\_frame).

Since \*\*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_frames is set to the number of frames in the data. This is needed
to

properly reach and/or free the data.

\*n\_particles is set to the number of particles in the returned data.
This is

needed to properly reach and/or free the data.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_particle\_data\_get(tng\_trajectory\_t
tng\_data,

int64\_t block\_id,

union data\_values \*\*\*\*values,

int64\_t \*n\_frames,

int64\_t \*n\_particles,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Read and retrieve particle data, in a specific interval, from the
trajectory.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the particle data block to read.

start\_frame\_nr is the index number of the first frame to read.

end\_frame\_nr is the index number of the last frame to read.

first\_particle\_number is the number of the first particle, for which
to

retrive data.

last\_particle\_number is the number of the last particle, for which to

retrieve data.

\*\*\*\*values is a pointer to a 3-dimensional array (memory
unallocated), which

will be filled with data. The array will be sized

(n\_frames \* n\_particles \* n\_values\_per\_frame).

Since \*\*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status
tng\_particle\_data\_interval\_get(tng\_trajectory\_t tng\_data,

int64\_t block\_id,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr,

int64\_t first\_particle\_number,

int64\_t last\_particle\_number,

union data\_values \*\*\*\*values,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Get the date and time of initial file creation in ISO format
(string).

tng\_data is a trajectory data container.

\*time is a pointer to the string in which the date will be stored.
Memory

must be reserved beforehand.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_time\_get\_str(tng\_trajectory\_t tng\_data,
char \*time);

Full API Header file:

(WIP)

\#ifndef \_TNGIO\_H

\#define \_TNGIO\_H 1

\#include <stdio.h\>

\#include <inttypes.h\>

\#define TNG\_VERSION 1

\#define TNG\_PARTICLE\_DEPENDENT 1

\#define TNG\_FRAME\_DEPENDENT 2

// \#define TNG\_MAX\_BLOCK\_PARTICLES 1000

\#define TNG\_MAX\_DATE\_STR\_LEN 24

\#define TNG\_HASH\_LEN 16

\#define TNG\_MAX\_STR\_LEN 1024

typedef enum {TNG\_BIG\_ENDIAN\_32,

TNG\_LITTLE\_ENDIAN\_32,

TNG\_BYTE\_PAIR\_SWAP\_32} tng\_endianness\_32;

typedef enum {TNG\_BIG\_ENDIAN\_64,

TNG\_LITTLE\_ENDIAN\_64,

TNG\_QUAD\_SWAP\_64,

TNG\_BYTE\_PAIR\_SWAP\_64,

TNG\_BYTE\_SWAP\_64} tng\_endianness\_64;

\#define min(a,b) \\

({ \_\_typeof\_\_ (a) \_a = (a); \\

\_\_typeof\_\_ (b) \_b = (b); \\

\_a < \_b ? \_a : \_b; })

\#define max(a,b) \\

({ \_\_typeof\_\_ (a) \_a = (a); \\

\_\_typeof\_\_ (b) \_b = (b); \\

\_a \> \_b ? \_a : \_b; })

typedef enum {TNG\_UNCOMPRESSED,

TNG\_XTC\_COMPRESSION,

TNG\_TNG\_COMPRESSION} tng\_compression;

typedef enum {TNG\_NON\_TRAJECTORY\_BLOCK, TNG\_TRAJECTORY\_BLOCK}
tng\_block\_type;

typedef enum {TNG\_ENDIANNESS\_AND\_STRING\_LENGTH,

TNG\_GENERAL\_INFO,

TNG\_MOLECULES,

TNG\_TRAJECTORY\_IDS\_AND\_NAMES,

TNG\_TRAJECTORY\_FRAME\_SET,

TNG\_BLOCK\_TABLE\_OF\_CONTENTS,

TNG\_PARTICLE\_MAPPING} tng\_non\_trajectory\_block\_ids;

typedef enum {TNG\_TRAJ\_BOX\_SHAPE = 10000,

TNG\_TRAJ\_POSITIONS,

TNG\_TRAJ\_VELOCITIES,

TNG\_TRAJ\_FORCES} tng\_trajectory\_block\_ids;

typedef enum {TNG\_NON\_PARTICLE\_BLOCK\_DATA,

TNG\_PARTICLE\_BLOCK\_DATA} tng\_particle\_block\_data;

/\*typedef enum {TNG\_NO\_HASH, TNG\_OTHER\_HASH,

TNG\_MD5\_HASH, TNG\_MD6\_HASH,

TNG\_SHA0\_HASH, TNG\_SHA1\_HASH,

TNG\_SHA256\_HASH, TNG\_SHA512\_HASH} tng\_hash\_type;\*/

typedef enum {FALSE, TRUE} tng\_bool;

typedef enum {TNG\_CONSTANT\_N\_ATOMS, TNG\_VARIABLE\_N\_ATOMS}

tng\_variable\_n\_atoms\_flag;

typedef enum {TNG\_SUCCESS, TNG\_FAILURE, TNG\_CRITICAL}
tng\_function\_status;

typedef enum {TNG\_NORMAL\_WRITE, TNG\_COPY\_EXISTING} tng\_write\_mode;

typedef enum {TNG\_SKIP\_HASH, TNG\_USE\_HASH} tng\_hash\_mode;

typedef enum {TNG\_CHAR\_DATA,

TNG\_INT\_DATA,

TNG\_FLOAT\_DATA,

TNG\_DOUBLE\_DATA} tng\_data\_type;

typedef struct tng\_trajectory \*tng\_trajectory\_t;

typedef struct tng\_molecule \*tng\_molecule\_t;

typedef struct tng\_chain \*tng\_chain\_t;

typedef struct tng\_residue \*tng\_residue\_t;

typedef struct tng\_atom \*tng\_atom\_t;

typedef struct tng\_bond \*tng\_bond\_t;

\#ifdef \_\_cplusplus

extern "C"

{

\#endif

struct tng\_bond {

int64\_t from\_atom\_id; /\* One of the atoms of the bond \*/

int64\_t to\_atom\_id; /\* The other atom of the bond \*/

};

struct tng\_atom {

tng\_residue\_t residue; /\* The molecule containing this atom \*/

int64\_t id; /\* A unique (per molecule) ID number of the atom \*/

char \*atom\_type; /\* The atom\_type (depending on the forcefield) \*/

char \*name; /\* The name of the atom \*/

};

struct tng\_residue {

tng\_chain\_t chain; /\* The chain containing this residue \*/

int64\_t id; /\* A unique (per chain) ID number of the residue \*/

char \*name; /\* The name of the residue \*/

int64\_t n\_atoms; /\* The number of atoms in the residue \*/

tng\_atom\_t atoms; /\* A list of atoms in the residue \*/

};

struct tng\_chain {

tng\_molecule\_t molecule; /\* The molecule containing this chain \*/

int64\_t id; /\* A unique (per molecule) ID number of the

chain \*/

char \*name; /\* The name of the chain \*/

int64\_t n\_residues; /\* The number of residues in the chain \*/

tng\_residue\_t residues; /\* A list of residues in the chain \*/

};

struct tng\_molecule {

int64\_t id; /\* A unique ID number of the molecule \*/

int64\_t quaternary\_str; /\* Quaternary structure of the molecule.

1 =\> monomeric

2 =\> dimeric

3 =\> trimeric

etc \*/

int64\_t n\_chains; /\* The number of chains in the molecule \*/

int64\_t n\_residues; /\* The number of residues in the molecule \*/

int64\_t n\_atoms; /\* The number of atoms in the molecule \*/

int64\_t n\_bonds; /\* The number of bonds in the molecule.

If the bonds are not specified this

value can be 0. \*/

char \*name; /\* The name of the molecule \*/

tng\_chain\_t chains; /\* A list of chains in the molecule \*/

tng\_residue\_t residues; /\* A list of residues in the molecule \*/

tng\_atom\_t atoms; /\* A list of the atoms in the molecule \*/

tng\_bond\_t bonds; /\* A list of the bonds in the molecule \*/

};

struct tng\_gen\_block {

int64\_t header\_contents\_size; /\* The size of the block header in
bytes \*/

int64\_t block\_contents\_size; /\* The size of the block contents in
bytes \*/

int64\_t id; /\* The ID of the block to determine its

type \*/

char hash[TNG\_HASH\_LEN]; /\* The MD5 hash of the block to verify

integrity \*/

char \*name; /\* The name of the block \*/

int64\_t block\_version; /\* The library version used to write the block
\*/

char \*header\_contents; /\* The full block header contents \*/

char \*block\_contents; /\* The full block contents \*/

};

struct tng\_frame\_set\_toc {

int64\_t n\_blocks; /\* The number of blocks listed in this

table of contents \*/

char \*\*block\_names; /\* A list of block names \*/

};

struct tng\_particle\_mapping {

int64\_t num\_first\_particle; /\* The index number of the first
particle

in this mapping block \*/

int64\_t n\_particles; /\* The number of particles list in this

mapping block \*/

int64\_t \*real\_particle\_numbers; /\* the mapping of index numbers to
the

real particle numbers in the

trajectory \*/

};

struct tng\_trajectory\_frame\_set {

struct tng\_frame\_set\_toc contents; /\* The table of contents of this
frame

set \*/

int64\_t n\_mapping\_blocks; /\* The number of different particle

mapping blocks present. \*/

struct tng\_particle\_mapping \*mappings; /\* The atom mappings of this
frame

set \*/

int64\_t first\_frame; /\* The first frame of this frame set \*/

int64\_t n\_frames; /\* The number of frames in this frame

set \*/

int64\_t \*molecule\_cnt\_list; /\* A list of the number of each
molecule

type - only used when using variable

number of atoms \*/

int64\_t n\_particles; /\* The number of particles/atoms - only

used when using variable number of

atoms \*/

int64\_t next\_frame\_set\_file\_pos; /\* The file position of the next
frame

set \*/

int64\_t prev\_frame\_set\_file\_pos; /\* The file position of the
previous

frame set \*/

int64\_t long\_stride\_next\_frame\_set\_file\_pos; /\* The file
position of the

frame set one stride step

ahead \*/

int64\_t long\_stride\_prev\_frame\_set\_file\_pos; /\* The file
position of the

frame set one stride step

behind \*/

/\* The data blocks in a frame set are trajectory data blocks \*/

int n\_particle\_data\_blocks; /\* The number of data blocks

of particle dependent data \*/

struct tng\_particle\_data \*tr\_particle\_data; /\* A list of data
blocks

containing particle

dependent data \*/

int n\_data\_blocks;

struct tng\_data \*tr\_data;

};

/\* Data can be either double, float, int or a string \*/

union data\_values {

double d;

float f;

int i;

char \*c;

};

/\* FIXME: Should there be a pointer to a tng\_gen\_block from each data
block? \*/

struct tng\_particle\_data {

int64\_t block\_id; /\* The block ID of the data block

containing this particle data.

This is used to determine the

kind of data that is stored \*/

char \*block\_name; /\* The name of the data block.

This is used to determine the

kind of data that is stored \*/

tng\_data\_type datatype; /\* The type of data stored. \*/

int64\_t first\_frame\_with\_data; /\* The first frame number of the

first data point \*/

int64\_t n\_frames; /\* The number of frames in this

frame set \*/

int64\_t n\_values\_per\_frame; /\* The number of values stored per

frame \*/

int64\_t stride\_length; /\* The number of frames between

each data point - e.g. when

storing sparse data. \*/

int64\_t codec\_id; /\* ID of the CODEC used for compression

0 == no compression. \*/

double compression\_multiplier; /\* The multiplier used for getting

integer values for compression \*/

union data\_values \*\*\*values; /\* A 3-dimensional array of values,

sized n\_frames x n\_particles x

n\_values\_per\_frame \*/

};

struct tng\_data {

int64\_t block\_id;

char \*block\_name; /\* The name of the data block.

This is used to determine the

kind of data that is stored \*/

tng\_data\_type datatype; /\* The type of data stored. \*/

int64\_t first\_frame\_with\_data;

int64\_t n\_frames;

int64\_t n\_values\_per\_frame;

int64\_t stride\_length;

int64\_t codec\_id; /\* ID of the CODEC used for compression

0 == no compression. \*/

double compression\_multiplier;

union data\_values \*\*values; /\* A 2-dimensional array of values,

sized n\_frames x n\_values\_per\_frame \*/

};

struct tng\_trajectory {

char \*input\_file\_path; /\* The path of the input trajectory file \*/

FILE \*input\_file; /\* A handle to the input file \*/

long int input\_file\_pos; /\* The reading position of the file \*/

long int input\_file\_len; /\* The length of the input file \*/

char \*output\_file\_path; /\* The path of the output trajectory

file \*/

FILE \*output\_file; /\* A handle to the output file \*/

long int output\_file\_pos; /\* The writing position of the file \*/

tng\_endianness\_32 endianness\_32; /\* The endianness of 32 bit values
of

the current computer \*/

tng\_endianness\_64 endianness\_64; /\* The endianness of 64 bit values
of

the current computer \*/

char \*program\_name; /\* The name of the program producing

this trajectory \*/

char \*forcefield\_name; /\* The forcefield used in the

simulations \*/

char \*user\_name; /\* The name of the user running the

simulations \*/

int64\_t time; /\* The time (n seconds since 1970) when

the file was created \*/

char \*computer\_name; /\* The name of the computer on which the

simulations were performed \*/

char \*pgp\_signature; /\* The PGP signature of the user creating

the file. \*/

char var\_num\_atoms\_flag; /\* A flag indicating if the number of atoms

can vary throughout the simulation, e.g.

using a grand canonical ensemble \*/

int64\_t frame\_set\_n\_frames; /\* The number of frames in a frame set.

It is allowed to have frame sets with

fewer frames, but this will help

searching for specific frames \*/

int64\_t stride\_length; /\* The number of frame sets in a long

stride step \*/

int64\_t n\_molecules; /\* The number of different kinds of

molecules in the trajectory \*/

tng\_molecule\_t molecules; /\* A list of molecules in the trajectory
\*/

int64\_t \*molecule\_cnt\_list; /\* A list of the count of each molecule
-

if using variable number of particles

this will be specified in each frame

set \*/

int64\_t n\_particles; /\* The total number of particles/atoms. If

using variable number of particles this

will be specified in each frame set \*/

int64\_t n\_id\_name\_pairs; /\* The number of ID-name

pairs \*/

struct tng\_block\_id\_name\_pair \*id\_name\_pairs; /\* A list of
ID-name

pairs \*/

int64\_t first\_trajectory\_frame\_set\_input\_file\_pos; /\* The pos in
the src

file of the first

frame set \*/

int64\_t first\_trajectory\_frame\_set\_output\_file\_pos; /\* The pos
in the dest

file of the first

frame set \*/

int64\_t last\_trajectory\_frame\_set\_input\_file\_pos; /\* The pos in
the src

file of the last

frame set \*/

int64\_t last\_trajectory\_frame\_set\_output\_file\_pos; /\* The pos in
the dest

file of the last

frame set \*/

struct tng\_trajectory\_frame\_set current\_trajectory\_frame\_set;/\*
The current-

ly active

frame set \*/

long int current\_trajectory\_frame\_set\_input\_file\_pos; /\* The pos
in the src

file of the current

frame set \*/

long int current\_trajectory\_frame\_set\_output\_file\_pos;/\* The pos
in the dest

file of the current

frame set \*/

int64\_t n\_trajectory\_frame\_sets; /\* The number of frame

sets in the

trajectory \*/

int64\_t n\_trajectory\_blocks; /\* The number of trajectory blocks

in the file \*/

int64\_t n\_non\_trajectory\_blocks; /\* The number of non-trajectory
blocks

in the file \*/

struct tng\_gen\_block non\_trajectory\_blocks[32]; /\* A list of
non-trajectory

blocks \*/

/\* These data blocks are non-trajectory data blocks \*/

int n\_particle\_data\_blocks; /\* The number of non-frame

dependent particle

dependent data blocks \*/

struct tng\_particle\_data \*non\_tr\_particle\_data; /\* A list of data
blocks

containing particle

dependent data \*/

int n\_data\_blocks; /\* The number of frame and particle

independent data blocks \*/

struct tng\_data \*non\_tr\_data; /\* A list of frame and particle

indepdendent data blocks \*/

};

/\* Setup a trajectory data container.

\*tng\_data is a pointer to pre-allocated memory.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_trajectory\_init(tng\_trajectory\_t
tng\_data);

/\* Clean up a trajectory data container.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

All allocated memory in the data structure is freed.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_trajectory\_destroy(tng\_trajectory\_t
tng\_data);

/\* Set the name of the input file.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_input\_file\_set(tng\_trajectory\_t
tng\_data,

const char \*file\_name);

/\* Set the name of the output file.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_output\_file\_set(tng\_trajectory\_t
tng\_data,

const char \*file\_name);

/\* Set the program name used when creating the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_program\_name\_set(tng\_trajectory\_t
tng\_data,

const char \*new\_name);

/\* Set the name of the forcefield used in the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_forcefield\_name\_set(tng\_trajectory\_t
tng\_data,

const char \*new\_name);

/\* Set the name of the user creating the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_user\_name\_set(tng\_trajectory\_t tng\_data,

const char \*new\_name);

/\* Set the name of the computer used when creating the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*new\_name is a pointer to the string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_computer\_name\_set(tng\_trajectory\_t
tng\_data,

const char \*new\_name);

/\* Set the pgp\_signature of the trajectory.

\*tng\_data is a pointer to pre-allocated memory containing trajectory
data.

\*signature is a pointer to the string containing the pgp\_signature.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_signature\_set(tng\_trajectory\_t tng\_data,

const char \*signature);

/\* Setup a molecule container.

\*molecule is a pointer to pre-allocated memory.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_init(tng\_molecule\_t molecule);

/\* Clean up a molecule container.

\*molecule is a pointer to pre-allocated memory containing a molecule.

All allocated memory in the data structure is freed.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_destroy(tng\_molecule\_t molecule);

/\* Setup a data block.

\*block is a pointer to pre-allocated memory.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_block\_init(struct tng\_gen\_block \*block);

/\* Clean up a data block.

\*block is a pointer to pre-allocated memory.

All allocated memory in the data structure is freed.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_block\_destroy(struct tng\_gen\_block
\*block);

/\* Set the name of a data block.

tng\_data is the trajectory data container containing the block..

\*block is a pointer to the block to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_block\_name\_set(tng\_trajectory\_t
tng\_data,

struct tng\_gen\_block \*block,

const char \*new\_name);

/\* Add a molecule to the trajectory.

tng\_data is the trajectory data container containing the block..

\*name is a pointer to the string containing the name of the new
molecule.

\*\*molecule is a pointer to a pointer to the newly created molecule.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_add(tng\_trajectory\_t tng\_data,

const char \*name,

tng\_molecule\_t \*molecule);

/\* Set the name of a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_name\_set(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

const char \*new\_name);

/\* Set the count of a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to rename.

\*cnt is a pointer to the variable to be populated with the count.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_cnt\_get(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

int64\_t \*cnt);

/\* Set the count of a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to rename.

cnt is the number of instances of this molecule.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_cnt\_set(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

int64\_t cnt);

/\* Add a chain to a molecule.

tng\_data is the trajectory data container containing the molecule..

\*molecule is a pointer to the molecule to add a chain to.

\*name is a pointer to a string containing the name of the chain.

\*\*chain is a pointer to a pointer to the newly created chain.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_molecule\_chain\_add(tng\_trajectory\_t
tng\_data,

tng\_molecule\_t molecule,

const char \*name,

tng\_chain\_t \*chain);

/\* Set the name of a chain.

tng\_data is the trajectory data container containing the atom..

\*chain is a pointer to the chain to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_chain\_name\_set(tng\_trajectory\_t
tng\_data,

tng\_chain\_t chain,

const char \*new\_name);

/\* Add a residue to a chain.

tng\_data is the trajectory data container containing the chain..

\*chain is a pointer to the chain to add a residue to.

\*name is a pointer to a string containing the name of the residue.

\*\*residue is a pointer to a pointer to the newly created residue.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_chain\_residue\_add(tng\_trajectory\_t
tng\_data,

tng\_chain\_t chain,

const char \*name,

tng\_residue\_t \*residue);

/\* Set the name of a residue.

tng\_data is the trajectory data container containing the residue.due.

\*residue is a pointer to the residue to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_residue\_name\_set(tng\_trajectory\_t
tng\_data,

tng\_residue\_t residue,

const char \*new\_name);

/\* Add an atom to a residue.

\*tng\_data is a pointer to the trajectory containing the residue.

\*residue is a pointer to the residue to add an atom to.

\*atom\_name is a pointer to a string containing the name of the atom.

\*atom\_type is a pointer to a string containing the atom type of the
atom.

\*\*atom is a pointer to a pointer to the newly created atom.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_residue\_atom\_add(tng\_trajectory\_t
tng\_data,

tng\_residue\_t residue,

const char \*atom\_name,

const char \*atom\_type,

tng\_atom\_t \*atom);

/\* Set the name of an atom.

tng\_data is the trajectory data container containing the atom..

\*atom is a pointer to the atom to rename.

\*new\_name is a string containing the wanted name.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_atom\_name\_set(tng\_trajectory\_t tng\_data,

tng\_atom\_t atom,

const char \*new\_name);

/\* Set the atom type of an atom.

tng\_data is the trajectory data container containing the atom..

\*atom is a pointer to the atom to change.

\*new\_name is a string containing the atom type.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_atom\_type\_set(tng\_trajectory\_t tng\_data,

tng\_atom\_t atom,

const char \*new\_type);

/\* Read the header blocks from the input\_file of tng\_data.

The trajectory blocks must be read separately and iteratively in chunks

to fit in memory.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

If hash\_mode == TNG\_USE\_HASH the written md5 hash in the file will be

compared to the md5 hash of the read contents to ensure valid data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_file\_headers\_read(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

/\* Write the header blocks to the output\_file of tng\_data.

The trajectory blocks must be written separately and iteratively in
chunks

to fit in memory.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path

specifies which file to write to. If the file (output\_file) is not open
it

will be opened.

If hash\_mode == TNG\_USE\_HASH an md5 hash for each header block will
be generated.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_file\_headers\_write(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

/\* Read one (the next) block (of any kind) from the input\_file of
tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

\*block\_data is a pointer to the struct which will be populated with
the

data. If block\_data-\>input\_file\_pos \> 0 it is the position from
where the

reading starts otherwise it starts from the current position.

If hash\_mode == TNG\_USE\_HASH the written md5 hash in the file will be

compared to the md5 hash of the read contents to ensure valid data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor

error has occurred or TNG\_CRITICAL (2) if a major error has occured.
\*/

tng\_function\_status tng\_block\_read\_next(tng\_trajectory\_t
tng\_data,

struct tng\_gen\_block \*block\_data,

const tng\_hash\_mode hash\_mode);

/\* Read one (the next) frame set, including toc, mapping and related
data blocks

from the input\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

If hash\_mode == TNG\_USE\_HASH the written md5 hash in the file will be

compared to the md5 hash of the read contents to ensure valid data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_set\_read\_next(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

/\* Write one frame set, including toc, mapping and related data blocks

to the output\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path
specifies

which file to write to. If the file (output\_file) is not open it will
be

opened.

If hash\_mode == TNG\_USE\_HASH an md5 hash for each block will be
generated.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_set\_write(tng\_trajectory\_t
tng\_data,

const tng\_hash\_mode hash\_mode);

tng\_function\_status tng\_frame\_set\_new(tng\_trajectory\_t tng\_data,

const int64\_t first\_frame,

const int64\_t n\_frames);

tng\_function\_status tng\_particle\_data\_block\_add(tng\_trajectory\_t
tng\_data,

const int64\_t id,

const char \*block\_name,

const char datatype,

const int64\_t n\_frames,

const int64\_t n\_values\_per\_frame,

const int64\_t stride\_length,

const int64\_t first\_particle\_number,

const int64\_t n\_particles,

int64\_t codec\_id,

void \*new\_data);

/\* Read one trajectory block from the input\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the ID of the block to read.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_traj\_block\_read(tng\_trajectory\_t
tng\_data,

int64\_t block\_id);

/\* Write one trajectory block to the output\_file of tng\_data.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path
specifies

which file to write to. If the file (output\_file) is not open it will
be

opened.

block\_id is the ID of the block to write.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_traj\_block\_write(tng\_trajectory\_t
tng\_data,

int64\_t block\_id);

/\* Read a requested frame set.

tng\_data is a trajectory data container.

tng\_data-\>current\_trajectory\_frame\_set will be the read frame set.

frame\_set\_nr is the number of the frame set to return (starting from
0).

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_set\_read\_nr(tng\_trajectory\_t
tng\_data,

int64\_t frame\_set\_nr);

/\* Read a number of consecutive trajectory frames from the input\_file
of

tng\_data.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

start\_frame\_nr is the index number of the first frame to read.

end\_frame\_nr is the index number of the last frame to read.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_read\_interval(tng\_trajectory\_t
tng\_data,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr);

/\* Write a number of consecutive trajectory frames to the output\_file
of

tng\_data.

tng\_data is a trajectory data container. tng\_data-\>output\_file\_path
specifies

which file to write to. If the file (output\_file) is not open it will
be

opened.

start\_frame\_nr is the index number of the first frame to write.

end\_frame\_nr is the index number of the last frame to write.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_frame\_write\_interval(tng\_trajectory\_t
tng\_data,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr);

/\* Free data is an array of values (2D).

\*\*values is the array to free and will be set to 0 afterwards.

n\_frames is the number of frames in the data array.

n\_values\_per\_frame is the number of values per frame in the data
array.

type is the data type of the data in the array.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_data\_values\_free(union data\_values
\*\*values,

int64\_t n\_frames,

int64\_t n\_values\_per\_frame,

tng\_data\_type type);

/\* Free data is an array of values (3D).

\*\*\*values is the array to free and will be set to 0 afterwards.

n\_frames is the number of frames in the data array.

n\_particles is the number of particles in the data array.

n\_values\_per\_frame is the number of values per frame in the data
array.

type is the data type of the data in the array.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_particle\_data\_values\_free(union
data\_values \*\*\*values,

int64\_t n\_frames,

int64\_t n\_particles,

int64\_t n\_values\_per\_frame,

tng\_data\_type type);

/\* Retrieve non-particle data from the last read frame set.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the data block to read.

\*\*\*values is a pointer to a 2-dimensional array (memory unallocated),
which

will point to the data of the requested data block. The array will be
sized

(n\_frames \* n\_values\_per\_frame).

Since \*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_frames is set to the number of frames in the data. This is needed
to

properly reach and/or free the data.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_data\_get(tng\_trajectory\_t tng\_data,

int64\_t block\_id,

union data\_values \*\*\*values,

int64\_t \*n\_frames,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Read and retrieve non-particle data, in a specific interval, from
the

trajectory.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the data block to read.

start\_frame\_nr is the index number of the first frame to read.

end\_frame\_nr is the index number of the last frame to read.

\*\*\*values is a pointer to a 2-dimensional array (memory unallocated),
which

will be filled with data. The array will be sized

(n\_frames \* n\_values\_per\_frame).

Since \*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_data\_interval\_get(tng\_trajectory\_t
tng\_data,

int64\_t block\_id,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr,

union data\_values \*\*\*values,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Retrieve particle data, from the last read frame set.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the particle data block to read.

\*\*\*\*values is a pointer to a 3-dimensional array (memory
unallocated), which

will point to the data of the requested data block. The array will be
sized

(n\_frames \* n\_particles \* n\_values\_per\_frame).

Since \*\*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_frames is set to the number of frames in the data. This is needed
to

properly reach and/or free the data.

\*n\_particles is set to the number of particles in the returned data.
This is

needed to properly reach and/or free the data.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_particle\_data\_get(tng\_trajectory\_t
tng\_data,

int64\_t block\_id,

union data\_values \*\*\*\*values,

int64\_t \*n\_frames,

int64\_t \*n\_particles,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Read and retrieve particle data, in a specific interval, from the
trajectory.

tng\_data is a trajectory data container. tng\_data-\>input\_file\_path
specifies

which file to read from. If the file (input\_file) is not open it will
be

opened.

block\_id is the id number of the particle data block to read.

start\_frame\_nr is the index number of the first frame to read.

end\_frame\_nr is the index number of the last frame to read.

first\_particle\_number is the number of the first particle, for which
to

retrive data.

last\_particle\_number is the number of the last particle, for which to

retrieve data.

\*\*\*\*values is a pointer to a 3-dimensional array (memory
unallocated), which

will be filled with data. The array will be sized

(n\_frames \* n\_particles \* n\_values\_per\_frame).

Since \*\*\*\*values is allocated in this function it is the callers

responsibility to free the memory.

\*n\_values\_per\_frame is set to the number of values per frame in the
data.

This is needed to properly reach and/or free the data.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status
tng\_particle\_data\_interval\_get(tng\_trajectory\_t tng\_data,

int64\_t block\_id,

int64\_t start\_frame\_nr,

int64\_t end\_frame\_nr,

int64\_t first\_particle\_number,

int64\_t last\_particle\_number,

union data\_values \*\*\*\*values,

int64\_t \*n\_values\_per\_frame,

tng\_data\_type \*type);

/\* Get the date and time of initial file creation in ISO format
(string).

tng\_data is a trajectory data container.

\*time is a pointer to the string in which the date will be stored.
Memory

must be reserved beforehand.

Returns TNG\_SUCCESS (0) if successful, TNG\_FAILURE (1) if a minor
error

has occurred or TNG\_CRITICAL (2) if a major error has occured. \*/

tng\_function\_status tng\_time\_get\_str(tng\_trajectory\_t tng\_data,
char \*time);

\#ifdef \_\_cplusplus

} /\* end extern "C" \*/

\#endif

\#endif /\* \_TNGIO\_H \*/

^[[g]](#cmnt7)^

[[a]](#cmnt_ref1)magnus.lundborg:

Currently sizes and offsets are not in the TOC block. I think it needs
further testing to decide if it is good or not.

* * * * *

Sander Pronk:

So how can you find out where the block is?

* * * * *

magnus.lundborg:

If the offsets are not listed in the TOC block you would have to read
the whole frame set, or at least all the block headers in the frame set,
which shouldn't be too bad.

[[b]](#cmnt_ref2)Roland Schulz:

I suggest not to use a version number. This is already a problem with
the tpx version and branches. Instead I suggest to have a bitvector
where each bit says whether a certain feature is present in this file.
Given a central registry of meaning of the bits, this allows different
groups/branches/software to add features. Which would be difficult with
a version number approach which has an inherently linear ordering. The
last bit probably should be reserved to signify whether the 64bit
bitvector is extended by a another 64bit. The data should be stored in
the order of the bitvectors. A reader which doesn't support a certain
bit, cannot read any of the following data if the bit is on.

[[c]](#cmnt_ref3)magnus.lundborg:

We will have a problem if we want to add data to a frame set in a file.
All subsequent frame sets will need to be rewritten. One alternative
would be to have a list of pointers to each block of each block type in
the frame set table of contents block. But we will have a problem adding
rows to that block as well, which in turn could be fixed by having a
pointer from the frame set block to the "current" table of contents
block and just let the old one remain. We could actually have a flag in
block headers to show if the block is "up-to-date". But there is a risk
that these pointers will be slow - especially when it comes to writing.

[[d]](#cmnt_ref4)Rossen Apostolov:

somehow the name of the FF doesn't fit naturally with the rest of the
info here :)

How about including the simulation setup in the file in a separate
block? That will be needed if the file can be used for restarts too.

[[e]](#cmnt_ref5)magnus.lundborg:

This introduces a bit of redundancy, but helps keeping track of the
data.

[[f]](#cmnt_ref6)Roland Schulz:

might be good to make this optional. And if it isn't given then the
numbering is consecutive. The would still give the flexibility that one
can specify the first and no of particles which isn't possible without
index block.

* * * * *

Daniel Spångberg:

if this is made optional, the comment below the section can be removed.
and particle mapping blocks required, since it will not cost much extra
to have it.

[[g]](#cmnt_ref7)Rossen Apostolov:

We should think of a different name for the traj. group blcok, it's
confusing
