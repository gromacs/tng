#include <stdlib.h>
#include <stdio.h>
#include <tng_io.h>
#include <omp.h>

int main(int argc, char **argv)
{
    tng_trajectory_t traj, local_traj = 0;
    union data_values ***positions = 0; // A 3-dimensional array to be populated
    int64_t n_particles, n_values_per_frame, n_frame_sets, n_frames;
    tng_data_type data_type;
    int i;
    int64_t particle = 0;
    char atom_name[64], res_name[64];

    if(argc <= 1)
    {
        printf("No file specified\n");
        printf("Usage:\n");
        printf("tng_io_read_pos <tng_file> [particle number = %"PRId64"]\n",
               particle);
        exit(1);
    }

    // A reference must be passed to allocate memory
    if(tng_trajectory_init(&traj) != TNG_SUCCESS)
    {
        tng_trajectory_destroy(&traj);
        exit(1);
    }
    tng_input_file_set(traj, argv[1]);

    // Read the file headers
    tng_file_headers_read(traj, TNG_USE_HASH);

    if(argc >= 3)
    {
        particle = strtoll(argv[2], 0, 10);
    }

    tng_num_frame_sets_get(traj, &n_frame_sets);

    printf("%"PRId64" frame sets\n", n_frame_sets);

    if(tng_atom_name_of_particle_nr_get(traj, particle, atom_name,
                                        sizeof(atom_name)) ==
       TNG_SUCCESS &&
       tng_residue_name_of_particle_nr_get(traj, particle, res_name,
                                           sizeof(res_name)) ==
       TNG_SUCCESS)
    {
        printf("Particle: %s (%s)\n", atom_name, res_name);
    }
    else
    {
        printf("Particle name not found\n");
    }

#pragma omp parallel \
private (n_frames, n_particles, n_values_per_frame, \
         data_type, positions) \
firstprivate (local_traj)
{
    positions = 0;
    tng_trajectory_copy(traj, &local_traj);
#pragma omp for
    for(i = 0; i < n_frame_sets; i++)
    {
        if(tng_frame_set_nr_find(local_traj, i) != TNG_SUCCESS)
        {
            printf("FAILED finding frame set %d!\n", i);
        }
        if(tng_particle_data_get(local_traj, TNG_TRAJ_POSITIONS, &positions,
                                 &n_frames, &n_particles, &n_values_per_frame,
                                 &data_type) != TNG_SUCCESS)
        {
            printf("FAILED getting particle data\n");
        }
//         printf("%"PRId64" %"PRId64" %"PRId64"\n", n_frames, n_particles, n_values_per_frame);
    }
    // Free memory
    if(positions)
    {
        tng_particle_data_values_free(local_traj, positions, n_frames, n_particles,
                                      n_values_per_frame, data_type);
    }
    tng_trajectory_destroy(&local_traj);
}
    tng_trajectory_destroy(&traj);

    return(0);
}
