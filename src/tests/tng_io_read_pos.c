#include <stdlib.h>
#include <stdio.h>
#include <tng_io.h>

int main(int argc, char **argv)
{
    tng_trajectory_t traj;
    union data_values ***positions = 0; // A 3-dimensional array to be populated
    int64_t n_particles, n_values_per_frame, n_frames;
    tng_data_type data_type;
    int i, j;
    int64_t particle = 0;
    // Set a default frame range
    int first_frame = 0, last_frame = 50;

    if(argc <= 1)
    {
        printf("No file specified\n");
        printf("Usage:\n");
        printf("tng_io_read_pos <tng_file> [particle number = %"PRId64"] "
               "[first_frame = %d] [last_frame = %d]\n",
               particle, first_frame, last_frame);
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
        if(argc >= 4)
        {
            first_frame = strtoll(argv[3], 0, 10);
            if(argc >= 5)
            {
                last_frame = strtoll(argv[4], 0, 10);
            }
        }
    }

    n_frames = last_frame - first_frame + 1;

    // Get the positions of all particles in the requested frame range.
    // The positions are stored in the positions array.
    // N.B. No proper error checks.
    if(tng_particle_data_interval_get(traj, TNG_TRAJ_POSITIONS, first_frame,
       last_frame, TNG_USE_HASH, &positions, &n_particles, &n_values_per_frame,
       &data_type) == TNG_SUCCESS)
    {
        // Print the positions of the wanted particle (zero based)
        for(i=0; i<n_frames; i++)
        {
            printf("%d", first_frame + i);
            for(j=0; j<n_values_per_frame; j++)
            {
                switch(data_type)
                {
                case TNG_INT_DATA:
                    printf("\t%"PRId64"", positions[i][particle][j].i);
                    break;
                case TNG_FLOAT_DATA:
                    printf("\t%f", positions[i][particle][j].f);
                    break;
                case TNG_DOUBLE_DATA:
                    printf("\t%f", positions[i][particle][j].d);
                    break;
                default:
                    break;
                }
                printf("\n");
            }
        }
    }
    else
    {
        printf("Cannot read positions\n");
    }

    // Free memory
    if(positions)
    {
        tng_particle_data_values_free(traj, positions, n_frames, n_particles,
                                      n_values_per_frame, data_type);
    }
    tng_trajectory_destroy(&traj);

    return(0);
}
