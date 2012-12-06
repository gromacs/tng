#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "tng_io.h"




static tng_function_status tng_setup_test_molecules(tng_trajectory_t traj)
{
    struct tng_molecule *molecule;
    struct tng_chain *chain;
    struct tng_residue *residue;
    struct tng_atom *atom;
    int64_t cnt;
//     int i;

    tng_molecule_add(traj, "water", &molecule);
    tng_molecule_chain_add(traj, molecule, "W", &chain);
    tng_chain_residue_add(traj, chain, "WAT", &residue);
    if(tng_residue_atom_add(traj, residue, "O", "O", &atom) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    if(tng_residue_atom_add(traj, residue, "HO1", "H", &atom) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    if(tng_residue_atom_add(traj, residue, "HO2", "H", &atom) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    tng_molecule_cnt_set(traj, molecule, 200);
    tng_molecule_cnt_get(traj, molecule, &cnt);
    printf("Created %"PRId64" %s molecules.\n", cnt, molecule->name);
//     traj->molecule_cnt_list[traj->n_molecules-1] = 5;
//     tng_molecule_name_set(traj, &traj->molecules[1], "ligand");
//     tng_molecule_name_set(traj, &traj->molecules[2], "water");
//     tng_molecule_name_set(traj, &traj->molecules[3], "dummy");
//     traj->molecules[0].id = 0;
//     traj->molecules[1].id = 1;
//     traj->molecules[2].id = 2;
//     traj->molecules[3].id = 3;

//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom1", "type1") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom2", "type1") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom3", "type1") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom4", "type2") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom5", "type2") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom6", "type2") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom7", "type3") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "C1", "C") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "O1", "O") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H11", "H") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H12", "H") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H13", "H") == TNG_CRITICAL)
//     {
//         return(TNG_CRITICAL);
//     }

    return(TNG_SUCCESS);
}

static tng_function_status tng_test_read_and_write_file
                (tng_trajectory_t traj)
{
    tng_function_status stat;

    stat = tng_file_headers_read(traj, TNG_USE_HASH);
    if(stat == TNG_CRITICAL)
    {
        return(stat);
    }
    stat = tng_file_headers_write(traj, TNG_USE_HASH);
    if(stat == TNG_CRITICAL)
    {
        return(stat);
    }

    while(stat != TNG_CRITICAL && traj->input_file_pos < traj->input_file_len &&
          traj->current_trajectory_frame_set.next_frame_set_file_pos != -1UL)
    {
        stat = tng_frame_set_read_next(traj, TNG_USE_HASH);
        if(stat == TNG_CRITICAL)
        {
            return(stat);
        }
        stat = tng_frame_set_write(traj, TNG_USE_HASH);
    }
    
    return(stat);
}

static tng_function_status tng_test_write_and_read_traj(tng_trajectory_t traj)
{
    int i, j, k, nr, tot_n_mols, cnt;
    float *data, *molpos;
    tng_function_status stat;
    
    /* Create molecules */
    if(tng_setup_test_molecules(traj) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }

    if(tng_block_init(&traj->non_trajectory_blocks[traj->n_non_trajectory_blocks]) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    traj->non_trajectory_blocks[traj->n_non_trajectory_blocks].id = TNG_MOLECULES;
    if(tng_block_name_set(traj,
                          &traj->non_trajectory_blocks[traj->n_non_trajectory_blocks++],
                          "MOLECULES") == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }

    if(tng_file_headers_write(traj, TNG_SKIP_HASH) == TNG_CRITICAL)
    {
        return(TNG_CRITICAL);
    }
    
    data = malloc(sizeof(float) * traj->n_particles *
                  traj->frame_set_n_frames * 3);
    if(!data)
    {
        printf("Cannot allocate memory. %s: %d\n", __FILE__, __LINE__);
        return(TNG_CRITICAL);
    }

    tot_n_mols = 0;
    for(i = 0; i < traj->n_molecules; i++)
    {
        tot_n_mols += traj->molecule_cnt_list[i];
    }
    molpos = malloc(sizeof(float) * tot_n_mols * 3);

    /* Set initial coordinates */
    for(i = 0; i < tot_n_mols; i++)
    {
        nr = i * 3;
        /* Somewhat random coordinates (between 0 and 100),
        * but not specifying a random seed */
        molpos[nr] = 100.0 * rand() / (RAND_MAX + 1.0);
        molpos[nr+1] = 100.0 * rand() / (RAND_MAX + 1.0);
        molpos[nr+2] = 100.0 * rand() / (RAND_MAX + 1.0);
    }

    /* Generate 200 frame sets - each with 100 frames (by default) */
    for(i = 0; i < 200; i++)
    {
        cnt = 0;
        for(j = 0; j < traj->frame_set_n_frames; j++)
        {
            for(k = 0; k < tot_n_mols; k++)
            {
                nr = k * 3;
                /* Move -1 to 1 */
                molpos[nr] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;
                molpos[nr+1] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;
                molpos[nr+2] += 2 * (rand() / (RAND_MAX + 1.0)) - 1;

                data[cnt++] = molpos[nr];
                data[cnt++] = molpos[nr + 1];
                data[cnt++] = molpos[nr + 2];
                data[cnt++] = molpos[nr] + 1;
                data[cnt++] = molpos[nr + 1] + 1;
                data[cnt++] = molpos[nr + 2] + 1;
                data[cnt++] = molpos[nr] - 1;
                data[cnt++] = molpos[nr + 1] - 1;
                data[cnt++] = molpos[nr + 2] - 1;
            }
        }
        if(tng_frame_set_new(traj, i * traj->frame_set_n_frames,
            traj->frame_set_n_frames) != TNG_SUCCESS)
        {
            printf("Error creating frame set %d. %s: %d\n",
                   i, __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }

        if(tng_particle_data_block_add(traj, TNG_TRAJ_POSITIONS,
                                       "POSITIONS",
                                       TNG_FLOAT_DATA,
                                       traj->frame_set_n_frames, 3,
                                       1, 0, traj->n_particles,
                                       TNG_UNCOMPRESSED,
                                       data) != TNG_SUCCESS)
        {
            printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }
        if(tng_frame_set_write(traj, TNG_SKIP_HASH) != TNG_SUCCESS)
        {
            printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TNG_CRITICAL);
        }
    }

//     tng_add_ids_names_pair(traj, TNG_TRAJ_VELOCITIES, "BOX SHAPE");
//     tng_add_ids_names_pair(traj, TNG_TRAJ_POSITIONS, "TRAJECTORY POSITIONS");
//     tng_add_ids_names_pair(traj, TNG_TRAJ_VELOCITIES, "TRAJECTORY VELOCITIES");
//     tng_add_ids_names_pair(traj, TNG_TRAJ_FORCES, "TRAJECTORY FORCES");
//     tng_add_ids_names_pair(traj, 11000, "TEST DATA");
        
    free(molpos);
    free(data);

    tng_trajectory_destroy(traj);
    tng_input_file_set(traj, "/tmp/tng_test.tng");
    
    stat = tng_file_headers_read(traj, TNG_SKIP_HASH);
    
    while(stat != TNG_CRITICAL && traj->input_file_pos < traj->input_file_len &&
          traj->current_trajectory_frame_set.next_frame_set_file_pos != -1ULL)
    {
        stat = tng_frame_set_read_next(traj, TNG_SKIP_HASH);
        if(stat == TNG_CRITICAL)
        {
            return(stat);
        }
    }

    return(stat);
}

int main()
{
    struct tng_trajectory traj;
    char time_str[TNG_MAX_DATE_STR_LEN];
    
    if(tng_trajectory_init(&traj) != TNG_SUCCESS)
    {
        tng_trajectory_destroy(&traj);
        printf("Test Init trajectory:\t\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    printf("Test Init trajectory:\t\t\t\tSucceeded.\n");

    tng_time_get_str(&traj, time_str);

    printf("Creation time: %s\n", time_str);
    
    tng_input_file_set(&traj, "tng_example.tng");
    tng_output_file_set(&traj, "/tmp/tng_example_out.tng");

//     if(tng_test_endianness(&traj) != TNG_SUCCESS)
//     {
//         printf("Test failed: Endianness. %s: %d\n", __FILE__, __LINE__);
//     }

    if(tng_test_read_and_write_file(&traj) == TNG_CRITICAL)
    {
        printf("Test Read and write file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Read and write file:\t\t\tSucceeded.\n");
    }
    
    if(tng_trajectory_destroy(&traj) == TNG_CRITICAL ||
       tng_trajectory_init(&traj) == TNG_CRITICAL)
    {
        printf("Test Destroy and init trajectory:\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Destroy and init trajectory:\t\tSucceeded.\n");
    }
    

    tng_output_file_set(&traj, "/tmp/tng_test.tng");
    
    if(tng_test_write_and_read_traj(&traj) == TNG_CRITICAL)
    {
        printf("Test Write and read file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Write and read file:\t\t\tSucceeded.\n");
    }
    
    if(tng_trajectory_destroy(&traj) == TNG_CRITICAL)
    {
        printf("Test Destroy trajectory:\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    else
    {
        printf("Test Destroy trajectory:\t\t\tSucceeded.\n");
    }
    

    printf("Tests finished\n");
    
    exit(0);
}