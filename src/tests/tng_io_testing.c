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

    tng_add_molecule(traj, "water", &molecule);
    tng_add_chain_to_molecule(traj, molecule, "W", &chain);
    tng_add_residue_to_chain(traj, chain, "WAT", &residue);
    if(tng_add_atom_to_residue(traj, residue, "O", "O", &atom) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    if(tng_add_atom_to_residue(traj, residue, "HO1", "H", &atom) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    if(tng_add_atom_to_residue(traj, residue, "HO2", "H", &atom) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    tng_set_molecule_cnt(traj, molecule, 200);
    tng_get_molecule_cnt(traj, molecule, &cnt);
    printf("Created %"PRId64" %s molecules.\n", cnt, molecule->name);
//     traj->molecule_cnt_list[traj->n_molecules-1] = 5;
//     tng_set_molecule_name(traj, &traj->molecules[1], "ligand");
//     tng_set_molecule_name(traj, &traj->molecules[2], "water");
//     tng_set_molecule_name(traj, &traj->molecules[3], "dummy");
//     traj->molecules[0].id = 0;
//     traj->molecules[1].id = 1;
//     traj->molecules[2].id = 2;
//     traj->molecules[3].id = 3;

//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom1", "type1") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom2", "type1") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom3", "type1") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom4", "type2") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom5", "type2") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom6", "type2") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[0], "atom7", "type3") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "C1", "C") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "O1", "O") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H11", "H") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H12", "H") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(tng_add_atom_to_molecule(traj, &traj->molecules[1], "H13", "H") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }

    return(TRG_SUCCESS);
}

static tng_function_status tng_test_read_and_write_file
                (tng_trajectory_t traj)
{
    tng_function_status stat;

    stat = tng_read_file_headers(traj, TRG_KEEP_FILE_OPEN);
    if(stat == TRG_CRITICAL)
    {
        return(stat);
    }
    stat = tng_write_file_headers(traj, TRG_KEEP_FILE_OPEN);
    if(stat == TRG_CRITICAL)
    {
        return(stat);
    }

    while(stat != TRG_CRITICAL && traj->input_file_pos < traj->input_file_len &&
          traj->current_trajectory_frame_set.next_frame_set_file_pos != -1UL)
    {
        stat = tng_read_next_frame_set(traj, TRG_KEEP_FILE_OPEN);
        if(stat == TRG_CRITICAL)
        {
            return(stat);
        }
        stat = tng_write_frame_set(traj, TRG_KEEP_FILE_OPEN);
    }
    
    return(stat);
}

static tng_function_status tng_test_write_and_read_traj(tng_trajectory_t traj)
{
    int i, j, k, nr, tot_n_mols, cnt;
    float *data, *molpos;
    tng_function_status stat;
    
    /* Create molecules */
    if(tng_setup_test_molecules(traj) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }

    if(tng_init_block(&traj->non_trajectory_blocks[traj->n_non_trajectory_blocks]) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    traj->non_trajectory_blocks[traj->n_non_trajectory_blocks].id = TRG_MOLECULES;
    if(tng_set_block_name(traj,
                          &traj->non_trajectory_blocks[traj->n_non_trajectory_blocks++],
                          "MOLECULES") == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }

    if(tng_write_file_headers(traj, TRG_KEEP_FILE_OPEN) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    
    data = malloc(sizeof(float) * traj->n_particles *
                  traj->frame_set_n_frames * 3);
    if(!data)
    {
        printf("Cannot allocate memory. %s: %d\n", __FILE__, __LINE__);
        return(TRG_CRITICAL);
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
        if(tng_new_frame_set(traj, i * traj->frame_set_n_frames,
            traj->frame_set_n_frames) != TRG_SUCCESS)
        {
            printf("Error creating frame set %d. %s: %d\n",
                   i, __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TRG_CRITICAL);
        }

        if(tng_add_particle_data_block(traj, TRG_TRAJ_POSITIONS,
                                       "POSITIONS",
                                       TRG_FLOAT_DATA,
                                       traj->frame_set_n_frames, 3,
                                       1, 0, traj->n_particles,
                                       TRG_UNCOMPRESSED,
                                       data) != TRG_SUCCESS)
        {
            printf("Error adding data. %s: %d\n", __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TRG_CRITICAL);
        }
        if(tng_write_frame_set(traj, TRG_KEEP_FILE_OPEN) != TRG_SUCCESS)
        {
            printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TRG_CRITICAL);
        }
    }

//     tng_add_ids_names_pair(traj, TRG_TRAJ_VELOCITIES, "BOX SHAPE");
//     tng_add_ids_names_pair(traj, TRG_TRAJ_POSITIONS, "TRAJECTORY POSITIONS");
//     tng_add_ids_names_pair(traj, TRG_TRAJ_VELOCITIES, "TRAJECTORY VELOCITIES");
//     tng_add_ids_names_pair(traj, TRG_TRAJ_FORCES, "TRAJECTORY FORCES");
//     tng_add_ids_names_pair(traj, 11000, "TEST DATA");
        
    free(molpos);
    free(data);

    tng_destroy_trajectory(traj);
    tng_set_input_file(traj, "/tmp/tng_test.tng");
    
    stat = tng_read_file_headers(traj, TRG_KEEP_FILE_OPEN);
    
    while(stat != TRG_CRITICAL && traj->input_file_pos < traj->input_file_len &&
          traj->current_trajectory_frame_set.next_frame_set_file_pos != -1ULL)
    {
        stat = tng_read_next_frame_set(traj, TRG_KEEP_FILE_OPEN);
        if(stat == TRG_CRITICAL)
        {
            return(stat);
        }
    }

    return(stat);
}

int main()
{
    struct tng_trajectory traj;
    char time_str[TRG_MAX_DATE_STR_LEN];
    
    if(tng_init_trajectory(&traj) != TRG_SUCCESS)
    {
        tng_destroy_trajectory(&traj);
        printf("Test Init trajectory:\t\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    printf("Test Init trajectory:\t\t\t\tSucceeded.\n");

    tng_get_time_str(&traj, time_str);

    printf("Creation time: %s\n", time_str);
    
    tng_set_input_file(&traj, "tng_example.tng");
    tng_set_output_file(&traj, "/tmp/tng_example_test.tng");

//     if(tng_test_endianness(&traj) != TRG_SUCCESS)
//     {
//         printf("Test failed: Endianness. %s: %d\n", __FILE__, __LINE__);
//     }

    if(tng_test_read_and_write_file(&traj) == TRG_CRITICAL)
    {
        printf("Test Read and write file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Read and write file:\t\t\tSucceeded.\n");
    }
    
    if(tng_destroy_trajectory(&traj) == TRG_CRITICAL ||
       tng_init_trajectory(&traj) == TRG_CRITICAL)
    {
        printf("Test Destroy and init trajectory:\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Destroy and init trajectory:\t\tSucceeded.\n");
    }
    

    tng_set_output_file(&traj, "/tmp/tng_test.tng");
    
    if(tng_test_write_and_read_traj(&traj) == TRG_CRITICAL)
    {
        printf("Test Write and read file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Write and read file:\t\t\tSucceeded.\n");
    }
    
    if(tng_destroy_trajectory(&traj) == TRG_CRITICAL)
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