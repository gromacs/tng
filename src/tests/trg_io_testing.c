#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "trg_io.h"




static trg_function_status trg_setup_test_molecules(struct trg_trajectory *traj)
{
    struct trg_molecule *molecule;
    struct trg_chain *chain;
    struct trg_residue *residue;
    struct trg_atom *atom;
    int64_t cnt;
//     int i;

    trg_add_molecule(traj, "water", &molecule);
    trg_add_chain_to_molecule(traj, molecule, "W", &chain);
    trg_add_residue_to_chain(traj, chain, "WAT", &residue);
    if(trg_add_atom_to_residue(traj, residue, "O", "O", &atom) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    if(trg_add_atom_to_residue(traj, residue, "HO1", "H", &atom) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    if(trg_add_atom_to_residue(traj, residue, "HO2", "H", &atom) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    trg_set_molecule_cnt(traj, molecule, 200);
    trg_get_molecule_cnt(traj, molecule, &cnt);
    printf("Created %"PRId64" %s molecules.\n", cnt, molecule->name);
//     traj->molecule_cnt_list[traj->n_molecules-1] = 5;
//     trg_set_molecule_name(traj, &traj->molecules[1], "ligand");
//     trg_set_molecule_name(traj, &traj->molecules[2], "water");
//     trg_set_molecule_name(traj, &traj->molecules[3], "dummy");
//     traj->molecules[0].id = 0;
//     traj->molecules[1].id = 1;
//     traj->molecules[2].id = 2;
//     traj->molecules[3].id = 3;

//     if(trg_add_atom_to_molecule(traj, &traj->molecules[0], "atom1", "type1") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[0], "atom2", "type1") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[0], "atom3", "type1") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[0], "atom4", "type2") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[0], "atom5", "type2") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[0], "atom6", "type2") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[0], "atom7", "type3") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[1], "C1", "C") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[1], "O1", "O") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[1], "H11", "H") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[1], "H12", "H") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }
//     if(trg_add_atom_to_molecule(traj, &traj->molecules[1], "H13", "H") == TRG_CRITICAL)
//     {
//         return(TRG_CRITICAL);
//     }

    return(TRG_SUCCESS);
}

static trg_function_status trg_test_read_and_write_file
                (struct trg_trajectory *traj)
{
    trg_function_status stat;

    stat = trg_read_file_headers(traj, TRG_KEEP_FILE_OPEN);
    if(stat == TRG_CRITICAL)
    {
        return(stat);
    }
    stat = trg_write_file_headers(traj, TRG_KEEP_FILE_OPEN);
    if(stat == TRG_CRITICAL)
    {
        return(stat);
    }

    while(stat != TRG_CRITICAL && traj->input_file_pos < traj->input_file_len &&
          traj->current_trajectory_frame_set.next_frame_set_file_pos != -1UL)
    {
        stat = trg_read_next_frame_set(traj, TRG_KEEP_FILE_OPEN);
        if(stat == TRG_CRITICAL)
        {
            return(stat);
        }
        stat = trg_write_frame_set(traj, TRG_KEEP_FILE_OPEN);
    }
    
    return(stat);
}

static trg_function_status trg_test_write_and_read_traj(struct trg_trajectory *traj)
{
    int i, j, k, nr, tot_n_mols, cnt;
    float *data, *molpos;
    trg_function_status stat;
    
    /* Create molecules */
    if(trg_setup_test_molecules(traj) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }

    if(trg_init_block(&traj->non_trajectory_blocks[traj->n_non_trajectory_blocks]) == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }
    traj->non_trajectory_blocks[traj->n_non_trajectory_blocks].id = TRG_MOLECULES;
    if(trg_set_block_name(traj,
                          &traj->non_trajectory_blocks[traj->n_non_trajectory_blocks++],
                          "MOLECULES") == TRG_CRITICAL)
    {
        return(TRG_CRITICAL);
    }

    if(trg_write_file_headers(traj, TRG_KEEP_FILE_OPEN) == TRG_CRITICAL)
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
        if(trg_new_frame_set(traj, i * traj->frame_set_n_frames,
            traj->frame_set_n_frames) != TRG_SUCCESS)
        {
            printf("Error creating frame set %d. %s: %d\n",
                   i, __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TRG_CRITICAL);
        }

        if(trg_add_particle_data_block(traj, TRG_TRAJ_POSITIONS,
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
        if(trg_write_frame_set(traj, TRG_KEEP_FILE_OPEN) != TRG_SUCCESS)
        {
            printf("Error writing frame set. %s: %d\n", __FILE__, __LINE__);
            free(molpos);
            free(data);
            return(TRG_CRITICAL);
        }
    }

//     trg_add_ids_names_pair(traj, TRG_TRAJ_VELOCITIES, "BOX SHAPE");
//     trg_add_ids_names_pair(traj, TRG_TRAJ_POSITIONS, "TRAJECTORY POSITIONS");
//     trg_add_ids_names_pair(traj, TRG_TRAJ_VELOCITIES, "TRAJECTORY VELOCITIES");
//     trg_add_ids_names_pair(traj, TRG_TRAJ_FORCES, "TRAJECTORY FORCES");
//     trg_add_ids_names_pair(traj, 11000, "TEST DATA");
        
    free(molpos);
    free(data);

    trg_destroy_trajectory(traj);
    trg_set_input_file(traj, "/tmp/trg_test.trg");
    
    stat = trg_read_file_headers(traj, TRG_KEEP_FILE_OPEN);
    
    while(stat != TRG_CRITICAL && traj->input_file_pos < traj->input_file_len &&
          traj->current_trajectory_frame_set.next_frame_set_file_pos != -1ULL)
    {
        stat = trg_read_next_frame_set(traj, TRG_KEEP_FILE_OPEN);
        if(stat == TRG_CRITICAL)
        {
            return(stat);
        }
    }

    return(stat);
}

int main()
{
    struct trg_trajectory traj;
    char time_str[TRG_MAX_DATE_STR_LEN];
    
    if(trg_init_trajectory(&traj) != TRG_SUCCESS)
    {
        trg_destroy_trajectory(&traj);
        printf("Test Init trajectory:\t\t\t\tFailed. %s: %d.\n",
               __FILE__, __LINE__);
        exit(1);
    }
    printf("Test Init trajectory:\t\t\t\tSucceeded.\n");

    trg_get_time_str(&traj, time_str);

    printf("Creation time: %s\n", time_str);
    
    trg_set_input_file(&traj, "trg_example.trg");
    trg_set_output_file(&traj, "/tmp/trg_example_test.trg");

//     if(trg_test_endianness(&traj) != TRG_SUCCESS)
//     {
//         printf("Test failed: Endianness. %s: %d\n", __FILE__, __LINE__);
//     }

    if(trg_test_read_and_write_file(&traj) == TRG_CRITICAL)
    {
        printf("Test Read and write file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Read and write file:\t\t\tSucceeded.\n");
    }
    
    if(trg_destroy_trajectory(&traj) == TRG_CRITICAL ||
       trg_init_trajectory(&traj) == TRG_CRITICAL)
    {
        printf("Test Destroy and init trajectory:\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Destroy and init trajectory:\t\tSucceeded.\n");
    }
    

    trg_set_output_file(&traj, "/tmp/trg_test.trg");
    
    if(trg_test_write_and_read_traj(&traj) == TRG_CRITICAL)
    {
        printf("Test Write and read file:\t\t\tFailed. %s: %d\n",
               __FILE__, __LINE__);
    }
    else
    {
        printf("Test Write and read file:\t\t\tSucceeded.\n");
    }
    
    if(trg_destroy_trajectory(&traj) == TRG_CRITICAL)
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