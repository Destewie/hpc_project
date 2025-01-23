#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#include <time.h>
#include <omp.h>
#include <mpi.h>

#define N_FISHES 30 // Numero di pesci
#define DIMENSIONS 2 // Dimensione dello spazio
#define BOUNDS_MIN 5.0   // Minimum bound of the search space
#define BOUNDS_MAX 30.0    // Maximum bound of the search space
#define BOUNDS_MIN_W 0.1   // Minimum bound of the search space
#define BOUNDS_MAX_W 10.0    // Maximum bound of the search space
#define MAX_ITER 100
#define MAX_INDIVIDUAL_STEP 1.5 // Passo massimo del movimento individuale
#define MAX_VOLITIVE_STEP 0.2 // Passo massimo del movimento volitivo
#define W_SCALE_MIN 1.0
#define W_SCALE_MAX 10.0
#define BREEDING_THRESHOLD 7.0 // minimus threshold of weight to breedh new fishes
#define FUNCTION "min_sphere"   //TODO: Capire se, al posto di fare un controllo su una stringa, possiamo passare alle funzioni direttamente un puntatore ad una funzione (in modo comodo, se no lasciamo perdere)
#define MULTIPLIER -1   // 1 in case of maximization, -1 in case of minimization
#define A 10.0 //rastrigin param

typedef struct {
    double position[DIMENSIONS];
    double new_position[DIMENSIONS];

    double previous_cycle_weight;
    double weight;

    double fitness;
    double new_fitness;

    double max_individual_step;
    double max_volitive_step;
} Fish;

//-------------------------------------------------------------------------------------------
//---------------------------- MATH FUNCTIONS -----------------------------------------------
//-------------------------------------------------------------------------------------------

double rosenbrok(double *x) {
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS - 1; i++) {
        double term1 = 100.0 * pow(x[i + 1] - x[i] * x[i], 2);
        double term2 = pow(1.0 - x[i], 2);
        sum += term1 + term2;
    }
    return sum;
}

double rastrigin(double *x) {
    double sum = A * DIMENSIONS; // Start with A * DIM
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += x[i] * x[i] - A * cos(2 * M_PI * x[i]);
    }
    return sum;
}

double min_sphere(double *x) {
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

double max_sphere(double *x) {
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += -(x[i] * x[i]);
    }
    return sum;
}

double ackley(double *x){
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum1 += x[i] * x[i];
        sum2 += cos(2 * M_PI * x[i]);
    }
    return -20 * exp(-0.2 * sqrt(sum1 / DIMENSIONS)) - exp(sum2 / DIMENSIONS) + 20 + M_E;
}

double objective_function(double *x) {
    if (strcmp(FUNCTION, "min_rosenbrok") == 0) {
        return rosenbrok(x);
    } else if (strcmp(FUNCTION, "min_sphere") == 0) {
        return min_sphere(x);
    } else if (strcmp(FUNCTION, "max_sphere") == 0) {
        return max_sphere(x);
    } else if (strcmp(FUNCTION, "min_rastrigin") == 0) {
        return rastrigin(x);
    } else if (strcmp(FUNCTION, "min_ackley") == 0) {
        return ackley(x);
    }else{
        return 0.0;
    }
}

//-------------------------------------------------------------------------------------------
//---------------------------- TESTING ---------------------------------------------------------
//-------------------------------------------------------------------------------------------

void print_fish(Fish *fish) {
    printf("Fish: ");
    for(int i=0; i<DIMENSIONS; i++){
        printf("pos: %f ", fish->position[i]);
    }
    printf("\tweight: %f \tfitness: %f\n", fish->weight, fish->fitness);
}


//-------------------------------------------------------------------------------------------
//---------------------------- FISH ---------------------------------------------------------
//-------------------------------------------------------------------------------------------

void initFish(Fish *fish) {
    for (int i = 0; i < DIMENSIONS; i++) {
        fish->position[i] = ((double)rand() / RAND_MAX) * (BOUNDS_MAX - BOUNDS_MIN) + BOUNDS_MIN;
        fish->new_position[i] = fish->position[i];
    }

    fish->weight = W_SCALE_MAX / 2;   // Peso iniziale
    fish->previous_cycle_weight = fish->weight;

    fish->fitness = objective_function(fish->position)*MULTIPLIER;        // Fitness iniziale //TODO: capire qual è il valore migliore di inizializzazione
    fish->new_fitness = fish->fitness;

    fish->max_individual_step = MAX_INDIVIDUAL_STEP; //TODO: capire qual è il valore migliore di inizializzazione e come aggiornarlo dinamicamente
    fish->max_volitive_step = MAX_VOLITIVE_STEP; //TODO: capire qual è il valore migliore di inizializzazione e come aggiornarlo dinamicamente
}

// Funzione per inizializzare un array di pesci
void initFishArray(Fish* fishArray, int n_fishes) {
    for (int i = 0; i < n_fishes; i++) {
        initFish(&fishArray[i]);  // Inizializza ciascun pesce
        // print_fish(fishArray[i]);
    }
}

// Per resettare le variabili all'inizio di ogni epoca
void variablesReset(float *local_tot_fitness, float *local_weighted_tot_fitness, float *local_max_improvement) {
    *local_tot_fitness = 0.0;

    #pragma omp parallel for
    for (int i = 0; i<DIMENSIONS; i++ ){
        local_weighted_tot_fitness[i] = 0.0;
    }
    *local_max_improvement = 0.0;
}

void individualMovement(Fish *fish, float *local_tot_delta_fitness, float *local_weighted_tot_delta_fitness, float *local_max_delta_fitness_improvement) {
    // Movimento casuale per ogni dimensione
    #pragma omp parallel for
    for (int d = 0; d < DIMENSIONS; d++)
    {
        double normalized_movement = (rand() / (double)RAND_MAX) * 2 - 1; //valore qualsiasi tra -1 e 1
        double individual_step = normalized_movement * fish->max_individual_step;

        fish->new_position[d] = fish->position[d] + individual_step;
    }

    // Aggiorno la fitness
    fish->new_fitness = objective_function(fish->new_position)*MULTIPLIER;

    // -------------- Update the collective variables
    double delta_fitness = fish->new_fitness - fish->fitness;

    // CI INTERESSANO SOLO I MOVIMENTI CHE VANNO AD AUMENTARE LA FITNESS DEL PESCE
    if (delta_fitness <= 0) {
        delta_fitness = 0.0;
    }

    //TODO: da rimuovere in favore del calcolo parallelo
    *local_tot_delta_fitness += delta_fitness;

    #pragma omp parallel for
    for (int d = 0; d < DIMENSIONS; d++)
    {
        // La delta_fitness farà in modo che pesci che si muovono "meglio" influenzino il movimento collettivo più degli altri
        local_weighted_tot_delta_fitness[d] += (fish->new_position[d] - fish->position[d]) * delta_fitness;
    }

    if (fabs(delta_fitness) > *local_max_delta_fitness_improvement) {
        *local_max_delta_fitness_improvement = delta_fitness;
    }

    // -------------- Finish update the collective variables

    // Update fish position considering only its individual movement
    #pragma omp parallel for
    for (int d = 0; d < DIMENSIONS; d++)
    {
        if (delta_fitness > 0) {
            // printf("Update for individual movement of %f, because of delta fitness %f  ", fish->new_position[d]-fish->position[d], delta_fitness);
            // printf("new_fitness: %f , old_fitness: %f\n ", fish->new_fitness, fish->fitness);
            fish->position[d] = fish->new_position[d];
        } else {
            // per sicurezza, se la delta_fitness non è positiva, non aggiorniamo la posizione
            fish->new_position[d] = fish->position[d];
        }
    }
}

void individualMovementArray (Fish *fishArray, int n_fishes, float *local_tot_delta_fitness, float *global_tot_delta_fitness, float *local_weighted_tot_delta_fitness, float *global_weighted_tot_delta_fitness, float *local_max_delta_fitness_improvement, float *global_max_delta_fitness_improvement) {
    for (int i = 0; i < n_fishes; i++) {
        individualMovement(&fishArray[i], local_tot_delta_fitness, local_weighted_tot_delta_fitness, local_max_delta_fitness_improvement);  // Inizializza ciascun pesce
    }

    //aggiornamento parallelo
    MPI_Allreduce(local_tot_delta_fitness, global_tot_delta_fitness, 1 , MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_weighted_tot_delta_fitness, global_weighted_tot_delta_fitness, DIMENSIONS, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_max_delta_fitness_improvement, global_max_delta_fitness_improvement, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

}


void collectiveMovement(Fish *fish, float *global_tot_delta_fitness, float *global_weighted_tot_delta_fitness) {
    if (*global_tot_delta_fitness == 0.0) {
        *global_tot_delta_fitness = 1;
    }

    #pragma omp parallel for
    for (int d =0;d<DIMENSIONS; d++){
        fish->new_position[d] = fish->position[d] + global_weighted_tot_delta_fitness[d] / *global_tot_delta_fitness;
        // printf("Update for collective movement of %f\n", fish->new_position[d]-fish->position[d]);
        fish->position[d] = fish->new_position[d]; //TODO: fa schifo, ma segue la logica dell'aggiornare prima la new position e poi quella current
    }
    fish->new_fitness = objective_function(fish->position) * MULTIPLIER; // questo va fatto per forza!
}

void collectiveMovementArray(Fish *fishArray, int n_fishes, float *global_tot_delta_fitness, float *global_weighted_tot_delta_fitness) {
    #pragma omp parallel for
    for (int i = 0; i < n_fishes; i++) {
        collectiveMovement(&fishArray[i], global_tot_delta_fitness, global_weighted_tot_delta_fitness);  // Inizializza ciascun pesce
    }
}

void updateWeights(Fish *fish, float *global_max_delta_fitness_improvement) {
    if (*global_max_delta_fitness_improvement != 0.0) { // Avoid division by zero
        fish->weight += (fish->new_fitness - fish->fitness)/ *global_max_delta_fitness_improvement;
    }    // fish->weight += (fish->new_fitness - fish->fitness);

    if (fish->weight<=W_SCALE_MIN) {
        fish->weight = W_SCALE_MIN; //TODO: non siamo sicure di questa cosa...
    } else if (fish->weight>W_SCALE_MAX) {
        fish->weight = W_SCALE_MAX;
    }

    //che qui la delta fitness sia positiva, non ci interessa...
    //a noi interessa che la delta fitness sia positiva prima di fare il movimento singolo
    fish->fitness = fish->new_fitness;
}

void updateWeightsArray(Fish *fishArray, int n_fishes, float *global_max_delta_fitness_improvement) {
    #pragma omp parallel for
    for (int i = 0; i < n_fishes; i++) {
        updateWeights(&fishArray[i], global_max_delta_fitness_improvement);
        // print_fish(*fishArray[i]);
    }
}


//QUA


void calculateBarycenter(Fish *fishArray, int n_fishes, float *global_barycenter, int rank){

    float local_numerator[DIMENSIONS];
    float global_numerator[DIMENSIONS];
    float local_denominator = 0.0;
    float global_denominator = 0.0;

    for (int d = 0; d<DIMENSIONS; d++){
        global_barycenter[d] = 0.0;
        local_numerator[d] = 0.0;
        global_numerator[d] = 0.0;
    }

    
    for (int i = 0; i < n_fishes; i++) {
        local_denominator += fishArray[i].weight;

        for (int d = 0; d<DIMENSIONS; d++){
            local_numerator[d] += fishArray[i].position[d] * fishArray[i].weight;
        }
    }

    MPI_Allreduce(local_numerator, global_numerator, DIMENSIONS, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_denominator, &global_denominator, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    for (int d = 0; d<DIMENSIONS; d++){
        if (global_denominator != 0.0) {
            global_barycenter[d] = global_numerator[d] / global_denominator;
        }
    }

    printf("[%d] Barycenter: %f %f\n", rank, global_barycenter[0], global_barycenter[1]);
}

// void calculateSumWeights(Fish *fishArray, float *old_sum, float *new_sum){
//     *old_sum = 0.0;
//     *new_sum = 0.0;

//     for (int i = 0; i < N_FISHES; i++) {
//         *old_sum += fishArray[i].previous_cycle_weight;
//         *new_sum += fishArray[i].weight;
//     }
// }

// void volitivePositionUpdateArray(Fish *fishArray, int shrink, float* global_barycenter){
//     double rand_mult = 0.0;

//     // questo codice si può ottimizare mettendo shrink -1,1
//     if (shrink==1) {
//         for (int i = 0; i < N_FISHES; i++) {
//             for (int d = 0; d < DIMENSIONS; d++) { // TODO: change max individual step with another step in order to have the possibility to tune it
//                 rand_mult= fmin(((double)rand() / (double)RAND_MAX) + 0.1, 1.0); //valore qualsiasi tra 0.1 e 1

//                 double temp= fishArray[i].position[d];
//                 fishArray[i].position[d] -= fishArray[i].max_volitive_step * rand_mult * (fishArray[i].position[d] - global_barycenter[d]);

//                 if (fishArray[i].position[d] > 1000.0 || fishArray[i].position[d] < -1000.0) {
//                     // printf("LAST STRANGE FISH: rand_mult: %f\n", rand_mult);
//                     // printf("dim= %d, pesce= %d\n", d, i);
//                     // printf("fishArray[%d].max_volitive_step: %f\n",i, fishArray[i].max_volitive_step);
//                     // printf("fishArray[%d].position[%d] before update: %f\n",i, d, temp);
//                     // printf("fishArray[%d].position[%d] after update: %f\n",i, d, fishArray[i].position[d]);
//                     // printf("global_barycenter: %f\n", global_barycenter[d]);
//                     printf("porco cane\n");
//                     exit(1);
//                 }
//             }
//         }
//     } else {
//         for (int i = 0; i < N_FISHES; i++) {
//             for (int d = 0; d < DIMENSIONS; d++){ // TODO: change max individual step with another step in order to have the possibility to tune it
//                 rand_mult= fmin(((double)rand() / (double)RAND_MAX) + 0.1, 1.0); //valore qualsiasi tra 0.1 e 1

//                 double temp= fishArray[i].position[d];
//                 fishArray[i].position[d] += fishArray[i].max_volitive_step * rand_mult * (fishArray[i].position[d] - global_barycenter[d]);

//                 if (fishArray[i].position[d] > 1000.0 || fishArray[i].position[d] < -1000.0) {
//                     // printf("LAST STRANGE FISH: rand_mult: %f\n", rand_mult);
//                     // printf("dim= %d, pesce= %d\n", d, i);
//                     // printf("fishArray[%d].max_volitive_step: %f\n",i, fishArray[i].max_volitive_step);
//                     // printf("fishArray[%d].position[%d] before update: %f\n",i, d, temp);
//                     // printf("fishArray[%d].position[%d] after update: %f\n",i, d, fishArray[i].position[d]);
//                     // printf("global_barycenter: %f\n", global_barycenter[d]);
//                     exit(1);
//                 }
//             }
//         }
//     }
// }

void collectiveVolitiveArray(Fish *fishes, int n_fishes, int rank) {
    float global_barycenter[DIMENSIONS];
    calculateBarycenter(fishes, n_fishes, global_barycenter, rank);

    // float old_sum_weights;
    // float new_sum_weights;
    // calculateSumWeights(fishes, &old_sum_weights, &new_sum_weights);

    // if (old_sum_weights < new_sum_weights) {
    //     //shrink = 1 -> il banco ha guadagnato peso quindi si deve avvicinare al baricentro
    //     // printf("MOVIMENTO YEAH -> IL BANCO SI È AVVICINATO\n");
    //     volitivePositionUpdateArray(fishes, 1, global_barycenter);
    // } else if (old_sum_weights > new_sum_weights) {
    //     //TODO: shrink = 0 -> il banco ha perso peso quindi si deve allargare in cerca di cibo
    //     // printf("MOVIMENTO BLEAH -> IL BANCO SI È ALLONTANATO\n");
    //     volitivePositionUpdateArray(fishes, 0, global_barycenter);
    // }else{
    //     // printf("EQUAL WEIGHTS, do nothing");
    // }

    // // for (int i = 0; i < N_FISHES; i++) {
    // //     print_fish(fishes[i]);
    // // }

    // // update previous_cycle_weight
    // for (int i = 0; i < N_FISHES; i++) {
    //     fishes[i].previous_cycle_weight = fishes[i].weight;
    // }
}

//-------------------------------------------------------------------------------------------
//---------------------------- MAIN ---------------------------------------------------------
//-------------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {

    //variabili MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //TODO: pensiamo ad un modo per implementare un timer

    //variabili locali al sottogruppo di pesci
    float local_best_fitness = -2000.0;
    float global_best_fitness = -2000.0;
    float local_total_fitness = 0.0;
    float global_total_fitness = 0.0;
    float local_weighted_total_fitness[DIMENSIONS];
    float global_weighted_total_fitness[DIMENSIONS];
    float local_max_improvement = 0.0;
    float global_max_improvement = 0.0;
    srand(time(NULL)+rank);  // Seed for random number generation have a different value for each process at each iteration

    //l'effettivo sottogruppo di pesci
    int local_n = N_FISHES / size;
    Fish *local_school = malloc(local_n * sizeof(Fish));
    initFishArray(local_school, local_n);

    for (int iter = 0; iter < MAX_ITER; iter++) {
        variablesReset(&local_total_fitness, local_weighted_total_fitness, &local_max_improvement);
        individualMovementArray(local_school, local_n, &local_total_fitness, &global_total_fitness, local_weighted_total_fitness, global_weighted_total_fitness, &local_max_improvement, &global_max_improvement);
        updateWeightsArray(local_school, local_n, &global_max_improvement);
        collectiveMovementArray(local_school, local_n, &global_total_fitness, global_weighted_total_fitness);
        collectiveVolitiveArray(local_school, local_n, rank);
    }

    free(local_school);
    MPI_Finalize();
    return 0;
}
