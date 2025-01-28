#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#include <time.h>
#include <string.h> // Include for strcmp

#define N_SCHOOLS 3
#define N_FISHES_PER_SCHOOL 15
#define N_FISHES (N_SCHOOLS*N_FISHES_PER_SCHOOL)
#define DIMENSIONS 2
#define MAX_ITER 100
#define BOUNDS_MIN 30.0   // Minimum bound of the search space
#define BOUNDS_MAX -30.0    // Maximum bound of the search space
#define BOUNDS_MIN_W 0.1   // Minimum bound of the search space
#define BOUNDS_MAX_W 10.0    // Maximum bound of the search space
#define MAX_INDIVIDUAL_STEP 1.7 // Maximum step for individual movement
#define MAX_VOLITIVE_STEP 0.2 // Maximum step for collective movement
#define W_SCALE_MIN 1.0
#define W_SCALE_MAX 10.0
#define BREEDING_THRESHOLD 7.0 // minimus threshold of weight to breedh new fishes
#define FUNCTION "min_ackley"   //TODO: Capire se, al posto di fare un controllo su una stringa, possiamo passare alle funzioni direttamente un puntatore ad una funzione (in modo comodo, se no lasciamo perdere)
#define MULTIPLIER -1   // 1 in case of maximization, -1 in case of minimization
#define A 10.0 //rastrigin param
#define LOG 1 // 1 to log the results, 0 otherwise

//10 very different colors that will be used by a python script to plot the results
const char *COLORS[] = {"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"};

typedef struct {
    double position[DIMENSIONS];
    double new_position[DIMENSIONS];

    double previous_cycle_weight;
    double weight;

    double fitness;
    double new_fitness;

    double max_individual_step;
    double max_volitive_step;
}Fish;

//-------------------------------------------------------------------------------------------
//----------------------------- UTILS -------------------------------------------------------
//-------------------------------------------------------------------------------------------

void WriteFishesToJson(Fish *fishes, FILE *file, int first, int last) {

    if (first) {
        // Scrive l'apertura dell'array principale solo se è la prima chiamata
        fprintf(file, "[\n");
    }

    fprintf(file, "\t[\n");

    for (int i = 0; i < N_FISHES; i++) {
        int school_index = i / N_FISHES_PER_SCHOOL;
        fprintf(file, "\t\t{\"x\": [");
        for (int d = 0; d < DIMENSIONS; d++) {
            fprintf(file, "%.6f", fishes[i].position[d]);
            if (d < DIMENSIONS - 1) {
                fprintf(file, ", ");
            }
        }
        fprintf(file, "], \"weight\": %.6f, \"color\": \"%s\"}", fishes[i].weight, COLORS[school_index]);

        if (i < N_FISHES - 1) {
            fprintf(file, ",\n");
        } else {
            fprintf(file, "\n");
        }
    }

    if (last) {
        // Chiude l'array principale se è l'ultima chiamata
        fprintf(file, "\t]\n]\n");
    } else {
        fprintf(file, "\t],\n");
    }
}

// Clamp positions to within bounds
double clamp(double value, double min, double max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

// Per resettare le variabili all'inizio di ogni epoca
void variablesReset(float *tot_fitness, float weighted_tot_fitness[N_SCHOOLS][DIMENSIONS], float *max_improvement) {
    for (int i = 0; i < N_SCHOOLS; i++) {
        tot_fitness[i] = 0.0;

        for (int d = 0; d<DIMENSIONS; d++ ){
            weighted_tot_fitness[i][d] = 0.0;
        }
        max_improvement[i] = 0.0;
    }
}

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

double objectiveFunction(double *x) {
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
//---------------------------- FISH ---------------------------------------------------------
//-------------------------------------------------------------------------------------------



void printFish(Fish fish){
    printf("Fish: ");
    for(int i=0; i<DIMENSIONS; i++){
        printf("pos: %f - ", fish.position[i]);
    }
    printf("\tweight: %f \tfitness: %f\n", fish.weight, fish.fitness);
}

// Funzione per inizializzare un singolo pesce
void initFish(Fish *fish, int fish_index) {

    // Posizioni iniziali divise per banco
    int school_index = fish_index / N_FISHES_PER_SCHOOL; // Calcola il banco
    double portion_bounds = (BOUNDS_MAX - BOUNDS_MIN) / N_SCHOOLS; // Calcola la porzione corretta per ciascun banco
    double lower_bound = BOUNDS_MIN + school_index * portion_bounds; // Limite inferiore per il banco
    double upper_bound = lower_bound + portion_bounds; // Limite superiore per il banco

    for (int d = 0; d < DIMENSIONS; d++) {
        // // Posizioni iniziali random
        // fish->position[d] = ((double)rand() / RAND_MAX) * (BOUNDS_MAX - BOUNDS_MIN) + BOUNDS_MIN;

        // Posizioni iniziali divise per banco
        if (d == 0) {
            fish->position[d] = ((double)rand() / RAND_MAX) * (upper_bound - lower_bound) + lower_bound;
            printf("[D0] lower_bound: %f, upper_bound: %f\n, x: %f", lower_bound, upper_bound, fish->position[d]);
        } else {
            fish->position[d] = ((double)rand() / RAND_MAX) * (BOUNDS_MAX - BOUNDS_MIN) + BOUNDS_MIN;
        }

        fish->new_position[d] = fish->position[d];
    }

    fish->weight = W_SCALE_MAX / 2;   // Peso iniziale
    fish->previous_cycle_weight = fish->weight;

    fish->fitness = objectiveFunction(fish->position)*MULTIPLIER;        // Fitness iniziale //TODO: capire qual è il valore migliore di inizializzazione
    fish->new_fitness = fish->fitness;

    fish->max_individual_step = MAX_INDIVIDUAL_STEP; //TODO: capire qual è il valore migliore di inizializzazione e come aggiornarlo dinamicamente
    fish->max_volitive_step = MAX_VOLITIVE_STEP; //TODO: capire qual è il valore migliore di inizializzazione e come aggiornarlo dinamicamente
}

// Funzione per inizializzare un array di pesci
void initFishArray(Fish* fishArray) {
    for (int i = 0; i < N_FISHES; i++) {
        initFish(&fishArray[i], i);  // Inizializza ciascun pesce
        // printFish(fishArray[i]);
    }
}

// Movimento individuale
void individualMovement(Fish *fish, float *tot_delta_fitness, float *weighted_tot_delta_fitness, float *max_delta_fitness_improvement) {

    // Movimento casuale per ogni dimensione
    for (int d = 0; d < DIMENSIONS; d++)
    {
        double normalized_movement = (rand() / (double)RAND_MAX) * 2 - 1; //valore qualsiasi tra -1 e 1
        double individual_step = normalized_movement * fish->max_individual_step;

        fish->new_position[d] = fish->position[d] + individual_step;
    }

    // Aggiorno la fitness
    fish->new_fitness = objectiveFunction(fish->new_position)*MULTIPLIER;



    // -------------- Update the collective variables
    double delta_fitness = fish->new_fitness - fish->fitness;

    // CI INTERESSANO SOLO I MOVIMENTI CHE VANNO AD AUMENTARE LA FITNESS DEL PESCE
    if (delta_fitness > 0) {

        *tot_delta_fitness += delta_fitness;

        for (int d = 0; d < DIMENSIONS; d++)
        {
            // La delta_fitness farà in modo che pesci che si muovono "meglio" influenzino il movimento collettivo più degli altri
            weighted_tot_delta_fitness[d] += (fish->new_position[d] - fish->position[d]) * delta_fitness;
        }


        if (fabs(delta_fitness) > *max_delta_fitness_improvement) {
            *max_delta_fitness_improvement = delta_fitness;
        }
    }
    // -------------- Finish update the collective variables

    // Update fish position considering only its individual movement
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


void individualMovementArray (Fish *fishArray, float *tot_delta_fitness, float weighted_tot_delta_fitness[N_SCHOOLS][DIMENSIONS], float *max_delta_fitness_improvement) {
    for (int s = 0; s < N_SCHOOLS; s++) {
        for (int i = 0; i < N_FISHES_PER_SCHOOL; i++) {
            individualMovement(&fishArray[s*N_FISHES_PER_SCHOOL+i], &tot_delta_fitness[s], weighted_tot_delta_fitness[s], &max_delta_fitness_improvement[s]);  // Inizializza ciascun pesce
        }
    }
}


void updateWeights(Fish *fish, float *max_delta_fitness_improvement) {
    if (*max_delta_fitness_improvement != 0.0) { // Avoid division by zero
        fish->weight += (fish->new_fitness - fish->fitness)/ *max_delta_fitness_improvement;
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

void updateWeightsArray(Fish *fishArray,  float *max_delta_fitness_improvement) {
    for (int s = 0; s < N_SCHOOLS; s++) {
        for (int i = 0; i < N_FISHES_PER_SCHOOL; i++) {
            updateWeights(&fishArray[s*N_FISHES_PER_SCHOOL+i], &max_delta_fitness_improvement[s]);
        }
    }
}


void collectiveMovement(Fish *fish, float *tot_delta_fitness, float *weighted_tot_delta_fitness) {
    if (*tot_delta_fitness == 0.0) {
        *tot_delta_fitness = 1;
    }

    for (int d =0;d<DIMENSIONS; d++){
        fish->new_position[d] = fish->position[d] + weighted_tot_delta_fitness[d] / *tot_delta_fitness;
        // printf("Update for collective movement of %f\n", fish->new_position[d]-fish->position[d]);
        fish->position[d] = fish->new_position[d]; //TODO: fa schifo, ma segue la logica dell'aggiornare prima la new position e poi quella current
    }
    fish->new_fitness = objectiveFunction(fish->position) * MULTIPLIER; // questo va fatto per forza!
}

void collectiveMovementArray(Fish *fishArray, float *tot_delta_fitness, float weighted_tot_delta_fitness[N_SCHOOLS][DIMENSIONS]) {
    for (int s = 0; s < N_SCHOOLS; s++) {
        for (int i = 0; i < N_FISHES_PER_SCHOOL; i++) {
            collectiveMovement(&fishArray[s*N_FISHES_PER_SCHOOL+i], &tot_delta_fitness[s], weighted_tot_delta_fitness[s]);  // Inizializza ciascun pesce
        }
    }
}

void calculateBarycenters(Fish *fishArray, float barycenter[N_SCHOOLS][DIMENSIONS]){
    float numerator[N_SCHOOLS][DIMENSIONS];
    float denominator[N_SCHOOLS][DIMENSIONS];

    for (int s = 0; s<N_SCHOOLS; s++) {
        for (int d = 0; d<DIMENSIONS; d++){
            numerator[s][d] = 0.0;
            denominator[s][d] = 0.0;
        }
    }

    for (int s = 0; s < N_SCHOOLS; s++) {
        for (int d = 0; d < DIMENSIONS; d++) {
            for (int i = 0; i < N_FISHES_PER_SCHOOL; i++) {
                numerator[s][d] += fishArray[s * N_FISHES_PER_SCHOOL + i].position[d] * fishArray[s * N_FISHES_PER_SCHOOL + i].weight;
                denominator[s][d] += fishArray[s * N_FISHES_PER_SCHOOL + i].weight;
            }

            if (denominator[s][d] != 0.0) {
                barycenter[s][d] = numerator[s][d] / denominator[s][d];
            }
        }
    }
}

void calculateSumWeights(Fish *fishArray, float *old_sum, float *new_sum){
    for (int s=0; s<N_SCHOOLS; s++) {
        old_sum[s] = 0.0;
        new_sum[s] = 0.0;
    }

    for (int s=0; s<N_SCHOOLS; s++) {
        for (int i = 0; i < N_FISHES_PER_SCHOOL; i++) {
            old_sum[s] += fishArray[s * N_FISHES_PER_SCHOOL + i].previous_cycle_weight;
            new_sum[s] += fishArray[s * N_FISHES_PER_SCHOOL + i].weight;
        }
    }
}

void volitivePositionUpdateArray(Fish *fishArray, int school_index, int shrink, float *barycenter){
    double rand_mult = 0.0;

    // questo codice si può ottimizare mettendo shrink -1,1
    if (shrink==1) {
        for (int i = school_index*N_FISHES_PER_SCHOOL; i < school_index*N_FISHES_PER_SCHOOL + N_FISHES_PER_SCHOOL; i++) {
            for (int d = 0; d < DIMENSIONS; d++) { // TODO: change max individual step with another step in order to have the possibility to tune it
                rand_mult= fmin(((double)rand() / (double)RAND_MAX) + 0.1, 1.0); //valore qualsiasi tra 0.1 e 1

                double temp= fishArray[i].position[d];
                fishArray[i].position[d] -= fishArray[i].max_volitive_step * rand_mult * (fishArray[i].position[d] - barycenter[d]);

                if (fishArray[i].position[d] > 1000.0 || fishArray[i].position[d] < -1000.0) {
                    // printf("LAST STRANGE FISH: rand_mult: %f\n", rand_mult);
                    // printf("dim= %d, pesce= %d\n", d, i);
                    // printf("fishArray[%d].max_volitive_step: %f\n",i, fishArray[i].max_volitive_step);
                    // printf("fishArray[%d].position[%d] before update: %f\n",i, d, temp);
                    // printf("fishArray[%d].position[%d] after update: %f\n",i, d, fishArray[i].position[d]);
                    // printf("barycenter: %f\n", barycenter[d]);
                    printf("Error: Fish position out of bounds.\n");
                    exit(1);
                }
            }
        }
    } else {
        for (int i = school_index*N_FISHES_PER_SCHOOL; i < school_index*N_FISHES_PER_SCHOOL + N_FISHES_PER_SCHOOL; i++) {
            for (int d = 0; d < DIMENSIONS; d++){ // TODO: change max individual step with another step in order to have the possibility to tune it
                rand_mult= fmin(((double)rand() / (double)RAND_MAX) + 0.1, 1.0); //valore qualsiasi tra 0.1 e 1

                double temp= fishArray[i].position[d];
                fishArray[i].position[d] += fishArray[i].max_volitive_step * rand_mult * (fishArray[i].position[d] - barycenter[d]);

                if (fishArray[i].position[d] > 1000.0 || fishArray[i].position[d] < -1000.0) {
                    // printf("LAST STRANGE FISH: rand_mult: %f\n", rand_mult);
                    // printf("dim= %d, pesce= %d\n", d, i);
                    // printf("fishArray[%d].max_volitive_step: %f\n",i, fishArray[i].max_volitive_step);
                    // printf("fishArray[%d].position[%d] before update: %f\n",i, d, temp);
                    // printf("fishArray[%d].position[%d] after update: %f\n",i, d, fishArray[i].position[d]);
                    // printf("barycenter: %f\n", barycenter[d]);
                    exit(1);
                }
            }
        }
    }
}

void collectiveVolitiveArray(Fish *fishes) {
    float barycenter[N_SCHOOLS][DIMENSIONS];
    calculateBarycenters(fishes, barycenter);

    float old_sum_weights[N_SCHOOLS];
    float new_sum_weights[N_SCHOOLS];
    calculateSumWeights(fishes, old_sum_weights, new_sum_weights);

    for (int s = 0; s < N_SCHOOLS; s++) {
        if (old_sum_weights[s] < new_sum_weights[s]) {
            //shrink = 1 -> il banco ha guadagnato peso quindi si deve avvicinare al baricentro
            volitivePositionUpdateArray(fishes, s, 1, barycenter[s]);
        } else if (old_sum_weights[s] > new_sum_weights[s]) {
            //shrink = 0 -> il banco ha perso peso quindi si deve allargare in cerca di cibo
            volitivePositionUpdateArray(fishes, s, 0, barycenter[s]);
        }else{
            // printf("EQUAL WEIGHTS, do nothing");
        }

        // update previous_cycle_weight
        for (int i = s * N_FISHES_PER_SCHOOL; i < (s + 1) * N_FISHES_PER_SCHOOL; i++) {
            fishes[i].previous_cycle_weight = fishes[i].weight;
        }
    }
}

void breeding(Fish *fishes) {
    for (int s = 0; s < N_SCHOOLS; s++) {
        // Indici del miglior, secondo miglior e peggiore pesce del banco
        int first_index = s * N_FISHES_PER_SCHOOL;
        int second_index = -1; // Inizializza con un valore non valido
        int worst_index = s * N_FISHES_PER_SCHOOL;

        for (int i = 1; i < N_FISHES_PER_SCHOOL; i++) {
            int current_index = s * N_FISHES_PER_SCHOOL + i;

            // Controlla se il peso corrente è maggiore del primo
            if (fishes[current_index].weight > fishes[first_index].weight) {
                second_index = first_index; // Aggiorna il secondo con il vecchio primo
                first_index = current_index;
            } 
            // Controlla se il peso corrente è maggiore del secondo (solo se valido)
            else if (second_index == -1 || fishes[current_index].weight > fishes[second_index].weight) {
                second_index = current_index;
            }

            // Aggiorna il peggiore
            if (fishes[current_index].weight < fishes[worst_index].weight) {
                worst_index = current_index;
            }
        }

        // Esegui il "breeding" solo se entrambi superano la threshold
        if (fishes[first_index].weight > BREEDING_THRESHOLD && fishes[second_index].weight > BREEDING_THRESHOLD) {
            for (int d = 0; d < DIMENSIONS; d++) {
                fishes[worst_index].position[d] = (fishes[first_index].position[d] + fishes[second_index].position[d]) / 2;
            }
            fishes[worst_index].weight = (fishes[first_index].weight + fishes[second_index].weight) / 2;
            fishes[worst_index].fitness = objectiveFunction(fishes[worst_index].position) * MULTIPLIER;

        }
    }
}


//-------------------------------------------------------------------------------------------
//------------------------------- MAIN ------------------------------------------------------
//-------------------------------------------------------------------------------------------

int main() {

    //create a timer
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    char filename[50];
    sprintf(filename, "../../evolution_logs/%s_%dd_log.json",FUNCTION, DIMENSIONS);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // float best_fitness[N_SCHOOLS];
    float total_fitness[N_SCHOOLS];
    float weighted_total_fitness[N_SCHOOLS][DIMENSIONS];
    float max_improvement[N_SCHOOLS];
    srand(time(NULL));  // Seed for random number generation

    // INITIALIZATION
    Fish fishes[N_FISHES]; //creiamo un vettore unico che sarà diviso in banchi di pesci in base agli indici
    initFishArray(fishes);
    if (DIMENSIONS <= 2 && LOG) {
        WriteFishesToJson(fishes, file, 1, 0);
    }

    // MAIN LOOP
    for (int iter = 0; iter < MAX_ITER; iter++) {

        variablesReset(total_fitness, weighted_total_fitness, max_improvement);

        // INDIVIDUAL MOVEMENT
        individualMovementArray(fishes, total_fitness, weighted_total_fitness, max_improvement);

        // UPDATE WEIGHTS
        updateWeightsArray(fishes, max_improvement);

        // COLLECTIVE MOVEMENT
        collectiveMovementArray(fishes, total_fitness, weighted_total_fitness);

        // COLLECTIVE VOLITIVE MOVEMENT
        collectiveVolitiveArray(fishes);

        // BREEDING
        breeding(fishes);

        // SAVE ON FILE
        if (DIMENSIONS <= 2 && LOG) {
            WriteFishesToJson(fishes, file, 0, iter==MAX_ITER-1?1:0);
        }

        // //calcolo la best fitness 
        // for (int i = 0; i < N_FISHES; i++) {
        //     // printFish(fishes[i]);
        //     if (best_fitness<fishes[i].fitness) {
        //         best_fitness = fishes[i].fitness;
        //     }
        // }
    }


    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    fclose(file);


    // print all the fishes
    for (int s = 0; s < N_SCHOOLS; s++) {
        printf("----------- School %d\n", s);
        for (int i = 0; i < N_FISHES_PER_SCHOOL; i++) {
            printFish(fishes[s*N_FISHES_PER_SCHOOL+i]);
        }
    }

    printf("Number of schools: %d\n", N_SCHOOLS);
    printf("Number of fishes per school: %d\n", N_FISHES_PER_SCHOOL);
    printf("Number of fishes: %d\n", N_FISHES);
    printf("Dimensions: %d\n", DIMENSIONS);
    printf("Epochs: %d\n", MAX_ITER);
    printf("Execution time: %f\n", cpu_time_used);
    // printf("Best fitness: %f\n", best_fitness/DIMENSIONS);

    return 0;
}
