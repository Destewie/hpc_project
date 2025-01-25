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

#define N_FISHES 15
#define DIMENSIONS 2
#define MAX_ITER 150
#define BOUNDS_MIN -30.0   // Minimum bound of the search space
#define BOUNDS_MAX 30.0    // Maximum bound of the search space
#define BOUNDS_MIN_W 0.1   // Minimum bound of the search space
#define BOUNDS_MAX_W 10.0    // Maximum bound of the search space
#define MAX_INDIVIDUAL_STEP 1.5 // Maximum step for individual movement
#define MAX_VOLITIVE_STEP 0.2 // Maximum step for collective movement
#define W_SCALE_MIN 1.0
#define W_SCALE_MAX 10.0
#define BREEDING_THRESHOLD 7.0 // minimus threshold of weight to breedh new fishes
#define FUNCTION "min_sphere"   //TODO: Capire se, al posto di fare un controllo su una stringa, possiamo passare alle funzioni direttamente un puntatore ad una funzione (in modo comodo, se no lasciamo perdere)
#define MULTIPLIER -1   // 1 in case of maximization, -1 in case of minimization
#define A 10.0 //rastrigin param
#define LOG 1 // 1 to log in the json, 0 not to log

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

void write_fishes_to_json(Fish *fishes, FILE *file, int first, int last) {
    if (first) {
        // Scrive l'apertura dell'array principale solo se è la prima chiamata
        fprintf(file, "[\n");
    }

    fprintf(file, "\t[\n");

    for (int i = 0; i < N_FISHES; i++) {
        if (DIMENSIONS == 1) {
            fprintf(file, "\t\t{\"x\": [%.6f],", fishes[i].position[0]);
        } else if (DIMENSIONS == 2) {
            fprintf(file, "\t\t{\"x\": [%.6f, %.6f],", fishes[i].position[0], fishes[i].position[1]);
        } else {
            fprintf(file, "\t\t{\"x\": [");
            for (int d = 0; d < DIMENSIONS; d++) {
                fprintf(file, "%.6f", fishes[i].position[d]);
                if (d < DIMENSIONS - 1) {
                    fprintf(file, ", ");
                }
            }
            fprintf(file, "],");
        }
        fprintf(file, "\"weight\": %.6f}", fishes[i].weight);

        if (i < N_FISHES - 1) {
            fprintf(file, ",\n");
        } else {
            fprintf(file, "\n");
        }
    }

    if (last) {
        // Chiude l'array principale se è l'ultima chiamata
        fprintf(file, "\t]\n");
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
void variables_reset(float *tot_fitness, float *weighted_tot_fitness, float *max_improvement) {
    *tot_fitness = 0.0;

    for (int i = 0; i<DIMENSIONS; i++ ){
        weighted_tot_fitness[i] = 0.0;
    }
    *max_improvement = 0.0;
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
//---------------------------- FISH ---------------------------------------------------------
//-------------------------------------------------------------------------------------------



void print_fish(Fish fish){
    printf("Fish: ");
    // for(int i=0; i<DIMENSIONS; i++){
    //     printf("pos: %f ", fish.position[i]);
    // }
    printf("\tweight: %f \tfitness: %f\n", fish.weight, fish.fitness);
}

// Funzione per inizializzare un singolo pesce
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
void initFishArray(Fish* fishArray) {
    for (int i = 0; i < N_FISHES; i++) {
        initFish(&fishArray[i]);  // Inizializza ciascun pesce
        // print_fish(fishArray[i]);
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
    fish->new_fitness = objective_function(fish->new_position)*MULTIPLIER;



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


void individualMovementArray (Fish *fishArray, float *tot_delta_fitness, float *weighted_tot_delta_fitness, float *max_delta_fitness_improvement) {
    for (int i = 0; i < N_FISHES; i++) {
        individualMovement(&fishArray[i], tot_delta_fitness, weighted_tot_delta_fitness, max_delta_fitness_improvement);  // Inizializza ciascun pesce
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
    fish->new_fitness = objective_function(fish->position) * MULTIPLIER; // questo va fatto per forza!
}

void collectiveMovementArray(Fish *fishArray, float *tot_delta_fitness, float *weighted_tot_delta_fitness) {
    for (int i = 0; i < N_FISHES; i++) {
        collectiveMovement(&fishArray[i], tot_delta_fitness, weighted_tot_delta_fitness);  // Inizializza ciascun pesce
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
    for (int i = 0; i < N_FISHES; i++) {
        updateWeights(&fishArray[i], max_delta_fitness_improvement);
        // print_fish(fishArray[i]);
    }
}

void calculateBarycenter(Fish *fishArray, float *barycenter){
    float numerator[DIMENSIONS];
    float denominator[DIMENSIONS];

    for (int d = 0; d<DIMENSIONS; d++){
        barycenter[d] = 0.0;
        numerator[d] = 0.0;
        denominator[d] = 0.0;
    }

    for (int d = 0; d<DIMENSIONS; d++){
        for (int i = 0; i < N_FISHES; i++) {
            numerator[d] += fishArray[i].position[d] * fishArray[i].weight;
            denominator[d] += fishArray[i].weight;
        }

        if (denominator[d] != 0.0) {
            barycenter[d] = numerator[d] / denominator[d];
        }
    }
}

void calculateSumWeights(Fish *fishArray, float *old_sum, float *new_sum){
    *old_sum = 0.0;
    *new_sum = 0.0;

    for (int i = 0; i < N_FISHES; i++) {
        *old_sum += fishArray[i].previous_cycle_weight;
        *new_sum += fishArray[i].weight;
    }
}

void volitivePositionUpdateArray(Fish *fishArray, int shrink, float* barycenter){
    double rand_mult = 0.0;

    // questo codice si può ottimizare mettendo shrink -1,1
    if (shrink==1) {
        for (int i = 0; i < N_FISHES; i++) {
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
                    printf("porco cane\n");
                    exit(1);
                }
            }
        }
    } else {
        for (int i = 0; i < N_FISHES; i++) {
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
    float barycenter[DIMENSIONS];
    calculateBarycenter(fishes, barycenter);

    float old_sum_weights;
    float new_sum_weights;
    calculateSumWeights(fishes, &old_sum_weights, &new_sum_weights);

    if (old_sum_weights < new_sum_weights) {
        //shrink = 1 -> il banco ha guadagnato peso quindi si deve avvicinare al baricentro
        volitivePositionUpdateArray(fishes, 1, barycenter);
    } else if (old_sum_weights > new_sum_weights) {
        //shrink = 0 -> il banco ha perso peso quindi si deve allargare in cerca di cibo
        volitivePositionUpdateArray(fishes, 0, barycenter);
    }else{
        // printf("EQUAL WEIGHTS, do nothing");
    }

    // update previous_cycle_weight
    for (int i = 0; i < N_FISHES; i++) {
        fishes[i].previous_cycle_weight = fishes[i].weight;
    }
}

void breeding(Fish *fishes){
    //mi salvo l'indice del pesce con weight maggiore e il secondo in ordine di weight
    //se la weigh supera la threshold costante, allora figlio -> creo un nuovo pesce con posizione media tra i due e peso medio tra i due
    int first_index = 0;
    int second_index = 0;
    int worst_index = 0;

    for (int i = 0; i < N_FISHES; i++) {
        if (fishes[i].weight > fishes[first_index].weight) {
            second_index = first_index;
            first_index = i;
        } else if (fishes[i].weight > fishes[second_index].weight) {
            second_index = i;
        }

        if (fishes[i].weight < fishes[worst_index].weight) {
            worst_index = i;
        }
    }

    //sopprimo il pesce peggiore e lo rimpiazzo con il figlio dei due migliori
    if (fishes[first_index].weight > BREEDING_THRESHOLD && fishes[second_index].weight > BREEDING_THRESHOLD) {
        for (int d = 0; d < DIMENSIONS; d++) {
            fishes[worst_index].position[d] = (fishes[first_index].position[d] + fishes[second_index].position[d]) / 2;
        }
        fishes[worst_index].weight = (fishes[first_index].weight + fishes[second_index].weight) / 2;
        fishes[worst_index].fitness = objective_function(fishes[worst_index].position)*MULTIPLIER;
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

    FILE *file; 

    // File opening
    if (DIMENSIONS <= 2 && LOG) {
        char filename[50];
        sprintf(filename, "../../evolution_logs/%s_%dd_log.json",FUNCTION, DIMENSIONS);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return 1;
        }
    }

    float best_fitness = -2000.0;
    float total_fitness = 0.0;
    float weighted_total_fitness[DIMENSIONS];
    float max_improvement = 0.0;
    srand(time(NULL));  // Seed for random number generation

    // INITIALIZATION
    Fish fishes[N_FISHES];
    initFishArray(fishes);

    if (DIMENSIONS <= 2 && LOG) {
        write_fishes_to_json(fishes, file, 1, 0);
    }

    // MAIN LOOP
    for (int iter = 0; iter < MAX_ITER; iter++) {

        variables_reset(&total_fitness, weighted_total_fitness, &max_improvement);

        // INDIVIDUAL MOVEMENT
        individualMovementArray(fishes, &total_fitness, weighted_total_fitness, &max_improvement);

        // UPDATE WEIGHTS
        updateWeightsArray(fishes, &max_improvement);

        // COLLECTIVE MOVEMENT
        collectiveMovementArray(fishes, &total_fitness, weighted_total_fitness);

        // COLLECTIVE VOLITIVE MOVEMENT
        collectiveVolitiveArray(fishes);

        // BREEDING
        breeding(fishes);

        // SAVE ON FILE
        if (DIMENSIONS <= 2 && LOG) {
            write_fishes_to_json(fishes, file, 0, iter==MAX_ITER-1?1:0);
        }

        for (int i = 0; i < N_FISHES; i++) {
            if (best_fitness<fishes[i].fitness) {
                best_fitness = fishes[i].fitness;
            }
        }
    }


    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    fprintf(file, "\n]");
    fclose(file);


    // print all the fishes
    for (int i = 0; i < N_FISHES; i++) {
        print_fish(fishes[i]);
    }

    printf("Number of fishes: %d\n", N_FISHES);
    printf("Dimensions: %d\n", DIMENSIONS);
    printf("Epochs: %d\n", MAX_ITER);
    printf("Execution time: %f\n", cpu_time_used);
    printf("Best fitness: %f\n", best_fitness/DIMENSIONS);

    return 0;
}
