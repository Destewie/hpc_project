#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> // Include for strcmp

#define DIMENSIONS 2
#define N_FISHES 10
#define MAX_ITER 100
#define BOUNDS_MIN -10.0   // Minimum bound of the search space
#define BOUNDS_MAX 10.0    // Maximum bound of the search space
#define BOUNDS_MIN_W 0.1   // Minimum bound of the search space
#define BOUNDS_MAX_W 10.0    // Maximum bound of the search space
#define W_SCALE 10.0
#define FUNCTION "min_sphere"   //TODO: Capire se, al posto di fare un controllo su una stringa, possiamo passare alle funzioni direttamente un puntatore ad una funzione (in modo comodo, se no lasciamo perdere)
#define MULTIPLIER -1   // 1 in case of maximization, -1 in case of minimization
#define A 10.0

typedef struct {
    double position[DIMENSIONS];
    double new_position[DIMENSIONS];
    
    double weight;

    double fitness;
    double new_fitness;

    double max_individual_step;
}Fish;

//-------------------------------------------------------------------------------------------
//----------------------------- UTILS -------------------------------------------------------
//-------------------------------------------------------------------------------------------

void write_fishes_to_json(Fish *fishes,FILE *file, int last) {

    fprintf(file, "\t[\n");

    for (int i = 0; i < N_FISHES; i++) {
        if (DIMENSIONS==1){
            fprintf(file, "\t\t{\"x\": [%.6f ],", fishes[i].position[0]);
        }else if(DIMENSIONS==2){
            fprintf(file, "\t\t{\"x\": [%.6f , %.6f],", fishes[i].position[0],fishes[i].position[1]);
        }
        fprintf(file, "\"weight\": %.6f }", fishes[i].weight);
        if (i == N_FISHES - 1 && last==0) {
            fprintf(file, "\n],");
        } else if(i == N_FISHES - 1 && last==1){
            fprintf(file, "\n]");
        }else{
            fprintf(file, ",\n");
        }
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

double objective_function(double *x) {
    if (strcmp(FUNCTION, "rosenbrok") == 0) {
        return rosenbrok(x);
    } else if (strcmp(FUNCTION, "min_sphere") == 0) {
        return min_sphere(x);
    } else if (strcmp(FUNCTION, "max_sphere") == 0) {
        return max_sphere(x);
    } else if (strcmp(FUNCTION, "min_rastrigin") == 0) {
        return max_sphere(x);
    } else {
        return 0.0;
    }
}

//-------------------------------------------------------------------------------------------
//---------------------------- FISH ---------------------------------------------------------
//-------------------------------------------------------------------------------------------



void print_fish(Fish fish){
    printf("Fish: ");
    for(int i=0; i<DIMENSIONS; i++){
        printf("%f ", fish.position[i]);
    }
    printf("weight: %f fitness: %f\n", fish.weight, fish.fitness);
}

// Funzione per inizializzare un singolo pesce
void initFish(Fish *fish) {
    for (int i = 0; i < DIMENSIONS; i++) {
        fish->position[i] = ((double)rand() / RAND_MAX) * (BOUNDS_MAX - BOUNDS_MIN) + BOUNDS_MIN;  
        fish->new_position[i] = fish->position[i];
    }

    fish->weight = W_SCALE / 2;   // Peso iniziale

    fish->fitness = objective_function(fish->position)*MULTIPLIER;        // Fitness iniziale //TODO: capire qual è il valore migliore di inizializzazione
    fish->new_fitness = fish->fitness;     // Fitness iniziale //TODO: capire qual è il valore migliore di inizializzazione

    fish->max_individual_step = 0.5; //TODO: capire qual è il valore migliore di inizializzazione
}

// Funzione per inizializzare un array di pesci
void initFishArray(Fish* fishArray) {
    for (int i = 0; i < N_FISHES; i++) {
        initFish(&fishArray[i]);  // Inizializza ciascun pesce
        print_fish(fishArray[i]);
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


        if ((delta_fitness) > *max_delta_fitness_improvement) {
            *max_delta_fitness_improvement = delta_fitness;
        }
    }
    // -------------- Finish update the collective variables

    // Update fish position considering only its individual movement
    for (int d = 0; d < DIMENSIONS; d++)
    {
        if (delta_fitness > 0) {
            printf("Update for individual movement of %f, because of delta fitness %f  ", fish->new_position[d]-fish->position[d], delta_fitness);
            printf("new_fitness: %f , old_fitness: %f\n ", fish->new_fitness, fish->fitness);
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
        printf("Update for collective movement of %f\n", fish->new_position[d]-fish->position[d]);
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
    // if (*max_delta_fitness_improvement > 0.0) { // Avoid division by zero
    //     fish->weight += (fish->new_fitness - fish->fitness)/ *max_delta_fitness_improvement;
    // }else{
    //     fish->weight += (fish->new_fitness - fish->fitness)/ 100; //tentativo di dividere per qualcosa
    // }
    fish->weight += (fish->new_fitness - fish->fitness);
    if (fish->weight<=0.0) {
        fish->weight = 0.1; //TODO: non siamo sicure di questa cosa...
    } else if (fish->weight>W_SCALE) {
        fish->weight = W_SCALE;
    }

    //che qui la delta fitness sia positiva, non ci interessa...
    //a noi interessa che la delta fitness sia positiva prima di fare il movimento singolo
    fish->fitness = fish->new_fitness; 
}

void updateWeightsArray(Fish *fishArray,  float *max_delta_fitness_improvement) {
    for (int i = 0; i < N_FISHES; i++) {
        updateWeights(&fishArray[i], max_delta_fitness_improvement);  
        print_fish(fishArray[i]);
    }
}




//-------------------------------------------------------------------------------------------
//------------------------------- MAIN ------------------------------------------------------
//-------------------------------------------------------------------------------------------

int main() {
    char filename[50];
    sprintf(filename, "../evolution_logs/%s_%dd_log.json",FUNCTION, DIMENSIONS);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // float volitive_step = 0.2; //TODO: Per ora non implementato. Forse va messo all'interno della struct del pesce come il max_individual_step
    
    float total_fitness = 0.0;
    float weighted_total_fitness[DIMENSIONS];
    float max_improvement = 0.0;
    srand(time(NULL));  // Seed for random number generation

    // INITIALIZATION
    Fish fishes[N_FISHES];
    fprintf(file, "[\n");
    initFishArray(fishes);
    write_fishes_to_json(fishes, file, 0);

    // MAIN LOOP
    for (int iter = 0; iter < MAX_ITER; iter++) {
        
        variables_reset(&total_fitness, weighted_total_fitness, &max_improvement);

        // INDIVIDUAL MOVEMENT
        printf("\n-------------------------ITER %d-------------------\n", iter);
        individualMovementArray(fishes, &total_fitness, weighted_total_fitness, &max_improvement);
        printf("total fitness: %f   ", total_fitness);
        
        for (int d = 0; d<DIMENSIONS; d++){
            printf("wtf: %f", weighted_total_fitness[d]);
        }
        printf("  max improvement: %f\n", max_improvement);

        // COLLECTIVE MOVEMENT
        collectiveMovementArray(fishes, &total_fitness, weighted_total_fitness);

        // UPDATE WEIGHTS
        updateWeightsArray(fishes, &max_improvement);

        // SAVE ON FILE
        write_fishes_to_json(fishes, file, iter==MAX_ITER-1?1:0);

        
    }

    fprintf(file, "\n]");
    fclose(file);

    return 0;
}
