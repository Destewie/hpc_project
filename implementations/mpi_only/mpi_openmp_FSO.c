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
#include <unistd.h>
#include <sys/time.h>
#include <string.h> // Include for strcmp
#include <mpi.h>

#define FUNCTION "min_ackley"   //TODO: Capire se, al posto di fare un controllo su una stringa, possiamo passare alle funzioni direttamente un puntatore ad una funzione (in modo comodo, se no lasciamo perdere)
#define MULTIPLIER -1   // 1 in case of maximization, -1 in case of minimization

#define SPAWN_BOUNDS_MIN -70.0   // Minimum bound of the spawn space
#define SPAWN_BOUNDS_MAX 70.0    // Maximum bound of the spawn space
#define MAX_INDIVIDUAL_STEP 2 // Maximum step for individual movement
#define MAX_VOLITIVE_STEP 0.2 // Maximum step for collective movement

#define W_SCALE_MIN 1.0         // Minimum Bounds of the weight
#define W_SCALE_MAX 10.0        // Maximum bount of the weight
#define BREEDING_THRESHOLD 7.0 // minimus threshold of weight to breedh new fishes
#define A 10.0 //rastrigin param

#define LOG 0 // 1 to log the results, 0 otherwise

//10 very different colors that will be used by a python script to plot the results
const char *COLORS[] = {"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"};

typedef struct {
    double *position;
    double *new_position;

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

void WriteFishesToJson(Fish *fishes, FILE *file, int first_call, int last_call, const int N_FISHES_PER_PROCESS, const int N_SCHOOLS, const int DIMENSIONS) {

    if (first_call) {
        // Scrive l'apertura dell'array principale solo se è la prima chiamata
        fprintf(file, "[\n");
    }

    fprintf(file, "\t[\n");

    for (int i = 0; i < N_FISHES_PER_PROCESS*N_SCHOOLS; i++) {
        int school_index = i / N_FISHES_PER_PROCESS;
        fprintf(file, "\t\t{\"x\": [");
        for (int d = 0; d < DIMENSIONS; d++) {
            fprintf(file, "%.6f", fishes[i].position[d]);
            if (d < DIMENSIONS - 1) {
                fprintf(file, ", ");
            }
        }
        fprintf(file, "], \"weight\": %.6f, \"color\": \"%s\"}", fishes[i].weight, COLORS[school_index]);

        if (i < N_FISHES_PER_PROCESS*N_SCHOOLS - 1) {
            fprintf(file, ",\n");
        } else {
            fprintf(file, "\n");
        }
    }

    if (last_call) {
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
void variablesReset(float* tot_fitness, float* weighted_tot_fitness,  float* max_improvement, const int DIMENSIONS) {

    for (int d = 0; d < DIMENSIONS; d++) {
        weighted_tot_fitness[d] = 0.0;
    }

    *tot_fitness = 0.0;
    *max_improvement = 0.0;
}


//-------------------------------------------------------------------------------------------
//---------------------------- MATH FUNCTIONS -----------------------------------------------
//-------------------------------------------------------------------------------------------

double rosenbrok(double *x,int DIMENSIONS) {
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS - 1; i++) {
        double term1 = 100.0 * pow(x[i + 1] - x[i] * x[i], 2);
        double term2 = pow(1.0 - x[i], 2);
        sum += term1 + term2;
    }
    return sum;
}

double rastrigin(double *x,const int DIMENSIONS) {
    double sum = A * DIMENSIONS; // Start with A * DIM
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += x[i] * x[i] - A * cos(2 * M_PI * x[i]);
    }
    return sum;
}

double min_sphere(double *x,const int DIMENSIONS) {
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

double max_sphere(double *x,const int DIMENSIONS) {
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += -(x[i] * x[i]);
    }
    return sum;
}

double ackley(double *x,const int DIMENSIONS){
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum1 += x[i] * x[i];
        sum2 += cos(2 * M_PI * x[i]);
    }
    return -20 * exp(-0.2 * sqrt(sum1 / DIMENSIONS)) - exp(sum2 / DIMENSIONS) + 20 + M_E;
}

double min_schwefel(double *x,const int DIMENSIONS) {
    double sum = 0.0;
    for (int i = 0; i < DIMENSIONS; i++) {
        sum += x[i] * sin(sqrt(fabs(x[i])));
    }
    return 418.9829 * DIMENSIONS - sum;
}


double objectiveFunction(double *x,const int DIMENSIONS) {
    if (strcmp(FUNCTION, "min_rosenbrock") == 0) {
        return rosenbrok(x, DIMENSIONS);
    } else if (strcmp(FUNCTION, "min_sphere") == 0) {
        return min_sphere(x, DIMENSIONS);
    } else if (strcmp(FUNCTION, "max_sphere") == 0) {
        return max_sphere(x, DIMENSIONS);
    } else if (strcmp(FUNCTION, "min_rastrigin") == 0) {
        return rastrigin(x, DIMENSIONS);
    } else if (strcmp(FUNCTION, "min_ackley") == 0) {
        return ackley(x, DIMENSIONS);
    } else if (strcmp(FUNCTION, "min_schwefel") == 0) {
        return min_schwefel(x, DIMENSIONS);
    }else{
        return 0.0;
    }
}

//-------------------------------------------------------------------------------------------
//---------------------------- FISH ---------------------------------------------------------
//-------------------------------------------------------------------------------------------



void printFish(Fish fish,const int DIMENSIONS){
    printf("Fish: ");
    for(int i=0; i<DIMENSIONS; i++){
        printf("pos: %f - ", fish.position[i]);
    }
    printf("\tweight: %f \tfitness: %f\n", fish.weight, fish.fitness);
}

// Funzione per inizializzare un singolo pesce
void initFish(Fish *fish, int process_rank, const int DIMENSIONS, const int N_FISHES_PER_PROCESS, const int N_SCHOOLS ) { // the number of schools is needed to divide the space

    // Posizioni iniziali divise per banco
    double portion_bounds = (SPAWN_BOUNDS_MAX - SPAWN_BOUNDS_MIN) / N_SCHOOLS; // Calcola la porzione corretta per ciascun banco
    double lower_bound = SPAWN_BOUNDS_MIN + process_rank * portion_bounds; // Limite inferiore per il banco
    double upper_bound = lower_bound + portion_bounds; // Limite superiore per il banco
    fish->position = malloc(DIMENSIONS*sizeof(double));
    fish->new_position = malloc(DIMENSIONS*sizeof(double));

    for (int d = 0; d < DIMENSIONS; d++) {
        // // Posizioni iniziali random
        // fish->position[d] = ((double)rand() / RAND_MAX) * (BOUNDS_MAX - SPAWN_BOUNDS_MIN) + SPAWN_BOUNDS_MIN;

        // Posizioni iniziali divise per banco
        if (d == 0) {
            fish->position[d] = ((double)rand() / RAND_MAX) * (upper_bound - lower_bound) + lower_bound;
            // printf("[D0] lower_bound: %f, upper_bound: %f\n, x: %f", lower_bound, upper_bound, fish->position[d]);
        } else {
            fish->position[d] = ((double)rand() / RAND_MAX) * (SPAWN_BOUNDS_MAX - SPAWN_BOUNDS_MIN) + SPAWN_BOUNDS_MIN;
        }

        fish->new_position[d] = fish->position[d];
    }

    fish->weight = W_SCALE_MAX / 2;   // Peso iniziale
    fish->previous_cycle_weight = fish->weight;

    fish->fitness = objectiveFunction(fish->position, DIMENSIONS)*MULTIPLIER;        // Fitness iniziale //TODO: capire qual è il valore migliore di inizializzazione
    fish->new_fitness = fish->fitness;

    fish->max_individual_step = MAX_INDIVIDUAL_STEP; //TODO: capire qual è il valore migliore di inizializzazione e come aggiornarlo dinamicamente
    fish->max_volitive_step = MAX_VOLITIVE_STEP; //TODO: capire qual è il valore migliore di inizializzazione e come aggiornarlo dinamicamente
}

// Funzione per inizializzare un array di pesci
void initFishArray(Fish* fishArray, const int DIMENSIONS, const int N_FISHES_PER_PROCESS, const int N_SCHOOLS, int rank) {
    for (int i = 0; i < N_FISHES_PER_PROCESS; i++){
        initFish(&fishArray[i], rank, DIMENSIONS, N_FISHES_PER_PROCESS, N_SCHOOLS);  // Inizializza ciascun pesce
        // printFish(fishArray[i]);
    }
}

void individualMovement(Fish *fish,
                        float *delta_fitness_out,
                        float *weighted_move_out,
                        float *max_improve_out,
                        int DIMENSIONS) {

    // Local new position buffer
    double *new_pos = (double *)malloc(DIMENSIONS * sizeof(double));
    if (!new_pos) return; // allocation failure

    // Compute candidate position and fitness
    for (int d = 0; d < DIMENSIONS; ++d) {
        double norm_move = ((double)rand() / (double)RAND_MAX) * 2.0 - 1.0;
        new_pos[d] = fish->position[d] + norm_move * fish->max_individual_step;
    }
    double new_fit = objectiveFunction(new_pos, DIMENSIONS) * MULTIPLIER;
    double delta = new_fit - fish->fitness;

    // Initialize outputs
    *delta_fitness_out = 0.0f;
    *max_improve_out   = 0.0f;
    for (int d = 0; d < DIMENSIONS; ++d) {
        weighted_move_out[d] = 0.0f;
    }

    fish->new_fitness = (float)new_fit; // the new fitness is updated in every case

    // If fitness improves, update outputs and fish state
    if (delta > 0.0) {
        *delta_fitness_out = (float)delta;
        *max_improve_out   = (float)delta;
        for (int d = 0; d < DIMENSIONS; ++d) {
            weighted_move_out[d] = (float)((new_pos[d] - fish->position[d]) * delta);
            fish->position[d] = new_pos[d];
        }
    }

    free(new_pos);
}

void individualMovementArray(Fish* fishArray,
                             float* tot_delta_fitness,
                             float* weighted_tot_delta_fitness,
                             float* max_delta_fitness_improvement,
                             int current_iter,
                             int N_FISHES_PER_PROCESS,
                             int DIMENSIONS,
                             int UPDATE_FREQUENCY) {


        // Work distribution: each fish
        for (int idx = 0; idx < N_FISHES_PER_PROCESS; ++idx) {
            // int s = idx / N_FISHES_PER_PROCESS;
            Fish *fish = &fishArray[idx];
            float dfit;
            float *wmove = (float *)malloc(DIMENSIONS * sizeof(float));
            float improve;

            individualMovement(fish, &dfit, wmove, &improve, DIMENSIONS);

            // Accumulate into thread-local buffers
            *tot_delta_fitness += dfit;
            *max_delta_fitness_improvement = fmaxf(*max_delta_fitness_improvement, improve);
            for (int d = 0; d < DIMENSIONS; ++d) {
                weighted_tot_delta_fitness[d] += wmove[d];
            }

            free(wmove);
        }


    // Periodic global update broadcasting complete totals
    if (current_iter % UPDATE_FREQUENCY == 0) {

        // Append tot_delta_fitness to weighted_tot_delta_fitness creating a new array
        float *temp_array = (float *)malloc((DIMENSIONS + 1) * sizeof(float));
        temp_array[0] = *tot_delta_fitness;
        for (int d = 0; d < DIMENSIONS; d++) {
            temp_array[d + 1] = weighted_tot_delta_fitness[d];
        }

        //TODO: SCRIVERE A COSA SERVONO LE REQUEST E GLI STATUS IN MPI
        MPI_Request request0;
        MPI_Status status0;
        MPI_Request request1;
        MPI_Status status1;


        MPI_Iallreduce(MPI_IN_PLACE, temp_array, DIMENSIONS + 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, &request1);
        MPI_Iallreduce(MPI_IN_PLACE, &max_delta_fitness_improvement, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD,  &request0);

        // MPI_Allreduce(MPI_IN_PLACE, temp_array, DIMENSIONS + 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        // MPI_Allreduce(MPI_IN_PLACE, &max_delta_fitness_improvement, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Update the original arrays
        MPI_Wait(&request1, &status1);
        *tot_delta_fitness = temp_array[0];
        for (int d = 0; d < DIMENSIONS; d++) {
            weighted_tot_delta_fitness[d] = temp_array[d + 1];
        }

        MPI_Wait(&request0, &status0);
    }

}


void updateWeights(Fish *fish, float max_delta_fitness_improvement) {
    if (max_delta_fitness_improvement != 0.0) { // Avoid division by zero
        fish->weight += (fish->new_fitness - fish->fitness)/ max_delta_fitness_improvement;
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

void updateWeightsArray(Fish *fishArray,  float max_delta_fitness_improvement, const int N_FISHES_PER_PROCESS) {
    for (int i = 0; i < N_FISHES_PER_PROCESS; i++) {
        updateWeights(&fishArray[i], max_delta_fitness_improvement);
    }
}


// Assumes Fish defined elsewhere, and objectiveFunction, MULTIPLIER available

// Local-version of collective movement: uses thread-local tot and weighted arrays
static inline void collectiveMovementLocal(Fish *fish,
                                           float *tot_delta_fitness,
                                           const float *weighted_delta,
                                           int DIMENSIONS) {
    // Avoid division by zero
    if (*tot_delta_fitness == 0.0f) {
        *tot_delta_fitness = 1.0f;
    }
    // Update each dimension
    for (int d = 0; d < DIMENSIONS; ++d) {
        float move = weighted_delta[d] / *tot_delta_fitness;
        fish->new_position[d] = fish->position[d] + move;
        fish->position[d] = fish->new_position[d];
    }
    // Recompute fitness after collective move
    fish->new_fitness = objectiveFunction(fish->position, DIMENSIONS) * MULTIPLIER;
}

void collectiveMovementArray(Fish *fishArray,
                             float *tot_delta_fitness,
                             float *weighted_tot_delta_fitness,
                             int N_FISHES_PER_PROCESS,
                             int DIMENSIONS) {
    for (int i = 0; i < N_FISHES_PER_PROCESS; ++i) {
        Fish *fish = &fishArray[i];
        collectiveMovementLocal(fish,
                                tot_delta_fitness,
                                weighted_tot_delta_fitness,
                                DIMENSIONS);
    }
}


void calculateBarycentersAndSumWeights(Fish *fishArray, float* barycenter, float *old_sum, float *new_sum, int current_iter, const int UPDATE_FREQUENCY, const int DIMENSIONS, const int N_FISHES_PER_PROCESS) {

    //common part for all processes where they calcuate their own data
    // BARYCENTER
    float numerator[DIMENSIONS];
    // float denominator[DIMENSIONS]; it is the same for everyone

    for (int d = 0; d < DIMENSIONS; d++){
        numerator[d] = 0.0;
        // denominator[d] = 0.0;
    }
    *new_sum = 0.0;
    *old_sum = 0.0;

    for (int d = 0; d < DIMENSIONS; d++) {
        for (int i = 0; i < N_FISHES_PER_PROCESS; i++) {
            numerator[d] += fishArray[i].position[d] * fishArray[i].weight;
            *new_sum += fishArray[i].weight;
            *old_sum += fishArray[i].previous_cycle_weight;
        }

        if (*new_sum != 0.0) {
            barycenter[d] = numerator[d] / *new_sum;
        } else {
            // printf("Denominator is zero...\n");
        }
    }

    // exchange of information among processes
    if (current_iter % UPDATE_FREQUENCY == 0) {
        // creare un unico array per il banco con tutte le cose dentro
        // tutte le cose per UN SOLO BANCO sono: il suo peso totale all'iterazione corrente, il suo peso totale all'iterazione precedente e la moltiplicazione del suo baricentro(d valori) per il suo peso attuale
        float* temp_array = (float *)malloc((DIMENSIONS + 2) * sizeof(float));
        temp_array[0] = *old_sum;
        temp_array[1] = *new_sum;

        for (int d = 0; d < DIMENSIONS; d++) {
            temp_array[d + 2] = numerator[d];
        }

        MPI_Allreduce(MPI_IN_PLACE, temp_array, DIMENSIONS + 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        // Update the original arrays
        *old_sum = temp_array[0];
        *new_sum = temp_array[1];

        if (*new_sum != 0.0) {
            for (int d = 0; d < DIMENSIONS; d++) {
                barycenter[d] = temp_array[d + 2] / *new_sum;
            }
        }else {
            // printf("Denominator is zero...\n");
        }

    }
}


// Assumes Fish, calculateBarycenters, calculateSumWeights defined elsewhere

void volitivePositionUpdateArray(Fish *fishArray,
                                 int shrink,
                                 float *barycenter,
                                 const int N_FISHES_PER_PROCESS,
                                 const int DIMENSIONS) {

    for (int idx = 0; idx < N_FISHES_PER_PROCESS; ++idx) {
        Fish *fish = &fishArray[idx];
        for (int d = 0; d < DIMENSIONS; ++d) {
            double rand_mult = fmin(((double)rand() / RAND_MAX) + 0.1, 1.0);

            double delta = fish->position[d] - barycenter[d];

            if (shrink == 1) { //attracted
                fish->position[d] -= fish->max_volitive_step * rand_mult * delta;
            } else { // repelled
                fish->position[d] += fish->max_volitive_step * rand_mult * delta;
            }
            if (fish->position[d] > 1000.0 || fish->position[d] < -1000.0) {
                shrink = 1;
                // printf("Error: Fish position out of bounds. Problematic error %f\n", fish->position[d]);

                fish->position[d] = fmax(fmin(fish->position[d], 1000.0), -1000.0);
            }
        }
    }
}

void collectiveVolitiveArray(Fish *fishes,
                             int current_iter,
                             const int DIMENSIONS,
                             const int N_FISHES_PER_PROCESS,
                             const int UPDATE_FREQUENCY) {
    float *barycenter = malloc(DIMENSIONS * sizeof(float));
    float old_weights;
    float new_weights;

    calculateBarycentersAndSumWeights(fishes, barycenter, &old_weights, &new_weights, current_iter, UPDATE_FREQUENCY, DIMENSIONS, N_FISHES_PER_PROCESS);



    if (old_weights == new_weights) return;
    int shrink = (old_weights < new_weights) ? 1 : 0;


    // Update positions in parallel per fish
    volitivePositionUpdateArray(fishes,
                                shrink,
                                barycenter,
                                N_FISHES_PER_PROCESS,
                                DIMENSIONS);

    // Update previous_cycle_weight for each fish
    // int base = s * N_FISHES_PER_PROCESS;
    for (int i = 0; i < N_FISHES_PER_PROCESS; ++i) {
        fishes[i].previous_cycle_weight = fishes[i].weight;
    }

    free(barycenter);
}



void breeding(Fish *fishes, int current_iter, const int N_FISHES_PER_PROCESS, const int DIMENSIONS) {
    // Indici del miglior, secondo miglior e peggiore pesce del banco
    int first_index = 0;
    int second_index = -1; // Inizializza con un valore non valido
    int worst_index = 0;

    for (int i = 1; i < N_FISHES_PER_PROCESS; i++) {
        int current_index = i;

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
        fishes[worst_index].fitness = objectiveFunction(fishes[worst_index].position, DIMENSIONS) * MULTIPLIER;
    }
}


//-------------------------------------------------------------------------------------------
//------------------------------- MAIN ------------------------------------------------------
//-------------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {

    if (argc != 6) {
        printf("Usage: %s N_FISHES_PER_PROCESS DIMENSIONS MAX_ITER UPDATE_FREQUENCY\n", argv[0]);
        return 1;
    }

    // MPI initialization
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //prendiamo in input n_fishes_per_process
    const int N_FISHES_PER_PROCESS = atoi(argv[1]);
    const int DIMENSIONS = atoi(argv[2]);
    const int MAX_ITER = atoi(argv[3]);
    const int UPDATE_FREQUENCY = atoi(argv[4]);
    char *place = argv[5];


    if (rank == 0) {
        printf("\nRUNNING WITH: N-PROCESSES %d - N_FISHES_PER_PROCESS %d - DIMENSIONS %d - MAX_ITER %d - UPDATE_FREQUENCY %d - PLACE %s\n", size, N_FISHES_PER_PROCESS, DIMENSIONS, MAX_ITER, UPDATE_FREQUENCY, place);
    }

    // Simulo la chiamata a whoami per ottenere il nome dell'utente corrente
    char *user;
    if (rank==0) {
        user = getenv("USER");
        if (user == NULL) {
            //return errror
            fprintf(stderr, "Error: USER environment variable not set.\n");
            MPI_Finalize();
            return 1;
        }
    }


    //create a timer
    double start = MPI_Wtime();
    double end = 0.0;

    FILE *file;
    if (rank==0 && LOG) {
        char filename[50];
        sprintf(filename, "/home/%s/hpc_project/our_versions/evolution_logs/%s_%dd_log.json", user, FUNCTION, DIMENSIONS);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return 1;
        }
    }

    // float best_fitness[N_SCHOOLS];
    float total_fitness;
    float* weighted_total_fitness = malloc(DIMENSIONS*sizeof(float));
    if (weighted_total_fitness==NULL){
        exit(2);
    }

    float max_improvement;


    // Seed for random number generation (different for different processes)
    const int seed = (unsigned int)time(NULL) + rank;
    srand(seed);


    // INITIALIZATION
    Fish *fishes = malloc(N_FISHES_PER_PROCESS*sizeof(Fish));//creiamo un vettore per ogni processo diverso

    initFishArray(fishes, DIMENSIONS, N_FISHES_PER_PROCESS, size, rank);
    if (DIMENSIONS <= 2 && LOG) {
        WriteFishesToJson(fishes, file, 1, 0, N_FISHES_PER_PROCESS, size, DIMENSIONS);
    }

    // MAIN LOOP
    double a, b, c, d, e, f, g, h, i, l, m, n;
    // le iterazioni devono essere sequenziali quindi non le possiamo parallelizzare
    for (int iter = 1; iter < MAX_ITER; iter++) {

        a = MPI_Wtime();
        variablesReset(&total_fitness, weighted_total_fitness, &max_improvement, DIMENSIONS);
        b = MPI_Wtime();

        // INDIVIDUAL MOVEMENT
        c = MPI_Wtime();
        individualMovementArray(fishes, &total_fitness, weighted_total_fitness, &max_improvement, iter, N_FISHES_PER_PROCESS, DIMENSIONS, UPDATE_FREQUENCY);
        d = MPI_Wtime();

        // UPDATE WEIGHTS
        e = MPI_Wtime();
        updateWeightsArray(fishes, max_improvement, N_FISHES_PER_PROCESS);
        f = MPI_Wtime();

        // COLLECTIVE MOVEMENT
        g = MPI_Wtime();
        collectiveMovementArray(fishes, &total_fitness, weighted_total_fitness, N_FISHES_PER_PROCESS, DIMENSIONS);
        h = MPI_Wtime();

        // COLLECTIVE VOLITIVE MOVEMENT
        i = MPI_Wtime();
        collectiveVolitiveArray(fishes, iter, DIMENSIONS, N_FISHES_PER_PROCESS, UPDATE_FREQUENCY);
        l = MPI_Wtime();

        // BREEDING
        m = MPI_Wtime();
        breeding(fishes, iter, N_FISHES_PER_PROCESS, DIMENSIONS);
        n = MPI_Wtime();


        // SAVE ON FILE
        if (DIMENSIONS <= 2 && LOG) {
            WriteFishesToJson(fishes, file, 0, iter==MAX_ITER-1?1:0,  N_FISHES_PER_PROCESS, size, DIMENSIONS);
        }
    }

    //timer stop
    end = MPI_Wtime();

    if (rank==0 && LOG) {
        fclose(file);
    }

    MPI_Finalize();

    printf("Variablee reset- %f\n", b-a);
    printf("Individual movement- %f\n", d-c);
    printf("Update weights- %f\n", f-e);
    printf("Collective movement- %f\n", h-g);
    printf("Collective volitive- %f\n", l-i);
    printf("Breeding- %f\n", n-m);

    printf("END: %f\n", end-start);


    return 0;
}
