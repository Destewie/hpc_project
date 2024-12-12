#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIM 2
#define N_FISHES 3
#define MAX_ITER 10
#define BOUNDS_MIN -10.0   // Minimum bound of the search space
#define BOUNDS_MAX 10.0    // Maximum bound of the search space
#define BOUNDS_MIN_W 0.1   // Minimum bound of the search space
#define BOUNDS_MAX_W 10.0    // Maximum bound of the search space
#define W_SCALE 10.0

typedef struct{
    double position[DIM];
    double weight;
    double fitness;
}Fish;

void print_fish(Fish fish){
    printf("Fish: ");
    for(int i=0; i<DIM; i++){
        printf("%f ", fish.position[i]);
    }
    printf("weight: %f fitness: %f\n", fish.weight, fish.fitness);
}

double rosenbrok(double *x) {
    double sum = 0.0;
    for (int i = 0; i < DIM - 1; i++) {
        double term1 = 100.0 * pow(x[i + 1] - x[i] * x[i], 2);
        double term2 = pow(1.0 - x[i], 2);
        sum += term1 + term2;
    }
    return sum;
}

double sphere(double *x) {
    double sum = 0.0;
    for (int i = 0; i < DIM; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}

double objective_function(char* function, double *x) {
    if (function == "rosenbrok") {
        return rosenbrok(x);
    } else if (function == "sphere") {
        return sphere(x);
    } else {
        return 0.0;
    }
}

void write_fish_to_json(Fish *fishes,FILE *file, int last) {

    fprintf(file, "\t[\n");

    for (int i = 0; i < N_FISHES; i++) {
        if (DIM==1){
            fprintf(file, "\t\t{\"x\": [%.6f ],", fishes[i].position[0]);
        }else if(DIM==2){
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

int main() {
    FILE *file = fopen("../evolution_logs/spherical_1d_log.json", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    char* function = "sphere";
    int multiplicator = -1; // 1 in case of maximization, -1 in case of minimization
    float individual_step = 0.5;
    float volitive_step = 0.2;
    srand(time(NULL));  // Seed for random number generation

    // // INITIALIZATION
    Fish fishes[N_FISHES];// Last element is the weight of the fish
    fprintf(file, "[\n");

    // // Initialize the fishes with random position and size
    for (int i = 0; i < N_FISHES; i++) {
        for (int j = 0; j < DIM; j++) {
            fishes[i].position[j] = (double)(BOUNDS_MIN + (BOUNDS_MAX - BOUNDS_MIN) * rand() / (RAND_MAX + 1.0));
        }
        fishes[i].weight = W_SCALE/2;
        fishes[i].fitness = objective_function(function, fishes[i].position);
        print_fish(fishes[i]);
    }
    write_fish_to_json(fishes, file, 0);

    // MAIN LOOP
    for (int iter = 0; iter < MAX_ITER; iter++) {
        // printf("Iteration %d -> FISH[0] %f %f %f \n ", iter, fishes[0].position[0], fishes[0].weight, fishes[0].fitness);

        // Individial movement
        float total_fitness = 0.0;
        float weighted_total_fitness = 0.0;
        float max_imp = 0.0;
        for(int i=0; i<N_FISHES; i++){
            float angle = ((float)rand()/RAND_MAX) *2* 3.14; // with 2 dimensions there are just 2 possibilities
            double new_pos[DIM];
            new_pos[0]= fishes[i].position[0] + individual_step*cos(angle);

            float new_fitness = objective_function(function, new_pos);
            //ci interessano solo i movimenti "buoni" -> che si avvicinano al nostro goal
            if ((new_fitness-fishes[i].fitness)*multiplicator>0){
                fishes[i].position[0] = new_pos[0];
                weighted_total_fitness += fabs(individual_step*cos(angle))* (new_fitness - fishes[i].fitness)*multiplicator;
                total_fitness += (new_fitness - fishes[i].fitness)*multiplicator;
            }
            if (max_imp< fabs((new_fitness - fishes[i].fitness)*multiplicator)){
                max_imp = fabs((new_fitness - fishes[i].fitness)*multiplicator);
            }

        }
        printf("Total fitness %f %f\n", total_fitness, weighted_total_fitness);

        // Collective movement
        for(int i=0; i<N_FISHES; i++){
                fishes[i].position[0] += weighted_total_fitness/total_fitness;
        }

        printf("Max imp %f\n", max_imp);
        // Update fish
        for(int i=0; i<N_FISHES; i++){
            float new_fitness = objective_function(function, fishes[i].position);
            fishes[i].weight += (new_fitness - fishes[i].fitness)*multiplicator/max_imp; // qua potrebbe dare nan nel caso in cui max_imp->0
            if (fishes[i].weight<1.0){
                printf("INFO weight<1");
                fishes[i].weight = 1.0;
            }else if(fishes[i].weight>W_SCALE){
                fishes[i].weight = W_SCALE;
            }
            fishes[i].fitness = new_fitness;
        }
        write_fish_to_json(fishes, file, iter==MAX_ITER-1?1:0);
    }

    printf("Final positions of the fishes\n");
    int index_best = 0;
    for (int i = 0; i < N_FISHES; i++) {
        print_fish(fishes[i]);
        if(fishes[i].fitness<fishes[index_best].fitness){
            index_best = i;
        }
    }

    printf("Best fish-> ");
    print_fish(fishes[index_best]);
    fprintf(file, "\n]");
    fclose(file);

    return 0;
}
