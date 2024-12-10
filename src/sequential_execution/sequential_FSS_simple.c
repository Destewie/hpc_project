#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DIM 2              
#define N_FISHES 30       
#define MAX_ITER 400     
#define BOUNDS_MIN -10.0   // Minimum bound of the search space
#define BOUNDS_MAX 10.0    // Maximum bound of the search space
#define BOUNDS_MIN_W 0.1   // Minimum bound of the search space
#define BOUNDS_MAX_W 10.0    // Maximum bound of the search space

typedef struct{
    double position[DIM];
    double weight;
    double fitness;
}Fish;


double rosenbrok(double *x) {
    double sum = 0.0;
    for (int i = 0; i < DIM - 1; i++) {
        double term1 = 100.0 * pow(x[i + 1] - x[i] * x[i], 2);
        double term2 = pow(1.0 - x[i], 2);
        sum += term1 + term2;
    }
    return sum;
}


int main() {
    float individual_step = 0.5;
    float volitive_step = 0.2;
    srand(time(NULL));  // Seed for random number generation

    // // INITIALIZATION
    Fish fishes[N_FISHES];// Last element is the weight of the fish
    
    // // Initialize the fishes with random position and size
    for (int i = 0; i < N_FISHES; i++) {
        for (int j = 0; j < DIM; j++) {
            fishes[i].position[j] = (double)(BOUNDS_MIN + (BOUNDS_MAX - BOUNDS_MIN) * rand() / (RAND_MAX + 1.0));
        }
        fishes[i].weight = 0.1;
        fishes[i].fitness = rosenbrok(fishes[i].position);
        printf("Fish %d: %f %f %f %f\n", i, fishes[i].position[0],fishes[i].position[1], fishes[i].weight, fishes[i].fitness);
    }

    // MAIN LOOP
    for (int iter = 0; iter < MAX_ITER; iter++) {
        printf("Iteration %d -> FISH[0] %f %f %f %f \n ", iter, fishes[0].position[0],fishes[0].position[1], fishes[0].weight, fishes[0].fitness);

        // Individial movement
        float total_fitness= 0.0;
        float weighted_total_fitness = 0.0;
        for(int i=0; i<N_FISHES; i++){
            float angle = ((float)rand()/RAND_MAX) *2* 3.14; // with 2 dimensions there are just 2 possibilities
            float new_pos[DIM];
            new_pos[0]= fishes[i].position[0] + individual_step*cos(angle);
            new_pos[1]= fishes[i].position[1] + individual_step*sin(angle); 
            float new_fitness = rosenbrok(new_pos);      
            total_fitness += fishes[i].fitness - new_fitness;  
            if (new_fitness<fishes[i].fitness){
                fishes[i].position[0] = new_pos[0];
                fishes[i].position[1] = new_pos[1];
                fishes[i].fitness = rosenbrok(fishes[i].position);
                weighted_total_fitness += individual_step* ( fishes[i].fitness - new_fitness);
            }
        }
        printf("Total fitness %f %f\n", total_fitness, weighted_total_fitness);  

        // Collective movement
        for(int i=0; i<N_FISHES; i++){
                fishes[i].position[0] += weighted_total_fitness/total_fitness;
                fishes[i].position[1] += weighted_total_fitness/total_fitness;
        }


        // Update the fish
        float max_imp = 0.0;
        for(int i=0; i<N_FISHES; i++){
            if(max_imp<(fishes[i].fitness-rosenbrok(fishes[i].position))){
                max_imp = fishes[i].fitness-rosenbrok(fishes[i].position);
            }
        }
        printf("Max imp %f", max_imp);

        for(int i=0; i<N_FISHES; i++){
            float fitness1 = rosenbrok(fishes[i].position);
            fishes[i].weight += (fishes[i].fitness-fitness1)/max_imp;
            fishes[i].fitness = fitness1;
        }
    }

    printf("Final positions of the fishes\n");
    int index_best = 0;
    for (int i = 0; i < N_FISHES; i++) {
        printf("Fish %d: %f %f %f %f\n", i, fishes[i].position[0],fishes[i].position[1], fishes[i].weight, fishes[i].fitness);
        if(fishes[i].fitness<fishes[index_best].fitness){
            index_best = i;
        }
    }

    
    printf("Best fish: {%d} %f %f %f %f\n", index_best, fishes[index_best].position[0],fishes[index_best].position[1], fishes[index_best].weight, fishes[index_best].fitness);

   
    return 0;
}