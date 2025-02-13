#!/bin/bash

# Compilazione del programma C
echo "Compilazione del programma..."
gcc sequential_FSS_more_schools.c -o main -lm
if [ $? -ne 0 ]; then
    echo "Errore durante la compilazione."
    exit 1
fi

# Esecuzione del programma
echo "Esecuzione del programma..."
# first param   = number of schools
# second param  = number of fishes per school
# third param   = number of dimensions
# fourth param  = number of iterations
# fifth param  = update frequency 

./main 10 10 3 100 1
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione di main (10, 10, 3, 100, 1)."
    exit 1
fi

./main 20 10 3 100 1
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione di main (20, 10, 3, 100, 1)."
    exit 1
fi

./main 40 10 3 100 1
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione di main (40, 10, 3, 100, 1)."
    exit 1
fi

./main 80 10 3 100 1
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione di main (80, 10, 3, 100, 1)."
    exit 1
fi

./main 160 10 3 100 1
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione di main (160, 10, 3, 100, 1)."
    exit 1
fi

