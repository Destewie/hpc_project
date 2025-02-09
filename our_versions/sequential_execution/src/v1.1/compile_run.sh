#!/bin/bash

# Compilazione del programma C
echo "Compilazione del programma..."
gcc sequential_FSS_more_schools.c -o main -lm
if [ $? -ne 0 ]; then
    echo "Errore durante la compilazione."
    exit 1
fi

# Esecuzione del programma
echo "Esecuzione del programma...\n"
./main 100 100 100 1
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione di main."
    exit 1
fi