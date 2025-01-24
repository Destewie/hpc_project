#!/bin/bash

# Compilazione del programma C
echo "Compilazione del programma..."
gcc sequential_FSS_simple.c -o main -lm
if [ $? -ne 0 ]; then
    echo "Errore durante la compilazione."
    exit 1
fi

# Esecuzione del programma
echo "Esecuzione del programma..."
./main
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione di main."
    exit 1
fi

# Avvio dello script di visualizzazione Python
echo "Avvio della visualizzazione..."
python3 ../../visualization/visuals_1d.py
if [ $? -ne 0 ]; then
    echo "Errore durante l'esecuzione dello script Python."
    exit 1
fi

echo "Operazione completata con successo!"

