# Cancella i file solo se esistono
for file in multithreading_FSO.sh.e* multithreading_FSO.sh.o* multithreading_FSO; do
    if [ -e "$file" ]; then
        rm "$file"
    fi
done

# Sottometti il job
qsub multithreading_FSO.sh