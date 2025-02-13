# Cancella i file solo se esistono
for file in sequential_FSO.sh.e* sequential_FSO.sh.o* sequential_FSO; do
    if [ -e "$file" ]; then
        rm "$file"
    fi
done

# Sottometti il job
qsub sequential_FSO.sh