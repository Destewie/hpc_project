# Cancella i file solo se esistono
for file in parallel_FSO.sh.e* parallel_FSO.sh.o* parallel_FSO; do
    if [ -e "$file" ]; then
        rm "$file"
    fi
done



# Sottometti il job
qsub parallel_FSO.sh