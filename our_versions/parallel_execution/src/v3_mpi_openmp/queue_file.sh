# Cancella i file solo se esistono
for file in mpi_openmp_FSO.sh.e* mpi_openmp_FSO.sh.o* mpi_openmp_FSO; do
    if [ -e "$file" ]; then
        rm "$file"
    fi
done

# Sottometti il job
qsub mpi_openmp_FSO.sh