# hpc_project
Io e Annina che facciamo muovere i pesci eheh



### FONTI

- https://en.wikipedia.org/wiki/Fish_School_Search
- https://it.wikipedia.org/wiki/Swarm_intelligence -> FSS è un algoritmo di Swarm Intelligence
- https://www.mql5.com/en/articles/11841 -> c'è dello pseudocodice in C (non ho letto bene tutto)
- https://www.youtube.com/watch?v=kSoA_zzrqsc -> Video spiegazione


## CRITICITÀ
- [ ] Gestione dell'aggiornamento del peso -> divisione per il max_improvement: abbiamo visto che si comporta meglio senza quella divisione. Spesso infatti sono valori piccoli e quando quel valore è 0 non dividiamo per nulla.
- [ ] Collective movement -> siamo sicur3 che stia funzionando correttamente? controlliamo i paper