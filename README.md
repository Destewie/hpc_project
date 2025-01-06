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
- [ ] Capiamo perché dobbiamo multiplicator con la rosenbrock va messo ad 1 e non a -1 nonostante noi vogliamo minimizzarla (quando invece guardiamo la min_sphere mettiamo -1 per minimizzare (crediamo))

## TODO

- [ ] Controllare che la scrittura del json funzioni per n dimensioni
- [ ] Capire ed implementare un criterio di fine 
- [ ] Una volta che abbiamo una versione funzionante (con le cose base), provare ad aggiungere anche altre funzionalità

## NOTE IMPORTANTI SULLE NOSTRE SCELTE IMPLEMENTATIVE 
- Abbiamo deciso di aggiornare i pesi dei pesci solo dopo il loro individual movement. In ogni caso, questo non influisce negativamente perché ad ogni individual movement verrà incluso anche il collective movement del ciclo prima.

- Il peso di un pesce non è direttamente correlabile alla sua posizione nella funzione (non è vero che un pesce più vicino al minimo abbia sempre un peso più grande di quello di un pesce distante), ma tuttavia è direttamente correlabile al suo miglioramento all'interno della funzione. 
Se un pesce migliora costantemente la sua fitness, il suo peso aumenta. Tuttavia:
Se la sua fitness si stabilizza o peggiora, il peso rimane invariato o diminuisce.
Pesci con miglioramenti recenti ma con fitness assoluta più bassa possono avere pesi superiori.