# hpc_project
Io e Annina che facciamo muovere i pesci eheh

### FONTI

- https://en.wikipedia.org/wiki/Fish_School_Search
- https://it.wikipedia.org/wiki/Swarm_intelligence -> FSS è un algoritmo di Swarm Intelligence
- https://www.mql5.com/en/articles/11841 -> c'è dello pseudocodice in C (non ho letto bene tutto)
- https://www.youtube.com/watch?v=kSoA_zzrqsc -> Video spiegazione

## CRITICITÀ

- [X] Gestione dell'aggiornamento del peso -> divisione per il max_improvement: abbiamo visto che si comporta meglio senza quella divisione. Spesso infatti sono valori piccoli e quando quel valore è 0 non dividiamo per nulla.
- [X] Collective movement -> siamo sicur3 che stia funzionando correttamente? controlliamo i paper

## TODO

- [X] Controllare che la scrittura del json funzioni per n dimensioni
- [ ] Una volta che abbiamo una versione funzionante (con le cose base), provare ad aggiungere anche altre funzionalità
- [ ] cambiamo la v1 parallela per usare mpi_reduce e mpi_broadcast
- [ ] cambiamo la v1 sequenziale per far comunicare i banchi tra di loro ogni n iterazioni

## NOTE IMPORTANTI SULLE NOSTRE SCELTE IMPLEMENTATIVE 
- Abbiamo deciso di aggiornare i pesi dei pesci solo dopo il loro individual movement. In ogni caso, questo non influisce negativamente perché ad ogni individual movement verrà incluso anche il collective movement del ciclo prima.
- Il peso di un pesce non è direttamente correlabile alla sua posizione nella funzione (non è vero che un pesce più vicino al minimo abbia sempre un peso più grande di quello di un pesce distante), ma tuttavia è direttamente correlabile al suo miglioramento all'interno della funzione. 
  Se un pesce migliora costantemente la sua fitness, il suo peso aumenta. Tuttavia:
  Se la sua fitness si stabilizza o peggiora, il peso rimane invariato o diminuisce.
  Pesci con miglioramenti recenti ma con fitness assoluta più bassa possono avere pesi superiori.
- Abbiamo trovato in varie altre implementazioni che il baricentro considera la media delle posizioni pesate sulla weight dei pesci, quindi va diviso per il peso e non la posizione (come invece scritto nel paper)



## PARALLELIZZAZIONE

Diverse possibilità (by chatgpt):

- parallelizzare il calcolo della fitness (cosa più pesante da fare), secondo me un po troppa comunicazione

- proposta di fare update non-sincronizzati, non aspettare tutti (molto time consuming)

- parallelizzare l'individual movement (il collective poi va fatto insieme)

- dividere lo spazio in vari domini, parallelizzare il calcolo per dominio (dubbio: come ci si assicura che siano sempre tutti separati?)

  Ci sono vari modi per fare questo, considerando *subset di pesci*- *population based,* quindi ogni processore mantiene un set di pesci, senza considerare dove si trovano. *Geometric domain decomposition* nel caso in cui si consideri una divisione in celle, in base allo spazio. In questo caso i pesci potrebbero cambiare zona e in quel caso la loro gestione passerebbe al processore che si occupa di quella zona. Le forme più probabili per la divisione geometrica sono la griglia uniforme, oppure li alberi ottali (octree) o quadramentali(quadree). Le ultime due strutture sono struttue gerarchiche che servono a partizionare lo spazio in modo ricorsivo. Attenzione, se i pesci cambiano dominio serve che siano sincronizzati. Quadtree è la versione in 2D, octtree è la versione in 3D.

- ibrido, unione di più soluzioni

