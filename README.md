# hpc_project
Io e Annina che facciamo muovere i pesci eheh

### FONTI

- https://en.wikipedia.org/wiki/Fish_School_Search
- https://it.wikipedia.org/wiki/Swarm_intelligence -> FSS è un algoritmo di Swarm Intelligence
- https://www.mql5.com/en/articles/11841 -> c'è dello pseudocodice in C (non ho letto bene tutto)
- https://www.youtube.com/watch?v=kSoA_zzrqsc -> Video spiegazione

## CRITICITÀ

- Gestione dell'aggiornamento del peso -> divisione per il max_improvement: abbiamo visto che si comporta meglio senza quella divisione. Spesso infatti sono valori piccoli e quando quel valore è 0 non dividiamo per nulla. [SOLVED]

## TODO
- [X] Controllare che la scrittura del json funzioni per n dimensioni
- [X] cambiamo la v1 sequenziale per far comunichre i banchi tra di loro ogni n iterazioni

## NOTE IMPORTANTI SULLE NOSTRE SCELTE IMPLEMENTATIVE 
- Abbiamo deciso di aggiornare i pesi dei pesci solo dopo il loro individual movement.
  In ogni caso, questo non influisce negativamente perché ad ogni individual movement verrà incluso anche il collective movement del ciclo prima.

- Il peso di un pesce non è direttamente correlabile alla sua posizione nella funzione (non è vero che un pesce più vicino al minimo abbia sempre un peso più grande di quello di un pesce distante), ma tuttavia è direttamente correlabile al suo miglioramento all'interno della funzione. 
  Se un pesce migliora costantemente la sua fitness, il suo peso aumenta. Tuttavia:
  Se la sua fitness si stabilizza o peggiora, il peso rimane invariato o diminuisce.
  Pesci con miglioramenti recenti ma con fitness assoluta più bassa possono avere pesi superiori.

- Abbiamo trovato in varie altre implementazioni che il baricentro considera la media delle posizioni pesate sulla weight dei pesci, quindi va diviso per il peso e non la posizione (come invece scritto nel paper2)

- Abbiamo notato che, per cercare di arginare l'overhead di comunicazione tra processi con MPI, possiamo dividere il baco in sottobanchi autonomi.
In ogni caso, vogliamo mantenere la comunicazione tra i vari sottobanchi: proprio a causa dell'overhead della comunicazione MPI, ogni N epoche ogni sottobanco condividerà il proprio stato affinchè sia comunque possibile il collective ed il volitive movement tra tutti i pesci.
ATTENZIONE però: se la comunicazione avviene troppo sporadicamente, l'aggiornamento causato dal collective e dal volitive movement generale non avrà un impatto significativo. iPensiamo sia il caso di svolgere la comunicazione collettiva ogni 5 epoche circa.

- Mettere sempre l'end del MPI_Wtime DOPO MPI_finalize() altrimenti da 0

## COSA POTER METTERE NEL REPORT
- Differenze tra le nostre versioni (grafici)
- Evoluzione delle tempistiche di ogni iterazione in un ambiente parallelo (vogliamo far vedere i picchi di tempo quando ci sono le mpi_allreduce)


## PARALLELIZZAZIONE

Diverse possibilità (by chatgpt):

- parallelizzare il calcolo della fitness (cosa più pesante da fare), secondo me un po troppa comunicazione

- proposta di fare update non-sincronizzati, non aspettare tutti (molto time consuming)

- parallelizzare l'individual movement (il collective poi va fatto insieme)

- dividere lo spazio in vari domini, parallelizzare il calcolo per dominio (dubbio: come ci si assicura che siano sempre tutti separati?)

  Ci sono vari modi per fare questo, considerando *subset di pesci*- *population based,* quindi ogni processore mantiene un set di pesci, senza considerare dove si trovano. *Geometric domain decomposition* nel caso in cui si consideri una divisione in celle, in base allo spazio. In questo caso i pesci potrebbero cambiare zona e in quel caso la loro gestione passerebbe al processore che si occupa di quella zona. Le forme più probabili per la divisione geometrica sono la griglia uniforme, oppure li alberi ottali (octree) o quadramentali(quadree). Le ultime due strutture sono struttue gerarchiche che servono a partizionare lo spazio in modo ricorsivo. Attenzione, se i pesci cambiano dominio serve che siano sincronizzati. Quadtree è la versione in 2D, octtree è la versione in 3D.

- ibrido, unione di più soluzioni


# Nostre versioni


## Sequenziali
### v0
La versione più semplice possibile.


Abbiamo un solo banco e tutti i pesci si aggiornano dopo ogni epoca.

### v1
ATTENZIONE: calcolo tempo non con MPI_Wtime
Versione più semplice che ha l'intento di essere corrispondente alle v1 parallela, quindi il banco non è unico ma ci sono più banchi definibili. Ciò implica che ogni banco si aggiorna internamente ad ongi iterazione e tutti i banchi si aggiornano insieme secondo l'UPDATE_FREQUENCY.

### v1.1
Equivalente alla v1, ma il tempo viene misurato correttamete e la versione è fatta per accettare parametri in input dal runner dello scheuduler del cluster, in modo da poter runnare le metriche in modo veloce e pratico direttamente da li. Per questo stesso motivo viene scritto su file il risultato in questo formato RUNNING WITH: N-SCHOOLS 20 - N_FISHES_PER_SCHOOL 200 - DIMENSIONS 1000 - MAX_ITER 100 - UPDATE_FREQUENCY 1
56903.189000

## Parallele
### v0
Come v0 sequenziale, con qualche MPI_allReduce e OMP parallel for un pò a caso. 
ATTENZIONE: tempo calcolato male non con MPI_Wtime ma con gettimeofday()
### v1
come V1 sequanziale, idea di affidare ogni banco a un core diverso e farli comunicare ogni UPDATE_FREQUENCY iterazioni
ATTENZIONE: tempo calcolato male non con MPI_Wtime ma con gettimeofday()

### v1.1
Stessa idea di V1.1 con tempo calcolato correttamente

### v2_multithreaded only
Idea di usare SOLO multithreading per parallelizzare le cose
DUBBIO: sarebbe forse da fare un unica versione in cui si può fare sia multithreading sia parallelizzazione ibrida?

Nei pragma, non serve specificare le variabili const e globali come shared.
Tuttavia, abbiamo avuto dei problemi quando dichiaravamo delle variabili come globali per poi inizializzarle nel main. 
Per ovviare a questo problema, abbiamo dovuto ripiegare sul passare le variabili precedentemente globali come parametri delle funzioni.


# TEST 28/04/2025
Vogliamo capire più o meno quanto i tempi di esecuzione si rimpiccioliscano all'aumentare dei thread (senza mulutiprocessing... per ora).
Per fare questi test, usiamo unicamente la versione solo multithreading andando a tweakare il numero di thread ed il numero di core che richiediamo.

STIAMO SEMPRE ANDANDO AD OPERARE CON UN SOLO PROCESSO
N_SCHOOLS - N_FISHES_PER_SCHOOL - DIMENSIONS - MAX_ITER - UPDATE-_FREQUENCY
1 - 5000 - 1000 - 100 - 1

Con 1 core e 1 thread (sequenziale)
91.32s

Con 2 core e 2 thread
83.02s

Con 4 core e 4 thread
156.00s

Con 8 core e 8 thread
247.77s

Qui c'è qualcosa che ci puzza. Noi ci aspettavamo questi tempi proporzionalmente decrescenti al numero di core che usiamo.
Non stiamo ben capendo qual è il motivo dietro a tutto questo, ma troviamo che la combinazione tra omp_places=cores e omp_proc_bind=close sembra dare buoni risultati.

Abbiamo poi aggiunto omp_scheule="dynamic" e "dynamic,1"

Forse ci starebbe usare la clausola di omp "collapse(2)" per parallelizzare due loop nestati in modo furbo.


# TEST 29/04/2025
Abbiamo visto che l'individual_movement è la parte più lenta di tutto il processo. A seguire collective movement e volitive movement. 

C'è la possibilità che le chiamate a rand() che facciamo serializzino le chiamate in quanto vanno ad agire su una variabile condivisa.
Per ovviare a questo, pensiamo sia il caso di usare rand_r() con un seed diverso per ogni thread.

Poi, per velocizzare, potremmo usare delle strutture dati locali e poi, solo alla fine, andare a scriverle su quelle condivise tra i thread.

Un'idea potrebbe essere controllare bene dove parallelizziamo e in caso capire quanti threads ha senso far spawnare.

Forse, in modo un po' più radicale, stiamo pensando di riscrivere la logica in modo che sia tutto più efficiente (evitando mille loop che ora servono solo ad avere codice più leggibile)



# TEST 05/05/2025
Il rand_t() ora sembra funzionare ma, tra i test con un processo e 1/2/4/8/12 thread, vince sempre il test ccon 1 thread. Qui qualche risultato.

(1)
Variablee reset- 0.000004
Individual movement- 0.312136
Update weights- 0.000066
Collective movement- 0.238642
Collective volitive- 0.170157
Breeding- 0.000067
END: 72.790821

(2)
Variablee reset- 0.000005
Individual movement- 0.349881
Update weights- 0.000099
Collective movement- 0.275390
Collective volitive- 0.269453
Breeding- 0.000078
END: 137.138517

(4)
Variablee reset- 0.000003
Individual movement- 0.340076
Update weights- 0.000071
Collective movement- 0.256104
Collective volitive- 0.218467
Breeding- 0.000069
END: 80.549758

(8)
Variablee reset- 0.000004
Individual movement- 0.425280
Update weights- 0.000339
Collective movement- 0.326726
Collective volitive- 0.290037
Breeding- 0.000080
END: 98.749366

(12)
Variablee reset- 0.000006
Individual movement- 0.773514
Update weights- 0.000111
Collective movement- 0.608344
Collective volitive- 0.554615
Breeding- 0.000115
END: 145.290174


With usign a local structure in the individual movement we reached a better result:
(1)  Individual movement- 0.325134
(2)  Individual movement- 0.192522
(4)  Individual movement- 0.199572
(8)  Individual movement- 0.070397
(12) Individual movement- 0.033476
(16) Individual movement- 0.026439
(20) Individual movement- 0.021122
(24) Individual movement- 0.023564 <- inizia a crescere leggermente, COME MAI?
(30) Individual movement- 0.027301


Ora proviamo a rendere locali anche le variabili del collective:



(1)  
Individual movement- 0.322119
Collective movement- 0.264698

(2)  
Individual movement- 0.374777
Collective movement- 0.297142
-----------------------------
Individual movement- 0.172758
Collective movement- 0.137967

(4)  
Individual movement- 0.087109
Collective movement- 0.068994

(8) 
Individual movement- 0.056808
Collective movement- 0.076849
-----------------------------
Individual movement- 0.106309
Collective movement- 0.133865
-----------------------------
Individual movement- 0.100382
Collective movement- 0.127348

(12)
Individual movement- 0.067111
Collective movement- 0.102759

(16)
Individual movement- 0.030886
Collective movement- 0.045614

(20)
Individual movement- 0.026262
Collective movement- 0.014225

(24)
Individual movement- 0.016126
Collective movement- 0.018476

... e del collective volitive:
(1)  
Individual movement- 0.338221
Collective movement- 0.276688
Collective volitive- 0.089173

(2)  
Individual movement- 0.162689
Collective movement- 0.134987
Collective volitive- 0.079497

(4)  

(8) 
Individual movement- 0.045884
Collective movement- 0.035155
Collective volitive- 0.092915

(12)

(16)
