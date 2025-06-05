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
Individual movement- 0.083644
Collective movement- 0.068647
Collective volitive- 0.096520

(8) 
Individual movement- 0.045884
Collective movement- 0.035155
Collective volitive- 0.092915

(12)
Individual movement- 0.037443
Collective movement- 0.050722
Collective volitive- 0.081791

(16)
Individual movement- 0.029870
Collective movement- 0.032265
Collective volitive- 0.084550

ORA FORSE SAREBBE IL CASO DI TESTARE NUOVAMENTE IL FUNZIONAMENTO DELL'ALGORITMO TRAMITE LE VISUALIZZAZIONI FATTE IN PYTHON.

# TEST 08/05/2025
Ora proviamo a vedere se funziona ancora tutto con l'aiuto delle visualizzazioni grafiche.
Siamo arrivati alla conclusione che in generale funziona, ma abbiamo un problemino con certe funzioni per cui la tendenza è quella di esplorare uno spazio troppo ampio e non capiamo se questo comportamento è corretto o se è causato da qualche problema nel nostro codice. Per provare a risolvere abbiamo normalizzato il vettore "delta" dentro "volitivePositionUpdateArray" con una funzione che non lo faccia in modo soft, ma purtoppo non ha funzionato. Per questo motivo abbiamo rinunicato a questa modifica e abbiamo costretto i pesci a stare dentro dei boundaries, se vanno fuori li rimettiamo al confine. Questo può essere giustificato dal fatto che vogliamo esplorare solo l'area interna a quei confini.


- prossima volta: 
Capiamo se vogliamo intraprendere la strada che ci porterà a sviluppare una versione poco leggibile, ma altamente parallelizzabile.
Se si, in bocca al lupo.
Se no, capiamo come immettere MPI.

# TEST 10/05/2025
Proviamo a inserire MPI nel codice, nuova directory v3_mpi_openmp. 
Per come la pensavamo molto tempo fa, noi definiamo quanti pesci totali ci saranno e poi il programma se li smazza automaticamente in schools facendo npesci / nprocessi.

Dato che ogni processo si gestisce un solo banco di pesci, dobbiamo ridefinire tutte le strutture che erano state fatte apposta per gestire n_schools (visto che ora ogni processo ne gestisce una sola).
Ora la versione funziona per l'individual movement, nonostante questo il tempo schizza su di molto molto e inoltre aumenta all'aumentare dei processi che si usano.
Alcuni dati:
RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 100 - DIMENSIONS 2 - MAX_ITER 200 - UPDATE_FREQUENCY 1 + 8 threads per processo ->>>>>>>Individual movement- 0.078004
RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 100 - DIMENSIONS 2 - MAX_ITER 200 - UPDATE_FREQUENCY 1 + 8 threads per processo ->>>>>>>Individual movement- 0.005673
RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 100 - DIMENSIONS 2 - MAX_ITER 200 - UPDATE_FREQUENCY 1 + 8 threads per processo ->>>>>>>Individual movement- 0.000015

Prossima volta: andare avanti con parallelizzazione e capire meglio quando usare scatter e pack e fare test

# TEST 15/05/2025
Abbiamo fatto un altro branch per provare a ridurre da 3 a 2 MPI_AllReduce in individual movement andando ad unire i due MPI_AllReduce che sommano le proprie componenti. Per farlo, però, dobbiamo andare a creare un nuovo vettore che comprende tutte le variabili.
Questo è stato fatto! Tendenzialmente riscontriamo un miglioramento sull'individual movement di circa 0.02 (ma non sempre, è un po' fluttuante... forse dipenderà dallo stato del cluster??)

Siamo al punto di adattare il collectiveVolitiveMovement, per farlo, dato che la funzione richiede comunicazione tra tutti i processi, vogliamo collassare insieme la funzione che si occupa di calcolare il baricentro e quella che si occupa di calcolare la somma totale dei pesi, in modo da ridurre la comunicazione totale. 
Con questa modifica funzionante abbiamo dei buoni risultati sulle tempistiche (negli esempi il numero di pesci totale è lo stesso 16000, con 8 core per processo)
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 20.770748 s
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 23.437921 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 21.707344 s
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 19.164604 s
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 1000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 36.393537 s (individual qua è molto alto, non so se sia per quello)


Nella parte di breeding, togliamo proprio la comunicazione tra tutti i banchi. Tanto il pesce peggiore tra tutti, sarà comunque il pesce peggiore del suo banco. Teniamo che ogni banco toglie il suo pesce peggiore e ne spawna uno nuovo dai suoi due pesci migliori.
Tempi di un altro test con breeding e avendo tolto le stampe (il numero di pesci totale è lo stesso 16000, con 8 core per processo):
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 68.865700 s
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 9.346894 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 18.798985 s
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 21.397372 s
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 1000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 53.934437 s

Consideriamo di usare MPI_IAllreduce, che è la versione non bloccante della stessa chiamata. Implementato qui male e da implementare meglio.
(il numero di pesci totale è lo stesso 16000, con 8 core per processo)
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 13.186210 s
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 21.615947 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 10.558105 s
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 14.17845 s
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 1000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 18.210253 s

-------------------------------------------------------------------------------------------------------------------------------
con 64000 pesci e 8 cores
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 64000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 120.599103 s
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 32000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 151.590752 s
  - place=pack ->  69.527637 s 
  - place=scatter ->  30.672575 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 103.626230 s
  - place=pack  ->   116.725785 s  
  - place=scatter -> 16.352849 s 
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 92.038892 s
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 117.446751 s
- RUNNING WITH: N-PROCESSES 32 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 64.285917 s

con 64000 pesci e 16 cores
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 64000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 94.166842 s
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 32000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 25.843306 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 24.979806 s
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 29.232648 s
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 31.953581 s
- RUNNING WITH: N-PROCESSES 32 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 39.903282 s

con 64000 pesci e 32 cores
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 64000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 92.584848 s
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 32000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 25.048108 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 23.050224 s
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 11.362145 s
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 8.297737 s
- RUNNING WITH: N-PROCESSES 32 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 100 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 16.806095 s
--------------------------------------------------------------------------------------------------------------------------------

con 64000 pesci e 16 cores (la differenza è che ci sono 1000 dimensioni invece che 100)
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 64000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 468.284195 s / 503.820818 s
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 32000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 358.015409 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 1000 - MAX_ITER 200  - UPDATE_FREQUENCY 1 => 264.799026 s
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 253.069131 s
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 201.931118 s
- RUNNING WITH: N-PROCESSES 32 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 97.150885 s

con 64000 pesci e 32 cores (la differenza è che ci sono 1000 dimensioni invece che 100)
- RUNNING WITH: N-PROCESSES 1 - N_FISHES_PER_PROCESS 64000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 442.861523 s               
- RUNNING WITH: N-PROCESSES 2 - N_FISHES_PER_PROCESS 32000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 288.869945 s
- RUNNING WITH: N-PROCESSES 4 - N_FISHES_PER_PROCESS 16000 - DIMENSIONS 1000 - MAX_ITER 200  - UPDATE_FREQUENCY 1 => 128.899041 s
- RUNNING WITH: N-PROCESSES 8 - N_FISHES_PER_PROCESS 8000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 57.821893 s 
- RUNNING WITH: N-PROCESSES 16 - N_FISHES_PER_PROCESS 4000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 54.816038 s 
  - place=pack ->  s
  - place=scatter -> 28.384023 s 
- RUNNING WITH: N-PROCESSES 32 - N_FISHES_PER_PROCESS 2000 - DIMENSIONS 1000 - MAX_ITER 200 - UPDATE_FREQUENCY 1 => 36.660738 s
  - place=pack    (3114370)
  - place=scatter 60.489892 s 
----------------------------------------------------------------------------------------------------------------------------------


Forse fare il mappazzone di codice non porterebbe a grandi risultati perché non ci permetterebbe di sfruttare appieno la funzionalita asincrona dell'mpi_Iallreduce... però controlliamo meglio la prossima volta.

Tenendo il numero di pesci fisso e modificando il numero di processi, vediamo che i tempi comunque si alzano... STRANO, uffi.

Proviamo anche ad aumentare i processi tenendo fisso il numero di pesci per processo. In questo caso guardiamo a come evolvono i tempi.
Con quelle premesse, i pesci crescono linearmente rispetto al numero di processi. Abbiamo vinto se i tempi non crescono anche loro linearmente :)

# 20/05/25
DA FARE:
- Se funziona lo script automatico per lanciare i job, fare la stessa cosa con i numeri corretti (Annina)
- Mergiare il branch (Des)
- Parser dei vari risulati ->
    cercare riga END per il tempo (uno a caso va bene)
    cercare rifa RUNNING WITH per tutte le info
    cercare riga per il numero di core

# 21/05/2025
Fatto uno script per lanciare job in massa:

- quando lo si lancia, vengono creati diversi file sh, uno per ogni job

- viene creato un file results.json che inizialmente serve ad avere una mappatura tra l'id del job e i suoi parametri

- una volta che si è contentx con i file ".o" e ".e" che sono stati prodotti dai job, puoi lanciare il parser.py che andrà ad aggiornare results.json con i tempi impiegati da ogni job (parser.py può essere sempre lanciato, anche se non tutti i file ".o" e ".e" sono stati prodotti. Ad ogni lancio di parser.py, il json viene aggiornato se ci dovessero essere nuovi file di output)

- Una volta che si sono raccolti i dati desiderati, basta lanciare lo script "remove_garbage.sh" per eliminare tutto quello che è stato prodotto nei precedenti punti

----------

Purtroppo, lanciare i job in massa produce dei risultati di merdaaaaaaaaaaaa [uffa uffa].
Credo che questo sia dovuto al fatto che si va a sovraccaricare il cluster.
Per ovviare a questo problema, potremmo fare in modo che ogni job venga lanciato solo quando il precedente ha prodotto i suoi file ".o" e ".e".
L'unica cosa che mi fa paura di questo approccio è che potrei essere bannato dal cluster perché ci sarebbe quello script che continua ad andare sul nodo di partenza... capiamola...

# 22/05/2025
Mi sono messo via lo script che lancia in automatico i job: non so per quale motivo ma, anche se lancio i job uno dopo l'altro, il tempo impiegato dal job è vergognosamente più alto rispetto a lanciarlo a mano uno per uno. 
Detto questo, augurateci buona fortuna: li lanceremo tutti a mano.

Teniamo fissi 64000 pesci totali

Abbiamo fatto un nuovo parser apposta per i file .sh.o{job_id} che aggiorna un file results.json con tutti i risultati che non sono ancora stati segnati.

# 26/05/25
Stiamo sviluppando i grafici per tirare le fila della nostra versione parallela.
Per ora abbiamo:

- un grafico che rappresenta la speedup all'aumentare dei processi, i dati sono raggrupati per numero di thread (visuals_2d_proc_speedup.py)

- un grafico che rappresenta il tempo in base all'aumentare dei core, i dati sono rappresentati per processo
Da notare che questi grafici servono per la strong scalability perchè aumenta il numero di processori ma i dati da usare sono sempre quelli.  (visuals_2d_tempo_core.py)

- un grafico che rappresenta l'efficiency in base al numero di core. I dati sono raggruppati per numero di processi. (visuals_2d_proc_efficiency) -> vediamo che l'efficiency decresce esponenzialmente all'aumentare dei thread: capiamo come cambia la curva all'aumentare dei pesci (tipo 128000)? Anche perché questo è di fatto un'analisi sulla strong scalability: il nostro algoritmo non sembra essere strongly scalable (se no avrebbe tenuto l'efficiency quasi costante all'aumentare delle risorse)


Cose da tenera a mente: scrivere tutte le scelte fatte sul report tipo come decidiamo di calcolare l'efficiency se per processo o per #thread*#processi. Poi fare studio con scatter e pack, fare studio sulla strong e weak scalability.

Attenzione: abbiamo modificato il calcolo della efficiency per tenere conto del numero di thread*processi così da avere una misura valida dell'occupazione di ogni core. Inoltre abbiamo trovato che per una versione ibrida, un aumento superlineare è normale.

RICORDIAMOCELOOO
Ora tutti i dati che abbiamo sono con "scatter", ma dovremmo anche prendere i dati con "pack".

TODO:
- sistemare pack e scatter (zio pera speriamo bene) -> capire perchè non parte pack con tante risorse
- finire di lanciare cose da 128000 per tirare delle conclusioni sensare -> bot
- REPORT e capire che cosa metterci dentro 

# 27/05/25
Abbiamo cambiato l'openMP scheduling da "dynamic,1" a "static,1" perchè il nostro carico è molto bilanciato

Ricordiamoci che tutti i test fatti finora sono per la strong scalability (che, spoiler, non abbiamo eheh).
Rifacciamo tutti i test tenendo fissi i fishes per CORE!
Per farlo, questi sono i parametri dei test:
- select: [1,2,4,8,16,32,64] -> da variare ad ogni esecuzione
- core: 2 (o 4, ma se lo facciamo con 2 possiamo motivare la scelta dicendo che, guardando i grafici dell'efficienza, con 2 thread questa è la migliore)
- 2000 pesci PER CORE (con due thread saranno 4000 pesci e con 4 thread saranno 8000. Questo numero è da lasciare fisso!)
- dimensioni, iterazioni e update frequency le lasciamo come al solito

E' già stato fatto uno script che devrebbe plottare le efficiency per la weak scalability, ma ritestiamolo una volta raccolti i dati qui sopra

Ricordiamoci anche di raccogliere gli stessi dati che abbiamo raccolto con scatter e 64000 pesci, ma ora con pack e sempre 64000 pesci.



# 28/05/25
Completati i test per la weak scalability, per fare una cosa sensata e dare un tempo per il calcolo dello speedup ho preso fatto 1 core 1 processo e 2000 pesci, spero abbia senso altrimenti si può cambiare

# 03/06/2025
Nella weak scalability cambia un po' il modo in cui ci rapportiamo alle metriche:
- Lo speedup non ha più senso (al massimo si può calcolare uno speeddown)
- La efficiency va calcolata in un altro modo: tenendo sempre la dimensione proporzionale al numero di core, si fa la frazione tra il tempo registrato con (1,1) e quello con (nProc, nThread).

L'idea è che se i tempi rimangono più o meno costanti al variare delle risorse e della dimensione del problema, allora si ha weak scalability.
A noi pare di averne un accenno eheh. si gode (un pochino, non troppo :)

Iniziamo a mettere sul latex la tabella con i dati e i vari plot per openmp e la versione ibrida

# 04/05/2025
TODO:
- finire di raccogliere i dati (o controllare che ci siano già)
  - weak scalability MPI
  - strong scalability MPI
- script per weak openMP (2 gruppi poi scegliamo quello che va meglio)
- overleaf -> leggi commenti ovunque e correggi quello che vuoi
- rimpinguare gli altri dati
- analisi dei risultati che abbiamo

# 05/05/2025
Grande refactory della repository (hybrid e openmp_only funzionano)
Cose da sistemare del refactory:
- colleganti grafi corretti
- mpi_only per seeds e tutto (è un po un casino)