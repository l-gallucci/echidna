> Release preparation and subsequent checks: [release-validation.json](release-validation.json) and [release guide](../releasing.md). The previous results below describe the first correction pass.

# Rivalidazione e correzioni — 3 novembre 2026

Il precedente AUDIT.md è una fotografia del commit `de465bcd`, non un giudizio su ogni run effettuata. Un risultato Clean con segnali realmente misurati può essere corretto. I difetti riprodotti riguardano trigger specifici; non dimostrano che i MAG delle analisi precedenti siano sbagliati.

Le riproduzioni storiche restano allegate con i loro esiti sul vecchio codice. Dopo le correzioni non devono continuare a “passare”: cercano deliberatamente il comportamento difettoso. La regressione corrente è in `tests/test_pipeline_contracts.py`.

| Rilievo | Rivalidazione e intervento |
|---|---|
| F01 | Confermato su FASTA + report CheckM2 con hash + clipping nativo. Le decisioni traducono gli ID e scrivono `contig_name_map.tsv`. |
| F02 | Confermato per assembly senza contig lunghi. Vengono prodotti assembly residuo e directory HQ anche in quel caso. |
| F03 | Confermato con report parziale. I mancanti restano nell'assembly come non valutati; report estranei/duplicati sono rifiutati. |
| F04 | Confermato nel piano con tassonomia. Il wrapper inoltra l'assembly, le tabelle finestre esistono anche senza righe. |
| F05 | Confermato con producer/consumer effettivi. Il lettore accetta `is_taxon_shift` e il vecchio alias `is_shift`. |
| F06 | Confermato su CIGAR di read lunga. Il calcolo usa overlap dei blocchi allineati e profondità per base; filtri coerenti con la media. |
| F07 | Confermato su reads flag=0/16. Statistiche paired disponibili solo per reads appaiate. |
| F08 | Confermato simulando exit BLAST 2. L'errore interrompe l'analisi, mentre zero hit dopo successo conserva la regola biologica precedente. Controllati anche return code samtools. |
| F09 | Comportamento precedente intenzionale, ma incoerente con un gate RNA richiesto. Risultati mancanti declassano HQ→MQ; non sono trasformati in RNA presenti. Soglie qualità invariate. |
| F10 | Confermato rilanciando l'organizzazione. Una cartella Selected popolata viene rifiutata con indicazione di usare un output nuovo; nessuna cancellazione automatica. Anche il riuso dell'assembly già processato viene rifiutato. |
| F11 | Confermato con zero contig risolti. Il segnale diventa assente; senza GUNC il verdetto è Unknown. Frazioni malformate causano errore. Vecchie tabelle senza colonna di disponibilità restano supportate. |
| F12 | Confermato con taxid=-1. Le reads non corrispondenti non contano come classificate; il rischio richiede almeno 10 reads classificate. Aggiunto roll-up opzionale con `--taxonomy-rank` e `--nodes-dmp`. |
| F13 | Confermato usando MAG multi-contig. La rimozione usa gli ID esatti letti dai FASTA, non il numero del bin; le liste devono contenere nomi esatti. |
| F14 | La colonna tier CheckM2 era una scelta documentata, non un bug aritmetico. Rimane compatibile; `--organized-quality` aggiunge `organized_tier` per evidenziare la tier effettiva. |
| F15 | Confermato con regex `\d` e virgolette. Replacement callable e serializzazione JSON compatibile YAML conservano il valore. |
| F16 | Confermato con `.fna`. Estensione inoltrata a MIMAG e controllo chimera. |
| F17 | Incompatibilità confermata, ma non basta tradurre Chimeric→High: un bin non è un contig. Il filtro ora richiede esplicitamente il report per-contig delle reads e rifiuta quello per-MAG invece di ignorarlo. |
| F18 | Confermato con file header-only. Il filtro continua e la classificazione produce un report vuoto valido. File assenti e formati illeggibili restano errori. |
| F19/F20 | Confermati contro lo schema ufficiale anvi’o 9. Parser per intestazioni, numero clipping corretto e range nativi zero-based; mantenuta la compatibilità legacy. |
| F21 | Limite di distribuzione, non algoritmo errato. Si mantiene il requisito di checkout per il setup Snakemake, ora esplicito in installazione e wiki; non si promette un workflow incluso nel package normale. |
| D01 | Scelta mantenuta: screening chimera e classificazione clipping report-only. Non viene aggiunto un filtro scientifico automatico. |

## Scelte metodologiche concordate

1. Soglie e nomi delle tier restano interni al progetto, distinti da MIMAG rigoroso. Nessun cambiamento alle bande 90/70/50 o ai confini inclusivi.
2. Regola di consenso GUNC + tassonomia contig invariata. Suggeriti database della stessa tassonomia/release, preparati nei formati di ciascun tool. Confronti alternativi vanno eseguiti in run separate.
3. Screening read-level mantenuto. Aggiunti rango comune opzionale e controllo sulla quantità di evidenza; suggerita verifica della giunzione con mapping e SA.
4. Annotazioni `circular=true/false` conservate e lette. Il report distingue circolarità annotata da candidato quasi completo (≥99% e ≤5% contaminazione). Lunghezza resta fallback nel filtro, non prova di circolarità.
5. Regola BLAST del workflow citato mantenuta. Un match locale non certifica ogni base, ma il software non cambia metodo: distingue gli errori tecnici dall'esito negativo.
6. Suffissi GTDB preservati. È un cambiamento intenzionale del confronto SSU; può far emergere disaccordi prima nascosti, senza cambiare nomi dei genomi. SILVA/GTDB richiedono cautela nella corrispondenza dei taxa.
7. Dopo splitting rinominare non basta: BAM/depth devono essere rigenerati sul FASTA modificato per il nuovo binning. Corretto l'esempio della guida.

## Verifica e limiti

La suite include controlli che mantengono Clean quando GUNC e tassonomia sono misurati e negativi, e mantiene Chimeric quando entrambi sono positivi. I test dei subprocess usano mock dichiarati, più uno smoke test opzionale con samtools reale. In questo ambiente samtools 1.21 è disponibile: su un BAM sintetico a copertura nota, il codice originale misura erroneamente profondità locale 0, il codice corretto misura 10× e riconosce 10 clipping finali con SA. Il confronto è in `real-bam-comparison.json`; il log finale è in `validation-current.json`.

Non è stata eseguita una validazione end-to-end con i database pesanti o sul cluster. I risultati storici dei propri dataset possono essere rivalutati selettivamente: schema dei file in input, corrispondenza degli ID, disponibilità dei segnali e warning. I test non consentono di certificare retroattivamente tutti i MAG.
