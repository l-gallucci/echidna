# Echidna — wiki

Guida al toolkit per controllo di assembly metagenomici, qualità dei MAG e screening delle chimere.
Questa pagina descrive le correzioni successive al checkout 0.0.4 esaminato il 3 novembre 2026; non annuncia una nuova release pubblicata.

## Due percorsi

| Percorso | Input | Scopo |
|---|---|---|
| Assembly | Assembly e reads HiFi rimappate sullo stesso assembly | Cercare misassemblies, separare contig lunghi, preparare il binning |
| MAG | Un FASTA per bin | Completezza, contaminazione, tassonomia, RNA e screening delle chimere |

Il core usa solo la libreria standard Python. Anvi’o, CheckM2, GTDB-Tk, GUNC, MMseqs2, barrnap, tRNAscan-SE, samtools e BLAST richiedono installazioni e database separati per i comandi che li utilizzano. Gli ambienti conda del progetto aiutano a separare le loro dipendenze.

Installare Echidna nell'ambiente che esegue ciascun wrapper. Da un checkout:

```bash
python -m pip install -e .
echidna --version
echidna --help
```

Il requisito dichiarato del core è Python ≥3.8; la compatibilità dei tool esterni va verificata separatamente. Le correzioni sono testate con l'interprete Miniforge disponibile; non implicano una certificazione di tutti i sistemi operativi o ambienti.

## Identità e formati

Usare nomi dei MAG e contig univoci. BAM e FASTA devono avere gli stessi identificativi e lo stesso riferimento. I read ID devono corrispondere al file del classificatore. Prima di rinominare i MAG completare i controlli e conservare gli ID originali dei report.

I parser riconoscono sia gli schemi nativi di anvi’o 9 sia le vecchie colonne Echidna:

```text
clipping: contig length pos relative_pos cov clipping clipping_ratio
zero-cov: contig length range range_size
```

Nel formato nativo, `clipping` è il numero di reads clippate; `relative_pos` non è una copertura. `range` è un intervallo zero-based half-open. File con sola intestazione sono risultati validi senza eventi; file mancanti o malformati sono errori. [Formato ufficiale anvi’o](https://anvio.org/help/9/programs/anvi-script-find-misassemblies/).

`process-large-contigs` scrive `contig_name_map.tsv` per collegare i nomi CheckM2 con hash ai contig originali. Un contig senza riga CheckM2 viene conservato come non valutato. `remove-hq-contigs --hq-dir` legge i contig dai FASTA dei MAG; i numeri nei nomi dei file non sono identificativi dei contig. In `--hq-list` usare i nomi FASTA esatti.

## Ramo assembly

Ordine delle dipendenze:

```text
assembly originale + BAM delle stesse reads
  → anvi’o: clipping e zero-coverage
  → process-large-contigs
  → check-read-chimeras, se è disponibile la classificazione delle reads
  → classify-clipping: annotazione report-only
  → filter-misassemblies
  → nuovo mapping, depth e binning esterno
```

Anteprima del wrapper, senza eseguire i tool:

```bash
echidna misassembly-wf --assembly assembly.fa --bam hifi.bam \
  --output-dir assembly_qc --dry-run
```

Il wrapper inoltra ora il FASTA al controllo read-level e il produttore scrive una tabella delle finestre anche quando vuota. `classify-clipping` legge il campo `is_taxon_shift` prodotto da quel controllo. Il suo report non decide automaticamente dove tagliare: il filtro usa i clipping originali, eventualmente selezionati con `--min-clip-coverage`.

La profondità locale è ora una media per base nella finestra, calcolata dai blocchi CIGAR allineati, anziché un conteggio delle partenze. Reads lunghe che iniziano fuori dalla finestra sono incluse se la attraversano. Le statistiche delle coppie sono NA per reads non appaiate come HiFi; i tag SA delle primarie forniscono il segnale di join senza contare due volte le supplementary.

Dopo lo splitting, per esempio `contig_A → contig_A_split1 + contig_A_split2`, il vecchio BAM usa ancora `contig_A`. Rimappare le reads sul FASTA modificato prima di produrre la depth per il nuovo binning. Il BAM originale continua a servire ai controlli sul riferimento originale.

`HQ_extracted/` e `assembly_for_filtering.fa` sono output complementari: i contig estratti non si trovano più nell'assembly residuo. Conservare e valutare entrambi. Un assembly senza contig lunghi produce comunque l'assembly per il passo successivo.

## Ramo MAG

```bash
echidna mag-quality-wf --bins-dir bins --genome-ext fa \
  --output-dir mag_qc --tiers HQ,MQ \
  --gtdbtk-data /data/gtdbtk_release \
  --gunc-db /data/gunc/reference.dmnd \
  --mmseqs-db /data/mmseqs_gtdb/gtdb \
  --split-by-domain --dry-run
```

Togliere `--dry-run` solo dopo aver predisposto ambienti, input e database. Il piano è:
`checkm2 → organize1 → gtdbtk → gunc → [mmseqs → contigtax] → mimag → chimeras → organize2`.

LQ non viene portato avanti con il default `--tiers HQ,MQ`. Specificare i tier desiderati. Le estensioni vengono inoltrate anche ai passi MIMAG e chimera. I nomi rimangono quelli originali fino all'ultima organizzazione.

I controlli chimera sono **report-only**: un MAG segnalato può comunque comparire nelle cartelle di qualità. Scegliere esplicitamente quali escludere e quali rivedere. `check-ssu-taxonomy` e `final-report` restano comandi separati; il DAG Snakemake fornito non comprende MMseqs/classify-contigs e non equivale al wrapper con tutte le opzioni.

## Database concordanti e confronto delle reads

Per confrontare GUNC, MMseqs2, GTDB-Tk e SSU, preferire riferimenti derivati dalla **stessa tassonomia e dalla stessa release**, quando disponibili. Usare una base tassonomica comune riduce differenze di nomenclatura e copertura; non rende indipendenti i due metodi.

Non si può passare lo stesso file `.dmnd` a MMseqs2: i tool richiedono formati e preparazioni differenti. Usare una release GUNC appropriata e il corrispondente database MMseqs2; registrare release, provenienza e checksum. Eseguire entrambi i passi produce i due segnali: non serve eseguire due volte lo stesso classificatore per ottenere concordanza. Database alternativi possono essere usati in run separate per un'analisi di sensibilità.

Per le reads si può aggiungere un confronto a rango comune:

```bash
echidna check-read-chimeras --bam hifi.bam --assembly assembly.fa \
  --kraken2-output reads.kraken2 --kraken2-report reads.report \
  --taxonomy-rank order --nodes-dmp /data/kraken/taxonomy/nodes.dmp \
  --output-dir reads_qc
```

`nodes.dmp` deve provenire dal database del classificatore. Senza `--taxonomy-rank` rimane il confronto dei taxid originali. Il roll-up riduce falsi conflitti fra un taxid di specie e uno del suo genere/ordine. Sono richieste almeno 10 reads classificate per il rischio read-level; reads mancanti dal classificatore non diventano un taxon artificiale. Le finestre valutano reads classificate.

Ulteriore evidenza utile: supporto SA coerente, identità e qualità del mapping, supporto dei due lati della giunzione e verifica in un viewer. Una tassonomia mista è il livello di screening disponibile, non una dimostrazione definitiva di chimera.

## Qualità, circolarità e suffissi GTDB

Le soglie interne rimangono intenzionalmente:

| Tier | Completezza | Contaminazione |
|---|---|---|
| HQ | ≥90% | ≤5% |
| MQ | ≥70% | ≤10% |
| LQ | ≥50% | ≤10% |

La selezione “HQ” di `process-large-contigs` usa invece 50/10 per default. Queste etichette operative restano distinte dallo standard MIMAG. Se si richiede `--mimag-rna-dir`, un HQ senza evidenza RNA completa, compreso un risultato mancante, viene declassato a MQ; il valore RNA rimane vuoto quando non valutato.

Il report distingue:

- `circular_annotated`: annotazione esplicita nel FASTA, per MAG a singolo contig; NA se assente. Una descrizione `circular=true` viene conservata nell'estrazione dei contig lunghi. `circular=false` prevale sull'euristica di lunghezza nel filtro.
- `near_complete_candidate`: completezza ≥99% e contaminazione ≤5%, indipendentemente dalla circolarità. È un'evidenziazione di qualità, non la dichiarazione di genoma completo.
- `n_contigs`: numero di contig nel FASTA. Queste informazioni strutturali richiedono `final-report --genomes-dir`.

Nel filtro assembly l'euristica di lunghezza resta un fallback per l'opzione `--split-hq-circular` quando manca l'annotazione. Non viene presentata come prova di circolarità; le definizioni degli assembler possono differire.

I confronti SSU **conservano `_A`, `_B` e gli altri suffissi GTDB**. Nessun nome FASTA viene cambiato da questa correzione. Un confronto prima positivo solo grazie alla rimozione del suffisso può ora diventare negativo: va valutato con la release del riferimento. Per SILVA vs GTDB non assumere una corrispondenza automatica fra nomi diversi; preferire riferimenti coerenti o una mappa tassonomica curata.

## Verdetti e report finale

| Evidenza | Verdetto MAG |
|---|---|
| GUNC fail e frazione estranea ≥0.05 | Chimeric |
| Un solo segnale anomalo | Check |
| Nessuna anomalia nei segnali misurati | Clean |
| Nessun segnale misurato | Unknown |

La frazione estranea è pesata sulla lunghezza dei contig **risolti**. Zero contig risolti significa segnale mancante, non frazione misurata uguale a zero. Per compatibilità sono accettate anche vecchie tabelle senza `n_contigs_resolved`, ma queste non permettono il controllo della disponibilità: preferire output aggiornati di `classify-contigs`.

GC, contaminazione CheckM2 e warning GTDB sono informativi. Un solo segnale misurato e negativo può dare Clean: leggerlo come “nessuna anomalia in quel segnale”. La regola di consenso è invariata.

```bash
echidna final-report --genomes-dir bins \
  --checkm2-report mag_qc/01_checkm2/quality_report.tsv \
  --gtdbtk-dir mag_qc/03_taxonomy --mimag-rna-dir mag_qc/05_mimag \
  --chimera-report mag_qc/06_chimera/chimera_report.tsv \
  --organized-quality mag_qc/07_final/merged_quality.tsv \
  --output mag_qc/final_report.tsv
```

`tier` mantiene il significato storico: soglie CheckM2 predefinite. `organized_tier` riporta la selezione effettiva, comprese soglie personalizzate e gate RNA, quando si passa `--organized-quality`. I genomi del report vanno identificati con i nomi originali, non con i nomi rinominati dell'output finale.

Le colonne aggiunte non rinominano le precedenti; i consumatori TSV devono leggere le intestazioni. `concerns` vuoto non certifica la qualità di un MAG. BLAST mantiene il criterio del workflow Trigodet/anvi’o: un'esecuzione riuscita senza hit può dare `unsupported`, un errore di esecuzione ora interrompe il comando.

## Ripresa e verifiche

Usare directory di output nuove per run con input o soglie cambiati. Per evitare tier obsolete, `organize-mags` rifiuta una selezione già popolata; non elimina file esistenti. `--from` presuppone che gli output dei passi precedenti siano validi. Per rifare l'organizzazione finale usare un output nuovo o spostare prima la precedente selezione.

La distribuzione si chiama **`echidnaxomics`**; comando e import restano **`echidna`**.
Fino alla pubblicazione installare dal checkout, non dal pacchetto omonimo `echidna`.
I wheel includono ora i template Snakemake: `echidna setup-workflow --export-workflow ./my-workflow`
li copia in una directory nuova e modificabile. `--repo /percorso/del/checkout` resta disponibile.
Con un workflow esportato passare a `fetch-databases` il percorso esplicito `--config`.
Gli ambienti forniti sono esempi da verificare con la propria release dei tool.

Il DAG traccia i singoli FASTA e le directory di output dei tool; nelle riesecuzioni Snakemake
ricrea gli output dichiarati, compresa `organized/Selected`, per evitare tier obsolete.
Conservare in queste directory solo gli output del workflow. Le prove del DAG usano report
a monte simulati e i comandi Echidna reali a valle: non sono una validazione con database reali.
Metadati autore, DOI e ricetta Bioconda sono descritti nella
[guida alla release](https://github.com/l-gallucci/echidna/blob/main/docs/releasing.md).

```bash
PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider
```

I test includono casi puliti di controllo, collegamenti fra formati reali e uno smoke test opzionale con un BAM sintetico e samtools reale. Non sostituiscono una run completa delle pipeline con i database installati. Il registro delle correzioni e delle scelte mantenute è nella [relazione del repository](https://github.com/l-gallucci/echidna/blob/main/docs/audit-2026-11-03/FIXES.md).
