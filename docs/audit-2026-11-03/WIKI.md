> Historical snapshot of commit de465bcd. For the revalidation and current fixes see [FIXES.md](FIXES.md); current wiki source: [Home.md](../../wiki/Home.md).

# Echidna — wiki tecnica e guida operativa

Pagina indipendente dal README, pronta da importare in una wiki Markdown.
Aggiornamento: **3 novembre 2026**. Riferimento: Echidna **0.0.4**, commit `de465bcd7dc47cba7be1e157f1af5e97f38bc0ba`.

> Stato verificato: 179 test del repository passano, ma sono presenti problemi di integrazione e interpretazione degli esiti. Le pipeline complete con tool e database reali non sono state validate in questa revisione. Consultare [l'audit](AUDIT.md) prima di usare gli output per selezioni definitive.

## A cosa serve

Echidna coordina controlli di qualità su assembly metagenomici e MAG. Il package usa la libreria standard Python; i calcoli biologici più pesanti sono affidati a tool esterni. L'installazione del core non installa automaticamente questi tool né i relativi database.

Ci sono due percorsi distinti:

| Percorso | Input | Risultato previsto |
|---|---|---|
| A — misassemblies | Assembly, BAM delle stesse reads HiFi usate per assemblarlo, eventuale tassonomia delle reads | Assembly per il binning, contig estratti e report delle anomalie |
| B — qualità MAG | Un FASTA per bin, database di qualità/tassonomia | Metriche CheckM2, GTDB-Tk, GUNC, RNA e selezione organizzata per tier |

Binning, mapping e validazione finale dei genomi restano attività esterne. Le cartelle di MAG organizzati possono contenere genomi segnalati come chimere: lo screening non li esclude automaticamente.

## Struttura del progetto

| Percorso nel repository | Responsabilità |
|---|---|
| `echidna/__main__.py` | Parser principale, registrazione dei 22 comandi, gestione errori CLI |
| `echidna/commands/` | Parser dei sottocomandi e coordinamento delle operazioni |
| `echidna/utils/` | FASTA, parser dei report, identità, statistiche e regole decisionali |
| `workflow/Snakefile` | DAG Snakemake del ramo MAG, con ambienti separati |
| `workflow/envs/` | Ricette conda fornite dal progetto |
| `tests/` | Suite di 179 casi nella revisione esaminata |
| `docs/` | Installazione, workflow, HPC, reference e fonti |

`filter-checkm2` è l'alias deprecato di `organize-mags`; il modulo Python si chiama ancora `filter_checkm2.py`.

## Installazione e ambienti

Il metadata dichiara Python ≥3.8, ma le guide e `check-envs` non sono ancora allineati. In questa revisione i test sono stati eseguiti su Python 3.13.13, gli help su 3.14.7. La compatibilità runtime con 3.8 non è stata provata. Scegliere per ogni tool un ambiente supportato dal tool stesso; evitare di dedurre la compatibilità dei wrapper dal solo import del package.

Da un checkout del progetto:

```bash
python -m pip install -e .
python -m echidna --version
python -m echidna --help
```

Per gli ambienti dei tool usare le ricette del repository e la guida `docs/installation.md`. Echidna deve essere installato anche nell'ambiente che esegue un suo wrapper. I comandi `*-wf` possono invocare `conda run -n <ambiente>`, mentre quelli standalone richiedono il tool nel proprio PATH.

| Ruolo | Tool / risorse richieste |
|---|---|
| Operazioni FASTA e report | Solo core Python |
| Misassemblies da mapping | anvi’o; samtools per le statistiche Echidna |
| Completezza/contaminazione | CheckM2 e database DIAMOND |
| Tassonomia del MAG | GTDB-Tk e directory dati compatibile |
| Coerenza genica | GUNC, DIAMOND, Prodigal, database GUNC |
| rRNA/tRNA | barrnap e tRNAscan-SE |
| Tassonomia per contig | MMseqs2 + database tassonomico, oppure output Kraken2/Kaiju e tassonomia NCBI |
| Supporto di sequenza / 16S | BLAST+, database appropriato; Meryl opzionale per k-mer |

Tipo dei percorsi database: GTDB-Tk = directory dati; CheckM2 e GUNC = file `.dmnd`; MMseqs2 e BLAST = prefissi dei rispettivi database. Registrare release, checksum e provenienza. `fetch-databases` aggiorna la configurazione Snakemake: i comandi standalone non la caricano automaticamente.

Finché le risorse workflow non vengono incluse nella distribuzione, `setup-workflow` richiede un checkout accessibile. Una normale installazione del package non garantisce la presenza di `workflow/`.

## Contratti degli input

- Il BAM deve riferirsi esattamente al FASTA usato nel controllo: stessi nomi e stesse coordinate. Per il ramo A seguire il requisito HiFi del progetto, conservando i tag SA e le informazioni supplementari necessarie.
- Usare identificativi dei contig globalmente unici fra i bin e nomi dei MAG univoci. I join di molti report usano lo stem del file. Alcuni join ignorano maiuscole/minuscole, altri no: conservare la stessa grafia.
- I read ID nel BAM devono corrispondere a quelli del classificatore. Una read assente dal file tassonomico non è una read classificata; questa distinzione è attualmente difettosa nel ramo read-level.
- I nomi rinominati da `organize-mags --rename` differiscono dagli ID dei report originali. Eseguire i controlli prima della rinomina e conservare una mappa esplicita originale→output; il package non produce un manifest completo di questa trasformazione.
- Il lettore FASTA usa solo la prima parola dell'header. Una descrizione contenente annotazioni biologiche non viene conservata come parte del nome.
- Gzip è gestito dalle utility FASTA, ma non è garantito uniformemente in tutti i wrapper e join. La procedura sotto usa file `.fa` non compressi.

Gli output anvi’o 9 reali hanno questi schemi:

```text
clipping: contig  length  pos  relative_pos  cov  clipping  clipping_ratio
zero-cov: contig  length  range  range_size
```

Il parser corrente di Echidna non li interpreta entrambi correttamente: vedere F19/F20 nell'audit. Non adattare silenziosamente i file senza conservare l'originale e documentare la trasformazione. La sorgente primaria degli schemi è la [guida ufficiale anvi’o 9](https://anvio.org/help/9/programs/anvi-script-find-misassemblies/).

## Ramo A — pulizia delle misassemblies

Ordine delle dipendenze:

```text
assembly + BAM sul riferimento originale
  → anvi’o: clipping / zero-coverage
  → process-large-contigs: CheckM2 dei contig lunghi e assembly residuo
  → check-read-chimeras: report e finestre, se disponibili tassonomie delle reads
  → classify-clipping: annotazione delle posizioni, report-only
  → filter-misassemblies: splitting/rimozione usando i clipping
  → nuovo mapping sul FASTA modificato
  → depth aggiornata e binning esterno
```

Il wrapper si può ispezionare senza eseguire i tool:

```bash
echidna misassembly-wf \
  --assembly assembly.fa --bam mapping.hifi.bam \
  --output-dir misassembly_preview --dry-run
```

Questo mostra il piano, non verifica la correttezza biologica né l'esistenza degli output intermedi.

**Limiti attuali:** F01/F03 possono gestire male contig lunghi; F04 impedisce il percorso con tassonomia delle reads; F05 perde gli shift; F06/F07 distorcono le statistiche HiFi; F02/F18 non gestiscono correttamente i casi senza eventi. Prima di affidare l'intero ramo al wrapper occorre correggere e verificare questi casi.

`classify-clipping` non controlla automaticamente quali posizioni tagliare. Il filtro usa tutti i clipping in modalità `--safe`, oppure quelli che superano `--min-clip-coverage`. “Safe” è il nome della modalità: non è una garanzia contro tagli non necessari. Il filtro può rimuovere gli HQ puliti per lasciarli come output separato; tenere quindi un inventario dell'assembly residuo **e** dei contig estratti.

Dopo aver cambiato nomi o sequenze dei contig, rigenerare BAM e depth per il binning successivo. La depth del riferimento originale non descrive automaticamente il FASTA splittato.

## Ramo B — metriche e tassonomia dei MAG

Il wrapper esegue `checkm2 → organize1 → gtdbtk → gunc → [mmseqs → contigtax] → mimag → chimeras → organize2`. MMseqs è abilitato da `--mmseqs-db`. Per default porta avanti **HQ,MQ** secondo le soglie interne; i MAG della tier LQ non ricevono questi controlli downstream nel wrapper.

Per ispezionarlo:

```bash
echidna mag-quality-wf \
  --bins-dir bins --output-dir mag_preview \
  --genome-ext fa --tiers HQ,MQ \
  --mmseqs-db /data/mmseqs_gtdb/gtdb \
  --split-by-domain --dry-run
```

`--from` riparte da un passo e `--skip` lo omette. Non ricostruiscono gli output mancanti né verificano che quelli esistenti corrispondano agli input correnti. F10 rende rischioso riusare le cartelle dopo variazioni di soglie o qualità.

### Esempio standalone con percorsi distinti

Esempio di coordinamento da adattare ai propri ambienti e database già installati. Non è una pipeline certificata end-to-end: applicare i limiti dell'audit e controllare gli esiti mancanti.

```bash
set -euo pipefail
E_BINS=/percorso/assoluto/bins
E_OUT=/percorso/assoluto/run_nuova
E_GTDB=/data/gtdbtk_release
E_GUNC=/data/gunc/reference.dmnd
E_CHECKM2=/data/checkm2/reference.dmnd
E_MMSEQS=/data/mmseqs_gtdb/gtdb
E_THREADS=16
mkdir -p "$E_OUT"

conda run -n checkm2 echidna run-checkm2 \
  --genome-dir "$E_BINS" --extension fa --threads "$E_THREADS" \
  --database "$E_CHECKM2" -o "$E_OUT/checkm2"

conda run -n gtdbtk env GTDBTK_DATA_PATH="$E_GTDB" \
  echidna run-gtdbtk --genome-dir "$E_BINS" --extension fa \
  --cpus "$E_THREADS" -o "$E_OUT/gtdbtk"

conda run -n gunc echidna run-gunc \
  --genome-dir "$E_BINS" --file-suffix .fa --threads "$E_THREADS" \
  --db-file "$E_GUNC" -o "$E_OUT/gunc"

conda run -n mimag echidna run-mimag-rna \
  --genome-dir "$E_BINS" --extension fa --threads "$E_THREADS" \
  --split-by-domain --gtdbtk-dir "$E_OUT/gtdbtk" --outseq \
  -o "$E_OUT/mimag"

conda run -n mmseqs2 echidna run-mmseqs-contigs \
  --bins-dir "$E_BINS" --db "$E_MMSEQS" --extensions fa \
  --threads "$E_THREADS" -o "$E_OUT/contigs.mmseqs.tsv"

echidna classify-contigs --bins-dir "$E_BINS" --extensions fa \
  --mmseqs-output "$E_OUT/contigs.mmseqs.tsv" \
  --rank order -o "$E_OUT/contigs.classified.tsv"

echidna check-chimeras --bins-dir "$E_BINS" --extensions fa \
  --gunc-dir "$E_OUT/gunc" \
  --contig-taxonomy "$E_OUT/contigs.classified.tsv.bins.tsv" \
  --checkm2-report "$E_OUT/checkm2/quality_report.tsv" \
  --gtdbtk-dir "$E_OUT/gtdbtk" -o "$E_OUT/chimera_report.tsv"

echidna final-report --genomes-dir "$E_BINS" --extensions fa \
  --checkm2-report "$E_OUT/checkm2/quality_report.tsv" \
  --gtdbtk-dir "$E_OUT/gtdbtk" --mimag-rna-dir "$E_OUT/mimag" \
  --chimera-report "$E_OUT/chimera_report.tsv" \
  --contig-taxonomy "$E_OUT/contigs.classified.tsv.bins.tsv" \
  -o "$E_OUT/final_report.tsv"
```

Scegliere una directory nuova per ogni run modificata. Verificare che ci sia esattamente una riga per ciascun MAG atteso in ogni report richiesto, controllare le basi tassonomicamente risolte e i warning, quindi definire gli ID da accettare o rivedere. `organize-mags` organizza per qualità ma non applica questa policy chimera.

Per organizzare i MAG mantenendo separata questa decisione:

```bash
echidna organize-mags \
  --checkm2-root "$E_OUT/checkm2" --genomes-dir "$E_BINS" \
  --gtdbtk-dir "$E_OUT/gtdbtk" --mimag-rna-dir "$E_OUT/mimag" \
  --tiers HQ,MQ,LQ --extensions fa --rename --copy \
  -o "$E_OUT/organized_by_quality"
```

Le copie risultanti sono organizzate per qualità; i report vanno ancora usati per la selezione scientifica. F09 lascia HQ alcuni esiti RNA mancanti; `final-report.tier` usa le soglie CheckM2 predefinite e non la tier declassata o personalizzata.

Il controllo SSU è separato: i FASTA di rRNA dello standalone si trovano in `mimag/rrna_fasta/`. Dare quel percorso a `check-ssu-taxonomy` insieme a un database BLAST SSU con lineage compatibili e ai risultati GTDB-Tk. Non basta assumere che qualsiasi FASTA SSU contenga la tassonomia nell'header. Passare poi `--ssu-report` a `final-report`.

## Interpretazione delle tier e dei verdetti

### Tier interne del codice corrente

Il primo criterio soddisfatto prevale:

| Tier | Completezza | Contaminazione |
|---|---|---|
| HQ | ≥90% | ≤5% |
| MQ | ≥70% | ≤10% |
| LQ | ≥50% | ≤10% |
| Fail | Nessuna delle precedenti | |

Queste sono soglie interne. Non chiamarle automaticamente MIMAG: cambiano bande e confini dello standard. L'RNA è un requisito aggiuntivo per HQ MIMAG e un risultato mancante non dimostra un pass. [Standard originale MIMAG](https://swifter.embl.de/publication/pdf/28787424.pdf).

`process-large-contigs` usa invece soglie predefinite 50/10 per la propria estrazione detta “HQ”: l'etichetta non è equivalente alla tier HQ della tabella.

### Chimere per MAG

| Segnali disponibili | Risultato implementato |
|---|---|
| GUNC negativo al controllo e frazione estranea ≥0.05 | `Chimeric` |
| Uno dei segnali indica anomalia | `Check` |
| Nessuna anomalia nei segnali disponibili | `Clean` |
| Nessun segnale disponibile | `Unknown` |

GUNC negativo al controllo significa `pass.GUNC=False`. La frazione estranea è pesata per lunghezza **tra i contig risolti**, non fra tutte le basi del bin. GC, contaminazione CheckM2 e warning GTDB sono colonne informative e non cambiano il verdetto.

`Clean` non prova assenza di chimere; `Unknown` non è Clean. F11 può oggi far sembrare disponibile il segnale MMseqs anche quando nulla è risolto. I due segnali non costituiscono una probabilità calibrata.

`check-read-chimeras` usa invece `High/Medium/Low/Clean/Insufficient` in `chimera_risk`. Non confondere questo schema con `chimera_verdict` del controllo MAG. Il filtro assembly legge ancora soltanto il vecchio campo di rischio.

### Tabella finale

`final-report` unisce gli ID originali dei MAG alle evidenze disponibili. `NA`, `N/A`, campi vuoti e Unknown indicano forme di assenza di informazione che devono essere esaminate; il formato non è ancora uniforme. Una colonna `concerns` vuota significa soltanto che nessuna regola di allarme implementata è scattata.

Concerns correnti: `contam>10`, `chimeric`, `chimera_check`, `contig_chimera`, `ssu_disagree`, `mimag_rna_incomplete`. Non sono un criterio formale di accettazione. Il report non sostituisce `merged_quality.tsv` e non conserva le soglie personalizzate.

## Snakemake e HPC

`workflow/Snakefile` orchestra CheckM2, GTDB-Tk, GUNC, RNA opzionale e organizzazione/report chimera. **Non contiene** il ramo assembly, MMseqs/classify-contigs, SSU o final-report. Il report chimera del DAG è quindi basato su GUNC soltanto, salvo estensioni del workflow.

Con `use_named_envs: false` usare `--use-conda` e preparare i file env affinché installino il checkout corretto. Con `use_named_envs: true` configurare tutti gli ambienti richiesti e usare il conda dell'utente. Il kingdom MIMAG nel DAG è fisso: un insieme misto batteri/archei richiede gestione appropriata.

CPU, RAM e tempo dipendono da tool, dataset e release dei database. Le indicazioni HPC della guida sono esempi, non misure effettuate durante questa revisione. Su filesystem condivisi tutti i nodi devono vedere input, checkout, database e ambienti.

## Catalogo rapido dei comandi

| Area | Comandi |
|---|---|
| Assembly | `run-anvi-misassemblies`, `process-large-contigs`, `classify-clipping`, `filter-misassemblies`, `remove-hq-contigs`, `check-zero-coverage` |
| Chimere e tassonomia | `check-read-chimeras`, `check-chimeras`, `run-gunc`, `run-mmseqs-contigs`, `classify-contigs`, `check-ssu-taxonomy` |
| Qualità e report | `run-checkm2`, `run-gtdbtk`, `run-mimag-rna`, `organize-mags`, `final-report` |
| Workflow e setup | `misassembly-wf`, `mag-quality-wf`, `setup-workflow`, `fetch-databases`, `check-envs` |

Ogni comando espone `--help`. La reference descrittiva va riallineata in alcuni punti al parser corrente.

## Troubleshooting essenziale

| Sintomo | Controllo |
|---|---|
| Nessun `chimera_read_windows.tsv` | Specificare il FASTA nel comando read-level; verificare lunghezze e reads per finestra. Il wrapper corrente lo omette. |
| Filtro coverage perde tutti i clipping | Controllare schema anvi’o e F19, prima di ridurre la soglia. |
| Nessuna regione zero-cov da un file non vuoto | Controllare F20: l'intervallo nativo è nella colonna `range`. |
| HQ senza evidenze RNA | Confrontare set degli ID e fallimenti dei tool; vedere F09. |
| MAG presente in due tier | Output riutilizzato: F10. Rigenerare in una directory nuova e confrontare i manifest. |
| Tier differente nel report finale | `tier` è ricalcolata da CheckM2; verificare `merged_quality.tsv` e gate RNA. |
| Nessun genoma MIMAG con input `.fna` | Nel wrapper manca la propagazione dell'estensione: F16. |
| Errore con `--sample-regex` contenente `\d` | F15: problema del writer di configurazione, non della regex biologica. |
| File workflow assenti dopo installazione | Usare il checkout esplicito fino alla correzione del packaging, F21. |

## Test, riproducibilità e manutenzione

Nel checkout, usando un interprete con pytest:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m pytest -p no:cacheprovider
```

Per gli esperimenti dell'audit, usando il percorso della cartella originale del progetto:

```bash
python reproduce_findings.py /percorso/echidna
```

Conservare input originali, hash del codice, comandi effettivi, versioni dei tool, release dei DB, log ed elenco dei genomi prima/dopo ogni selezione. Le riproduzioni dell'audit usano mock e fixture sintetiche: non eseguono i tool pesanti.

La manutenzione prioritaria comprende formati nativi versionati, preservazione degli ID, stati mancanti espliciti, test dei collegamenti producer→consumer e due smoke test reali separati. L'elenco completo con trigger e rimedi è in [AUDIT.md](AUDIT.md).
