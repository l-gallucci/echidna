> Historical snapshot of commit de465bcd. For the revalidation and current fixes see [FIXES.md](FIXES.md); current wiki source: [Home.md](../../wiki/Home.md).

# Echidna — revisione tecnica e documentale

Data: 3 novembre 2026. Versione: **0.0.4**. Checkout: `/Users/lgallucc/Desktop/LocalBioinf/echidna`.
Commit esaminato: `de465bcd7dc47cba7be1e157f1af5e97f38bc0ba`.

## Esito

**Non è possibile considerare la pipeline interamente affidabile nello stato esaminato.** La suite esistente passa, ma non intercetta diversi difetti nei collegamenti tra comandi, nei formati dei tool e nella gestione degli esiti mancanti. Alcuni casi possono estrarre contig problematici, perdere contig dagli output di lavoro, rimuovere sequenze sbagliate o produrre verdetti rassicuranti senza evidenza sufficiente.

La separazione fra logica Python e wrapper esterni rende il progetto leggibile e verificabile. Sono utili anche i piani `--dry-run`, le chiamate subprocess generalmente costruite come liste di argomenti e i test dei parser. Questi aspetti non sostituiscono però una validazione delle due pipeline su dati reali e formati nativi.

Questa revisione **non modifica il codice applicativo, il README o i risultati scientifici esistenti**. La pagina wiki è un documento separato. Le priorità sotto sono valutazioni di revisione: P1 = rischio sui risultati o blocco di un percorso documentato; P2 = robustezza, compatibilità o semantica/documentazione da correggere.

## Verifiche eseguite e limiti

| Verifica | Esito |
|---|---|
| Suite del repository, Python 3.13.13 di Miniforge | **179 passed**, 0.16 s |
| `--version`, help principale, help dei 22 sottocomandi | **24/24** terminano con codice 0, Python 3.14.7 |
| Parsing sintattico dei file Python tracciati | **61/61**, inclusi 42 file del package; grammatica Python 3.8 |
| `bash -n examples/complete_workflow.sh` | Passa; verifica solo la sintassi |
| Esperimenti mirati stdlib su dati sintetici e mock dichiarati | **21/21 comportamenti confermati**: F01–F20 e D01 |
| Creazione sdist in copia temporanea del checkout | Passa; **workflow non incluso** (F21) |
| Coverage quantitativa | Non misurata: `pytest-cov` non disponibile nell'interprete di test |
| Tool esterni, database, Snakemake/SLURM su dati reali | **Non eseguiti end-to-end** |

Comando della suite:

```bash
PYTHONDONTWRITEBYTECODE=1 /Users/lgallucc/software/miniforge3/bin/python \
  -m pytest -p no:cacheprovider
```

Il `python3` predefinito era Python 3.14.7 e non aveva pytest; è stato usato l'interprete Miniforge già presente, senza installare pacchetti. Il parsing con grammatica 3.8 non certifica il funzionamento su Python 3.8. L'analisi riguarda il checkout corrente tracciato, non le vecchie copie MaMISA in `.claude/worktrees/`. Nessuna percentuale di copertura o accuratezza biologica è deducibile dai test passati.

Gli esperimenti sono in [reproduce_findings.py](reproduce_findings.py); gli esiti in [reproduction-results.json](reproduction-results.json). Un esperimento che passa **riproduce il difetto attuale**, non certifica il risultato come corretto. È deliberatamente separato dai test applicativi.

## Risultati con impatto maggiore

### F01 — P1 — Identità dei contig persa tra CheckM2 e misassemblies

Riferimenti: `echidna/commands/process_large_contigs.py:77`, `:160`, `:369`.

I file individuali sono denominati con `safe_filename(nome)` e hash finale. CheckM2 restituisce quel nome, mentre i clipping sono indicizzati dal nome FASTA originale. Il confronto `contig_name in misasm_contigs` non converte le identità. Nell'esperimento `contig_7`, con clipping a 200 e qualità 95/1, diventa `contig_7_93105b`, viene marcato `EXTRACT_HQ` e sparisce dall'assembly destinato allo splitting. La sequenza rimane in `HQ_extracted`, ma viene sottratta al controllo previsto.

Correzione proposta: persistere una mappa biunivoca ID originale ↔ file CheckM2 e usarla per decisioni, report e ricostruzione. Test necessario: input FASTA + clipping nativo + report CheckM2 con hash, verifica dell'output finale.

### F03 — P1 — Contig senza riga CheckM2 persi dagli output successivi

Riferimenti: `process_large_contigs.py:369`, `:392` e `create_updated_assembly()`.

Le decisioni percorrono solo le righe presenti nel report. Un contig lungo estratto ma assente da un report parziale, precedente o non corrispondente non finisce né in `HQ_extracted` né in `assembly_for_filtering.fa`. I file intermedi conservano la sequenza, ma il successivo binning la perde. La fixture con due contig e una sola riga CheckM2 conferma il caso.

Correzione: controllare che tutti i contig estratti abbiano un esito univoco; fermarsi con errore esplicito oppure mantenere i mancanti come qualità sconosciuta. Verificare anche righe duplicate ed estranee.

### F04 — P1 — `misassembly-wf` con tassonomia richiede un file che non produce

Riferimenti: `echidna/commands/misassembly_wf.py:61`, `:76`; `check_read_chimeras.py:450` e fase di scrittura.

Il piano non inoltra `--assembly` a `check-read-chimeras`. Senza lunghezze i profili non memorizzano posizioni e non eseguono l'analisi a finestre. Il file `chimera_read_windows.tsv` non viene scritto; il passo `classify` lo richiede comunque e si ferma per file mancante. La stessa omissione appare negli esempi della guida.

Correzione: inoltrare l'assembly originale e produrre sempre una tabella con intestazione anche quando non ci sono finestre utilizzabili. Testare il contratto di output, non solo la presenza degli argomenti nel piano.

### F06 — P1 — Le statistiche locali BAM contano partenze, non copertura

Riferimenti: `echidna/utils/bam_stats.py:256`, `:412`.

Gli accumulatori sono aggiornati solo quando la **partenza** dell'allineamento è entro ±`window` dalla posizione di clipping. Una read lunga che attraversa la posizione o termina lì, ma parte più lontano, è ignorata. La fixture SAM `9000M100S`, POS=1000, clipping=10000 perde sia la read sia il tag SA del partner.

Inoltre `local_depth` è un numero di read con partenza nella finestra, confrontato con una profondità media per base di `samtools coverage`: non sono misure equivalenti. Il rapporto varia con lunghezza delle reads e finestra anche a copertura uniforme. Le etichette repeat/deletion possono quindi risultare distorte.

Correzione: calcolare overlap e basi allineate tramite CIGAR, adottare coordinate uniche e filtri coerenti, normalizzare la profondità locale per le basi della regione. Validare il rapporto su BAM sintetici con copertura nota e reads di lunghezze diverse.

### F07 — P1 — Le reads HiFi non appaiate diventano tutte “discordanti”

Riferimenti: `bam_stats.py:190`; `echidna/commands/classify_clipping.py:239`.

Una read primaria senza flag proper-pair incrementa `n_discordant` anche se non è paired-end. Dieci reads non appaiate, senza altri segnali anomali e con rapporto di profondità 1, producono `discordant_fraction=1` e `chimera_candidate`. Questo contraddice l'uso HiFi del ramo A.

Correzione: statistiche di coppia solo per allineamenti con flag paired; per HiFi riportare NA e non contribuire ai punteggi. I test attuali degli accumulatori usano flag=0, ma non verificano questa proprietà.

### F08 — P1 — Errore BLAST trasformato in evidenza di sequenza senza supporto

Riferimenti: `echidna/utils/zero_coverage.py:299`; `echidna/commands/check_zero_coverage.py:162`.

`blast_regions_vs_reads()` restituisce `{}` sia quando non ci sono hit sia quando BLAST termina con errore. Il chiamante usa `hits.get(id, 0)` e assegna `unsupported`. Un errore tecnico simulato con exit 2 produce quindi un verdetto biologico negativo e può indurre a eliminare sequenza valida.

Correzione: distinguere esecuzione fallita da esecuzione riuscita senza hit; propagare errore e log del tool. Lo stesso principio va applicato ai flussi samtools: alcuni percorsi ignorano il return code o degradano a risultati vuoti.

### F09 — P1 — Il controllo MIMAG richiesto lascia HQ gli esiti mancanti

Riferimenti: `echidna/commands/filter_checkm2.py:211`, `:241`; `run_mimag_rna.py:_process_dir()`.

Un MAG HQ senza riga rRNA/tRNA rimane HQ. Se una directory non contiene summary validi, il gate non viene applicato. I fallimenti per singolo genoma in barrnap/tRNAscan vengono omessi dal summary mentre altri genomi possono completare la run. Il comportamento è intenzionalmente codificato da `test_gate_missing_record_keeps_hq`, ma **non supporta la promessa di un gate MIMAG completo**.

Correzione: mantenere uno stato “non valutato/errore”, verificare la corrispondenza completa dei MAG e impedire una certificazione HQ MIMAG senza tutte le evidenze richieste. Separare `tier_checkm2` da `mimag_status`.

### F10 — P1 — Rilanciare l'organizzazione conserva selezioni obsolete

Riferimento: `filter_checkm2.py:327` e ciclo di copia/link.

La seconda esecuzione crea o sostituisce i file nella tier nuova senza rimuovere quelli precedenti. Un MAG inizialmente HQ e poi MQ compare in entrambe le cartelle. Un MAG che diventa Fail può restare nella vecchia selezione. Il rischio riguarda anche `mag-quality-wf --from`.

Correzione: manifest degli output posseduti dalla run e pubblicazione atomica di una nuova selezione, oppure rifiutare output già popolati senza opzione esplicita. Non cancellare indiscriminatamente file dell'utente.

### F11 — P1 — Nessun contig classificato può diventare `Clean`

Riferimenti: `echidna/utils/contig_taxonomy.py:_summarize_bins()`; `echidna/commands/check_chimeras.py:87`, `:232`.

Un bin con zero contig risolti genera `foreign_frac=0`. Il consumatore carica solo la frazione, ignora `n_contigs_resolved` e considera il segnale tassonomico disponibile e negativo. Senza GUNC l'esperimento produce **Clean con zero contig classificati**. Anche un valore malformato della frazione viene convertito a 0.

Correzione: preservare disponibilità e quantità dell'evidenza; usare NA/Unknown sotto una soglia minima di basi risolte. Valori malformati devono generare errore, non assenza di contaminazione. Documentare che il denominatore della frazione corrente comprende solo le sequenze risolte.

### F12 — P1 — Le reads assenti dal classificatore risultano classificate

Riferimenti: `echidna/commands/check_read_chimeras.py:112`; `echidna/utils/read_taxonomy.py:314`, `:317`.

Il collegamento BAM assegna taxid=-1 alle reads non presenti nel file tassonomico. Le statistiche escludono solo taxid=0: -1 diventa un taxon reale, entra nelle finestre e incrementa il numero di reads classificate. Cento reads assenti producono classificazione 100%, taxon dominante `taxid:-1` e un possibile esito Clean.

Correzione: ammettere come classificati solo taxid positivi; mantenere separati unclassified e unmatched; richiedere una copertura minima del join e del classificatore. Non confrontare inoltre taxid a ranghi diversi come se fossero organismi diversi: manca un roll-up uniforme del ramo read-level.

### F13 — P1 — `remove-hq-contigs --hq-dir` può rimuovere contig estranei ai MAG

Riferimenti: `echidna/commands/remove_hq_contigs.py:29`; `filter_misassemblies.py:138`.

La directory HQ viene interpretata attraverso numeri nei **nomi dei file**, senza leggere i contig presenti nei FASTA. Un `bin.42.fa` contenente `contig_7` e `contig_9` fa rimuovere `contig_42` dall'assembly, lasciando invece entrambi i contig del MAG HQ. L'esempio completo usa proprio una cartella di MAG selezionati, quindi rende il problema concretamente raggiungibile. L'estrazione di soli ID numerici può anche collidere fra campioni/assembler diversi e con i suffissi hash.

Correzione: leggere gli identificativi FASTA esatti dai MAG e confrontarli con quelli dell'assembly; rilevare duplicati e ID non corrispondenti. La lista esplicita deve dichiarare se contiene identificativi interi o nomi completi.

### F19 — P1 — Colonna di clipping incompatibile con l'output anvi’o 9

Riferimento: `echidna/utils/misassembly.py:40–80`.

Il parser presume `contig_name, sample_name, position, clipping_coverage`. Lo schema ufficiale è `contig, length, pos, relative_pos, cov, clipping, clipping_ratio`. Il campo letto come intero è quindi una frazione; il fallback imposta 1. Nella fixture ufficiale 120 reads clippate diventano 1 e `--min-clip-coverage 50` elimina l'evento.

Correzione: leggere le colonne per nome e riconoscere esplicitamente le versioni del formato. La sorgente primaria e gli esempi ufficiali confermano lo schema: [guida anvi’o 9](https://anvio.org/help/9/programs/anvi-script-find-misassemblies/), [sorgente tag v9](https://github.com/merenlab/anvio/blob/v9/anvio/cli/find_misassemblies.py).

### F20 — P1 — Gli intervalli zero-coverage reali vengono tutti scartati

Riferimento: `echidna/utils/zero_coverage.py:55–99`.

Il parser si aspetta start/end numerici separati, mentre anvi’o 9 scrive `range` come `100500-101000`, seguito dalla dimensione dell'intervallo. `int(range)` fallisce, la riga viene saltata e il comando può concludere che non ci siano regioni da verificare. Il produttore usa intervalli zero-based half-open; il parser dichiara una conversione da coordinate 1-based.

Correzione: supportare il formato nativo `range`, mantenere coordinate corrette e non trasformare un file interamente non interpretabile in un risultato negativo valido. Fonti: [documentazione ufficiale](https://anvio.org/help/9/programs/anvi-script-find-misassemblies/), [implementazione v9](https://github.com/merenlab/anvio/blob/v9/anvio/cli/find_misassemblies.py).

## Altre criticità riprodotte

| ID | Priorità | Evidenza e intervento |
|---|---|---|
| F02 | P2 | `process_large_contigs.py:329`: nessun contig >300 kbp termina normalmente senza `assembly_for_filtering.fa`; il workflow lo richiede comunque. Scrivere sempre gli output contrattuali, anche nel caso vuoto. |
| F05 | P2 | `check_read_chimeras.py:216` scrive `is_taxon_shift`; `classify_clipping.py:116` legge `is_shift`. Le finestre di shift sono ignorate anche quando presenti. Usare schema condiviso e test producer→consumer. |
| F14 | P2 | `utils/final_report.py:183` ricalcola le tier CheckM2 con soglie predefinite: un MAG declassato a MQ dal gate RNA appare HQ nel report, e le soglie personalizzate sono perse. La reference dice che è una tier CheckM2, ma il nome generico `tier` facilita errori: rinominare la colonna o importare la selezione effettiva. |
| F15 | P2 | `setup_workflow.py:36–39`: `re.sub` interpreta la regex campione come replacement; `^(S\d+)_` genera `bad escape \d`. Servono replacement callable e escaping YAML corretto; verificare anche virgolette e backslash. |
| F16 | P2 | `mag_quality_wf.py:110–120`: `--genome-ext fna` non viene inoltrato a MIMAG, che cerca `.fa`. L'estensione non arriva neppure a `check-chimeras`, rilevante per estensioni fuori dai suoi default. Propagare il parametro a tutti i passi. |
| F17 | P2 | `filter_misassemblies.py:81`: è pubblicizzato input da `check-chimeras`, ma si legge solo `chimera_risk`; l'output 0.0.4 usa `chimera_verdict`. Un `Chimeric` viene ignorato. Occorre anche una mappa bin→contig, non una semplice traduzione di etichette. |
| F18 | P2 | `filter_misassemblies.py:289` e `classify_clipping.py:420`: un file clipping valido senza eventi causa exit 1. Distinguere “nessun evento” da file assente o schema non valido; gestire la run pulita. |
| F21 | P2 | `MANIFEST.in` / `setup.py`: la sdist costruita da una copia pulita contiene il package, ma nessun file `workflow/`. `setup-workflow` presume invece `<package>/../workflow`. L'uso con checkout editable è diverso da una normale installazione. Distribuire le risorse e trovarle con `importlib.resources`, oppure dichiarare esplicitamente il vincolo al checkout. |

## D01 — Comportamento progettuale verificato: screening ≠ esclusione

`mag_quality_wf.py:124` costruisce `organize2` senza il report chimera. Anche il DAG Snakemake organizza i MAG e calcola il report come rami separati. Quindi un MAG `Chimeric` può essere copiato nelle cartelle finali HQ/MQ se supera CheckM2/RNA. `classify-clipping` produce inoltre un report che `filter-misassemblies` **non legge**: la classificazione delle posizioni non guida automaticamente i tagli.

Non è necessariamente un bug: mantenere report-only è una scelta valida. Va però dichiarato nel README, nei diagrammi e nella wiki; “MAG finali filtrati su tutti i segnali” sarebbe una promessa falsa. Serve una policy esplicita di esclusione/revisione e una lista tracciabile degli ID selezionati. `mag-quality-wf` non comprende neppure `check-ssu-taxonomy` e `final-report`; si invocano separatamente.

## Fallacie e limiti metodologici

1. **“HQ” non equivale sempre a MIMAG HQ.** `process-large-contigs` chiama HQ una selezione che parte da completezza 50 e contaminazione 10. `organize-mags` usa HQ ≥90/≤5, MQ ≥70/≤10, LQ ≥50/≤10. È una classificazione interna e non coincide con lo standard. MIMAG usa HQ >90 e <5 più RNA, MQ ≥50 e <10, LQ <50 e <10. Anche la ricetta “strict MIMAG” della guida mantiene confronti inclusivi e non risolve i limiti di confine. [Bowers et al., tabella 1 e testo](https://swifter.embl.de/publication/pdf/28787424.pdf).
2. **Concordanza di due classificatori non prova indipendenza statistica.** GUNC e MMseqs possono condividere riferimenti e bias. Il consenso è una regola euristica implementata, non una probabilità calibrata di chimera. I test verificano la regola, non sensibilità/specificità. `Clean` con un solo segnale disponibile significa soltanto “nessuna anomalia in quel segnale”.
3. **Tassonomia eterogenea non dimostra da sola un join chimerico.** Reads ambiguamente assegnate, cross-mapping, regioni conservate e livelli tassonomici diversi richiedono interpretazione. Nel ramo read-level i taxid non sono riportati a un rango comune. La frase nel codice “the only reliable signal” è troppo assoluta.
4. **Lunghezza non dimostra circolarità.** Le euristiche usano >100 kbp nel GC e ≥200 kbp in un altro comando. Inoltre il lettore FASTA conserva solo la prima parola dell'header: annotazioni `circular=true` nella descrizione vengono perse. Il candidato è una priorità di ispezione, non evidenza di replicone completo.
5. **Un hit BLAST non certifica supporto per tutta una regione.** Il controllo zero-coverage verifica presenza di hit, senza percentuale minima di copertura della regione/identità. Un match locale non giustifica l'intera sequenza. Analogamente `check-ssu-taxonomy` sceglie il best bitscore; segnala bassa identità ma non la esclude dal confronto e non impone copertura allineata minima.
6. **Rimuovere i suffissi GTDB può cancellare differenze biologiche.** `normalize_taxon()` fa coincidere `Firmicutes_A` e `Firmicutes_B`. Il confronto SSU è una concordanza normalizzata, non identità rigorosa fra cladi GTDB.
7. **Cambiare l'assembly cambia il riferimento del BAM.** Dopo splitting/rimozione occorre riallineare prima di calcolare la depth per il nuovo binning. L'esempio della guida calcola invece depth dal BAM originale e lo associa all'assembly rinominato: contig e coordinate possono non corrispondere.

## Lacune del README e documentazione da riallineare

Il README è un buon indice dei comandi, ma mancano informazioni essenziali per interpretare una run. Alcune sono già presenti nei documenti collegati: vanno riassunte nel README senza duplicare tutta la reference.

| Tema | Situazione attuale | Aggiunta/correzione richiesta |
|---|---|---|
| Stato validazione | Nessuna distinzione test unitari vs pipeline reali | Dichiarare 179 test in questa revisione, limiti end-to-end e casi noti. La guida dice ancora 96. |
| Python | `INSTALL.md` ≥3.7, `setup.py`/docs ≥3.8, `check-envs` ≥3.9 | Unificare requisito e matrice CI; dichiarare gli interpreti effettivamente verificati. |
| Due rami | Nel README non è esplicito il limite HiFi del ramo A | Input, ordine corretto anvi→process, condizioni d'uso; binning/remapping sono esterni. |
| “Dependency-free” | Vero per il package, non per la pipeline | Separare dipendenze Python core, tool esterni, database e orchestratore. Gli env sono una strategia fornita, non una necessità logica universale. |
| Qualità | Non riassume le soglie né la selezione `HQ,MQ` del wrapper | Mostrare soglie interne, LQ escluso di default dal wrapper e differenza rispetto a MIMAG. |
| Gate RNA | Promessa di HQ completo senza politica sugli esiti mancanti | Distinguere pass, fail, not evaluated; spiegare il comportamento F09 fino alla correzione. |
| Chimere | Descrizione breve, nessuna policy finale | Definire Chimeric/Check/Clean/Unknown e dire che i MAG segnalati non vengono esclusi automaticamente. |
| Input e identità | Mancano schemi minimi e vincoli globali | FASTA/BAM sullo stesso riferimento, read IDs esatti, nomi dei bin e dei contig unici, mapping originale→rinominato. |
| File output | Nessuna lista dei risultati chiave | Tabelle, cartelle `Selected`, conservazione di HQ estratti, risultati vuoti, reinserimento finale o inventario separato. |
| Database | Link alla guida, ma nessun registro di release | Registrare release, checksum e compatibilità dei DB; `fetch-databases` aggiorna config Snakemake, non tutti i wrapper standalone. |
| Ripresa | `--from/--skip` senza precondizioni | Output intermedi necessari, directory nuova per run modificata, rischio risultati obsoleti e limiti di `--skip-existing`. |
| Distribuzione | `pip install .` e helpers presentati insieme | Spiegare checkout richiesto dagli helper finché F21 non è corretto. |
| Piattaforme | “macOS out of the box” non validato per tutti i tool | Separare portabilità core da disponibilità dei tool e risorse HPC. |

Altri punti documentali concreti:

- `docs/pipeline.md` e `docs/workflows.md` mostrano in alcuni schemi `process-large-contigs` prima di anvi’o, benché consumi i clipping di anvi’o.
- Alcuni comandi multilinea hanno un commento dopo `\`: la continuazione di shell non funziona perché il backslash non è l'ultimo carattere della riga.
- La guida riutilizza lo stesso percorso per il TSV raw MMseqs e per quello classificato: sovrascrive l'evidenza originaria. Usare due nomi distinti.
- La reference sostiene che `--assembly` serva se le lunghezze non sono nel BAM; il ramo read-level non le ricava dal BAM, quindi il FASTA serve comunque per le finestre.
- `final-report` elenca ancora concerns `chimera_High/chimera_Medium`; il codice corrente emette `chimeric/chimera_check`.
- `INSTALL.md` conserva tool/versioni e dipendenze associate a vecchi nomi; `filter-misassemblies` non esegue anvi’o e `organize-mags` non esegue CheckM2.
- Il DAG Snakemake non contiene MMseqs/classify-contigs, SSU o final-report: non equivale al wrapper MAG con tutti i segnali opzionali. Nei set misti usa un kingdom MIMAG fisso.
- Non sono presenti workflow CI tracciati. L'esempio shell contiene il binning come placeholder e non costituisce uno smoke test completo.
- `docs/command-reference.md` non elenca tutti i flag, per esempio `process-large-contigs --checkm2-env` e le soglie self-BLAST; `--help` è la fonte più vicina al parser.

## Robustezza da verificare dopo i problemi principali

Lettura statica, senza misurazione di incidenza su dataset reali: collisioni di basename e canonicalizzazione case-insensitive; contig con stesso nome in bin diversi; gzip ed estensioni composte non trattati uniformemente; coordinate clipping fuori limite o duplicate; input FASTQ troncati; directory riutilizzate da MMseqs; quote dei percorsi nelle regole shell Snakemake; risorse condivise sul cluster; timeout e return code dei processi; log stderr soppressi; summary GTDB riusati senza controllo sul set di genomi. `check-envs` verifica soprattutto presenza/versioni e alcuni vincoli GUNC, non la compatibilità completa delle API dei tool.

Diversi percorsi caricano sequenze o tassonomie intere in memoria nonostante l'uso di lettori streaming. La profilazione a finestre riesamina tutte le posizioni per ogni finestra. Servono benchmark prima di promettere scalabilità su metagenomi grandi. Gli argomenti di risorse e le soglie numeriche hanno validatori disponibili ma spesso non li utilizzano.

## Ordine di intervento suggerito

1. Formati nativi e identità: F19/F20, F01/F03/F13; preservare tutte le sequenze e i mapping.
2. Affidabilità degli esiti: F08/F09/F11/F12; stato mancante distinto da negativo.
3. Ramo HiFi: F04/F05/F06/F07 e casi senza eventi F02/F18.
4. Selezione ripetibile e report coerenti: F10/F14/F16/F17; policy esplicita D01.
5. Packaging/config, README e wiki; CI su Python supportati e test producer→consumer usando fixture native versionate.
6. Due smoke test reali separati, uno per ramo, con genomi/reads noti e risultati attesi; poi benchmark di accuratezza e risorse su dati rappresentativi. Eseguire anche ripresa, output vuoti e fallimento dei tool.

## Materiale consegnato

- [WIKI.md](WIKI.md): pagina autonoma in italiano con uso, contratti e limiti della revisione corrente.
- [README-ADDENDUM.md](README-ADDENDUM.md): testo proposto per integrare il README senza riscriverlo.
- [reproduce_findings.py](reproduce_findings.py), [reproduction-results.json](reproduction-results.json), [static-checks.json](static-checks.json), [packaging-check.json](packaging-check.json): evidenze ripetibili e risultati delle verifiche.

I riferimenti `file:linea` nel testo sono relativi alla radice del checkout al commit sopra. Le righe possono cambiare dopo correzioni. Le verifiche esterne sono limitate alle fonti primarie citate per formati anvi’o e standard MIMAG; gli altri rilievi derivano dal codice locale e dagli esperimenti allegati.
