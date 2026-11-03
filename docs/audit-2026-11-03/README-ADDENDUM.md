> Historical snapshot of commit de465bcd. For the revalidation and current fixes see [FIXES.md](FIXES.md); current wiki source: [Home.md](../../wiki/Home.md).

# Integrazione proposta per il README

Bozza da adattare alla lingua inglese del README attuale. Non applicata automaticamente: descrive il comportamento del commit `de465bcd7dc47cba7be1e157f1af5e97f38bc0ba` e non presume che i difetti dell'audit siano già corretti.

## Da aggiungere dopo il Quickstart

Echidna comprende due rami: pulizia dell'assembly, documentata per PacBio HiFi, e qualità/tassonomia di MAG già ottenuti da un binner esterno. Il primo richiede un BAM delle stesse reads usate per generare l'assembly. Dopo splitting o rinomina dei contig è necessario rifare il mapping per calcolare la depth del nuovo riferimento.

Il core non richiede librerie Python aggiuntive. I wrapper richiedono tool bioinformatici e database installati separatamente. Gli ambienti conda forniti sono la strategia di installazione del progetto; Echidna deve essere disponibile negli ambienti che eseguono i suoi comandi.

## Da aggiungere prima della lista della documentazione

Le tier predefinite di `organize-mags` sono HQ ≥90%/≤5%, MQ ≥70%/≤10%, LQ ≥50%/≤10% per completezza/contaminazione. Sono tier interne, non una corrispondenza rigorosa allo standard MIMAG. Per valutare i requisiti RNA utilizzare `run-mimag-rna`; un esito mancante non equivale a un esito positivo. `process-large-contigs` usa un'altra selezione, chiamata HQ, con default 50%/10%.

`check-chimeras` produce `chimera_verdict`: Chimeric se GUNC e tassonomia per contig indicano entrambi un'anomalia, Check se la indica un solo segnale, Clean se nessuno dei segnali disponibili la indica, Unknown se non sono disponibili segnali. GC, contaminazione CheckM2 e warning GTDB sono informativi. La concordanza è una regola di screening, non una probabilità di correttezza del MAG.

`mag-quality-wf` seleziona HQ,MQ per default. Il report chimera non esclude automaticamente MAG dalle cartelle finali. Anche `classify-clipping` è report-only: la sua classificazione non viene consumata dal filtro assembly. Definire una policy esplicita di revisione/esclusione dei MAG segnalati.

`final-report` e `check-ssu-taxonomy` sono passaggi separati dai wrapper. Il DAG Snakemake copre un sottoinsieme del ramo MAG e non include MMseqs/classify-contigs. I nomi dei genomi devono rimanere coerenti durante i controlli; conservare una mappa quando si rinomina.

## Stato della verifica da rendere visibile

Alla revisione del 3 novembre 2026 passano 179 test del repository e gli help dei 22 comandi. La validazione end-to-end dei tool esterni non è stata eseguita. Sono stati riprodotti difetti di integrazione rilevanti, inclusi parser dei formati anvi’o 9, matching degli ID dei contig lunghi, gestione degli esiti mancanti e statistiche BAM HiFi. Consultare audit e wiki prima di usare i risultati per una selezione definitiva.

Questo paragrafo deve essere sostituito da uno stato aggiornato quando le correzioni e i test reali saranno disponibili, con riferimento al nuovo commit e alle versioni dei tool.

## Rimandi da aggiungere

- Wiki separata: `docs/audit-2026-11-03/WIKI.md`.
- Relazione di revisione: `docs/audit-2026-11-03/AUDIT.md`.
- Requisiti di runtime uniformati fra metadata, INSTALL e `check-envs`.
- Schemi input/output minimi, database/release e inventario degli output.
- Precondizioni di `--from/--skip`, directory nuove per run modificate e gestione dei risultati vuoti.
- Chiarezza sul checkout richiesto da `setup-workflow` finché le risorse non sono distribuite col package.

## Correzioni ai documenti collegati

Riallineare il numero di test, l'ordine anvi→process, i nomi dei campi `is_taxon_shift`/`chimera_verdict`, la lista dei concerns e i requisiti Python. Togliere i commenti dopo i backslash di continuazione. Separare il TSV MMseqs originale dal TSV di classificazione. Correggere l'esempio del binning per rigenerare mapping e depth dopo lo splitting. Evitare “strict MIMAG” e “full MIMAG” quando confini o evidenze mancanti non sono gestiti in modo coerente.
