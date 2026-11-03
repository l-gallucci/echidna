# Wiki di Echidna

La pagina mantenuta nel repository è [wiki/Home.md](../wiki/Home.md). La wiki GitHub, dopo l'attivazione, sarà disponibile su [Echidna Wiki](https://github.com/l-gallucci/echidna/wiki).

Per attivarla: nel repository GitHub aprire **Settings → General → Features → Wikis**, abilitare Wikis e creare/salvare una prima pagina Home dalla scheda Wiki. GitHub usa un repository Git separato per la wiki; il push del codice principale non la pubblica automaticamente.

Dopo aver creato la prima pagina, eseguire personalmente:

```bash
git clone git@github.com:l-gallucci/echidna.wiki.git /Users/lgallucc/Desktop/LocalBioinf/echidna-wiki
cp /Users/lgallucc/Desktop/LocalBioinf/echidna/wiki/Home.md /Users/lgallucc/Desktop/LocalBioinf/echidna-wiki/Home.md
cd /Users/lgallucc/Desktop/LocalBioinf/echidna-wiki
git diff -- Home.md
git add Home.md
git commit -m "docs: publish Echidna wiki"
git push
```

Se il clone wiki esiste già, aggiornarlo con `git pull --ff-only` invece di clonarlo di nuovo. Verificare prima eventuali modifiche locali e incorporare gli aggiornamenti fatti dal browser prima di sostituire Home.md.

Il progetto non pubblica automaticamente un sito GitHub Pages: la wiki integrata è un servizio separato. Fonti ufficiali: [abilitare la wiki](https://docs.github.com/en/communities/documenting-your-project-with-wikis/disabling-wikis), [aggiungere o modificare pagine](https://docs.github.com/en/communities/documenting-your-project-with-wikis/adding-or-editing-wiki-pages).
