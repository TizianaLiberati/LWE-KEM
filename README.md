# LWE-KEM
Master Thesis project

TODO:

1) Dal confronto con FrodoKEM (altro KEM basato su LWE presentato alla competizione NIST) questo risulta meno "evoluto" quindi bisognerebbe presentarlo non come un'innovazione dal punto di vista delle performance o della sicurezza ma come uno studio accademico con valutazione delle performance (? - da discutere)
2) Analizzare passaggio su GPU (il bottleneck sono le funzioni di hash)

DONE:

Generali:
1) implementa XOF (funzione che preso un seed genera in maniera deterministica una stringa di lunghezza variabile)
2) implementa SHA3-256 da sostituire con SHA256
3) cambiare il generatore random mt19937 con uno CSPRNG (rispetta le guide NIST) preso dalla libreria Botan
4) sostituite funzioni di hash con quelle presenti nella libreria OpenSSL (più veloci)
5) aggiunta parallelizzazione con openMP (generazione di N chiavi)

PKE:
- KeyGen: 
1) generare una stringa casuale z lunga 256 e concatenarla con sk (serve per l'implicit rejection)

- Encrypt:
1) calcolare i noise r, e1, e2 non più da stringhe casuali ma da seeds che verranno passate da Encaps


KEM:
- Encaps: 
1) generare un plaintext random lungo 256
2) calcolare:
   pk_h = Hash(A || t);
   seed = Hash( pk_h || m);
   XOF( seed ) (questa serve per calcolare le stringhe da cui vengono generati i noise in encrypt legandoli al messaggio così che nella decaps si possa tornare allo stesso ciphertext, non so quanto deve essere lungo l'output, dipende da come decidiamo di generare i noise)

- Decaps: 
1) come per Encaps vanno rigenerati i noise
2) implementare FO con implicit rejection
