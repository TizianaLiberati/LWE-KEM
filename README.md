# LWE-KEM
Master Thesis project

TODO:

Generali:
1) valuta se esistono funzioni "migliori" per fare i sample dalle distribuzioni (da questo dipende la lunghezza della stringa in output della XOF) 

PKE:

KEM:

 
Avevo provato a vedere delle ottimizzazioni a livello di aritmetica come suggeriva Matteo (RNS ad esempio) però molte sono efficienti con un modulo molto più grande o quando il numero di prodotti/addizioni è maggiore del numero di conversioni da fare (un numero va convertito in un vettore particolare ed è in generale costoso). 
C'era anche l'idea di fare una cifratura (encrypt) a blocchi, non bit a bit, ma avremmo comunque un ciphertext molto grande e non so se ci sarebbe una reale efficienza


DONE:

Generali:
1) implementa XOF (funzione che preso un seed genera in maniera deterministica una stringa di lunghezza variabile)
2) implementa SHA3-256 da sostituire con SHA256
3) cambiare il generatore random mt19937 con uno CSPRNG (rispetta le guide NIST) preso dalla libreria Botan

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
