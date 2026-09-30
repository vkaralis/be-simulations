# Προσομοιώσεις βιοϊσοδυναμίας σε Python

Καμπύλες πιθανότητας αποδοχής βιοϊσοδυναμίας (%) ως προς τον **πραγματικό GMR T/R**. Συγγραφέας: **Vangelis D. Karalis**.

Python simulations of bioequivalence studies and Monte Carlo power curves for 2x2, 2x3, 2x4 and 3x3 crossover designs.

**Tags:** Bioequivalence, Simulations, Monte Carlo, Power Curves, 2x2, 2x3, 2x4, 3x3, Crossover Designs.

**GitHub topics:** `bioequivalence`, `simulations`, `monte-carlo`, `power-curves`, `2x2`, `2x3`, `2x4`, `3x3`, `crossover-designs`, `python`.

## Εγκατάσταση

Python 3.10 ή νεότερη:

```bash
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -e .
```

## Εκτέλεση

Κάθε αρχείο ρωτά συνολικό αριθμό ατόμων, επαναλήψεις Monte Carlo **ανά GMR**, residual CV (%) και τιμή ή εύρος GMR:

```bash
python BE2x2.py
python BE2x3.py
python BE2x4.py
python BE3x3.py
```

| Αρχείο | Ακολουθίες | Προϋπόθεση N |
| --- | --- | --- |
| BE2x2.py | TR / RT | άρτιος, τουλάχιστον 4 |
| BE2x3.py | RTR / TRT | άρτιος, τουλάχιστον 4 |
| BE2x4.py | TRTR / RTRT | άρτιος, τουλάχιστον 4 |
| BE3x3.py | TRR / RTR / RRT | πολλαπλάσιο του 3, τουλάχιστον 6 |

Το `3x3` εδώ σημαίνει **3 ακολουθίες × 3 περιόδους με δύο σκευάσματα**, όπως στον αρχικό κώδικα. Δεν σημαίνει τρία διαφορετικά σκευάσματα. Ίσος αριθμός ατόμων ανά ακολουθία· μη έγκυρο N απορρίπτεται και δεν αυξάνεται αυτόματα.

Παράδειγμα χωρίς διαδραστικές ερωτήσεις:

```bash
python BE2x2.py --subjects 24 --trials 100000 --cv 30 --gmr 0.70:1.40:0.01 --seed 2026
python -m besim --design 3x3 --subjects 24 --trials 10000 --cv 30 --gmr 0.95 --seed 42 --show
```

`--cv 30` σημαίνει CV=30%, **όχι** 0.30%. Το εύρος γράφεται `αρχή:τέλος:βήμα` και το τέλος περιλαμβάνεται μόνο όταν ανήκει στο πλέγμα. Μία τιμή GMR δίνει ένα σημείο στο γράφημα.

Στον φάκελο `results/` δημιουργούνται:

- `power_<design>.png`: καμπύλη ισχύος και σημειακά 95% Wilson διαστήματα αβεβαιότητας Monte Carlo.
- `power_<design>.csv`: GMR, πλήθος αποδοχών, ποσοστό αποδοχής και όρια αβεβαιότητας.
- `power_<design>.json`: παράμετροι, seed, ακολουθίες, βαθμοί ελευθερίας και μοντέλο.

Τα ίδια ονόματα αντικαθίστανται σε επόμενη εκτέλεση. Χρησιμοποιήστε `--output results/run_02` για ξεχωριστό φάκελο. Το `--show` ανοίγει επίσης το γράφημα. Το seed παράγεται και αποθηκεύεται αν δεν δοθεί· για αναπαραγωγή χρησιμοποιήστε το ίδιο seed και τις ίδιες εκδόσεις εξαρτήσεων.

## Στατιστικό μοντέλο

Για ένα log-PK endpoint:

`log(PK_ij) = subject_i + period_j + log(GMR) * I(T) + error_ij`

Ανεξάρτητα κανονικά residuals με κοινή διακύμανση για T και R, `sigma² = log(1 + (CV%/100)²)`. Το μοντέλο περιλαμβάνει ξεχωριστή σταθερά ανά άτομο, άρα απορροφά και την ακολουθία. Η ενδοατομική αντίθεση θεραπείας και το residual MSE εκτιμώνται όπως σε μοντέλο subject + period + treatment. Για πλήρη δεδομένα με κοινό residual και τυχαία σταθερά ατόμου, αυτή η ενδοατομική ανάλυση εξαλείφει την επίδραση ατόμου.

Αποδοχή όταν ολόκληρο το 90% CI του T/R βρίσκεται στο [0.80, 1.25], ισοδύναμα TOST με α=0.05. Το CI με t κατανομή χρησιμοποιεί **εκτιμώμενη**, όχι γνωστή, residual διακύμανση.

Ο κοινός κινητήρας κατασκευάζει τον πίνακα X και υπολογίζει `df = observations - rank(X)` και `v = [(X'X)^(-1)]_treatment,treatment`. Στους ισορροπημένους σχεδιασμούς:

| Σχεδιασμός | Residual df | Var(εκτίμησης log GMR) / sigma² |
| --- | --- | --- |
| 2x2 | N−2 | 2/N |
| 2x3 | 2N−3 | 1.5/N |
| 2x4 | 3N−4 | 1/N |
| 3x3 | 2N−3 | 1.5/N |

Προσομοιώνονται απευθείας τα επαρκή στατιστικά κάθε μελέτης: η αντίθεση από `Normal(log(GMR), sigma²*v)` και ανεξάρτητα το MSE από `sigma²*ChiSquare(df)/df`. Αυτή είναι ακριβής προσομοίωση της κατανομής των αποτελεσμάτων του παραπάνω κανονικού μοντέλου, χωρίς την αργή παραγωγή και προσαρμογή όλων των PK εγγραφών. Δεν παράγονται ατομικά δεδομένα. Κοινοί τυχαίοι αριθμοί μεταξύ GMR μειώνουν τον θόρυβο της καμπύλης· οι εκτιμήσεις διαφορετικών σημείων συσχετίζονται και η ζώνη είναι σημειακή, όχι ταυτόχρονη.

Υποθέσεις: πλήρεις μελέτες, χωρίς dropout, carryover ή subject-by-treatment interaction. Υποστηρίζεται **κλασική average BE για ένα endpoint** με κοινό CV. Τα replicate designs δεν ενεργοποιούν αυτόματα scaled BE, ούτε διαφορετικό CV για T/R. Η πιθανότητα αφορά το συγκεκριμένο endpoint, όχι κοινή αποδοχή AUC και Cmax. Πρόκειται για εργαλείο προσομοίωσης υπό αυτές τις υποθέσεις.

Πηγές κριτηρίου και πλαισίου: [EMA guideline](https://www.ema.europa.eu/system/files/documents/scientific-guideline/wc500070039_en.pdf), [EMA: Investigation of bioequivalence](https://www.ema.europa.eu/en/investigation-bioequivalence-scientific-guideline). Η EMA επισημαίνει ότι το ICH M13A έχει αντικαταστήσει σχετικά τμήματα για non-replicate designs από 25/01/2025. Οι προσομοιώσεις εδώ δεν υλοποιούν κάθε ειδική κανονιστική μέθοδο.

## Διορθώσεις και έλεγχοι

Δείτε [CHANGELOG.md](CHANGELOG.md) για τα σφάλματα των MATLAB αρχείων και τις αλλαγές μεθόδου.

```bash
python -m unittest discover -s tests -v
```

Έλεγχοι: ανεξάρτητοι τύποι df/διακύμανσης, προσαρμογή πραγματικών προσομοιωμένων log-PK δεδομένων, σύγκριση ισχύος με αριθμητική ολοκλήρωση της ακριβούς κατανομής, σφάλμα τύπου Ι στα όρια, επαναληψιμότητα και μη έγκυρες είσοδοι. Το GitHub Actions εκτελεί τους ελέγχους.

## Συγγραφέας

Vangelis D. Karalis. Η αρχική εφαρμογή αναπτύχθηκε τον Μάρτιο 2011.


## Δομή repository

```text
BEsimALL/
├── BE2x2.py
├── BE2x3.py
├── BE2x4.py
├── BE3x3.py
├── besim/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py
│   └── core.py
├── tests/
│   ├── test_core.py
│   └── test_cli.py
├── .github/workflows/tests.yml
├── .gitignore
├── pyproject.toml
├── requirements.txt
├── README.md
└── CHANGELOG.md
```
