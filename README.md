# Transposition de matrices creuses

- [Transposition de matrices creuses](#transposition-de-matrices-creuses)
  - [Installation](#installation)
  - [Utilisation](#utilisation)
    - [Matrices](#matrices)
    - [Intel Math Kernel Library (MKL)](#intel-math-kernel-library-mkl)
    - [Intel Threading Building Block (TBB)](#intel-threading-building-block-tbb)
  - [TODO](#todo)
  - [Licence](#licence)

## Installation

Testé avec:
Python 3.8.3
matplotlib 3.2.1
numpy 1.18.5
pandas 1.0.3
Intel Math Kernel Library (MKL) 2020.0.166
Intel Threading Building Block 2020.0.166
ICC TODO
GCC 6.3.0 20170516 (Debian 6.3.0-18+deb9u1)

## Utilisation

### Matrices

Le jeu de matrices utilisé par *Wang et. al* est disponible sur la machine ppti-gpu-1 dans le dossier ou dans le fichier `matrices.txt` :

`/Infos/lmd/2019/master/ue/MU4IN903-2020fev`

Taille : 5Go.

Un autre jeu de matrices (le mien) est disponible dans le dossier :

`/Infos/lmd/2019/master/ue/MU4IN903-2020fev/RSA`

Taille : 14Go. Ces matrices sont plus rectangulaires.

### Intel Math Kernel Library (MKL)

Pour utiliser la MKL sur ppti-gpu-1 :

`source /usr/intel/mkl/bin/mklvars.sh intel64 [ilp64]`

La dernière option (qui est facultative) active l'usage d'entiers 64 bits. On peut essayer pour voir si ça résout nos problèmes...

Il faut l'activer dans le Makefile (`USE_TBB=yes`), puis faire `make clean` et recompiler.

Après, il y a plusieurs options (séquentiel, parallèle avec OpenMP, etc). Choisir la bonne dans le Makefile. J'ai utilisé le [Intel® Math Kernel Library Link Line Advisor](https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html) pour les déterminer...

### Intel Threading Building Block (TBB)

Pour utiliser TBB sur ppti-gpu-1 :

`source /usr/intel/tbb/bin/tbbvars.sh intel64`

Ensuite, `make clean` puis recompiler.

## TODO

- [x] afficher la version de MKL
- [x] afficher la version de TBB
- [x] afficher le nombre de threads utilisés avec MKL
- [x] afficher le nombre de threads utilisés avec TBB
- [x] trouver le remplaçant de task_scheduler_init
- [x] créer un nouveau makefile
- [x] changer les noms des fichiers
- [x] changer la hierarchie des fichiers
- [x] initialiser time et Total
- [x] créer fonction triplet_transpose
- [x] fixer le nombre de threads de omp for avec num_threads
- [x] donner num_threads à MKL
- [x] gérer les benchmarks
- [x] créer un csv
- [x] fonctionne avec gcc/g++ ? OUI
- [x] fonctionne avec gcc/g++ et tbb ? OUI, mais fuites mémoire non gérables dues à TBB
- [x] fonctionne avec gcc/g++ et MKL ? OUI
- [x] fonctionne avec icc/icpc? OUI
- [x] fonctionne avec icc/icpc et tbb ? OUI
- [x] fonctionne avec icpc et mkl ? OUI
- [x] mettre dans un sous-dossier, makefile qui créer deux executables.
- [x] créer regle make benchmarks
- [x] modifications des entetes fallacieux des matrices RSA.
- [x] utiliser valgrind pour corriger les fuites mémoire et les erreurs de driver_wang : sptrans -> segfault, vec -> segfault?
- [x] mettre un assert dans la fonction load en supposant que l'entête est faux à la place de spasm_max
- [x] afficher dans le csv lire use_mkl,compiler_flags
- [x] mettre Have_mkl dans le name dans le csv
- [x] enlever le chemin d'accès des matrices dans le csv
- [x] prendre les algo parallel, moyenner les algo parallel, parcourir les matrices, parcourir les moyennes d'algo parallel (choisir un algo et parcourir le même sur les threads), tracer transpose en fonction de threads, idem pour transpose_tr, récupérer thread = argmax
- [x] parcourir les matrices, parcourir les algo, tracer les boites à moustaches
- [x] améliorer la ligne de compilation (les sauts de lignes dus à -D)
- [x] traiter le csv et afficher des graphes (boxplot pour la variance, pour parallèle : courbe de l'accélération), dataframe pandas
- [x] créer un sous-dossier dans charts automatiquement
- [x] pour la durée : regarder le min sur la moyenne des deux ?
- [x] sauvegarder les fichiers avec des _ au lieux des  
- [x] mettre une ligne avec plus de pointillés et des symboles plus petits
- [x] pour résumer : prendre l'addition des medianes par algo sur toutes les matrices
- [x] graphique accélération, calculer l'accélération avec un algo séquentiel pour chaque matrice : Gustavson -O2
- [x] refactor code
- [x] utiliser typedef u32 pour n et m et u32/u64 pour nz
- [x] finalize,spasm_add_entry, spasm_human_format mis dans .h et sans static
- [ ] essayer avec O2 et O3 pour chaque algo, O3 defavorables ? Gustavon OK, std::sort OK, tbb::sort
- [ ] commenter le code **driver et python**
- [ ] checker les todo
- [ ] COMPARER AVEC LE CODE SUR PPTI AVANT DE COMMIT (Makefile...)
- [ ] à remettre: la vérification avec la MKL, citation/licence SpaSM ou refactor ?
- [ ] la seule utilité de matrix_entry_t est pour finalize, utiliser spasm_triplet ?
- [ ] taille de W : n, nnz, max(n,m)+1 ?
- [x] reformater le code
- [x] décrire le fonctionnement de tbb:parallel_sort
- [x] décrire le fonctionnement de std::sort

GCC, MKL sequentielle, MKL iomp

keep_better_parallel ne regarde que les threads de la première transpose, semble OK

std::sort est plus rapide que classical sur pre-transpose[6,7,8,9,10,11,12]. Son écart de durée est d'environ 5 ms (< 10ms). La variabilité inter-algo est alors plus petite que la variablité intra-algo.

Lorsque classical transpose et transpose_tr sont relativement proches, std transpose est plus lent que std transpose_tr

## Licence

Code from Wang: GNU LGPL 2.1
Code from SpaSM: GNU GPL 3
