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
- [x] commenter le code
- [x] reformater le code
- [x] à remettre: la vérification avec la MKL, 
- [ ] je dois trouver le min de la moyenne des deux transpositions dans find_minima ?
- [ ] dans boxplot, actualiser le bon numéro de thread
- [ ] revérifier qu'il n'y a pas de fuites mémoires avec wang+MKL
- [ ] dans wang, les asserts ne fonctionne pas avec scantrans^2 et les pre-transpose car celles-ci sont triées par ligne et non par ligne puis par colonne. Ainsi, la sortie colIdx de scantrans qui est triée par ligne puis par colonne n'est pas dans le même ordre que dans le fichier .mtx
- [ ] essayer avec O2 et O3 pour chaque algo, O3 defavorables ? Gustavon OK, std::sort OK, tbb::sort OK, MKL , scan, merge
- [ ] comparer le temps de scan et de merge avec leur version et la mienne
- [ ] comparer la MKL séquentielle et la parallèle avec 1 threads
- [ ] citation/licence SpaSM ou refactor ?
- [ ] comparer matrice par matrice quel algo est le meilleur
- [ ] taille de W : n, nnz, max(n,m)+1 ?
- [ ] checker les todo
- [x] décrire le fonctionnement de tbb:parallel_sort
- [x] décrire le fonctionnement de std::sort
- [ ] comparer la compression ?

GCC, (MKL iomp O3 AVX2), ScanTrans omp O3 AVX2, MergeTrans omp O3 AVX2

## Remarques

keep_better_parallel ne regarde que les threads de la première transpose, semble OK

std::sort est plus rapide que classical sur pre-transpose[6,7,8,9,10,11,12]. Son écart de durée est d'environ 5 ms (< 10ms). La variabilité inter-algo est alors plus petite que la variablité intra-algo.

Lorsque classical transpose et transpose_tr sont relativement proches, std transpose est plus lent que std transpose_tr

MKL iomp:
 s'arrete au nombre de coeurs physique (44)

MKL iomp 2PV:
 est idéal entre 8 et 24 threads (voire souvent entre 8 et 15, se confirme sur le graphique grobal). Pour les matrices RSA: vers 8 pour les premières puis vers 12 pour les dernières. Pour les matrices de Wang, la courbe atteint parfois un minimum plus tard vers 20-24.
 grande différence selon la forme du rectangle
 plus la matrice est (petite ou dense (première de pre-transpose)) plus le passage à l'echelle est mauvais

MergeTrans:
	est idéal entre 12 et 32 (souvent vers 16, plutôt entre 20 et 32 sur le graphique global)
	peu de différence selon la forme du rectangle sur les matrices de Wang, très grandes différence sur les 10 premières de pre-transpose (pas le même algo ?)

ScanTrans:
	ressemble à MergeTrans, mais le passage à l'échelle semble meilleur car le minimum est atteint un peu plus loin (la décroissance est plus applati jusqu'à 40 et la croissance est moins semble forte)
	sur les 10 premières de pre-transpose les deux transpositions utilisent le même algo
	encore moins de différences entre les deux transpositions que MergeTrans

transient est spéciale ?

TBB:
	sur les 10 premières de pre-transpose, il se passe quelque chose à partir de 80 et grand écart entre les deux transpositions
	courbe très plate, minimum atteint vers 16 (ou 44 sur les dernières)

ScanTrans est globalement le meilleur algo (overall duration et speed up)

## Licence

Code from pandas: BSD
Code from Wang: GNU LGPL 2.1
Code from SpaSM: GNU GPL 3
