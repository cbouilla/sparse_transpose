# Algorithmes de tri

- [Algorithmes de tri](#algorithmes-de-tri)
  - [Liens utiles](#liens-utiles)
  - [Algorithmes en cours d'étude](#algorithmes-en-cours-détude)
  - [Algorithmes par comparaison](#algorithmes-par-comparaison)
    - [Tri rapide (*Quicksort*)](#tri-rapide-quicksort)
      - [Caractéristique](#caractéristique)
      - [Complexité en temps](#complexité-en-temps)
      - [Complexité en espace](#complexité-en-espace)
      - [Optimisations](#optimisations)
    - [*Sledgesort*](#sledgesort)
      - [Caractéristiques](#caractéristiques)
    - [*Introsort* (*introspective sort*)](#introsort-introspective-sort)
      - [Caractéristiques](#caractéristiques-1)
      - [Complexité en temps](#complexité-en-temps-1)
    - [*Quick radix sort*](#quick-radix-sort)
    - [*BlockQuicksort*](#blockquicksort)
    - [Tri par tas *heapsort*](#tri-par-tas-heapsort)
      - [Caractéristiques](#caractéristiques-2)
      - [Complexité en temps](#complexité-en-temps-2)
      - [Complexité en espace](#complexité-en-espace-1)
      - [Optimisations](#optimisations-1)
    - [*Smoothsort*](#smoothsort)
      - [Complexité en temps](#complexité-en-temps-3)
      - [Complexité en espace](#complexité-en-espace-2)
    - [Tri fusion (*merge sort*)](#tri-fusion-merge-sort)
      - [Caractéristiques](#caractéristiques-3)
      - [Complexité en temps](#complexité-en-temps-4)
      - [Complexité en espace](#complexité-en-espace-3)
      - [Optimisations](#optimisations-2)
    - [*Block merge sort*](#block-merge-sort)
      - [Caractéristiques](#caractéristiques-4)
      - [Complexité en temps](#complexité-en-temps-5)
      - [Complexité en espace](#complexité-en-espace-4)
    - [*Timsort*](#timsort)
      - [Caractérisitques](#caractérisitques)
      - [Complexité en temps](#complexité-en-temps-6)
      - [Complexité en espace](#complexité-en-espace-5)
    - [*Batcher odd-even mergesort*](#batcher-odd-even-mergesort)
  - [Algorithmes sans comparaison](#algorithmes-sans-comparaison)
    - [*American flag sort*](#american-flag-sort)
      - [Caractéristiques](#caractéristiques-5)
    - [Tri par base (*Radix sort*)](#tri-par-base-radix-sort)
      - [Complexité en temps](#complexité-en-temps-7)
      - [Complexité en espace](#complexité-en-espace-6)
    - [*Burstsort*](#burstsort)
    - [Tri comptage (*compting sort*)](#tri-comptage-compting-sort)
      - [Complexité en temps](#complexité-en-temps-8)
      - [Complexité en espace](#complexité-en-espace-7)

## Liens utiles

https://pub.phyks.me/sdz/sdz/ri-rapide-ameliorations.html
https://openclassrooms.com/forum/sujet/tri-decouvrez-timsort-en-c-96826
https://en.cppreference.com/w/cpp/algorithm/sort

## Algorithmes en cours d'étude

- Intel TBB parallel_sort : std::sort parallel ? + std::sort sequentiel avec seuil à 500
- std::sort : introsort (quicksort avec median-of-3 pivot + heap sort avec seuil à 2 log_2(n)?) + insertion sort avec seuil à 16 (Knuth ?).
- séparer les nnz en p blocs, trier les p blocs (l'algo dépend de la taille), [fusionner les p blocs](https://en.wikipedia.org/wiki/K-way_merge_algorithm) et utiliser le fait que l'on sache approximativement où chercher les entrées

## Algorithmes par comparaison

[wikipedia en](https://en.wikipedia.org/wiki/Comparison_sort)

### Tri rapide (*Quicksort*)

[wikipedia en](https://en.wikipedia.org/wiki/Quicksort)
[wikipedia fr](https://fr.wikipedia.org/wiki/Tri_rapide)

#### Caractéristique

- en place
- non stable
- ne tire pas avantage de l'état de l'entrée
- parallèlisable

#### Complexité en temps

- Pire cas $O(n^{2})$
- Moyenne $O(n\log n)$
- Meilleur cas $O(n\log n)$

#### Complexité en espace

- Pire cas $O(n)$
- Moyenne $O(\log n)$

#### Optimisations

- choix du pivot : prendre une approximation de la mediane des éléments
- [double pivot](https://en.wikipedia.org/wiki/Quicksort#cite_note-:0-10)
- partitionnement
- récursivité terminale
- multi-key quicksort/three-way radix sort

### *Sledgesort*

[wikipedia en](https://en.wikipedia.org/wiki/Quicksort)

#### Caractéristiques

- utiliser le *quicksort* et s'arrêter à des petites sous-listes
- utiliser le tri par insertion soit sur le tableau entier (cause pleins de fautes de caches sur les grands tableaux) soit directement sur les petites sous-listes

### *Introsort* (*introspective sort*)

[wikipedia en](https://en.wikipedia.org/wiki/Quicksort)
[wikipedia en](https://en.wikipedia.org/wiki/Introsort)
[wikipedia fr](https://fr.wikipedia.org/wiki/Introsort)

#### Caractéristiques

- en place
- non stable
- utiliser le *quicksort*
- utiliser le tri par tas ou *smoothsort* (complexité en $O(n\log n)$ dans le pire cas) lorsque la profondeur de la récursion dépasse un certain seuil pour limiter la complexité dans le pire cas
- utiliser le insertion sort lorsque le nombre d'éléments est petit

#### Complexité en temps

- Pire cas $O(n\log n)$
- Moyenne $O(n\log n)$

### *Quick radix sort*

[wikipedia en](https://en.wikipedia.org/wiki/Quicksort)

### *BlockQuicksort*

[wikipedia en](https://en.wikipedia.org/wiki/Quicksort)

### Tri par tas *heapsort*

[wikipedia en](https://en.wikipedia.org/wiki/Heapsort)
[wikipedia fr](https://fr.wikipedia.org/wiki/Tri_par_tas)

#### Caractéristiques

- en place
- non stable
- parallèlisable ?

#### Complexité en temps

- Pire cas $O(n\log n)$
- Moyenne $O(n\log n)$
- Meilleur cas $O(n\log n)$ (distinct keys) or $O(n)$ (equal keys)

#### Complexité en espace

- Pire cas $O(n)$ total $O(1)$ auxiliary

#### Optimisations

- construction de Floyd
- construction en profondeur dès que possible pour les grandes grandes données [](https://en.wikipedia.org/wiki/Heapsort#cite_note-Bojesen00-6) [](https://en.wikipedia.org/wiki/Heapsort#cite_note-7)
- bottom up si la comparaison est couteuse

### *Smoothsort*

#### Complexité en temps

- Pire cas $O(n\log n)$
- Moyenne $O(n\log n)$
- Meilleur cas $O(n)$

#### Complexité en espace

- Pire cas $O(n)$ total $O(1)$ auxiliary

### Tri fusion (*merge sort*)

[wikipedia en](https://en.wikipedia.org/wiki/Merge_sort)
[wikipedia fr](https://fr.wikipedia.org/wiki/Tri_fusion)

#### Caractéristiques

- convient au cache/bonne localité des données
- stable
- pas en place
- parallèlisable

#### Complexité en temps

- Pire cas $O(n\log n)$
- Moyenne $O(n\log n)$
- Meilleur cas $O(n\log n)$, $O(n)$ natural variant

#### Complexité en espace

- Pire cas $O(n)$ total avec $O(n)$ auxialiary

#### Optimisations

- en place
- *tiled merge sort* utiliser un autre tri lorsque les sous-listes tiennent dans le cache [](https://en.wikipedia.org/wiki/Merge_sort#CITEREFLaMarcaLadner1997)
- mémoire en $O(1)$ [](https://en.wikipedia.org/wiki/Merge_sort#CITEREFKatajainenPasanenTeuhola1996)

### *Block merge sort*

[wikipedia en](https://en.wikipedia.org/wiki/Block_sort)

#### Caractéristiques

- stable
- en place

#### Complexité en temps

- Pire cas $O(n\log n)$
- Moyenne $O(n\log n)$
- Meilleur cas $O(n)$

#### Complexité en espace

- Pire cas $O(1)$

### *Timsort*

[wikipedia fr](https://fr.wikipedia.org/wiki/Timsort)

#### Caractérisitques

- stable
- utilie le fait que la liste soit déjà triée

#### Complexité en temps

- Pire cas $O(n\log n)$
- Moyenne $O(n\log n)$
- Meilleur cas $O(n)$

#### Complexité en espace

- Pire cas $O(n)$

### *Batcher odd-even mergesort*

[wikipedia en](https://en.wikipedia.org/wiki/Batcher_odd%E2%80%93even_mergesort)

## Algorithmes sans comparaison

### *American flag sort*

[wikipedia en](https://en.wikipedia.org/wiki/American_flag_sort)

#### Caractéristiques

- mauvaise gestion du cache
- en place

### Tri par base (*Radix sort*)

[wikipedia en](https://en.wikipedia.org/wiki/Radix_sort)
[wikipedia fr](https://fr.wikipedia.org/wiki/Tri_par_base)

#### Complexité en temps

- Pire cas $O(w n)$, where w is the number of bits required to store each key.

#### Complexité en espace

- Pire cas $O(w n)$, where w is the number of bits required to store each key.

### *Burstsort*

[wikipedia en](https://en.wikipedia.org/wiki/Burstsort)

### Tri comptage (*compting sort*)

[wikipedia en](https://en.wikipedia.org/wiki/Counting_sort)
[wikipedia fr](https://fr.wikipedia.org/wiki/Tri_comptage)

#### Complexité en temps

- Pire cas $O(k + n)$, where k is the range of the non-negative key values.

#### Complexité en espace

- Pire cas $O(k + n)$, where k is the range of the non-negative key values.
