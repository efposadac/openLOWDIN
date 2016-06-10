/*
file: Iterators.h
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "Iterators.h"
#include <algorithm>

IntraIntsIt::IntraIntsIt(const int &n1, const int &n2, const int &n3,
                         const int &n4, const int &fi1, const int &fi2,
                         const int &fi3, const int &fi4)

    : ni(n1), nj(n2), nk(n3), nl(n4), fii(fi1), fij(fi2), fik(fi3), fil(fi4) {

  done = false;

  // printf("Quartet info: n: (%2d %2d | %2d %2d), fi: (%2d %2d | %2d %2d)\n",
  // ni, nj, nk, nl, fii, fij, fik, fil);

  iimax = ni - 1;
  if (fii == fij && fik == fil && fii == fik) {
    kkmax = 0;
    llmax = 0;
    jjmax = 0;
  } else if (fii == fik && fij == fil) {
    kkmax = 0;
    llmax = 0;
    jjmax = nj - 1;
  } else {
    kkmax = nk - 1;
    jjmax = (fii == fij) ? 0 : nj - 1;
    llmax = (fik == fil) ? 0 : nl - 1;
  }

  ii = 0;
  jj = 0;
  kk = 0;
  ll = 0;
}

void IntraIntsIt::first() {
  current.i = 0 + fii;
  current.j = 0 + fij;
  current.k = 0 + fik;
  current.l = 0 + fil;
  current.index = 0;
  if (fii == fij && fik == fil && fii == fik) { // (aa|aa) case
  } else if (fii == fik && fij == fil) {
    if (current.i < current.j) {
      std::swap(current.i, current.j);
      std::swap(current.k, current.l);
    }
    if (current.i < current.k) {
      std::swap(current.i, current.k);
      std::swap(current.j, current.l);
    }
  } else {
    if (current.i < current.j) {
      std::swap(current.i, current.j);
    }
    if (current.k < current.l) {
      std::swap(current.k, current.l);
    }
    if ((current.i < current.k) ||
        (current.i == current.k && current.j < current.l)) {
      std::swap(current.i, current.k);
      std::swap(current.j, current.l);
    }
  }

  // printf("first: (%2d %2d | %2d %2d) \n", current.i, current.j, current.k,
  //        current.l);
}

void IntraIntsIt::next() {
  if (fii == fij && fik == fil && fii == fik) {
    ++ll;
    if (ll > llmax) {
      ++kk;
      ll = 0;
      if (kk > kkmax) {
        kk = 0;
        ++jj;
        if (jj > jjmax) {
          jj = 0;
          ++ii;
          if (ii > iimax) {
            done = true;
          }
          jjmax = ii;
        }
        kkmax = ii;
      }
      llmax = (kk == ii) ? jj : kk;
    }
    current.i = ii + fii;
    current.j = jj + fij;
    current.k = kk + fik;
    current.l = ll + fil;
    current.index = ll + nl * (kk + nk * (jj + nj * ii));

    // printf("case 1: (%2d %2d | %2d %2d) \n", current.i, current.j, current.k,
    //        current.l);

  } else if (fii == fik && fij == fil) { //(ab|ab)
    ++ll;
    if (ll > llmax) {
      ++kk;
      ll = 0;
      if (kk > kkmax) {
        kk = 0;
        ++jj;
        if (jj > jjmax) {
          jj = 0;
          ++ii;
          if (ii > iimax) {
            done = true;
          }
        }
        kkmax = ii;
      }
      llmax = (kk == ii) ? jj : nl - 1;
    }
    current.i = ii + fii;
    current.j = jj + fij;
    current.k = kk + fik;
    current.l = ll + fil;
    current.index = ll + nl * (kk + nk * (jj + nj * ii));
    if (current.i < current.j) {
      std::swap(current.i, current.j);
      std::swap(current.k, current.l);
    }
    if (current.i < current.k) {
      std::swap(current.i, current.k);
      std::swap(current.j, current.l);
    }
    // printf("case 2: (%2d %2d | %2d %2d) \n", current.i, current.j, current.k,
    //        current.l);
  } else {
    ++ll;
    if (ll > llmax) {
      ++kk;
      ll = 0;
      if (kk > kkmax) {
        kk = 0;
        ++jj;
        if (jj > jjmax) {
          jj = 0;
          ++ii;
          if (ii > iimax) {
            done = true;
          }
          jjmax = (fii == fij) ? ii : nj - 1;
        }
      }
      llmax = (fik == fil) ? kk : nl - 1;
    }
    current.i = ii + fii;
    current.j = jj + fij;
    current.k = kk + fik;
    current.l = ll + fil;
    current.index = ll + nl * (kk + nk * (jj + nj * ii));
    if (current.i < current.j) {
      std::swap(current.i, current.j);
    }
    if (current.k < current.l) {
      std::swap(current.k, current.l);
    }
    if ((current.i < current.k) ||
        (current.i == current.k && current.j < current.l)) {
      std::swap(current.i, current.k);
      std::swap(current.j, current.l);
    }
    // printf("case 3: (%2d %2d | %2d %2d) \n", current.i, current.j, current.k,
    //        current.l);
  }
}
