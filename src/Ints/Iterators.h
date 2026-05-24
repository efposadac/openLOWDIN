/*
file: Iterators.h
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#ifndef ITERATORS_H
#define ITERATORS_H

class IntraIntsIt {
 private:
  struct Integral {
    int i;
    int j;
    int k;
    int l;
    unsigned int index;
  };

  Integral current;

  bool done;

  int ii, iimax, jj, jjmax, kk, kkmax, ll, llmax;
  const int &ni, &nj, &nk, &nl, &fii, &fij, &fik, &fil;

 public:
  IntraIntsIt(const int &n1, const int &n2, const int &n3, const int &n4,
              const int &fi1, const int &fi2, const int &fi3, const int &fi4);

  void first();
  void next();
  bool is_done() { return done; }

  int i() const { return current.i; }
  int j() const { return current.j; }
  int k() const { return current.k; }
  int l() const { return current.l; }
  int index() const { return current.index; }
};

#endif