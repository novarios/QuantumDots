#include "CCfunctions.hpp"
#include "MATHfunctions.hpp"

//   Function to return Hash index for 2 indices
int Hash2(const int &p, const int &q, const int &size)
{
  return size*p + q;
}

//   Function to return Hash index for 3 indices
int Hash3(const int &p, const int &q, const int &r, const int &size)
{
  return size*size*p + q*size + r;
}

//   Function to search for index1 of p in vec1 (size num1) and index2 of q in vec2 (size num2)
//   Returns index of pq in matrix <p|q>. Index = index1*num2 + index2.
int Index11(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[i] == p){ ind1 = i; break; }
    else if(i == num1 - 1){
      for(int j = 0; j < num1; ++j){ std::cout << vec1[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index11 for " << p << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[i] == q){ ind2 = i; break; }
    else if(i == num2 - 1){
      for(int j = 0; j < num2; ++j){ std::cout << vec2[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index11 for " << q << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;
}

//   Function to search for index of pq in vec1 (size num1)
//   Returns index of pq in matrix <p|q>.
int Index2(const int *vec1, const int &num1, const int &p, const int &q)
{
  int ind1;
  for(int i = 0; i < num1; ++i){
    if(vec1[2*i] == p && vec1[2*i + 1] == q){ ind1 = i; break; }
    else if(i == num1 - 1){ return -1; }
  }
  return ind1;
}


//   Returns index of p in vector <p>.
int Index1(const int *vec1, const int &num1, const int &p)
{
  int ind1;
  for(int i = 0; i < num1; ++i){
    if(vec1[i] == p){ ind1 = i; break; }
    else if(i == num1 - 1){
      for(int j = 0; j < num1; ++j){ std::cout << vec1[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index1 for " << p << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1;
}

//   Function to search for index1 of pq in vec1 (size num1) and index2 of rs in vec2 (size num2)
//   Returns index of pqrs in matrix <pq|rs>. Index = index1*num2 + index2.
int Index22(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[2*i] == p && vec1[2*i + 1] == q){
      ind1 = i;
      break;
    }
    else if(i == num1 - 1){
      for(int j = 0; j < num1; ++j){ std::cout << vec1[2*j] << "," << vec1[2*j + 1] << " "; }
      std::cout << std::endl;
      std::cerr << "Index22 for " << p << " " << q << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[2*i] == r && vec2[2*i + 1] == s){
      ind2 = i;
      break;
    }
    else if(i == num2 - 1){
      for(int j = 0; j < num2; ++j){ std::cout << vec2[2*j] << "," << vec2[2*j + 1] << " "; }
      std::cout << std::endl;
      std::cerr << "Index22 for " << r << " " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

//   Function to search for index1 of p in vec1 (size num1) and index2 of qrs in vec2 (size num2)
//   Returns index of pqrs in matrix <p|qrs>. Index = index1*num2 + index2.
int Index13(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[i] == p){ ind1 = i; break; }
    else if(i == num1 - 1){ 
      for(int j = 0; j < num1; ++j){ std::cout << vec1[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index13 for " << p << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[3*i] == q && vec2[3*i + 1] == r && vec2[3*i + 2] == s){ ind2 = i; break; }
    else if(i == num2 - 1){ 
      for(int j = 0; j < num2; ++j){ std::cout << vec2[3*j] << "," << vec2[3*j + 1] << "," << vec2[3*j + 2] << " "; }
      std::cout << std::endl;
      std::cerr << "Index13 for " << q << " " << r << " " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

//   Function to search for index1 of pqr in vec1 (size num1) and index2 of s in vec2 (size num2)
//   Returns index of pqrs in matrix <pqr|s>. Index = index1*num2 + index2.
int Index31(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[3*i] == p && vec1[3*i + 1] == q && vec1[3*i + 2] == r){ ind1 = i; break; }
    else if(i == num1 - 1){ 
      for(int j = 0; j < num1; ++j){ std::cout << vec1[3*j] << "," << vec1[3*j + 1] << "," << vec1[3*j + 2] << " "; }
      std::cout << std::endl;
      std::cerr << "Index31 for " << p << " " << q << " " << r << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[i] == s){ ind2 = i; break; }
    else if(i == num2 - 1){ 
      for(int j = 0; j < num2; ++j){ std::cout << vec2[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index31 for " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

void plus(State &S, const State &S1, const State &S2){
  S.t = S1.t + S2.t;
  S.m = S1.m + S2.m;
  S.nx = S1.nx + S2.nx;
  S.ny = S1.ny + S2.ny;
  S.nz = S1.nz + S2.nz;
  S.ml = S1.ml + S2.ml;
  S.par = S1.par * S2.par;
  S.j = S1.j + S2.j;
}

void minus(State &S, const State &S1, const State &S2){
  S.t = S1.t - S2.t;
  S.m = S1.m - S2.m;
  S.nx = S1.nx - S2.nx;
  S.ny = S1.ny - S2.ny;
  S.nz = S1.nz - S2.nz;
  S.ml = S1.ml - S2.ml;
  S.par = S1.par * S2.par;
  S.j = S1.j + S2.j;
}

bool equal(const State &S1, const State &S2){
  return (S1.t == S2.t && S1.m == S2.m && S1.nx == S2.nx && S1.ny == S2.ny && S1.nz == S2.nz && S1.ml == S2.ml && S1.par == S2.par && S1.j == S2.j);
}

double Amplitudes::get_energy(const Input_Parameters &Parameters, const Channels &Chan, const Interactions &Ints)
{
  double energy = 0.0;
  int nhh, npp;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    
    #pragma omp parallel
    {
      #pragma omp for schedule(static) reduction(+:energy)
      for(int hh = 0; hh < Chan.nhh[chan]; ++hh){
	for(int pp = 0; pp < Chan.npp[chan]; ++pp){
	  energy += 0.25 * D1.T1[chan][hh * Chan.npp[chan] + pp] * Ints.D_ME1.V4[chan][pp * Chan.nhh[chan] + hh];
	}
      }
    }
  }
  
  if(Parameters.approx == "singles" || Parameters.approx == "triples"){

    #pragma omp parallel
    {
      #pragma omp for schedule(static) reduction(+:energy)
      for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
	for(int hp2 = 0; hp2 < Chan.nhp1[Chan.ind0]; ++hp2){
	  energy += 0.5 * S1.T1[hp1] * S1.T1[hp2] * Ints.D_ME1.V9[Chan.ind0][hp1 * Chan.nhp1[Chan.ind0] + hp2];
	}
      }
    }
  }

  return energy;
}

Amplitudes::Amplitudes(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1 = Doubles_1(Parameters, Space, Chan);
  }
  else if(Parameters.approx == "singles"){
    D1 = Doubles_1(Parameters, Space, Chan);
    S1 = Singles_1(Parameters, Space, Chan);
  }
  else if(Parameters.approx == "triples"){
    D1 = Doubles_1(Parameters, Space, Chan);
    S1 = Singles_1(Parameters, Space, Chan);
  }
}

void Amplitudes::delete_struct(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1.delete_struct(Chan);
  }
  else if(Parameters.approx == "singles"){
    D1.delete_struct(Chan);
    S1.delete_struct(Chan);
  }
  else if(Parameters.approx == "triples"){
    D1.delete_struct(Chan);
    S1.delete_struct(Chan);
  }
}

void Amplitudes::zero(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1.zero(Chan);
  }
  else if(Parameters.approx == "singles"){
    D1.zero(Chan);
    S1.zero(Chan);
  }
  else if(Parameters.approx == "triples"){
    D1.zero(Chan);
    S1.zero(Chan);
  }
}

Doubles_1::Doubles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  int length0, length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  Tmap = new int*[Chan.size1];
  Evec = new double*[Chan.size1];
  T1 = new double*[Chan.size1];
  T2 = new double*[Chan.size2];
  T3 = new double*[Chan.size2];
  T4 = new double*[Chan.size2]; 
  T5 = new double*[Chan.size2];
  T6 = new double*[Chan.size3];
  T7 = new double*[Chan.size3];
  T8 = new double*[Chan.size3];
  T9 = new double*[Chan.size3];
  S1 = new double*[Chan.size1];
  S2 = new double*[Chan.size3];
  S3 = new double*[Chan.size3];
  S4 = new double*[Chan.size3];
  S5 = new double*[Chan.size3];
  S6 = new double*[Chan.size2];
  S7 = new double*[Chan.size2];

  Q11 = new double*[Chan.size1];
  Q21 = new double*[Chan.size1];
  Qmap1 = new int*[Chan.size1];
  Qmap2 = new int*[Chan.size1];
  Q12 = new double*[Chan.size3];
  Q22 = new double*[Chan.size3];

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      Tmap[i] = new int[16 * length];
      Evec[i] = new double[length];
      T1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	for(int k = 0; k < 16; ++k){ Tmap[i][16*j + k] = -1; }
	Evec[i][j] = 0.0;
	T1[i][j] = 0.0;
      }
    }
    length = nhh * nhh;
    if(length != 0){
      S1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S1[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      Q11[i] = new double[length];
      Qmap1[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
	Qmap1[i][2*j] = -1;
	Qmap1[i][2*j + 1] = -1;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      Q21[i] = new double[length];
      Qmap2[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
	Qmap2[i][2*j] = -1;
	Qmap2[i][2*j + 1] = -1;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    length = nhpp * nh;
    if(length != 0){
      T6[i] = new double[length];
      T7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T6[i][j] = 0.0;
	T7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      T8[i] = new double[length];
      T9[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T8[i][j] = 0.0;
	T9[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      Q12[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(nhhp * nh != 0){
      Q22[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q22[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      S2[i] = new double[length];
      S3[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S2[i][j] = 0.0;
	S3[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      S4[i] = new double[length];
      S5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S4[i][j] = 0.0;
	S5[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      T2[i] = new double[length];
      T3[i] = new double[length];
      T4[i] = new double[length];
      T5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
	T4[i][j] = 0.0;
	T5[i][j] = 0.0;
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      S6[i] = new double[length];
      S7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S6[i][j] = 0.0;
	S7[i][j] = 0.0;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, pind1, pind2;
      int hh1, pp1;
      length = Chan.nhh[i] * Chan.npp[i];
      length0 = Chan.npp[i];
      #pragma omp for schedule(static)
      for(int hhpp = 0; hhpp < length; ++hhpp){
	pp1 = int(hhpp%length0);
	hh1 = int((hhpp - pp1)/length0);
	hind1 = Chan.hhvec[i][2*hh1];
	hind2 = Chan.hhvec[i][2*hh1 + 1];
	pind1 = Chan.ppvec[i][2*pp1];
	pind2 = Chan.ppvec[i][2*pp1 + 1];
	//T2 = ((ia)(jb)')
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind1, pind1, hind2, pind2);
	Tmap[i][16 * hhpp] = ind;
	Tmap[i][16 * hhpp + 1] = ind1;
	//T3 = ((jb)(ia)')
	minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind2, pind2, hind1, pind1);
	Tmap[i][16 * hhpp + 2] = ind;
	Tmap[i][16 * hhpp + 3] = ind1;
	//T4 = ((ib)(ja)')
	minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind1, pind2, hind2, pind1);
	Tmap[i][16 * hhpp + 4] = ind;
	Tmap[i][16 * hhpp + 5] = ind1;
	//T5 = ((ja)(ib)')
	minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind2, pind1, hind1, pind2);
	Tmap[i][16 * hhpp + 6] = ind;
	Tmap[i][16 * hhpp + 7] = ind1;
 	//T6 = ((jab)(i))
	ind = Chan.indvec[hind1];
	key1 = Chan.hpp_map[ind][Hash3(hind2, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind1];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hppvec[ind], Chan.hvec[ind], Chan.nhpp[ind], Chan.nh[ind], hind2, pind1, pind2, hind1);
	Tmap[i][16 * hhpp + 8] = ind;
	Tmap[i][16 * hhpp + 9] = ind1;
	//T7 = ((iab)(j))
	ind = Chan.indvec[hind2];
	key1 = Chan.hpp_map[ind][Hash3(hind1, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind2];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hppvec[ind], Chan.hvec[ind], Chan.nhpp[ind], Chan.nh[ind], hind1, pind1, pind2, hind2);
	Tmap[i][16 * hhpp + 10] = ind;
	Tmap[i][16 * hhpp + 11] = ind1;
	//T8 = ((ijb)(a))
	ind = Chan.indvec[pind1];
	key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind2, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hhpvec[ind], Chan.pvec[ind], Chan.nhhp[ind], Chan.np[ind], hind1, hind2, pind2, pind1);
	Tmap[i][16 * hhpp + 12] = ind;
	Tmap[i][16 * hhpp + 13] = ind1;
	//T9 = ((ija)(b))
	ind = Chan.indvec[pind2];
	key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind1, Space.indtot)];
	key2 = Chan.p_map[ind][pind2];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hhpvec[ind], Chan.pvec[ind], Chan.nhhp[ind], Chan.np[ind], hind1, hind2, pind1, pind2);
	Tmap[i][16 * hhpp + 14] = ind;
	Tmap[i][16 * hhpp + 15] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, hind3, pind1, pind2, pind3;
      int hp1, pp1, hh1;
      length = Chan.nhp[i] * Chan.npp[i];
      length0 = Chan.npp[i];
      #pragma omp for schedule(static)
      for(int hppp = 0; hppp < length; ++hppp){
	pp1 = int(hppp%length0);
	hp1 = int((hppp - pp1)/length0);
	hind1 = Chan.hpvec[i][2*hp1];
	pind1 = Chan.hpvec[i][2*hp1 + 1];
	pind2 = Chan.ppvec[i][2*pp1];
	pind3 = Chan.ppvec[i][2*pp1 + 1];
	ind = Chan.indvec[pind1];
	key1 = Chan.hpp_map[ind][Hash3(hind1, pind2, pind3, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hppvec[ind], Chan.pvec[ind], Chan.nhpp[ind], Chan.np[ind], hind1, pind2, pind3, pind1);
	Qmap1[i][2 * hppp] = ind;
	Qmap1[i][2 * hppp + 1] = ind1;
      }
      length = Chan.nhh[i] * Chan.nhp[i];
      length0 = Chan.nhp[i];
      #pragma omp for schedule(static)
      for(int hhhp = 0; hhhp < length; ++hhhp){
	hp1 = int(hhhp%length0);
	hh1 = int((hhhp - hp1)/length0);
	hind1 = Chan.hhvec[i][2*hh1];
	hind2 = Chan.hhvec[i][2*hh1 + 1];
	hind3 = Chan.hpvec[i][2*hp1];
	pind1 = Chan.hpvec[i][2*hp1 + 1];
	ind = Chan.indvec[hind3];
	key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind1, Space.indtot)];
	key2 = Chan.h_map[ind][hind3];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hhpvec[ind], Chan.hvec[ind], Chan.nhhp[ind], Chan.nh[ind], hind1, hind2, pind1, hind3);
	Qmap2[i][2 * hhhp] = ind;
	Qmap2[i][2 * hhhp + 1] = ind1;
      }
    }
  }

  /*std::ofstream mapfile;
  mapfile.open ("mapfileN4.txt");
  for(int i = 0; i < Chan.size1; ++i){ mapfile << Chan.hh[i] * Chan.pp[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size2; ++i){ mapfile << Chan.hp1[i] * Chan.hp2[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size3; ++i){ mapfile << Chan.h[i] * Chan.hpp[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size3; ++i){ mapfile << Chan.p[i] * Chan.hhp[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.hh[i]*Chan.pp[i]; ++j){
      mapfile << i << " ";
      for(int k = 0; k < 15; ++k){ mapfile << Tmap[i][15*j + k] << " "; }
      mapfile << std::endl;
    }
  }
  mapfile.close();*/
}

void Doubles_1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(nhh * npp != 0){
      delete[] Tmap[i];
      delete[] Evec[i];
      delete[] T1[i];
    }
    if(nhh != 0){
      delete[] S1[i];
    }
    if(nhp * npp != 0){
      delete[] Q11[i];
      delete[] Qmap1[i];
    }
    if(nhh * nhp != 0){
      delete[] Q21[i];
      delete[] Qmap2[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nhpp * nh != 0){
      delete[] T6[i];
      delete[] T7[i];
    }
    if(nhhp * np != 0){
      delete[] T8[i];
      delete[] T9[i];
    }
    if(nhpp * np != 0){
      delete[] Q12[i];
    }
    if(nhhp * nh != 0){
      delete[] Q22[i];
    }
    if(nh != 0){
      delete[] S2[i];
      delete[] S3[i];
    }
    if(np != 0){
      delete[] S4[i];
      delete[] S5[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp1 * nhp2 != 0){
      delete[] T2[i];
      delete[] T3[i];
      delete[] T4[i];
      delete[] T5[i];
    }
    if(nhp2 != 0){
      delete[] S6[i];
      delete[] S7[i];
    }
  }
  delete[] Tmap;
  delete[] Evec;
  delete[] T1;
  delete[] T2;
  delete[] T3;
  delete[] T4;
  delete[] T5;
  delete[] T6;
  delete[] T7;
  delete[] T8;
  delete[] T9;
  delete[] S1;
  delete[] S2;
  delete[] S3;
  delete[] S4;
  delete[] S5;
  delete[] S6;
  delete[] S7;

  delete[] Q11;
  delete[] Q21;
  delete[] Qmap1;
  delete[] Qmap2;
  delete[] Q12;
  delete[] Q22;
}

void Doubles_1::zero(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T1[i][j] = 0.0;
      }
    }
    length = nhh * nhh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S1[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T6[i][j] = 0.0;
	T7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T8[i][j] = 0.0;
	T9[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(nhhp * nh != 0){
      for(int j = 0; j < length; ++j){
	Q22[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S2[i][j] = 0.0;
	S3[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S4[i][j] = 0.0;
	S5[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
	T4[i][j] = 0.0;
	T5[i][j] = 0.0;
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S6[i][j] = 0.0;
	S7[i][j] = 0.0;
      }
    }
  }
}

void Doubles_1::set_T(int i, int j, double T)
{
  T1[i][j] = T;
  T2[Tmap[i][16*j]][Tmap[i][16*j + 1]] = T;
  T3[Tmap[i][16*j + 2]][Tmap[i][16*j + 3]] = T;
  T4[Tmap[i][16*j + 4]][Tmap[i][16*j + 5]] = T;
  T5[Tmap[i][16*j + 6]][Tmap[i][16*j + 7]] = T;
  T6[Tmap[i][16*j + 8]][Tmap[i][16*j + 9]] = T;
  T7[Tmap[i][16*j + 10]][Tmap[i][16*j + 11]] = T;
  T8[Tmap[i][16*j + 12]][Tmap[i][16*j + 13]] = T;
  T9[Tmap[i][16*j + 14]][Tmap[i][16*j + 15]] = T;
}

double Doubles_1::get_T(int i, int j) const
{
  double tempt;
  tempt = T1[i][j];
  tempt += T2[Tmap[i][16*j]][Tmap[i][16*j + 1]];    
  tempt += T3[Tmap[i][16*j + 2]][Tmap[i][16*j + 3]];
  tempt += T4[Tmap[i][16*j + 4]][Tmap[i][16*j + 5]];
  tempt += T5[Tmap[i][16*j + 6]][Tmap[i][16*j + 7]];
  tempt += T6[Tmap[i][16*j + 8]][Tmap[i][16*j + 9]];
  tempt += T7[Tmap[i][16*j + 10]][Tmap[i][16*j + 11]];
  tempt += T8[Tmap[i][16*j + 12]][Tmap[i][16*j + 13]];
  tempt += T9[Tmap[i][16*j + 14]][Tmap[i][16*j + 15]];
  return tempt;
}

void Doubles_1::set_T_2(const Channels &Chan, Interactions &Ints)
{
  double fac1 = 1.0, fac2 = 0.0;
  char N = 'N';
  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.nhh[i];
    int pp = Chan.npp[i];
    int hp = Chan.nhp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    dgemm_NN(Ints.S_ME1.V19[i], T1[i], Q11[i], &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    dgemm_NN(T1[i], Ints.S_ME1.V20[i], Q21[i], &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	int ind0 = hp1 * Chan.npp[i] + pp1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
	int ind0 = hh1 * Chan.nhp[i] + hp1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }
}

Singles_1::Singles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  int length0, length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    Tmap = new int[3 * nhp1];
    Evec = new double[nhp1];
    T1 = new double[nhp1];
    Tmap2 = new int[18 * nhp1 * nhp1];
    S3 = new double[nhp1 * nhp1];
    S4 = new double[nhp1];
    for(int i = 0; i < nhp1; ++i){
      for(int j = 0; j < 3; ++j){ Tmap[3*i + j] = -1; }
      Evec[i] = 0.0;
      T1[i] = 0.0;
      S4[i] = 0.0;
      for(int j = 0; j < nhp1; ++j){
	for(int k = 0; k < 18; ++k){ Tmap2[18 * (nhp1 * i + j) + k] = -1; }
	S3[nhp1 * i +j] = 0.0;
      }
    }
  }
  if(nhh1 != 0){
    Q31 = new double[nhh1];
    Qmap3 = new int[2 * nhh1];
    for(int i = 0; i < nhh1; ++i){
      Q31[i] = 0.0;
      Qmap3[2*i] = -1;
      Qmap3[2*i + 1] = -1;
    }
  }
  if(npp1 != 0){
    Q41 = new double[npp1];
    Qmap4 = new int[2 * npp1];
    for(int i = 0; i < npp1; ++i){
      Q41[i] = 0.0;
      Qmap4[2*i] = -1;
      Qmap4[2*i + 1] = -1;
    }
  }

  T2 = new double*[Chan.size3];
  T3 = new double*[Chan.size3];

  E1 = new double*[Chan.size1];
  E2 = new double*[Chan.size2];
  E3 = new double*[Chan.size2];
  E4 = new double*[Chan.size2];
  E5 = new double*[Chan.size2];
  E6 = new double*[Chan.size3];
  E7 = new double*[Chan.size3];
  E8 = new double*[Chan.size3];
  E9 = new double*[Chan.size3];

  S1 = new double*[Chan.size3];
  S2 = new double*[Chan.size3];

  Q11 = new double*[Chan.size3];
  Q21 = new double*[Chan.size3];
  Qmap1 = new int*[Chan.size3];
  Qmap2 = new int*[Chan.size3];
  Q12 = new double*[Chan.size2];
  Q22 = new double*[Chan.size2];

  Q32 = new double*[Chan.size3];
  Q42 = new double*[Chan.size3];

  Q51 = new double*[Chan.size1];
  Q61 = new double*[Chan.size1];
  Qmap5 = new int*[Chan.size1];
  Qmap6 = new int*[Chan.size1];
  Q52 = new double*[Chan.size3];
  Q62 = new double*[Chan.size3];

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      E1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E1[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      Q61[i] = new double[length];
      Qmap6[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q61[i][j] = 0.0;
	Qmap6[i][2*j] = -1;
	Qmap6[i][2*j + 1] = -1;
      }
    }
    length = nhp * npp;
    if(length != 0){
      Q51[i] = new double[length];
      Qmap5[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q51[i][j] = 0.0;
	Qmap5[i][2*j] = -1;
	Qmap5[i][2*j + 1] = -1;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    length = nh * np;
    if(length != 0){
      T2[i] = new double[length];
      T3[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
      }
    }
    length = nhpp * nh;
    if(length != 0){
      E6[i] = new double[length];
      E7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E6[i][j] = 0.0;
	E7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      E8[i] = new double[length];
      E9[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E8[i][j] = 0.0;
	E9[i][j] = 0.0;
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      Q11[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      Q21[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      Q32[i] = new double[length];
      S2[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q32[i][j] = 0.0;
	S2[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      Q42[i] = new double[length];
      S1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q42[i][j] = 0.0;
	S1[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(length != 0){
      Q62[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q62[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      Q52[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q52[i][j] = 0.0;
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      Qmap1[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Qmap1[i][2*j] = -1;
	Qmap1[i][2*j + 1] = -1;
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      Qmap2[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Qmap2[i][2*j] = -1;
	Qmap2[i][2*j + 1] = -1;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      E2[i] = new double[length];
      E3[i] = new double[length];
      E4[i] = new double[length];
      E5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E2[i][j] = 0.0;
	E3[i][j] = 0.0;
	E4[i][j] = 0.0;
	E5[i][j] = 0.0;
      }
    }
    length = nhp1 * nhp1;
    if(length != 0){
      Q12[i] = new double[length];
      Q22[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
	Q22[i][j] = 0.0;
      }
    }
  }

  #pragma omp parallel
  {
    State tb;
    int ind, ind1, key1, key2;
    int hind1, pind1;
    #pragma omp for schedule(static)
    for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
      hind1 = Chan.hp1vec[Chan.ind0][2*hp1];
      pind1 = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      ind = Chan.indvec[hind1];
      Tmap[3 * hp1] = ind;
      key1 = Chan.p_map[ind][pind1];
      key2 = Chan.h_map[ind][hind1];
      ind1 = key1 * Chan.nh[ind] + key2;
      //ind1 = Index11(Chan.pvec[ind], Chan.hvec[ind], Chan.np[ind], Chan.nh[ind], pind1, hind1);
      Tmap[3 * hp1 + 1] = ind1;
      ind1 = key2 * Chan.np[ind] + key1;
      //ind1 = Index11(Chan.hvec[ind], Chan.pvec[ind], Chan.nh[ind], Chan.np[ind], hind1, pind1);
      Tmap[3 * hp1 + 2] = ind1;
    }
  }

  #pragma omp parallel
  {
    State tb;
    int ind, ind1, key1, key2;
    int hind1, hind2, pind1, pind2;
    int hp1, hp2;
    length = Chan.nhp1[Chan.ind0] * Chan.nhp1[Chan.ind0];
    length0 = Chan.nhp1[Chan.ind0];
    #pragma omp for schedule(static)
    for(int hphp = 0; hphp < length; ++hphp){
      hp2 = int(hphp%length0);
      hp1 = int((hphp - hp2)/length0);
      hind1 = Chan.hp1vec[Chan.ind0][2*hp1];
      pind2 = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      hind2 = Chan.hp1vec[Chan.ind0][2*hp2];
      pind1 = Chan.hp1vec[Chan.ind0][2*hp2 + 1];

      //E1 = ((ij)(ab))
      plus(tb, Space.qnums[hind1], Space.qnums[hind2]);
      ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
      key1 = Chan.hh_map[ind][Hash2(hind1, hind2, Space.indtot)];
      key2 = Chan.pp_map[ind][Hash2(pind1, pind2, Space.indtot)];
      ind1 = key1 * Chan.npp[ind] + key2;
      //ind1 = Index22(Chan.hhvec[ind], Chan.ppvec[ind], Chan.nhh[ind], Chan.npp[ind], hind1, hind2, pind1, pind2);
      Tmap2[18 * hphp] = ind;
      Tmap2[18 * hphp + 1] = ind1;
      //E2 = ((ia)(jb)')
      minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
      key2 = Chan.hp2_map[ind][Hash2(hind2, pind2, Space.indtot)];
      ind1 = key1 * Chan.nhp2[ind] + key2;
      //ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind1, pind1, hind2, pind2);
      Tmap2[18 * hphp + 2] = ind;
      Tmap2[18 * hphp + 3] = ind1;
      //E3 = ((jb)(ia)')
      minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
      key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
      ind1 = key1 * Chan.nhp2[ind] + key2;
      //ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind2, pind2, hind1, pind1);
      Tmap2[18 * hphp + 4] = ind;
      Tmap2[18 * hphp + 5] = ind1;
      //E4 = ((ib)(ja)')
      minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
      key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
      ind1 = key1 * Chan.nhp2[ind] + key2;
      //ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind1, pind2, hind2, pind1);
      Tmap2[18 * hphp + 6] = ind;
      Tmap2[18 * hphp + 7] = ind1;
      //E5 = ((ja)(ib)')
      minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
      key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
      ind1 = key1 * Chan.nhp2[ind] + key2;
      //ind1 = Index22(Chan.hp1vec[ind], Chan.hp2vec[ind], Chan.nhp1[ind], Chan.nhp2[ind], hind2, pind1, hind1, pind2);
      Tmap2[18 * hphp + 8] = ind;
      Tmap2[18 * hphp + 9] = ind1;
      //E6 = ((jab)(i))
      ind = Chan.indvec[hind1];
      key1 = Chan.hpp_map[ind][Hash3(hind2, pind1, pind2, Space.indtot)];
      key2 = Chan.h_map[ind][hind1];
      ind1 = key1 * Chan.nh[ind] + key2;
      //ind1 = Index31(Chan.hppvec[ind], Chan.hvec[ind], Chan.nhpp[ind], Chan.nh[ind], hind2, pind1, pind2, hind1);
      Tmap2[18 * hphp + 10] = ind;
      Tmap2[18 * hphp + 11] = ind1;
      //E7 = ((iab)(j))
      ind = Chan.indvec[hind2];
      key1 = Chan.hpp_map[ind][Hash3(hind1, pind1, pind2, Space.indtot)];
      key2 = Chan.h_map[ind][hind2];
      ind1 = key1 * Chan.nh[ind] + key2;
      //ind1 = Index31(Chan.hppvec[ind], Chan.hvec[ind], Chan.nhpp[ind], Chan.nh[ind], hind1, pind1, pind2, hind2);
      Tmap2[18 * hphp + 12] = ind;
      Tmap2[18 * hphp + 13] = ind1;
      //E8 = ((ijb)(a))
      ind = Chan.indvec[pind1];
      key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind2, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      //ind1 = Index31(Chan.hhpvec[ind], Chan.pvec[ind], Chan.nhhp[ind], Chan.np[ind], hind1, hind2, pind2, pind1);
      Tmap2[18 * hphp + 14] = ind;
      Tmap2[18 * hphp + 15] = ind1;
      //E9 = ((ija)(b))
      ind = Chan.indvec[pind2];
      key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind1, Space.indtot)];
      key2 = Chan.p_map[ind][pind2];
      ind1 = key1 * Chan.np[ind] + key2;
      //ind1 = Index31(Chan.hhpvec[ind], Chan.pvec[ind], Chan.nhhp[ind], Chan.np[ind], hind1, hind2, pind1, pind2);
      Tmap2[18 * hphp + 16] = ind;
      Tmap2[18 * hphp + 17] = ind1;
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, pind1, pind2;
      int h1, hpp1, p1, hhp1;
      length = Chan.nhpp1[i] * Chan.nh[i];
      length0 = Chan.nh[i];
      #pragma omp for schedule(static)
      for(int hpph = 0; hpph < length; ++hpph){
	h1 = int(hpph%length0);
	hpp1 = int((hpph - h1)/length0);
	hind2 = Chan.hpp1vec[i][3*hpp1];
	pind1 = Chan.hpp1vec[i][3*hpp1 + 1];
	pind2 = Chan.hpp1vec[i][3*hpp1 + 2];
	hind1 = Chan.hvec[i][h1];
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	//ind1 = Index22(Chan.hp1vec[ind], Chan.hp1vec[ind], Chan.nhp1[ind], Chan.nhp1[ind], hind1, pind1, hind2, pind2);
	Qmap1[i][2 * hpph] = ind;
	Qmap1[i][2 * hpph + 1] = ind1;
      }
      length = Chan.nhhp1[i] * Chan.np[i];
      length0 = Chan.np[i];
      #pragma omp for schedule(static)
      for(int hhpp = 0; hhpp < length; ++hhpp){
	p1 = int(hhpp%length0);
	hhp1 = int((hhpp - p1)/length0);
	hind2 = Chan.hhp1vec[i][3*hhp1];
	hind1 = Chan.hhp1vec[i][3*hhp1 + 1];
	pind2 = Chan.hhp1vec[i][3*hhp1 + 2];
	pind1 = Chan.pvec[i][p1];
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	//ind1 = Index22(Chan.hp1vec[ind], Chan.hp1vec[ind], Chan.nhp1[ind], Chan.nhp1[ind], hind1, pind1, hind2, pind2);
	Qmap2[i][2 * hhpp] = ind;
	Qmap2[i][2 * hhpp + 1] = ind1;
      }
    }
  }

  #pragma omp parallel
  {
    State tb;
    int ind, ind1, key1, key2;
    int hind1, hind2, pind1, pind2;
    #pragma omp for schedule(static)
    for(int hh1 = 0; hh1 < Chan.nhh1[Chan.ind0]; ++hh1){
      hind1 = Chan.hh1vec[Chan.ind0][2*hh1];
      hind2 = Chan.hh1vec[Chan.ind0][2*hh1 + 1];
      ind = Chan.indvec[hind1];
      key1 = Chan.h_map[ind][hind2];
      key2 = Chan.h_map[ind][hind1];
      ind1 = key1 * Chan.nh[ind] + key2;
      //ind1 = Index11(Chan.hvec[ind], Chan.hvec[ind], Chan.nh[ind], Chan.nh[ind], hind2, hind1);
      Qmap3[2 * hh1] = ind;
      Qmap3[2 * hh1 + 1] = ind1;
    }
    #pragma omp for schedule(static)
    for(int pp1 = 0; pp1 < Chan.npp1[Chan.ind0]; ++pp1){
      pind1 = Chan.pp1vec[Chan.ind0][2*pp1];
      pind2 = Chan.pp1vec[Chan.ind0][2*pp1 + 1];
      ind = Chan.indvec[pind1];
      key1 = Chan.p_map[ind][pind1];
      key2 = Chan.p_map[ind][pind2];
      ind1 = key1 * Chan.np[ind] + key2;
      //ind1 = Index11(Chan.pvec[ind], Chan.pvec[ind], Chan.np[ind], Chan.np[ind], pind1, pind2);
      Qmap4[2 * pp1] = ind;
      Qmap4[2 * pp1 + 1] = ind1;
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, hind3, pind1, pind2, pind3;
      int hp1, pp1, hh1;
      length = Chan.nhp[i] * Chan.npp[i];
      length0 = Chan.npp[i];
      #pragma omp for schedule(static)
      for(int hppp = 0; hppp < length; ++hppp){
	pp1 = int(hppp%length0);
	hp1 = int((hppp - pp1)/length0);
	hind1 = Chan.hpvec[i][2*hp1];
	pind1 = Chan.hpvec[i][2*hp1 + 1];
	pind2 = Chan.ppvec[i][2*pp1];
	pind3 = Chan.ppvec[i][2*pp1 + 1];
	ind = Chan.indvec[pind1];
	key1 = Chan.hpp_map[ind][Hash3(hind1, pind2, pind3, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hppvec[ind], Chan.pvec[ind], Chan.nhpp[ind], Chan.np[ind], hind1, pind2, pind3, pind1);
	Qmap5[i][2 * hppp] = ind;
	Qmap5[i][2 * hppp + 1] = ind1;
      }
      length = Chan.nhh[i] * Chan.nhp[i];
      length0 = Chan.nhp[i];
      #pragma omp for schedule(static)
      for(int hhhp = 0; hhhp < length; ++hhhp){
	hp1 = int(hhhp%length0);
	hh1 = int((hhhp - hp1)/length0);
	hind1 = Chan.hhvec[i][2*hh1];
	hind2 = Chan.hhvec[i][2*hh1 + 1];
	hind3 = Chan.hpvec[i][2*hp1];
	pind1 = Chan.hpvec[i][2*hp1 + 1];
	ind = Chan.indvec[hind3];
	key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind1, Space.indtot)];
	key2 = Chan.h_map[ind][hind3];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hhpvec[ind], Chan.hvec[ind], Chan.nhhp[ind], Chan.nh[ind], hind1, hind2, pind1, hind3);
	Qmap6[i][2 * hhhp] = ind;
	Qmap6[i][2 * hhhp + 1] = ind1;
      }
    }
  }
}

void Singles_1::zero(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    for(int i = 0; i < nhp1; ++i){
      T1[i] = 0.0;
      S4[i] = 0.0;
      for(int j = 0; j < nhp1; ++j){
	S3[nhp1 * i +j] = 0.0;
      }
    }
  }
  if(nhh1 != 0){
    for(int i = 0; i < nhh1; ++i){
      Q31[i] = 0.0;
    }
  }
  if(npp1 != 0){
    for(int i = 0; i < npp1; ++i){
      Q41[i] = 0.0;
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E1[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q61[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q51[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    length = nh * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
      }
    }
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E6[i][j] = 0.0;
	E7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E8[i][j] = 0.0;
	E9[i][j] = 0.0;
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q32[i][j] = 0.0;
	S2[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q42[i][j] = 0.0;
	S1[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q52[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q62[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E2[i][j] = 0.0;
	E3[i][j] = 0.0;
	E4[i][j] = 0.0;
	E5[i][j] = 0.0;
      }
    }
    length = nhp1 * nhp1;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
	Q22[i][j] = 0.0;
      }
    }
  }
}

void Singles_1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    delete[] Tmap;
    delete[] Evec;
    delete[] T1;
    delete[] Tmap2;
    delete[] S3;
    delete[] S4;
  }
  if(nhh1 != 0){
    delete[] Q31;
    delete[] Qmap3;
  }
  if(npp1 != 0){
    delete[] Q41;
    delete[] Qmap4;
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(nhh * npp != 0){
      delete[] E1[i];
    }
    if(nhh * nhp != 0){
      delete[] Q61[i];
      delete[] Qmap6[i];
    }
    if(nhp * npp != 0){
      delete[] Q51[i];
      delete[] Qmap5[i];
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(nh * np != 0){
      delete[] T2[i];
      delete[] T3[i];
    }
    if(nhpp * nh != 0){
      delete[] E6[i];
      delete[] E7[i];
    }
    if(nhhp * np != 0){
      delete[] E8[i];
      delete[] E9[i];
    }
    if(nhpp1 * nh != 0){
      delete[] Q11[i];
    }
    if(nhhp1 * np != 0){
      delete[] Q21[i];
    }
    if(nh != 0){
      delete[] Q32[i];
      delete[] S2[i];
    }
    if(np != 0){
      delete[] Q42[i];
      delete[] S1[i];
    }
    if(nhpp * np != 0){
      delete[] Q52[i];
    }
    if(nhhp * nh != 0){
      delete[] Q62[i];
    }
    if(nhpp1 * nh != 0){
      delete[] Qmap1[i];
    }
    if(nhhp1 * np != 0){
      delete[] Qmap2[i];
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp1 * nhp2 != 0){
      delete[] E2[i];
      delete[] E3[i];
      delete[] E4[i];
      delete[] E5[i];
    }
    if(nhp1 != 0){
      delete[] Q12[i];
      delete[] Q22[i];
    }
  }

  delete[] T2;
  delete[] T3;
  delete[] E1;
  delete[] E2;
  delete[] E3;
  delete[] E4;
  delete[] E5;
  delete[] E6;
  delete[] E7;
  delete[] E8;
  delete[] E9;
  delete[] S1;
  delete[] S2;
  delete[] Q11;
  delete[] Q21;
  delete[] Qmap1;
  delete[] Qmap2;
  delete[] Q12;
  delete[] Q22;
  delete[] Q32;
  delete[] Q42;
  delete[] Q51;
  delete[] Q61;
  delete[] Qmap5;
  delete[] Qmap6;
  delete[] Q52;
  delete[] Q62;
}

void Singles_1::set_T(int i, double T)
{
  T1[i] = T;
  T2[Tmap[3*i]][Tmap[3*i + 1]] = T;
  T3[Tmap[3*i]][Tmap[3*i + 2]] = T;
}

double Singles_1::get_T(int i) const
{
  double tempt = T1[i];
  tempt += T2[Tmap[3*i]][Tmap[3*i + 1]];
  tempt += T3[Tmap[3*i]][Tmap[3*i + 2]];
  return tempt;
}

void Singles_1::set_T_2(const Channels &Chan, Interactions &Ints)
{
  for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
    for(int hp2 = 0; hp2 < Chan.nhp1[Chan.ind0]; ++hp2){
      int ind0 = hp1 * Chan.nhp1[Chan.ind0] + hp2;
      double T = T1[hp1] * T1[hp2];
      E1[Tmap2[18 * ind0]][Tmap2[18 * ind0 + 1]] = T;
      E2[Tmap2[18 * ind0 + 2]][Tmap2[18 * ind0 + 3]] = T;
      E3[Tmap2[18 * ind0 + 4]][Tmap2[18 * ind0 + 5]] = T;
      E4[Tmap2[18 * ind0 + 6]][Tmap2[18 * ind0 + 7]] = T;
      E5[Tmap2[18 * ind0 + 8]][Tmap2[18 * ind0 + 9]] = T;
      E6[Tmap2[18 * ind0 + 10]][Tmap2[18 * ind0 + 11]] = T;
      E7[Tmap2[18 * ind0 + 12]][Tmap2[18 * ind0 + 13]] = T;
      E8[Tmap2[18 * ind0 + 14]][Tmap2[18 * ind0 + 15]] = T;
      E9[Tmap2[18 * ind0 + 16]][Tmap2[18 * ind0 + 17]] = T;
    }
  }

  double fac1 = 1.0, fac2 = 0.0;
  int one = 1;
  char N = 'N';
  for(int i = 0; i < Chan.size3; ++i){
    int h = Chan.nh[i];
    int p = Chan.np[i];
    int hpp1 = Chan.nhpp1[i];
    int hhp1 = Chan.nhhp1[i];
    if(h == 0 || p == 0){ continue; }
    if(hpp1 != 0){ dgemm_NN(Ints.S_ME1.V13[i], T2[i], Q11[i], &hpp1, &h, &p, &fac1, &fac2, &N, &N); }
    if(hhp1 != 0){ dgemm_NN(Ints.S_ME1.V14[i], T3[i], Q21[i], &hhp1, &p, &h, &fac1, &fac2, &N, &N); }
    for(int hpp2 = 0; hpp2 < hpp1; ++hpp2){
      for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
	int ind0 = hpp2 * Chan.nh[i] + h1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hhp2 = 0; hhp2 < hhp1; ++hhp2){
      for(int p1 = 0; p1 < Chan.np[i]; ++p1){
	int ind0 = hhp2 * Chan.np[i] + p1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }

  int hp = Chan.nhp1[Chan.ind0];
  int hh = Chan.nhh1[Chan.ind0];
  int pp = Chan.npp1[Chan.ind0];
  if(hp != 0 && hh != 0){ dgemm_NN(Ints.S_ME1.V15[Chan.ind0], T1, Q31, &hh, &one, &hp, &fac1, &fac2, &N, &N); }
  if(hp != 0 && pp != 0){ dgemm_NN(Ints.S_ME1.V16[Chan.ind0], T1, Q41, &pp, &one, &hp, &fac1, &fac2, &N, &N); }
  for(int hh1 = 0; hh1 < Chan.nhh1[Chan.ind0]; ++hh1){ Q32[Qmap3[2 * hh1]][Qmap3[2 * hh1 + 1]] = Q31[hh1]; }
  for(int pp1 = 0; pp1 < Chan.npp1[Chan.ind0]; ++pp1){ Q42[Qmap4[2 * pp1]][Qmap4[2 * pp1 + 1]] = Q41[pp1]; }

  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.nhh[i];
    int pp = Chan.npp[i];
    int hp = Chan.nhp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    dgemm_NN(Ints.S_ME1.V19[i], E1[i], Q51[i], &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    dgemm_NN(E1[i], Ints.S_ME1.V20[i], Q61[i], &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	int ind0 = hp1 * Chan.npp[i] + pp1;
	Q52[Qmap5[i][2 * ind0]][Qmap5[i][2 * ind0 + 1]] = Q51[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
	int ind0 = hh1 * Chan.nhp[i] + hp1;
	Q62[Qmap6[i][2 * ind0]][Qmap6[i][2 * ind0 + 1]] = Q61[i][ind0];
      }
    }
  }
}


Interactions::Interactions(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D_ME1 = Doubles_ME1(Chan);
  }
  else if(Parameters.approx == "singles"){
    D_ME1 = Doubles_ME1(Chan);
    S_ME1 = Singles_ME1(Chan);
  }
  else if(Parameters.approx == "triples"){
    D_ME1 = Doubles_ME1(Chan);
    S_ME1 = Singles_ME1(Chan);
  }
}

void Interactions::delete_struct(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D_ME1.delete_struct(Chan);;
  }
  else if(Parameters.approx == "singles"){
    D_ME1.delete_struct(Chan);
    S_ME1.delete_struct(Chan);
  }
  else if(Parameters.approx == "triples"){
    D_ME1.delete_struct(Chan);
    S_ME1.delete_struct(Chan);
  }
}

Doubles_ME1::Doubles_ME1(const Channels &Chan)
{
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  V1 = new double*[Chan.size1];
  V2 = new double*[Chan.size1];
  V3 = new double*[Chan.size2];
  V4 = new double*[Chan.size1];
  V5 = new double*[Chan.size3];
  V6 = new double*[Chan.size3];
  V7 = new double*[Chan.size3];
  V8 = new double*[Chan.size3];
  V9 = new double*[Chan.size2];
  V10 = new double*[Chan.size2];
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    if(npp != 0){
      V1[i] = new double[npp * npp];
      for(int j = 0; j < npp * npp; ++j){ V1[i][j] = 0.0; }
    }
    if(nhh != 0){
      V2[i] = new double[nhh * nhh];
      for(int j = 0; j < nhh * nhh; ++j){ V2[i][j] = 0.0; }
    }
    if(npp * nhh != 0){
      V4[i] = new double[npp * nhh];
      for(int j = 0; j < npp * nhh; ++j){ V4[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nh * nhpp != 0){
      V5[i] = new double[nh * nhpp];
      V6[i] = new double[nh * nhpp];
      for(int j = 0; j < nh * nhpp; ++j){ V5[i][j] = 0.0; V6[i][j] = 0.0; }
    }
    if(np * nhhp != 0){
      V7[i] = new double[np * nhhp];
      V8[i] = new double[np * nhhp];
      for(int j = 0; j < np * nhhp; ++j){ V7[i][j] = 0.0; V8[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp2 != 0){
      V3[i] = new double[nhp2 * nhp2];
      for(int j = 0; j < nhp2 * nhp2; ++j){ V3[i][j] = 0.0; }
    }
    if(nhp2 * nhp1){
      V9[i] = new double[nhp2 * nhp1];
      V10[i] = new double[nhp2 * nhp1];
      for(int j = 0; j < nhp2 * nhp1; ++j){ V9[i][j] = 0.0; V10[i][j] = 0.0; }
    }
  }
}

void Doubles_ME1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    if(npp != 0){
      delete[] V1[i];
    }
    if(nhh != 0){
      delete[] V2[i];
    }
    if(npp * nhh != 0){
      delete[] V4[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nh * nhpp != 0){
      delete[] V5[i];
      delete[] V6[i];
    }
    if(np * nhhp != 0){
      delete[] V7[i];
      delete[] V8[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp2 != 0){
      delete[] V3[i];
    }
    if(nhp2 * nhp1){
      delete[] V9[i];
      delete[] V10[i];
    }
  }
  delete[] V1;
  delete[] V2;
  delete[] V3;
  delete[] V4;
  delete[] V5;
  delete[] V6;
  delete[] V7;
  delete[] V8;
  delete[] V9;
  delete[] V10;
}

Singles_ME1::Singles_ME1(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhh1, npp1, nhpp1, nhhp1;
  V11 = new double*[Chan.size3];
  V12 = new double*[Chan.size3];
  V17 = new double*[Chan.size3];
  V18 = new double*[Chan.size3];
  V13 = new double*[Chan.size3];
  V14 = new double*[Chan.size3];
  V20 = new double*[Chan.size1];
  V19 = new double*[Chan.size1];
  V16 = new double*[Chan.size2];
  V15 = new double*[Chan.size2];
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(npp * nhp != 0){
      V20[i] = new double[npp * nhp];
      for(int j = 0; j < npp * nhp; ++j){ V20[i][j] = 0.0; }
    }
    if(nhp * nhh != 0){
      V19[i] = new double[nhp * nhh];
      for(int j = 0; j < nhp * nhh; ++j){ V19[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(np * nhpp != 0){
      V11[i] = new double[np * nhpp];
      V17[i] = new double[nhpp * np];
      for(int j = 0; j < np * nhpp; ++j){ V11[i][j] = 0.0; V17[i][j] = 0.0; }
    }
    if(nh * nhhp != 0){
      V12[i] = new double[nh * nhhp];
      V18[i] = new double[nhhp * nh];
      for(int j = 0; j < nh * nhhp; ++j){ V12[i][j] = 0.0; V18[i][j] = 0.0; }
    }
    if(nhpp1 * np != 0){
      V13[i] = new double[nhpp1 * np];
      for(int j = 0; j < nhpp1 * np; ++j){ V13[i][j] = 0.0; }
    }
    if(nhhp1 * nh != 0){
      V14[i] = new double[nhhp1 * nh];
      for(int j = 0; j < nhhp1 * nh; ++j){ V14[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    if(npp1 * nhp1 != 0){
      V16[i] = new double[npp1 * nhp1];
      for(int j = 0; j < npp1 * nhp1; ++j){ V16[i][j] = 0.0; }
    }
    if(nhh1 * nhp1 != 0){
      V15[i] = new double[nhh1 * nhp1];
      for(int j = 0; j < nhh1 * nhp1; ++j){ V15[i][j] = 0.0; }
    }
  }
}

void Singles_ME1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhh1, npp1, nhpp1, nhhp1;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(npp * nhp != 0){
      delete[] V20[i];
    }
    if(nhp * nhh != 0){
      delete[] V19[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(np * nhpp != 0){
      delete[] V11[i];
      delete[] V17[i];
    }
    if(nh * nhhp != 0){
      delete[] V12[i];
      delete[] V18[i];
    }
    if(nhpp1 * np != 0){
      delete[] V13[i];
    }
    if(nhhp1 * nh != 0){
      delete[] V14[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    if(npp1 * nhp1 != 0){
      delete[] V16[i];
    }
    if(nhh1 * nhp1 != 0){
      delete[] V15[i];
    }
  }
  delete[] V11;
  delete[] V12;
  delete[] V17;
  delete[] V18;
  delete[] V13;
  delete[] V14;
  delete[] V20;
  delete[] V19;
  delete[] V16;
  delete[] V15;
}


//Function to setup Channels
Channels::Channels(const Input_Parameters &Parameters, const Model_Space &Space)
{
  std::cout << "Building Channels ...  ";

  State state;
  int key;
  int count0, ind1, ind2, h, p, hh, pp, hp, hp1, hp2, hh1, pp1, hpp, hhp, hpp1, hhp1, hhh, ppp;
  State *qnumstemp = new State[Space.indtot]; // max number of qnums groups
  int *hnum = new int[Space.indtot]; // count holes in each qnums group
  int *pnum = new int[Space.indtot]; // count particles in each qnums group
  indvec = new int[Space.indtot]; // index of qnums group for each state
  for(int i = 0; i < Space.indtot; ++i){ hnum[i] = 0; pnum[i] = 0; }

  qnumstemp[0] = Space.qnums[0];
  indvec[0] = 0;
  if( Space.qnums[0].type == "hole" ){ ++hnum[0]; }
  else{ ++pnum[0]; }

  // count # of qnums groups and # of hs and ps in each qnums group, and fill indvec
  count0 = 1;
  for(int i = 1; i < Space.indtot; ++i){
    state = Space.qnums[i];
    for(int k = 0; k < count0; ++k){
      if( equal(state, qnumstemp[k]) ){
	indvec[i] = k;
	if(Space.qnums[i].type == "hole"){ ++hnum[k]; }
	else{ ++pnum[k]; }
	goto stop;
      }
      if(k == count0 - 1){
	qnumstemp[count0] = Space.qnums[i];
	indvec[i] = count0;
	if(Space.qnums[i].type == "hole"){ ++hnum[count0]; }
	else{ ++pnum[count0]; }
	++count0;
	break;
      }
    }
  stop:;
  }
  size3 = count0;

  // allocate memory for Hvec and Pvec, reset hnum and pnum
  qnums3 = new State[size3];
  hvec = new int*[size3];
  pvec = new int*[size3];
  nh = new int[size3];
  np = new int[size3];
  h_map = new std::unordered_map<int,int>[size3];
  p_map = new std::unordered_map<int,int>[size3];
  for(int i = 0; i < size3; ++i){
    h = hnum[i];
    p = pnum[i];
    qnums3[i] = qnumstemp[i];
    if(h != 0){ hvec[i] = new int[h]; }
    if(p != 0){ pvec[i] = new int[p]; }
    nh[i] = 0;
    np[i] = 0;
  }
  delete[] qnumstemp;
  delete[] hnum;
  delete[] pnum;

  // place states in appropriate Hvec or Pvec position
  for(int i = 0; i < Space.indtot; ++i){
    for(int k = 0; k < size3; ++k){
      if( equal(Space.qnums[i], qnums3[k]) ){
	if(Space.qnums[i].type == "hole"){
	  hvec[k][nh[k]] = i;
	  h_map[k][i] = nh[k];
	  ++nh[k];
	}
	else{
	  pvec[k][np[k]] = i;
	  p_map[k][i] = np[k];
	  ++np[k];
	}
	break;
      }
    }
  }

  /*for(int i = 0; i < size3; ++i){
    std::cout << "Chan3 = " << i << ", " << qnums3[i].par << " " << qnums3[i].ml << " " << qnums3[i].m << std::endl;
    for(int j = 0; j < nh[i]; ++j){ std::cout << hvec[i][j] << " "; }
    std::cout << std::endl;
    for(int j = 0; j < np[i]; ++j){ std::cout << pvec[i][j] << " "; }
    std::cout << std::endl << std::endl;
    }*/

  size1 = Space.size_2b;
  size2 = Space.size_2b;
  std::cout << " Size1 = " << size1 << ", Size2 = " << size2 << ", Size3 = " << size3 << std::endl;

  qnums1 = new State[size1];
  qnums2 = new State[size2];

  hhvec = new int*[size1];
  ppvec = new int*[size1];
  hpvec = new int*[size1];
  hp1vec = new int*[size2];
  hp2vec = new int*[size2];
  hh1vec = new int*[size2];
  pp1vec = new int*[size2];

  hh_map = new std::unordered_map<int,int>[size1];
  pp_map = new std::unordered_map<int,int>[size1];
  hp_map = new std::unordered_map<int,int>[size1];
  hp1_map = new std::unordered_map<int,int>[size2];
  hp2_map = new std::unordered_map<int,int>[size2];
  hh1_map = new std::unordered_map<int,int>[size2];
  pp1_map = new std::unordered_map<int,int>[size2];

  nhh = new int[size1];
  npp = new int[size1];
  nhp = new int[size1];
  nhp1 = new int[size2];
  nhp2 = new int[size2];
  nhh1 = new int[size2];
  npp1 = new int[size2];

  for(int chan1 = 0; chan1 < size1; ++chan1){
    nhh[chan1] = 0;
    npp[chan1] = 0;
    nhp[chan1] = 0;
  }
  for(int chan2 = 0; chan2 < size2; ++chan2){
    nhp1[chan2] = 0;
    nhp2[chan2] = 0;
    nhh1[chan2] = 0;
    npp1[chan2] = 0;
  }

  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      plus(state, Space.qnums[i], Space.qnums[j]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      qnums1[ind1] = state;
      ++nhh[ind1];
      minus(state, Space.qnums[i], Space.qnums[j]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      qnums2[ind2] = state;
      ++nhh1[ind2];
    }
  }

  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      plus(state, Space.qnums[a], Space.qnums[b]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      qnums1[ind1] = state;
      ++npp[ind1];
      minus(state, Space.qnums[a], Space.qnums[b]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      qnums2[ind2] = state;
      ++npp1[ind2];
    }
  }

  for(int i = 0; i < Space.indhol; ++i){
    for(int a = Space.indhol; a < Space.indtot; ++a){
      plus(state, Space.qnums[i], Space.qnums[a]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      qnums1[ind1] = state;
      ++nhp[ind1];
      minus(state, Space.qnums[i], Space.qnums[a]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      qnums2[ind2] = state;
      ++nhp1[ind2];
      minus(state, Space.qnums[a], Space.qnums[i]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      qnums2[ind2] = state;
      ++nhp2[ind2];
    }
  }

  for(int i = 0; i < size1; ++i){
    hh = nhh[i];
    pp = npp[i];
    hp = nhp[i];
    if(hh != 0){ hhvec[i] = new int[2 * hh]; }
    if(pp != 0){ ppvec[i] = new int[2 * pp]; }
    if(hp != 0){ hpvec[i] = new int[2 * hp]; }
    nhh[i] = 0;
    npp[i] = 0;
    nhp[i] = 0;
  }
  for(int i = 0; i < size2; ++i){
    hp1 = nhp1[i];
    hp2 = nhp2[i];
    hh1 = nhh1[i];
    pp1 = npp1[i];
    if(hp1 != 0){ hp1vec[i] = new int[2 * hp1]; }
    if(hp2 != 0){ hp2vec[i] = new int[2 * hp2]; }
    if(hh1 != 0){ hh1vec[i] = new int[2 * hh1]; }
    if(pp1 != 0){ pp1vec[i] = new int[2 * pp1]; }
    nhp1[i] = 0;
    nhp2[i] = 0;
    nhh1[i] = 0;
    npp1[i] = 0;
  }

  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      key = Hash2(i, j, Space.indtot);
      plus(state, Space.qnums[i], Space.qnums[j]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      hhvec[ind1][2 * nhh[ind1]] = i;
      hhvec[ind1][2 * nhh[ind1] + 1] = j;
      hh_map[ind1][key] = nhh[ind1];
      ++nhh[ind1];
      minus(state, Space.qnums[i], Space.qnums[j]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      hh1vec[ind2][2 * nhh1[ind2]] = i;
      hh1vec[ind2][2 * nhh1[ind2] + 1] = j;
      hh1_map[ind2][key] = nhh1[ind2];
      ++nhh1[ind2];
    }
  }

  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      key = Hash2(a, b, Space.indtot);
      plus(state, Space.qnums[a], Space.qnums[b]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      ppvec[ind1][2 * npp[ind1]] = a;
      ppvec[ind1][2 * npp[ind1] + 1] = b;
      pp_map[ind1][key] = npp[ind1];
      ++npp[ind1];
      minus(state, Space.qnums[a], Space.qnums[b]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      pp1vec[ind2][2 * npp1[ind2]] = a;
      pp1vec[ind2][2 * npp1[ind2] + 1] = b;
      pp1_map[ind2][key] = npp1[ind2];
      ++npp1[ind2];
    }
  }

  for(int i = 0; i < Space.indhol; ++i){
    for(int a = Space.indhol; a < Space.indtot; ++a){
      key = Hash2(i, a, Space.indtot);
      plus(state, Space.qnums[i], Space.qnums[a]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      hpvec[ind1][2 * nhp[ind1]] = i;
      hpvec[ind1][2 * nhp[ind1] + 1] = a;
      hp_map[ind1][key] = nhp[ind1];
      ++nhp[ind1];
      minus(state, Space.qnums[i], Space.qnums[a]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      hp1vec[ind2][2 * nhp1[ind2]] = i;
      hp1vec[ind2][2 * nhp1[ind2] + 1] = a;
      hp1_map[ind2][key] = nhp1[ind2];
      ++nhp1[ind2];
      minus(state, Space.qnums[a], Space.qnums[i]);
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
      hp2vec[ind2][2 * nhp2[ind2]] = i;
      hp2vec[ind2][2 * nhp2[ind2] + 1] = a;
      hp2_map[ind2][key] = nhp2[ind2];
      ++nhp2[ind2];
    }
  }
  
  // for singles case
  if(Parameters.approx == "singles" || Parameters.approx == "triples"){
    state.t = 0; //tz
    state.m = 0; //jz
    state.par = 1; //par
    state.ml = 0; //lz
    state.j = 0; //j
    ind0 = ChanInd_2b_cross(Parameters.basis, Space, state);
  }

  hppvec = new int*[size3];
  hhpvec = new int*[size3];
  hpp1vec = new int*[size3];
  hhp1vec = new int*[size3];
  hhhvec = new int*[size3];
  pppvec = new int*[size3];

  hpp_map = new std::unordered_map<int,int>[size3];
  hhp_map = new std::unordered_map<int,int>[size3];
  hpp1_map = new std::unordered_map<int,int>[size3];
  hhp1_map = new std::unordered_map<int,int>[size3];
  hhh_map = new std::unordered_map<int,int>[size3];
  ppp_map = new std::unordered_map<int,int>[size3];

  nhpp = new int[size3];
  nhhp = new int[size3];
  nhpp1 = new int[size3];
  nhhp1 = new int[size3];
  nhhh = new int[size3];
  nppp = new int[size3];

  for(int i = 0; i < size3; ++i){
    nhpp[i] = 0;
    nhhp[i] = 0;
    nhpp1[i] = 0;
    nhhp1[i] = 0;
    nhhh[i] = 0;
    nppp[i] = 0;
  }

  for(int i = 0; i < size3; ++i){
    for(int j = 0; j < size3; ++j){
      h = nh[j];
      p = np[j];
      plus(state, qnums3[i], qnums3[j]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      hh = nhh[ind1];
      pp = npp[ind1];
      hp = nhp[ind1];
      nhhp[i] += p * hh; // <p1p2|h1h2> -> <p1|h1h2p2>
      nhpp[i] += h * pp; // <h1h2|p1p2> -> <h1|h2p1p2>
      nhhp1[i] += h * hp; // <h1h2|h3p1> -> <h1|h2h3p1>
      nhpp1[i] += p * hp; // <p1p2|h1p3> -> <p1|h1p3p2>
      nhhh[i] += h * hh; // <p1h1|h2h3> -> <p1|h1h2h3>
      nppp[i] += p * pp; // <h1p1|p2p3> -> <h1|p1p2p3>
    }
  }

  for(int i = 0; i < size3; ++i){
    hpp = nhpp[i];
    hhp = nhhp[i];
    hpp1 = nhpp1[i];
    hhp1 = nhhp1[i];
    hhh = nhhh[i];
    ppp = nppp[i];
    if(hpp != 0){ hppvec[i] = new int[3 * hpp]; }
    if(hhp != 0){ hhpvec[i] = new int[3 * hhp]; }
    if(hpp1 != 0){ hpp1vec[i] = new int[3 * hpp1]; }
    if(hhp1 != 0){ hhp1vec[i] = new int[3 * hhp1]; }
    if(hhh != 0){ hhhvec[i] = new int[3 * hhh]; }
    if(ppp != 0){ pppvec[i] = new int[3 * ppp]; }

    nhpp[i] = 0;
    nhhp[i] = 0;
    nhpp1[i] = 0;
    nhhp1[i] = 0;
    nhhh[i] = 0;
    nppp[i] = 0;
  }

  for(int i = 0; i < size3; ++i){
    for(int j = 0; j < size3; ++j){
      plus(state, qnums3[i], qnums3[j]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      for(int p1 = 0; p1 < np[j]; ++p1){ // <p1p2|h1h2> -> <p1|h1h2p2>
	for(int hh1 = 0; hh1 < nhh[ind1]; ++hh1){
	  key = Hash3(hhvec[ind1][2 * hh1], hhvec[ind1][2 * hh1 + 1], pvec[j][p1], Space.indtot);
	  hhpvec[i][3 * nhhp[i]] = hhvec[ind1][2 * hh1];
	  hhpvec[i][3 * nhhp[i] + 1] = hhvec[ind1][2 * hh1 + 1];
	  hhpvec[i][3 * nhhp[i] + 2] = pvec[j][p1];
	  hhp_map[i][key] = nhhp[i];
	  ++nhhp[i];
	}
      }
      for(int h1 = 0; h1 < nh[j]; ++h1){ // <h1h2|p1p2> -> <h1|h2p1p2>
	for(int pp1 = 0; pp1 < npp[ind1]; ++pp1){
	  key = Hash3(hvec[j][h1], ppvec[ind1][2 * pp1], ppvec[ind1][2 * pp1 + 1], Space.indtot);
	  hppvec[i][3 * nhpp[i]] = hvec[j][h1];
	  hppvec[i][3 * nhpp[i] + 1] = ppvec[ind1][2 * pp1];
	  hppvec[i][3 * nhpp[i] + 2] = ppvec[ind1][2 * pp1 + 1];
	  hpp_map[i][key] = nhpp[i];
	  ++nhpp[i];
	}
      }
      for(int h1 = 0; h1 < nh[j]; ++h1){ // <h1h2|h3p1> -> <h1|h2h3p1>
	for(int hp1 = 0; hp1 < nhp[ind1]; ++hp1){
	  key = Hash3(hvec[j][h1], hpvec[ind1][2 * hp1], hpvec[ind1][2 * hp1 + 1], Space.indtot);
	  hhp1vec[i][3 * nhhp1[i]] = hvec[j][h1];
	  hhp1vec[i][3 * nhhp1[i] + 1] = hpvec[ind1][2 * hp1];
	  hhp1vec[i][3 * nhhp1[i] + 2] = hpvec[ind1][2 * hp1 + 1];
	  hhp1_map[i][key] = nhhp1[i];
	  ++nhhp1[i];
	}
      }
      for(int p1 = 0; p1 < np[j]; ++p1){ // <p1p2|h1p3> -> <p1|h1p3p2>
	for(int hp1 = 0; hp1 < nhp[ind1]; ++hp1){
	  key = Hash3(hpvec[ind1][2 * hp1], hpvec[ind1][2 * hp1 + 1], pvec[j][p1], Space.indtot);
	  hpp1vec[i][3 * nhpp1[i]] = hpvec[ind1][2 * hp1];
	  hpp1vec[i][3 * nhpp1[i] + 1] = hpvec[ind1][2 * hp1 + 1];
	  hpp1vec[i][3 * nhpp1[i] + 2] = pvec[j][p1];
	  hpp1_map[i][key] = nhpp1[i];
	  ++nhpp1[i];
	}
      }
      for(int h1 = 0; h1 < nh[j]; ++h1){ // <p1h1|h2h3> -> <p1|h1h2h3>
	for(int hh1 = 0; hh1 < nhh[ind1]; ++hh1){
	  key = Hash3(hvec[j][h1], hhvec[ind1][2 * hh1], hhvec[ind1][2 * hh1 + 1], Space.indtot);
	  hhhvec[i][3 * nhhh[i]] = hvec[j][h1];
	  hhhvec[i][3 * nhhh[i] + 1] = hhvec[ind1][2 * hh1];
	  hhhvec[i][3 * nhhh[i] + 2] = hhvec[ind1][2 * hh1 + 1];
	  hhh_map[i][key] = nhhh[i];
	  ++nhhh[i];
	}
      }
      for(int p1 = 0; p1 < np[j]; ++p1){ // <h1p1|p2p3> -> <h1|p1p2p3>
	for(int pp1 = 0; pp1 < npp[ind1]; ++pp1){
	  key = Hash3(pvec[j][p1], ppvec[ind1][2 * pp1], ppvec[ind1][2 * pp1 + 1], Space.indtot);
	  pppvec[i][3 * nppp[i]] = pvec[j][p1];
	  pppvec[i][3 * nppp[i] + 1] = ppvec[ind1][2 * pp1];
	  pppvec[i][3 * nppp[i] + 2] = ppvec[ind1][2 * pp1 + 1];
	  ppp_map[i][key] = nppp[i];
	  ++nppp[i];
	}
      }
    }
  }

  double memory = 0.0;
  int intsize = sizeof(int);
  int doubsize = sizeof(double);
  for(int i = 0; i < size1; ++i){
    memory += (16 * intsize + 3 * doubsize) * nhh[i] * npp[i]; // Tmap, Evec, T1, V4
    memory += 2 * doubsize * nhh[i] * nhh[i]; // S1, V2
    memory += doubsize * npp[i] * npp[i]; // V1
  }
  for(int i = 0; i < size3; ++i){
    memory += 4 * doubsize * nhpp[i] * nh[i]; // T6, T7, V5, V6
    memory += 4 * doubsize * nhhp[i] * np[i]; // T8, T9, V7, V8
    memory += 2 * doubsize * nh[i] * nh[i]; // S2, S3
    memory += 2 * doubsize * np[i] * np[i]; // S4, S5
  }
  for(int i = 0; i < size2; ++i){
    memory += 6 * doubsize * nhp1[i] * nhp2[i]; // T2, T3, T4, T5, V9, V10
    memory += 3 * doubsize * nhp2[i] * nhp2[i]; // S6, S7, V3
  }

  if(Parameters.approx == "singles"){
    memory += (3 * intsize + 3 * doubsize) * nhp1[ind0]; // Tmap, Evec, T1, S4
    memory += (18 * intsize + doubsize) * nhp1[ind0] * nhp1[ind0]; // Tmap2, S3
    memory += (doubsize + 2 * intsize) * nhh1[ind0]; // Q31, Qmap3
    memory += (doubsize + 2 * intsize) * npp1[ind0]; // Q41, Qmap4
    for(int i = 0; i < size1; ++i){
      memory += (doubsize + 2 * intsize) * nhp[i] * npp[i]; // Q11, Qmap1
      memory += (doubsize + 2 * intsize) * nhh[i] * nhp[i]; // Q21, Qmap2
      memory += doubsize * nhh[i] * npp[i]; // E1
      memory += (2 * doubsize + 2 * intsize) * nhh[i] * nhp[i]; // Q61, Qmap6, V19
      memory += (2 * doubsize + 2 * intsize) * nhp[i] * npp[i]; // Q51, Qmap5, V20
    }
    for(int i = 0; i < size3; ++i){
      memory += 3 * doubsize * nhpp[i] * np[i]; // Q12, V11, V17
      memory += 3 * doubsize * nhhp[i] * nh[i]; // Q22, V12, V18
      memory += 2 * doubsize * nh[i] * np[i]; // T2, T3
      memory += 2 * doubsize * nhpp[i] * nh[i]; // E6, E7
      memory += 2 * doubsize * nhhp[i] * np[i]; // E8, E9
      memory += doubsize * nhpp1[i] * nh[i]; // Q11
      memory += doubsize * nhhp1[i] * np[i]; // Q21
      memory += 2 * doubsize * nh[i] * nh[i]; // Q32, S2
      memory += 2 * doubsize * np[i] * np[i]; // Q42, S1
      memory += doubsize * nhhp[i] * nh[i]; // Q62
      memory += doubsize * nhpp[i] * np[i]; // Q52
      memory += 2 * intsize * nhpp1[i] * nh[i]; // Qmap1
      memory += 2 * intsize * nhhp1[i] * np[i]; // Qmap2
      memory += doubsize * nhpp1[i] * np[i]; // V13
      memory += doubsize * nhhp1[i] * nh[i]; // V14
    }
    for(int i = 0; i < size2; ++i){
      memory += 4 * doubsize * nhp1[i] * nhp2[i]; // E2, E3, E4, E5
      memory += 2 * doubsize * nhp1[i] * nhp1[i]; // Q12, Q22
      memory += doubsize * npp1[i] * nhp1[i]; // V16
      memory += doubsize * nhh1[i] * nhp1[i]; // V15
    }
  }

  if(Parameters.extra == -1 || Parameters.extra == 0 || Parameters.extra == 1){
    memory += (2 * doubsize + 6 * intsize) * nhp1[ind0]; // X_ia1, Map_ia, X_ai1, Map_ai
    memory += (doubsize + 3 * intsize) * npp1[ind0]; // X_ab1, Map_ab
    memory += (2 * doubsize + 3 * intsize) * nhh1[ind0]; // X_ij1, X1_ij1, Map_ij
    memory += (doubsize + 3 * intsize) * nhp1[ind0]; // X_ai1, Map_ai
    for(int i = 0; i < size3; ++i){
      memory += 4 * doubsize * nh[i] * np[i]; // X_ia2, X_ia3, X_ai2, X_ai3
      memory += 2 * doubsize * np[i] * np[i]; // X_ab2, X_ab3
      memory += 4 * doubsize * nh[i] * nh[i]; // X_ij2, X_ij3, X1_ij2, X1_ij3
      memory += (3 * doubsize + 20 * intsize) * np[i] * nhpp[i]; // X1_iabc1, X_iabc1, Map_iabc, X_abic1, Map_abic
      memory += (4 * doubsize + 18 * intsize) * nh[i] * nhhp[i]; // X1_ijka1, X_ijka1, Map_ijka, X2_iajk1, X_iajk1, Map_iajk
      memory += 2 * doubsize * nh[i] * nppp[i]; // X1_iabc2, X_abic2
      memory += 2 * doubsize * np[i] * nhhh[i]; // X1_ijka2, X2_iajk2
      memory += 4 * doubsize * np[i] * nhpp1[i]; // X1_iabc3, X_iabc3, X_abic3, X_abic4
      memory += 2 * doubsize * nh[i] * nhhp1[i]; // X2_iajk3, X2_iajk4
      memory += 3 * doubsize * np[i] * nppp[i]; // X1_abcd2, X1_abcd3, V_abcd
      memory += 3 * doubsize * nh[i] * nhhh[i]; // X_ijkl2, X_ijkl3, V_ijkl
      memory += 4 * doubsize * nh[i] * nhpp1[i]; // X1_iajb3, X1_iajb4, X3_iajb3, X_iajb3
      memory += 3 * doubsize * np[i] * nhhp1[i]; // X1_iajb2, X3_iajb2, X3_iajb5
    }
    for(int i = 0; i < size1; ++i){
      memory += doubsize * nhh[i] * npp[i]; // X_ijab1
      memory += (2 * doubsize + 8 * intsize) * nhh[i] * nhh[i]; // X_ijkl1, X_ijkl4, Map_ijkl
      memory += (2 * doubsize + 6 * intsize) * npp[i] * npp[i]; // X1_abcd1, X_abcd1, Map_abcd
      memory += 2 * doubsize * npp[i] * nhp[i]; // X_iabc5, X_abic7
      memory += 2 * doubsize * nhh[i] * nhp[i]; // X_ijka5, X2_iajk7
    }    
    for(int i = 0; i < size2; ++i){
      memory += doubsize * npp1[i] * nhp1[i]; // X_iabc4
      memory += doubsize * nhh1[i] * nhp1[i]; // X_ijka4
      memory += (3 * doubsize + 8 * intsize) * nhp2[i] * nhp2[i]; // X1_iajb1, X3_iajb1, X_iajb1, Map_iajb
      memory += 2 * doubsize * npp1[i] * nhp2[i]; // X_abic5, X_abic6
      memory += 2 * doubsize * nhh1[i] * nhp2[i]; // X2_iajk5, X2_iajk6
    }
  }
  std::cout << std::endl << "Estimated Memory = " << memory/1000000.0 << " MB" << std::endl << std::endl;
}
  
void Channels::delete_struct()
{
  int h, p, hh, pp, hp, hp1, hp2, hh1, pp1, hpp, hhp, hpp1, hhp1, hhh, ppp;
  for(int i = 0; i < size3; ++i){
    h = nh[i];
    p = np[i];
    if(h != 0){ delete[] hvec[i]; }
    if(p != 0){ delete[] pvec[i]; }
  }
  delete[] indvec;
  delete[] qnums3;
  delete[] hvec;
  delete[] pvec;
  delete[] h_map;
  delete[] p_map;
  delete[] nh;
  delete[] np;

  delete[] qnums1;
  delete[] qnums2;

  for(int i = 0; i < size1; ++i){
    hh = nhh[i];
    pp = npp[i];
    hp = nhp[i];
    if(hh != 0){ delete[] hhvec[i]; }
    if(pp != 0){ delete[] ppvec[i]; }
    if(hp != 0){ delete[] hpvec[i]; }
  }
  delete[] hhvec;
  delete[] ppvec;
  delete[] hpvec;
  delete[] hh_map;
  delete[] pp_map;
  delete[] hp_map;
  delete[] nhh;
  delete[] npp;
  delete[] nhp;

  for(int i = 0; i < size2; ++i){
    hp1 = nhp1[i];
    hp2 = nhp2[i];
    hh1 = nhh1[i];
    pp1 = npp1[i];
    if(hp1 != 0){ delete[] hp1vec[i]; }
    if(hp2 != 0){ delete[] hp2vec[i]; }
    if(hh1 != 0){ delete[] hh1vec[i]; }
    if(pp1 != 0){ delete[] pp1vec[i]; }
  }
  delete[] hp1vec;
  delete[] hp2vec;
  delete[] hh1vec;
  delete[] pp1vec;
  delete[] hp1_map;
  delete[] hp2_map;
  delete[] hh1_map;
  delete[] pp1_map;
  delete[] nhp1;
  delete[] nhp2;
  delete[] nhh1;
  delete[] npp1;

  for(int i = 0; i < size3; ++i){
    hpp = nhpp[i];
    hhp = nhhp[i];
    hpp1 = nhpp1[i];
    hhp1 = nhhp1[i];
    hhh = nhhh[i];
    ppp = nppp[i];
    if(hpp != 0){ delete[] hppvec[i]; }
    if(hhp != 0){ delete[] hhpvec[i]; }
    if(hpp1 != 0){ delete[] hpp1vec[i]; }
    if(hhp1 != 0){ delete[] hhp1vec[i]; }
    if(hhh != 0){ delete[] hhhvec[i]; }
    if(ppp != 0){ delete[] pppvec[i]; }
  }
  delete[] hppvec;
  delete[] hhpvec;
  delete[] hpp1vec;
  delete[] hhp1vec;
  delete[] hhhvec;
  delete[] pppvec;
  delete[] hpp_map;
  delete[] hhp_map;
  delete[] hpp1_map;
  delete[] hhp1_map;
  delete[] hhh_map;
  delete[] ppp_map;
  delete[] nhpp;
  delete[] nhhp;
  delete[] nhpp1;
  delete[] nhhp1;
  delete[] nhhh;
  delete[] nppp;
}


CC_Eff::CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  std::cout << "Building Effective Hamiltonian ..." << std::endl;
  int length0, length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1, nhhh, nppp;

  nhp1 = Chan.nhp1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  nhp1 = Chan.nhp1[Chan.ind0];
  length = nhp1;
  if(length != 0){
    X_ia1 = new double[length]; // (ia)
    Map_ia = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ia1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ia[3*j + k] = -1; }
    }
  }
  length = npp1;
  if(length != 0){
    X_ab1 = new double[length]; // (ba)
    Map_ab = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ab1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ab[3*j + k] = -1; }
    }
  }
  length = nhh1;
  if(nhh1 != 0){
    X_ij1 = new double[length]; // (ji)
    X1_ij1 = new double[length]; // (ji)
    Map_ij = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ij1[j] = 0.0;
      X1_ij1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ij[3*j + k] = -1; }
    }
  }
  length = nhp1;
  if(nhp1 != 0){
    X_ai1 = new double[length]; // (ia)
    Map_ai = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ai1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ai[3*j + k] = -1; }
    }
  }

  X_ia2 = new double*[Chan.size3];
  X_ia3 = new double*[Chan.size3];

  X_ab2 = new double*[Chan.size3];
  X_ab3 = new double*[Chan.size3];

  X_ij2 = new double*[Chan.size3];
  X_ij3 = new double*[Chan.size3];
  X1_ij2 = new double*[Chan.size3];
  X1_ij3 = new double*[Chan.size3];

  X_ai2 = new double*[Chan.size3];
  X_ai3 = new double*[Chan.size3];

  X_ijab1 = new double*[Chan.size1];

  X1_iabc1 = new double*[Chan.size3];
  X1_iabc2 = new double*[Chan.size3];
  X1_iabc3 = new double*[Chan.size3];
  X_iabc1 = new double*[Chan.size3];
  X_iabc3 = new double*[Chan.size3];
  X_iabc4 = new double*[Chan.size2];
  X_iabc5 = new double*[Chan.size1];
  Map_iabc = new int*[Chan.size3];

  X1_ijka1 = new double*[Chan.size3];
  X1_ijka2 = new double*[Chan.size3];
  X_ijka1 = new double*[Chan.size3];
  X_ijka4 = new double*[Chan.size2];
  X_ijka5 = new double*[Chan.size1];
  Map_ijka = new int*[Chan.size3];

  X1_abcd1 = new double*[Chan.size1];
  X1_abcd2 = new double*[Chan.size3];
  X1_abcd3 = new double*[Chan.size3];
  X_abcd1 = new double*[Chan.size1];
  V_abcd = new double*[Chan.size3];
  Map_abcd = new int*[Chan.size1];

  X_ijkl1 = new double*[Chan.size1];
  X_ijkl2 = new double*[Chan.size3];
  X_ijkl3 = new double*[Chan.size3];
  X_ijkl4 = new double*[Chan.size1];
  V_ijkl = new double*[Chan.size3];
  Map_ijkl = new int*[Chan.size1];

  X1_iajb1 = new double*[Chan.size2];
  X1_iajb2 = new double*[Chan.size3];
  X1_iajb3 = new double*[Chan.size3];
  X1_iajb4 = new double*[Chan.size3];
  X3_iajb1 = new double*[Chan.size2];
  X3_iajb2 = new double*[Chan.size3];
  X3_iajb3 = new double*[Chan.size3];
  X3_iajb5 = new double*[Chan.size3];
  X_iajb1 = new double*[Chan.size2];
  X_iajb3 = new double*[Chan.size3];
  Map_iajb = new int*[Chan.size2];

  X_abic1 = new double*[Chan.size3];
  X_abic2 = new double*[Chan.size3];
  X_abic3 = new double*[Chan.size3];
  X_abic4 = new double*[Chan.size3];
  X_abic5 = new double*[Chan.size2];
  X_abic6 = new double*[Chan.size2];
  X_abic7 = new double*[Chan.size1];
  Map_abic = new int*[Chan.size3];

  X2_iajk1 = new double*[Chan.size3];
  X2_iajk2 = new double*[Chan.size3];
  X2_iajk3 = new double*[Chan.size3];
  X2_iajk4 = new double*[Chan.size3];
  X2_iajk5 = new double*[Chan.size2];
  X2_iajk6 = new double*[Chan.size2];
  X2_iajk7 = new double*[Chan.size1];
  X_iajk1 = new double*[Chan.size3];
  Map_iajk = new int*[Chan.size3];

  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    nhhh = Chan.nhhh[i];
    nppp = Chan.nppp[i];
    length = nh * np;
    if(length != 0){
      X_ia2[i] = new double[length]; // (i,a)
      X_ia3[i] = new double[length]; // (a,i)
      X_ai2[i] = new double[length]; // (a,i)
      X_ai3[i] = new double[length]; // (i,a)
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ia2[i][j] = 0.0;
	X_ia3[i][j] = 0.0;
	X_ai2[i][j] = 0.0;
	X_ai3[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      X_ab2[i] = new double[length]; // (a,b)
      X_ab3[i] = new double[length]; // (b,a)
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ab2[i][j] = 0.0;
	X_ab3[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      X_ij2[i] = new double[length]; // (i,j)
      X_ij3[i] = new double[length]; // (j,i)
      X1_ij2[i] = new double[length]; // (i,j)
      X1_ij3[i] = new double[length]; // (j,i)
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ij2[i][j] = 0.0;
	X_ij3[i][j] = 0.0;
	X1_ij2[i][j] = 0.0;
	X1_ij3[i][j] = 0.0;
      }
    }
    length = np * nhpp;
    if(length != 0){
      X1_iabc1[i] = new double[length];
      X_iabc1[i] = new double[length];
      Map_iabc[i] = new int[8 * length];
      X_abic1[i] = new double[length];
      Map_abic[i] = new int[12 * length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_iabc1[i][j] = 0.0;
	X_iabc1[i][j] = 0.0;
	X_abic1[i][j] = 0.0;
	for(int k = 0; k < 8; ++k){ Map_iabc[i][8*j + k] = -1; }
	for(int k = 0; k < 12; ++k){ Map_abic[i][12*j + k] = -1; }
      }
    }
    length = nh * nhhp;
    if(length != 0){
      X1_ijka1[i] = new double[length];
      X_ijka1[i] = new double[length];
      Map_ijka[i] = new int[6 * length];
      X2_iajk1[i] = new double[length];
      X_iajk1[i] = new double[length];
      Map_iajk[i] = new int[12 * length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_ijka1[i][j] = 0.0;
	X_ijka1[i][j] = 0.0;
	X2_iajk1[i][j] = 0.0;
	X_iajk1[i][j] = 0.0;
	for(int k = 0; k < 6; ++k){ Map_ijka[i][6*j + k] = -1; }
	for(int k = 0; k < 12; ++k){ Map_iajk[i][12*j + k] = -1; }
      }
    }
    length = nh * nppp;
    if(length != 0){
      X1_iabc2[i] = new double[length];
      X_abic2[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_iabc2[i][j] = 0.0;
	X_abic2[i][j] = 0.0;
      }
    }
    length = np * nhhh;
    if(length != 0){
      X1_ijka2[i] = new double[length];
      X2_iajk2[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_ijka2[i][j] = 0.0;
	X2_iajk2[i][j] = 0.0;
      }
    }
    length = np * nhpp1;
    if(length != 0){
      X1_iabc3[i] = new double[length];
      X_iabc3[i] = new double[length];
      X_abic3[i] = new double[length];
      X_abic4[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_iabc3[i][j] = 0.0;      
	X_iabc3[i][j] = 0.0;
	X_abic3[i][j] = 0.0;
	X_abic4[i][j] = 0.0;
      }
    }
    length = nh * nhhp1;
    if(length != 0){
      X2_iajk3[i] = new double[length];
      X2_iajk4[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X2_iajk3[i][j] = 0.0;
	X2_iajk4[i][j] = 0.0;
      }
    }
    length = np * nppp;
    if(length != 0){
      X1_abcd2[i] = new double[length];
      X1_abcd3[i] = new double[length];
      V_abcd[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_abcd2[i][j] = 0.0;
	X1_abcd3[i][j] = 0.0;
	V_abcd[i][j] = 0.0;
      }
    }
    length = nh * nhhh;
    if(length != 0){
      X_ijkl2[i] = new double[length];
      X_ijkl3[i] = new double[length];
      V_ijkl[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ijkl2[i][j] = 0.0;
	X_ijkl3[i][j] = 0.0;
	V_ijkl[i][j] = 0.0;
      }
    }
    length = nh * nhpp1;
    if(length != 0){
      X1_iajb3[i] = new double[length];
      X1_iajb4[i] = new double[length];
      X3_iajb3[i] = new double[length];
      X_iajb3[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_iajb3[i][j] = 0.0;
	X1_iajb4[i][j] = 0.0;
	X3_iajb3[i][j] = 0.0;
	X_iajb3[i][j] = 0.0;
      }
    }
    length = np * nhhp1;
    if(length != 0){
      X1_iajb2[i] = new double[length];
      X3_iajb2[i] = new double[length];
      X3_iajb5[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_iajb2[i][j] = 0.0;
	X3_iajb2[i][j] = 0.0;
	X3_iajb5[i][j] = 0.0;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      X_ijab1[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ijab1[i][j] = 0.0;
      }
    }
    length = nhh * nhh;
    if(length != 0){
      X_ijkl1[i] = new double[length];
      X_ijkl4[i] = new double[length];
      Map_ijkl[i] = new int[8 * length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ijkl1[i][j] = 0.0;
	X_ijkl4[i][j] = 0.0;
	for(int k = 0; k < 8; ++k){ Map_ijkl[i][8*j + k] = -1; }
      }
    }
    length = npp * npp;
    if(length != 0){
      X1_abcd1[i] = new double[length];
      X_abcd1[i] = new double[length];
      Map_abcd[i] = new int[6 * length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_abcd1[i][j] = 0.0;
	X_abcd1[i][j] = 0.0;
	for(int k = 0; k < 6; ++k){ Map_abcd[i][6*j + k] = -1; }
      }
    }
    length = npp * nhp;
    if(length != 0){
      X_iabc5[i] = new double[length];
      X_abic7[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_iabc5[i][j] = 0.0;
	X_abic7[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      X_ijka5[i] = new double[length];
      X2_iajk7[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ijka5[i][j] = 0.0;
	X2_iajk7[i][j] = 0.0;
      }
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    length = npp1 * nhp1;
    if(length != 0){
      X_iabc4[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_iabc4[i][j] = 0.0;
      }
    }
    length = nhh1 * nhp1;
    if(length != 0){
      X_ijka4[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_ijka4[i][j] = 0.0;
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      X1_iajb1[i] = new double[length];
      X3_iajb1[i] = new double[length];
      X_iajb1[i] = new double[length];
      Map_iajb[i] = new int[8 * length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X1_iajb1[i][j] = 0.0;
	X3_iajb1[i][j] = 0.0;
	X_iajb1[i][j] = 0.0;
	for(int k = 0; k < 8; ++k){ Map_iajb[i][8*j + k] = -1; }
      }
    }
    length = npp1 * nhp2;
    if(length != 0){
      X_abic5[i] = new double[length];
      X_abic6[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X_abic5[i][j] = 0.0;
	X_abic6[i][j] = 0.0;
      }
    }
    length = nhh1 * nhp2;
    if(length != 0){
      X2_iajk5[i] = new double[length];
      X2_iajk6[i] = new double[length];
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < length; ++j){
	X2_iajk5[i][j] = 0.0;
	X2_iajk6[i][j] = 0.0;
      }
    }
  }
 
  #pragma omp parallel
  {
    int ind, ind1, key1, key2;
    int hind1, hind2, pind1, pind2;
    length = Chan.nhp1[Chan.ind0];
    #pragma omp for schedule(static)
    for(int i = 0; i < length; ++i){
      hind1 = Chan.hp1vec[Chan.ind0][2*i];
      pind1 = Chan.hp1vec[Chan.ind0][2*i + 1];
      ind = Chan.indvec[hind1];
      key1 = Chan.h_map[ind][hind1];
      key2 = Chan.p_map[ind][pind1];
      Map_ia[3*i] = ind;
      ind1 = key1 * Chan.np[ind] + key2;
      //ind1 = Index11(Chan.hvec[ind], Chan.pvec[ind], Chan.nh[ind], Chan.np[ind], hind1, pind1);
      Map_ia[3*i + 1] = ind1;
      ind1 = key2 * Chan.nh[ind] + key1;
      //ind1 = Index11(Chan.pvec[ind], Chan.hvec[ind], Chan.np[ind], Chan.nh[ind], pind1, hind1);
      Map_ia[3*i + 2] = ind1;
    }
    length = Chan.npp1[Chan.ind0];
    #pragma omp for schedule(static)
    for(int i = 0; i < length; ++i){
      pind2 = Chan.pp1vec[Chan.ind0][2*i];
      pind1 = Chan.pp1vec[Chan.ind0][2*i + 1];
      ind = Chan.indvec[pind1];
      key1 = Chan.p_map[ind][pind1];
      key2 = Chan.p_map[ind][pind2];
      Map_ab[3*i] = ind;
      ind1 = key1 * Chan.np[ind] + key2;
      //ind1 = Index11(Chan.pvec[ind], Chan.pvec[ind], Chan.np[ind], Chan.np[ind], pind1, pind2);
      Map_ab[3*i + 1] = ind1;
      ind1 = key2 * Chan.np[ind] + key1;
      //ind1 = Index11(Chan.pvec[ind], Chan.pvec[ind], Chan.np[ind], Chan.np[ind], pind2, pind1);
      Map_ab[3*i + 2] = ind1;
    }
    length = Chan.nhh1[Chan.ind0];
    #pragma omp for schedule(static)
    for(int i = 0; i < length; ++i){
      hind2 = Chan.hh1vec[Chan.ind0][2*i];
      hind1 = Chan.hh1vec[Chan.ind0][2*i + 1];
      ind = Chan.indvec[hind1];
      key1 = Chan.h_map[ind][hind1];
      key2 = Chan.h_map[ind][hind2];
      Map_ij[3*i] = ind;
      ind1 = key1 * Chan.nh[ind] + key2;
      //ind1 = Index11(Chan.hvec[ind], Chan.hvec[ind], Chan.nh[ind], Chan.nh[ind], hind1, hind2);
      Map_ij[3*i + 1] = ind1;
      ind1 = key2 * Chan.nh[ind] + key1;
      //ind1 = Index11(Chan.hvec[ind], Chan.hvec[ind], Chan.nh[ind], Chan.nh[ind], hind2, hind1);
      Map_ij[3*i + 2] = ind1;
    }
    length = Chan.nhp1[Chan.ind0];
    #pragma omp for schedule(static)
    for(int i = 0; i < length; ++i){
      hind1 = Chan.hp1vec[Chan.ind0][2*i];
      pind1 = Chan.hp1vec[Chan.ind0][2*i + 1];
      ind = Chan.indvec[hind1];
      key1 = Chan.p_map[ind][pind1];
      key2 = Chan.h_map[ind][hind1];
      Map_ai[3*i] = ind;
      ind1 = key1 * Chan.nh[ind] + key2;
      //ind1 = Index11(Chan.pvec[ind], Chan.hvec[ind], Chan.np[ind], Chan.nh[ind], pind1, hind1);
      Map_ai[3*i + 1] = ind1;
      ind1 = key2 * Chan.np[ind] + key1;
      //ind1 = Index11(Chan.hvec[ind], Chan.pvec[ind], Chan.nh[ind], Chan.np[ind], hind1, pind1);
      Map_ai[3*i + 2] = ind1;
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, pind1, pind2, pind3;
      int p1, hpp1;
      length = Chan.np[i] * Chan.nhpp[i];
      length0 = Chan.nhpp[i];
      #pragma omp for schedule(static)
      for(int phpp = 0; phpp < length; ++phpp){
	hpp1 = int(phpp%length0);
	p1 = int((phpp - hpp1)/length0);
	pind1 = Chan.pvec[i][p1];
	hind1 = Chan.hppvec[i][3*hpp1];
	pind2 = Chan.hppvec[i][3*hpp1 + 1];
	pind3 = Chan.hppvec[i][3*hpp1 + 2];

	ind = Chan.indvec[hind1];
	Map_iabc[i][8*phpp] = ind;
	key1 = Chan.ppp_map[ind][Hash3(pind1, pind2, pind3, Space.indtot)];
	key2 = Chan.h_map[ind][hind1];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.pppvec[ind], Chan.hvec[ind], Chan.nppp[ind], Chan.nh[ind], pind1, pind2, pind3, hind1);
	Map_iabc[i][8*phpp + 1] = ind1;
	ind = Chan.indvec[pind2];
	Map_iabc[i][8*phpp + 2] = ind;
	key1 = Chan.hpp1_map[ind][Hash3(hind1, pind1, pind3, Space.indtot)];
	key2 = Chan.p_map[ind][pind2];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hpp1vec[ind], Chan.pvec[ind], Chan.nhpp1[ind], Chan.np[ind], hind1, pind1, pind3, pind2);
	Map_iabc[i][8*phpp + 3] = ind1;
	minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iabc[i][8*phpp + 4] = ind;
	key1 = Chan.pp1_map[ind][Hash2(pind3, pind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	//ind1 = Index22(Chan.pp1vec[ind], Chan.hp1vec[ind], Chan.npp1[ind], Chan.nhp1[ind], pind3, pind1, hind1, pind2);
	Map_iabc[i][8*phpp + 5] = ind1;
	plus(tb, Space.qnums[pind2], Space.qnums[pind3]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_iabc[i][8*phpp + 6] = ind;
	key1 = Chan.pp_map[ind][Hash2(pind2, pind3, Space.indtot)];
	key2 = Chan.hp_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp[ind] + key2;
	//ind1 = Index22(Chan.ppvec[ind], Chan.hpvec[ind], Chan.npp[ind], Chan.nhp[ind], pind2, pind3, hind1, pind1);
	Map_iabc[i][8*phpp + 7] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, hind3, pind1;
      int h1, hhp1;
      length = Chan.nh[i] * Chan.nhhp[i];
      length0 = Chan.nhhp[i];
      #pragma omp for schedule(static)
      for(int hhhp = 0; hhhp < length; ++hhhp){
	hhp1 = int(hhhp%length0);
	h1 = int((hhhp - hhp1)/length0);
	hind3 = Chan.hvec[i][h1];
	hind1 = Chan.hhpvec[i][3*hhp1];
	hind2 = Chan.hhpvec[i][3*hhp1 + 1];
	pind1 = Chan.hhpvec[i][3*hhp1 + 2];
	
	ind = Chan.indvec[pind1];
	Map_ijka[i][6*hhhp] = ind;
	key1 = Chan.hhh_map[ind][Hash3(hind3, hind1, hind2, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hhhvec[ind], Chan.pvec[ind], Chan.nhhh[ind], Chan.np[ind], hind3, hind1, hind2, pind1);
	Map_ijka[i][6*hhhp + 1] = ind1;
	minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_ijka[i][6*hhhp + 2] = ind;
	key1 = Chan.hh1_map[ind][Hash2(hind3, hind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	//ind1 = Index22(Chan.hh1vec[ind], Chan.hp1vec[ind], Chan.nhh1[ind], Chan.nhp1[ind], hind3, hind1, hind2, pind1);
	Map_ijka[i][6*hhhp + 3] = ind1;
	plus(tb, Space.qnums[hind1], Space.qnums[hind2]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_ijka[i][6*hhhp + 4] = ind;
	key1 = Chan.hp_map[ind][Hash2(hind3, pind1, Space.indtot)];
	key2 = Chan.hh_map[ind][Hash2(hind1, hind2, Space.indtot)];
	ind1 = key1 * Chan.nhh[ind] + key2;
	//ind1 = Index22(Chan.hpvec[ind], Chan.hhvec[ind], Chan.nhp[ind], Chan.nhh[ind], hind3, pind1, hind1, hind2);
	Map_ijka[i][6*hhhp + 5] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int pind1, pind2, pind3, pind4;
      int pp1, pp2;
      length = Chan.npp[i] * Chan.npp[i];
      length0 = Chan.npp[i];
      #pragma omp for schedule(static)
      for(int pppp = 0; pppp < length; ++pppp){
	pp2 = int(pppp%length0);
	pp1 = int((pppp - pp2)/length0);
	pind3 = Chan.ppvec[i][2*pp1];
	pind4 = Chan.ppvec[i][2*pp1 + 1];
	pind1 = Chan.ppvec[i][2*pp2];
	pind2 = Chan.ppvec[i][2*pp2 + 1];

	ind = Chan.indvec[pind2];
	Map_abcd[i][6*pppp] = ind;
	key1 = Chan.ppp_map[ind][Hash3(pind1, pind3, pind4, Space.indtot)];
	key2 = Chan.p_map[ind][pind2];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.pppvec[ind], Chan.pvec[ind], Chan.nppp[ind], Chan.np[ind], pind1, pind3, pind4, pind2);
	Map_abcd[i][6*pppp + 1] = ind1;
	ind = Chan.indvec[pind1];
	Map_abcd[i][6*pppp + 2] = ind;
	key1 = Chan.ppp_map[ind][Hash3(pind2, pind3, pind4, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.pppvec[ind], Chan.pvec[ind], Chan.nppp[ind], Chan.np[ind], pind2, pind3, pind4, pind1);
	Map_abcd[i][6*pppp + 3] = ind1;
	ind = Chan.indvec[pind3];
	Map_abcd[i][6*pppp + 4] = ind;
	key1 = Chan.ppp_map[ind][Hash3(pind4, pind1, pind2, Space.indtot)];
	key2 = Chan.p_map[ind][pind3];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.pppvec[ind], Chan.pvec[ind], Chan.nppp[ind], Chan.np[ind], pind4, pind1, pind2, pind3);
	Map_abcd[i][6*pppp + 5] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, hind3, hind4;
      int hh1, hh2;
      length = Chan.nhh[i] * Chan.nhh[i];
      length0 = Chan.nhh[i];
      #pragma omp for schedule(static)
      for(int hhhh = 0; hhhh < length; ++hhhh){
	hh2 = int(hhhh%length0);
	hh1 = int((hhhh - hh2)/length0);
	hind1 = Chan.hhvec[i][2*hh1];
	hind2 = Chan.hhvec[i][2*hh1 + 1];
	hind3 = Chan.hhvec[i][2*hh2];
	hind4 = Chan.hhvec[i][2*hh2 + 1];

	ind = Chan.indvec[hind4];
	Map_ijkl[i][8*hhhh] = ind;
	key1 = Chan.hhh_map[ind][Hash3(hind3, hind1, hind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind4];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hhhvec[ind], Chan.hvec[ind], Chan.nhhh[ind], Chan.nh[ind], hind3, hind1, hind2, hind4);
	Map_ijkl[i][8*hhhh + 1] = ind1;
	ind = Chan.indvec[hind3];
	Map_ijkl[i][8*hhhh + 2] = ind;
	key1 = Chan.hhh_map[ind][Hash3(hind4, hind1, hind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind3];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hhhvec[ind], Chan.hvec[ind], Chan.nhhh[ind], Chan.nh[ind], hind4, hind1, hind2, hind3);
	Map_ijkl[i][8*hhhh + 3] = ind1;
	ind = i;
	Map_ijkl[i][8*hhhh + 4] = ind;
	Map_ijkl[i][8*hhhh + 5] = Chan.nhh[i]*hh2 + hh1;
	ind = Chan.indvec[hind2];
	Map_ijkl[i][8*hhhh + 6] = ind;
	key1 = Chan.hhh_map[ind][Hash3(hind1, hind3, hind4, Space.indtot)];
	key2 = Chan.h_map[ind][hind2];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hhhvec[ind], Chan.hvec[ind], Chan.nhhh[ind], Chan.nh[ind], hind1, hind3, hind4, hind2);
	Map_ijkl[i][8*hhhh + 7] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, pind1, pind2;
      int hp1, hp2;
      length = Chan.nhp2[i] * Chan.nhp2[i];
      length0 = Chan.nhp2[i];
      #pragma omp for schedule(static)
      for(int hphp = 0; hphp < length; ++hphp){
	hp2 = int(hphp%length0);
	hp1 = int((hphp - hp2)/length0);
	hind1 = Chan.hp2vec[i][2*hp1];
	pind2 = Chan.hp2vec[i][2*hp1 + 1];
	hind2 = Chan.hp2vec[i][2*hp2];
	pind1 = Chan.hp2vec[i][2*hp2 + 1];

	ind = Chan.indvec[pind1];
	Map_iajb[i][8*hphp] = ind;
	key1 = Chan.hhp1_map[ind][Hash3(hind1, hind2, pind2, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hhp1vec[ind], Chan.pvec[ind], Chan.nhhp1[ind], Chan.np[ind], hind1, hind2, pind2, pind1);
	Map_iajb[i][8*hphp + 1] = ind1;
	ind = Chan.indvec[hind2];
	Map_iajb[i][8*hphp + 2] = ind;
	key1 = Chan.hpp1_map[ind][Hash3(hind1, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind2];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hpp1vec[ind], Chan.hvec[ind], Chan.nhpp1[ind], Chan.nh[ind], hind1, pind1, pind2, hind2);
	Map_iajb[i][8*hphp + 3] = ind1;
	ind = Chan.indvec[hind1];
	Map_iajb[i][8*hphp + 4] = ind;
	key1 = Chan.hpp1_map[ind][Hash3(hind2, pind2, pind1, Space.indtot)];
	key2 = Chan.h_map[ind][hind1];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hpp1vec[ind], Chan.hvec[ind], Chan.nhpp1[ind], Chan.nh[ind], hind2, pind2, pind1, hind1);
	Map_iajb[i][8*hphp + 5] = ind1;
	ind = Chan.indvec[pind2];
	Map_iajb[i][8*hphp + 6] = ind;
	key1 = Chan.hhp1_map[ind][Hash3(hind2, hind1, pind1, Space.indtot)];
	key2 = Chan.p_map[ind][pind2];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hhp1vec[ind], Chan.pvec[ind], Chan.nhhp1[ind], Chan.np[ind], hind2, hind1, pind1, pind2);
	Map_iajb[i][8*hphp + 7] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, pind1, pind2, pind3;
      int hpp1, p1;
      length = Chan.nhpp[i] * Chan.np[i];
      length0 = Chan.np[i];
      #pragma omp for schedule(static)
      for(int hppp = 0; hppp < length; ++hppp){
	p1 = int(hppp%length0);
	hpp1 = int((hppp - p1)/length0);
	hind1 = Chan.hppvec[i][3*hpp1];
	pind1 = Chan.hppvec[i][3*hpp1 + 1];
	pind2 = Chan.hppvec[i][3*hpp1 + 2];
	pind3 = Chan.pvec[i][p1];
	
	ind = Chan.indvec[hind1];
	Map_abic[i][12*hppp] = ind;
	key1 = Chan.ppp_map[ind][Hash3(pind3, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind1];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.pppvec[ind], Chan.hvec[ind], Chan.nppp[ind], Chan.nh[ind], pind3, pind1, pind2, hind1);
	Map_abic[i][12*hppp + 1] = ind1;
	ind = Chan.indvec[pind1];
	Map_abic[i][12*hppp + 2] = ind;
	key1 = Chan.hpp1_map[ind][Hash3(hind1, pind3, pind2, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hpp1vec[ind], Chan.pvec[ind], Chan.nhpp1[ind], Chan.np[ind], hind1, pind3, pind2, pind1);
	Map_abic[i][12*hppp + 3] = ind1;
	ind = Chan.indvec[pind2];
	Map_abic[i][12*hppp + 4] = ind;
	key1 = Chan.hpp1_map[ind][Hash3(hind1, pind3, pind1, Space.indtot)];
	key2 = Chan.p_map[ind][pind2];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hpp1vec[ind], Chan.pvec[ind], Chan.nhpp1[ind], Chan.np[ind], hind1, pind3, pind1, pind2);
	Map_abic[i][12*hppp + 5] = ind1;
	minus(tb, Space.qnums[pind3], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_abic[i][12*hppp + 6] = ind;
	key1 = Chan.pp1_map[ind][Hash2(pind3, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.pp1vec[ind], Chan.hp2vec[ind], Chan.npp1[ind], Chan.nhp2[ind], pind3, pind2, hind1, pind1);
	Map_abic[i][12*hppp + 7] = ind1;
	minus(tb, Space.qnums[pind3], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_abic[i][12*hppp + 8] = ind;
	key1 = Chan.pp1_map[ind][Hash2(pind3, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.pp1vec[ind], Chan.hp2vec[ind], Chan.npp1[ind], Chan.nhp2[ind], pind3, pind1, hind1, pind2);
	Map_abic[i][12*hppp + 9] = ind1;
	plus(tb, Space.qnums[pind1], Space.qnums[pind2]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_abic[i][12*hppp + 10] = ind;
	key1 = Chan.hp_map[ind][Hash2(hind1, pind3, Space.indtot)];
	key2 = Chan.pp_map[ind][Hash2(pind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.npp[ind] + key2;
	//ind1 = Index22(Chan.hpvec[ind], Chan.ppvec[ind], Chan.nhp[ind], Chan.npp[ind], hind1, pind3, pind1, pind2);
	Map_abic[i][12*hppp + 11] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel
    {
      State tb;
      int ind, ind1, key1, key2;
      int hind1, hind2, hind3, pind1;
      int hhp1, h1;
      length = Chan.nhhp[i] * Chan.nh[i];
      length0 = Chan.nh[i];
      #pragma omp for schedule(static)
      for(int hhph = 0; hhph < length; ++hhph){
	h1 = int(hhph%length0);
	hhp1 = int((hhph - h1)/length0);
	hind2 = Chan.hhpvec[i][3*hhp1];
	hind3 = Chan.hhpvec[i][3*hhp1 + 1];
	pind1 = Chan.hhpvec[i][3*hhp1 + 2];
	hind1 = Chan.hvec[i][h1];
	
	ind = Chan.indvec[pind1];
	Map_iajk[i][12*hhph] = ind;
	key1 = Chan.hhh_map[ind][Hash3(hind1, hind2, hind3, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	//ind1 = Index31(Chan.hhhvec[ind], Chan.pvec[ind], Chan.nhhh[ind], Chan.np[ind], hind1, hind2, hind3, pind1);
	Map_iajk[i][12*hhph + 1] = ind1;
	ind = Chan.indvec[hind3];
	Map_iajk[i][12*hhph + 2] = ind;
	key1 = Chan.hhp1_map[ind][Hash3(hind2, hind1, pind1, Space.indtot)];
	key2 = Chan.h_map[ind][hind3];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hhp1vec[ind], Chan.hvec[ind], Chan.nhhp1[ind], Chan.nh[ind], hind2, hind1, pind1, hind3);
	Map_iajk[i][12*hhph + 3] = ind1;
	ind = Chan.indvec[hind2];
	Map_iajk[i][12*hhph + 4] = ind;
	key1 = Chan.hhp1_map[ind][Hash3(hind3, hind1, pind1, Space.indtot)];
	key2 = Chan.h_map[ind][hind2];
	ind1 = key1 * Chan.nh[ind] + key2;
	//ind1 = Index31(Chan.hhp1vec[ind], Chan.hvec[ind], Chan.nhhp1[ind], Chan.nh[ind], hind3, hind1, pind1, hind2);
	Map_iajk[i][12*hhph + 5] = ind1;
	minus(tb, Space.qnums[hind2], Space.qnums[hind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iajk[i][12*hhph + 6] = ind;
	key1 = Chan.hh1_map[ind][Hash2(hind2, hind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind3, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.hh1vec[ind], Chan.hp2vec[ind], Chan.nhh1[ind], Chan.nhp2[ind], hind2, hind1, hind3, pind1);
	Map_iajk[i][12*hhph + 7] = ind1;
	minus(tb, Space.qnums[hind3], Space.qnums[hind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iajk[i][12*hhph + 8] = ind;
	key1 = Chan.hh1_map[ind][Hash2(hind3, hind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	//ind1 = Index22(Chan.hh1vec[ind], Chan.hp2vec[ind], Chan.nhh1[ind], Chan.nhp2[ind], hind3, hind1, hind2, pind1);
	Map_iajk[i][12*hhph + 9] = ind1;
	plus(tb, Space.qnums[hind2], Space.qnums[hind3]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_iajk[i][12*hhph + 10] = ind;
	key1 = Chan.hh_map[ind][Hash2(hind2, hind3, Space.indtot)];
	key2 = Chan.hp_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp[ind] + key2;
	//ind1 = Index22(Chan.hhvec[ind], Chan.hpvec[ind], Chan.nhh[ind], Chan.nhp[ind], hind2, hind3, hind1, pind1);
	Map_iajk[i][12*hhph + 11] = ind1;
      }
    }
  }
}

void CC_Eff::delete_struct(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1, nhhh, nppp;

  nhp1 = Chan.nhp1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  nhp1 = Chan.nhp1[Chan.ind0];
  length = nhp1;
  if(length != 0){
    delete[] X_ia1; // (ia)
    delete[] Map_ia;
  }
  length = npp1;
  if(length != 0){
    delete[] X_ab1; // (ba)
    delete[] Map_ab;
  }
  length = nhh1;
  if(nhh1 != 0){
    delete[] X_ij1; // (ji)
    delete[] X1_ij1; // (ji)
    delete[] Map_ij;
  }
  length = nhp1;
  if(nhp1 != 0){
    delete[] X_ai1; // (ia)
    delete[] Map_ai;
  }

  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    nhhh = Chan.nhhh[i];
    nppp = Chan.nppp[i];
    length = nh * np;
    if(length != 0){
      delete[] X_ia2[i]; // (i,a)
      delete[] X_ia3[i]; // (a,i)
      delete[] X_ai2[i]; // (a,i)
      delete[] X_ai3[i]; // (i,a)
    }
    length = np * np;
    if(length != 0){
      delete[] X_ab2[i]; // (a,b)
      delete[] X_ab3[i]; // (b,a)
    }
    length = nh * nh;
    if(length != 0){
      delete[] X_ij2[i]; // (i,j)
      delete[] X_ij3[i]; // (j,i)
      delete[] X1_ij2[i]; // (i,j)
      delete[] X1_ij3[i]; // (j,i)
    }
    length = np * nhpp;
    if(length != 0){
      delete[] X1_iabc1[i];
      delete[] X_iabc1[i];
      delete[] Map_iabc[i];
      delete[] X_abic1[i];
      delete[] Map_abic[i];
    }
    length = nh * nhhp;
    if(length != 0){
      delete[] X1_ijka1[i];
      delete[] X_ijka1[i];
      delete[] Map_ijka[i];
      delete[] X2_iajk1[i];
      delete[] X_iajk1[i];
      delete[] Map_iajk[i];
    }
    length = nh * nppp;
    if(length != 0){
      delete[] X1_iabc2[i];
      delete[] X_abic2[i];
    }
    length = np * nhhh;
    if(length != 0){
      delete[] X1_ijka2[i];
      delete[] X2_iajk2[i];
    }
    length = np * nhpp1;
    if(length != 0){
      delete[] X1_iabc3[i];
      delete[] X_iabc3[i];
      delete[] X_abic3[i];
      delete[] X_abic4[i];
    }
    length = nh * nhhp1;
    if(length != 0){
      delete[] X2_iajk3[i];
      delete[] X2_iajk4[i];
    }
    length = np * nppp;
    if(length != 0){
      delete[] X1_abcd2[i];
      delete[] X1_abcd3[i];
      delete[] V_abcd[i];
    }
    length = nh * nhhh;
    if(length != 0){
      delete[] X_ijkl2[i];
      delete[] X_ijkl3[i];
      delete[] V_ijkl[i];
    }
    length = nh * nhpp1;
    if(length != 0){
      delete[] X1_iajb3[i];
      delete[] X1_iajb4[i];
      delete[] X3_iajb3[i];
      delete[] X_iajb3[i];
    }
    length = np * nhhp1;
    if(length != 0){
      delete[] X1_iajb2[i];
      delete[] X3_iajb2[i];
      delete[] X3_iajb5[i];
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      delete[] X_ijab1[i];
    }
    length = nhh * nhh;
    if(length != 0){
      delete[] X_ijkl1[i];
      delete[] X_ijkl4[i];
      delete[] Map_ijkl[i];
    }
    length = npp * npp;
    if(length != 0){
      delete[] X1_abcd1[i];
      delete[] X_abcd1[i];
      delete[] Map_abcd[i];
    }
    length = npp * nhp;
    if(length != 0){
      delete[] X_iabc5[i];
      delete[] X_abic7[i];
    }
    length = nhh * nhp;
    if(length != 0){
      delete[] X_ijka5[i];
      delete[] X2_iajk7[i];
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    length = npp1 * nhp1;
    if(length != 0){
      delete[] X_iabc4[i];
    }
    length = nhh1 * nhp1;
    if(length != 0){
      delete[] X_ijka4[i];
    }
    length = nhp2 * nhp2;
    if(length != 0){
      delete[] X1_iajb1[i];
      delete[] X3_iajb1[i];
      delete[] X_iajb1[i];
      delete[] Map_iajb[i];
    }
    length = npp1 * nhp2;
    if(length != 0){
      delete[] X_abic5[i];
      delete[] X_abic6[i];
    }
    length = nhh1 * nhp2;
    if(length != 0){
      delete[] X2_iajk5[i];
      delete[] X2_iajk6[i];
    }
  }
  delete[] X_ia2;
  delete[] X_ia3;
  delete[] X_ab2;
  delete[] X_ab3;
  delete[] X_ij2;
  delete[] X_ij3;
  delete[] X1_ij2;
  delete[] X1_ij3;
  delete[] X_ai2;
  delete[] X_ai3;
  delete[] X_ijab1;
  delete[] X1_iabc1;
  delete[] X1_iabc2;
  delete[] X1_iabc3;
  delete[] X_iabc1;
  delete[] X_iabc3;
  delete[] X_iabc4;
  delete[] X_iabc5;
  delete[] Map_iabc;
  delete[] X1_ijka1;
  delete[] X1_ijka2;
  delete[] X_ijka1;
  delete[] X_ijka4;
  delete[] X_ijka5;
  delete[] Map_ijka;
  delete[] X1_abcd1;
  delete[] X1_abcd2;
  delete[] X1_abcd3;
  delete[] X_abcd1;
  delete[] V_abcd;
  delete[] Map_abcd;
  delete[] X_ijkl1;
  delete[] X_ijkl2;
  delete[] X_ijkl3;
  delete[] X_ijkl4;
  delete[] V_ijkl;
  delete[] Map_ijkl;
  delete[] X1_iajb1;
  delete[] X1_iajb2;
  delete[] X1_iajb3;
  delete[] X1_iajb4;
  delete[] X3_iajb1;
  delete[] X3_iajb2;
  delete[] X3_iajb3;
  delete[] X3_iajb5;
  delete[] X_iajb1;
  delete[] X_iajb3;
  delete[] Map_iajb;
  delete[] X_abic1;
  delete[] X_abic2;
  delete[] X_abic3;
  delete[] X_abic4;
  delete[] X_abic5;
  delete[] X_abic6;
  delete[] X_abic7;
  delete[] Map_abic;
  delete[] X2_iajk1;
  delete[] X2_iajk2;
  delete[] X2_iajk3;
  delete[] X2_iajk4;
  delete[] X2_iajk5;
  delete[] X2_iajk6;
  delete[] X2_iajk7;
  delete[] X_iajk1;
  delete[] Map_iajk;
}

void CC_Eff::set_X_ia(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhp1[Chan.ind0]; ++i){
    x = X_ia1[i];
    X_ia2[Map_ia[3*i]][Map_ia[3*i + 1]] = x;
    X_ia3[Map_ia[3*i]][Map_ia[3*i + 2]] = x;
  }
}

void CC_Eff::set_X_ab(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.npp1[Chan.ind0]; ++i){
    x = 0.0;
    x += X_ab1[i];
    x += X_ab2[Map_ab[3*i]][Map_ab[3*i + 1]];
    x += X_ab3[Map_ab[3*i]][Map_ab[3*i + 2]];
    X_ab1[i] = x;
    X_ab2[Map_ab[3*i]][Map_ab[3*i + 1]] = x;
    X_ab3[Map_ab[3*i]][Map_ab[3*i + 2]] = x;
  }
}

void CC_Eff::set_X_ij(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhh1[Chan.ind0]; ++i){
    x = 0.0;
    x += X_ij1[i];
    x += X_ij2[Map_ij[3*i]][Map_ij[3*i + 1]];
    x += X_ij3[Map_ij[3*i]][Map_ij[3*i + 2]];
    X_ij1[i] = x;
    X_ij2[Map_ij[3*i]][Map_ij[3*i + 1]] = x;
    X_ij3[Map_ij[3*i]][Map_ij[3*i + 2]] = x;
  }
}

void CC_Eff::set_X1_ij(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhh1[Chan.ind0]; ++i){
    x = 0.0;
    x += X1_ij1[i];
    x += X1_ij2[Map_ij[3*i]][Map_ij[3*i + 1]];
    x += X1_ij3[Map_ij[3*i]][Map_ij[3*i + 2]];
    X1_ij1[i] = x;
    X1_ij2[Map_ij[3*i]][Map_ij[3*i + 1]] = x;
    X1_ij3[Map_ij[3*i]][Map_ij[3*i + 2]] = x;
  }
}

void CC_Eff::set_X_ai(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhp1[Chan.ind0]; ++i){
    x = 0.0;
    x += X_ai1[i];
    x += X_ai2[Map_ai[3*i]][Map_ai[3*i + 1]];
    x += X_ai3[Map_ai[3*i]][Map_ai[3*i + 2]];
    X_ai1[i] = x;
    X_ai2[Map_ai[3*i]][Map_ai[3*i + 1]] = x;
    X_ai3[Map_ai[3*i]][Map_ai[3*i + 2]] = x;
  }
}

void CC_Eff::set_X1_iabc(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int p1 = 0; p1 < Chan.np[i]; ++p1){
      for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
	ind0 = p1*Chan.nhpp[i] + hpp1;
	X1_iabc2[Map_iabc[i][8*ind0]][Map_iabc[i][8*ind0 + 1]] = X1_iabc1[i][ind0];
	X1_iabc3[Map_iabc[i][8*ind0 + 2]][Map_iabc[i][8*ind0 + 3]] = X1_iabc1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X_iabc(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int p1 = 0; p1 < Chan.np[i]; ++p1){
      for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
	ind0 = p1*Chan.nhpp[i] + hpp1;
	X_iabc3[Map_iabc[i][8*ind0 + 2]][Map_iabc[i][8*ind0 + 3]] = X_iabc1[i][ind0];
	X_iabc4[Map_iabc[i][8*ind0 + 4]][Map_iabc[i][8*ind0 + 5]] = X_iabc1[i][ind0];
	X_iabc5[Map_iabc[i][8*ind0 + 6]][Map_iabc[i][8*ind0 + 7]] = X_iabc1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X1_ijka(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
      for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
	ind0 = h1*Chan.nhhp[i] + hhp1;
	X1_ijka2[Map_ijka[i][6*ind0]][Map_ijka[i][6*ind0 + 1]] = X1_ijka1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X_ijka(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
      for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
	ind0 = h1*Chan.nhhp[i] + hhp1;
	X_ijka4[Map_ijka[i][6*ind0 + 2]][Map_ijka[i][6*ind0 + 3]] = X_ijka1[i][ind0];
	X_ijka5[Map_ijka[i][6*ind0 + 4]][Map_ijka[i][6*ind0 + 5]] = X_ijka1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X1_abcd(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size1; ++i){
    for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
      for(int pp2 = 0; pp2 < Chan.npp[i]; ++pp2){
	x = 0.0;
	ind0 = Chan.npp[i]*pp1 + pp2;
	V_abcd[Map_abcd[i][6*ind0 + 4]][Map_abcd[i][6*ind0 + 5]] = X1_abcd1[i][ind0];
	x += X1_abcd1[i][ind0];
	x += X1_abcd2[Map_abcd[i][6*ind0]][Map_abcd[i][6*ind0 + 1]];
	x += X1_abcd3[Map_abcd[i][6*ind0 + 2]][Map_abcd[i][6*ind0 + 3]];
	X1_abcd1[i][ind0] = x;
	X1_abcd2[Map_abcd[i][6*ind0]][Map_abcd[i][6*ind0 + 1]] = x;
	X1_abcd3[Map_abcd[i][6*ind0 + 2]][Map_abcd[i][6*ind0 + 3]] = x;
      }
    }
  }
}

void CC_Eff::set_X_ijkl(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size1; ++i){
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hh2 = 0; hh2 < Chan.nhh[i]; ++hh2){
	x = 0.0;
	ind0 = Chan.nhh[i]*hh1 + hh2;
	V_ijkl[Map_ijkl[i][8*ind0 + 6]][Map_ijkl[i][8*ind0 + 7]] = X_ijkl1[i][ind0];
	x += X_ijkl1[i][ind0];
	x += X_ijkl2[Map_ijkl[i][8*ind0]][Map_ijkl[i][8*ind0 + 1]];
	x += X_ijkl3[Map_ijkl[i][8*ind0 + 2]][Map_ijkl[i][8*ind0 + 3]];
	x += X_ijkl4[Map_ijkl[i][8*ind0 + 4]][Map_ijkl[i][8*ind0 + 5]];
	X_ijkl1[i][ind0] = x;
	X_ijkl2[Map_ijkl[i][8*ind0]][Map_ijkl[i][8*ind0 + 1]] = x;
	X_ijkl3[Map_ijkl[i][8*ind0 + 2]][Map_ijkl[i][8*ind0 + 3]] = x;
	X_ijkl4[Map_ijkl[i][8*ind0 + 4]][Map_ijkl[i][8*ind0 + 5]] = x;
      }
    }
  }
}

void CC_Eff::set_X1_iajb(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size2; ++i){
    for(int hp1 = 0; hp1 < Chan.nhp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.nhp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.nhp2[i]*hp1 + hp2;
	x += X1_iajb1[i][ind0];
	x += X1_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]];
	x += X1_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]];
	X1_iajb1[i][ind0] = x;
	X1_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]] = x;
	X1_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]] = x;
	X1_iajb4[Map_iajb[i][8*ind0 + 4]][Map_iajb[i][8*ind0 + 5]] = x;
      }
    }
  }
}

void CC_Eff::set_X3_iajb(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size2; ++i){
    for(int hp1 = 0; hp1 < Chan.nhp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.nhp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.nhp2[i]*hp1 + hp2;
	x += X3_iajb1[i][ind0];
	x += X3_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]];
	x += X3_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]];
	X3_iajb1[i][ind0] = x;
	X3_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]] = x;
	X3_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]] = x;
	X3_iajb5[Map_iajb[i][8*ind0 + 6]][Map_iajb[i][8*ind0 + 7]] = x;
      }
    }
  }
}

void CC_Eff::set_X_iajb(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size2; ++i){
    for(int hp1 = 0; hp1 < Chan.nhp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.nhp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.nhp2[i]*hp1 + hp2;
	x += X_iajb1[i][ind0];
	x += X_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]];
	X_iajb1[i][ind0] = x;
	X_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]] = x;
      }
    }
  }
}

void CC_Eff::set_X_abic(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size3; ++i){
    for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
      for(int p1 = 0; p1 < Chan.np[i]; ++p1){
	x = 0.0;
	ind0 = Chan.np[i]*hpp1 + p1;
	x += X_abic1[i][ind0];
	x += X_abic2[Map_abic[i][12*ind0]][Map_abic[i][12*ind0 + 1]];
	x += X_abic3[Map_abic[i][12*ind0 + 2]][Map_abic[i][12*ind0 + 3]];
	x += X_abic4[Map_abic[i][12*ind0 + 4]][Map_abic[i][12*ind0 + 5]];
	x += X_abic5[Map_abic[i][12*ind0 + 6]][Map_abic[i][12*ind0 + 7]];
	x += X_abic6[Map_abic[i][12*ind0 + 8]][Map_abic[i][12*ind0 + 9]];
	x += X_abic7[Map_abic[i][12*ind0 + 10]][Map_abic[i][12*ind0 + 11]];
	X_abic1[i][ind0] = x;
	X_abic2[Map_abic[i][12*ind0]][Map_abic[i][12*ind0 + 1]] = x;
	X_abic3[Map_abic[i][12*ind0 + 2]][Map_abic[i][12*ind0 + 3]] = x;
	X_abic4[Map_abic[i][12*ind0 + 4]][Map_abic[i][12*ind0 + 5]] = x;
	X_abic5[Map_abic[i][12*ind0 + 6]][Map_abic[i][12*ind0 + 7]] = x;
	X_abic6[Map_abic[i][12*ind0 + 8]][Map_abic[i][12*ind0 + 9]] = x;
	X_abic7[Map_abic[i][12*ind0 + 10]][Map_abic[i][12*ind0 + 11]] = x;
      }
    }
  }
}

void CC_Eff::set_X2_iajk(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size3; ++i){
    for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
      for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
	x = 0.0;
	ind0 = Chan.nh[i]*hhp1 + h1;
	x += X2_iajk1[i][ind0];
	x += X2_iajk2[Map_iajk[i][12*ind0]][Map_iajk[i][12*ind0 + 1]];
	x += X2_iajk3[Map_iajk[i][12*ind0 + 2]][Map_iajk[i][12*ind0 + 3]];
	x += X2_iajk4[Map_iajk[i][12*ind0 + 4]][Map_iajk[i][12*ind0 + 5]];
	x += X2_iajk5[Map_iajk[i][12*ind0 + 6]][Map_iajk[i][12*ind0 + 7]];
	x += X2_iajk6[Map_iajk[i][12*ind0 + 8]][Map_iajk[i][12*ind0 + 9]];
	x += X2_iajk7[Map_iajk[i][12*ind0 + 10]][Map_iajk[i][12*ind0 + 11]];
	X2_iajk1[i][ind0] = x;
	X2_iajk2[Map_iajk[i][12*ind0]][Map_iajk[i][12*ind0 + 1]] = x;
	X2_iajk3[Map_iajk[i][12*ind0 + 2]][Map_iajk[i][12*ind0 + 3]] = x;
	X2_iajk4[Map_iajk[i][12*ind0 + 4]][Map_iajk[i][12*ind0 + 5]] = x;
	X2_iajk5[Map_iajk[i][12*ind0 + 6]][Map_iajk[i][12*ind0 + 7]] = x;
	X2_iajk6[Map_iajk[i][12*ind0 + 8]][Map_iajk[i][12*ind0 + 9]] = x;
	X2_iajk7[Map_iajk[i][12*ind0 + 10]][Map_iajk[i][12*ind0 + 11]] = x;
      }
    }
  }
}


//Initialize program from input file
void Get_Input_Parameters(std::string &infile, Input_Parameters &Input)
{ 
  std::string path; // Input File Path
  std::string line; // Line from input file
  std::ifstream filestream; // File stream
  int index; // Count input line
  size_t colon; // Size for finding ':' which proceeds every input
  std::string substr; // Substd::string for extracting input

  path = PATH + infile;
  filestream.open(path.c_str());
  if (!filestream.is_open()){ std::cerr << "Input file, " << path << ", does not exist" << std::endl; exit(1); }

  //find lines that start with '\*'
  index = 0;
  while (getline(filestream, line)){ 
    if(line[0] == '\\' && line[1] == '*'){  
      ++index;
      colon = line.find(':');
      if( colon == line.size() - 1 ){ continue; };
      substr = line.substr(colon + 2, line.size());
      switch(index){
      case 1:
	Input.calc_case = substr;
	break;
      case 2:
	Input.basis = substr;
	break;
      case 3:
	Input.approx = substr;
	break;
      case 4:
	Input.obstrength = atof(substr.c_str());
	break;
      case 5:
	Input.tbstrength = atof(substr.c_str());
	break;
      case 6:
	Input.Pshells = atoi(substr.c_str());
	break;
      case 7:
	Input.Nshells = atoi(substr.c_str());
	break;
      case 8:
	Input.LevelScheme = substr;
	break;
      case 9:
	Input.MatrixElements = substr;
	break;
      case 10:
	Input.extra = atoi(substr.c_str());
	break;
      } 
    }
    else{ continue; };
  }
}


void Print_Parameters(const Input_Parameters &Parameters, const Model_Space &Space)
{
  std::cout << std::endl;
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << "Case = " << Parameters.calc_case << ", Basis = " << Parameters.basis << ", Approximation = " << Parameters.approx << std::endl;
  if(Parameters.LevelScheme.size() > 0){ 
    std::cout << "Levels Scheme = " << Parameters.LevelScheme << std::endl;
    if(Parameters.MatrixElements.size() > 0){ std::cout << "Interaction = " << Parameters.MatrixElements << std::endl; }
  }
  if(Parameters.calc_case == "nuclear"){
    if(Parameters.basis == "finite_J"){ std::cout << "Total States = " << Space.indtot << std::endl; }
    else{ std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.indtot << std::endl; }
    std::cout << "Proton Shells = " << Parameters.Pshells << ", Neutron Shells = " << Parameters.Nshells << std::endl;
    std::cout << "Protons = " << Parameters.P << ", Neutrons = " << Parameters.N << std::endl;
    if(Parameters.calc_case == "infinite"){ std::cout << "Density = " << Parameters.density << std::endl; }
  }
  else if(Parameters.calc_case == "electronic"){
    std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.indtot << std::endl;
    std::cout << "Electron Shells = " << Parameters.Pshells << ", Electrons = " << Parameters.P << std::endl;
    std::cout << "Density = " << Parameters.density << std::endl;
  }
  else if(Parameters.calc_case == "quantum_dot"){
    std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.indtot << std::endl;
    std::cout << "Electron Shells = " << Parameters.Pshells << ", Electrons = " << Parameters.P << std::endl;
    std::cout << "Oscillator Energy = " << Parameters.density << std::endl;
  }
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;
}

void Model_Space::delete_struct(Input_Parameters &Parameters)
{
  delete[] qnums;
  if(Parameters.basis == "infinite"){ delete[] map_2b; }
  if(Parameters.basis == "finite_J"){ delete[] shellsm; }
}

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  std::string fullpath; // Model Space file path
  std::string phline; // Std::String for each file line
  std::ifstream splevels; // Model space file

  int ind; // state index
  int l; // orbital angular momentum quantum number
  double energy; // state energy

  int pcount = 0; // proton state count
  int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  int nhcount = 0; // neutron hole state count
  
  fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open()){ std::cerr << "Level Scheme file does not exist" << std::endl; exit(1); };

  getline(splevels, phline);
  std::stringstream(phline) >> Space.indtot;

  // allocate memory for quntum numbers for each state
  Space.qnums = new State[Space.indtot];

  // initialize mins and maxs for each quantum number
  Space.qmins.t = 1000;
  Space.qmins.m = 1000;
  Space.qmins.nx = 1000;
  Space.qmins.ny = 1000;
  Space.qmins.nz = 1000;
  Space.qmins.ml = 1000;
  Space.qmins.j = 1000;
  Space.qmins.par = 1000;
  Space.qmaxs.t = -1000;
  Space.qmaxs.m = -1000;
  Space.qmaxs.nx = -1000;
  Space.qmaxs.ny = -1000;
  Space.qmaxs.nz = -1000;
  Space.qmaxs.ml = -1000;
  Space.qmaxs.j = -1000;
  Space.qmaxs.par = -1000;

  // States must be ordered by energy from low to high
  for(int i = 0; i < Space.indtot; ++i){
    getline(splevels, phline);
    if(Parameters.basis == "infinite"){ // ind, nx, ny, nz, sz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].nx >> Space.qnums[i].ny >> Space.qnums[i].nz >> Space.qnums[i].m >> Space.qnums[i].t >> energy;
      Space.qnums[i].ml = 0;
      Space.qnums[i].par = 1;
      Space.qnums[i].j = 0;
      if(Space.qnums[i].nx < Space.qmins.nx){ Space.qmins.nx = Space.qnums[i].nx; }
      if(Space.qnums[i].nx > Space.qmaxs.nx){ Space.qmaxs.nx = Space.qnums[i].nx; }
      if(Space.qnums[i].ny < Space.qmins.ny){ Space.qmins.ny = Space.qnums[i].ny; }
      if(Space.qnums[i].ny > Space.qmaxs.ny){ Space.qmaxs.ny = Space.qnums[i].ny; }
      if(Space.qnums[i].nz < Space.qmins.nz){ Space.qmins.nz = Space.qnums[i].nz; }
      if(Space.qnums[i].nz > Space.qmaxs.nz){ Space.qmaxs.nz = Space.qnums[i].nz; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_M"){ // ind, n, l, j, jz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> l >> Space.qnums[i].j >> Space.qnums[i].m >> Space.qnums[i].t >> energy;
      Space.qnums[i].par = -2*(l%2) + 1;
      Space.qnums[i].ml = 0;
      Space.qnums[i].nx = 0;
      Space.qnums[i].ny = 0;
      Space.qnums[i].nz = 0;
      if(Space.qnums[i].par < Space.qmins.par){ Space.qmins.par = Space.qnums[i].par; }
      if(Space.qnums[i].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[i].par; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_HO"){ // ind, n, l, lz, sz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> l >> Space.qnums[i].ml >> Space.qnums[i].m >> Space.qnums[i].t >> energy;
      Space.qnums[i].par = -2*(l%2) + 1;
      Space.qnums[i].j = 0;
      Space.qnums[i].nx = 0;
      Space.qnums[i].ny = 0;
      Space.qnums[i].nz = 0;
      if(Space.qnums[i].ml < Space.qmins.ml){ Space.qmins.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].ml > Space.qmaxs.ml){ Space.qmaxs.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_J"){ // ind, n, l, j, tz, l2n
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> l >> Space.qnums[i].j >> Space.qnums[i].t >> ind >> energy;
      Space.qnums[i].par = -2*(l%2) + 1;
      Space.qnums[i].m = 0;
      Space.qnums[i].ml = 0;
      Space.qnums[i].nx = 0;
      Space.qnums[i].ny = 0;
      Space.qnums[i].nz = 0;
      if(Space.qnums[i].par < Space.qmins.par){ Space.qmins.par = Space.qnums[i].par; }
      if(Space.qnums[i].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[i].par; }
      if(Space.qnums[i].j < Space.qmins.j){ Space.qmins.j = Space.qnums[i].j; }
      if(Space.qnums[i].j > Space.qmaxs.j){ Space.qmaxs.j = Space.qnums[i].j; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    Space.qnums[i].energy = energy * Parameters.obstrength;
  }
  splevels.close();

  /*for(int i = 0; i < Space.indtot; ++i){
    std::cout << i+1 << " " << Space.qnums[i].ml << " " << Space.qnums[i].par << " " << Space.qnums[i].t << "   " << Space.qnums[i].energy << std::endl;
    }*/

  // Determine shell structure and assign hole/particle labels
  int *pshell;
  int *nshell;
  int pshell_num = 0;
  int nshell_num = 0;
  double p_en = -1000.0;
  double n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }
  pshell = new int[pshell_num];
  nshell = new int[nshell_num];
  pshell_num = 0;
  nshell_num = 0;
  p_en = -1000.0;
  n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      pshell[pshell_num] = i;
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      nshell[nshell_num] = i;
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }

  if(Parameters.Pshells > pshell_num){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells > nshell_num){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      ++pcount;
      if(i < pshell[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++parcount; }
    }
    else if(Space.qnums[i].t == 1){
      ++ncount;
      if(i < nshell[Parameters.Nshells]){ Space.qnums[i].type = "hole"; ++holcount; ++nhcount; }
      else{ Space.qnums[i].type = "particle"; ++parcount; }
    }
  }
  delete[] pshell;
  delete[] nshell;
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  // Find range of 2body quantum numbers from mins and maxs
  Space.qsizes.t = Space.qmaxs.t - Space.qmins.t + 1;
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.nx = 2*(Space.qmaxs.nx - Space.qmins.nx) + 1;
  Space.qsizes.ny = 2*(Space.qmaxs.ny - Space.qmins.ny) + 1;
  Space.qsizes.nz = 2*(Space.qmaxs.nz - Space.qmins.nz) + 1;
  Space.qsizes.ml = 2*(Space.qmaxs.ml - Space.qmins.ml) + 1;
  Space.qsizes.j = Space.qmaxs.j + 1;
  Space.qsizes.par = int((Space.qmaxs.par - Space.qmins.par)/2) + 1;

  // Get total number of 2body states (and setup index array for infinite case)
  if(Parameters.basis == "infinite"){
    int N2max = 0;
    int N2maxtemp;
    for(int i = 0; i < Space.indtot; ++i){
      for(int j = i+1; j < Space.indtot; ++j){
	N2maxtemp = (Space.qnums[i].nx + Space.qnums[j].nx)*(Space.qnums[i].nx + Space.qnums[j].nx) + 
	  (Space.qnums[i].ny + Space.qnums[j].ny)*(Space.qnums[i].ny + Space.qnums[j].ny) + 
	  (Space.qnums[i].nz + Space.qnums[j].nz)*(Space.qnums[i].nz + Space.qnums[j].nz);
	if(N2maxtemp > N2max){ N2max = N2maxtemp; }
      }
    }
    Space.map_2b = new int[Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz];
    int count1 = 0;
    for(int nx = 2*Space.qmins.nx; nx <= 2*Space.qmaxs.nx; ++nx){
      for(int ny = 2*Space.qmins.ny; ny <= 2*Space.qmaxs.ny; ++ny){
	for(int nz = 2*Space.qmins.nz; nz <= 2*Space.qmaxs.nz; ++nz){
	  if(nx*nx + ny*ny + nz*nz <= N2max){
	    Space.map_2b[(nx - 2*Space.qmins.nx)*(Space.qsizes.ny*Space.qsizes.nz) + (ny - 2*Space.qmins.ny)*Space.qsizes.nz + (nz - 2*Space.qmins.nz)] = count1;
	    ++count1;
	  }
	  else{
	    Space.map_2b[(nx - 2*Space.qmins.nx)*(Space.qsizes.ny*Space.qsizes.nz) + (ny - 2*Space.qmins.ny)*Space.qsizes.nz + (nz - 2*Space.qmins.nz)] = -1;
	  }
	}
      }
    }
    Space.size_2b = count1 * Space.qsizes.t * Space.qsizes.m;
  }
  else if(Parameters.basis == "finite_M"){
    Space.size_2b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
  }
  else if(Parameters.basis == "finite_HO"){
    Space.size_2b = Space.qsizes.ml * Space.qsizes.m * Space.qsizes.t;
  }
  else if(Parameters.basis == "finite_J"){
    Space.size_2b = Space.qsizes.par * Space.qsizes.j * Space.qsizes.t;
  }
}

void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space &Space)
{
  int ind; // state index
  int m; // angular momenum projection
  int shelllength; // number of single particle states for each shell

  int pcount = 0; // proton state count
  int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  int nhcount = 0; // neutron hole state count

  State *states = Space.qnums;
  int indtotj = Space.indtot;
  Space.shellsm = new int*[Space.indtot];
  // reset Space.indtot with degeneracies 2j + 1
  Space.indtot = 0;
  for(int i = 0; i < indtotj; ++i){ Space.indtot += Space.qnums[i].j + 1; }

  // allocate memory for quntum numbers for each state
  delete Space.qnums;
  Space.qnums = new State[Space.indtot];

  // initialize mins and maxs for each quantum number
  Space.qmins.t = 1000;
  Space.qmins.m = 1000;
  Space.qmins.j = 1000;
  Space.qmins.par = 1000;
  Space.qmaxs.t = -1000;
  Space.qmaxs.m = -1000;
  Space.qmaxs.j = -1000;
  Space.qmaxs.par = -1000;

  ind = 0;
  for(int i = 0; i < indtotj; ++i){
    shelllength = states[i].j + 1;
    Space.shellsm[i] = new int[shelllength];
    for(int j = 0; j < shelllength; ++j){
      m = -states[i].j + 2*j;
      Space.qnums[ind].par = states[i].par;
      Space.qnums[ind].j = states[i].j;
      Space.qnums[ind].m = m;
      Space.qnums[ind].t = states[i].t;
      Space.qnums[ind].energy = states[i].energy;
      Space.qnums[ind].type = states[i].type;
      Space.qnums[ind].ml = 0;
      Space.qnums[ind].nx = 0;
      Space.qnums[ind].ny = 0;
      Space.qnums[ind].nz = 0;
      Space.shellsm[i][j] = ind;
      if(Space.qnums[ind].par < Space.qmins.par){ Space.qmins.par = Space.qnums[ind].par; }
      if(Space.qnums[ind].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[ind].par; }
      if(Space.qnums[ind].m < Space.qmins.m){ Space.qmins.m = Space.qnums[ind].m; }
      if(Space.qnums[ind].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[ind].m; }
      if(Space.qnums[ind].t < Space.qmins.t){ Space.qmins.t = Space.qnums[ind].t; }
      if(Space.qnums[ind].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[ind].t; }
      ++ind;
    }
  }

  // Determine shell structure and assign hole/particle labels
  int *pshell;
  int *nshell;
  int pshell_num = 0;
  int nshell_num = 0;
  double p_en = -1000.0;
  double n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }
  pshell = new int[pshell_num];
  nshell = new int[nshell_num];
  pshell_num = 0;
  nshell_num = 0;
  p_en = -1000.0;
  n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      pshell[pshell_num] = i;
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      nshell[nshell_num] = i;
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }

  if(Parameters.Pshells > pshell_num){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells > nshell_num){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      if(i < pshell[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.qnums[i].t == 1){
      if(i < nshell[Parameters.Nshells]){ Space.qnums[i].type = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.qnums[i].type = "particle"; ++ncount; ++parcount; }
    }
  }
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;
  delete[] pshell;
  delete[] nshell;

  Space.qsizes.t = Space.qmaxs.t - Space.qmins.t + 1;
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.par = int((Space.qmaxs.par - Space.qmins.par)/2) + 1;
  Space.size_2b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
}


void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;
  double electron_prefac = hbarc_HartA * hbarc_HartA / (2.0 * m_electronc2_Hart);
  double neutron_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_neutronc2);
  double proton_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_protonc2);

  int count = 0; // total state count
  int pcount = 0; // proton state count
  int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  int nhcount = 0; // neutron hole state count

  int shellnums [] = {1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46};

  // Find appropriate Nmax for the number of shells using shellnums
  if(Parameters.Shells > 40){ std::cerr << "Nmax too big!" << std::endl; exit(1); }
  Space.Nmax = shellnums[Parameters.Shells - 1];

  if(Parameters.calc_case == "electronic"){
    Parameters.Nshells = 0;
    double r_b = hbarc_HartA / (m_electronc2_Hart * fine_struct);
    Parameters.density = 3.0/(4.0 * PI * pow(Parameters.density * r_b, 3));
  }

  // Find maximum number of states (indtot) depending on Nmax and whether or not there are protons/neutrons
  if(Parameters.Pshells != 0 && Parameters.Nshells != 0){ Space.qmins.t = -1, Space.qmaxs.t = 1, Space.qsizes.t = 2; }
  else if(Parameters.Pshells != 0 && Parameters.Nshells == 0){ Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1; }
  else if(Parameters.Pshells == 0 && Parameters.Nshells != 0){ Space.qmins.t = 1, Space.qmaxs.t = 1, Space.qsizes.t = 1; }
  else if(Parameters.calc_case == "nuclear"){ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }
  else if(Parameters.calc_case == "electronic"){ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  // Find total number of states
  Space.nmax = std::floor(std::sqrt(Space.Nmax));
  for(int shell = 0; shell <= Space.Nmax; ++shell){
    for(int nx = -Space.nmax; nx <= Space.nmax; ++nx){
      for(int ny = -Space.nmax; ny <= Space.nmax; ++ny){
	for(int nz = -Space.nmax; nz <= Space.nmax; ++nz){
	  if(shell != nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    for( int tz = -1; tz <= 1; tz = tz+2){
	      if(tz == -1 && Parameters.Pshells == 0){ continue; }
	      if(tz == 1 && Parameters.Nshells == 0){ continue; }
	      count++;
	    }
	  }
	}
      }
    }
  }
  Space.indtot = count;
  count = 0;

  // allocate memory for quantum numbers for each state
  Space.qnums = new State[Space.indtot];
  for(int shell = 0; shell <= Space.Nmax; ++shell){
    for(int nx = -Space.nmax; nx <= Space.nmax; ++nx){    
      for(int ny = -Space.nmax; ny <= Space.nmax; ++ny){	
	for(int nz = -Space.nmax; nz <= Space.nmax; ++nz){	  
	  if(shell != nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    for( int tz = -1; tz <= 1; tz = tz+2){
	      if(tz == -1){
		if(Parameters.Pshells == 0){ continue; }
		++pcount;
		E = 4.0*(nx*nx + ny*ny + nz*nz);
		Space.qnums[count].energy = E;
		if(shell < Parameters.Pshells){ Space.qnums[count].type = "hole"; ++holcount; ++phcount; }
		else{ Space.qnums[count].type = "particle"; ++parcount; }
	      }
	      if(tz == 1){
		if(Parameters.Nshells == 0){ continue; }
		++ncount;
		E = 4.0*(nx*nx + ny*ny + nz*nz);
		Space.qnums[count].energy = E;
		if(shell < Parameters.Nshells){ Space.qnums[count].type = "hole"; ++holcount; ++nhcount; }
		else{ Space.qnums[count].type = "particle"; ++parcount; }
	      }
	      Space.qnums[count].nx = nx;
	      Space.qnums[count].ny = ny;
	      Space.qnums[count].nz = nz;
	      Space.qnums[count].m = sz;
	      Space.qnums[count].t = tz;
	      Space.qnums[count].ml = 0;
	      Space.qnums[count].par = 1;
	      Space.qnums[count].j = 0;
	      count++;
	    }
	  }   
	} 
      }
    }
  }
  Space.qsizes.m = 3; // -2, 0, 2
  Space.qsizes.nx = 4*Space.nmax + 1;
  Space.qsizes.ny = 4*Space.nmax + 1;
  Space.qsizes.nz = 4*Space.nmax + 1;
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  // With the number of hole states, the length scale and state energies can be calculated
  double L = pow(holcount/Parameters.density, 1.0/3.0);
  if(Parameters.calc_case == "nuclear"){
    for(int i = 0; i < Space.indtot; ++i){
      if(Space.qnums[i].t == -1){ Space.qnums[i].energy *= proton_prefac * M_PI*M_PI / (L*L); }
      else if(Space.qnums[i].t == 1){ Space.qnums[i].energy *= neutron_prefac * M_PI*M_PI / (L*L); }
    }
    // Change energies to Hartree-Fock energies, E_p = E_p + 2*V_pipi
    for(int p = 0; p < Space.indtot; ++p){
      for(int i = 0; i < Space.indhol; ++i){
	if(p == i){ continue; }
	Space.qnums[p].energy += 2*vint_Minnesota_Momentum(Space, p, i, p, i, L);
      }
    }
  }
  else if(Parameters.calc_case == "electronic"){
    for(int i = 0; i < Space.indtot; ++i){ Space.qnums[i].energy *= electron_prefac * M_PI*M_PI / (L*L); }
    // Change energies to Hartree-Fock energies, E_p = E_p + 2*V_pipi
    for(int p = 0; p < Space.indtot; ++p){
      for(int i = 0; i < Space.indhol; ++i){
	if(p == i){ continue; }
	Space.qnums[p].energy += 2*Coulomb_Inf(Space, p, i, p, i, L);;
      }
    }
  }

  // Order two-body states with respect to Nx, Ny, Nz for channel index function
  int count1 = 0;
  int N2max = 4*Space.Nmax;
  Space.map_2b = new int[Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz];
  for(int nx = -2 * Space.nmax; nx <= 2 * Space.nmax; ++nx){
    for(int ny = -2 * Space.nmax; ny <= 2 * Space.nmax; ++ny){
      for(int nz = -2 * Space.nmax; nz <= 2 * Space.nmax; ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max){
	  Space.map_2b[(nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz + (ny + 2*Space.nmax) * Space.qsizes.nz + (nz + 2*Space.nmax)] = count1;
	  ++count1;
	}
	  else{
	    Space.map_2b[(nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz + (ny + 2*Space.nmax) * Space.qsizes.nz + (nz + 2*Space.nmax)] = -1;
	  }
      }
    }
  }
  Space.size_2b = count1 * Space.qsizes.t * Space.qsizes.m;
}


void QD_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;

  int count = 0; // total state count
  int pcount = 0; // proton state count
  //int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  //int nhcount = 0; // neutron hole state count

  Parameters.Nshells = 0;

  Space.qmins.ml = 1000;
  Space.qmins.par = 1000;
  Space.qmaxs.ml = -1000;
  Space.qmaxs.par = -1000;
  Space.nmax = -1000;

  if(Parameters.Pshells != 0){ Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1, Space.qsizes0.t = 1; }
  else{ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  // Find total number of states
  for(int shell = 0; shell < Parameters.Shells; ++shell){
    for(int n = 0; n < Parameters.Shells; ++n){
      for(int ml = -1 * int(shell - 2*n); ml <= int(shell - 2*n); ++ml){
	if(2*n + abs(ml) != shell){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  if(n > Space.nmax){ Space.nmax = n; }
	  if(ml < Space.qmins.ml){ Space.qmins.ml = ml; }
	  if(ml > Space.qmaxs.ml){ Space.qmaxs.ml = ml; }
	  count++;
	}
      }
    }
  }
  Space.indtot = count;
  count = 0;

  // Allocate memory for quantum numbers for each state
  Space.qnums = new State[Space.indtot];
  Space.nmax = 0;
  for(int shell = 0; shell < Parameters.Shells; ++shell){
    for(int n = 0; n < Parameters.Shells; ++n){
      for(int ml = -1 * int(shell - 2*n); ml <= int(shell - 2*n); ++ml){
	if(2*n + abs(ml) != shell){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  E = Parameters.density*(shell + 1.0);
	  Space.qnums[count].energy = E;
	  Space.qnums[count].n = n;
	  Space.qnums[count].par = -2*(abs(ml)%2) + 1;
	  Space.qnums[count].ml = ml;
	  Space.qnums[count].m = sz;
	  Space.qnums[count].t = -1;
	  Space.qnums[count].nx = 0;
	  Space.qnums[count].ny = 0;
	  Space.qnums[count].nz = 0;
	  Space.qnums[count].j = 0;
	  if(shell < Parameters.Pshells){ Space.qnums[count].type = "hole"; ++holcount; ++phcount; }
	  else{ Space.qnums[count].type = "particle"; ++parcount; }
	  //std::cout << count+1 << ": " << n << " " << Space.qnums[count].par << " " << ml << " " << sz << " " << Space.qnums[count].type << " : " << E << std::endl;
	  count++;
	}
      }
    }
  }
  Space.qsizes.m = 3; // -2, 0, +2
  Space.qsizes.ml = 2*(Space.qmaxs.ml - Space.qmins.ml) + 1;
  Space.qsizes.par = 2; // -1, +1
  Space.size_2b = Space.qsizes.ml * Space.qsizes.m * Space.qsizes.t;

  Space.qsizes0.m = 2; // -1, +1
  Space.qsizes0.ml = Space.qmaxs.ml - Space.qmins.ml + 1;

  int key;
  for(int i = 0; i < Space.indtot; ++i){
    key = ChanInd_1b(Parameters.basis, Space, Space.qnums[i]);
    Space.map_1b[key] = i;
    //std::cout << "! " << i << ", " << key << " : " << Space.map_1b[key] << std::endl;
  }

  Space.indp = pcount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = 0;
}

/*Model_Space CART_Build_Model_Space_Twist(Input_Parameters &Parameters, const double &thx, const double &thy, const double &thz)
{
  Model_Space Space; // Model Space information
  double E;
  int count, pcount, ncount, holcount, parcount, phcount, nhcount;
  
  double hbarc = 197.3269788; // MeVfm
  double m_neutronc2 = 939.565378; // MeV
  //double m_protonc2 = 938.272046; // MeV
  double m_protonc2 = 939.565378; // MeV
  double neutron_prefac = hbarc*hbarc/(2.0*m_neutronc2);
  double proton_prefac = hbarc*hbarc/(2.0*m_protonc2);
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1.0/3.0);

  if(Parameters.P != 0 && Parameters.N != 0){ Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2*2; }
  else if(Parameters.P != 0 || Parameters.N != 0){ Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2; }
  else{ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }

  Space.levelsind.resize(Space.indtot);
  Space.levelsm.resize(Space.indtot);
  Space.levelst.resize(Space.indtot);
  Space.levelstype.resize(Space.indtot);
  Space.levelsen.resize(Space.indtot);
  Space.levelsnx.resize(Space.indtot);
  Space.levelsny.resize(Space.indtot);
  Space.levelsnz.resize(Space.indtot);

  count = 0;
  pcount = 0;
  ncount = 0;
  parcount = 0;
  holcount = 0;
  phcount = 0;
  nhcount = 0;

  int nmax = std::ceil(std::sqrt(Parameters.Nmax));
  for(int nx = -nmax; nx <= nmax; ++nx){    
    for(int ny = -nmax; ny <= nmax; ++ny){	
      for(int nz = -nmax; nz <= nmax; ++nz){	  
	if(Parameters.Nmax < nx*nx + ny*ny + nz*nz){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  for( int tz = -1; tz <= 1; tz = tz+2){
	    if(tz == -1){
	      if(Parameters.P == 0){ continue; }
	      E = 4.0*(nx*nx + ny*ny + nz*nz) + (4.0/M_PI)*(nx*thx + ny*thy + nz*thz) + (1.0/(M_PI*M_PI))*(thx*thx + thy*thy + thz*thz);
	      Space.levelsen[count] = E;
	    }
	    if(tz == 1){
	      if(Parameters.N == 0){ continue; }
	      E = 4.0*(nx*nx + ny*ny + nz*nz) + (4.0/M_PI)*(nx*thx + ny*thy + nz*thz) + (1.0/(M_PI*M_PI))*(thx*thx + thy*thy + thz*thz);
	      Space.levelsen[count] = E;
	    }
	    Space.levelsind[count] = count + 1;
	    Space.levelsnx[count] = nx;
	    Space.levelsny[count] = ny;
	    Space.levelsnz[count] = nz;
	    Space.levelsm[count] = sz;
	    Space.levelst[count] = tz;
	    count++;
	  }
	}   
      } 
    }
  }

  //order states by energy
  double tempE;
  int tempind;
  int nx1, ny1, nz1;
  double sz1, tz1, en1;
  for(int i = 0; i < count; ++i){
    tempE = Space.levelsen[i];
    tempind = i;
    for(int j = i+1; j < count; ++j){
      if(Space.levelsen[j] < tempE){ tempE = Space.levelsen[j]; tempind = j; }
    }
    if(tempind != i){
      nx1 = Space.levelsnx[tempind];
      ny1 = Space.levelsny[tempind];
      nz1 = Space.levelsnz[tempind];
      sz1 = Space.levelsm[tempind];
      tz1 = Space.levelst[tempind];
      en1 = Space.levelsen[tempind];
      for(int j = tempind; j > i; --j){
	Space.levelsnx[j] = Space.levelsnx[j-1];
	Space.levelsny[j] = Space.levelsny[j-1];
	Space.levelsnz[j] = Space.levelsnz[j-1];
	Space.levelsm[j] = Space.levelsm[j-1];
	Space.levelst[j] = Space.levelst[j-1];
	Space.levelsen[j] = Space.levelsen[j-1];
      }
      Space.levelsnx[i] = nx1;
      Space.levelsny[i] = ny1;
      Space.levelsnz[i] = nz1;
      Space.levelsm[i] = sz1;
      Space.levelst[i] = tz1;
      Space.levelsen[i] = en1;
    }
  }

  for(int i = 0; i < count; ++i){
    if(Space.levelst[i] == -1){
      if(phcount < Parameters.P){ Space.levelstype[i] = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.levelstype[i] = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.levelst[i] == 1){
      if(nhcount < Parameters.N){ Space.levelstype[i] = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.levelstype[i] = "particle"; ++ncount; ++parcount; }
    }
  }

  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;

  for(int i = 0; i < count; ++i){
    if(Space.levelst[i] == -1){ Space.levelsen[i] *= proton_prefac*M_PI*M_PI/(L*L); }
    else if(Space.levelst[i] == 1){ Space.levelsen[i] *= neutron_prefac*M_PI*M_PI/(L*L); }
  }

  Space.levelsind.resize(count);
  Space.levelsm.resize(count);
  Space.levelst.resize(count);
  Space.levelstype.resize(count);
  Space.levelsen.resize(count);
  Space.levelsnx.resize(count);
  Space.levelsny.resize(count);
  Space.levelsnz.resize(count);
  Space.indtot = count;

  Space.nmax = int(std::floor(std::sqrt(Parameters.Nmax)));

  int count1 = 0;
  int Nxsize = 4*Space.nmax + 1;
  int Nysize, Nzsize;
  Space.Nxmin = -2*Space.nmax;
  Space.Nymin.resize(Nxsize);
  Space.Nzmin.resize(Nxsize);
  Space.CART_tb1Indvec.resize(Nxsize);
  for(int i = 0; i < Nxsize; ++i){
    Space.Nymin[i] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i))));
    Nysize = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i)))) + 1;
    Space.Nzmin[i].resize(Nysize);
    Space.CART_tb1Indvec[i].resize(Nysize);
    for(int j = 0; j < Nysize; ++j){
      Space.Nzmin[i][j] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i) - (Space.Nymin[i]+j)*(Space.Nymin[i]+j))));
      Nzsize = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i) - (Space.Nymin[i]+j)*(Space.Nymin[i]+j)))) + 1;
      Space.CART_tb1Indvec[i][j].resize(Nzsize);
      for(int k = 0; k < Nzsize; ++k){
	Space.CART_tb1Indvec[i][j][k] = count1;
	++count1;
      }
    }
  }
  Space.CART_tb1Indsize = count1;

  int count2 = 0;
  int Nx2size = 4*Space.nmax + 1;
  int Ny2size, Nz2size;
  Space.Nx2min = -2*Space.nmax;
  Space.Ny2min.resize(Nx2size);
  Space.Nz2min.resize(Nx2size);
  Space.CART_tb2Indvec.resize(Nx2size);
  for(int i = 0; i < Nx2size; ++i){
    Space.Ny2min[i] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i))));
    Ny2size = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i)))) + 1;
    Space.Nz2min[i].resize(Ny2size);
    Space.CART_tb2Indvec[i].resize(Ny2size);
    for(int j = 0; j < Ny2size; ++j){
      Space.Nz2min[i][j] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i) - (Space.Ny2min[i]+j)*(Space.Ny2min[i]+j))));
      Nz2size = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i) - (Space.Ny2min[i]+j)*(Space.Ny2min[i]+j)))) + 1;
      Space.CART_tb2Indvec[i][j].resize(Nz2size);
      for(int k = 0; k < Nz2size; ++k){
	Space.CART_tb2Indvec[i][j][k] = count2;
	++count2;
      }
    }
  }
  Space.CART_tb2Indsize = count2;
  
  Space.Mmin = -2;
  Space.Msize = 3;
  Space.M2min = -2;
  Space.M2size = 3;
  if(Parameters.P != 0 && Parameters.N != 0){
    Space.Tmin = -2;
    Space.Tsize = 3;
    Space.T2min = -2;
    Space.T2size = 3;
  }
  else if(Parameters.P != 0 && Parameters.N == 0){
    Space.Tmin = -2;
    Space.Tsize = 1;
    Space.T2min = 0;
    Space.T2size = 1;
  }
  else if(Parameters.P == 0 && Parameters.N != 0){
    Space.Tmin = 2;
    Space.Tsize = 1;
    Space.T2min = 0;
    Space.T2size = 1;
  }

  return Space;
}*/

int ChanInd_1b(const std::string &basis, const Model_Space &Space, const State &State)
{
  // for tz = -1 only;
  if(basis == "finite_HO"){
    return State.n*Space.qsizes0.ml*2 + (State.ml - Space.qmins.ml)*2 + int(0.5*(State.m + 1));
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const State &State)
{
  if(basis == "infinite"){
    return Space.map_2b[(State.nx + 2*Space.nmax)*Space.qsizes.ny*Space.qsizes.nz + (State.ny + 2*Space.nmax)*Space.qsizes.nz +
			(State.nz + 2*Space.nmax)]*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t + 
                         int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_M"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.m*Space.qsizes.t + int((State.m - 2*Space.qmins.m)/2)*Space.qsizes.t + 
           int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_HO"){
    return int(State.ml - 2*Space.qmins.ml)*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t +
           int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_J"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.j*Space.qsizes.t + int(State.j/2)*Space.qsizes.t +
           int((State.t - 2*Space.qmins.t)/2);
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}


int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const State &State)
{
  if(basis == "infinite"){
    return Space.map_2b[(State.nx + 2*Space.nmax)*Space.qsizes.ny*Space.qsizes.nz + (State.ny + 2*Space.nmax)*Space.qsizes.nz +
			(State.nz + 2*Space.nmax)]*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t + 
                         int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_M"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.m*Space.qsizes.t + int((State.m - Space.qmins.m+Space.qmaxs.m)/2)*Space.qsizes.t + 
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_HO"){
    return int(State.ml - Space.qmins.ml+Space.qmaxs.ml)*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t +
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_J"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.j*Space.qsizes.t + int(State.j/2)*Space.qsizes.t +
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}


void Coulomb_Inf_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  
  double L = pow(Parameters.P/Parameters.density, 1.0/3.0);
  int i, j, k, l, a, b, c, d;
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  double TBME;

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    
    if(npp == 0){ goto stop1; }
    #pragma omp parallel private(a, b, c, d, TBME)
    {
      #pragma omp for schedule(static)
      for(int cd = 0; cd < npp; ++cd){
	c = Chan.ppvec[chan][2*cd];
	d = Chan.ppvec[chan][2*cd + 1];
	if(c == d){ continue; }
	for(int ab = cd; ab < npp; ++ab){
	  a = Chan.ppvec[chan][2*ab];
	  b = Chan.ppvec[chan][2*ab + 1];
	  if(a == b){ continue; }
	  TBME = Coulomb_Inf(Space, a, b, c, d, L);
	  Ints.D_ME1.V1[chan][npp*cd + ab] = TBME;
	  Ints.D_ME1.V1[chan][npp*ab + cd] = TBME;
	}
      }
    }
  stop1:;
    if(nhh == 0){ goto stop2; }
    #pragma omp parallel private(i, j, k, l, TBME)
    {
      #pragma omp for schedule(static)
      for(int ij = 0; ij < nhh; ++ij){
	i = Chan.hhvec[i][2*ij];
	j = Chan.hhvec[i][2*ij + 1];
	if(i == j){ continue; }
	for(int kl = ij; kl < nhh; ++kl){
	  k = Chan.hhvec[i][2*kl];
	  l = Chan.hhvec[i][2*kl + 1];
	  if(k == l){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, k, l, L);
	  Ints.D_ME1.V2[i][nhh*ij + kl] = TBME;
	  Ints.D_ME1.V2[i][nhh*kl + ij] = TBME;
	}
      }
    }
  stop2:;
    if(nhh * npp == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int ab = 0; ab < npp; ++ab){
	a = Chan.ppvec[i][2*ab];
	b = Chan.ppvec[i][2*ab + 1];
	if(a == b){ continue; }
	for(int ij = 0; ij < nhh; ++ij){
	  i = Chan.hhvec[i][2*ij];
	  j = Chan.hhvec[i][2*ij + 1];
	  if(i == j){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V4[i][nhh*ab + ij] = TBME;
	}
      }
    }
  }
  
  for(int chan = 0; chan < Chan.size3; ++chan){
    nh = Chan.nh[chan];
    np = Chan.np[chan];
    nhpp = Chan.nhpp[chan];
    nhhp = Chan.nhhp[chan];

    if(nh * nhpp == 0){ goto stop3; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int h = 0; h < nh; ++h){
	j = Chan.hvec[chan][h];
	for(int hpp = 0; hpp < nhpp; ++hpp){
	  i = Chan.hppvec[chan][3*hpp];
	  a = Chan.hppvec[chan][3*hpp + 1];
	  b = Chan.hppvec[chan][3*hpp + 2];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V5[chan][nhpp*h + hpp] = TBME;
	  Ints.D_ME1.V6[chan][nhpp*h + hpp] = -1.0 * TBME;
	}
      }
    }
  stop3:;
    if(np * nhhp == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int p = 0; p < np; ++p){
	b = Chan.pvec[chan][p];
	for(int hhp = 0; hhp < nhhp; ++hhp){
	  i = Chan.hhpvec[chan][3*hhp];
	  j = Chan.hhpvec[chan][3*hhp + 1];
	  a = Chan.hhpvec[chan][3*hhp + 2];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V7[chan][nhhp*p + hhp] = TBME;
	  Ints.D_ME1.V8[chan][nhhp*p + hhp] = -1.0 * TBME;
	}
      }
    }
  }
  
  for(int chan = 0; chan < Chan.size2; ++chan){
    nhp1 = Chan.nhp1[chan];
    nhp2 = Chan.nhp2[chan];

    if(nhp2 == 0){ goto stop4; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int hp21 = 0; hp21 < nhp2; ++hp21){
	i = Chan.hp2vec[chan][2*hp21];
	b = Chan.hp2vec[chan][2*hp21 + 1];
	for(int hp22 = hp21; hp22 < nhp2; ++hp22){
	  j = Chan.hp2vec[chan][2*hp22];
	  a = Chan.hp2vec[chan][2*hp22 + 1];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, a, j, b, L);
	  Ints.D_ME1.V3[chan][nhp2*hp21 + hp22] = TBME;
	  Ints.D_ME1.V3[chan][nhp2*hp22 + hp21] = TBME;
	}
      }
    }
  stop4:;
    if(nhp1 * nhp2 == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int hp2 = 0; hp2 < nhp2; ++hp2){
	i = Chan.hp2vec[chan][2*hp2];
	a = Chan.hp2vec[chan][2*hp2 + 1];
	for(int hp1 = 0; hp1 < Chan.nhp1[chan]; ++hp1){
	  j = Chan.hp1vec[chan][2*hp1];
	  b = Chan.hp1vec[chan][2*hp1 + 1];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V9[chan][nhp1*hp2 + hp1] = TBME;
	  Ints.D_ME1.V10[chan][nhp1*hp2 + hp1] = -1.0 * TBME;
	}
      }
    }
  }
}


// Minnesota Potential for momentum basis
double Coulomb_Inf(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
{
  double term = 0.0;
  double e_sq = hbarc_HartA * fine_struct;
  double prefactor = e_sq/(L*L*L);
  double qSquared1;
  double kX1, kY1, kZ1;

  if(Space.qnums[qi].nx + Space.qnums[qj].nx != Space.qnums[qk].nx + Space.qnums[ql].nx || 
     Space.qnums[qi].ny + Space.qnums[qj].ny != Space.qnums[qk].ny + Space.qnums[ql].ny || 
     Space.qnums[qi].nz + Space.qnums[qj].nz != Space.qnums[qk].nz + Space.qnums[ql].nz){ return 0.0; }
  if(Space.qnums[qi].m == Space.qnums[qk].m && Space.qnums[qj].m == Space.qnums[ql].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qk].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qk].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qk].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-8){ term += 0.0; }
    else{ term += prefactor*4.0*M_PI/qSquared1; }
  }
  if(Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[ql].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[ql].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[ql].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-8){ term -= 0.0; }
    else{ term -= prefactor*4.0*M_PI/qSquared1; }
  }
  return term;
}


// Minnesota Potential for momentum basis
double Coulomb_HO(const Input_Parameters &Parameters, const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql)
{
  int g1, g2, g3, g4, G, L;
  double dir = 0.0;
  double exch = 0.0;
  int n1, m1, n2, m2, n3, m3, n4, m4;
  double M11, M21, M31, M41, M12, M22, M32, M42;
  double temp;
  double LogRatio1, LogProd2, LogRatio2, ExProd, ExTerm;
  double ExTerm2, LogProd3, LGam1, LGam2;
  n1 = Space.qnums[qi].n; //1
  m1 = Space.qnums[qi].ml;
  M11 = 0.5*(abs(m1) + m1);
  M12 = 0.5*(abs(m1) - m1);
  n2 = Space.qnums[qj].n; //2
  m2 = Space.qnums[qj].ml;
  M21 = 0.5*(abs(m2) + m2);
  M22 = 0.5*(abs(m2) - m2);
  n3 = Space.qnums[ql].n; //4
  m3 = Space.qnums[ql].ml;
  M31 = 0.5*(abs(m3) + m3);
  M32 = 0.5*(abs(m3) - m3);
  n4 = Space.qnums[qk].n; //3
  m4 = Space.qnums[qk].ml;
  M41 = 0.5*(abs(m4) + m4);
  M42 = 0.5*(abs(m4) - m4);
  if(Space.qnums[qi].m == Space.qnums[qk].m && Space.qnums[qj].m == Space.qnums[ql].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j4 = 0; j4 <= n4; ++j4){
	g1 = int(j1 + j4 + M11 + M42);
	g4 = int(j4 + j1 + M41 + M12);
	for(int j2 = 0; j2 <= n2; ++j2){
	  for(int j3 = 0; j3 <= n3; ++j3){
	    g2 = int(j2 + j3 + M21 + M32);
	    g3 = int(j3 + j2 + M31 + M22);
	    G = g1 + g2 + g3 + g4;
	    LogRatio1 = logratio1(j1, j2, j3, j4);
	    LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    LogRatio2 = logratio2(G);
	    ExTerm = (-2*((j1 + j2 + j3 + j4)%2) + 1);
	    ExProd = exp(LogRatio1 + LogProd2 + LogRatio2);
	    temp = 0.0;
	    for(int l2 = 0; l2 <= g2; ++l2){
	      for(int l3 = 0; l3 <= g3; ++l3){
		ExTerm2 = (-2*((g2 + g3 - l2 - l3)%2) + 1);
		for(int l1 = 0; l1 <= g1; ++l1){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 == l3 + l4){
		      L = l1 + l2 + l3 + l4;
		      LogProd3 = logproduct3(l1, l2, l3, l4, g1, g2, g3, g4);
		      LGam1 = lgamma(1.0 + 0.5*L);
		      LGam2 = lgamma(0.5*(G - L + 1.0));
		      temp += ExTerm2 * exp(LogProd3 + LGam1 + LGam2);
		    }
		  }
		}
	      }
	    }
	    dir += ExTerm * ExProd * temp;
	  }
	}
      }
    }
    dir *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }

  n3 = Space.qnums[qk].n; //3
  m3 = Space.qnums[qk].ml;
  M31 = 0.5*(abs(m3) + m3);
  M32 = 0.5*(abs(m3) - m3);
  n4 = Space.qnums[ql].n; //4
  m4 = Space.qnums[ql].ml;
  M41 = 0.5*(abs(m4) + m4);
  M42 = 0.5*(abs(m4) - m4);
  if(Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j4 = 0; j4 <= n4; ++j4){
	g1 = int(j1 + j4 + M11 + M42);
	g4 = int(j4 + j1 + M41 + M12);
	for(int j2 = 0; j2 <= n2; ++j2){
	  for(int j3 = 0; j3 <= n3; ++j3){
	    g2 = int(j2 + j3 + M21 + M32);
	    g3 = int(j3 + j2 + M31 + M22);
	    G = g1 + g2 + g3 + g4;
	    LogRatio1 = logratio1(j1, j2, j3, j4);
	    LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    LogRatio2 = logratio2(G);
	    ExTerm = (-2*((j1 + j2 + j3 + j4)%2) + 1);
	    ExProd = exp(LogRatio1 + LogProd2 + LogRatio2);
	    temp = 0.0;
	    for(int l2 = 0; l2 <= g2; ++l2){
	      for(int l3 = 0; l3 <= g3; ++l3){
		ExTerm2 = (-2*((g2 + g3 - l2 - l3)%2) + 1);
		for(int l1 = 0; l1 <= g1; ++l1){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 == l3 + l4){
		      L = l1 + l2 + l3 + l4;
		      LogProd3 = logproduct3(l1, l2, l3, l4, g1, g2, g3, g4);
		      LGam1 = lgamma(1.0 + 0.5*L);
		      LGam2 = lgamma(0.5*(G - L + 1.0));
		      temp += ExTerm2 * exp(LogProd3 + LGam1 + LGam2);
		    }
		  }
		}
	      }
	    }
	    exch += ExTerm * ExProd * temp;
	  }
	}
      }
    }
    exch *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }
  return std::sqrt(Parameters.density)*(dir - exch);
}


void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);

  for(int i = 0; i < Chan.size1; ++i){
    for(int pq = 0; pq < Chan.npp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec[i][2*pq];
      shell2 = Chan.ppvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = pq; rs < Chan.npp[i]; ++rs){
	shell3 = Chan.ppvec[i][2*rs];
	shell4 = Chan.ppvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V1[i][Chan.npp[i]*pq + rs] = TBME;
	Ints.D_ME1.V1[i][Chan.npp[i]*rs + pq] = TBME;
      }
    }
    for(int pq = 0; pq < Chan.nhh[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hhvec[i][2*pq];
      shell2 = Chan.hhvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = pq; rs < Chan.nhh[i]; ++rs){
	shell3 = Chan.hhvec[i][2*rs];
	shell4 = Chan.hhvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V2[i][Chan.nhh[i]*pq + rs] = TBME;
	Ints.D_ME1.V2[i][Chan.nhh[i]*rs + pq] = TBME;
      }
    }
    for(int pq = 0; pq < Chan.npp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec[i][2*pq];
      shell2 = Chan.ppvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = 0; rs < Chan.nhh[i]; ++rs){
	shell3 = Chan.hhvec[i][2*rs];
	shell4 = Chan.hhvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V4[i][Chan.nhh[i]*pq + rs] = TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    for(int q = 0; q < Chan.nh[i]; ++q){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell2 = Chan.hvec[i][q];
      for(int prs = 0; prs < Chan.nhpp[i]; ++prs){
	shell1 = Chan.hppvec[i][3*prs];
	shell3 = Chan.hppvec[i][3*prs + 1];
	shell4 = Chan.hppvec[i][3*prs + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V5[i][Chan.nhpp[i]*q + prs] = TBME;
	Ints.D_ME1.V6[i][Chan.nhpp[i]*q + prs] = -1.0 * TBME;
      }
    }
    for(int s = 0; s < Chan.np[i]; ++s){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell4 = Chan.pvec[i][s];
      for(int pqr = 0; pqr < Chan.nhhp[i]; ++pqr){
	shell1 = Chan.hhpvec[i][3*pqr];
	shell2 = Chan.hhpvec[i][3*pqr + 1];
	shell3 = Chan.hhpvec[i][3*pqr + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V7[i][Chan.nhhp[i]*s + pqr] = TBME;
	Ints.D_ME1.V8[i][Chan.nhhp[i]*s + pqr] = -1.0 * TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size2; ++i){
    for(int ps = 0; ps < Chan.nhp2[i]; ++ps){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec[i][2*ps];
      shell4 = Chan.hp2vec[i][2*ps + 1];
      for(int qr = ps; qr < Chan.nhp2[i]; ++qr){
	shell3 = Chan.hp2vec[i][2*qr];
	shell2 = Chan.hp2vec[i][2*qr + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V3[i][Chan.nhp2[i]*ps + qr] = TBME;
	Ints.D_ME1.V3[i][Chan.nhp2[i]*qr + ps] = TBME;
      }
    }
    for(int pr = 0; pr < Chan.nhp2[i]; ++pr){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec[i][2*pr];
      shell3 = Chan.hp2vec[i][2*pr + 1];
      for(int qs = 0; qs < Chan.nhp1[i]; ++qs){
	shell2 = Chan.hp1vec[i][2*qs];
	shell4 = Chan.hp1vec[i][2*qs + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell4, shell3, L);
	Ints.D_ME1.V9[i][Chan.nhp1[i]*pr + qs] = TBME;
	Ints.D_ME1.V10[i][Chan.nhp1[i]*pr + qs] = -1.0 * TBME;
      }
    }
  }
}

int kron_del(const int &i, const int &j)
{
  if(i != j){ return 0; }
  return 1;
}

// Minnesota Potential for momentum basis
double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
{
  double V_R1, V_T1, V_S1, V_R2, V_T2, V_S2;
  double V_0R, V_0T, V_0S;
  double kappa_R, kappa_T, kappa_S;
  double kX1, kY1, kZ1, kX2, kY2, kZ2;
  double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
  double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;
  V_0R = 200; //MeV
  V_0T = 178; //MeV
  V_0S = 91.85; //MeV
  kappa_R = 1.487; //fm^-2
  kappa_T = 0.639; //fm^-2
  kappa_S = 0.465; //fm^-2

  if(Space.qnums[qi].nx + Space.qnums[qj].nx != Space.qnums[qk].nx + Space.qnums[ql].nx){ return 0.0; }
  if(Space.qnums[qi].ny + Space.qnums[qj].ny != Space.qnums[qk].ny + Space.qnums[ql].ny){ return 0.0; }
  if(Space.qnums[qi].nz + Space.qnums[qj].nz != Space.qnums[qk].nz + Space.qnums[ql].nz){ return 0.0; }
  if(Space.qnums[qi].m + Space.qnums[qj].m != Space.qnums[qk].m + Space.qnums[ql].m){ return 0.0; }
  if(Space.qnums[qi].t + Space.qnums[qj].t != Space.qnums[qk].t + Space.qnums[ql].t){ return 0.0; }

  kX1 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[qk].nx + Space.qnums[ql].nx);
  kY1 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[qk].ny + Space.qnums[ql].ny);
  kZ1 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[qk].nz + Space.qnums[ql].nz);

  kX2 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[ql].nx + Space.qnums[qk].nx);
  kY2 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[ql].ny + Space.qnums[qk].ny);
  kZ2 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[ql].nz + Space.qnums[qk].nz);

  qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
  qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;
  
  V_R1 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared1/(4*kappa_R));
  V_T1 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared1/(4*kappa_T));
  V_S1 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared1/(4*kappa_S));

  V_R2 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared2/(4*kappa_R));
  V_T2 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared2/(4*kappa_T));
  V_S2 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared2/(4*kappa_S));
  
  spinEx1 = spinExchangeMtxEle(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[qk].m, Space.qnums[ql].m);
  isoSpinEx1 = spinExchangeMtxEle(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[qk].t, Space.qnums[ql].t);

  spinEx2 = spinExchangeMtxEle(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[ql].m, Space.qnums[qk].m);
  isoSpinEx2 = spinExchangeMtxEle(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[ql].t, Space.qnums[qk].t);
  
  IsIt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m) * kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsIt1 = spinEx1 * kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m)*kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * isoSpinEx1;

  IsIt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsIt2 = spinEx2 * kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * isoSpinEx2;

  return 0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * IsIt1 + 
    0.25 * (V_T1 - V_S1) * PsIt1 - 
    0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * PsPt1 - 
    0.25 * (V_T1 - V_S1) * IsPt1 -
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * IsIt2 - 
    0.25 * (V_T2 - V_S2) * PsIt2 +
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * PsPt2 + 
    0.25 * (V_T2 - V_S2) * IsPt2;
}

int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l)
{
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}


/*double vint_Conversion(const Model_Space &Space, const std::vector<std::vector<double> > &ME, const int &qi, const int &qj, const int &qk, const int &ql)
{
  int kx1, ky1, kz1, kx2, ky2, kz2, Tz, Sz, ind1, ind2, ind3, MEind;
  double TBME = 0.0, CGC10, CGC11, CGC20, CGC21;
  kx1 = Space.levelsnx[qi] - Space.levelsnx[qj];
  ky1 = Space.levelsny[qi] - Space.levelsny[qj];
  kz1 = Space.levelsnz[qi] - Space.levelsnz[qj];
  kx2 = Space.levelsnx[qk] - Space.levelsnx[ql];
  ky2 = Space.levelsny[qk] - Space.levelsny[ql];
  kz2 = Space.levelsnz[qk] - Space.levelsnz[ql];
  Tz = int(0.5*Space.levelst[qi] + 0.5*Space.levelst[qj]);
  Sz = int(0.5*Space.levelsm[qi] + 0.5*Space.levelsm[qj]);
  MEind = (Sz+1)*3 + (Tz+1);
  ind1 = CART_tbInd2_rel(Space, kx1, ky1, kz1);
  ind2 = CART_tbInd2_rel(Space, kx2, ky2, kz2);
  CGC10 = CGC(0.5,0.5*Space.levelsm[qi],0.5,0.5*Space.levelsm[qj],0,Sz);
  CGC11 = CGC(0.5,0.5*Space.levelsm[qi],0.5,0.5*Space.levelsm[qj],1,Sz);
  CGC20 = CGC(0.5,0.5*Space.levelsm[qk],0.5,0.5*Space.levelsm[ql],0,Sz);
  CGC21 = CGC(0.5,0.5*Space.levelsm[qk],0.5,0.5*Space.levelsm[ql],1,Sz);
  if(Sz == 0){
    ind3 = (2*Space.CART_tb2Indsize)*ind1 + ind2;
    TBME += ME[MEind][ind3]*CGC10*CGC20; //S=00
    ind3 = (2*Space.CART_tb2Indsize)*ind1 + (Space.CART_tb2Indsize+ind2);
    TBME += ME[MEind][ind3]*CGC10*CGC21; //S=01
    ind3 = (2*Space.CART_tb2Indsize)*(Space.CART_tb2Indsize+ind1) + ind2;
    TBME += ME[MEind][ind3]*CGC11*CGC20; //S=10
    ind3 = (2*Space.CART_tb2Indsize)*(Space.CART_tb2Indsize+ind1) + (Space.CART_tb2Indsize+ind2);
    TBME += ME[MEind][ind3]*CGC11*CGC21; //S=11
  }
  else{
    ind3 = Space.CART_tb2Indsize*ind1 + ind2;
    TBME += ME[MEind][ind3]*CGC11*CGC21; //S=11
  }
  std::cout << "** " << kx1 << " " << ky1 << " " << kz1 << " " << kx2 << " " << ky2 << " " << kz2 << " " << Sz << " " << Tz << " " << TBME << std::endl;
  return TBME;
  }*/

void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps)
{
  std::cout << "Calculating Coupled-Cluster Amplitudes ... " << std::endl;
  
  int ind = 0;
  int i, j, a, b;
  int nhh, npp;
  double tempen, tempt;
  double error0 = 1000.0;
  double error1;
  double error2;
  double total1;
  double total2;
  Amplitudes Amps2 = Amplitudes(Parameters, Space, Chan);

  double CCinE, CCoutE;
  int Stot = 0;
  int Dtot = 0;

  ////////////////////////////////////////////////////////////////

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }

    #pragma omp parallel private(i, j, a, b, tempen, tempt)
    {
      #pragma omp for schedule(static)
      for(int pp = 0; pp < npp; ++pp){
	a = Chan.ppvec[chan][2*pp];
	b = Chan.ppvec[chan][2*pp + 1];
	for(int hh = 0; hh < nhh; ++hh){
	  i = Chan.hhvec[chan][2*hh];
	  j = Chan.hhvec[chan][2*hh + 1];
	  if(i == j || a == b){ continue; }
	  ++Dtot;
	  tempen = Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy;
	  Amps.D1.Evec[chan][hh * npp + pp] = tempen;
	  Amps2.D1.Evec[chan][hh * npp + pp] = tempen;
	  tempt = Ints.D_ME1.V4[chan][pp * nhh + hh] / tempen;
	  Amps.D1.set_T(chan, hh * npp + pp, tempt);
	  //std::cout << "T: " << i << " " << j << " " << a << " " << b << " = " << tempt << " , " << tempen << std::endl;
	  /*if(ind1 == ind2 || ind3 == ind4){ tempt = 0.0; }
	    else if(ind1 < ind2 && ind3 < ind4){ tempt = 1000*ind3 + 100*ind4 + 10*ind1 + ind2; }
	    else if(ind1 > ind2 && ind3 > ind4){ tempt = 1000*ind4 + 100*ind3 + 10*ind2 + ind1; }
	    else if(ind1 < ind2 && ind3 > ind4){ tempt = -1.0*(1000*ind4 + 100*ind3 + 10*ind1 + ind2); }
	    else if(ind1 > ind2 && ind3 < ind4){ tempt = -1.0*(1000*ind3 + 100*ind4 + 10*ind2 + ind1); }
	    Amps.D1.set_T(chan, hind * Chan.pp[chan] + pind, tempt);*/
	}
      }
    }
  }

  if(Parameters.approx == "singles"){
    for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
      ++Stot;
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempen = Space.qnums[i].energy - Space.qnums[a].energy;
      tempt = 0.0;
      //tempt = 10*ind2 + ind1;
      //tempt = 0.01;
      //std::cout << "t: " << i << " " << a << " = " << tempt << " " << tempen << std::endl;
      Amps.S1.Evec[hp1] = tempen;
      Amps2.S1.Evec[hp1] = tempen;
      Amps.S1.set_T(hp1, tempt);
    }
    Amps.S1.set_T_2(Chan, Ints);
    Amps.D1.set_T_2(Chan, Ints);
  }

  CCoutE = Amps.get_energy(Parameters, Chan, Ints);
  if(Parameters.approx == "singles"){  
    std::cout << "Iteration Number = " << ind << ", CCSD Energy = " << CCoutE << std::endl;
  }
  else{
    std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCoutE << std::endl;
  }

  //std::cout << std::endl;
  while((error0 > 1e-12 && ind < 5000) || ind < 40){
    Doubles_Step(Space, Chan, Ints, Amps, Amps2);
    if(Parameters.approx == "singles"){
      Doubles_Step_2(Space, Chan, Ints, Amps, Amps2);
      Singles_Step(Space, Chan, Ints, Amps, Amps2);
    }

    error1 = 0.0;
    error2 = 0.0;
    total1 = 0.0;
    total2 = 0.0;
    for(int chan = 0; chan < Chan.size1; ++chan){
      nhh = Chan.nhh[chan];
      npp = Chan.npp[chan];
      if(nhh * npp == 0){ continue; }

      #pragma omp parallel private(i, j, a, b, tempt)
      {
        #pragma omp for schedule(static) reduction(+:error2, total2)
	for(int pp = 0; pp < npp; ++pp){
	  a = Chan.ppvec[chan][2*pp];
	  b = Chan.ppvec[chan][2*pp + 1];
	  for(int hh = 0; hh < nhh; ++hh){
	    i = Chan.hhvec[chan][2*hh];
	    j = Chan.hhvec[chan][2*hh + 1];
	    if(i == j || a == b){ continue; }
	    tempt = Amps2.D1.get_T(chan, hh * npp + pp);
	    tempt += Ints.D_ME1.V4[chan][pp * nhh + hh];
	    tempt /= Amps2.D1.Evec[chan][hh * npp + pp];
	    if(fabs(0.45*tempt + 0.55*Amps.D1.T1[chan][hh * npp + pp]) > 1e-10){
	      error2 += fabs((tempt - Amps.D1.T1[chan][hh * npp + pp])/(0.45*tempt + 0.55*Amps.D1.T1[chan][hh * npp + pp]));
	    }
	    Amps.D1.set_T(chan, hh * npp + pp, 0.45*tempt + 0.55*Amps.D1.T1[chan][hh * npp + pp]);
	    //std::cout << "T: " << i << " " << j << " " << a << " " << b << " = " << tempt << std::endl;
	  }
	}
      }
    }
    error2 /= Dtot;

    if(Parameters.approx == "singles" && Stot != 0){
      //std::cout << std::endl;
      for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
	i = Chan.hp1vec[Chan.ind0][2*hp1];
	a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	tempt = Amps2.S1.get_T(hp1);
	tempt /= Amps2.S1.Evec[hp1];
	if(fabs(0.45*tempt + 0.55*Amps.S1.T1[hp1]) > 1e-10){
	  error1 += fabs((tempt - Amps.S1.T1[hp1])/(0.45*tempt + 0.55*Amps.S1.T1[hp1]));
	}
	Amps.S1.set_T(hp1, 0.45*tempt + 0.55*Amps.S1.T1[hp1]);
	//std::cout << "t: " << i << " " << a << " = " << tempt << std::endl;
      }
      Amps.S1.set_T_2(Chan, Ints);
      Amps.D1.set_T_2(Chan, Ints);
      error1 /= Stot;
    }
    
    error0 = error1 + error2;
    CCoutE = Amps.get_energy(Parameters, Chan, Ints);
    CCinE = CCoutE;
    Amps2.zero(Parameters, Chan);

    //std::cout << std::endl;
    if(Parameters.approx == "singles"){
      std::cout << "Iteration Number = " << ind+1 << ", CCSD Energy = " << CCoutE << ", error = " << error0 << "\r"; // << std::endl;
    }
    else{
      std::cout << "Iteration Number = " << ind+1 << ", CCD Energy = " << CCoutE << ", error = " << error0 << "\r"; // << std::endl;      
    }
    ++ind;
  }

  std::cout << std::endl << std::endl;
  Amps2.delete_struct(Parameters, Chan);
}

double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Ints)
{
  double energy = 0.0;
  int ind, ind1;
  State tb;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].type != "hole"){ continue; }
    energy += Space.qnums[i].energy;
    for(int j = 0; j < Space.indtot; ++j){
      if(Space.qnums[j].type != "hole" || i == j){ continue; }
      plus(tb, Space.qnums[i], Space.qnums[j]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
      ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], i, j, i, j);
      energy -= 0.5 * Ints.D_ME1.V2[ind1][ind];
    }
  }
  return energy;
}

void Doubles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac4 = 0.25, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.nhh[chan];
    int pp = Chan.npp[chan];
    if(hh == 0 || pp == 0){ continue; }
    /*std::cout << "hh" << std::endl;
    for(int i = 0; i < hh; ++i){ std::cout << Chan.hhvec[chan][2*i] << Chan.hhvec[chan][2*i + 1] << " "; }
    std::cout << std::endl;
    std::cout << "pp" << std::endl;
    for(int i = 0; i < pp; ++i){ std::cout << Chan.ppvec[chan][2*i] << Chan.ppvec[chan][2*i + 1] << " "; }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = 0.5 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V1[chan], Amps2.D1.T1[chan], &hh, &pp, &pp, &fac3, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (1)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = 0.5 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
    dgemm_NN(Ints.D_ME1.V2[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac3, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (2)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = 0.25 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac4, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (7)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan];
    int hp2 = Chan.nhp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    /*std::cout << "hp1" << std::endl;
    for(int i = 0; i < hp1; ++i){ std::cout << Chan.hp1vec[chan][2*i] << Chan.hp1vec[chan][2*i + 1] << " "; }
    std::cout << std::endl;
    std::cout << "hp2" << std::endl;
    for(int i = 0; i < hp2; ++i){ std::cout << Chan.hp2vec[chan][2*i] << Chan.hp2vec[chan][2*i + 1] << " "; }
    std::cout << std::endl;*/
    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3)
    dgemm_NN(Amps1.D1.T2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = -T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4)
    //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5)
    //T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6)
    for(int i = 0; i < hp1 * hp2; ++i){
      Amps2.D1.T3[chan][i] = Amps2.D1.T2[chan][i];
      Amps2.D1.T4[chan][i] = -1.0 * Amps2.D1.T2[chan][i];
      Amps2.D1.T5[chan][i] = -1.0 * Amps2.D1.T2[chan][i];
    }
    /*std::cout << "T2 (3)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T2[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "T3 (4)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T3[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "T4 (5)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T4[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "T5 (6)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T5[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.D1.T2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T2 (12)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T2[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.D1.T2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T4 (13)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T4[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    /*std::cout << "h" << std::endl;
    for(int i = 0; i < h; ++i){ std::cout << Chan.hvec[chan][i] << " "; }
    std::cout << std::endl;
    std::cout << "hpp" << std::endl;
    for(int i = 0; i < hpp; ++i){ std::cout << Chan.hppvec[chan][3*i] << Chan.hppvec[chan][3*i + 1] << Chan.hppvec[chan][3*i + 2] << " "; }
    std::cout << std::endl;*/
    //T6(ab|ij){jab,i} = -0.5 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8)
    dgemm_NN(Ints.D_ME1.V5[chan], Amps1.D1.T6[chan], Amps2.D1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S2[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac6, &fac1, &N, &N); //fac1
    /*std::cout << "T6 (8)" << std::endl;
    for(int i = 0; i < hpp; ++i){
      for(int j = 0; j < h; ++j){
	std::cout << Amps2.D1.T6[chan][h*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T7(ab|ij){iab,j} = -0.5 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9)
    dgemm_NN(Ints.D_ME1.V6[chan], Amps1.D1.T6[chan], Amps2.D1.S3[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S3[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac6, &fac1, &N, &N); //fac1
    /*std::cout << "T7 (9)" << std::endl;
    for(int i = 0; i < hpp; ++i){
      for(int j = 0; j < h; ++j){
	std::cout << Amps2.D1.T7[chan][h*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    /*std::cout << "p" << std::endl;
    for(int i = 0; i < p; ++i){ std::cout << Chan.pvec[chan][i] << " "; }
    std::cout << std::endl;
    std::cout << "hhp" << std::endl;
    for(int i = 0; i < hhp; ++i){ std::cout << Chan.hhpvec[chan][3*i] << Chan.hhpvec[chan][3*i + 1] << Chan.hhpvec[chan][3*i + 2] << " "; }
    std::cout << std::endl;*/
    //T8(ab|ij){ijb,a} = -0.5 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10)
    dgemm_NN(Ints.D_ME1.V7[chan], Amps1.D1.T8[chan], Amps2.D1.S4[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S4[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac6, &fac1, &N, &N); //fac1
    /*std::cout << "T8 (10)" << std::endl;
    for(int i = 0; i < hhp; ++i){
      for(int j = 0; j < p; ++j){
	std::cout << Amps2.D1.T8[chan][p*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T9(ab|ij){ija,b} = -0.5 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11)
    dgemm_NN(Ints.D_ME1.V8[chan], Amps1.D1.T8[chan], Amps2.D1.S5[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S5[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac6, &fac1, &N, &N); //fac1
    /*std::cout << "T9 (11)" << std::endl;
    for(int i = 0; i < hhp; ++i){
      for(int j = 0; j < p; ++j){
	std::cout << Amps2.D1.T9[chan][p*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
  }
}

void Singles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = -1.0, fac4 = -0.5;
  char N = 'N';
  int hp = Chan.nhp1[Chan.ind0];
  int one = 1;
  if(hp != 0){
    /*std::cout << "hp(ind0)" << std::endl;
    for(int i = 0; i < hp; ++i){ std::cout << Chan.hp1vec[Chan.ind0][2*i] << Chan.hp1vec[Chan.ind0][2*i + 1] << " "; }
    std::cout << std::endl;*/
    //t1(ia){ia} = -V3(ka|ic){ia,kc}.t1(kc){kc} (1)
    dgemm_NN(Ints.D_ME1.V3[Chan.ind0], Amps1.S1.T1, Amps2.S1.T1, &hp, &one, &hp, &fac3, &fac1, &N, &N); //fac1
    /*std::cout << "t1 (1)" << std::endl;
    for(int i = 0; i < hp; ++i){
      std::cout << Amps2.S1.T1[i] << " ";
    }
    std::cout << std::endl;*/
    //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6)
    dgemm_NN(Ints.D_ME1.V9[Chan.ind0], Amps1.D1.T2[Chan.ind0], Amps2.S1.S3, &hp, &hp, &hp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.T1, Amps2.S1.S3, Amps2.S1.T1, &one, &hp, &hp, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "t1 (6)" << std::endl;
    for(int i = 0; i < hp; ++i){
      std::cout << Amps2.S1.T1[i] << " ";
    }
    std::cout << std::endl;*/
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    /*std::cout << "h" << std::endl;
    for(int i = 0; i < h; ++i){ std::cout << Chan.hvec[chan][i] << " "; }
    std::cout << std::endl;
    std::cout << "hpp" << std::endl;
    for(int i = 0; i < hpp; ++i){ std::cout << Chan.hppvec[chan][3*i] << Chan.hppvec[chan][3*i + 1] << Chan.hppvec[chan][3*i + 2] << " "; }
    std::cout << std::endl;
    std::cout << "p" << std::endl;
    for(int i = 0; i < p; ++i){ std::cout << Chan.pvec[chan][i] << " "; }
    std::cout << std::endl;
    std::cout << "hhp" << std::endl;
    for(int i = 0; i < hhp; ++i){ std::cout << Chan.hhpvec[chan][3*i] << Chan.hhpvec[chan][3*i + 1] << Chan.hhpvec[chan][3*i + 2] << " "; }
    std::cout << std::endl;*/
    if(p != 0 && h != 0){
      if(hpp != 0){
	//t2(a|i){a,i} = -0.5 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2)
	dgemm_NN(Ints.S_ME1.V11[chan], Amps1.D1.T6[chan], Amps2.S1.T2[chan], &p, &h, &hpp, &fac4, &fac1, &N, &N); //fac1
	/*std::cout << "t2 (2)" << std::endl;
	for(int i = 0; i < p; ++i){
	  for(int j = 0; j < h; ++j){
	    std::cout << Amps2.S1.T2[chan][h*i + j] << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;*/
	//t2(a|i){a,i} = V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} (3)
	dgemm_NN(Ints.S_ME1.V11[chan], Amps1.S1.E6[chan], Amps2.S1.T2[chan], &p, &h, &hpp, &fac1, &fac1, &N, &N); //fac1
	/*std::cout << "t2 (3)" << std::endl;
	for(int i = 0; i < p; ++i){
	  for(int j = 0; j < h; ++j){
	    std::cout << Amps2.S1.T2[chan][h*i + j] << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;*/
	//t2(a|i){a,i} = -0.5 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7)
	dgemm_NN(Ints.D_ME1.V6[chan], Amps1.D1.T6[chan], Amps2.S1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T2[chan], Amps2.S1.S2[chan], Amps2.S1.T2[chan], &p, &h, &h, &fac4, &fac1, &N, &N); //fac1
	/*std::cout << "t2 (7)" << std::endl;
	for(int i = 0; i < p; ++i){
	  for(int j = 0; j < h; ++j){
	    std::cout << Amps2.S1.T2[chan][h*i + j] << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;*/
      }
      if(hhp != 0){
	//t3(i|a){i,a} = -0.5 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4)
	dgemm_NN(Ints.S_ME1.V12[chan], Amps1.D1.T8[chan], Amps2.S1.T3[chan], &h, &p, &hhp, &fac4, &fac1, &N, &N); //fac1
 	/*std::cout << "t3 (4)" << std::endl;
	for(int i = 0; i < h; ++i){
	  for(int j = 0; j < p; ++j){
	    std::cout << Amps2.S1.T3[chan][p*i + j] << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;*/
	//t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} (5)
	dgemm_NN(Ints.S_ME1.V12[chan], Amps1.S1.E8[chan], Amps2.S1.T3[chan], &h, &p, &hhp, &fac1, &fac1, &N, &N); //fac1
 	/*std::cout << "t3 (5)" << std::endl;
	for(int i = 0; i < h; ++i){
	  for(int j = 0; j < p; ++j){
	    std::cout << Amps2.S1.T3[chan][p*i + j] << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;*/
	//t3(i|a){i,a} = -0.5 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8)
	dgemm_NN(Ints.D_ME1.V8[chan], Amps1.D1.T8[chan], Amps2.S1.S1[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T3[chan], Amps2.S1.S1[chan], Amps2.S1.T3[chan], &h, &p, &p, &fac4, &fac1, &N, &N); //fac1
 	/*std::cout << "t3 (8)" << std::endl;
	for(int i = 0; i < h; ++i){
	  for(int j = 0; j < p; ++j){
	    std::cout << Amps2.S1.T3[chan][p*i + j] << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;*/
	//t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} (9)
	dgemm_NN(Ints.D_ME1.V8[chan], Amps1.S1.E8[chan], Amps2.S1.S1[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T3[chan], Amps2.S1.S1[chan], Amps2.S1.T3[chan], &h, &p, &p, &fac1, &fac1, &N, &N); //fac1
 	/*std::cout << "t3 (9)" << std::endl;
	for(int i = 0; i < h; ++i){
	  for(int j = 0; j < p; ++j){
	    std::cout << Amps2.S1.T3[chan][p*i + j] << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;*/
      }
    }
  }
}
 
void Doubles_Step_2(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.nhh[chan];
    int pp = Chan.npp[chan];
    if(hh * pp == 0){ continue; }
    /*std::cout << "hh" << std::endl;
    for(int i = 0; i < hh; ++i){ std::cout << Chan.hhvec[chan][2*i] << Chan.hhvec[chan][2*i + 1] << " "; }
    std::cout << std::endl;
    std::cout << "pp" << std::endl;
    for(int i = 0; i < pp; ++i){ std::cout << Chan.ppvec[chan][2*i] << Chan.ppvec[chan][2*i + 1] << " "; }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = -E1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
    dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V1[chan], Amps2.D1.T1[chan], &hh, &pp, &pp, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (1)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = -V2(ij|kl){ij,kl}.E1(ab|kl){kl,ab} (2)
    dgemm_NN(Ints.D_ME1.V2[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (2)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = -0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (3)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (3)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = -0.5 * E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} (4)
    dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (4)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T1(ab|ij){ij,ab} = E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (5)
    dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T1 (5)" << std::endl;
    for(int i = 0; i < hh; ++i){
      for(int j = 0; j < pp; ++j){
	std::cout << Amps2.D1.T1[chan][pp*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan];
    int hp2 = Chan.nhp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    /*std::cout << "hp1" << std::endl;
    for(int i = 0; i < hp1; ++i){ std::cout << Chan.hp1vec[chan][2*i] << Chan.hp1vec[chan][2*i + 1] << " "; }
    std::cout << std::endl;
    std::cout << "hp2" << std::endl;
    for(int i = 0; i < hp2; ++i){ std::cout << Chan.hp2vec[chan][2*i] << Chan.hp2vec[chan][2*i + 1] << " "; }
    std::cout << std::endl;*/
    //T2(ab|ij){ia,jb} = E2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (6)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T2 (6)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T2[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T3(ab|ij){jb,ia} = E2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (7)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T3[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T3 (7)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T3[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T4(ab|ij){ib,ja} = -E2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (8)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T4 (8)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T4[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T5(ab|ij){ja,ib} = -E2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (9)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T5[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T5 (9)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T5[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T2(ab|ij){ia,jb} = -E2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (10)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.D1.T2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T2 (10)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T2[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.E2(db|lj){ld,jb} (11)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.S1.E2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T2 (11)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T2[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T4(ab|ij){ib,ja} = -E2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (12)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.D1.T2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T4 (12)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T4[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T4(ab|ij){ib,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.E2(ca|lj){lc,ja} (13)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.S1.E2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T4 (13)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T4[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T2(ab|ij){ia,jb} = -Q12(ad|ik){ia,kd}.T2(db|kj){kd,jb} (14)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T2 (14)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T2[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T3(ab|ij){jb,ia} = -Q12(bd|jk){jb,kd}.T2(da|ki){kd,ia} (15)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T3 (15)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T3[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T4(ab|ij){ib,ja} = Q12(bd|ik){ib,kd}.T2(da|kj){kd,ja} (16)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T4 (16)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T4[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T5(ab|ij){ja,ib} = Q12(ad|jk){ja,kd}.T2(db|ki){kd,ib} (17)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T5 (17)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T5[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T2(ab|ij){ia,jb} = -Q22(ac|il){ia,kc}.T2(cb|lj){kc,jb} (18)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T2 (18)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T2[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T3(ab|ij){jb,ia} = -Q22(bc|jl){jb,kc}.T2(ca|li){kc,ia} (19)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T3 (19)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T3[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T4(ab|ij){ib,ja} = Q22(bc|il){ib,kc}.T2(ca|lj){kc,ja} (20)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T4 (20)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T4[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T5(ab|ij){ja,ib} = Q22(ac|jl){ja,kc}.T2(cb|li){kc,ib} (21)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T5 (21)" << std::endl;
    for(int i = 0; i < hp1; ++i){
      for(int j = 0; j < hp2; ++j){
	std::cout << Amps2.D1.T5[chan][hp2*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.np[chan];
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    /*std::cout << "h" << std::endl;
    for(int i = 0; i < h; ++i){ std::cout << Chan.hvec[chan][i] << " "; }
    std::cout << std::endl;
    std::cout << "hpp" << std::endl;
    for(int i = 0; i < hpp; ++i){ std::cout << Chan.hppvec[chan][3*i] << Chan.hppvec[chan][3*i + 1] << Chan.hppvec[chan][3*i + 2] << " "; }
    std::cout << std::endl;
    std::cout << "p" << std::endl;
    for(int i = 0; i < p; ++i){ std::cout << Chan.pvec[chan][i] << " "; }
    std::cout << std::endl;*/
    if(h == 0 || hpp == 0){ continue; }
    //T6(ab|ij){jab,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.E6(cd|ik){kcd,i} (22)
    dgemm_NN(Ints.D_ME1.V5[chan], Amps1.S1.E6[chan], Amps2.D1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S2[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T6 (22)" << std::endl;
    for(int i = 0; i < hpp; ++i){
      for(int j = 0; j < h; ++j){
	std::cout << Amps2.D1.T6[chan][h*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T7(ab|ij){iab,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.E6(cd|jl){lcd,j} (23)
    dgemm_NN(Ints.D_ME1.V6[chan], Amps1.S1.E6[chan], Amps2.D1.S3[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S3[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T7 (23)" << std::endl;
    for(int i = 0; i < hpp; ++i){
      for(int j = 0; j < h; ++j){
	std::cout << Amps2.D1.T7[chan][h*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T6(ab|ij){jab,i} = T6(ab|lj){jab,l}.Q32(l|i){l,i} (26)
    dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T6 (26)" << std::endl;
    for(int i = 0; i < hpp; ++i){
      for(int j = 0; j < h; ++j){
	std::cout << Amps2.D1.T6[chan][h*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T7(ab|ij){iab,j} = -T6(ab|li){iab,l}.Q32(l|j){l,j} (27)
    dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T7 (27)" << std::endl;
    for(int i = 0; i < hpp; ++i){
      for(int j = 0; j < h; ++j){
	std::cout << Amps2.D1.T7[chan][h*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    if(p != 0){
      //T6(ab|ij){jab,i} = -V17(jc|ab){jab,c}.t2(c|i){c,i} (30)
      dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
      /*std::cout << "T6 (30)" << std::endl;
      for(int i = 0; i < hpp; ++i){
	for(int j = 0; j < h; ++j){
	  std::cout << Amps2.D1.T6[chan][h*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T7(ab|ij){iab,j} = V17(ic|ab){iab,c}.t2(c|j){c,j} (31)
      dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
      /*std::cout << "T7 (31)" << std::endl;
      for(int i = 0; i < hpp; ++i){
	for(int j = 0; j < h; ++j){
	  std::cout << Amps2.D1.T7[chan][h*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T6(ab|ij){jab,i} = -0.5 * DQ12(ab|jc){jab,c}.t2(c|i){c,i} (34)
      dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac6, &fac1, &N, &N); //fac1
      /*std::cout << "T6 (34)" << std::endl;
      for(int i = 0; i < hpp; ++i){
	for(int j = 0; j < h; ++j){
	  std::cout << Amps2.D1.T6[chan][h*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T7(ab|ij){iab,j} = 0.5 * DQ12(ab|ic){iab,c}.t2(c|j){c,j} (35)
      dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac3, &fac1, &N, &N); //fac1
      /*std::cout << "T7 (35)" << std::endl;
      for(int i = 0; i < hpp; ++i){
	for(int j = 0; j < h; ++j){
	  std::cout << Amps2.D1.T7[chan][h*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T6(ab|ij){jab,i} = Q52(jc|ab){jab,c}.t2(c|i){c,i} (38)
      dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
      /*std::cout << "T6 (38)" << std::endl;
      for(int i = 0; i < hpp; ++i){
	for(int j = 0; j < h; ++j){
	  std::cout << Amps2.D1.T6[chan][h*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T7(ab|ij){iab,j} = -Q52(ic|ab){iab,c}.t2(c|j){c,j} (39)
      dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
      /*std::cout << "T7 (39)" << std::endl;
      for(int i = 0; i < hpp; ++i){
	for(int j = 0; j < h; ++j){
	  std::cout << Amps2.D1.T7[chan][h*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    /*std::cout << "h" << std::endl;
    for(int i = 0; i < h; ++i){ std::cout << Chan.hvec[chan][i] << " "; }
    std::cout << std::endl;
    std::cout << "p" << std::endl;
    for(int i = 0; i < p; ++i){ std::cout << Chan.pvec[chan][i] << " "; }
    std::cout << std::endl;
    std::cout << "hhp" << std::endl;
    for(int i = 0; i < hhp; ++i){ std::cout << Chan.hhpvec[chan][3*i] << Chan.hhpvec[chan][3*i + 1] << Chan.hhpvec[chan][3*i + 2] << " "; }
    std::cout << std::endl;*/
    if(p == 0 || hhp == 0){ continue; }
    //T8(ab|ij){ijb,a} = T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.E8(ac|kl){klc,a} (24)
    dgemm_NN(Ints.D_ME1.V7[chan], Amps1.S1.E8[chan], Amps2.D1.S4[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S4[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T8 (24)" << std::endl;
    for(int i = 0; i < hhp; ++i){
      for(int j = 0; j < p; ++j){
	std::cout << Amps2.D1.T8[chan][p*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T9(ab|ij){ija,b} = T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.E8(bd|kl){kld,b} (25)
    dgemm_NN(Ints.D_ME1.V8[chan], Amps1.S1.E8[chan], Amps2.D1.S5[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S5[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T9 (25)" << std::endl;
    for(int i = 0; i < hhp; ++i){
      for(int j = 0; j < p; ++j){
	std::cout << Amps2.D1.T9[chan][p*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T8(ab|ij){ijb,a} = T8(db|ij){ijb,d}*Q42(d|a){d,a} (28)
    dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    /*std::cout << "T8 (28)" << std::endl;
    for(int i = 0; i < hhp; ++i){
      for(int j = 0; j < p; ++j){
	std::cout << Amps2.D1.T8[chan][p*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    //T9(ab|ij){ija,b} = -T8(da|ij){ija,d}*Q42(d|b){d,b} (29)
    dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac5, &fac1, &N, &N); //fac1
    /*std::cout << "T9 (29)" << std::endl;
    for(int i = 0; i < hhp; ++i){
      for(int j = 0; j < p; ++j){
	std::cout << Amps2.D1.T9[chan][p*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/
    if(h != 0){
      //T8(ab|ij){ijb,a} = -V18(ij|kb){ijb,k}.t3(k|a){k,a} (32)
      dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
      /*std::cout << "T8 (32)" << std::endl;
      for(int i = 0; i < hhp; ++i){
	for(int j = 0; j < p; ++j){
	  std::cout << Amps2.D1.T8[chan][p*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T9(ab|ij){ija,b} = V18(ij|ka){ija,k}.t3(k|b){k,b} (33)
      dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
      /*std::cout << "T9 (33)" << std::endl;
      for(int i = 0; i < hhp; ++i){
	for(int j = 0; j < p; ++j){
	  std::cout << Amps2.D1.T9[chan][p*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T8(ab|ij){ijb,a} = -0.5 * DQ22(kb|ij){ijb,k}.t3(k|a){k,a} (36)
      dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac6, &fac1, &N, &N); //fac1
      /*std::cout << "T8 (36)" << std::endl;
      for(int i = 0; i < hhp; ++i){
	for(int j = 0; j < p; ++j){
	  std::cout << Amps2.D1.T8[chan][p*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T9(ab|ij){ija,b} = 0.5 * DQ22(ka|ij){ija,k}.t3(k|b){k,b} (37)
      dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac3, &fac1, &N, &N); //fac1
      /*std::cout << "T9 (37)" << std::endl;
      for(int i = 0; i < hhp; ++i){
	for(int j = 0; j < p; ++j){
	  std::cout << Amps2.D1.T9[chan][p*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T8(ab|ij){ijb,a} = -Q62(ij|kb){ijb,k}.t3(k|a){k,a} (40)
      dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
      /*std::cout << "T8 (40)" << std::endl;
      for(int i = 0; i < hhp; ++i){
	for(int j = 0; j < p; ++j){
	  std::cout << Amps2.D1.T8[chan][p*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      //T9(ab|ij){ija,b} = Q62(ij|ka){ija,k}.t3(k|b){k,b} (41)
      dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
      /*std::cout << "T9 (41)" << std::endl;
      for(int i = 0; i < hhp; ++i){
	for(int j = 0; j < p; ++j){
	  std::cout << Amps2.D1.T9[chan][p*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
    }
  }
}

void Build_CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, CC_Eff &V_Eff)
{
  double fac1 = 1.0, fac5 = -1.0, fac4 = -0.5, fac3 = 0.5;
  int one = 1, hp0 = Chan.nhp1[Chan.ind0], pp0 = Chan.npp1[Chan.ind0], hh0 = Chan.nhh1[Chan.ind0];
  char N = 'N';

  if(hp0 != 0){
    // X_ia1(i|a){ia} = V9(ik|ac){ia,kc}.t1(c|k){kc}
    dgemm_NN(Ints.D_ME1.V9[Chan.ind0], Amps.S1.T1, V_Eff.X_ia1, &hp0, &one, &hp0, &fac1, &fac1, &N, &N);
    V_Eff.set_X_ia(Chan);
  }

  if(pp0 != 0){
    // X_ab1(a|b){ba} = f_ab.delta(a,b)
    for(int pp = 0; pp < pp0; ++pp){
      if(Chan.pp1vec[Chan.ind0][2*pp] == Chan.pp1vec[Chan.ind0][2*pp + 1]){
	V_Eff.X_ab1[pp] += Space.qnums[Chan.pp1vec[Chan.ind0][2*pp]].energy;
      }
    }
    if(hp0 != 0){
      // X_ab1(a|b){ba} = V16(ka|cb){ba,kc}.t1(c|k){kc}
      dgemm_NN(Ints.S_ME1.V16[Chan.ind0], Amps.S1.T1, V_Eff.X_ab1, &pp0, &one, &hp0, &fac1, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp[chan];
    if(h1 != 0 && p1 != 0){
      // X_ab3(a|b){b,a} = -X_ia3(k|b){b,k}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X_ia3[chan], Amps.S1.T3[chan], V_Eff.X_ab3[chan], &p1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(p1 != 0 && hhp1 != 0){
      // X_ab3(a|b){b,a} = -(1/2).V8(kl|bc){b,klc}.T8(ac|kl){klc,a}
      dgemm_NN(Ints.D_ME1.V8[chan], Amps.D1.T8[chan], V_Eff.X_ab3[chan], &p1, &p1, &hhp1, &fac4, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ab(Chan);

  if(hh0 != 0){
    // X1_ij1(i|j){ji} = f_ij.delta(i,j)
    for(int hh = 0; hh < hh0; ++hh){
      if(Chan.hh1vec[Chan.ind0][2*hh] == Chan.hh1vec[Chan.ind0][2*hh + 1]){
	V_Eff.X1_ij1[hh] += Space.qnums[Chan.hh1vec[Chan.ind0][2*hh]].energy;
      }
    }
    if(hp0 != 0){
      // X1_ij1(i|j){ji} = -V15(ki|jc){ji,kc}.t1(c|k){kc}
      dgemm_NN(Ints.S_ME1.V15[Chan.ind0], Amps.S1.T1, V_Eff.X1_ij1, &hh0, &one, &hp0, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], hpp1 = Chan.nhpp[chan];
    if(h1 != 0 && hpp1 != 0){
      // X1_ij2(i|j){i,j} = (1/2).V6(ik|cd){i,kcd}.T6(cd|jk){kcd,j}
      dgemm_NN(Ints.D_ME1.V6[chan], Amps.D1.T6[chan], V_Eff.X1_ij2[chan], &h1, &h1, &hpp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_ij(Chan);

  // X_ij1(i|j){ji} = X1_ij1(i|j){ji}
  if(hh0 != 0){
    for(int hh = 0; hh < hh0; ++hh){
      V_Eff.X_ij1[hh] = V_Eff.X1_ij1[hh];
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan];
    if(h1 != 0 && p1 != 0){
      // X_ij2(i|j){i,j} = X_ia2(i|d){i,d}.t2(d|j){d,j}
      dgemm_NN(V_Eff.X_ia2[chan], Amps.S1.T2[chan], V_Eff.X_ij2[chan], &h1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ij(Chan);

  if(hp0 != 0){
    // X_ai1(a|i){ia} = -V3(ka|ic){ia,kc}.t1(c|k){kc}
    dgemm_NN(Ints.D_ME1.V3[Chan.ind0], Amps.S1.T1, V_Eff.X_ai1, &hp0, &one, &hp0, &fac5, &fac1, &N, &N);
    // X_ai1(a|i){ia} = T2(ac|ik){ia,kc}.X_ia1(k|c){kc}
    dgemm_NN(Amps.D1.T2[Chan.ind0], V_Eff.X_ia1, V_Eff.X_ai1, &hp0, &one, &hp0, &fac1, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp[chan], hpp1 = Chan.nhpp[chan];
    if(h1 != 0 && p1 != 0){
      // X_ai2(a|i){a,i} = X_ab2(a|c){a,c}.t2(c|i){c,i}
      dgemm_NN(V_Eff.X_ab2[chan], Amps.S1.T2[chan], V_Eff.X_ai2[chan], &p1, &h1, &p1, &fac1, &fac1, &N, &N);
      if(hpp1 != 0){
	// X_ai2(a|i){a,i} = -(1/2).V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i}
	dgemm_NN(Ints.S_ME1.V11[chan], Amps.D1.T6[chan], V_Eff.X_ai2[chan], &p1, &h1, &hpp1, &fac4, &fac1, &N, &N);
      }
      // X_ai3(a|i){i,a} = -X1_ij3(k|i){i,k}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X1_ij3[chan], Amps.S1.T3[chan], V_Eff.X_ai3[chan], &h1, &p1, &h1, &fac5, &fac1, &N, &N);
      if(hhp1 != 0){
	// X_ai3(a|i){i,a} = -(1/2).V12(kl|ic){i,klc}.T8(ac|kl){klc,a}
	dgemm_NN(Ints.S_ME1.V12[chan], Amps.D1.T8[chan], V_Eff.X_ai3[chan], &h1, &p1, &hhp1, &fac4, &fac1, &N, &N);
      }
    }
  }
  V_Eff.set_X_ai(Chan);

  // X_ijab1(ij|ab){ab,ij} = V4(ij|ab){ab,ij}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.npp[chan] * Chan.nhh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_ijab1[chan][i] = Ints.D_ME1.V4[chan][i];
      }
    }
  }

  // X(1)_iabc(ia|bc){a,ibc} = V11(ia|bc){a,ibc}
  // X(1)_ijka(ij|ka){k,ija} = V12(ij|ka){k,ija}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.np[chan] * Chan.nhpp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_iabc1[chan][i] = Ints.S_ME1.V11[chan][i];
	V_Eff.X_iabc1[chan][i] = Ints.S_ME1.V11[chan][i];
      }
    }
    length = Chan.nh[chan] * Chan.nhhp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_ijka1[chan][i] = Ints.S_ME1.V12[chan][i];
	V_Eff.X_ijka1[chan][i] = Ints.S_ME1.V12[chan][i];
      }
    }
  }

  for(int chan = 0; chan < Chan.size3; ++chan){
    int p1 = Chan.np[chan], h1 = Chan.nh[chan], hpp1 = Chan.nhpp[chan], hhp1 = Chan.nhhp[chan];
    if(p1 != 0 && h1 != 0 && hpp1 != 0){
      // X(1)_iabc(ia|bc){a,ibc} = -((1/2)).t2(a|k){a,k}.V5(ik|bc){k,ibc}
      dgemm_NN(Amps.S1.T2[chan], Ints.D_ME1.V5[chan], V_Eff.X1_iabc1[chan], &p1, &hpp1, &h1, &fac4, &fac1, &N, &N);
      dgemm_NN(Amps.S1.T2[chan], Ints.D_ME1.V5[chan], V_Eff.X_iabc1[chan], &p1, &hpp1, &h1, &fac5, &fac1, &N, &N);
    }
    if(p1 != 0 && h1 != 0 && hhp1 != 0){
      // X(1)_ijka(ij|ka){k,ija} = -((1/2)).t3(c|k){k,c}.V7(ij|ac){c,ija}
      dgemm_NN(Amps.S1.T3[chan], Ints.D_ME1.V7[chan], V_Eff.X1_ijka1[chan], &h1, &hhp1, &p1, &fac4, &fac1, &N, &N);
      dgemm_NN(Amps.S1.T3[chan], Ints.D_ME1.V7[chan], V_Eff.X_ijka1[chan], &h1, &hhp1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_iabc(Chan);
  V_Eff.set_X1_iabc(Chan);
  V_Eff.set_X_ijka(Chan);
  V_Eff.set_X1_ijka(Chan);

  // X1_abcd(ab|cd){cd,ab} = V1(ab|cd){cd,ab}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.npp[chan] * Chan.npp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_abcd1[chan][i] = Ints.D_ME1.V1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int ppp1 = Chan.nppp[chan], p1 = Chan.np[chan], h1 = Chan.nh[chan];
    if(ppp1 != 0 && p1 != 0 && h1 != 0){
      // X1_abcd(ab|cd){acd,b} = X1_iabc(ka|cd}{acd,k}.t3(b|k){k,b}
      dgemm_NN(V_Eff.X1_iabc2[chan], Amps.S1.T3[chan], V_Eff.X1_abcd2[chan], &ppp1, &p1, &h1, &fac1, &fac1, &N, &N);
      // X1_abcd(ab|cd){bcd,a} = -X1_iabc(kb|cd}{bcd,k}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X1_iabc2[chan], Amps.S1.T3[chan], V_Eff.X1_abcd3[chan], &ppp1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_abcd(Chan);

  // X_abcd(ab|cd){cd,ab} = X1_abcd(ab|cd){cd,ab}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.npp[chan] * Chan.npp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_abcd1[chan][i] = V_Eff.X1_abcd1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int pp1 = Chan.npp[chan], hh1 = Chan.nhh[chan];
    if(pp1 != 0 && hh1 != 0){
      // X_abcd(ab|cd){cd,ab} = (1/2).V4(kl|cd}{cd,kl}.T1(ab|kl){kl,ab}
      dgemm_NN(Ints.D_ME1.V4[chan], Amps.D1.T1[chan], V_Eff.X_abcd1[chan], &pp1, &pp1, &hh1, &fac3, &fac1, &N, &N);
    }
  }

  // X_ijkl(ij|kl){ij,kl} = V2(ij|kl){ij,kl}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.nhh[chan] * Chan.nhh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_ijkl1[chan][i] = Ints.D_ME1.V2[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int hhh1 = Chan.nhhh[chan], h1 = Chan.nh[chan], p1 = Chan.np[chan];
    if(hhh1 != 0 && h1 != 0 && p1 != 0){
      // X_ijkl(ij|kl){ijk,l} = X1_ijkl(ij|kc}{ijk,c}.t3(c|l){c,l}
      dgemm_NN(V_Eff.X1_ijka2[chan], Amps.S1.T2[chan], V_Eff.X_ijkl2[chan], &hhh1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X_ijkl(ij|kl){ijl,k} = -X1_ijkl(ij|lc}{ijl,c}.t3(c|k){c,k}
      dgemm_NN(V_Eff.X1_ijka2[chan], Amps.S1.T2[chan], V_Eff.X_ijkl3[chan], &hhh1, &h1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh1 = Chan.nhh[chan], pp1 = Chan.npp[chan];
    if(hh1 != 0 && pp1 != 0){
      // X_ijkl(ij|kl){kl,ij} = (1/2).T1(kl|cd}{kl,cd}.V4(ij|cd){cd,ij}
      dgemm_NN(Amps.D1.T1[chan], Ints.D_ME1.V4[chan], V_Eff.X_ijkl4[chan], &hh1, &hh1, &pp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ijkl(Chan);

  // X1_iajb(ia|jb){ib,ja} = V(ia|jb){ib,ja}
  // X3_iajb(ia|jb){ib,ja} = V(ia|jb){ib,ja}
  for(int chan = 0; chan < Chan.size2; ++chan){
    int length = Chan.nhp2[chan] * Chan.nhp2[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_iajb1[chan][i] = Ints.D_ME1.V3[chan][i];
	V_Eff.X3_iajb1[chan][i] = Ints.D_ME1.V3[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp1[chan], hpp1 = Chan.nhpp1[chan];
    if(h1 * p1 * hhp1 != 0){
      // X1_iajb(ia|jb){ijb,a} = (1/2).V14(ki|jb}{ijb,k}.t3(a|k){k,a}
      dgemm_NN(Ints.S_ME1.V14[chan], Amps.S1.T3[chan], V_Eff.X1_iajb2[chan], &hhp1, &p1, &h1, &fac3, &fac1, &N, &N);
      // X3_iajb(ia|jb){ijb,a} = V14(ki|jb}{ijb,k}.t3(a|k){k,a}
      dgemm_NN(Ints.S_ME1.V14[chan], Amps.S1.T3[chan], V_Eff.X3_iajb2[chan], &hhp1, &p1, &h1, &fac1, &fac1, &N, &N);
    }
    if(h1 * p1 * hpp1 != 0){
      // X1_iajb(ia|jb){iab,j} = X1_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      dgemm_NN(V_Eff.X1_iabc3[chan], Amps.S1.T2[chan], V_Eff.X1_iajb3[chan], &hpp1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X3_iajb(ia|jb){iab,j} = (1/2).X_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      dgemm_NN(V_Eff.X_iabc3[chan], Amps.S1.T2[chan], V_Eff.X3_iajb3[chan], &hpp1, &h1, &p1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_iajb(Chan);
  V_Eff.set_X3_iajb(Chan);

  // X_iajb(ia|jb){ib,ja} = X3_iajb(ia|jb){ib,ja}
  for(int chan = 0; chan < Chan.size2; ++chan){
    int length = Chan.nhp2[chan] * Chan.nhp2[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_iajb1[chan][i] = V_Eff.X3_iajb1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan], hp2 = Chan.nhp2[chan];
    if(hp1 != 0 && hp2 != 0){
      // X_iajb(ia|jb){ib,ja} = -V10(ik|cb}{ib,kc}.T5(ca|jk){kc,ja}
      dgemm_NN(Ints.D_ME1.V10[chan], Amps.D1.T5[chan], V_Eff.X_iajb1[chan], &hp2, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hpp1 = Chan.nhpp1[chan];
    if(h1 != 0 && p1 != 0 && hpp1 != 0){
      // X_iajb(ia|jb){iab,j} = (1/2).X_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      dgemm_NN(V_Eff.X_iabc3[chan], Amps.S1.T2[chan], V_Eff.X_iajb3[chan], &hpp1, &h1, &p1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_iajb(Chan);

  // X_abic(ab|ic){iab,c} = V17(ab|ic){iab,c}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.nhpp[chan] * Chan.np[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_abic1[chan][i] = Ints.S_ME1.V17[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], ppp1 = Chan.nppp[chan], hpp1 = Chan.nhpp[chan], hpp2 = Chan.nhpp1[chan];
    if(h1 != 0 && p1 != 0 && hpp1 != 0){
      // X_abic(ab|ic){iab,c} = -T7(ab|ik}{iab,k}.X_ia(k|c){k,c}  
      dgemm_NN(Amps.D1.T7[chan], V_Eff.X_ia2[chan], V_Eff.X_abic1[chan], &hpp1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && ppp1 != 0){
      // X_abic(ab|ic){abc,i} = V*(ab|dc){abc,d}.t2(d|i){d,i}
      dgemm_NN(V_Eff.V_abcd[chan], Amps.S1.T2[chan], V_Eff.X_abic2[chan], &ppp1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && hpp2 != 0){
      // X_abic(ab|ic){icb',a'} = -X1_iajb(kb|ic){icb',k'}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X1_iajb4[chan], Amps.S1.T3[chan], V_Eff.X_abic3[chan], &hpp2, &p1, &h1, &fac5, &fac1, &N, &N);
      // X_abic(ab|ic){ica',b'} = X1_iajb(ka|ic){ica',k'}.t3(b|k){k,b}
      dgemm_NN(V_Eff.X1_iajb4[chan], Amps.S1.T3[chan], V_Eff.X_abic4[chan], &hpp2, &p1, &h1, &fac1, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int pp1 = Chan.npp1[chan], hp1 = Chan.nhp1[chan], hp2 = Chan.nhp2[chan];
    if(pp1 != 0 && hp1 != 0 && hp2 != 0){
      // X_abic(ab|ic){bc',ia'} = X_iabc(kb|dc){bc',kd}.T3(ad|ik){kd,ia'}
      dgemm_NN(V_Eff.X_iabc4[chan], Amps.D1.T3[chan], V_Eff.X_abic5[chan], &pp1, &hp2, &hp1, &fac1, &fac1, &N, &N);
      // X_abic(ab|ic){ac',ib'} = -X_iabc(ka|dc){ac',kd}.T3(bd|ik){kd,ib'}
      dgemm_NN(V_Eff.X_iabc4[chan], Amps.D1.T3[chan], V_Eff.X_abic6[chan], &pp1, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int pp1 = Chan.npp[chan], hh1 = Chan.nhh[chan], hp1 = Chan.nhp[chan];
    if(pp1 != 0 && hh1 != 0 && hp1 != 0){
      // X_abic(ab|ic){ic,ab} = (1/2).X_ijka(kl|ic){ic,kl}.T1(ab|kl){kl,ab}
      dgemm_NN(V_Eff.X_ijka5[chan], Amps.D1.T1[chan], V_Eff.X_abic7[chan], &hp1, &pp1, &hh1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_abic(Chan);

  // X2_iajk(ia|jk){jka,i} = V18(ia|jk){jka,i}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.nhhp[chan] * Chan.nh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X2_iajk1[chan][i] = Ints.S_ME1.V18[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhh1 = Chan.nhhh[chan], hhp1 = Chan.nhhp1[chan];
    if(h1 != 0 && p1 != 0 && hhh1 != 0){
      // X2_iajk(ia|jk){ijk,a} = -V*(il|jk){ijk,l}.t3(a|l){l,a}
      dgemm_NN(V_Eff.V_ijkl[chan], Amps.S1.T3[chan], V_Eff.X2_iajk2[chan], &hhh1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && hhp1 != 0){
      // X2_iajk(ia|jk){jia',k'} = X3_iajb(ia|jd){jia',d'}.t2(d|k){d,k}
      dgemm_NN(V_Eff.X3_iajb5[chan], Amps.S1.T2[chan], V_Eff.X2_iajk3[chan], &hhp1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X2_iajk(ia|jk){kia',j'} = -X3_iajb(ia|kd){kia',d'}.t2(d|j){d,j}
      dgemm_NN(V_Eff.X3_iajb5[chan], Amps.S1.T2[chan], V_Eff.X2_iajk4[chan], &hhp1, &h1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hh1 = Chan.nhh1[chan], hp1 = Chan.nhp1[chan], hp2 = Chan.nhp2[chan];
    if(hh1 != 0 && hp1 != 0 && hp2 != 0){
      // X2_iajk(ia|jk){ij',ka'} = X_ijka(il|jc){ij',lc}.T3(ac|kl){lc,ka'}
      dgemm_NN(V_Eff.X_ijka4[chan], Amps.D1.T3[chan], V_Eff.X2_iajk5[chan], &hh1, &hp2, &hp1, &fac1, &fac1, &N, &N);
      // X2_iajk(ia|jk){ik',ja'} = -X_ijka(il|kc){ik',lc}.T3(ac|jl){lc,ja'}
      dgemm_NN(V_Eff.X_ijka4[chan], Amps.D1.T3[chan], V_Eff.X2_iajk6[chan], &hh1, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hp1 = Chan.nhp[chan], hh1 = Chan.nhh[chan], pp1 = Chan.npp[chan];
    if(hh1 != 0 && pp1 != 0 && hp1 != 0){
      // X2_iajk(ia|jk){jk,ia} = (1/2).T1(cd|jk){jk,cd}.X_iabc(ia|cd){cd,ia}
      dgemm_NN(Amps.D1.T1[chan], V_Eff.X_iabc5[chan], V_Eff.X2_iajk7[chan], &hh1, &hp1, &pp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X2_iajk(Chan);

  // X_iajk(ia|jk){jka,i} = X2_iajk(ia|jk){jka,i}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.nhhp[chan] * Chan.nh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_iajk1[chan][i] = V_Eff.X2_iajk1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp[chan];
    if(h1 != 0 && p1 != 0 && hhp1 != 0){
      // X_iajk(ia|jk){jka,i} = T8(ca|jk){jka,c}.X_ia(i|c){c,i}
      dgemm_NN(Amps.D1.T8[chan], V_Eff.X_ia3[chan], V_Eff.X_iajk1[chan], &hhp1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
  }

  /*std::cout << "X_ia" << std::endl;
  for(int i = 0; i < Chan.nhp1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.hp1vec[Chan.ind0][2*i] << "|" << Chan.hp1vec[Chan.ind0][2*i + 1] << "> = " << V_Eff.X_ia1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X_ab" << std::endl;
  for(int i = 0; i < Chan.npp1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.pp1vec[Chan.ind0][2*i + 1] << "|" << Chan.pp1vec[Chan.ind0][2*i] << "> = " << V_Eff.X_ab1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X1_ij" << std::endl;
  for(int i = 0; i < Chan.nhh1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.hh1vec[Chan.ind0][2*i + 1] << "|" << Chan.hh1vec[Chan.ind0][2*i] << "> = " << V_Eff.X1_ij1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X_ij" << std::endl;
  for(int i = 0; i < Chan.nhh1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.hh1vec[Chan.ind0][2*i + 1] << "|" << Chan.hh1vec[Chan.ind0][2*i] << "> = " << V_Eff.X_ij1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X_ai" << std::endl;
  for(int i = 0; i < Chan.nhp2[Chan.ind0]; ++i){
    std::cout << "<" << Chan.hp1vec[Chan.ind0][2*i + 1] << "|" << Chan.hp1vec[Chan.ind0][2*i] << "> = " << V_Eff.X_ai1[i] << std::endl;
  }
  std::cout << std::endl;*/
  /*std::cout << "X_ijab" << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.npp[i]; ++j){
      int pind1 = Chan.ppvec[i][2*j];
      int pind2 = Chan.ppvec[i][2*j + 1];
      for(int k = 0; k < Chan.nhh[i]; ++k){
	int hind1 = Chan.hhvec[i][2*k];
	int hind2 = Chan.hhvec[i][2*k + 1];
	std::cout << "<" << hind1 << " " << hind2 << "|" << pind1 << " " << pind2 << "> = " << V_Eff.X_ijab1[i][j*Chan.nhh[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X1_iabc" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.np[i]; ++j){
      int pind1 = Chan.pvec[i][j];
      for(int k = 0; k < Chan.nhpp[i]; ++k){
	int hind1 = Chan.hppvec[i][3*k];
	int pind2 = Chan.hppvec[i][3*k + 1];
	int pind3 = Chan.hppvec[i][3*k + 2];
	std::cout << "<" << hind1 << " " << pind1 << "|" << pind2 << " " << pind3 << "> = " << V_Eff.X1_iabc1[i][j*Chan.nhpp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_iabc" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.np[i]; ++j){
      int pind1 = Chan.pvec[i][j];
      for(int k = 0; k < Chan.nhpp[i]; ++k){
	int hind1 = Chan.hppvec[i][3*k];
	int pind2 = Chan.hppvec[i][3*k + 1];
	int pind3 = Chan.hppvec[i][3*k + 2];
	std::cout << "<" << hind1 << " " << pind1 << "|" << pind2 << " " << pind3 << "> = " << V_Eff.X_iabc1[i][j*Chan.nhpp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X1_ijka" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nh[i]; ++j){
      int hind3 = Chan.hvec[i][j];
      for(int k = 0; k < Chan.nhhp[i]; ++k){
	int hind1 = Chan.hhpvec[i][3*k];
	int hind2 = Chan.hhpvec[i][3*k + 1];
	int pind1 = Chan.hhpvec[i][3*k + 2];
	std::cout << "<" << hind1 << " " << hind2 << "|" << hind3 << " " << pind1 << "> = " << V_Eff.X1_ijka1[i][j*Chan.nhhp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_ijka" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nh[i]; ++j){
      int hind3 = Chan.hvec[i][j];
      for(int k = 0; k < Chan.nhhp[i]; ++k){
	int hind1 = Chan.hhpvec[i][3*k];
	int hind2 = Chan.hhpvec[i][3*k + 1];
	int pind1 = Chan.hhpvec[i][3*k + 2];
	std::cout << "<" << hind1 << " " << hind2 << "|" << hind3 << " " << pind1 << "> = " << V_Eff.X_ijka1[i][j*Chan.nhhp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X1_abcd" << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.npp[i]; ++j){
      int pind3 = Chan.ppvec[i][2*j];
      int pind4 = Chan.ppvec[i][2*j + 1];
      for(int k = 0; k < Chan.npp[i]; ++k){
	int pind1 = Chan.ppvec[i][2*k];
	int pind2 = Chan.ppvec[i][2*k + 1];
	std::cout << "<" << pind1 << " " << pind2 << "|" << pind3 << " " << pind4 << "> = " << V_Eff.X1_abcd1[i][j*Chan.npp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_abcd" << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.npp[i]; ++j){
      int pind3 = Chan.ppvec[i][2*j];
      int pind4 = Chan.ppvec[i][2*j + 1];
      for(int k = 0; k < Chan.npp[i]; ++k){
	int pind1 = Chan.ppvec[i][2*k];
	int pind2 = Chan.ppvec[i][2*k + 1];
	std::cout << "<" << pind1 << " " << pind2 << "|" << pind3 << " " << pind4 << "> = " << V_Eff.X_abcd1[i][j*Chan.npp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X_ijkl" << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.nhh[i]; ++j){
      int hind1 = Chan.hhvec[i][2*j];
      int hind2 = Chan.hhvec[i][2*j + 1];
      for(int k = 0; k < Chan.nhh[i]; ++k){
	int hind3 = Chan.hhvec[i][2*k];
	int hind4 = Chan.hhvec[i][2*k + 1];
	std::cout << "<" << hind1 << " " << hind2 << "|" << hind3 << " " << hind4 << "> = " << V_Eff.X_ijkl1[i][j*Chan.nhh[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X1_iajb" << std::endl;
  for(int i = 0; i < Chan.size2; ++i){
    for(int j = 0; j < Chan.nhp2[i]; ++j){
      int hind1 = Chan.hp2vec[i][2*j];
      int pind2 = Chan.hp2vec[i][2*j + 1];
      for(int k = 0; k < Chan.nhp2[i]; ++k){
	int hind2 = Chan.hp2vec[i][2*k];
	int pind1 = Chan.hp2vec[i][2*k + 1];
	std::cout << "<" << hind1 << " " << pind1 << "|" << hind2 << " " << pind2 << "> = " << V_Eff.X1_iajb1[i][j*Chan.nhp2[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X3_iajb" << std::endl;
  for(int i = 0; i < Chan.size2; ++i){
    for(int j = 0; j < Chan.nhp2[i]; ++j){
      int hind1 = Chan.hp2vec[i][2*j];
      int pind2 = Chan.hp2vec[i][2*j + 1];
      for(int k = 0; k < Chan.nhp2[i]; ++k){
	int hind2 = Chan.hp2vec[i][2*k];
	int pind1 = Chan.hp2vec[i][2*k + 1];
	std::cout << "<" << hind1 << " " << pind1 << "|" << hind2 << " " << pind2 << "> = " << V_Eff.X3_iajb1[i][j*Chan.nhp2[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_iajb" << std::endl;
  for(int i = 0; i < Chan.size2; ++i){
    for(int j = 0; j < Chan.nhp2[i]; ++j){
      int hind1 = Chan.hp2vec[i][2*j];
      int pind2 = Chan.hp2vec[i][2*j + 1];
      for(int k = 0; k < Chan.nhp2[i]; ++k){
	int hind2 = Chan.hp2vec[i][2*k];
	int pind1 = Chan.hp2vec[i][2*k + 1];
	std::cout << "<" << hind1 << " " << pind1 << "|" << hind2 << " " << pind2 << "> = " << V_Eff.X_iajb1[i][j*Chan.nhp2[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X_abic" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nhpp[i]; ++j){
      int hind1 = Chan.hppvec[i][3*j];
      int pind1 = Chan.hppvec[i][3*j + 1];
      int pind2 = Chan.hppvec[i][3*j + 2];
      for(int k = 0; k < Chan.np[i]; ++k){
	int pind3 = Chan.pvec[i][k];
	std::cout << "<" << pind1 << " " << pind2 << "|" << hind1 << " " << pind3 << "> = " << V_Eff.X_abic1[i][j*Chan.np[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X2_iajk" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nhhp[i]; ++j){
      int hind2 = Chan.hhpvec[i][3*j];
      int hind3 = Chan.hhpvec[i][3*j + 1];
      int pind1 = Chan.hhpvec[i][3*j + 2];
      for(int k = 0; k < Chan.nh[i]; ++k){
	int hind1 = Chan.hvec[i][k];
	std::cout << "<" << hind1 << " " << pind1 << "|" << hind2 << " " << hind3 << "> = " << V_Eff.X2_iajk1[i][j*Chan.nh[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
  /*std::cout << "X_iajk" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nhhp[i]; ++j){
      int hind2 = Chan.hhpvec[i][3*j];
      int hind3 = Chan.hhpvec[i][3*j + 1];
      int pind1 = Chan.hhpvec[i][3*j + 2];
      for(int k = 0; k < Chan.nh[i]; ++k){
	int hind1 = Chan.hvec[i][k];
	std::cout << "<" << hind1 << " " << pind1 << "|" << hind2 << " " << hind3 << "> = " << V_Eff.X_iajk1[i][j*Chan.nh[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;*/
}



void EE_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff)
{
  std::cout << "Computing Excited States ..." << std::endl << std::endl;

  double *Ham;
  unsigned long long *bra;
  unsigned long long *ket;
  int *p0 = new int[0];
  int *q0 = new int[0];
  int *p1 = new int[1];
  int *q1 = new int[1];
  int *p2 = new int[2];
  int *q2 = new int[2];
  int *p3 = new int[3];
  int *q3 = new int[3];
  int *q4 = new int[4];
  double ME, ME1, ME0;
  double tempen;
  int total;

  int Nbit = std::ceil(Space.indtot/64.0);
  int **chanvec1;
  int **chanvec2;
  chanvec1 = new int*[Chan.size2];
  chanvec2 = new int*[Chan.size2];

  State *Spectrum;

  State state;
  int *count1;
  int *count2;
  int ind;
  count1 = new int[Chan.size2];
  count2 = new int[Chan.size2];
  for(int i = 0; i < Chan.size2; ++i){
    count1[i] = 0;
    count2[i] = 0;
  }

  for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // p
    for(int chan2 = 0; chan2 < Chan.size3; ++chan2){ // h
      if(Chan.np[chan1] == 0 || Chan.nh[chan2] == 0){ continue; }
      minus(state, Chan.qnums3[chan1], Chan.qnums3[chan2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, state);
      ++count1[ind];
    }
  }
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){ // pp
    for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // hh
      if(Chan.npp[chan1] == 0 || Chan.nhh[chan2] == 0){ continue; }
      minus(state, Chan.qnums1[chan1], Chan.qnums1[chan2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, state);
      if(ind < 0 || ind >= Chan.size2){ continue; }
      ++count2[ind];
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    if(count1[i] != 0){ chanvec1[i] = new int[2 * count1[i]]; }
    if(count2[i] != 0){ chanvec2[i] = new int[2 * count2[i]]; }
    count1[i] = 0;
    count2[i] = 0;
  }

  for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // p
    for(int chan2 = 0; chan2 < Chan.size3; ++chan2){ // h
      if(Chan.np[chan1] == 0 || Chan.nh[chan2] == 0){ continue; }
      minus(state, Chan.qnums3[chan1], Chan.qnums3[chan2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, state);
      chanvec1[ind][2 * count1[ind]] = chan1;
      chanvec1[ind][2 * count1[ind] + 1] = chan2;
      ++count1[ind];
    }
  }
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){ // pp
    for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // hh
      if(Chan.npp[chan1] == 0 || Chan.nhh[chan2] == 0){ continue; }
      minus(state, Chan.qnums1[chan1], Chan.qnums1[chan2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, state);
      if(ind < 0 || ind >= Chan.size2){ continue; }
      chanvec2[ind][2 * count2[ind]] = chan1;
      chanvec2[ind][2 * count2[ind] + 1] = chan2;
      ++count2[ind];
    }
  }

  for(int chan = 0; chan < Chan.size2; ++chan){
    if(chan == Chan.ind0){ ++total; }
    for(int i = 0; i < count1[chan]; ++i){
      for(int j = 0; j < Chan.nh[chanvec1[chan][2*i + 1]]; ++j){
	for(int k = 0; k < Chan.np[chanvec1[chan][2*i]]; ++k){
	  ++total;
	}
      }
    }
    for(int i = 0; i < count2[chan]; ++i){
      for(int j = 0; j < Chan.nhh[chanvec2[chan][2*i + 1]]; ++j){
	if(Chan.hhvec[chanvec2[chan][2*i + 1]][2*j] >= Chan.hhvec[chanvec2[chan][2*i + 1]][2*j+1]){ continue; }
	for(int k = 0; k < Chan.npp[chanvec2[chan][2*i]]; ++k){
	  if(Chan.ppvec[chanvec2[chan][2*i]][2*k] >= Chan.ppvec[chanvec2[chan][2*i]][2*k+1]){ continue; }
	  ++total;
	}
      }
    }
  }
  //std::cout << "Total = " << total << std::endl;

  Spectrum = new State[total];
  total = 0;

  int *hpvec2;
  int *hhppvec2;
  for(int chan = 0; chan < Chan.size2; ++chan){
    int count01 = 0;
    int count02 = 0;

    for(int i = 0; i < count1[chan]; ++i){
      for(int j = 0; j < Chan.nh[chanvec1[chan][2*i + 1]]; ++j){
	for(int k = 0; k < Chan.np[chanvec1[chan][2*i]]; ++k){
	  ++count01;
	}
      }
    }
    for(int i = 0; i < count2[chan]; ++i){
      for(int j = 0; j < Chan.nhh[chanvec2[chan][2*i + 1]]; ++j){
	if(Chan.hhvec[chanvec2[chan][2*i + 1]][2*j] >= Chan.hhvec[chanvec2[chan][2*i + 1]][2*j+1]){ continue; }
	for(int k = 0; k < Chan.npp[chanvec2[chan][2*i]]; ++k){
	  if(Chan.ppvec[chanvec2[chan][2*i]][2*k] >= Chan.ppvec[chanvec2[chan][2*i]][2*k+1]){ continue; }
	  ++count02;
	}
      }
    }

    hpvec2 = new int[2 * count01];
    hhppvec2 = new int[4 * count02];
    count01 = 0;
    count02 = 0;
    for(int i = 0; i < count1[chan]; ++i){
      for(int j = 0; j < Chan.nh[chanvec1[chan][2*i + 1]]; ++j){
	for(int k = 0; k < Chan.np[chanvec1[chan][2*i]]; ++k){
	  hpvec2[2*count01] = Chan.hvec[chanvec1[chan][2*i + 1]][j];
	  hpvec2[2*count01 + 1] = Chan.pvec[chanvec1[chan][2*i]][k];
	  ++count01;
	}
      }
    }
    for(int i = 0; i < count2[chan]; ++i){
      for(int j = 0; j < Chan.nhh[chanvec2[chan][2*i + 1]]; ++j){
	if(Chan.hhvec[chanvec2[chan][2*i + 1]][2*j] >= Chan.hhvec[chanvec2[chan][2*i + 1]][2*j+1]){ continue; }
	for(int k = 0; k < Chan.npp[chanvec2[chan][2*i]]; ++k){
	  if(Chan.ppvec[chanvec2[chan][2*i]][2*k] >= Chan.ppvec[chanvec2[chan][2*i]][2*k+1]){ continue; }
	  hhppvec2[4*count02] = Chan.hhvec[chanvec2[chan][2*i + 1]][2*j];
	  hhppvec2[4*count02 + 1] = Chan.hhvec[chanvec2[chan][2*i + 1]][2*j+1];
	  hhppvec2[4*count02 + 2] = Chan.ppvec[chanvec2[chan][2*i]][2*k];
	  hhppvec2[4*count02 + 3] = Chan.ppvec[chanvec2[chan][2*i]][2*k+1];
	  ++count02;
	}
      }
    }

    int N = count01 + count02;
    if(chan == Chan.ind0){ N += 1; }
    if(N == 0){ continue; }
    Ham = new double[N*N];
    for(int i = 0; i < N*N; ++i){ Ham[i] = 0.0; }

    for(int col = 0; col < N; ++col){
      ket = new unsigned long long[Nbit];
      for(int i = 0; i < Nbit; ++i){ ket[i] = 0; }
      if(col < count01){ bitsetup(hpvec2, ket, col, 2); }
      else if(col < count01 + count02){ bitsetup(hhppvec2, ket, (col - count01), 4); }
      for(int row = 0; row < N; ++row){
	bra = new unsigned long long[Nbit];
	for(int i = 0; i < Nbit; ++i){ bra[i] = 0; }
	if(row < count01){ bitsetup(hpvec2, bra, row, 2); }
	else if(row < count01 + count02){ bitsetup(hhppvec2, bra, (row - count01), 4); }
	ME = 0.0;

	for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){ // {i+a} -> ia
	  ME0 = V_Eff.X_ia1[hp1];
	  if(ME0 == 0.0){ continue; }
	  q2[0] = Chan.hp1vec[Chan.ind0][2*hp1];
	  q2[1] = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	  ME1 = matrixe(bra, ket, Nbit, p0, 0, q2, 2, ME0);
	  ME += ME1;
	}

	for(int pp1 = 0; pp1 < Chan.npp1[Chan.ind0]; ++pp1){ // {a+b}
	  ME0 = V_Eff.X_ab1[pp1];
	  if(ME0 == 0.0){ continue; }
	  p1[0] = Chan.pp1vec[Chan.ind0][2*pp1 + 1];
	  q1[0] = Chan.pp1vec[Chan.ind0][2*pp1];
	  ME1 = matrixe(bra, ket, Nbit, p1, 1, q1, 1, ME0);
	  ME += ME1;
	}

	for(int hh1 = 0; hh1 < Chan.nhh1[Chan.ind0]; ++hh1){ // {i+j} -> ij+ = -1 * j+i
	  ME0 = V_Eff.X_ij1[hh1];
	  if(ME0 == 0.0){ continue; }
	  p1[0] = Chan.hh1vec[Chan.ind0][2*hh1];
	  q1[0] = Chan.hh1vec[Chan.ind0][2*hh1 + 1];
	  ME1 = matrixe(bra, ket, Nbit, p1, 1, q1, 1, ME0);
	  ME -= ME1;
	}

	for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){ // {a+i} -> a+i+
	  ME0 = V_Eff.X_ai1[hp1];
	  if(ME0 == 0.0){ continue; }
	  p2[0] = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	  p2[1] = Chan.hp1vec[Chan.ind0][2*hp1];
	  ME1 = matrixe(bra, ket, Nbit, p2, 2, q0, 0, ME0);
	  ME += ME1;
	}

	for(int i = 0; i < Chan.size3; ++i){ // {i+a+cb} -> ia+cb = a+ibc
	  for(int p = 0; p < Chan.np[i]; ++p){
	    for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
	      ME0 = V_Eff.X_iabc1[i][p*Chan.nhpp[i] + hpp1];
	      if(ME0 == 0.0){ continue; }
	      p1[0] = Chan.pvec[i][p];
	      q3[0] = Chan.hppvec[i][3*hpp1];
	      q3[1] = Chan.hppvec[i][3*hpp1 + 1];
	      q3[2] = Chan.hppvec[i][3*hpp1 + 2];
	      ME1 = matrixe(bra, ket, Nbit, p1, 1, q3, 3, ME0);
	      ME += 0.5 * ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size3; ++i){ // {i+j+ak} -> ijak+ = -1 * k+ija
	  for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
	    for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
	      ME0 = V_Eff.X_ijka1[i][h1*Chan.nhhp[i] + hhp1];
	      if(ME0 == 0.0){ continue; }
	      p1[0] = Chan.hvec[i][h1];
	      q3[0] = Chan.hhpvec[i][3*hhp1];
	      q3[1] = Chan.hhpvec[i][3*hhp1 + 1];
	      q3[2] = Chan.hhpvec[i][3*hhp1 + 2];
	      ME1 = matrixe(bra, ket, Nbit, p1, 1, q3, 3, ME0);
	      ME -= 0.5 * ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size3; ++i){ // {a+b+ci} -> a+b+ci+ = -1 * i+a+b+c
	  for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
	    for(int p = 0; p < Chan.np[i]; ++p){
	      ME0 = V_Eff.X_abic1[i][hpp1*Chan.np[i] + p];
	      if(ME0 == 0.0){ continue; }
	      p3[0] = Chan.hppvec[i][3*hpp1];
	      p3[1] = Chan.hppvec[i][3*hpp1 + 1];
	      p3[2] = Chan.hppvec[i][3*hpp1 + 2];
	      q1[0] = Chan.pvec[i][p];
	      ME1 = matrixe(bra, ket, Nbit, p3, 3, q1, 1, ME0);
	      ME -= 0.5 * ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size3; ++i){ // {i+a+kj} -> ia+k+j+ = j+k+a+i
	  for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
	    for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
	      ME0 = V_Eff.X_iajk1[i][hhp1*Chan.nh[i] + h1];
	      if(ME0 == 0.0){ continue; }
	      p3[0] = Chan.hhpvec[i][3*hhp1];
	      p3[1] = Chan.hhpvec[i][3*hhp1 + 1];
	      p3[2] = Chan.hhpvec[i][3*hhp1 + 2];
	      q1[0] = Chan.hvec[i][h1];
	      ME1 = matrixe(bra, ket, Nbit, p3, 3, q1, 1, ME0);
	      ME += 0.5 * ME1;
	    }
	  }
	}
	
	for(int i = 0; i < Chan.size1; ++i){ // {a+b+dc}
	  for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	    for(int pp2 = 0; pp2 < Chan.npp[i]; ++pp2){
	      ME0 = V_Eff.X_abcd1[i][pp1*Chan.npp[i] + pp2];
	      if(ME0 == 0.0){ continue; }
	      p2[0] = Chan.ppvec[i][2*pp2];
	      p2[1] = Chan.ppvec[i][2*pp2 + 1];
	      q2[0] = Chan.ppvec[i][2*pp1 + 1];
	      q2[1] = Chan.ppvec[i][2*pp1];
	      ME1 = matrixe(bra, ket, Nbit, p2, 2, q2, 2, ME0);
	      ME += 0.25 * ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size1; ++i){ // {i+j+lk} -> ijl+k+ = l+k+ij
	  for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
	    for(int hh2 = 0; hh2 < Chan.nhh[i]; ++hh2){
	      ME0 = V_Eff.X_ijkl1[i][hh1*Chan.nhh[i] + hh2];
	      if(ME0 == 0.0){ continue; }
	      p2[0] = Chan.hhvec[i][2*hh2 + 1];
	      p2[1] = Chan.hhvec[i][2*hh2];
	      q2[0] = Chan.hhvec[i][2*hh1];
	      q2[1] = Chan.hhvec[i][2*hh1 + 1];
	      ME1 = matrixe(bra, ket, Nbit, p2, 2, q2, 2, ME0);
	      ME += 0.25 * ME1;
	    }
	  }
	}
	
	for(int i = 0; i < Chan.size2; ++i){ // {i+a+bj} -> ia+bj+ = j+a+ib
	  for(int hp1 = 0; hp1 < Chan.nhp2[i]; ++hp1){
	    for(int hp2 = 0; hp2 < Chan.nhp2[i]; ++hp2){
	      ME0 = V_Eff.X_iajb1[i][hp1*Chan.nhp2[i] + hp2];
	      if(ME0 == 0.0){ continue; }
	      p2[0] = Chan.hp2vec[i][2*hp2];
	      p2[1] = Chan.hp2vec[i][2*hp2 + 1];
	      q2[0] = Chan.hp2vec[i][2*hp1];
	      q2[1] = Chan.hp2vec[i][2*hp1 + 1];
	      ME1 = matrixe(bra, ket, Nbit, p2, 2, q2, 2, ME0);
	      ME += ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size1; ++i){ // {i+j+ba} -> ijba
	  for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
	      ME0 = V_Eff.X_ijab1[i][pp1*Chan.nhh[i] + hh1];
	      if(ME0 == 0.0){ continue; }
	      q4[0] = Chan.hhvec[i][2*hh1];
	      q4[1] = Chan.hhvec[i][2*hh1 + 1];
	      q4[2] = Chan.ppvec[i][2*pp1 + 1];
	      q4[3] = Chan.ppvec[i][2*pp1];
	      ME1 = matrixe(bra, ket, Nbit, p0, 0, q4, 4, ME0);
	      ME += 0.25 * ME1;
	    }
	  }
	}

	// fill column major
	Ham[N*col + row] = ME;
	delete[] bra;
      }
      delete[] ket;
    }
  
    char job = 'V';
    int lwork = 5*N, info;
    double *Vl = new double[N * N];
    double *Vr = new double[N * N];
    double *wl = new double[N];
    double *wr = new double[N];
    double *work = new double[5*N];
    dgeev_(&job, &job, &N, Ham, &N, wr, wl, Vl, &N, Vr, &N, work, &lwork, &info);

    for(int j = 0; j < N; ++j){
      Spectrum[total] = Chan.qnums2[chan];
      Spectrum[total].energy = std::sqrt(wr[j]*wr[j] + wl[j]*wl[j]);
      ++total;
    }
   
    delete[] Vl;
    delete[] Vr;
    delete[] wl;
    delete[] wr;
    delete[] work;
    delete[] hpvec2;
    delete[] hhppvec2;
    delete[] Ham;
  }
  delete[] p0;
  delete[] q0;
  delete[] p1;
  delete[] q1;
  delete[] p2;
  delete[] q2;
  delete[] p3;
  delete[] q3;
  delete[] q4;

  delete[] count1;
  delete[] count2;

  for(int chan = 0; chan < Chan.size2; ++chan){
    delete[] chanvec1[chan];
    delete[] chanvec2[chan];
  }
  delete[] chanvec1;
  delete[] chanvec2;

  for(int i = 0; i < total - 1; ++i){
    ind = i;
    tempen = Spectrum[i].energy;
    for(int j = i + 1; j < total; ++j){
      if(tempen - Spectrum[j].energy >= 1e-8){
	tempen = Spectrum[j].energy;
	ind = j;
      }
    }
    state = Spectrum[i];
    Spectrum[i] = Spectrum[ind];
    Spectrum[ind] = state;
  }

  for(int i = 0; i < total; ++i){
    std::cout << Spectrum[i].par << " " << Spectrum[i].ml << " " << Spectrum[i].m << " : " << Spectrum[i].energy << std::endl;
  }

  delete[] Spectrum;
  
}


void PA_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff, State *states, double *nums)
{
  std::cout << "Computing Particle-Attached States ..." << std::endl << std::endl;

  double *Ham;
  double tempen;
  State state;
  int ind;
  int *pvec2;
  int *hppvec2;
  int length;

  for(int chan = 0; chan < Chan.size3; ++chan){
    if(Chan.qnums3[chan].m != 1){ continue; }
    if(Chan.qnums3[chan].ml != 0 && Chan.qnums3[chan].ml != 1 && Chan.qnums3[chan].ml != 2){ continue; }

    int count01 = 0;
    int count02 = 0;

    count01 = Chan.np[chan];
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // h
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // pp
	if(Chan.nh[chan1] * Chan.npp[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int pp = 0; pp < Chan.npp[chan2]; ++pp){
	      if(Chan.ppvec[chan2][2*pp] < Chan.ppvec[chan2][2*pp + 1]){
		count02 += Chan.nh[chan1];
	      }
	    }
	  }
	}
      }
    }

    pvec2 = new int[count01];
    hppvec2 = new int[3 * count02];
    count01 = 0;
    count02 = 0;
    for(int i = 0; i < Chan.np[chan]; ++i){
      pvec2[count01] = Chan.pvec[chan][i];
      ++count01;
    }
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // h
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // pp
	if(Chan.nh[chan1] * Chan.npp[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int pp = 0; pp < Chan.npp[chan2]; ++pp){
	      if(Chan.ppvec[chan2][2*pp] < Chan.ppvec[chan2][2*pp + 1]){
		for(int h = 0; h < Chan.nh[chan1]; ++h){
		  hppvec2[3*count02] = Chan.hvec[chan1][h];
		  hppvec2[3*count02 + 1] = Chan.ppvec[chan2][2*pp];
		  hppvec2[3*count02 + 2] = Chan.ppvec[chan2][2*pp + 1];
		  ++count02;
		}
	      }
	    }
	  }
	}
      }
    }

    int N = count01 + count02;
    if(N == 0){ continue; }
    Ham = new double[N*N];

    #pragma omp parallel
    {
      double ME;
      int b1, b2, b3, k1, k2, k3, ind0, ind;
      int bra, ket, key1, key2;
      State tb;
      length = count01 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + bra] = 0.0;
	b1 = pvec2[bra];
	k1 = pvec2[ket];
	ME = 0.0;
	ind = Chan.pp1_map[Chan.ind0][Hash2(k1, b1, Space.indtot)];
	ME += V_Eff.X_ab1[ind];
	Ham[N*ket + bra] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + (bra + count01)] = 0.0;
	b1 = hppvec2[3*bra];
	b2 = hppvec2[3*bra + 1];
	b3 = hppvec2[3*bra + 2];
	k1 = pvec2[ket];
	ME = 0.0;
	if(b2 == k1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(b1, b3, Space.indtot)];
	  ME += V_Eff.X_ai1[ind];
	}
	if(b3 == k1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(b1, b2, Space.indtot)];
	  ME -= V_Eff.X_ai1[ind];
	}
	ind0 = Chan.indvec[k1];
	key1 = Chan.hpp_map[ind0][Hash3(b1, b2, b3, Space.indtot)];
	key2 = Chan.p_map[ind0][k1];
	ind = key1 * Chan.np[ind0] + key2;
	ME -= V_Eff.X_abic1[ind0][ind];
	Ham[N*ket + (bra + count01)] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	bra = int(braket%count01);
	ket = int((braket - bra)/count01);
	Ham[N*(ket + count01) + bra] = 0.0;
	k1 = hppvec2[3*ket];
	k2 = hppvec2[3*ket + 1];
	k3 = hppvec2[3*ket + 2];
	b1 = pvec2[bra];
	ME = 0.0;
	if(k2 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k1, k3, Space.indtot)];
	  ME += V_Eff.X_ia1[ind];
	}
	if(k3 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k1, k2, Space.indtot)];
	  ME -= V_Eff.X_ia1[ind];
	}
	ind0 = Chan.indvec[b1];
	key1 = Chan.p_map[ind0][b1];
	key2 = Chan.hpp_map[ind0][Hash3(k1, k2, k3, Space.indtot)];
	ind = key1 * Chan.nhpp[ind0] + key2;
	ME -= V_Eff.X_iabc1[ind0][ind];
	Ham[N*(ket + count01) + bra] += ME;
      }
      length = count02 * count02;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count02);
	bra = int((braket - ket)/count02);
	Ham[N*(ket + count01) + (bra + count01)] = 0.0;
	b1 = hppvec2[3*bra];
	b2 = hppvec2[3*bra + 1];
	b3 = hppvec2[3*bra + 2];
	k1 = hppvec2[3*ket];
	k2 = hppvec2[3*ket + 1];
	k3 = hppvec2[3*ket + 2];
	ME = 0.0;
	if(b1 == k1){
	  if(b3 == k3){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k2, b2, Space.indtot)];
	    ME += V_Eff.X_ab1[ind];
	  }
	  if(b2 == k2){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k3, b3, Space.indtot)];
	    ME += V_Eff.X_ab1[ind];
	  }
	  if(b3 == k2){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k3, b2, Space.indtot)];
	    ME -= V_Eff.X_ab1[ind];
	  }
	  if(b2 == k3){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k2, b3, Space.indtot)];
	    ME -= V_Eff.X_ab1[ind];
	  }
	}
	if(b2 == k2 && b3 == k3){
	  ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k1, Space.indtot)];
	  ME -= V_Eff.X_ij1[ind];
	}
	if(b1 == k1){
	  plus(tb, Space.qnums[b2], Space.qnums[b3]);
	  ind0 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  plus(tb, Space.qnums[k2], Space.qnums[k3]);
	  if(ind0 == ChanInd_2b_dir(Parameters.basis, Space, tb)){
	    key1 = Chan.pp_map[ind0][Hash2(k2, k3, Space.indtot)];
	    key2 = Chan.pp_map[ind0][Hash2(b2, b3, Space.indtot)];
	    ind = key1 * Chan.npp[ind0] + key2;
	    ME += V_Eff.X_abcd1[ind0][ind];
	  }
	}
	if(b3 == k3){
	  minus(tb, Space.qnums[k2], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b2], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k2, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b2, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b2 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b3 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b2], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b2, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b2 == k3){
	  minus(tb, Space.qnums[k2], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k2, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	Ham[N*(ket + count01) + (bra + count01)] += ME;
      }
    }

    /*for(int i = 0; i < N; ++i){
      for(int j = 0; j < N; ++j){
	std::cout << Ham[N*j + i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/

    double norm1p;
    int info = 0;
    if(N <= 3){
      char job1 = 'N';
      char job2 = 'V';
      int lwork = 10*N;
      double *Vl = new double[N * N];
      double *Vr = new double[N * N];
      double *wl = new double[N];
      double *wr = new double[N];
      double *work = new double[lwork];
      info = 0;
      dgeev_(&job1, &job2, &N, Ham, &N, wr, wl, Vl, &N, Vr, &N, work, &lwork, &info);

      ind = 0;
      tempen = wr[0];
      for(int i = 1; i < N; ++i){
	if(tempen - wr[i] >= 1e-8){
	  tempen = wr[i];
	  ind = i;
	}
      }

      //for(int i = 0; i < N; ++i){ std::cout << Vr[N*ind + i] << " "; }
      //std::cout << std::endl;

      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += Vr[N*ind + i] * Vr[N*ind + i]; }
      norm1p = std::sqrt(norm1p);

      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = wr[ind];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;
      
      //std::cout << "! " << Chan.qnums3[chan].ml << " " << Chan.qnums3[chan].m << " : " << wr[ind] << " : " << norm1p << std::endl;
      
      delete[] Vl;
      delete[] Vr;
      delete[] wl;
      delete[] wr;
      delete[] work;
    }
    else{
      int ido = 0; // status integer is zero at start
      char bmat[] = "I"; // standard eigenvalue problem
      char which[] = "SM"; // smallest magnitude eigenvalues
      int nev = 1; // number of eigenvalues to calculate
      double tol = 1.0e-10; // error tolerance
      double *resid = new double[N];
      int ncv = 3*nev + 2; //5*nev; // number of lanczos vectors
      if(ncv > N){ ncv = N; }
      double *v = new double[N*ncv];
      for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
      int ldv = N;
      double *workd = new double[3*N];
      int lworkl = 4*ncv*(ncv + 2);
      double *workl = new double[lworkl];
      info = 0;
      int iparam[11];
      int ipntr[14];
      
      int ishift = 1;
      int mxiter = 5000;
      int mode = 1;
      iparam[0] = ishift;
      iparam[2] = mxiter;
      iparam[6] = mode;
      
      do{
	dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
	
	if(ido != -1 && ido != 1){ break; }

        #pragma omp parallel for
	for(int j = 0; j < N; ++j){
	  workd[ipntr[1]-1 + j] = 0.0;
	  for(int k = 0; k < N; ++k){
	    workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k];
	  }
	}
      } while(true);
      
      bool rvec = true;
      char howmny = 'A';
      int *select = new int[ncv];
      double *dr = new double[nev + 1];
      double *di = new double[nev + 1];
      double *z = new double[N * (nev + 1)];
      for(int i = 0; i < nev + 1; ++i){
	dr[i] = 0.0;
	di[i] = 0.0;
	for(int j = 0; j < N; ++j){
	  z[N*i + j] = 0.0;
	}
      }
      int ldz = N;
      double sigmar;
      double sigmai;
      double *workev = new double[3 * ncv];

      dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

      //for(int i = 0; i < N; ++i){ std::cout << z[i] << " "; }
      //std::cout << std::endl;

      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += z[i] * z[i]; }
      norm1p = std::sqrt(norm1p);

      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = dr[0];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;

      //std::cout << "! " << Chan.qnums3[chan].ml << " " << Chan.qnums3[chan].m << " : " << dr[0] << " : " << norm1p << std::endl;
   
      delete[] resid;
      delete[] v;
      delete[] workd;
      delete[] workl;
      delete[] select;
      delete[] dr;
      delete[] di;
      delete[] z;
      delete[] workev;
    }
    delete[] pvec2;
    delete[] hppvec2;
    delete[] Ham;
  }
}

void PR_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff, State *states, double *nums)
{
  std::cout << "Computing Particle-Removed States ..." << std::endl << std::endl;

  double *Ham;
  double tempen;
  State state;
  int ind;
  int *hvec2;
  int *hhpvec2;
  int length;

  for(int chan = 0; chan < Chan.size3; ++chan){
    if(Chan.qnums3[chan].m != 1){ continue; }
    if(Chan.qnums3[chan].ml != 0 && Chan.qnums3[chan].ml != 1 && Chan.qnums3[chan].ml != 2){ continue; }

    int count01 = 0;
    int count02 = 0;

    count01 = Chan.nh[chan];
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // p
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // hh
	if(Chan.np[chan1] * Chan.nhh[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int hh = 0; hh < Chan.nhh[chan2]; ++hh){
	      if(Chan.hhvec[chan2][2*hh] < Chan.hhvec[chan2][2*hh + 1]){
		count02 += Chan.np[chan1];
	      }
	    }
	  }
	}
      }
    }

    hvec2 = new int[count01];
    hhpvec2 = new int[3 * count02];
    count01 = 0;
    count02 = 0;
    for(int i = 0; i < Chan.nh[chan]; ++i){
      hvec2[count01] = Chan.hvec[chan][i];
      ++count01;
    }
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // p
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // hh
	if(Chan.np[chan1] * Chan.nhh[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int hh = 0; hh < Chan.nhh[chan2]; ++hh){
	      if(Chan.hhvec[chan2][2*hh] < Chan.hhvec[chan2][2*hh + 1]){
		for(int p = 0; p < Chan.np[chan1]; ++p){
		  hhpvec2[3*count02] = Chan.hhvec[chan2][2*hh];
		  hhpvec2[3*count02 + 1] = Chan.hhvec[chan2][2*hh + 1];
		  hhpvec2[3*count02 + 2] = Chan.pvec[chan1][p];
		  ++count02;
		}
	      }
	    }
	  }
	}
      }
    }

    int N = count01 + count02;
    if(N == 0){ continue; }
    Ham = new double[N*N];

    #pragma omp parallel
    {
      double ME;
      int b1, b2, b3, k1, k2, k3, ind0, ind;
      int bra, ket, key1, key2;
      State tb;
      length = count01 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + bra] = 0.0;
	b1 = hvec2[bra];
	k1 = hvec2[ket];
	ME = 0.0;
	ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k1, Space.indtot)];
	ME -= V_Eff.X_ij1[ind];
	Ham[N*ket + bra] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + (bra + count01)] = 0.0;
	b1 = hhpvec2[3*bra];
	b2 = hhpvec2[3*bra + 1];
	b3 = hhpvec2[3*bra + 2];
	k1 = hvec2[ket];
	ME = 0.0;
	if(b1 == k1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(b2, b3, Space.indtot)];
	  ME -= V_Eff.X_ai1[ind];
	}
	if(b2 == k1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(b1, b3, Space.indtot)];
	  ME += V_Eff.X_ai1[ind];
	}
	ind0 = Chan.indvec[k1];
	key1 = Chan.hhp_map[ind0][Hash3(b1, b2, b3, Space.indtot)];
	key2 = Chan.h_map[ind0][k1];
	ind = key1 * Chan.nh[ind0] + key2;
	ME += V_Eff.X_iajk1[ind0][ind];
	Ham[N*ket + (bra + count01)] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	bra = int(braket%count01);
	ket = int((braket - bra)/count01);
	Ham[N*(ket + count01) + bra] = 0.0;
	k1 = hhpvec2[3*ket];
	k2 = hhpvec2[3*ket + 1];
	k3 = hhpvec2[3*ket + 2];
	b1 = hvec2[bra];
	ME = 0.0;
	if(k1 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k2, k3, Space.indtot)];
	  ME -= V_Eff.X_ia1[ind];
	}
	if(k2 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k1, k3, Space.indtot)];
	  ME += V_Eff.X_ia1[ind];
	}
	ind0 = Chan.indvec[b1];
	key1 = Chan.h_map[ind0][b1];
	key2 = Chan.hhp_map[ind0][Hash3(k1, k2, k3, Space.indtot)];
	ind = key1 * Chan.nhhp[ind0] + key2;
	ME += V_Eff.X_ijka1[ind0][ind];
	Ham[N*(ket + count01) + bra] += ME;
      }
      length = count02 * count02;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count02);
	bra = int((braket - ket)/count02);
	Ham[N*(ket + count01) + (bra + count01)] = 0.0;
	b1 = hhpvec2[3*bra];
	b2 = hhpvec2[3*bra + 1];
	b3 = hhpvec2[3*bra + 2];
	k1 = hhpvec2[3*ket];
	k2 = hhpvec2[3*ket + 1];
	k3 = hhpvec2[3*ket + 2];
	ME = 0.0;
	if(b3 == k3){
	  if(b2 == k2){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k1, Space.indtot)];
	    ME -= V_Eff.X_ij1[ind];
	  }
	  if(b1 == k1){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b2, k2, Space.indtot)];
	    ME -= V_Eff.X_ij1[ind];
	  }
	  if(b2 == k1){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k2, Space.indtot)];
	    ME += V_Eff.X_ij1[ind];
	  }
	  if(b1 == k2){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b2, k1, Space.indtot)];
	    ME += V_Eff.X_ij1[ind];
	  }
	}
	if(b1 == k1 && b2 == k2){
	  ind = Chan.pp1_map[Chan.ind0][Hash2(k3, b3, Space.indtot)];
	  ME += V_Eff.X_ab1[ind];
	}
	if(b3 == k3){
	  plus(tb, Space.qnums[b1], Space.qnums[b2]);
	  ind0 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  plus(tb, Space.qnums[k1], Space.qnums[k2]);
	  if(ind0 == ChanInd_2b_dir(Parameters.basis, Space, tb)){
	    key1 = Chan.hh_map[ind0][Hash2(k1, k2, Space.indtot)];
	    key2 = Chan.hh_map[ind0][Hash2(b1, b2, Space.indtot)];
	    ind = key1 * Chan.nhh[ind0] + key2;
	    ME += V_Eff.X_ijkl1[ind0][ind];
	  }
	}
	if(b2 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b1 == k1){
	  minus(tb, Space.qnums[k3], Space.qnums[k2]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b2]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k2, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b2, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b2 == k1){
	  minus(tb, Space.qnums[k3], Space.qnums[k2]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k2, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b1 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b2]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b2, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	Ham[N*(ket + count01) + (bra + count01)] += ME;
      }
    }

    /*for(int i = 0; i < N; ++i){
      for(int j = 0; j < N; ++j){
	std::cout << Ham[N*j + i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/

    double norm1p;
    int info = 0;
    if(N <= 3){
      char job1 = 'N';
      char job2 = 'V';
      int lwork = 10*N;
      double *Vl = new double[N * N];
      double *Vr = new double[N * N];
      double *wl = new double[N];
      double *wr = new double[N];
      double *work = new double[lwork];
      info = 0;
      dgeev_(&job1, &job2, &N, Ham, &N, wr, wl, Vl, &N, Vr, &N, work, &lwork, &info);

      ind = 0;
      tempen = wr[0];
      for(int i = 1; i < N; ++i){
	if(tempen - wr[i] <= 1e-8){
	  tempen = wr[i];
	  ind = i;
	}
      }

      //for(int i = 0; i < N; ++i){ std::cout << Vr[N*ind + i] << " "; }
      //std::cout << std::endl;

      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += Vr[N*ind + i] * Vr[N*ind + i]; }
      norm1p = std::sqrt(norm1p);
      
      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = wr[ind];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;

      //std::cout << "! " << Chan.qnums3[chan].ml << " " << Chan.qnums3[chan].m << " : " << wr[ind] << " : " << norm1p << std::endl;
      
      delete[] Vl;
      delete[] Vr;
      delete[] wl;
      delete[] wr;
      delete[] work;
    }
    else{
      int ido = 0; // status integer is zero at start
      char bmat[] = "I"; // standard eigenvalue problem
      char which[] = "LM"; // smallest magnitude eigenvalues
      int nev = 1; // number of eigenvalues to calculate
      double tol = 1.0e-10; // error tolerance
      double *resid = new double[N];
      int ncv = 3*nev + 2; //5*nev; // number of lanczos vectors
      if(ncv > N){ ncv = N; }
      double *v = new double[N*ncv];
      for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
      int ldv = N;
      double *workd = new double[3*N];
      int lworkl = 4*ncv*(ncv + 2);
      double *workl = new double[lworkl];
      info = 0;
      int iparam[11];
      int ipntr[14];
      
      int ishift = 1;
      int mxiter = 5000;
      int mode = 1;
      iparam[0] = ishift;
      iparam[2] = mxiter;
      iparam[6] = mode;
      
      do{
	dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
	
	if(ido != -1 && ido != 1){ break; }

        #pragma omp parallel for
	for(int j = 0; j < N; ++j){
	  workd[ipntr[1]-1 + j] = 0.0;
	  for(int k = 0; k < N; ++k){
	    workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k];
	  }
	}
      } while(true);
      
      bool rvec = true;
      char howmny = 'A';
      int *select = new int[ncv];
      double *dr = new double[nev + 1];
      double *di = new double[nev + 1];
      double *z = new double[N * (nev + 1)];
      for(int i = 0; i < nev + 1; ++i){
	dr[i] = 0.0;
	di[i] = 0.0;
	for(int j = 0; j < N; ++j){
	  z[N*i + j] = 0.0;
	}
      }
      int ldz = N;
      double sigmar;
      double sigmai;
      double *workev = new double[3 * ncv];

      dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

      //for(int i = 0; i < N; ++i){ std::cout << z[i] << " "; }
      //std::cout << std::endl;

      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += z[i] * z[i]; }
      norm1p = std::sqrt(norm1p);

      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = dr[0];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;

      //std::cout << "! " << Chan.qnums3[chan].ml << " " << Chan.qnums3[chan].m << " : " << dr[0] << " : " << norm1p << std::endl;
   
      delete[] resid;
      delete[] v;
      delete[] workd;
      delete[] workl;
      delete[] select;
      delete[] dr;
      delete[] di;
      delete[] z;
      delete[] workev;
    }
    delete[] hvec2;
    delete[] hhpvec2;
    delete[] Ham;
  }
}

void bitsetup(int *vec, unsigned long long *state, const int &begin, const int &size)
{
  unsigned long long one = 1;
  for(int i = size*begin; i < size*(begin + 1); ++i){ state[int(std::floor(vec[i]/64.0))] += (one << (vec[i]%64)); }
}

double matrixe(unsigned long long *bra, unsigned long long *ket, const int &size, int *p, const int &psize, int *q, const int &qsize, const double &ME)
{
  unsigned long long one = 1, zero = 0;
  unsigned long long comp;
  unsigned long long *tempket = new unsigned long long[size];
  for(int i = 0; i < size; ++i){ tempket[i] = ket[i]; }

  unsigned long long *pbit = new unsigned long long[psize];
  unsigned long long *qbit = new unsigned long long[qsize];
  int *p1 = new int[psize];
  int *q1 = new int[qsize];
  double phase = 1.0;
  int bcount;
  int flag1 = 1;
  std::bitset<64> bket;

  for(int i = 0; i < psize; ++i){
    pbit[i] = one << (p[i]%64);
    p1[i] = int(std::floor(p[i]/64));
    if((pbit[i] & bra[p1[i]]) == 0){ flag1 = 0; goto stop1; }
  }
  for(int i = 0; i < qsize; ++i){
    qbit[i] = one << (q[i]%64);
    q1[i] = int(std::floor(q[i]/64));
    if((qbit[i] & tempket[q1[i]]) == 0){ flag1 = 0; goto stop1; }
  }

 stop1:;
  if(flag1 == 1){
    for(int i = qsize-1; i >= 0; --i){
      comp = tempket[q1[i]] & ~(~zero << (q[i]%64));
      for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
      for(int j = q1[i]-1; j >= 0; --j){ bket = tempket[j]; bcount += bket.count(); }
      phase = phase*pow(-1.0, bcount); tempket[q1[i]] = tempket[q1[i]]^qbit[i];
    }
    for(int i = psize-1; i >= 0; --i){
      comp = tempket[p1[i]] & ~(~zero << (p[i]%64));
      for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
      for(int j = p1[i]-1; j >= 0; --j){ bket = tempket[j]; bcount += bket.count(); }
      phase = phase*pow(-1.0, bcount); tempket[p1[i]] = tempket[p1[i]]^pbit[i];
    }
    for(int i = 0; i < size; ++i){
      if(bra[i] != tempket[i]){ phase = 0.0; break; }
    }
  }
  else{ phase = 0.0; }

  delete[] pbit;
  delete[] qbit;
  delete[] tempket;
  delete[] p1;
  delete[] q1;

  return phase * ME;
}


/*std::vector<std::vector<std::vector<double> > > p(1);
  std::vector<std::vector<std::vector<double> > > delp(1);
  std::vector<double> B;
  std::vector<double> B2;
  std::vector<int> ipiv;
  std::vector<double> work;
  int lwork;
  int info = 0;
  int maxl = 10;
  int N = 1;

  p[0] = CCin.T1;
  delp[0] = CCin.T1;
  CCout.CCDE = 0.0;
  error = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	tempen = Space.levelsen[ind1] + Space.levelsen[ind2] - Space.levelsen[ind3] - Space.levelsen[ind4];
	tempt = CCME.HHPP1[chan][pind * Chan.hh[chan] + hind] / tempen;
	CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	error += 1.0e10;
	CCin.Evec[chan][hind * Chan.pp[chan] + pind] = tempen;
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	p[0][chan][hind * Chan.pp[chan] + pind] = tempt;
	delp[0][chan][hind * Chan.pp[chan] + pind] = 1.0e10;
      }
    }
  }
  CCout = CCin;
  error = std::sqrt(error)/(Parameters.P + Parameters.N);

  while((error > 10e-10 && ind < 1000) || ind < 10){
    if(N != maxl){
      p.push_back(CCin.T1);
      delp.push_back(CCin.T1);
      ++N;
      B.assign((N + 1) * (N + 1), 0.0);
      for(int i = N - 1; i >= 0; --i){
	for(int j = N - 1; j >= 0; --j){
	  B[(N + 1) * i + j] = B[N * i + j];
	}
      }
    }
    else{
      p[0] = CCin.T1;
      delp[0] = CCin.T1;
      p.resize(1);
      delp.resize(1);
      N = 1;
    }
    Doubles_Step(Space, Chan, CCME, CCin, CCout);
    CCout.CCDE = 0.0;
    error = 0.0;
    for(int chan = 0; chan < Chan.size1; ++chan){
      for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	  tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	  tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	  tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	  CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	  error += (tempt - p[N - 2][chan][hind * Chan.pp[chan] + pind]) * (tempt - p[N - 2][chan][hind * Chan.pp[chan] + pind]);
	  p[N - 1][chan][hind * Chan.pp[chan] + pind] = tempt;
	  delp[N - 1][chan][hind * Chan.pp[chan] + pind] = tempt - p[N - 2][chan][hind * Chan.pp[chan] + pind];
	  for(int i = 0; i < N; ++i){
	    B[(N + 1) * i + N - 1] += delp[i][chan][hind * Chan.pp[chan] + pind] * delp[N - 1][chan][hind * Chan.pp[chan] + pind];
	    if(i == N - 1){ break; } // don't double count
	    B[(N + 1) * (N - 1) + i] += delp[N - 1][chan][hind * Chan.pp[chan] + pind] * delp[i][chan][hind * Chan.pp[chan] + pind];
	  }
	}
      }
    }
    for(int i = 0; i < N; ++i){ B[(N + 1) * i + N] = -1.0; B[(N + 1) * N + i] = -1.0; }
    B[(N + 1) * N + N] = 0.0;

    int P = N + 1;
    lwork = sizeof(double) * P;
    ipiv.resize(P);
    work.resize(sizeof(double) * P);
    B2 = B;
    dgetrf_(&P, &P, & *B2.begin(), &P, & *ipiv.begin(), &info);
    dgetri_(&P, & *B2.begin(), &P, & *ipiv.begin(), & *work.begin(), &lwork, &info);
    error = std::sqrt(error)/(Parameters.P + Parameters.N);
    
    for(int chan = 0; chan < Chan.size1; ++chan){
      for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	  double tempt = 0.0;
	  for(int i = 0; i < N; ++i){ tempt += B2[(N + 1) * i + N] * p[i][chan][hind * Chan.pp[chan] + pind]; }
	  CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	}
      }
    }
    std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCout.CCDE << ", error = " << error << std::endl;
    ++ind;
  }
  std::cout << std::endl << std::endl;
  
  return CCout;*/

  ////////////////////////////////////////

  /*std::vector<std::vector<std::vector<double> > > Vin(3); // 0, V/E
  std::vector<std::vector<std::vector<double> > > F(2); // Vout - Vin
  std::vector<std::vector<std::vector<double> > > delV(1);
  std::vector<std::vector<std::vector<double> > > delF(1);
  std::vector<std::vector<std::vector<double> > > u(1);
  std::vector<double> a(1, 0.0);
  std::vector<double> B(1, 0.0);
  std::vector<double> w(2, 1.0);
  std::vector<int> ipiv;
  std::vector<double> work;
  double tempt, tempen, norm;
  int lwork;
  int info = 0;
  int maxl = 10;
  int P = 1;
  int N = 2;

  double alpha = 0.9;
  double w0 = 0.01;

  Vin[0] = CCin.T1;
  Vin[1] = CCin.T1;
  F[0] = CCin.T1;
  F[1] = CCin.T1;
  delV[0] = CCin.T1;
  delF[0] = CCin.T1;
  u[0] = CCin.T1;
  norm = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    int ind1, ind2, ind3, ind4;
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	tempen = Space.levelsen[ind1] + Space.levelsen[ind2] - Space.levelsen[ind3] - Space.levelsen[ind4];
	tempt = alpha * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind] / tempen;
	CCin.Evec[chan][hind * Chan.pp[chan] + pind] = tempen;
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	Vin[1][chan][hind * Chan.pp[chan] + pind] = tempt;
	F[0][chan][hind * Chan.pp[chan] + pind] = tempt;
	norm += tempt*tempt;
      }
    }
  }
  w[0] = 0.0;
  norm = std::sqrt(norm);
  CCout.Evec = CCin.Evec;

  Doubles_Step(Space, Chan, CCME, CCin, CCout);

  norm = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	F[1][chan][hind * Chan.pp[chan] + pind] = tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind];
	norm += (F[1][chan][hind * Chan.pp[chan] + pind] - F[0][chan][hind * Chan.pp[chan] + pind]) *
	  (F[1][chan][hind * Chan.pp[chan] + pind] - F[0][chan][hind * Chan.pp[chan] + pind]);
      }
    }
  }
  w[1] = 1.0;
  norm = std::sqrt(norm);

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	delV[0][chan][hind * Chan.pp[chan] + pind] = (Vin[1][chan][hind * Chan.pp[chan] + pind] - Vin[0][chan][hind * Chan.pp[chan] + pind])/norm;
	delF[0][chan][hind * Chan.pp[chan] + pind] = (F[1][chan][hind * Chan.pp[chan] + pind] - F[0][chan][hind * Chan.pp[chan] + pind])/norm;
	u[0][chan][hind * Chan.pp[chan] + pind] = alpha * delF[0][chan][hind * Chan.pp[chan] + pind] + delV[0][chan][hind * Chan.pp[chan] + pind];
	a[0] += w[0] * delF[0][chan][hind * Chan.pp[chan] + pind] * w[0] * delF[0][chan][hind * Chan.pp[chan] + pind];
      }
    }
  }
  a[0] += w0*w0;

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	tempt = Vin[1][chan][hind * Chan.pp[chan] + pind] + alpha * F[1][chan][hind * Chan.pp[chan] + pind];
	tempt -= w[0] * delF[0][chan][hind * Chan.pp[chan] + pind] * F[1][chan][hind * Chan.pp[chan] + pind] * u[0][chan][hind * Chan.pp[chan] + pind] / a[0];
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
      }
    }
  }
  Vin[2] = CCin.T1;

  while((error > 10e-8 && ind < 1000) || ind < 10)
    {
      Doubles_Step(Space, Chan, CCME, CCin, CCout);
      P = int(delF.size());
      N = int(delF.size() + 1);
      if(P != maxl){
	F.resize(F.size() + 1);
	Vin.resize(Vin.size() + 1);
	delV.resize(delV.size() + 1);
	delF.resize(delF.size() + 1);
	u.resize(u.size() + 1);
	w.resize(w.size() + 1);
	Vin[Vin.size() - 1] = CCout.T1;
	F[F.size() - 1] = CCout.T1;
	delV[delV.size() - 1] = CCout.T1;
	delF[delF.size() - 1] = CCout.T1;
	u[u.size() - 1] = CCout.T1;
	w[w.size() - 1] = 1;
	P = int(delF.size());
	N = int(delF.size() + 1);
	a.resize(P * P);
	for(int i = P - 2; i >= 0; --i){
	  for(int j = P - 2; j >= 0; --j){
	    a[P * i + j] = a[(P - 1) * i + j];
	  }
	}
      }
      else{
	for(int i = 0; i < (P - 1); ++i){ delV[i] = delV[i+1]; delF[i] = delF[i+1]; u[i] = u[i+1]; }
	for(int i = 0; i < (N - 1); ++i){ F[i] = F[i+1]; w[i] = w[i+1]; }
	for(int i = 0; i < N; ++i){ Vin[i] = Vin[i+1]; }
	for(int i = 0; i < (P - 1); ++i){ 
	  for(int j = 0; j < (P - 1); ++j){
	    a[P * i + j] = a[P * (i + 1) + (j + 1)];
	  }
	}
      }

      for(int i = 0; i < P; ++i){
	a[P * i + (P - 1)] = 0.0;
	a[P * (P - 1) + i] = 0.0;
      }

      CCout.CCDE = 0.0;
      error = 0.0;
      norm = 0.0;
      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	    tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	    tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	    CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	    F[N - 1][chan][hind * Chan.pp[chan] + pind] = tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind];
	    error += (tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind]) * (tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind]);
	    norm += (F[N - 1][chan][hind * Chan.pp[chan] + pind] - F[N - 2][chan][hind * Chan.pp[chan] + pind]) *
	      (F[N - 1][chan][hind * Chan.pp[chan] + pind] - F[N - 2][chan][hind * Chan.pp[chan] + pind]);
	  }
	}
      }
      w[N - 1] = 1.0/error;
      norm = std::sqrt(norm);
      error = std::sqrt(error);

      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    delV[P - 1][chan][hind * Chan.pp[chan] + pind] = (Vin[N - 1][chan][hind * Chan.pp[chan] + pind] - Vin[N - 2][chan][hind * Chan.pp[chan] + pind])/norm;
	    delF[P - 1][chan][hind * Chan.pp[chan] + pind] = (F[N - 1][chan][hind * Chan.pp[chan] + pind] - F[N - 2][chan][hind * Chan.pp[chan] + pind])/norm;
	    u[P - 1][chan][hind * Chan.pp[chan] + pind] = alpha * delF[P - 1][chan][hind * Chan.pp[chan] + pind] + delV[P - 1][chan][hind * Chan.pp[chan] + pind];
	    for(int i = 0; i < P; ++i){
	      a[P * i + (P - 1)] += w[i] * delF[i][chan][hind * Chan.pp[chan] + pind] * w[P - 1] * delF[P - 1][chan][hind * Chan.pp[chan] + pind];
	      if(i == P - 1){ break; } // don't double count
	      a[P * (P - 1) + i] += w[P - 1] * delF[P - 1][chan][hind * Chan.pp[chan] + pind] * w[i] * delF[i][chan][hind * Chan.pp[chan] + pind];
	    }
	  }
	}
      }
      a[P * (P - 1) + (P - 1)] += w0 * w0;

      lwork = P;
      ipiv.resize(P);
      work.resize(sizeof(double) * P);
      B = a;
      dgetrf_(&P, &P, & *B.begin(), &P, & *ipiv.begin(), &info);
      dgetri_(&P, & *B.begin(), &P, & *ipiv.begin(), & *work.begin(), &lwork, &info);

      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    tempt = Vin[N - 1][chan][hind * Chan.pp[chan] + pind] + alpha * F[N - 1][chan][hind * Chan.pp[chan] + pind];
	    for(int i = 0; i < P; ++i){
	      for(int j = 0; j < P; ++j){
		tempt += w[i] * delF[i][chan][hind * Chan.pp[chan] + pind] * F[N - 1][chan][hind * Chan.pp[chan] + pind] * 
		  B[P * i + j] * u[j][chan][hind * Chan.pp[chan] + pind];
	      }
	    }
	    CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	  }
	}
      }
      Vin[N] = CCin.T1;
      
      std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCout.CCDE << ", error = " << error << std::endl;
      ++ind;
    }
  std::cout << std::endl << std::endl;
    
return CCout;*/
