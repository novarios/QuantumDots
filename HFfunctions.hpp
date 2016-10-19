#ifndef HFFUNCTIONS_H
#define HFFUNCTIONS_H

#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>
#include <bitset>
#include <iomanip>
#include <omp.h>
#include <cstdarg>
#include <cstdlib>
#include <stdlib.h>
#include <stdint.h>
#include <unordered_map>

extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

struct Single_Particle_States;
struct HF_Channels;
struct HF_Matrix_Elements;

//void Separate_Particles_Holes(Single_Particle_States &States, const HF_Channels &Chan);
//void Build_Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &States);
//void Setup_Channels_HF(const Input_Parameters &Parameters, const Model_Space &Space, HF_Channels &Chan);
void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Read_Matrix_Elements_M(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Hartree_Fock_States(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &States, const HF_Matrix_Elements &ME);
void Convert_To_HF_Matrix_Elements(const HF_Channels &Chan, const Single_Particle_States &States, HF_Matrix_Elements &ME);
void Setup_HF_Space(Model_Space &Space, const Single_Particle_States &States, const HF_Channels &Chan);
void Get_Matrix_Elements(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_J(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Read_Matrix_Elements_HO(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME);
void Read_HO_ME_From_File(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME);

struct Single_Particle_States{
  int hp;
  int hn;
  double*** holes; //list of sp hole states given as vectors of coefficients
  double*** particles; //list of sp particle states given as vectors of coefficients
  double** h_energies; //list of sp-hole energies
  double** pt_energies; //list of sp-particle energies
  double*** vectors; 
  double** energies;
  int* h; //number of holes in OBchan
  int* p; //number of particles in OBchan

  Single_Particle_States(){};
  Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan);
  void delete_struct(const HF_Channels &Chan);
  void Separate(const HF_Channels &Chan);
};

//Structure for holding channel information
struct HF_Channels{
  int size1;
  int size3;
  State *qnums1;
  State *qnums3;
  int *indvec;
  int *ntb;
  int **tbvec;
  //std::unordered_map<std::string, int> tbmap;
  int *nob;
  int **obvec;
  HF_Channels(){};
  HF_Channels(const Input_Parameters &Parameters, const Model_Space &Space);
  void delete_struct();
};

struct HF_Matrix_Elements{
  double **V;
  HF_Matrix_Elements(const HF_Channels &Chan);
  HF_Matrix_Elements(){};
  void delete_struct(const HF_Channels &Chan);
};

#endif
