#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"

HF_Matrix_Elements::HF_Matrix_Elements(const HF_Channels &Chan)
{
  int ntb1;
  V = new double*[Chan.size1];
  for(int i = 0; i < Chan.size1; ++i){
    ntb1 = Chan.ntb[i];
    if(ntb1 != 0){
      ntb1 *= ntb1;
      V[i] = new double[ntb1];
      for(int j = 0; j < ntb1; ++j){
	V[i][j] = 0.0;
      }
    }
  }
}

void HF_Matrix_Elements::delete_struct(const HF_Channels &Chan)
{
  int ntb1;
  for(int i = 0; i < Chan.size1; ++i){
    ntb1 = Chan.ntb[i];
    if(ntb1 != 0){
      delete[] V[i];
    }
  }
  delete[] V;
}

//Function to setup Channels
HF_Channels::HF_Channels(const Input_Parameters &Parameters, const Model_Space &Space)
{
  std::cout << "Building HF Channels ... " << std::endl;

  State state;
  int count0, ind1, nob1, ntb1;
  State *qnumstemp = new State[Space.indtot]; // max number of qnums groups
  int *obnum = new int[Space.indtot]; // count states in each qnums group
  indvec = new int[Space.indtot]; // index of qnums group for each state
  for(int i = 0; i < Space.indtot; ++i){ obnum[i] = 0; }

  qnumstemp[0] = Space.qnums[0];
  indvec[0] = 0;
  obnum[0] = 1;

  // count # of qnums groups and # of hs and ps in each qnums group, and fill indvec
  count0 = 1;
  for(int i = 1; i < Space.indtot; ++i){
    state = Space.qnums[i];
    for(int k = 0; k < count0; ++k){
      if( equal(state, qnumstemp[k]) ){
	indvec[i] = k;
	++obnum[k];
	goto stop;
      }
      else if(k == count0 - 1){
	qnumstemp[count0] = Space.qnums[i];
	indvec[i] = count0;
	++obnum[count0];
	++count0;
	break;
      }
    }
  stop:;
  }
  size3 = count0;

  // allocate memory for Hvec and Pvec, reset hnum and pnum
  qnums3 = new State[size3];
  obvec = new int*[size3];
  nob = new int[size3];
  for(int i = 0; i < size3; ++i){
    nob1 = obnum[i];
    qnums3[i] = qnumstemp[i];
    if(nob1 != 0){ obvec[i] = new int[nob1]; }
    nob[i] = 0;
  }
  delete[] qnumstemp;
  delete[] obnum;

  // place states in appropriate Hvec or Pvec position
  for(int i = 0; i < Space.indtot; ++i){
    for(int k = 0; k < size3; ++k){
      if( equal(Space.qnums[i], qnums3[k]) ){
	obvec[k][nob[k]] = i;
	++nob[k];
	break;
      }
    }
  }

  /*for(int i = 0; i < size3; ++i){
    std::cout << qnums3[i].ml << " " << qnums3[i].m << std::endl;
    for(int j = 0; j < nob[i]; ++j){
      std::cout << obvec[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/

  size1 = Space.size_2b;
  qnums1 = new State[size1];
  tbvec = new int*[size1];
  ntb = new int[size1];
  for(int i = 0; i < size1; ++i){
    ntb[i] = 0;
  }

  for(int p = 0; p < Space.indtot; ++p){
    for(int q = 0; q < Space.indtot; ++q){
      plus(state, Space.qnums[p], Space.qnums[q]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      qnums1[ind1] = state;
      ++ntb[ind1];
    }
  }

  for(int i = 0; i < size1; ++i){
    ntb1 = ntb[i];
    if(ntb1 != 0){
      tbvec[i] = new int[2 * ntb[i]];
    }
    ntb[i] = 0;
  }

  for(int p = 0; p < Space.indtot; ++p){
    for(int q = 0; q < Space.indtot; ++q){
      plus(state, Space.qnums[p], Space.qnums[q]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      tbvec[ind1][2 * ntb[ind1]] = p;
      tbvec[ind1][2 * ntb[ind1] + 1] = q;
      ++ntb[ind1];
    }
  }

  /*for(int i = 0; i < size1; ++i){
    std::cout << qnums1[i].ml << " " << qnums1[i].m << std::endl;
    for(int j = 0; j < ntb[i]; ++j){
      std::cout << tbvec[i][2*j] << "," << tbvec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/
}

void HF_Channels::delete_struct()
{
  int ntb1;
  for(int i = 0; i < size3; ++i){
    delete[] obvec[i];
  }
  delete[] indvec;
  delete[] qnums3;
  delete[] obvec;
  delete[] nob;

  for(int i = 0; i < size1; ++i){
    ntb1 = ntb[i];
    if(ntb1 != 0){
      delete[] tbvec[i];
    }
  }
  delete[] qnums1;
  delete[] tbvec;
  delete[] ntb;
}


Single_Particle_States::Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan)
{
  int ind; // count for filling states
  int nob1, h1, p1;
  double tempen, tempen1, vec1; // temp energy for ordering states
  if(Parameters.Pshells != 0){ hp = Parameters.P; }
  else{ hp = 0; }
  if(Parameters.Nshells != 0){ hn = Parameters.N; }
  else{ hn = 0; }

  vectors = new double**[Chan.size3];
  energies = new double*[Chan.size3];
  h = new int[Chan.size3];
  p = new int[Chan.size3];
  holes = new double**[Chan.size3];
  particles = new double**[Chan.size3];
  h_energies = new double*[Chan.size3];
  pt_energies = new double*[Chan.size3];

  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    h[i] = 0;
    p[i] = 0;
    if(nob1 != 0){
      vectors[i] = new double*[nob1];
      energies[i] = new double[nob1];
      for(int j = 0; j < nob1; ++j){
	energies[i][j] = Space.qnums[Chan.obvec[i][j]].energy;
	vectors[i][j] = new double[nob1];
	if(Space.qnums[Chan.obvec[i][j]].type == "hole"){ ++h[i]; }
	else{ ++p[i]; }
	for(int k = 0; k < nob1; ++k){
	  if(j == k){ vectors[i][j][k] = 1.0; }
	  else{ vectors[i][j][k] = 0.0; }
	}
      }
    }
  }

  // Order states by energy
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nob[i] - 1; ++j){
      ind = j;
      tempen = energies[i][j];
      for(int k = j + 1; k < Chan.nob[i]; ++k){
	if(energies[i][k] < tempen){ tempen = energies[i][k]; ind = k; }
      }
      tempen1 = energies[i][j];
      energies[i][j] = energies[i][ind];
      energies[i][ind] = tempen1;
      for(int k = 0; k < Chan.nob[i]; ++k){
	vec1 = vectors[i][j][k];
	vectors[i][j][k] = vectors[i][ind][k];
	vectors[i][ind][k] = vec1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    h1 = h[i];
    p1 = p[i];
    if(h1 != 0){
      holes[i] = new double*[h1];
      h_energies[i] = new double[h1];
      for(int j = 0; j < h1; ++j){
	holes[i][j] = new double[nob1];
	h_energies[i][j] = energies[i][j];
	for(int k = 0; k < nob1; ++k){
	  holes[i][j][k] = vectors[i][j][k];
	}
      }
    }
    if(p1 != 0){
      particles[i] = new double*[p1];
      pt_energies[i] = new double[p1];
      for(int j = 0; j < p1; ++j){
	particles[i][j] = new double[nob1];
	pt_energies[i][j] = energies[i][j + h1];
	for(int k = 0; k < nob1; ++k){
	  particles[i][j][k] = vectors[i][j + h1][k];
	}
      }
    }
  }

  //Separate States
  Separate(Chan);
}

void Single_Particle_States::delete_struct(const HF_Channels &Chan)
{
  int h1, p1, nob1;
  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    h1 = h[i];
    p1 = p[i];
    if(nob1 != 0){
      for(int j = 0; j < nob1; ++j){
	delete[] vectors[i][j];
      }
      delete[] vectors[i];
      delete[] energies[i];
    }
    if(h1 != 0){
      for(int j = 0; j < h1; ++j){
	delete[] holes[i][j];
      }
      delete[] holes[i];
      delete[] h_energies[i];
    }
    if(p1 != 0){
      for(int j = 0; j < p1; ++j){
	delete[] particles[i][j];
      }
      delete[] particles[i];
      delete[] pt_energies[i];
    }
  }
  delete[] vectors;
  delete[] energies;
  delete[] h;
  delete[] holes;
  delete[] h_energies;
  delete[] p;
  delete[] particles;
  delete[] pt_energies;
}


//ignores existing holes and particles and fills them from protons and neutrons
void Single_Particle_States::Separate(const HF_Channels &Chan)
{
  //define holes/particles
  int h1, p1, nob1;
  for(int i = 0; i < Chan.size3; ++i){
    h1 = h[i];
    p1 = p[i];
    nob1 = Chan.nob[i];
    if(h1 != 0){
      for(int j = 0; j < h1; ++j){
	delete[] holes[i][j];
      }
      delete[] holes[i];
      delete[] h_energies[i];
      holes[i] = new double*[h1];
      h_energies[i] = new double[h1];
      for(int j = 0; j < h1; ++j){
	holes[i][j] = new double[nob1];
	h_energies[i][j] = energies[i][j];
	for(int k = 0; k < nob1; ++k){
	  holes[i][j][k] = vectors[i][j][k];
	}
      }
    }
    if(p1 != 0){
      for(int j = 0; j < p1; ++j){
	delete[] particles[i][j];
      }
      delete[] particles[i];
      delete[] pt_energies[i];
      particles[i] = new double*[p1];
      pt_energies[i] = new double[p1];
      for(int j = 0; j < p1; ++j){
	particles[i][j] = new double[nob1];
	pt_energies[i][j] = energies[i][j + h1];
	for(int k = 0; k < nob1; ++k){
	  particles[i][j][k] = vectors[i][j + h1][k];
	}
      }
    }
  }
}


void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  size_t index1, index2; // indicies for finding parameters among file lines
  double TBME, hom, r2, p2; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4, coupJ, coupT, par; // interaction file contents
  int ind1, ind;
  State tb;

  ME = HF_Matrix_Elements(Chan);

  fullpath1 = PATH + Parameters.MatrixElements + ".int";

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << Parameters.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> number;
  while (number != "Total"){ 
    getline(interaction, interactionline);
    interactionstream.str(interactionline);
    interactionstream >> number;
  }
  index1 = interactionline.find_first_of("0123456789");
  index2 = interactionline.find_last_of("0123456789");
  NumElements = std::atoi( interactionline.substr(index1, index2 - index1 + 1).c_str() );

  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> number;
  while (number != "Tz"){
    getline(interaction, interactionline);
    interactionstream.str(interactionline);
    interactionstream >> number;
  }

  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> coupT >> par >> coupJ >> shell1 >> shell2 >> shell3 >> shell4 >> TBME >> hom >> r2 >> p2;
    //std::cout << coupT << " " << par << " " << coupJ << " " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << ", " << TBME << std::endl;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    if((shell1 == shell2 || shell3 == shell4) && coupJ%2 != 0){ continue; }
    if(shell1 == shell2){ TBME *= sqrt(2.0); } // !! check
    if(shell3 == shell4){ TBME *= sqrt(2.0); } // !! check
    //if(shell1 > shell2){ std::swap(shell1, shell2); TBME *= pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ))); }
    //if(shell3 > shell4){ std::swap(shell3, shell4); TBME *= pow(-1.0, int(0.5*(Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ))); }
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    tb.j = coupJ;
    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);

    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell1, shell2, shell3, shell4);
    ME.V[ind1][ind] = TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell2, shell1, shell3, shell4);
    ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ) + 1)) * TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell1, shell2, shell4, shell3);
    ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ) + 1)) * TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell2, shell1, shell4, shell3);
    ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j))) * TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell3, shell4, shell1, shell2);
      ME.V[ind1][ind] = TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell3, shell4, shell2, shell1);
      ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ) + 1)) * TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell4, shell3, shell1, shell2);
      ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ) + 1)) * TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell4, shell3, shell2, shell1);
      ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j))) * TBME;
    }
  }
  interaction.close();
}

void Read_Matrix_Elements_M(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind1, ind;
  State tb;
  ME = HF_Matrix_Elements(Chan);

  fullpath1 = PATH + Parameters.MatrixElements + ".int";

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << Parameters.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> NumElements;

  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
    TBME *= Parameters.tbstrength;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    if(shell1 == shell2 || shell3 == shell4){ continue; }
    if(shell1 > shell2){ std::swap(shell1, shell2); TBME *= -1.0; }
    if(shell3 > shell4){ std::swap(shell3, shell4); TBME *= -1.0; }
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);

    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell1, shell2, shell3, shell4);
    ME.V[ind1][ind] = TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell2, shell1, shell3, shell4);
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell1, shell2, shell4, shell3);
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell2, shell1, shell4, shell3);
    ME.V[ind1][ind] = TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell3, shell4, shell1, shell2);
      ME.V[ind1][ind] = TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell3, shell4, shell2, shell1);
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell4, shell3, shell1, shell2);
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell4, shell3, shell2, shell1);
      ME.V[ind1][ind] = TBME;
    }
  }
  interaction.close();
}

void Read_HO_ME_From_File(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  std::ifstream interaction;	// interaction file
  ME = HF_Matrix_Elements(Chan);

  /*int p, q, r, s, ntb;
  for(int chan = 0; chan < Chan.size1; ++chan){
    ntb = Chan.ntb[chan];
    if(ntb == 0){ continue; }
    
    for(int pq = 0; pq < ntb; ++pq){
      p = Chan.tbvec[chan][2*pq];
      q = Chan.tbvec[chan][2*pq + 1];
      for(int rs = 0; rs < ntb; ++rs){
	r = Chan.tbvec[chan][2*rs];
	s = Chan.tbvec[chan][2*rs + 1];
	std::cout << "! " << p << " " << q << " | " << r << " " << s << " : " << ME.V[chan][pq*ntb + rs] << std::endl;;
      }
    }
    }*/

  fullpath1 = PATH + "coulomb-ho2d-elements-20-shells.dat";

  interaction.open(fullpath1.c_str(), std::ios::binary);
  if(interaction.is_open()){
    //get length of file
    interaction.seekg(0, interaction.end);
    int length = interaction.tellg();
    interaction.seekg(0, interaction.beg);

    char *buffer = new char[length];
    interaction.read(buffer, length);
    interaction.close();

    float TBME;
    unsigned int neg1 = 128;
    unsigned int neg2 = 4294967040;
    unsigned int n1, n2, n3, n4;
    int ml1, ml2, ml3, ml4;
    int nmax = Parameters.Shells;
    int key;
    State tb;
    int ind, ind1;
    State statep, stateq, stater, states;
    int p, q, r, s;
    //int S1, S2;
    statep.t = -1, stateq.t = -1, stater.t = -1, states.t = -1;
    length /= 16;

    #pragma omp parallel private(n1, n2, n3, n4, ml1, ml2, ml3, ml4, key, p, q, r, s, statep, stateq, stater, states, TBME, ind, ind1, tb)
    {
      #pragma omp for schedule(static)
      for(int i = 0; i < length; ++i){
	n1 = 0, n2 = 0, n3 = 0, n4 = 0;
	ml1 = 0, ml2 = 0, ml3 = 0, ml4 = 0;
	
	n1 = *(buffer + 16*i);
	ml1 = *(buffer + 16*i + 1);
	if((neg1 & ml1) != 0){ ml1 = (ml1 | neg2); }// ml1 < 0
	n2 = *(buffer + 16*i + 2);
	ml2 = *(buffer + 16*i + 3);
	if((neg1 & ml2) != 0){ ml2 = (ml2 | neg2); }// ml2 < 0
	n3 = *(buffer + 16*i + 4);
	ml3 = *(buffer + 16*i + 5);
	if((neg1 & ml3) != 0){ ml3 = (ml3 | neg2); }// ml3 < 0
	n4 = *(buffer + 16*i + 6);
	ml4 = *(buffer + 16*i + 7);
	if((neg1 & ml4) != 0){ ml4 = (ml4 | neg2); }// ml4 < 0
	TBME = *(double*)(buffer + 16*i + 8);
	if(int(2*n1 + abs(ml1)) >= nmax || int(2*n2 + abs(ml2)) >= nmax || int(2*n3 + abs(ml3)) >= nmax || int(2*n4 + abs(ml4)) >= nmax){ continue; }
	//std::cout << "< " << n1 << "," << ml1 << " ; " << n2 << "," << ml2 << " |V| " << n3 << "," << ml3 << " ; " << n4 << "," << ml4 << " > = " << TBME << std::endl;

	TBME *= std::sqrt(Parameters.density);
	statep.n = n1;
	statep.ml = ml1;
	stateq.n = n2;
	stateq.ml = ml2;
	stater.n = n3;
	stater.ml = ml3;
	states.n = n4;
	states.ml = ml4;
	for(int s1 = -1; s1 <= 1; s1 += 2){
	  statep.m = s1;
	  key = ChanInd_1b(Parameters.basis, Space, statep);
	  p = Space.map_1b[key];
	  for(int s2 = -1; s2 <= 1; s2 += 2){
	    stateq.m = s2;
	    key = ChanInd_1b(Parameters.basis, Space, stateq);
	    q = Space.map_1b[key];
	    if(p == q){ continue; }
	    for(int s3 = -1; s3 <= 1; s3 += 2){
	      stater.m = s3;
	      key = ChanInd_1b(Parameters.basis, Space, stater);
	      r = Space.map_1b[key];
	      if(s3 != s1){ continue; }
	      for(int s4 = -1; s4 <= 1; s4 += 2){
		states.m = s4;
		key = ChanInd_1b(Parameters.basis, Space, states);
		s = Space.map_1b[key];
		if(r == s || s4 != s2){ continue; }
		
		plus(tb, Space.qnums[p], Space.qnums[q]);
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);

		//std::cout << "1--------------------------------------------------------------" << std::endl;
		//std::cout << "? <" << p << " " << q << "||" << r << " " << s << "> += <" << p << " " << q << "|" << r << " " << s << ">" << std::endl;
		//std::cout << "? <" << p << " " << q << "||" << s << " " << r << "> -= <" << p << " " << q << "|" << r << " " << s << ">" << std::endl;
		// C(p1q2r3s4) -> <p1q2 || r3s4>
		ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], p, q, r, s);
		ME.V[ind1][ind] += TBME;
		// C(p1q2r3s4) -> -<p1q2 || s4r3>
		ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], p, q, s, r);
		ME.V[ind1][ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  //std::cout << "2--------------------------------------------------------------" << std::endl;
		  //std::cout << "? <" << q << " " << p << "||" << s << " " << r << "> += <" << q << " " << p << "|" << s << " " << r << ">" << std::endl;
		  //std::cout << "? <" << q << " " << p << "||" << r << " " << s << "> -= <" << q << " " << p << "|" << s << " " << r << ">" << std::endl;
		  // C(p1q2r3s4) = C(q2p1s4r3) -> <q2p1 || s4r3>
		  ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], q, p, s, r);
		  ME.V[ind1][ind] += TBME;
		  // C(p1q2r3s4) = C(q2p1s4r3) -> -<q2p1 || r3s4>
		  ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], q, p, r, s);
		  ME.V[ind1][ind] -= TBME;
		}
		if(((n1 == n3 && ml1 == ml3) && (n2 == n4 && ml2 == ml4)) || ((n1 == n4 && ml1 == ml4) && (n2 == n3 && ml2 == ml3))){ continue; }
		//std::cout << "3--------------------------------------------------------------" << std::endl;
		//std::cout << "? <" << r << " " << s << "||" << p << " " << q << "> += <" << r << " " << s << "|" << p << " " << q << ">" << std::endl;
		//std::cout << "? <" << r << " " << s << "||" << q << " " << p << "> -= <" << r << " " << s << "|" << p << " " << q << ">" << std::endl;
		// C(p1q2r3s4) = C(r3s4p1q2) -> <r3s4 || p1q2>
		ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], r, s, p, q);
		ME.V[ind1][ind] += TBME;
		// C(p1q2r3s4) = C(r3s4p1q2) -> -<r3s4 || q2p1>
		ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], r, s, q, p);
		ME.V[ind1][ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  //std::cout << "4--------------------------------------------------------------" << std::endl;
		  //std::cout << "? <" << s << " " << r << "||" << q << " " << p << "> += <" << s << " " << r << "|" << q << " " << p << ">" << std::endl;
		  //std::cout << "? <" << s << " " << r << "||" << p << " " << q << "> -= <" << s << " " << r << "|" << q << " " << p << ">" << std::endl;
		  // C(p1q2r3s4) = C(s4r3q2p1) -> <s4r3 || q2p1>
		  ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], s, r, q, p);
		  ME.V[ind1][ind] += TBME;
		  // C(p1q2r3s4) = C(s4r3q2p1) -> -<s4r3 || p1q2>
		  ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], s, r, p, q);
		  ME.V[ind1][ind] -= TBME;
		  //std::cout << "--------------------------------------------------------------" << std::endl;
		}
	      }
	    }
	  }
	}
      }
    }
    delete[] buffer;
  }

  /*for(int chan = 0; chan < Chan.size1; ++chan){
    ntb = Chan.ntb[chan];
    if(ntb == 0){ continue; }
    
    for(int pq = 0; pq < ntb; ++pq){
      p = Chan.tbvec[chan][2*pq];
      q = Chan.tbvec[chan][2*pq + 1];
      for(int rs = 0; rs < ntb; ++rs){
	r = Chan.tbvec[chan][2*rs];
	s = Chan.tbvec[chan][2*rs + 1];
	std::cout << "! " << p << " " << q << " | " << r << " " << s << " : " << ME.V[chan][pq*ntb + rs] << std::endl;;
      }
    }
    }*/
}

void Read_Matrix_Elements_HO(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  struct timespec time1, time2;
  double elapsed0 = 0.0;

  std::cout << "Computing Coulomb Matrix Elements ..." << std::endl;
  ME = HF_Matrix_Elements(Chan);
  clock_gettime(CLOCK_MONOTONIC, &time1);

  int pq, rs, p, q, r, s;
  int ntb, ntb0;
  int ind1, ind2, ind3, ind4;
  int length, length0;
  int *tbvec0;
  double TBME;
  for(int chan = 0; chan < Chan.size1; ++chan){
    ntb = Chan.ntb[chan];
    if(ntb == 0){ continue; }
    tbvec0 = new int[ntb];
    ntb0 = 0;
    for(int i = 0; i < ntb; ++i){
      if(Chan.tbvec[chan][2*i] < Chan.tbvec[chan][2*i + 1]){
	tbvec0[2*ntb0] = Chan.tbvec[chan][2*i];
	tbvec0[2*ntb0 + 1] = Chan.tbvec[chan][2*i + 1];
	++ntb0;
      }
    }
    
    length = int(0.5 * ntb0 * (ntb0 + 1));
    #pragma omp parallel private(pq, rs, p, q, r, s, length0, ind1, ind2, ind3, ind4, TBME)
    {
      #pragma omp for schedule(static)
      for(int pqrs = 0; pqrs < length; ++pqrs){
	pq = std::floor((2*ntb0 - 1 - std::sqrt(1 + 4*ntb0 + 4*ntb0*ntb0 - 8*pqrs))/2) + 1;
	length0 = int(0.5 * pq * (2*ntb0 - pq + 1));
	rs = int(pq + pqrs - length0);
	p = tbvec0[2*pq];
	q = tbvec0[2*pq + 1];
	r = tbvec0[2*rs];
	s = tbvec0[2*rs + 1];
	ind1 = Index2(Chan.tbvec[chan], Chan.ntb[chan], p, q);
	ind2 = Index2(Chan.tbvec[chan], Chan.ntb[chan], r, s);
	ind3 = Index2(Chan.tbvec[chan], Chan.ntb[chan], q, p);
	ind4 = Index2(Chan.tbvec[chan], Chan.ntb[chan], s, r);
	TBME = Coulomb_HO(Parameters, Space, p, q, r, s);
	ME.V[chan][ind1 * ntb + ind2] = TBME;
	ME.V[chan][ind3 * ntb + ind2] = -1.0 * TBME;
	ME.V[chan][ind1 * ntb + ind4] = -1.0 * TBME;
	ME.V[chan][ind3 * ntb + ind4] = TBME;
	if(ind1 != ind2){
	  ME.V[chan][ind2 * ntb + ind1] = TBME;
	  ME.V[chan][ind4 * ntb + ind1] = -1.0 * TBME;
	  ME.V[chan][ind2 * ntb + ind3] = -1.0 * TBME;
	  ME.V[chan][ind4 * ntb + ind3] = TBME;
	}
      }
    }
  }

  clock_gettime(CLOCK_MONOTONIC, &time2);
  elapsed0 = (time2.tv_sec - time1.tv_sec);
  elapsed0 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
  std::cout << std::endl << "!! Runtime = " << elapsed0 << " sec. " << std::endl;

}

void Hartree_Fock_States(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &HF, const HF_Matrix_Elements &ME)
{
  std::cout << "Diagonalizing Hartree-Fock Matrix ..." << std::endl;
  double error; // Energy error between iterations
  double error1; // Energy error between iterations
  double Bshift; // level shift parameter
  int ind; // Index to keep track of iteration number
  char jobz, uplo; // Parameters for Diagonalization, Multiplication
  int lda; // Parameter for Diagonalization
  int lwork, info; // Parameters for Diagonaliztion
  double term;
  State tb;
  int nob1;
  double tempen2;
  double tempen3;
  double vec3;

  int *temph;
  double **tempen;
  double ***tempvec;
  temph = new int[Chan.size3];
  tempen = new double*[Chan.size3];
  tempvec = new double**[Chan.size3];
  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    temph[i] = HF.h[i];
    tempen[i] = new double[nob1];
    tempvec[i] = new double*[nob1];
    for(int j = 0; j < nob1; ++j){
      tempen[i][j] = HF.energies[i][j];
      tempvec[i][j] = new double[nob1];
      for(int k = 0; k < nob1; ++k){
	tempvec[i][j][k] = HF.vectors[i][j][k];
      }
    }
  }

  double *w;
  double *work;
  double *fock;

  jobz = 'V';
  uplo = 'U';
  Bshift = 50.0;

  ind = 0;
  error = 1000;

  while((error > 1e-10 && ind < 5000) || ind < 10){ //10
    ++ind;
    error = 0.0;

    int size1, size2;
    int ind1, ind2;
    int mind, nind, lind, kind;
    int m, l;
    int length, length0;
    //Make Fock Matrix
    for(int i = 0; i < Chan.size3; ++i){
      error1 = 0.0;
      size1 = Chan.nob[i];
      fock = new double[size1 * size1];
      length = int(0.5 * size1 * (size1 + 1));
      #pragma omp parallel private(mind, lind, size2, nind, tb, ind1, kind, ind2, term, m, l, length0)
      {
        #pragma omp for schedule(static)
	for(int ml = 0; ml < length; ++ml){
	  m = std::floor((2*size1 - 1 - std::sqrt(1 + 4*size1 + 4*size1*size1 - 8*ml))/2) + 1;
	  length0 = int(0.5 * m * (2*size1 - m + 1));
	  l = int(m + ml - length0);

	  mind = Chan.obvec[i][m];
	  lind = Chan.obvec[i][l];
	  if(m == l){ fock[size1 * m + m] = Space.qnums[mind].energy; }	// Add diagonal elements to fock matrices
	  else{ fock[size1 * m + l] = 0.0; }

	  for(int j = 0; j < Chan.size3; ++j){
	    size2 = Chan.nob[j];
	    for(int beta = 0; beta < temph[j]; ++beta){ // Sum over occupied levels
	      for(int n = 0; n < size2; ++n){
		nind = Chan.obvec[j][n];
		plus(tb, Space.qnums[mind], Space.qnums[nind]);
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		for(int k = 0; k < size2; ++k){
		  kind = Chan.obvec[j][k];
		  ind2 = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], mind, nind, lind, kind);
		  term = tempvec[j][beta][n] * tempvec[j][beta][k] * ME.V[ind1][ind2];
		  fock[size1 * m + l] += term;
		}
	      }
	    }
	  }
	  for(int beta = 0; beta < temph[i]; ++beta){ // Sum over occupied levels
	    fock[size1 * m + l] -= tempvec[i][beta][m] * tempvec[i][beta][l] * Bshift;
	  }
	  if(m != l){ fock[size1 * l + m] = fock[size1 * m + l]; }
	}
      }

      /*std::cout << "!! " << ind << ", " << i << " " << Chan.qnums3[i].ml << " " << Chan.qnums3[i].m << " " << Chan.qnums3[i].par << " " << Chan.qnums3[i].t << std::endl;
      for(int j = 0; j < size1; ++j){
	for(int k = 0; k < size1; ++k){
	  std::cout << fock[size1*j + k] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/

      lda = size1;
      lwork = (3+2)*size1;
      w = new double[lda];
      work = new double[lwork];
      for(int j = 0; j < size1; ++j){ w[j] = 0.0; }
      for(int j = 0; j < (3+2)*size1; ++j){ work[j] = 0.0; }
    
      if(size1 != 0){ dsyev_(&jobz, &uplo, &size1, fock, &lda, w, work, &lwork, &info); }
      for(int j = 0; j < size1*size1; ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}
      for(int j = 0; j < temph[i]; ++j){ w[j] += Bshift; } //Add back Level-shift parameter

      /*for(int j = 0; j < size1; ++j){
	for(int k = 0; k < size1; ++k){
	  std::cout << fock[size1*j + k] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/

      for(int j = 0; j < size1; ++j){
	for(int k = 0; k < size1; ++k){
	  if(fabs(fock[size1 * j + k]) > 1e-10){
	    error1 += fabs((HF.vectors[i][j][k] - fock[size1 * j + k])/fock[size1 * j + k]);
	  }
	  HF.vectors[i][j][k] = fock[size1 * j + k];
	}
      }

      // Order states by energy
      for(int j = 0; j < size1 - 1; ++j){
	ind2 = j;
	tempen2 = w[j];
	for(int k = j + 1; k < size1; ++k){
	  if(w[k] < tempen2){ tempen2 = w[k]; ind2 = k; }
	}
	tempen3 = w[j];
	w[j] = w[ind2];
	w[ind2] = tempen3;
	for(int k = 0; k < size1; ++k){
	  vec3 = HF.vectors[i][j][k];
	  HF.vectors[i][j][k] = HF.vectors[i][ind2][k];
	  HF.vectors[i][ind2][k] = vec3;
	}
      }

      for(int j = 0; j < size1; ++j){
	HF.energies[i][j] = w[j];
      }
      GramSchmidt(HF.vectors[i], size1);

      delete[] fock;
      delete[] w;
      delete[] work;

      error += (error1 / (size1*size1));
    }
    HF.Separate(Chan);

    // HFtemp = HF
    for(int i = 0; i < Chan.size3; ++i){
      temph[i] = HF.h[i];
      nob1 = Chan.nob[i];
      for(int j = 0; j < nob1; ++j){
	tempen[i][j] = HF.energies[i][j];
	for(int k = 0; k < nob1; ++k){
	  tempvec[i][j][k] = HF.vectors[i][j][k];
	}
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    for(int j = 0; j < nob1; ++j){
      delete[] tempvec[i][j];
    }
    delete[] tempen[i];
    delete[] tempvec[i];
  }
  delete[] temph;
  delete[] tempen;
  delete[] tempvec;

  ind = 0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < HF.h[i]; ++j){
      Space.qnums[Chan.obvec[i][j]] = Chan.qnums3[i];
      Space.qnums[Chan.obvec[i][j]].energy = HF.energies[i][j];
      Space.qnums[Chan.obvec[i][j]].type = "hole";
      ++ind;
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = HF.h[i]; j < HF.h[i] + HF.p[i]; ++j){
      Space.qnums[Chan.obvec[i][j]] = Chan.qnums3[i];
      Space.qnums[Chan.obvec[i][j]].energy = HF.energies[i][j];
      Space.qnums[Chan.obvec[i][j]].type = "particle";
      ++ind;
    }
  }

  /*for(int i = 0; i < Space.indtot; ++i){
    std::cout << Space.qnums[i].par << " " << Space.qnums[i].t << " " << Space.qnums[i].ml << " " << Space.qnums[i].m << " " << Space.qnums[i].j << " : " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
    }*/


  /*std::ofstream HFlevelfile;
  std::string filename = PATH + Parameters.LevelScheme + "_HF.sp";
  HFlevelfile.open(filename.c_str());
  HFlevelfile << "Mass number A of chosen nucleus (important for CoM corrections): \t" << Space.A << "\n";
  HFlevelfile << "Oscillator energy: \t" << Space.HOEnergy << "\n";
  HFlevelfile << "Total number of single-particle orbits: \t" << Space.shelltot << "\n";
  HFlevelfile << "Legend:   \tn \tl \t2j \ttz \t2n+l \tHO-energy \tevalence \tparticle/hole \tinside/outside \n";
  for(int i = 0; i < int(HF.protons.size()); ++i){
    HFlevelfile << "Number:   " << i+1 << "\t" << p_n[i] << "\t" << HF.p_l[i] << "\t" << int(2*HF.p_j[i]) << "\t";
    HFlevelfile << "-1" << "\t" << 2*p_n[i]+HF.p_l[i] << "\t" << std::setprecision(8) << HF.p_energies[i] << "\t" << "0.000000" << "\t";
    if(i < Space.Pocc){ HFlevelfile << "hole    " << "\t" << "inside" << "\n"; }
    else{ HFlevelfile << "particle" << "\t" << "inside" << "\n"; }
  }
  for(int i = 0; i < int(HF.neutrons.size()); ++i){
    HFlevelfile << "Number:   " << i+int(HF.protons.size())+1 << "\t" << n_n[i] << "\t" << HF.n_l[i] << "\t" << int(2*HF.n_j[i]) << "\t";
    HFlevelfile << "1" << "\t" << 2*n_n[i]+HF.n_l[i] << "\t" << std::setprecision(8) << HF.n_energies[i] << "\t" << "0.000000" << "\t";
    if(i < Space.Nocc){ HFlevelfile << "hole    " << "\t" << "inside" << "\n"; }
    else{ HFlevelfile << "particle" << "\t" << "inside" << "\n"; }
  }
  
  HFlevelfile.close();*/
  
}


void Convert_To_HF_Matrix_Elements(const HF_Channels &Chan, const Single_Particle_States &States, HF_Matrix_Elements &ME)
{
  std::cout << "Converting Matrix Elements to HF basis ..." << std::endl;

  int length, matlength; // max length of M-Scheme indicies, length of J_ME
  double tempel;
  double *M1;//, *M2; // Matrices of coefficients
  double *C;

  char transa, transb;
  double alpha1, beta1;
  std::ofstream jschemefile; // file to print M-Scheme matrix elements
  int p, q, a, g;
  int ind1, ind2;

  for(int chan = 0; chan < Chan.size1; ++chan){
    length = Chan.ntb[chan];
    if(length == 0){ continue; }
    matlength = pow(length, 2.0);

    M1 = new double[matlength];
    C = new double[matlength];
    for(int i = 0; i < matlength; ++i){
      M1[i] = 0.0;
      C[i] = 0.0;
    }

    #pragma omp parallel private(p, q, a, g, ind1, ind2, tempel)
    {
      #pragma omp for schedule(static)
      for(int pq = 0; pq < length; ++pq){
	for(int ag = 0; ag < length; ++ag){
	  p = Chan.tbvec[chan][2*pq];
	  q = Chan.tbvec[chan][2*pq + 1];
	  a = Chan.tbvec[chan][2*ag];
	  g = Chan.tbvec[chan][2*ag + 1];
	  if(Chan.indvec[p] != Chan.indvec[a] || Chan.indvec[q] != Chan.indvec[g]){ continue; }
	  ind1 = Chan.indvec[p];
	  ind2 = Chan.indvec[q];
	  p = Index1(Chan.obvec[ind1], Chan.nob[ind1], p);
	  q = Index1(Chan.obvec[ind2], Chan.nob[ind2], q);
	  a = Index1(Chan.obvec[ind1], Chan.nob[ind1], a);
	  g = Index1(Chan.obvec[ind2], Chan.nob[ind2], g);
	  tempel = States.vectors[ind1][a][p] * States.vectors[ind2][g][q];
	  M1[length * pq + ag] = tempel;
	}
      }
    }

    /*std::cout << "!! ";
    for(int i = 0; i < length; ++i){
      std::cout << Chan.tbvec[chan][2*i] << Chan.tbvec[chan][2*i + 1] << " ";
    }
    std::cout << std::endl << std::endl;

    for(int i = 0; i < length; ++i){
      for(int j = 0; j < length; ++j){
	std::cout << M1[length*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    for(int i = 0; i < length; ++i){
      for(int j = 0; j < length; ++j){
	std::cout << ME.V[chan][length*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/

    transa = 'N';
    transb = 'T';
    alpha1 = 1.0;
    beta1 = 0.0;
    
    dgemm_NN(ME.V[chan], M1, C, &length, &length, &length, &alpha1, &beta1, &transa, &transa);
    dgemm_TN(M1, C, ME.V[chan], &length, &length, &length, &alpha1, &beta1, &transb, &transa);

    /*for(int i = 0; i < length; ++i){
      for(int j = 0; j < length; ++j){
	std::cout << ME.V[chan][length*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/

    delete[] M1;
    delete[] C;
  }

  // print J_ME to file
  /*jschemefile.open((PATH + MatrixElements + "_HF.int").c_str());
  jschemefile << "Total number of twobody matx elements:" << "\t" << ind << "\n";
  jschemefile << "----> Interaction part\n";   
  jschemefile << "Nucleon-Nucleon interaction model:n3lo\n";            
  jschemefile << "Type of calculation: nocore\n";              
  jschemefile << "Number and value of starting energies:   1	0.000000E+00\n";
  jschemefile << "Total number of twobody matx elements:\t" << ind << "\n";
  jschemefile << "Tz      Par      2J      a      b      c      d      <ab|V|cd>\n";
  for(int J = 0; J <= Space.max2J; ++J){
    for(int m = 0; m < plength; ++m){
      for(int n = m; n < plength; ++n){
	for(int l = m; l < plength; ++l){
	  for(int k = l; k < plength; ++k){
	    tempel = HF_ME.get_ppJME(m, n, l, k, J);
	    if(fabs(tempel) < 1.0e-10){ continue; }
	    jschemefile << "\t -1\t" << pow(-1.0,States.p_l[m]+States.p_l[n]) << "\t" << 2*J << "\t"
			<< m + 1 << "\t" << n + 1 << "\t" << l + 1 << "\t" << k + 1
			<< "\t" << std::setprecision(8) << tempel << "\n";
	  }
	}
      }
    }
    for(int m = 0; m < plength; ++m){
      for(int n = 0; n < nlength; ++n){
	for(int l = m; l < plength; ++l){
	  for(int k = 0; k < nlength; ++k){
	    tempel = HF_ME.get_pnJME(m, n, l, k, J);
	    if(fabs(tempel) < 1.0e-10){ continue; }
	    jschemefile << "\t  0\t" << pow(-1.0,States.p_l[m]+States.n_l[n]) << "\t" << 2*J << "\t"
			<< m + 1 << "\t" << n + plength + 1 << "\t" << l + 1 << "\t" << k + plength + 1 
			<< "\t" << std::setprecision(8) << tempel << "\n";
	  }
	}
      }
    }
    for(int m = 0; m < nlength; ++m){
      for(int n = m; n < nlength; ++n){
	for(int l = m; l < nlength; ++l){
	  for(int k = l; k < nlength; ++k){
	    tempel = HF_ME.get_nnJME(m, n, l, k, J);
	    if(fabs(tempel) < 1.0e-10){ continue; }
	    jschemefile << "\t  1\t" << pow(-1.0,States.n_l[m]+States.n_l[n]) << "\t" << 2*J << "\t"
			<< m + plength + 1 << "\t" << n + plength + 1 << "\t" << l + plength + 1 << "\t"
			<< k + plength + 1 << "\t" << std::setprecision(8) << tempel << "\n";
	  }
	}
      }
    }
  }
  
  jschemefile.close();*/
}

void Get_Matrix_Elements(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind, ind1, ind2;
  std::string ptype, qtype, rtype, stype;
  State tb;

  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    #pragma omp parallel private(shell1, shell2, shell3, shell4, ptype, qtype, rtype, stype, ind, ind1, ind2, tb, TBME)
    {
      #pragma omp for schedule(static)
      for(int tb1 = 0; tb1 < HF_Chan.ntb[chan]; ++tb1){
	shell1 = HF_Chan.tbvec[chan][2*tb1];
	shell2 = HF_Chan.tbvec[chan][2*tb1 + 1];
	ptype = Space.qnums[shell1].type;
	qtype = Space.qnums[shell2].type;
	if(shell1 >= shell2){ continue; }
	for(int tb2 = 0; tb2 < HF_Chan.ntb[chan]; ++tb2){
	  shell3 = HF_Chan.tbvec[chan][2*tb2];
	  shell4 = HF_Chan.tbvec[chan][2*tb2 + 1];
	  rtype = Space.qnums[shell3].type;
	  stype = Space.qnums[shell4].type;
	  if(shell3 >= shell4){ continue; }
	  if(shell1 == shell3 && shell2 > shell4){ continue; }
	  
	  if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "hole"){ continue; }
	  if(ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "particle"){ continue; }
	  
	  TBME = HF_ME.V[chan][tb1*HF_Chan.ntb[chan] + tb2];
	  //std::cout << "V_hf: " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << " = " << TBME << std::endl;
	  if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V2[ind1][ind] = TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V2[ind1][ind] = TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell3, shell4, shell1, shell2);
	    Ints.D_ME1.V2[ind1][ind] = TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell4, shell3, shell1, shell2);
	    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell3, shell4, shell2, shell1);
	    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell4, shell3, shell2, shell1);
	    Ints.D_ME1.V2[ind1][ind] = TBME;
	  }
	  else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V1[ind1][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V1[ind1][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell3, shell4, shell1, shell2);
	    Ints.D_ME1.V1[ind1][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell4, shell3, shell1, shell2);
	    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell3, shell4, shell2, shell1);
	    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell4, shell3, shell2, shell1);
	    Ints.D_ME1.V1[ind1][ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp2vec[ind1], Chan.nhp2[ind1], Chan.nhp2[ind1], shell1, shell4, shell3, shell2);
	    Ints.D_ME1.V3[ind1][ind] = TBME;
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp2vec[ind1], Chan.nhp2[ind1], Chan.nhp2[ind1], shell3, shell2, shell1, shell4);
	    Ints.D_ME1.V3[ind1][ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell3, shell4, shell1, shell2);
	    Ints.D_ME1.V4[ind1][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell4, shell3, shell1, shell2);
	    Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell3, shell4, shell2, shell1);
	    Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell4, shell3, shell2, shell1);
	    Ints.D_ME1.V4[ind1][ind] = TBME;
	    
	    ind2 = Chan.indvec[shell2];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V5[ind2][ind] = TBME;
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
	    ind2 = Chan.indvec[shell1];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V5[ind2][ind] = TBME;
	    
	    ind2 = Chan.indvec[shell1];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V6[ind2][ind] = TBME;
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
	    ind2 = Chan.indvec[shell2];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V6[ind2][ind] = TBME;
	    
	    ind2 = Chan.indvec[shell4];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell1, shell2, shell3);
	    Ints.D_ME1.V7[ind2][ind] = TBME;
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell2, shell1, shell3);
	    Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
	    ind2 = Chan.indvec[shell3];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	    Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	    Ints.D_ME1.V7[ind2][ind] = TBME;
	    
	    ind2 = Chan.indvec[shell3];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	    Ints.D_ME1.V8[ind2][ind] = TBME;
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	    Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
	    ind2 = Chan.indvec[shell4];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell1, shell2, shell3);
	    Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell2, shell1, shell3);
	    Ints.D_ME1.V8[ind2][ind] = TBME;
	    
	    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell3, shell2, shell4);
	    Ints.D_ME1.V9[ind1][ind] = TBME;
	    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell3, shell1, shell4);
	    Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
	    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell4, shell2, shell3);
	    Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
	    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell4, shell1, shell3);
	    Ints.D_ME1.V9[ind1][ind] = TBME;
	    
	    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell4, shell2, shell3);
	    Ints.D_ME1.V10[ind1][ind] = TBME;
	    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell4, shell1, shell3);
	    Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
	    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell3, shell2, shell4);
	    Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
	    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell3, shell1, shell4);
	    Ints.D_ME1.V10[ind1][ind] = TBME;
	  }
	  if(Parameters.approx == "singles"){
	    if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	      ind2 = Chan.indvec[shell2];
	      ind = Index13(Chan.pvec[ind2], Chan.hppvec[ind2], Chan.np[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	      Ints.S_ME1.V11[ind2][ind] = TBME;
	      ind = Index13(Chan.pvec[ind2], Chan.hppvec[ind2], Chan.np[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	      Ints.S_ME1.V11[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell3];
	      ind = Index31(Chan.hpp1vec[ind2], Chan.pvec[ind2], Chan.nhpp1[ind2], Chan.np[ind2], shell1, shell2, shell4, shell3);
	      Ints.S_ME1.V13[ind2][ind] = TBME;
	      ind2 = Chan.indvec[shell4];
	      ind = Index31(Chan.hpp1vec[ind2], Chan.pvec[ind2], Chan.nhpp1[ind2], Chan.np[ind2], shell1, shell2, shell3, shell4);
	      Ints.S_ME1.V13[ind2][ind] = -1.0 * TBME;
	      
	      minus(tb, Space.qnums[shell1], Space.qnums[shell3]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.pp1vec[ind2], Chan.hp1vec[ind2], Chan.npp1[ind2], Chan.nhp1[ind2], shell4, shell2, shell1, shell3);
	      Ints.S_ME1.V16[ind2][ind] = TBME;
	      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.pp1vec[ind2], Chan.hp1vec[ind2], Chan.npp1[ind2], Chan.nhp1[ind2], shell3, shell2, shell1, shell4);
	      Ints.S_ME1.V16[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell2];
	      ind = Index31(Chan.hppvec[ind2], Chan.pvec[ind2], Chan.nhpp[ind2], Chan.np[ind2], shell1, shell3, shell4, shell2);
	      Ints.S_ME1.V17[ind2][ind] = TBME;
	      ind = Index31(Chan.hppvec[ind2], Chan.pvec[ind2], Chan.nhpp[ind2], Chan.np[ind2], shell1, shell4, shell3, shell2);
	      Ints.S_ME1.V17[ind2][ind] = -1.0 * TBME;
	      
	      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      ind = Index22(Chan.ppvec[ind2], Chan.hpvec[ind2], Chan.npp[ind2], Chan.nhp[ind2], shell3, shell4, shell1, shell2);
	      Ints.S_ME1.V20[ind2][ind] = TBME;
	      ind = Index22(Chan.ppvec[ind2], Chan.hpvec[ind2], Chan.npp[ind2], Chan.nhp[ind2], shell4, shell3, shell1, shell2);
	      Ints.S_ME1.V20[ind2][ind] = -1.0 * TBME;
	    }
	    else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
	      ind2 = Chan.indvec[shell3];
	      ind = Index13(Chan.hvec[ind2], Chan.hhpvec[ind2], Chan.nh[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	      Ints.S_ME1.V12[ind2][ind] = TBME;
	      ind = Index13(Chan.hvec[ind2], Chan.hhpvec[ind2], Chan.nh[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	      Ints.S_ME1.V12[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell1];
	      ind = Index31(Chan.hhp1vec[ind2], Chan.hvec[ind2], Chan.nhhp1[ind2], Chan.nh[ind2], shell2, shell3, shell4, shell1);
	      Ints.S_ME1.V14[ind2][ind] = TBME;
	      ind2 = Chan.indvec[shell2];
	      ind = Index31(Chan.hhp1vec[ind2], Chan.hvec[ind2], Chan.nhhp1[ind2], Chan.nh[ind2], shell1, shell3, shell4, shell2);
	      Ints.S_ME1.V14[ind2][ind] = -1.0 * TBME;
	      
	      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.hh1vec[ind2], Chan.hp1vec[ind2], Chan.nhh1[ind2], Chan.nhp1[ind2], shell3, shell2, shell1, shell4);
	      Ints.S_ME1.V15[ind2][ind] = TBME;
	      minus(tb, Space.qnums[shell2], Space.qnums[shell4]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.hh1vec[ind2], Chan.hp1vec[ind2], Chan.nhh1[ind2], Chan.nhp1[ind2], shell3, shell1, shell2, shell4);
	      Ints.S_ME1.V15[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell3];
	      ind = Index31(Chan.hhpvec[ind2], Chan.hvec[ind2], Chan.nhhp[ind2], Chan.nh[ind2], shell1, shell2, shell4, shell3);
	      Ints.S_ME1.V18[ind2][ind] = TBME;
	      ind = Index31(Chan.hhpvec[ind2], Chan.hvec[ind2], Chan.nhhp[ind2], Chan.nh[ind2], shell2, shell1, shell4, shell3);
	      Ints.S_ME1.V18[ind2][ind] = -1.0 * TBME;
	      
	      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      ind = Index22(Chan.hpvec[ind2], Chan.hhvec[ind2], Chan.nhp[ind2], Chan.nhh[ind2], shell3, shell4, shell1, shell2);
	      Ints.S_ME1.V19[ind2][ind] = TBME;
	      ind = Index22(Chan.hpvec[ind2], Chan.hhvec[ind2], Chan.nhp[ind2], Chan.nhh[ind2], shell3, shell4, shell2, shell1);
	      Ints.S_ME1.V19[ind2][ind] = -1.0 * TBME;
	    }
	  }
	}
      }
    }
  }
}

void Get_Matrix_Elements_J(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind, ind1, ind2;
  std::string ptype, qtype, rtype, stype;
  State tb;

  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    #pragma omp parallel private(shell1, shell2, shell3, shell4, ptype, qtype, rtype, stype, ind, ind1, ind2, tb, TBME)
    {
      #pragma omp for schedule(static)
      for(int tb1 = 0; tb1 < HF_Chan.ntb[chan]; ++tb1){
	shell1 = HF_Chan.tbvec[chan][2*tb1];
	shell2 = HF_Chan.tbvec[chan][2*tb1 + 1];
	ptype = Space.qnums[shell1].type;
	qtype = Space.qnums[shell2].type;
	if(shell1 >= shell2){ continue; }
	for(int tb2 = 0; tb2 < HF_Chan.ntb[chan]; ++tb2){
	  shell3 = HF_Chan.tbvec[chan][2*tb2];
	  shell4 = HF_Chan.tbvec[chan][2*tb2 + 1];
	  rtype = Space.qnums[shell3].type;
	  stype = Space.qnums[shell4].type;
	  if(shell3 >= shell4){ continue; }
	  if(shell1 == shell3 && shell2 > shell4){ continue; }
	  
	  if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "hole"){ continue; }
	  if(ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "particle"){ continue; }
	  
	  TBME = HF_ME.V[chan][tb1*HF_Chan.ntb[chan] + tb2];
	  if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V2[ind1][ind] = TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V2[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V2[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V2[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell3, shell4, shell1, shell2);
	    Ints.D_ME1.V2[ind1][ind] = TBME;
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell4, shell3, shell1, shell2);
	    Ints.D_ME1.V2[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell3, shell4, shell2, shell1);
	    Ints.D_ME1.V2[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell4, shell3, shell2, shell1);
	    Ints.D_ME1.V2[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	  }
	  else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V1[ind1][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V1[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V1[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V1[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell3, shell4, shell1, shell2);
	    Ints.D_ME1.V1[ind1][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell4, shell3, shell1, shell2);
	    Ints.D_ME1.V1[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell3, shell4, shell2, shell1);
	    Ints.D_ME1.V1[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell4, shell3, shell2, shell1);
	    Ints.D_ME1.V1[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	  }
	  else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp2vec[ind1], Chan.nhp2[ind1], Chan.nhp2[ind1], shell1, shell4, shell3, shell2);
	    Ints.D_ME1.V3[ind1][ind] = TBME;
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp2vec[ind1], Chan.nhp2[ind1], Chan.nhp2[ind1], shell3, shell2, shell1, shell4);
	    Ints.D_ME1.V3[ind1][ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell3, shell4, shell1, shell2);
	    Ints.D_ME1.V4[ind1][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell4, shell3, shell1, shell2);
	    Ints.D_ME1.V4[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell3, shell4, shell2, shell1);
	    Ints.D_ME1.V4[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell4, shell3, shell2, shell1);
	    Ints.D_ME1.V4[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    
	    ind2 = Chan.indvec[shell2];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V5[ind2][ind] = TBME;
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V5[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind2 = Chan.indvec[shell1];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V5[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V5[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    
	    ind2 = Chan.indvec[shell1];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell3, shell4);
	    Ints.D_ME1.V6[ind2][ind] = TBME;
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell4, shell3);
	    Ints.D_ME1.V6[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind2 = Chan.indvec[shell2];
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	    Ints.D_ME1.V6[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	    Ints.D_ME1.V6[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    
	    ind2 = Chan.indvec[shell4];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell1, shell2, shell3);
	    Ints.D_ME1.V7[ind2][ind] = TBME;
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell2, shell1, shell3);
	    Ints.D_ME1.V7[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind2 = Chan.indvec[shell3];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	    Ints.D_ME1.V7[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	    Ints.D_ME1.V7[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    
	    ind2 = Chan.indvec[shell3];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	    Ints.D_ME1.V8[ind2][ind] = TBME;
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	    Ints.D_ME1.V8[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    ind2 = Chan.indvec[shell4];
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell1, shell2, shell3);
	    Ints.D_ME1.V8[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell2, shell1, shell3);
	    Ints.D_ME1.V8[ind2][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    
	    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell3, shell2, shell4);
	    Ints.D_ME1.V9[ind1][ind] = TBME;
	    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell3, shell1, shell4);
	    Ints.D_ME1.V9[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell4, shell2, shell3);
	    Ints.D_ME1.V9[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell4, shell1, shell3);
	    Ints.D_ME1.V9[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	    
	    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell4, shell2, shell3);
	    Ints.D_ME1.V10[ind1][ind] = TBME;
	    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell4, shell1, shell3);
	    Ints.D_ME1.V10[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - Chan.qnums1[chan].j);
	    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell3, shell2, shell4);
	    Ints.D_ME1.V10[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - Chan.qnums1[chan].j);
	    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell3, shell1, shell4);
	    Ints.D_ME1.V10[ind1][ind] = TBME * pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j);
	  }
	  /*if(Parameters.approx == "singles"){
	    if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	      ind2 = Chan.indvec[shell2];
	      ind = Index13(Chan.pvec[ind2], Chan.hppvec[ind2], Chan.np[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	      Ints.S_ME1.V11[ind2][ind] = TBME;
	      ind = Index13(Chan.pvec[ind2], Chan.hppvec[ind2], Chan.np[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	      Ints.S_ME1.V11[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell3];
	      ind = Index31(Chan.hpp1vec[ind2], Chan.pvec[ind2], Chan.nhpp1[ind2], Chan.np[ind2], shell1, shell2, shell4, shell3);
	      Ints.S_ME1.V13[ind2][ind] = TBME;
	      ind2 = Chan.indvec[shell4];
	      ind = Index31(Chan.hpp1vec[ind2], Chan.pvec[ind2], Chan.nhpp1[ind2], Chan.np[ind2], shell1, shell2, shell3, shell4);
	      Ints.S_ME1.V13[ind2][ind] = -1.0 * TBME;
	      
	      minus(tb, Space.qnums[shell1], Space.qnums[shell3]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.pp1vec[ind2], Chan.hp1vec[ind2], Chan.npp1[ind2], Chan.nhp1[ind2], shell4, shell2, shell1, shell3);
	      Ints.S_ME1.V16[ind2][ind] = TBME;
	      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.pp1vec[ind2], Chan.hp1vec[ind2], Chan.npp1[ind2], Chan.nhp1[ind2], shell3, shell2, shell1, shell4);
	      Ints.S_ME1.V16[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell2];
	      ind = Index31(Chan.hppvec[ind2], Chan.pvec[ind2], Chan.nhpp[ind2], Chan.np[ind2], shell1, shell3, shell4, shell2);
	      Ints.S_ME1.V17[ind2][ind] = TBME;
	      ind = Index31(Chan.hppvec[ind2], Chan.pvec[ind2], Chan.nhpp[ind2], Chan.np[ind2], shell1, shell4, shell3, shell2);
	      Ints.S_ME1.V17[ind2][ind] = -1.0 * TBME;
	      
	      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      ind = Index22(Chan.ppvec[ind2], Chan.hpvec[ind2], Chan.npp[ind2], Chan.nhp[ind2], shell3, shell4, shell1, shell2);
	      Ints.S_ME1.V20[ind2][ind] = TBME;
	      ind = Index22(Chan.ppvec[ind2], Chan.hpvec[ind2], Chan.npp[ind2], Chan.nhp[ind2], shell4, shell3, shell1, shell2);
	      Ints.S_ME1.V20[ind2][ind] = -1.0 * TBME;
	    }
	    else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
	      ind2 = Chan.indvec[shell3];
	      ind = Index13(Chan.hvec[ind2], Chan.hhpvec[ind2], Chan.nh[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	      Ints.S_ME1.V12[ind2][ind] = TBME;
	      ind = Index13(Chan.hvec[ind2], Chan.hhpvec[ind2], Chan.nh[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	      Ints.S_ME1.V12[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell1];
	      ind = Index31(Chan.hhp1vec[ind2], Chan.hvec[ind2], Chan.nhhp1[ind2], Chan.nh[ind2], shell2, shell3, shell4, shell1);
	      Ints.S_ME1.V14[ind2][ind] = TBME;
	      ind2 = Chan.indvec[shell2];
	      ind = Index31(Chan.hhp1vec[ind2], Chan.hvec[ind2], Chan.nhhp1[ind2], Chan.nh[ind2], shell1, shell3, shell4, shell2);
	      Ints.S_ME1.V14[ind2][ind] = -1.0 * TBME;
	      
	      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.hh1vec[ind2], Chan.hp1vec[ind2], Chan.nhh1[ind2], Chan.nhp1[ind2], shell3, shell2, shell1, shell4);
	      Ints.S_ME1.V15[ind2][ind] = TBME;
	      minus(tb, Space.qnums[shell2], Space.qnums[shell4]);
	      ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      ind = Index22(Chan.hh1vec[ind2], Chan.hp1vec[ind2], Chan.nhh1[ind2], Chan.nhp1[ind2], shell3, shell1, shell2, shell4);
	      Ints.S_ME1.V15[ind2][ind] = -1.0 * TBME;
	      
	      ind2 = Chan.indvec[shell3];
	      ind = Index31(Chan.hhpvec[ind2], Chan.hvec[ind2], Chan.nhhp[ind2], Chan.nh[ind2], shell1, shell2, shell4, shell3);
	      Ints.S_ME1.V18[ind2][ind] = TBME;
	      ind = Index31(Chan.hhpvec[ind2], Chan.hvec[ind2], Chan.nhhp[ind2], Chan.nh[ind2], shell2, shell1, shell4, shell3);
	      Ints.S_ME1.V18[ind2][ind] = -1.0 * TBME;
	      
	      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      ind = Index22(Chan.hpvec[ind2], Chan.hhvec[ind2], Chan.nhp[ind2], Chan.nhh[ind2], shell3, shell4, shell1, shell2);
	      Ints.S_ME1.V19[ind2][ind] = TBME;
	      ind = Index22(Chan.hpvec[ind2], Chan.hhvec[ind2], Chan.nhp[ind2], Chan.nhh[ind2], shell3, shell4, shell2, shell1);
	      Ints.S_ME1.V19[ind2][ind] = -1.0 * TBME;
	    }
	    }*/
	}
      }
    }
  }
}

/*void Get_Matrix_Elements_J(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME0, TBME, m1, m2, t1, t2, CGC1, CGC2; // interaction two-body interaction ME and two-body COM ME
  int pind, qind, rind, sind;
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind, ind1, ind2;
  std::string ptype, qtype, rtype, stype;
  State tb;

  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    for(int tb1 = 0; tb1 < HF_Chan.tb[chan]; ++tb1){
      shell1 = HF_Chan.tbvec1[chan][2*tb1];
      shell2 = HF_Chan.tbvec1[chan][2*tb1 + 1];
      ptype = Space.qnums[shell1].type;
      qtype = Space.qnums[shell2].type;
      if(ptype == qtype && shell1 > shell2){ continue; }
      if(ptype == "particle" && qtype == "hole"){ continue; }
      for(int tb2 = tb1; tb2 < HF_Chan.tb[chan]; ++tb2){
	shell3 = HF_Chan.tbvec1[chan][2*tb2];
	shell4 = HF_Chan.tbvec1[chan][2*tb2 + 1];
	rtype = Space.qnums[shell3].type;
	stype = Space.qnums[shell4].type;
	if(rtype == stype && shell3 > shell4){ continue; }
	if(rtype == "particle" && stype == "hole"){ continue; }
	if(ptype == rtype && qtype == stype && shell1 == shell3 && shell2 > shell4){ continue; }
	TBME0 = HF_ME.V[chan][tb1*HF_Chan.tb[chan] + tb2];
	//if(shell1 == shell2){ TBME0 *= sqrt(2.0); } // !!check
	//if(shell3 == shell4){ TBME0 *= sqrt(2.0); } // !!check
	for(int jz = -HF_Chan.qnums1[chan].j; jz <= HF_Chan.qnums1[chan].j; jz+=2){
	  for(int p = 0; p < int(Space.shellsm[shell1].size()); ++p){
	    for(int q = 0; q < int(Space.shellsm[shell2].size()); ++q){
	      for(int r = 0; r < int(Space.shellsm[shell3].size()); ++r){
		for(int s = 0; s < int(Space.shellsm[shell4].size()); ++s){
		  pind = Space.shellsm[shell1][p];
		  qind = Space.shellsm[shell2][q];
		  rind = Space.shellsm[shell3][r];
		  sind = Space.shellsm[shell4][s];
		  m1 = Space.qnums[pind].m + Space.qnums[qind].m;
		  m2 = Space.qnums[rind].m + Space.qnums[sind].m;
		  t1 = Space.qnums[pind].t + Space.qnums[qind].t;
		  t2 = Space.qnums[rind].t + Space.qnums[sind].t;
		  if(t1 != t2 || m1 != jz || m2 != jz){ continue; }
		  if(pind >= qind && shell1 == shell2){ continue; }
		  if(rind >= sind && shell3 == shell4){ continue; }
		  if(pind > rind && shell1 == shell3 && shell2 == shell4){ continue; }
		  if(pind == rind && qind > sind && shell1 == shell3 && shell2 == shell4){ continue; }
		  CGC1 = CGC(0.5*Space.qnums[pind].j, 0.5*Space.qnums[pind].m, 0.5*Space.qnums[qind].j, 0.5*Space.qnums[qind].m, 0.5*HF_Chan.qnums1[chan].j, double(0.5*jz));
		  CGC2 = CGC(0.5*Space.qnums[rind].j, 0.5*Space.qnums[rind].m, 0.5*Space.qnums[sind].j, 0.5*Space.qnums[sind].m, 0.5*HF_Chan.qnums1[chan].j, double(0.5*jz));
		  TBME = TBME0 * CGC1 * CGC2;
		  if(fabs(TBME) < 1e-12){ continue; }
		  ptype = Space.qnums[pind].type;
		  qtype = Space.qnums[qind].type;
		  rtype = Space.qnums[rind].type;
		  stype = Space.qnums[sind].type;
		  if(ptype == "particle" && qtype == "hole"){ std::swap(pind, qind); std::swap(ptype, qtype); TBME *= -1.0; }
		  if(rtype == "particle" && stype == "hole"){ std::swap(rind, sind); std::swap(rtype, stype); TBME *= -1.0; }
		  if((ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "hole") || 
		     (ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "hole") ||
		     (ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "particle")){
		    std::swap(pind, rind);
		    std::swap(qind, sind);
		    std::swap(ptype, rtype);
		    std::swap(qtype, stype);
		  }
		  if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
		    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		  }
		  else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2); 
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell1, shell2);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell1, shell2);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell2, shell1);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell2, shell1);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
		    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell1, shell4, shell3, shell2);
		    Ints.D_ME1.V3[ind1][ind] = TBME;
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell3, shell2, shell1, shell4);
		    Ints.D_ME1.V3[ind1][ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
		    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
		    Ints.D_ME1.V4[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
		    Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
		    Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
		    Ints.D_ME1.V4[ind1][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell3];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
		    Ints.D_ME1.V5[ind2][ind] = TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
		    Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell4];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
		    Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
		    Ints.D_ME1.V5[ind2][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell4];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
		    Ints.D_ME1.V6[ind2][ind] = TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
		    Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell3];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
		    Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
		    Ints.D_ME1.V6[ind2][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell1];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V7[ind2][ind] = TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell2];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V7[ind2][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell2];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V8[ind2][ind] = TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell1];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V8[ind2][ind] = TBME;
		    
		    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
		    Ints.D_ME1.V9[ind1][ind] = TBME;
		    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
		    Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
		    Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
		    Ints.D_ME1.V9[ind1][ind] = TBME;
		    
		    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
		    Ints.D_ME1.V10[ind1][ind] = TBME;
		    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
		    Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
		    Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
		    Ints.D_ME1.V10[ind1][ind] = TBME;
		  }
		  if(Parameters.approx == "singles"){
		    if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		      ind2 = Chan.indvec[shell2];
		      ind = Index13(Chan.pvec1[ind2], Chan.hppvec1[ind2], Chan.p[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
		      Ints.S_ME1.HPPP[ind2][ind] = TBME;
		      ind = Index13(Chan.pvec1[ind2], Chan.hppvec1[ind2], Chan.p[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
		      Ints.S_ME1.HPPP[ind2][ind] = -1.0 * TBME;
		      
		      ind = Index31(Chan.hppvec1[ind2], Chan.pvec1[ind2], Chan.hpp[ind2], Chan.p[ind2], shell1, shell3, shell4, shell2);
		      Ints.S_ME1.HPPP2[ind2][ind] = TBME;
		      ind = Index31(Chan.hppvec1[ind2], Chan.pvec1[ind2], Chan.hpp[ind2], Chan.p[ind2], shell1, shell4, shell3, shell2);
		      Ints.S_ME1.HPPP2[ind2][ind] = -1.0 * TBME;
		      
		      ind2 = Chan.indvec[shell3];
		      ind = Index31(Chan.hpp2vec1[ind2], Chan.pvec1[ind2], Chan.hpp2[ind2], Chan.p[ind2], shell1, shell2, shell4, shell3);
		      Ints.S_ME1.V13[ind2][ind] = TBME;
		      ind2 = Chan.indvec[shell4];
		      ind = Index31(Chan.hpp2vec1[ind2], Chan.pvec1[ind2], Chan.hpp2[ind2], Chan.p[ind2], shell1, shell2, shell3, shell4);
		      Ints.S_ME1.V13[ind2][ind] = -1.0 * TBME;
		      
		      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		      ind = Index22(Chan.ppvec1[ind2], Chan.hpvec1[ind2], Chan.pp[ind2], Chan.hp[ind2], shell3, shell4, shell1, shell2);
		      Ints.S_ME1.HPPP4[ind2][ind] = TBME;
		      ind = Index22(Chan.ppvec1[ind2], Chan.hpvec1[ind2], Chan.pp[ind2], Chan.hp[ind2], shell4, shell3, shell1, shell2);
		      Ints.S_ME1.HPPP4[ind2][ind] = -1.0 * TBME;

		      minus(tb, Space.qnums[shell1], Space.qnums[shell3]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell1, shell3);
		      ind = Index22(Chan.pp1vec1[ind2], Chan.hp1vec1[ind2], Chan.pp1[ind2], Chan.hp1[ind2], shell4, shell2, shell1, shell3);
		      Ints.S_ME1.HPPP5[ind2][ind] = TBME;
		      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell1, shell4);
		      ind = Index22(Chan.pp1vec1[ind2], Chan.hp1vec1[ind2], Chan.pp1[ind2], Chan.hp1[ind2], shell3, shell2, shell1, shell4);
		      Ints.S_ME1.HPPP5[ind2][ind] = -1.0 * TBME;
		    }
		    else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
		      ind2 = Chan.indvec[shell3];
		      ind = Index13(Chan.hvec1[ind2], Chan.hhpvec1[ind2], Chan.h[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
		      Ints.S_ME1.HHHP[ind2][ind] = TBME;
		      ind = Index13(Chan.hvec1[ind2], Chan.hhpvec1[ind2], Chan.h[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
		      Ints.S_ME1.HHHP[ind2][ind] = -1.0 * TBME;
		      
		      ind = Index31(Chan.hhpvec1[ind2], Chan.hvec1[ind2], Chan.hhp[ind2], Chan.h[ind2], shell1, shell2, shell4, shell3);
		      Ints.S_ME1.HHHP2[ind2][ind] = TBME;
		      ind = Index31(Chan.hhpvec1[ind2], Chan.hvec1[ind2], Chan.hhp[ind2], Chan.h[ind2], shell2, shell1, shell4, shell3);
		      Ints.S_ME1.HHHP2[ind2][ind] = -1.0 * TBME;
		      
		      ind2 = Chan.indvec[shell1];
		      ind = Index31(Chan.hhp2vec1[ind2], Chan.hvec1[ind2], Chan.hhp2[ind2], Chan.h[ind2], shell2, shell3, shell4, shell1);
		      Ints.S_ME1.V14[ind2][ind] = TBME;
		      ind2 = Chan.indvec[shell2];
		      ind = Index31(Chan.hhp2vec1[ind2], Chan.hvec1[ind2], Chan.hhp2[ind2], Chan.h[ind2], shell1, shell3, shell4, shell2);
		      Ints.S_ME1.V14[ind2][ind] = -1.0 * TBME;

		      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);		      
		      //ind2 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		      ind = Index22(Chan.hpvec1[ind2], Chan.hhvec1[ind2], Chan.hp[ind2], Chan.hh[ind2], shell3, shell4, shell1, shell2);
		      Ints.S_ME1.HHHP4[ind2][ind] = TBME;
		      ind = Index22(Chan.hpvec1[ind2], Chan.hhvec1[ind2], Chan.hp[ind2], Chan.hh[ind2], shell3, shell4, shell2, shell1);
		      Ints.S_ME1.HHHP4[ind2][ind] = -1.0 * TBME;
		      
		      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell1, shell4);
		      ind = Index22(Chan.hh1vec1[ind2], Chan.hp1vec1[ind2], Chan.hh1[ind2], Chan.hp1[ind2], shell3, shell2, shell1, shell4);
		      Ints.S_ME1.HHHP5[ind2][ind] = TBME;
		      minus(tb, Space.qnums[shell2], Space.qnums[shell4]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell2, shell4);
		      ind = Index22(Chan.hh1vec1[ind2], Chan.hp1vec1[ind2], Chan.hh1[ind2], Chan.hp1[ind2], shell3, shell1, shell2, shell4);
		      Ints.S_ME1.HHHP5[ind2][ind] = -1.0 * TBME;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}*/
