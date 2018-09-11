#include "gen_torsions.h"
#include "residue_library.h"

extern "C" {
void gen_torsions_( int* lev, int* ftype, char* ifilename ){
   int level; level=(*lev);
   // Setup the residue library using the header file
   std::vector<resLibObj> residue_library; std::vector<residueAlias> aliases;
   setup_residue_library( level, residue_library, aliases );

   // Read in the gromacs input
   std::vector<atomObj> atoms;
   if( (*ftype)==1 ){ readGromacsInput( ifilename, atoms ); }
   else{ std::cerr<<"ERROR - filetype not implemented"<<std::endl; abort(); }

   // Change the names of any pesky alised residues and get the number of residues
   int nres; nres=sortResidues( atoms, residue_library, aliases );
   std::cerr<<"Found "<<nres<<" residues in protein"<<std::endl;

   // Split the atoms into their constituent residues
   std::vector< std::vector<atomObj> > residues(nres); std::vector<atomObj> dummy;
   for(unsigned long i=0;i<nres;++i){ 
     getResidue( i+1, atoms, residues[i] ); 
     std::cerr<<"RESIDUE "<<i+1<<" "<<residues[i][0].getName()<<" "<<libNumber( residues[i][0], residue_library )<<std::endl;
     //for(unsigned long j=0;j<residues[i].size();++j){std::cerr<<residues[i][j]<<std::endl; }
     //std::cerr<<std::endl;
  }

   // Create torsions for the N-terminus
   int resID; resID=libNumber( residues[0][0], residue_library );
   residue_library[resID].createTorsions( 1, nres, residues[0], dummy, residues[1] );

   // Create all other torsions
   for(int i=1;i<nres-1;++i){
      // Find out what residue this is
      resID=libNumber( residues[i][0], residue_library );
      // Create the torsions
      residue_library[resID].createTorsions( i+1, nres, residues[i], residues[i-1], residues[i+1] ); 
   }

   // Create torsions for the C-terminus
   resID=libNumber( residues[nres-1][0], residue_library );
   residue_library[resID].createTorsions( nres, nres, residues[nres-1], residues[nres-2], dummy );

   return ;
}
}

// Input and output for the atom object
std::istream& operator>>(std::istream& istr , atomObj& a ){
  istr>>a.resnum>>a.resname>>a.atomname>>a.atnum;
};

std::ostream& operator<<(std::ostream& ostr , const atomObj& a ){
  ostr<<a.resnum<<" "<<a.resname<<" "<<a.atomname<<" "<<a.atnum;
};

std::ostream& operator<<(std::ostream& ostr , const torsionObj& t ){
  ostr<<t.at1<<" "<<t.at2<<" "<<t.at3<<" "<<t.at4;
};

int inLibrary( const atomObj& a , const std::vector<resLibObj>& residue_library ){
   for(unsigned long i=0;i<residue_library.size();++i){
      if( residue_library[i].isRes( a.resname )==1 ){ return 1; }
   }
   return 0; 
}
int libNumber( const atomObj& a , const std::vector<resLibObj>& residue_library ){
   for(unsigned long i=0;i<residue_library.size();++i){
      if( residue_library[i].isRes( a.resname )==1 ){ return i; }
   }
   std::cerr<<"! ERROR - RESIDUE APPEARS NOT TO BE IN LIBRARY"<<std::endl; abort();
   return -1;
}

void resLibObj::createTorsions( const int& resno, const int& nres, const std::vector<atomObj>& residue, const std::vector<atomObj>& mresidue, const std::vector<atomObj>& presidue ){

  int a1, a2, a3, a4;
  for(unsigned long i=0;i<ntorsions;++i){
     if( resno==1 && i==1 ){ continue; }                // No Psi for residue at N-terminus
     if( resno==nres && i==0 ){ continue; }             // No Phi for residue at C-terminus
     torsions[i].findAtoms( residue, mresidue, presidue, a1, a2, a3, a4 );
     std::cout<<"TORSION LIST "<<a1<<" "<<a2<<" "<<a3<<" "<<a4<<"              # "<<name<<" atoms "<<torsions[i]<<std::endl;
  }

  return;
}

void torsionObj::findAtoms( const std::vector<atomObj>& residue, const std::vector<atomObj>& mresidue, const std::vector<atomObj>& presidue, int& a1, int& a2, int& a3, int& a4 ) const {

  a1=a2=a3=a4=-1;

  // Find the first atom in the torsion
  if( at1=="-C" ){ for(int i=0;i<mresidue.size();++i){ if(mresidue[i].AtomNameIs("C")==0){ a1=mresidue[i].AtomNumber(); break; } } }
  else{ for(int i=0;i<residue.size();++i){ if(residue[i].AtomNameIs(at1)==0){ a1=residue[i].AtomNumber(); break; } } }

  // Find the second atom in the torsion
  for(int i=0;i<residue.size();++i){ if(residue[i].AtomNameIs(at2)==0){ a2=residue[i].AtomNumber(); break; } }

  // Find the third atom in the torsion
  for(int i=0;i<residue.size();++i){ if(residue[i].AtomNameIs(at3)==0){ a3=residue[i].AtomNumber(); break; } }

  // Find the fourth atom in the torsion
  if( at4=="+N" ){ for(int i=0;i<presidue.size();++i){ if(presidue[i].AtomNameIs("N")==0){ a4=presidue[i].AtomNumber(); break; } } }
  else{ for(int i=0;i<residue.size();++i){ if(residue[i].AtomNameIs(at4)==0){ a4=residue[i].AtomNumber(); break; } } }
}

// This splits out each residue in the atom list
void getResidue( const int& resno, const std::vector<atomObj>& atoms, std::vector<atomObj>& residue ){
  // Count the number of atoms in this residue
  int natoms=0; for(unsigned long i=0;i<atoms.size();++i){ if( atoms[i].ResNo()==resno ){ natoms++; } }
  // Resize the residue
  residue.resize( natoms ); int n=0;
  // Put the atoms in this residue
  for(unsigned long i=0;i<atoms.size();++i){ if( atoms[i].ResNo()==resno ){ residue[n]=atoms[i]; n++; } }
  // And a little check
  if( n!=natoms ){ std::cerr<<"ERROR - something has gone wrong"<<std::endl; abort(); }
}

int sortResidues( std::vector<atomObj>& atoms, const std::vector<resLibObj>& residue_library, const std::vector<residueAlias>& aliases ){
   // First run through the atoms changing the names of any funky residues
   for(unsigned long i=0;i<atoms.size();++i){
      for(unsigned long j=0;j<aliases.size();++j){ atoms[i].checkAlias( aliases[j] ); }
      // std::cerr<<"ATOM CHECK "<<atoms[i]<<" "<<atoms[i].inLibrary( residue_library )<<std::endl;
   } 

   // This counts the number of residues
   int nres=1;
   for(unsigned long i=1;i<atoms.size();++i){
      if( sameResidue( atoms[i-1], atoms[i] )==0 ){
        if(i==1){std::cerr<<"ERROR - SERIOUSLY WEIRD PROTEIN INPUTTED"<<std::endl; abort(); }
        if( inLibrary( atoms[i], residue_library )==1 && inLibrary( atoms[i-1], residue_library )==1 ){ nres++; }
        else if( inLibrary( atoms[i], residue_library )==1 ){ std::cerr<<"ERROR - NON-CONTIGUOUS AMINO ACIDS IN PROTEIN"<<std::endl; abort(); }
      }
   }
   return nres;
}

void readGromacsInput( char* ifilename, std::vector<atomObj>& atoms ){
   // Open an ifstream to read in the file
   std::ifstream afile; afile.open(ifilename);

   // Read in the first line this is a comment
   std::string line; getline( afile, line );  

   // Read in the number of atoms
   unsigned long natoms; getline( afile, line ); std::istringstream inpt (line); inpt>>natoms;  
   atoms.resize( natoms ); std::cerr<<"Found "<<natoms<<" atoms"<<std::endl;
   
   // Read in the atoms
   for(unsigned long i=0;i<natoms;++i){ getline( afile, line ); std::istringstream atomIn (line); atomIn>>atoms[i]; }

   afile.close(); // Close the ifstream
}
