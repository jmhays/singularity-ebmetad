#include <vector>
#include <valarray> 
#include <iostream> 
#include <string>
#include <sstream>
#include <fstream>

class residueAlias{
private:
  std::string fake, real;
public:
  void set( const std::string& f, const std::string& r ) { fake=f; real=r; }
  int fakeName( const std::string& nam ) const { return fake.compare( nam ); }
  std::string getRealName() const { return real; }
};

class resLibObj;

class atomObj{
friend std::istream& operator>>(std::istream& , atomObj&);
friend std::ostream& operator<<(std::ostream& , const atomObj& );
friend int sameResidue( const atomObj& , const atomObj& );
friend int inLibrary( const atomObj& , const std::vector<resLibObj>& );
friend int libNumber( const atomObj& , const std::vector<resLibObj>& );
private:
  int resnum, atnum;
  std::string resname, atomname;
public:

  int AtomNameIs( const std::string& nam ) const { return atomname.compare( nam ); }
  int AtomNumber() const { return atnum; }
  int ResNo() const { return resnum; }
  std::string getName() const { return resname; }
  void checkAlias( const residueAlias& alias ){ if( alias.fakeName( resname )==0 ){ resname=alias.getRealName(); } }
};

class torsionObj{
friend std::ostream& operator<<(std::ostream&  , const torsionObj&  );
private:
  std::string at1,at2,at3,at4;
public:
  void set( const std::string& a1, const std::string& a2, const std::string& a3, const std::string& a4 ){ at1=a1; at2=a2; at3=a3; at4=a4; } 
  void findAtoms( const std::vector<atomObj>& residue, const std::vector<atomObj>& mresidue, const std::vector<atomObj>& presidue, int& a1, int& a2, int& a3, int& a4 ) const ;
};

class resLibObj{
private:
  std::string name;
  unsigned long ntorsions;
  std::vector< torsionObj > torsions;
public:
  void set( const std::string& n, const int& num ){ name=n; ntorsions=num; torsions.resize( ntorsions ); }
  void setTorsion( const int& n, const std::string& a1, const std::string& a2, const std::string& a3, const std::string& a4 ){
     torsions[n].set( a1, a2, a3, a4 );
  }
  int isRes( const std::string& nam ) const { if( name.compare( nam )==0 ){ return 1; } return 0; }
  void createTorsions( const int& resno, const int& nres, const std::vector<atomObj>& residue, const std::vector<atomObj>& mresidue, const std::vector<atomObj>& presidue );
};

// This compares the residue numbers of two amino acids
int sameResidue( const atomObj& a1, const atomObj& a2){ if( a1.resnum==a2.resnum ){ return 1; }  return 0; }
// This splits out each residue in the atom list
void getResidue( const int& resno, const std::vector<atomObj>& atoms, std::vector<atomObj>& residue );

void readGromacsInput( char* ifilename, std::vector<atomObj>& atoms );
int sortResidues( std::vector<atomObj>& atoms, const std::vector<resLibObj>& residue_library, const std::vector<residueAlias>& aliases );
