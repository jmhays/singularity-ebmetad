void setup_residue_library( const int& lev, std::vector<resLibObj>& residue_library, std::vector<residueAlias>& aliases ) {
  if( lev!=0 && lev!=1 && lev!=2 ){ std::cerr<<"!ERRROR - not sure what you are trying to do"<<std::endl; abort(); }
  // Adjust these when you add extra aliases or resdiues
  unsigned long nresidues=22, nalias=5; 

  // Resize the vectors that we store everything in
  residue_library.resize( nresidues ); aliases.resize( nalias );

  // Setup the aliases - amino acids that are named differently from their three letter codes
  aliases[0].set( "NASN", "ASN" );
  aliases[1].set( "CSER", "SER" );
  aliases[2].set( "LYP", "LYS" );  
  aliases[3].set( "LYSH", "LYS" );
  aliases[4].set( "CYS2", "CYS" );

  // ALANINE  -- ONLY BACKBONE ANGLES
  residue_library[0].set( "ALA",2 ); 
  residue_library[0].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[0].setTorsion( 1, "-C", "N", "CA", "C" );
  
  // GLYCINE -- ONLY BACKBONE ANGLES 
  residue_library[1].set( "GLY",2 );
  residue_library[1].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[1].setTorsion( 1, "-C", "N", "CA", "C" );
  
  // LEUCINE -- Side chain -CH(CH3)CH2CH3
  if( lev==0 ){ residue_library[2].set( "LEU",2 ); }
  else if( lev==1 ){ residue_library[2].set( "LEU",3 ); }
  else if( lev==2 ){ residue_library[2].set( "LEU",4 ); }
  residue_library[2].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[2].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[2].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){ residue_library[2].setTorsion( 3, "CA", "CB", "CG", "CD1" ); }
 
  // ARGENINE  -- Side chain -CH2CH2CH2NHC(NH2)2
  if( lev==0 ){ residue_library[3].set( "ARG",2 ); }
  else if( lev==1 ){ residue_library[3].set( "ARG",3 ); }
  else if( lev==2 ){ residue_library[3].set( "ARG",6 ); }
  residue_library[3].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[3].setTorsion( 1, "-C", "N", "CA", "C" ); 
  if(lev>0){ residue_library[3].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){
    residue_library[3].setTorsion( 3, "CA", "CB", "CG", "CD" );
    residue_library[3].setTorsion( 4, "CB", "CG", "CD", "NE" );
    residue_library[3].setTorsion( 5, "CG", "CD", "NE", "CZ" ); 
  }
 
  // ASPARTIC ACID -- Side chain -CH(NH2)CH2COO-
  if( lev==0 ){ residue_library[4].set( "ASP",2 ); }
  else if( lev==1 ){ residue_library[4].set( "ASP",3 ); }
  else if( lev==2 ){ residue_library[4].set( "ASP",4 ); }
  residue_library[4].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[4].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[4].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){ residue_library[4].setTorsion( 3, "CA", "CB", "CG", "OD1" ); }
  
  // SERINE -- Side chain -CH2OH
  if( lev==0 ){ residue_library[5].set( "SER",2 ); }
  else if( lev==1 || lev==2 ){ residue_library[5].set( "SER",3 ); }
  residue_library[5].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[5].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[5].setTorsion( 2, "N", "CA", "CB", "OG" ); }
  
  // PROLINE  -- Cyclic amino acid - ONLY BACKBONE ANGLES 
  residue_library[6].set( "PRO",2 );
  residue_library[6].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[6].setTorsion( 1, "-C", "N", "CA", "C" );
  
  // GLUTAMINE -- Side chain -CH2CH2CONH2
  if( lev==0 ){ residue_library[7].set( "GLN",2 ); }
  else if( lev==1 ){ residue_library[7].set( "GLN",3 ); }
  else if( lev==2 ){ residue_library[7].set( "GLN",5 ); }
  residue_library[7].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[7].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[7].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){
    residue_library[7].setTorsion( 3, "CA", "CB", "CG", "CD" );
    residue_library[7].setTorsion( 4, "CB", "CG", "CD", "OE1"  );
  }
  
  // TRYPTOPHAN -- Side chain -CH2(RingThing)
  if( lev==0 ){ residue_library[8].set( "TRP",2 ); }
  else if( lev==1 ){ residue_library[8].set( "TRP",3 ); }
  else if( lev==2 ){ residue_library[8].set( "TRP",4 ); }
  residue_library[8].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[8].setTorsion( 1, "-C", "N", "CA", "C" ); 
  if(lev>0){residue_library[8].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){residue_library[8].setTorsion( 3, "CA", "CB", "CG", "CD1" ); }
  
  // TYROSINE -- Side chain -CH2Ph-OH
  if( lev==0 ){ residue_library[9].set( "TYR",2 ); }
  else if( lev==1 ){ residue_library[9].set( "TYR",3 ); }
  else if( lev==2 ){ residue_library[9].set( "TYR",4 ); }
  residue_library[9].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[9].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){residue_library[9].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){residue_library[9].setTorsion( 3, "CA", "CB", "CG", "CD1" ); }  

  // ISOLEUCINE -- Side chain -CH(CH3)CH2CH3
  if(lev==0 ){ residue_library[10].set( "ILE",2 ); }
  else if(lev==1 ){ residue_library[10].set( "ILE",3 ); }
  else if(lev==2 ){ residue_library[10].set( "ILE",4 ); }
  residue_library[10].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[10].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[10].setTorsion( 2, "N", "CA", "CB", "CG1" ); }
  if(lev>1){ residue_library[10].setTorsion( 3, "CA", "CB", "CG1", "CD" ); }  

  //  LYSINE -- Side chain -CH2CH2CH2CH2NH3+
  if( lev==0 ){ residue_library[11].set( "LYS",2 ); }
  else if( lev==1 ){ residue_library[11].set( "LYS",3 ); }
  else if( lev==2 ){ residue_library[11].set( "LYS",6 ); }
  residue_library[11].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[11].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[11].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){
    residue_library[11].setTorsion( 3, "CA", "CB", "CG", "CD" );
    residue_library[11].setTorsion( 4,"CB", "CG", "CD", "CE" );
    residue_library[11].setTorsion( 5,"CG", "CD", "CE", "NZ"  );
  }

  // ASPARAGINE 
  if( lev==0 ){ residue_library[12].set( "ASN",2 ); }
  else if( lev==1 ){ residue_library[12].set( "ASN",3 ); }
  else if( lev==2 ){ residue_library[12].set( "ASN",4 ); } 
  residue_library[12].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[12].setTorsion( 1, "-C", "N", "CA", "C" ); 
  if(lev>0){ residue_library[12].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){ residue_library[12].setTorsion( 3, "CA", "CB", "CG", "ND2" ); }

  // THREONINE
  if( lev==0 ){ residue_library[13].set( "THR",2 ); }
  else if( lev==1 || lev==2 ){ residue_library[13].set( "THR",3 ); }
  residue_library[13].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[13].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[13].setTorsion( 2, "N", "CA", "CB", "CG2" ); }

  // HISTIDINE 
  if(lev==0){ residue_library[14].set( "HIS",2 ); }
  else if(lev==1){ residue_library[14].set( "HIS",3 ); }
  else if(lev==2){ residue_library[14].set( "HIS",4 ); }
  residue_library[14].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[14].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[14].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){ residue_library[14].setTorsion( 3, "CA", "CB", "CG", "CD2" ); }
  
  // GLUTAMIC ACID
  if(lev==0){ residue_library[15].set( "GLU",2 ); }
  else if(lev==1){ residue_library[15].set( "GLU",3 ); }
  else if(lev==2){ residue_library[15].set( "GLU",5 ); }
  residue_library[15].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[15].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[15].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){
     residue_library[15].setTorsion( 3, "CA", "CB", "CG", "CD");
     residue_library[15].setTorsion( 4, "CB", "CG", "CD", "OE1");
  } 

  // Phenylalanine
  if(lev==0){ residue_library[16].set( "PHE",2 ); }
  else if(lev==1){ residue_library[16].set( "PHE",3 ); }
  else if(lev==2){ residue_library[16].set( "PHE",4 ); }
  residue_library[16].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[16].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[16].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){ residue_library[16].setTorsion( 3, "CA", "CB", "CG", "CD1" ); }

  // Valine
  if(lev==0){ residue_library[17].set( "VAL",2 ); }
  else if(lev==1){ residue_library[17].set( "VAL",3 ); }
  else if(lev==2){ residue_library[17].set( "VAL",3 ); }
  residue_library[17].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[17].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[17].setTorsion( 2, "N", "CA", "CB", "CG1" ); }

  // Cysteine
  if(lev==0){ residue_library[18].set( "CYS",2 ); }
  else if(lev==1){ residue_library[18].set( "CYS",3 ); }
  else if(lev==2){ residue_library[18].set( "CYS",3 ); }
  residue_library[18].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[18].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[18].setTorsion( 2, "N", "CA", "CB", "SG" ); }

  // Methionine
  if(lev==0){ residue_library[19].set( "MET",2 ); }
  else if(lev==1){ residue_library[19].set( "MET",3 ); }
  else if(lev==2){ residue_library[19].set( "MET",5 ); }
  residue_library[19].setTorsion( 0, "N", "CA", "C", "+N" );
  residue_library[19].setTorsion( 1, "-C", "N", "CA", "C" );
  if(lev>0){ residue_library[19].setTorsion( 2, "N", "CA", "CB", "CG" ); }
  if(lev>1){
     residue_library[19].setTorsion( 3, "CA", "CB", "CG", "SD" );
     residue_library[19].setTorsion( 4, "CB", "CG", "SD", "CE" );  
  }

  // ACETYLENE TERMINAL GROUP - no torsions
  residue_library[20].set( "ACE",0 ); 

  // NAC TERMINAL GROUP - no torsions
  residue_library[21].set( "NAC",0 ); 

}
