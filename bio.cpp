#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// (My Function)
// Converts every T charecter into U
string T_to_U (string str){
  string aaa = "";
  for( auto x : str){
    if( x == 'T'){
      aaa += 'U';
      continue;
    }
    aaa += x;
  }
  str = aaa;
  return str;
}

// (My Function)
// Finds the position of the each element of vector from another vector, and adds the position (integer) to another vector. Returns vector <int>
vector <int> Find_i( const vector<string> & v_1, const vector<string> & v_2){
  vector <int> res;
  int z = 0;

  for( int i = 0; i < v_1.size(); i++){
    cout << v_1[i] << " ";
  }
  cout << endl;
  cout << endl;
  for( int i = 0; i < v_2.size(); i++){
    cout << v_2[i] << " ";
  }
  cout << endl;


  for( int i = 0; i < res.size(); i++){
    cout << res[i] << " ";
  }
  cout << endl;

  return res;
}

/*
This function should return true if and only if
every character in the input is one of ATCG.
*/
bool IsValidDNASequence(const std::string & input){
  bool b = true;
  for( int i = 0; i < input.length(); ++i){
    if( input[i] != 'T' && input[i] != 'A' && input[i] != 'C' && input[i] != 'G' ){
      b = false;
      return b;
      break;
    }
  }
  return b;
}


/*
This function should calculate the reverse complement DNA sequence.

The first argument is the sequence, the second argument is a pointer to
an empty string, which you should modify to store the result.

This is obtained by reversing the input sequence and swaping each
nucleotide/letter with it's complement:
A <-> T
C <-> G

Example:
input = AAATTCGGGG
reverse = GGGGCTTAAA
reverse complement = CCCCGAATTT
*/
void GetReverseComplementSequence(const std::string & input,  std::string * const output){
  string str = "";
  for( int i = input.length() - 1; i >= 0; --i){
    if( input[i] == 'A' ){
      str += 'T';
    } else if( input[i] == 'T' ){
      str += 'A';
    } else if( input[i] == 'C' ){
      str += 'G';
    } else if( input[i] == 'G' ){
      str += 'C';
    }
  }
  *output = str; 
}


/*
This function should return the RNA transcript from a DNA sequence.

A RNA transcript is the reverse complement of the DNA sequence, but RNA
has U (uracil) instead of T (thiamine).

Make sure you don't have any redundant code.
*/
std::string GetRNATranscript(const std::string & input){
  string str = "";
  string* ptr = &str;
  GetReverseComplementSequence(input, ptr);
  for( int i = 0; i < str.length(); ++i){
    if( str[i] == 'T'){
      str[i] = 'U';
    }
  }
  return str;
}


/*
This function should return a vector of vector of strings with each possible RNA
reading frame from the given DNA sequence.

There are three possible reading frames (because the genetic code has three
nucleotides per amino acid) in each direction (you can also transcribe DNA in
the reverse complement direction, called the antiparallel strand).

Order the sequences like so:
1: Original (0 offset)
2: Original (1 offset)
3: Original (2 offset)
4: Antiparallel (0 offset)
5: Antiparallel (1 offset)
6: Antiparallel (2 offset)

With in the input sequence of: AATTCCCGAAA
Original RNA transcript = AAUUCCCGAAA
Antiparallel RNA transcript = UUUGCCCAAUU

The offsets (starting at pos 0, 1, and 2) of the two RNA transcripts
UUUGCCCAAUU
UUGCCCAAUU
UGCCCAAUU
AAUUCCCGAAA
AUUCCCGAAA
UUCCCGAAA

Instead of returning a vector of 6 strings, break each string into a vector
of length 3 strings (called codons) These codons will be useful
for the next translation step.

UUUGCCCAAUU -> {"UUU", "GCC", "CAA"}
// drop any remaining letters that don't fill a codon
UUGCCCAAUU -> {"UUG", "CCC", "AAU"}
UGCCCAAUU -> {"UGC", "CCA", "AUU"}
AAUUCCCGAAA -> {"AAU", "UCC", "CGA"}
AUUCCCGAAA -> {"AUU", "CCC", "GAA"}
UUCCCGAAA -> {"UUC", "CCG", "AAA"}
*/
std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string & input){
  vector< vector<string> > v_1(6);
  string o_rna = "", a_rna = "";

  string line = "";     //For saving the three charecters and inputting to the vector after
  int l_count = 0;      //Counter for "line"

  o_rna = T_to_U(input);

  GetReverseComplementSequence(input, &a_rna);
  a_rna = T_to_U(a_rna);

  //cout << a_rna << endl;
  //cout << o_rna << endl;
  //cout << endl;

  for(int i = 0; i < 3; i++){
    for( int j = 0; j < a_rna.length(); j++){
      
      line += a_rna[j];
      l_count++;

      //Resetting the "line" and counter
      if( l_count == 3){
        //cout << line << " ";
        
        //Adding codon to a vector
        v_1[i].push_back(line);
        
        l_count = 0;
        line = "";
      }
    
      //cout << a_rna[j];
    }
    
    //Reseating "line" and counter before each new sequence
    l_count = 0;
    line = "";
    
    //Deliting fisrt element in string
    a_rna.erase( a_rna.begin() );
    //cout << endl;
  }

  for(int i = 3; i < 6; i++){
    for( int j = 0; j < o_rna.length(); j++){
      
      line += o_rna[j];
      l_count++;

      //Resetting the "line" and counter
      if( l_count == 3){
        //cout << line << " ";
        
        //Adding codon to a vector
        v_1[i].push_back(line);
        
        l_count = 0;
        line = "";
      }
    
      //cout << a_rna[j];
    }
    
    //Reseating "line" and counter before each new sequence
    l_count = 0;
    line = "";
    
    //Deliting fisrt element in string
    o_rna.erase( o_rna.begin() );
    //cout << endl;
  }

  
  //Deletes empty element in string
  for( int i = 0; i < v_1.size(); i++){
    for( int j  = 0 ; j < v_1[i].size(); j++){
      if( v_1[i][j] == "" ){    //Deletes empty element in string
        v_1[i].erase( v_1[i].begin() + j);
      }
    }
  }
  

  return v_1;
}

/*
This function translates/converts a vector<string> (vector of codons) into a
string of amino acids using the genetic code
(see https://en.wikipedia.org/wiki/Genetic_code).

For example, the codons:
{"UUU", "GCC", "CAA"}
translates to:
F (Phenylalanine), A (Alanine), Q (Glutamine)
abreviated:
FAQ

http://www.soc-bdr.org/rds/authors/unit_tables_conversions_and_genetic_dictionaries/genetic_code_tables/


To make your lives easier, here's a list of the possible codons:
"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG",
"AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
"GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
"UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
"CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
"ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
"GUG", "UAG", "UGA", "UAA"

And there corresponding amino acids ("*" represents STOP codons,
more on them later):

"A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D",
"C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
"I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
"P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
"Y", "V", "V", "V", "V", "*", "*", "*"
*/
std::string Translate(const std::vector<std::string> & codon_sequence){
  string res = "";
  string key = "";
  
  vector <string> v_codon = {"GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC", "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA", "GUG", "UAG", "UGA", "UAA"};
  vector <string> v_acid = {"A", "A", "A", "A", "R", "R", "R", "R", "R", "R", "N", "N", "D", "D", "C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I", "I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P", "P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y", "Y", "V", "V", "V", "V", "*", "*", "*"};
  

  for( int i = 0 ; i < codon_sequence.size(); i++){
    key = codon_sequence[i];

    vector<string>::iterator itr = find(v_codon.begin(), v_codon.end(), key);

    if ( itr != v_codon.cend()){
      res += v_acid[distance(v_codon.begin(), itr)];
    }
  }

  return res;
}

/*
This function takes a DNA sequence and returns the longest possible
amino acid sequence / protein that is encoded by that sequence
(open reading frame). A valid open reading frame begins with the
codon AUG (the amino acid, Methionine (M)) and runs until a stop codon (*)
is encountered. There may be multiple open reading frames in a sequence, and
you need to check all six reading frames in order given by
get_reading_frames_as_codons. If there are ties for longest, favor the first
one found.

Return the longest open reading frame as an amino acid sequence. It must start
with an 'M' and end with a '*' with no other '*''s within.
*/
std::string GetLongestOpenReadingFrame(const std::string & DNA_sequence){
  
  //Adding two string, one for saving string from "line". Second, is for getting longest one
  string s_max = "";
  string s_test = "";
  
  string line = "";
  
  vector< vector<string> > v_1;

  //Getting two integers, one is for getting the length of new sequence in line. Second, is maximum length of sequence in DNA.
  int l_max = -999999999;
  int l_test = 0;

  v_1 = GetReadingFramesAsCodons(DNA_sequence);

  for( int i = 0; i < v_1.size(); i++){
    line = Translate(v_1[i]);
    //cout << line << endl;
    
    //Looking for the first position of char 'M'
    size_t pos = line.find('M');
    
    //Looking every char in line, starting from 'M' position.
    for( int j = pos; j < line.size(); j++){
      
      l_test++;
      s_test += line[j];
     
      if( line[j] == '*' ){
        break;
      
      }
    }

    //Checking for maximum/bigger
    if( l_test > l_max ){
      l_max = l_test;
      s_max = s_test;
    }

    //Reseating back
    l_test = 0;
    s_test = "";
  }

  return s_max;
}
