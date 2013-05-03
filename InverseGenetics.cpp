/*
 * File: InverseGenetics.cpp
 * --------------------------
 * Section: SCPD, Aaron Broder <abroder@stanford.edu>
 * Copyright 2013 Eric Beach <ebeach@google.com>
 * Assignment 3 - Pt. 3 - Inverse Genetics
 */

#include <iostream>
#include <string>
#include <fstream>
#include "set.h"
#include "map.h"
#include "console.h"
using namespace std;

////////// FUNCTION PROTOTYPES //////////
/*
 * Given a protein and a map from amino acid codes to the codons for
 * that code, lists all possible RNA strands that could generate
 * that protein.
 */
void listAllRNAStrandsFor(string protein, Map<char, Set<string> >& codons);

/* Loads the codon mapping table from a file. */
Map<char, Set<string> > loadCodonMap();

/* Assemble */
void assembleAllRNAStrandsFor(string protein,
                              const Map<char, Set<string> >& codons,
                              Vector<string> rnaStrands);

////////// FUNCTION IMPLEMENTATION //////////
/*
 * Main function to run tests and list all RNA strands for a given protein.
 */
int main() {    
    /* Load the codon map. */
    Map<char, Set<string> > codons = loadCodonMap();
    
    /* Get protein */
    string protein = "SK";
    
    /* Assemble RNA Strands */
    listAllRNAStrandsFor(protein, codons);
    return 0;
}

/*
 * Assemble RNS strands that could represent a protein with a given amino
 *   acid sequence.
 */
void assembleAllRNAStrandsFor(string protein,
                              const Map<char, Set<string> >& codons,
                              Vector<string> rnaStrands) {
    if (protein.size() == 0) {
        cout << rnaStrands << endl;
        return;
    } else {
        char aminoAcid = protein[0];
        Set<string> particularCodons = codons.get(aminoAcid);
        for (string codon : particularCodons) {
            Vector<string> newRnaStrands = rnaStrands;
            newRnaStrands += codon;
            assembleAllRNAStrandsFor(protein.substr(1), codons, newRnaStrands);
        }
    }
}

/*
 * Print out all RNA strands for a given amino acid sequence.
 */
void listAllRNAStrandsFor(string protein, Map<char, Set<string> >& codons) {
    Vector<string> rnaStrands;
    assembleAllRNAStrandsFor(protein, codons, rnaStrands);
}

/* You do not need to change this function. */
Map<char, Set<string> > loadCodonMap() {
    ifstream input("codons.txt");
    Map<char, Set<string> > result;
    
    /* The current codon / protein combination. */
    string codon;
    char protein;
    
    /* Continuously pull data from the file until all data has been
     * read.
     */
    while (input >> codon >> protein) {
        result[protein] += codon;
    }
    
    return result;
}
