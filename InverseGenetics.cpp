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

/*
 * Nucleotides = letters A, C, U, and G
 * RNA Strands = Sequences of four nucleotides, represented as A, C, U, and G
 * Three Nucleotides = Codon
 * Codon = Encodes a specific amino acid
 * Protein = sequence of amino acids
 */

/* Function: listAllRNAStrandsFor(string protein,
 *                                Map<char, Set<string> >& codons);
 * Usage: listAllRNAStrandsFor("PARTY", codons);
 * ==================================================================
 * Given a protein and a map from amino acid codes to the codons for
 * that code, lists all possible RNA strands that could generate
 * that protein
 */
void listAllRNAStrandsFor(string protein, Map<char, Set<string> >& codons);


/* Function: loadCodonMap();
 * Usage: Map<char, Lexicon> codonMap = loadCodonMap();
 * ==================================================================
 * Loads the codon mapping table from a file.
 */
Map<char, Set<string> > loadCodonMap();

int main() {
    /* Load the codon map. */
    Map<char, Set<string> > codons = loadCodonMap();

    /* Get protein */
    string protein = "KWS";

    /* Assemble RNA Strands */
    listAllRNAStrandsFor(protein, codons);
    return 0;
}

void permuteCodons(Set<string>& particularCodons, Set<string>& rnaStrands) {
    if (particularCodons.size() == 0) {
        return;
    } else {
        string codon = particularCodons.first();
        // stuck.. need to figure out how to recursively
        // put the codon into each vector but not dup within the codon set
    }
}

void assembleAllRNAStrandsFor(string protein,
                           Map<char, Set<string> >& codons,
                           Set<string>& rnaStrands) {
    if (protein == "") {
        return;
    } else {
        char aminoAcid = protein[0];
        Set<string> particularCodons = codons.get(aminoAcid);
        // take each element in $particularCodons and do a cross-join
        //   with every existing element in $rnaStrands;
        //   in other words, we need to generate the permutations
        permuteCodons(particularCodons, rnaStrands);

        assembleAllRNAStrandsFor(protein.substr(1), codons, rnaStrands);
    }
}

void listAllRNAStrandsFor(string protein, Map<char, Set<string> >& codons) {
    Set<string> rnaStrands;
    assembleAllRNAStrandsFor(protein, codons, rnaStrands);
    for (string rnaStrand : rnaStrands) {
        cout << rnaStrands << endl;
    }
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