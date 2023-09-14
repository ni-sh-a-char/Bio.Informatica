markdown

# BioInformatica

BioInformatica is a JavaScript package that provides various bioinformatics functions for working with biological sequences and chemical data.

## Installation

You can install BioInformatica using npm:

```
npm install bioinformatica
```

## Usage

To use BioInformatica in your JavaScript/Node.js project, you can import the functions you need as follows:


```
const bioinformatica = require('bioinformatica');

// Example usage:
const seq1 = "ATCGATCG";
const seq2 = "ATGGATCG";
const k = 1;
const t = 1;

bioinformatica.dotplot(seq1, seq2, k, t);
const gc = bioinformatica.gcContent(seq1);
const at = bioinformatica.atContent(seq1);

console.log("GC Content:", gc);
console.log("AT Content:", at);

```

## Functions

BioInformatica provides the following functions:

    dotplot(seq1, seq2, k, t): Generates a dot plot for two sequences.
    gcContent(seq): Calculates the GC content of a DNA sequence.
    atContent(seq): Calculates the AT content of a DNA sequence.
    generate(smilesList): Generates molecular descriptors for a list of SMILES strings.
    dnaNucleotideCount(seq): Counts the occurrences of DNA nucleotides in a sequence.

## Dependencies

BioInformatica relies on the following dependencies:

    tablesaw-js: A JavaScript library for working with tabular data.
    rdkit-js: A JavaScript library for cheminformatics (for chemical data functions).

You can install these dependencies using npm as described in the Installation section.