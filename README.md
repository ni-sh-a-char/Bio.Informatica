# Bioinformatica

A Dart library that provides various bioinformatics functions, including sequence analysis and molecular descriptors calculation.

## Features

- Calculate sequence similarity and create dot plots.
- Calculate GC and AT content of DNA sequences.
- Calculate molecular descriptors for chemical compounds.
- Count DNA nucleotides.

## Installation

Add this to your Dart project's `pubspec.yaml`:

```
dependencies:
  bioinformatica: ^1.0.0
```

## Usage
```
import 'package:my_bioinformatics_library/my_bioinformatics_library.dart';

void main() {
  // Example usage:
  
  // Calculate sequence similarity and create a dot plot
  dotplot("ATCG", "ATGG");

  // Calculate GC and AT content of DNA sequences
  double gc = gcContent("ATCGATCG");
  double at = atContent("ATCGATCG");

  // Calculate molecular descriptors for chemical compounds
  List<String> smilesList = ["CCO", "CCC", "CCN"];
  List<MolecularDescriptor> descriptors = generate(smilesList);

  // Count DNA nucleotides
  Map<String, int> nucleotideCount = dnaNucleotideCount("ATCGATCG");
}
```