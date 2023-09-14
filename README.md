# Bio.Informatica now as library for Java

# Usage

Usage
Sequence Analysis

## Delta Function
```
char x = 'A';
char y = 'T';
int result = Bioinformatics.delta(x, y);
System.out.println("Delta: " + result); // Output: Delta: 1
```

## GC Content

```

String seq = "ATGCGTA";
double gc = Bioinformatics.gcContent(seq);
System.out.println("GC Content: " + gc + "%"); // Output: GC Content: 57.142857%
```

## DNA Nucleotide Count

```
String seq = "ATGCGTA";
Map<Character, Integer> countMap = Bioinformatics.dnaNucleotideCount(seq);
System.out.println("A: " + countMap.get('A')); // Output: A: 2
System.out.println("T: " + countMap.get('T')); // Output: T: 1
System.out.println("G: " + countMap.get('G')); // Output: G: 3
System.out.println("C: " + countMap.get('C')); // Output: C: 1
```
## Molecular Descriptor Generation
Generate Descriptors from SMILES

```
List<String> smilesList = Arrays.asList("CCO", "CCC", "CCN");
Table descriptors = Bioinformatics.generate(smilesList);
System.out.println(descriptors.toString());
```