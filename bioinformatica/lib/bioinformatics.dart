import 'dart:math';

int delta(String x, String y) {
  return x == y ? 0 : 1;
}

int M(String seq1, String seq2, int i, int j, int k) {
  int sum = 0;
  for (int index = 0; index < k; index++) {
    sum += delta(seq1[i + index], seq2[j + index]);
  }
  return sum;
}

List<List<int>> makeMatrix(String seq1, String seq2, int k) {
  int n = seq1.length;
  int m = seq2.length;
  List<List<int>> matrix = List.generate(n - k + 1, (i) {
    return List.generate(m - k + 1, (j) {
      return M(seq1, seq2, i, j, k);
    });
  });
  return matrix;
}

void plotMatrix(List<List<int>> M, int t, String seq1, String seq2,
    {String nonblank = '\u25A0', String blank = ' '}) {
  print(' |$seq2');
  print('-' * (2 + seq2.length));
  for (int i = 0; i < seq1.length; i++) {
    String line = seq1[i] +
        '|' +
        M[i].map((s) => s < t ? nonblank : blank).join('');
    print(line);
  }
}

void dotplot(String seq1, String seq2, {int k = 1, int t = 1}) {
  List<List<int>> M = makeMatrix(seq1, seq2, k);
  plotMatrix(M, t, seq1, seq2);
}

double gcContent(String seq) {
  int countG = seq.split('G').length - 1;
  int countC = seq.split('C').length - 1;
  return (countG + countC) / seq.length * 100;
}

double atContent(String seq) {
  int countA = seq.split('A').length - 1;
  int countT = seq.split('T').length - 1;
  return (countA + countT) / seq.length * 100;
}

double aromaticProportion(Chem.Mol mol) {
  List<bool> aromaticAtoms = List.generate(mol.getNumAtoms(), (i) {
    return mol.getAtomWithIdx(i).getIsAromatic();
  });
  int aaCount = aromaticAtoms.where((i) => i == true).length;
  int heavyAtom = Descriptors.HeavyAtomCount(mol);
  double AR = aaCount / heavyAtom;
  return AR;
}

class MolecularDescriptor {
  double molLogP;
  double molWt;
  double numRotatableBonds;
  double aromaticProportion;

  MolecularDescriptor(
      this.molLogP, this.molWt, this.numRotatableBonds, this.aromaticProportion);
}

List<MolecularDescriptor> generate(List<String> smilesList) {
  List<MolecularDescriptor> descriptors = [];

  for (String smiles in smilesList) {
    Chem.Mol? mol = Chem.MolFromSmiles(smiles);
    if (mol != null) {
      double molLogP = Descriptors.MolLogP(mol);
      double molWt = Descriptors.MolWt(mol);
      double numRotatableBonds = Descriptors.NumRotatableBonds(mol);
      double aromaticProp = aromaticProportion(mol);

      MolecularDescriptor descriptor =
          MolecularDescriptor(molLogP, molWt, numRotatableBonds, aromaticProp);
      descriptors.add(descriptor);
    }
  }

  return descriptors;
}

Map<String, int> dnaNucleotideCount(String seq) {
  Map<String, int> countMap = {
    'A': seq.split('A').length - 1,
    'T': seq.split('T').length - 1,
    'G': seq.split('G').length - 1,
    'C': seq.split('C').length - 1,
  };
  return countMap;
}
