import 'package:test/test.dart';
import 'package:bioinformatica/bioinformatics.dart';

void main() {
  group('Bioinformatics Library Tests', () {
    test('Delta Function', () {
      expect(delta('A', 'A'), equals(0));
      expect(delta('A', 'T'), equals(1));
      expect(delta('G', 'C'), equals(1));
    });

    test('GC Content', () {
      expect(gcContent('ATCGATCG'), equals(50.0));
    });

    test('Aromatic Proportion', () {
      // You may need to provide a valid Chem.Mol for testing this function.
      // Mock a valid Chem.Mol for testing purposes.
      final mockMol = MockChemMol(); // Replace MockChemMol with your mock.

      expect(aromaticProportion(mockMol), equals(0.25));
    });

    // Add more test cases for your other functions.
  });
}

class MockChemMol {
  // Mock this class to simulate Chem.Mol for testing aromaticProportion.
}
