import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;
import org.RDKit.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BioInformatics {

    // Define the delta function
    public static int delta(char x, char y) {
        return (x == y) ? 0 : 1;
    }

    // Define the M function
    public static int M(String seq1, String seq2, int i, int j, int k) {
        int sum = 0;
        for (int index = 0; index < k; index++) {
            sum += delta(seq1.charAt(i + index), seq2.charAt(j + index));
        }
        return sum;
    }

    // Define the makeMatrix function
    public static int[][] makeMatrix(String seq1, String seq2, int k) {
        int n = seq1.length();
        int m = seq2.length();
        int[][] matrix = new int[n - k + 1][m - k + 1];

        for (int i = 0; i < n - k + 1; i++) {
            for (int j = 0; j < m - k + 1; j++) {
                matrix[i][j] = M(seq1, seq2, i, j, k);
            }
        }

        return matrix;
    }

    // Define the plotMatrix function
    public static void plotMatrix(int[][] M, int t, String seq1, String seq2, char nonblank, char blank) {
        System.out.println(" |" + seq2);
        System.out.println("-".repeat(2 + seq2.length()));

        for (int i = 0; i < seq1.length(); i++) {
            System.out.print(seq1.charAt(i) + "|");
            for (int j = 0; j < seq2.length(); j++) {
                char s = M[i][j] < t ? nonblank : blank;
                System.out.print(s);
            }
            System.out.println();
        }
    }

    // Define the dotplot function
    public static void dotplot(String seq1, String seq2, int k, int t) {
        int[][] M = makeMatrix(seq1, seq2, k);
        plotMatrix(M, t, seq1, seq2, '\u25A0', ' ');
    }

    public static double gcContent(String seq) {
        long countG = seq.chars().filter(ch -> ch == 'G').count();
        long countC = seq.chars().filter(ch -> ch == 'C').count();
        return (double) (countG + countC) / seq.length() * 100;
    }

    public static double atContent(String seq) {
        long countA = seq.chars().filter(ch -> ch == 'A').count();
        long countT = seq.chars().filter(ch -> ch == 'T').count();
        return (double) (countA + countT) / seq.length() * 100;
    }

    public static double aromaticProportion(String smiles) {
        try {
            RWMol mol = RWMol.MolFromSmiles(smiles);
            int numAtoms = mol.getNumAtoms();
            int aromaticAtomCount = 0;

            for (int i = 0; i < numAtoms; i++) {
                Atom atom = mol.getAtomWithIdx(i);
                if (atom.getIsAromatic()) {
                    aromaticAtomCount++;
                }
            }

            return (double) aromaticAtomCount / numAtoms;
        } catch (Exception e) {
            e.printStackTrace();
            return 0.0;
        }
    }

    public static Table generate(List<String> smilesList) {
        Table descriptors = Table.create("Molecule Descriptors");

        // Create columns with the same number of rows as smilesList
        StringColumn smilesColumn = StringColumn.create("SMILES");
        DoubleColumn molLogPColumn = DoubleColumn.create("MolLogP");
        DoubleColumn molWtColumn = DoubleColumn.create("MolWt");
        DoubleColumn numRotatableBondsColumn = DoubleColumn.create("NumRotatableBonds");
        DoubleColumn aromaticPropColumn = DoubleColumn.create("AromaticProportion");

        for (String smiles : smilesList) {
            try {
                RWMol mol = RWMol.MolFromSmiles(smiles);
                double molLogP = Descriptors.MolLogP(mol);
                double molWt = Descriptors.MolWt(mol);
                double numRotatableBonds = Descriptors.NumRotatableBonds(mol);
                double aromaticProp = aromaticProportion(smiles);

                // Add data to the columns
                smilesColumn.append(smiles);
                molLogPColumn.append(molLogP);
                molWtColumn.append(molWt);
                numRotatableBondsColumn.append(numRotatableBonds);
                aromaticPropColumn.append(aromaticProp);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        // Add the columns to the table
        descriptors.addColumns(smilesColumn, molLogPColumn, molWtColumn, numRotatableBondsColumn, aromaticPropColumn);

        return descriptors;
    }

    public static Map<Character, Integer> dnaNucleotideCount(String seq) {
        Map<Character, Integer> countMap = new HashMap<>();
        countMap.put('A', 0);
        countMap.put('T', 0);
        countMap.put('G', 0);
        countMap.put('C', 0);

        for (char ch : seq.toCharArray()) {
            countMap.put(ch, countMap.get(ch) + 1);
        }

        return countMap;
    }
}
