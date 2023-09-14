import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import tech.tablesaw.api.Table;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class BioinformaticsTest {

    @Test
    void testDelta() {
        assertEquals(0, Bioinformatics.delta('A', 'A'));
        assertEquals(1, Bioinformatics.delta('A', 'T'));
        assertEquals(1, Bioinformatics.delta('G', 'C'));
    }

    @Test
    void testM() {
        assertEquals(0, Bioinformatics.M("ATG", "ATG", 0, 0, 3));
        assertEquals(2, Bioinformatics.M("ATG", "CTA", 0, 0, 3));
    }

    @Test
    void testMakeMatrix() {
        int[][] matrix = Bioinformatics.makeMatrix("ATG", "CTA", 3);
        assertEquals(2, matrix[0][0]);
        assertEquals(2, matrix[1][0]);
    }

    @Test
    void testGCContent() {
        double gc = Bioinformatics.gcContent("ATGCGTA");
        assertEquals(57.142857, gc, 0.000001);
    }

    @Test
    void testATContent() {
        double at = Bioinformatics.atContent("ATGCGTA");
        assertEquals(42.857142, at, 0.000001);
    }

    @Test
    void testGenerate() {
        List<String> smilesList = Arrays.asList("CCO", "CCC", "CCN");
        Table descriptors = Bioinformatics.generate(smilesList);

        assertNotNull(descriptors);
        assertEquals(3, descriptors.rowCount());
        assertEquals(5, descriptors.columnCount());
    }

    @Test
    void testDNANucleotideCount() {
        String seq = "ATGCGTA";
        Map<Character, Integer> countMap = Bioinformatics.dnaNucleotideCount(seq);

        assertNotNull(countMap);
        assertEquals(2, countMap.get('A'));
        assertEquals(1, countMap.get('T'));
        assertEquals(3, countMap.get('G'));
        assertEquals(1, countMap.get('C'));
    }
}
