class BioInformatica {
    static delta(x, y) {
      return x === y ? 0 : 1;
    }
  
    static M(seq1, seq2, i, j, k) {
      let sum = 0;
      for (let index = 0; index < k; index++) {
        sum += this.delta(seq1.charAt(i + index), seq2.charAt(j + index));
      }
      return sum;
    }
  
    static makeMatrix(seq1, seq2, k) {
      const n = seq1.length;
      const m = seq2.length;
      const matrix = new Array(n - k + 1);
  
      for (let i = 0; i < n - k + 1; i++) {
        matrix[i] = new Array(m - k + 1);
        for (let j = 0; j < m - k + 1; j++) {
          matrix[i][j] = this.M(seq1, seq2, i, j, k);
        }
      }
  
      return matrix;
    }
  
    static plotMatrix(M, t, seq1, seq2, nonblank, blank) {
      console.log(" |" + seq2);
      console.log("-".repeat(2 + seq2.length));
  
      for (let i = 0; i < seq1.length; i++) {
        process.stdout.write(seq1.charAt(i) + "|");
        for (let j = 0; j < seq2.length; j++) {
          const s = M[i][j] < t ? nonblank : blank;
          process.stdout.write(s);
        }
        console.log();
      }
    }
  
    static dotplot(seq1, seq2, k, t) {
      const M = this.makeMatrix(seq1, seq2, k);
      this.plotMatrix(M, t, seq1, seq2, '\u25A0', ' ');
    }
  
    static gcContent(seq) {
      const countG = (seq.match(/G/g) || []).length;
      const countC = (seq.match(/C/g) || []).length;
      return ((countG + countC) / seq.length) * 100;
    }
  
    static atContent(seq) {
      const countA = (seq.match(/A/g) || []).length;
      const countT = (seq.match(/T/g) || []).length;
      return ((countA + countT) / seq.length) * 100;
    }
  
    static aromaticProportion(smiles) {
      try {
        const mol = RDKit.MolFromSmiles(smiles);
        const numAtoms = mol.getNumAtoms();
        let aromaticAtomCount = 0;
  
        for (let i = 0; i < numAtoms; i++) {
          const atom = mol.getAtomWithIdx(i);
          if (atom.getIsAromatic()) {
            aromaticAtomCount++;
          }
        }
  
        return (aromaticAtomCount / numAtoms) * 100;
      } catch (e) {
        console.error(e);
        return 0.0;
      }
    }
  
    static generate(smilesList) {
      const descriptors = Table.create("Molecule Descriptors")
        .addColumns(
          StringColumn.create("SMILES"),
          DoubleColumn.create("MolLogP"),
          DoubleColumn.create("MolWt"),
          DoubleColumn.create("NumRotatableBonds"),
          DoubleColumn.create("AromaticProportion")
        );
  
      for (const smiles of smilesList) {
        try {
          const mol = RDKit.MolFromSmiles(smiles);
          const molLogP = Descriptors.MolLogP(mol);
          const molWt = Descriptors.MolWt(mol);
          const numRotatableBonds = Descriptors.NumRotatableBonds(mol);
          const aromaticProp = this.aromaticProportion(smiles);
  
          descriptors.append(
            smiles,
            molLogP,
            molWt,
            numRotatableBonds,
            aromaticProp
          );
        } catch (e) {
          console.error(e);
        }
      }
  
      return descriptors;
    }
  
    static dnaNucleotideCount(seq) {
      const countMap = {
        'A': 0,
        'T': 0,
        'G': 0,
        'C': 0
      };
  
      for (const ch of seq) {
        countMap[ch]++;
      }
  
      return countMap;
    }
  
    static main() {
      // Example usage
      const seq1 = "ATCGATCG";
      const seq2 = "ATGGATCG";
      const k = 1;
      const t = 1;
  
      this.dotplot(seq1, seq2, k, t);
      const gc = this.gcContent(seq1);
      const at = this.atContent(seq1);
  
      console.log("GC Content: " + gc);
      console.log("AT Content: " + at);
  
      // Example usage for generate function
      const smilesList = ["CCO", "CCC", "CCN"];
      const descriptors = this.generate(smilesList);
      console.log(descriptors.toString());
    }
  }
  
  // Example usage
  BioinformaticsFunctions.main();
  