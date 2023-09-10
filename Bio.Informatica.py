import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
import neatbio.sequtils as utils
from collections import Counter
import io
import Bio.Data.CodonTable

# Data Vis Pkgs
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import numpy as np

import pandas as pd
import pickle
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Descriptors

import altair as alt

import re

import subprocess
import os
import base64

import mols2grid
import streamlit.components.v1 as components
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors

def delta(x,y):
    return 0 if x == y else 1


def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]


def plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)


def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice


# Convert to Fxn
def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()


def gc_content(seq):
	result = float(str(seq).count('G') + str(seq).count('C'))/len(seq) * 100
	return result

def at_content(seq):
	result = float(str(seq).count('A') + str(seq).count('T'))/len(seq) * 100
	return result

def AromaticProportion(m):
  aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
  aa_count = []
  for i in aromatic_atoms:
    if i==True:
      aa_count.append(1)
  AromaticAtom = sum(aa_count)
  HeavyAtom = Descriptors.HeavyAtomCount(m)
  AR = AromaticAtom/HeavyAtom
  return AR

def generate(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData= np.arange(1,1)
    i=0
    for mol in moldata:

        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MolWt = Descriptors.MolWt(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
        desc_AromaticProportion = AromaticProportion(mol)

        row = np.array([desc_MolLogP,
                        desc_MolWt,
                        desc_NumRotatableBonds,
                        desc_AromaticProportion])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1

    columnNames=["MolLogP","MolWt","NumRotatableBonds","AromaticProportion"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors

def DNA_nucleotide_count(seq):
            d = dict([
                    ('A',seq.count('A')),
                    ('T',seq.count('T')),
                    ('G',seq.count('G')),
                    ('C',seq.count('C'))
                    ])
            return d



def main():
    
    # Logo image
    image = Image.open('logo.png')

    st.image(image, use_column_width=True)

    st.title("Bio.Informatica")

    menu = ["Intro", "DNA Sequence", "DNA Nucleotides", "Bioactivity Prediction", "DotPlot", "Solubility Prediction", "Molecular Descriptor Calculator", "FDA Approved Drugs", "Humsafar"]
    choice = st.sidebar.selectbox("Select Activity", menu)

    if choice == "Intro":
        st.subheader("A  Bundelkhand University, Jhansi  initiative to study Bioinformatics, Computational Biology and Computational Chemistry")
    elif choice == "DNA Sequence":
        st.subheader("DNA Sequence Analysis")

        seq_file = st.file_uploader("Upload FASTA File",type=["fasta", "fa"])

        if seq_file is not None:
            byte_str = seq_file.read()
            text_obj = byte_str.decode("UTF-8")
            dna_record = SeqIO.read(io.StringIO(text_obj),"fasta")
            #st.write(dna_record)
            dna_seq = dna_record.seq
        
            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record.description)
            elif details == "Sequence":
                st.write(dna_record.seq)


            # Nucleotide Frequencies
            st.subheader("Nucleotide Frequency")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_color = st.color_picker("Adenine Color")
            thymine_color = st.color_picker("Thymine Color")
            guanine_color = st.color_picker("Guanine Color")
            cytosine_color = st.color_picker("Cytosine Color")

            if st.button("Plot Freq"):
                barlist = plt.bar(dna_freq.keys(),dna_freq.values())
                barlist[2].set_color(adenine_color)
                barlist[3].set_color(thymine_color)
                barlist[1].set_color(guanine_color)
                barlist[0].set_color(cytosine_color)

                st.set_option('deprecation.showPyplotGlobalUse', False)
                st.pyplot()

            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content":gc_score,"AT Content":at_score})

            # Nucleotide Count
            nt_count = st.text_input("Type Nucleotide Alphabet To Find Frequency")
            st.header("Number of  {} Nucleotide is ::  {}".format((nt_count),str(dna_seq).count(nt_count)))


    elif choice == "DNA Nucleotides":

        st.write("""
        # Nucleotide composition of query DNA!
        ***
        """)


        ######################
        # Input Text Box
        ######################

        #st.sidebar.header('Enter DNA sequence')
        st.header('Enter DNA sequence')

        sequence_input = ">DNA Query \n"

        #sequence = st.sidebar.text_area("Sequence input", sequence_input, height=250)
        sequence = st.text_area("Sequence Input  (Remove stop codons if any)", sequence_input, height=250)
        sequence = sequence.splitlines()
        sequence = sequence[1:] # Skips the sequence name (first line)
        sequence = ''.join(sequence) # Concatenates list to string

        st.write("""
        ***
        """)

        ## Prints the input DNA sequence
        st.header('INPUT (DNA Query)')
        sequence

        ## DNA nucleotide count
        st.header('OUTPUT (DNA Nucleotide Count)')

        ### 1. Print dictionary
        st.subheader('1. Print dictionary')
        
        X = DNA_nucleotide_count(sequence)

        #X_label = list(X)
        #X_values = list(X.values())

        X

        ### 2. Print text
        st.subheader('2. Print text')
        st.write('There are  ' + str(X['A']) + ' adenine (A)')
        st.write('There are  ' + str(X['T']) + ' thymine (T)')
        st.write('There are  ' + str(X['G']) + ' guanine (G)')
        st.write('There are  ' + str(X['C']) + ' cytosine (C)')

        ### 3. Display DataFrame
        st.subheader('3. Display DataFrame')
        df = pd.DataFrame.from_dict(X, orient='index')
        df = df.rename({0: 'count'}, axis='columns')
        df.reset_index(inplace=True)
        df = df.rename(columns = {'index':'nucleotide'})
        st.write(df)

        ### 4. Display Bar Chart using Altair
        st.subheader('4. Display Bar chart')
        p = alt.Chart(df).mark_bar().encode(
            x='nucleotide',
            y='count'
        )
        p = p.properties(
            width=alt.Step(80)  # controls width of bar.
        )
        st.write(p)

        # Protein Synthesis
        st.subheader("Protein Synthesis")
        p1 = Bio.Seq.translate(sequence)
        aa_freq = Counter(str(p1))

        if st.checkbox("Transcription"):
            st.write(Bio.Seq.transcribe(sequence))

        if st.checkbox("Back-Transcription"):
            st.write(Bio.Seq.back_transcribe(sequence))

        elif st.checkbox("Translation"):
            st.write(Bio.Seq.translate(sequence))

        elif st.checkbox("Complement"):
            st.write(Bio.Seq.complement(sequence))

        elif st.checkbox("Reverse Complement"):
            st.write(Bio.Seq.reverse_complement(sequence))

        elif st.checkbox("AA Frequency"):
            st.write(aa_freq)

        elif st.checkbox("Plot AA Frequency"):
            st.set_option('deprecation.showPyplotGlobalUse', False)
            aa_color = st.color_picker("Pick An Amino Acid Color")
            # barlist = plt.bar(aa_freq.keys(),aa_freq.values())
            # barlist[2].set_color(aa_color)
            plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
            st.pyplot()

        elif st.checkbox("Full Amino Acid Name"):
            aa_name = str(p1).replace("*", "")
            aa3 = utils.convert_1to3(aa_name)
            st.write(aa_name)            
            st.write("======================")
            st.write(aa3)

            st.write("======================")
            st.write(utils.get_acid_name(aa3))
            # Top Most Common Amino

    elif choice == "DotPlot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file1 = st.file_uploader("Upload 1st FASTA File",type=["fasta", "fa"])
        seq_file2 = st.file_uploader("Upload 2nd FASTA File",type=["fasta", "fa"])

        if seq_file1 and seq_file2 is not None:
            byte_str1 = seq_file1.read()
            text_obj1 = byte_str1.decode("UTF-8")
            byte_str2 = seq_file2.read()
            text_obj2 = byte_str2.decode("UTF-8")
            dna_record1 = SeqIO.read(io.StringIO(text_obj1),"fasta")
            dna_record2 = SeqIO.read(io.StringIO(text_obj2),"fasta")
            #st.write(dna_record)
            dna_seq1 = dna_record1.seq
            dna_seq2 = dna_record2.seq
        
            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(dna_record1.description)
                st.write("=====================")
                st.write(dna_record2.description)
            elif details == "Sequence":
                st.write(dna_record1.seq)
                st.write("=====================")
                st.write(dna_record2.seq)

            cus_limit = st.number_input("Select Max number of Nucleotide",10,200,50)
            if st.button("Dot Plot"):
                st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                dotplotx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit])

                st.pyplot()

    elif choice == "Bioactivity Prediction":            

        # Molecular descriptor calculator
        def desc_calc():
            # Performs the descriptor calculation
            bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            os.remove('molecule.smi')

        # File download
        def filedownload(df):
            csv = df.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
            href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
            return href

        # Model building
        def build_model(input_data):
            # Reads in saved regression model
            load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
            # Apply model to make predictions
            prediction = load_model.predict(input_data)
            st.header('**Prediction output**')
            prediction_output = pd.Series(prediction, name='pIC50')
            molecule_name = pd.Series(load_data[1], name='molecule_name')
            df = pd.concat([molecule_name, prediction_output], axis=1)
            st.write(df)
            st.markdown(filedownload(df), unsafe_allow_html=True)

        # Logo image
        #image = Image.open('logo.png')

        #st.image(image, use_column_width=True)

        # Page title
        st.markdown("""
        # Bioactivity Prediction (Acetylcholinesterase)

        This section allows you to predict the bioactivity towards inhibting the `Acetylcholinesterase` enzyme. `Acetylcholinesterase` is a drug target for Alzheimer's disease.

        **About**
        - Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
        ---
        """)

        # Sidebar
        with st.sidebar.header('1. Upload your CSV data'):
            uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
            st.sidebar.markdown("""
        [Example input file](https://raw.githubusercontent.com/PIYUSH-MISHRA-00/data/main/example_acetylcholinesterase.txt)
        """)

        if st.sidebar.button('Predict'):
            load_data = pd.read_table(uploaded_file, sep=' ', header=None)
            load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

            st.header('**Original input data**')
            st.write(load_data)

            with st.spinner("Calculating descriptors..."):
                desc_calc()

            # Read in calculated descriptors and display the dataframe
            st.header('**Calculated molecular descriptors**')
            desc = pd.read_csv('descriptors_output.csv')
            st.write(desc)
            st.write(desc.shape)

            # Read descriptor list used in previously built model
            st.header('**Subset of descriptors from previously built models**')
            Xlist = list(pd.read_csv('descriptor_list.csv').columns)
            desc_subset = desc[Xlist]
            st.write(desc_subset)
            st.write(desc_subset.shape)

            # Apply trained model to make prediction on query compounds
            build_model(desc_subset)
        else:
            st.info('Upload input data in the sidebar to start!')

    elif choice == "Solubility Prediction":


        ######################
        # Page Title
        ######################

        #image = Image.open('solubility-logo.jpg')

        #st.image(image, use_column_width=True)

        st.write("""
        # Molecular Solubility Prediction
        This section predicts the **Solubility (LogS)** values of molecules!
        Data obtained from the John S. Delaney. [ESOL:  Estimating Aqueous Solubility Directly from Molecular Structure](https://pubs.acs.org/doi/10.1021/ci034243x). ***J. Chem. Inf. Comput. Sci.*** 2004, 44, 3, 1000-1005.
        ***
        """)

        ######################
        # Input molecules (Side Panel)
        ######################

        st.sidebar.header('User Input Features')

        ## Read SMILES input
        SMILES_input = "CCCCC\nCCC\nCN"

        SMILES = st.sidebar.text_area("SMILES input", SMILES_input)
        SMILES = "C\n" + SMILES #Adds C as a dummy, first item
        SMILES = SMILES.split('\n')

        st.header('Input SMILES')
        SMILES[1:] # Skips the dummy first item

        ## Calculate molecular descriptors
        st.header('Computed molecular descriptors')
        X = generate(SMILES)
        X[1:] # Skips the dummy first item

        ######################
        # Pre-built model
        ######################

        # Reads in saved model
        load_model = pickle.load(open('solubility_model.pkl', 'rb'))

        # Apply model to make predictions
        prediction = load_model.predict(X)
        #prediction_proba = load_model.predict_proba(X)

        st.header('Predicted LogS values')
        prediction[1:] # Skips the dummy first item

    elif choice == "Molecular Descriptor Calculator":

        # Molecular descriptor calculator
        def desc_calc():
            # Performs the descriptor calculation
            bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/%s -dir ./ -file descriptors_output.csv" % selected_fp
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            # Read in calculated descriptors and display the dataframe
            st.subheader('Calculated molecular descriptors')
            desc = pd.read_csv('descriptors_output.csv')
            st.write(desc)
            st.markdown(filedownload(desc), unsafe_allow_html=True)
            # Write the data dimension (number of molecules and descriptors)
            nmol = desc.shape[0]
            ndesc = desc.shape[1]
            st.info('Selected fingerprint: ' + user_fp)
            st.info('Number of molecules: ' + str(nmol))
            st.info('Number of descriptors: ' + str(ndesc-1))
            os.remove('molecule.smi')

        # File download
        def filedownload(df):
            csv = df.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
            href = f'<a href="data:file/csv;base64,{b64}" download="descriptor_{user_fp}.csv">Download CSV File</a>'
            return href

        # Page title
        st.markdown("""
        # Molecular Descriptor Calculator

        This section allows you to calculate descriptors of molecules (aka **molecular descriptors**) that you can use for computational drug discovery projects such as for the construction of quantitative structure-activity/property relationship (QSAR/QSPR) models.

        This section gives you access to 12 **molecular fingerprints** (`AtomPairs2D`, `AtomPairs2DCount`, `CDK`, `CDKextended`, `CDKgraphonly`, `EState`, `KlekotaRoth`, `KlekotaRothCount`, `MACCS`, `PubChem`, `Substructure` and `SubstructureCount`).

        **About**
        - Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) software.
        - Yap CW. [PaDEL‐descriptor: An open source software to calculate molecular descriptors and fingerprints](https://doi.org/10.1002/jcc.21707). ***J Comput Chem*** 32 (2011) 1466-1474.
        ---
        """)

        # Sidebar
        with st.sidebar.header('1. Upload your CSV data'):
            uploaded_file = st.sidebar.file_uploader("Upload your input CSV file", type=["csv"])
            st.sidebar.markdown("""
        [Example CSV input file](https://raw.githubusercontent.com/PIYUSH-MISHRA-00/data/main/acetylcholinesterase_04_bioactivity_data_3class_pIC50.csv)
        """)

        with st.sidebar.header('2. Enter column names for 1) Molecule ID and 2) SMILES'):
            name_mol = st.sidebar.text_input('Enter column name for Molecule ID', 'molecule_chembl_id')
            name_smiles = st.sidebar.text_input('Enter column name for SMILES', 'canonical_smiles')

        with st.sidebar.header('3. Set parameters'):
            # Select fingerprint
            fp_dict = {'AtomPairs2D':'AtomPairs2DFingerprinter.xml',
                    'AtomPairs2DCount':'AtomPairs2DFingerprintCount.xml',
                    'CDK':'Fingerprinter.xml',
                    'CDKextended':'ExtendedFingerprinter.xml',
                    'CDKgraphonly':'GraphOnlyFingerprinter.xml',
                    'EState':'EStateFingerprinter.xml',
                    'KlekotaRoth':'KlekotaRothFingerprinter.xml',
                    'KlekotaRothCount':'KlekotaRothFingerprintCount.xml',
                    'MACCS':'MACCSFingerprinter.xml',
                    'PubChem':'PubchemFingerprinter.xml',
                    'Substructure':'SubstructureFingerprinter.xml',
                    'SubstructureCount':'SubstructureFingerprintCount.xml'}
            user_fp = st.sidebar.selectbox('Choose fingerprint to calculate', list(fp_dict.keys()) )
            selected_fp = fp_dict[user_fp]

            # Set number of molecules to compute
            df0 = pd.read_csv('acetylcholinesterase_04_bioactivity_data_3class_pIC50.csv')
            all_mol = df0.shape[0]
            number2calc = st.sidebar.slider('How many molecules to compute?', min_value=10, max_value=all_mol, value=10, step=10)


        if uploaded_file is not None:
            # Read CSV data
            @st.cache
            def load_csv():
                csv = pd.read_csv(uploaded_file).iloc[:number2calc,1:]
                return csv
            df = load_csv()
            df2 = pd.concat([df[name_smiles], df[name_mol]], axis=1)
            # Write CSV data
            df2.to_csv('molecule.smi', sep = '\t', header = False, index = False)
            st.subheader('Initial data from CSV file')
            st.write(df)
            st.subheader('Formatted as PADEL input file')
            st.write(df2)
            with st.spinner("Calculating descriptors..."):
                desc_calc()

        else:
            st.info('Awaiting for CSV file to be uploaded.')
            if st.button('Press to use Example Dataset'):
                # Read CSV data
                @st.cache
                def load_data():
                    # number2calc specifies the number of molecules to compute
                    df = pd.read_csv('acetylcholinesterase_04_bioactivity_data_3class_pIC50.csv').iloc[:number2calc,1:]
                    return df
                df = load_data()
                df2 = pd.concat([df[name_smiles], df[name_mol]], axis=1)
                # Write CSV data
                df2.to_csv('molecule.smi', sep = '\t', header = False, index = False)
                st.subheader('Initial data from CSV file')
                st.write(df)
                st.subheader('Formatted as PADEL input file')
                st.write(df2)
                with st.spinner("Calculating descriptors..."):
                    desc_calc()

    elif choice == "FDA Approved Drugs":
        st.title("Filter FDA Approved Drugs by Lipinski's Rule-of-Five")

        @st.cache(allow_output_mutation=True)
        def download_dataset():
            """Loads once then cached for subsequent runs"""
            df = pd.read_csv(
                "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", sep="\t"
            ).dropna()
            return df

        # Calculate descriptors
        def calc_mw(smiles_string):
            """Given a smiles string (ex. C1CCCCC1), calculate and return the molecular weight"""
            mol = Chem.MolFromSmiles(smiles_string)
            return ExactMolWt(mol)

        def calc_logp(smiles_string):
            """Given a smiles string (ex. C1CCCCC1), calculate and return the LogP"""
            mol = Chem.MolFromSmiles(smiles_string)
            return MolLogP(mol)

        def calc_NumHDonors(smiles_string):
            """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHDonors"""
            mol = Chem.MolFromSmiles(smiles_string)
            return NumHDonors(mol)

        def calc_NumHAcceptors(smiles_string):
            """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHAcceptors"""
            mol = Chem.MolFromSmiles(smiles_string)
            return NumHAcceptors(mol)


        # Copy the dataset so any changes are not applied to the original cached version
        df = download_dataset().copy()
        df["MW"] = df.apply(lambda x: calc_mw(x["smiles"]), axis=1)
        df["LogP"] = df.apply(lambda x: calc_logp(x["smiles"]), axis=1)
        df["NumHDonors"] = df.apply(lambda x: calc_NumHDonors(x["smiles"]), axis=1)
        df["NumHAcceptors"] = df.apply(lambda x: calc_NumHAcceptors(x["smiles"]), axis=1)


        # Sidebar panel
        st.sidebar.header('Set parameters')
        st.sidebar.write('*Note: Display compounds having values less than the following thresholds*')
        weight_cutoff = st.sidebar.slider(
            label="Molecular weight",
            min_value=0,
            max_value=1000,
            value=500,
            step=10,
        )
        logp_cutoff = st.sidebar.slider(
            label="LogP",
            min_value=-10,
            max_value=10,
            value=5,
            step=1,
        )
        NumHDonors_cutoff = st.sidebar.slider(
            label="NumHDonors",
            min_value=0,
            max_value=15,
            value=5,
            step=1,
        )
        NumHAcceptors_cutoff = st.sidebar.slider(
            label="NumHAcceptors",
            min_value=0,
            max_value=20,
            value=10,
            step=1,
        )

        df_result = df[df["MW"] < weight_cutoff]
        df_result2 = df_result[df_result["LogP"] < logp_cutoff]
        df_result3 = df_result2[df_result2["NumHDonors"] < NumHDonors_cutoff]
        df_result4 = df_result3[df_result3["NumHAcceptors"] < NumHAcceptors_cutoff]

        st.write(df_result4.shape)
        st.write(df_result4)


        raw_html = mols2grid.display(df_result4,
                                    #subset=["Name", "img"],
                                    subset=["img", "Name", "MW", "LogP", "NumHDonors", "NumHAcceptors"],
                                    rename={"smiles": "SMILES", "generic_name": "Name"})._repr_html_()
        components.html(raw_html, width=900, height=1100, scrolling=False)

    elif choice == "Humsafar":
        st.subheader("Embedded Browser")
        process = subprocess.Popen(['python','Humsafar.py'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        components.iframe("process", width=800, height=600)


if __name__ == '__main__':
    main()