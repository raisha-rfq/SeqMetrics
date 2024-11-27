import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.RNASequence;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map;

public class BioJavaGUI extends JFrame {

    private JTextArea outputArea;
    private JFileChooser fileChooser;

    public BioJavaGUI() {
        setTitle("SeqMetrics");
        setSize(600, 500);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        // Create GUI Components
        outputArea = new JTextArea();
        outputArea.setEditable(false);
        JScrollPane scrollPane = new JScrollPane(outputArea);
        
        JButton loadButton = new JButton("Load FASTA File");
        JButton analyzeButton = new JButton("Analyze Sequences");

        JPanel panel = new JPanel();
        panel.setLayout(new FlowLayout());
        panel.add(loadButton);
        panel.add(analyzeButton);

        add(scrollPane, BorderLayout.CENTER);
        add(panel, BorderLayout.SOUTH);

        // Set up file chooser
        fileChooser = new JFileChooser();
        fileChooser.setCurrentDirectory(new File(System.getProperty("user.home")));

        // Load Button Action
        loadButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int result = fileChooser.showOpenDialog(null);
                if (result == JFileChooser.APPROVE_OPTION) {
                    File selectedFile = fileChooser.getSelectedFile();
                    outputArea.append("Selected file: " + selectedFile.getAbsolutePath() + "\n");
                }
            }
        });

        // Analyze Button Action
        analyzeButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                File selectedFile = fileChooser.getSelectedFile();
                if (selectedFile == null) {
                    outputArea.append("Please select a file first!\n");
                } else {
                    analyzeSequences(selectedFile);
                }
            }
        });
    }

    // Analyze sequences
    private void analyzeSequences(File file) {
        try {
            outputArea.append("Analyzing sequences...\n");

            // Step 1: Load FASTA file
            Map<String, DNASequence> dnaSequences = FastaReaderHelper.readFastaDNASequence(file);

            // Validate sequences
            if (!validateParameters((LinkedHashMap<String, DNASequence>) dnaSequences)) {
                outputArea.append("Error: Invalid sequences.\n");
                return;
            }

            // Iterate through sequences
            for (String key : dnaSequences.keySet()) {
                DNASequence dnaSequence = dnaSequences.get(key);
                outputArea.append("\nProcessing sequence: " + key + "\n");
                outputArea.append("Sequence: " + dnaSequence.getSequenceAsString().toUpperCase() + "\n");

                // Validate sequence format
                if (!validateSequenceFormat(dnaSequence)) {
                    outputArea.append("Error: Invalid sequence format for sequence: " + key + "\n");
                    continue;
                }

                // Step 2: Translate DNA to Protein
                ProteinSequence proteinSequence = translateDNAToProtein(dnaSequence);
                outputArea.append("Protein Sequence: " + proteinSequence.getSequenceAsString() + "\n");

                // Step 3: Calculate molecular weight
                float molecularWeight = calculateMolecularWeight(proteinSequence);
                outputArea.append("Molecular Weight: " + molecularWeight + " Da\n");

                // Step 4: Detect palindromes
                if (findPalindromes(dnaSequence)) {
                    outputArea.append("Palindrome detected.\n");
                } else {
                    outputArea.append("No palindrome detected.\n");
                }
            }

        } catch (Exception ex) {
            ex.printStackTrace();
            outputArea.append("Error: " + ex.getMessage() + "\n");
        }
    }

    // Function to check if the file loaded properly
    private static boolean validateParameters(LinkedHashMap<String, DNASequence> sequences) {
        return sequences != null && !sequences.isEmpty();
    }

    // Function to ensure it's a valid DNA sequence
    private static boolean validateSequenceFormat(DNASequence sequence) {
        String validBases = "ATGCNRYKMSWBDHV";
        String seq = sequence.getSequenceAsString().toUpperCase().replaceAll("[^ATGCNRYKMSWBDHV]", "");
        for (char base : seq.toCharArray()) {
            if (!validBases.contains(String.valueOf(base))) {
                return false;
            }
        }
        return true;
    }

    // Function to translate DNA to Protein Sequence
    private static ProteinSequence translateDNAToProtein(DNASequence dnaSequence) throws CompoundNotFoundException {
        String rnaSequenceStr = dnaSequence.getSequenceAsString().replace('T', 'U');
        RNASequence rnaSequence = new RNASequence(rnaSequenceStr);

        StringBuilder proteinStr = new StringBuilder();
        String rnaSeq = rnaSequence.getSequenceAsString();
        String[] codons = {"UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AUU", "AUC", "AUA", "AUG",
                            "GUU", "GUC", "GUA", "GUG", "UCU", "UCC", "UCA", "UCG", "CCU", "CCC", "CCA", "CCG",
                            "ACU", "ACC", "ACA", "ACG", "GCU", "GCC", "GCA", "GCG", "UAU", "UAC", "UAA", "UAG",
                            "CAU", "CAC", "CAA", "CAG", "AAU", "AAC", "AAA", "AAG", "GAU", "GAC", "GAA", "GAG",
                            "UGU", "UGC", "UGA", "UGG"};

        String[] aminoAcids = {"F", "F", "L", "L", "L", "L", "I", "I", "I", "M", "V", "V", "V", "V", "S", "S", "S", 
                               "S", "P", "P", "P", "P", "T", "T", "T", "T", "A", "A", "A", "A", "Y", "Y", "*", "*", 
                               "H", "H", "Q", "Q", "N", "N", "K", "K", "D", "D", "E", "E", "C", "C", "*", "W"};

        for (int i = 0; i < rnaSeq.length() - 2; i += 3) {
            String codon = rnaSeq.substring(i, i + 3);
            boolean codonFound = false;
            for (int j = 0; j < codons.length; j++) {
                if (codons[j].equals(codon)) {
                    if (j < aminoAcids.length) {
                        proteinStr.append(aminoAcids[j]);
                    }
                    codonFound = true;
                    break;
                }
            }
            if (!codonFound) {
                proteinStr.append("X");
            }
        }
        return new ProteinSequence(proteinStr.toString());
    }

    // Function to calculate molecular weight
    private static float calculateMolecularWeight(ProteinSequence proteinSequence) {
        float totalMolecularWeight = 0.0f;
        for (AminoAcidCompound aminoAcid : proteinSequence) {
            Float molecularWeight = aminoAcid.getMolecularWeight();
            if (molecularWeight != null) {
                totalMolecularWeight += molecularWeight;
            }
        }
        return totalMolecularWeight;
    }

    // Function to find palindromes in the DNA sequence
    private static boolean findPalindromes(DNASequence sequence) {
        String seq = sequence.getSequenceAsString();
        String revComp = sequence.getReverseComplement().getSequenceAsString();
        return seq.equalsIgnoreCase(revComp);
    }

    public static void main(String[] args) {
        // Launch the GUI
        SwingUtilities.invokeLater(() -> {
            BioJavaGUI gui = new BioJavaGUI();
            gui.setVisible(true);
        });
    }
}
