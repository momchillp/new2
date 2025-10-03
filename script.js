function translateDNA() {
  const codonTable = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_',
    'TGC':'C','TGT':'C','TGA':'_','TGG':'W',
  };

  let input = document.getElementById("dna-input").value;

  // Clean input: remove FASTA headers, whitespace, and invalid characters
  let dna = input
    .toUpperCase()
    .split('\n')                         // Split into lines
    .filter(line => !line.startsWith('>')) // Remove FASTA headers
    .join('')                             // Combine lines
    .replace(/[^ATCG]/g, '');             // Remove invalid characters

  let protein = "";

  for (let i = 0; i < dna.length - 2; i += 3) {
    const codon = dna.substring(i, i + 3);
    protein += codonTable[codon] || '?';
  }

  document.getElementById("protein-output").textContent = protein;
}
