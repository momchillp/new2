function reverseComplement() {
  const rawInput = document.getElementById("dna-input").value;

  // Remove FASTA headers (lines starting with ">") and join all lines
  const cleaned = rawInput
    .split("\n")
    .filter(line => !line.startsWith(">"))
    .join("")
    .toUpperCase()
    .replace(/[^ATCG]/g, ""); // Remove non-DNA characters

  if (cleaned.length === 0) {
    document.getElementById("reverse-output").textContent =
      "Please enter a valid DNA sequence (A, T, C, G only).";
    return;
  }

  const complementMap = { A: "T", T: "A", C: "G", G: "C" };

  const reverseComplement = cleaned
    .split("")
    .reverse()
    .map(base => complementMap[base] || "?")
    .join("");

  document.getElementById("reverse-output").textContent = reverseComplement;
}
