function parseDNAInput(nameInput, seqInput) {
  let name = nameInput.value.trim() || "Unnamed";
  let sequence = seqInput.value.replace(/[^ACGTacgt]/g, "").toUpperCase();
  return { name, sequence };
}

// Needleman-Wunsch Global Alignment (simplified)
function needlemanWunsch(seq1, seq2) {
  const gap = -1;
  const match = 1;
  const mismatch = -1;
  const m = seq1.length;
  const n = seq2.length;

  const score = Array.from({ length: m + 1 }, () => Array(n + 1).fill(0));
  for (let i = 0; i <= m; i++) score[i][0] = i * gap;
  for (let j = 0; j <= n; j++) score[0][j] = j * gap;

  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const matchScore = seq1[i - 1] === seq2[j - 1] ? match : mismatch;
      score[i][j] = Math.max(
        score[i - 1][j - 1] + matchScore,
        score[i - 1][j] + gap,
        score[i][j - 1] + gap
      );
    }
  }

  // Traceback not needed, just compute identity %
  let aligned = Math.max(m, n);
  let identical = (score[m][n] + aligned) / 2;
  return 1 - identical / aligned;
}

function computeDistanceMatrix(sequences) {
  const n = sequences.length;
  const matrix = Array.from({ length: n }, () => Array(n).fill(0));

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const dist = needlemanWunsch(sequences[i].sequence, sequences[j].sequence);
      matrix[i][j] = matrix[j][i] = dist;
    }
  }

  return matrix;
}

// UPGMA clustering
function buildUPGMATree(sequences) {
  let clusters = sequences.map((s, i) => ({
    name: s.name,
    sequence: s.sequence,
    children: [],
    indices: [i],
    distance: 0
  }));

  let matrix = computeDistanceMatrix(sequences);

  while (clusters.length > 1) {
    let minDist = Infinity;
    let pair = [0, 1];

    for (let i = 0; i < clusters.length; i++) {
      for (let j = i + 1; j < clusters.length; j++) {
        let d = 0;
        for (let a of clusters[i].indices) {
          for (let b of clusters[j].indices) {
            d += matrix[a][b];
          }
        }
        d /= clusters[i].indices.length * clusters[j].indices.length;
        if (d < minDist) {
          minDist = d;
          pair = [i, j];
        }
      }
    }

    let [i, j] = pair;
    let merged = {
      name: "",
      children: [clusters[i], clusters[j]],
      indices: [...clusters[i].indices, ...clusters[j].indices],
      distance: minDist / 2
    };

    if (i > j) [i, j] = [j, i];
    clusters.splice(j, 1);
    clusters.splice(i, 1);
    clusters.push(merged);
  }

  return clusters[0];
}

function renderTree(root) {
  const svg = document.getElementById("tree-svg");
  svg.innerHTML = "";

  const width = svg.clientWidth;
  const height = svg.clientHeight;

  let leaves = [];
  function collectLeaves(node) {
    if (!node.children || node.children.length === 0) {
      leaves.push(node);
    } else {
      node.children.forEach(collectLeaves);
    }
  }
  collectLeaves(root);

  leaves.forEach((leaf, i) => (leaf.y = (i + 1) * height / (leaves.length + 1)));

  function assignX(node) {
    if (!node.children || node.children.length === 0) {
      node.x = width/2;
    } else {
      node.children.forEach(assignX);
      node.x = Math.min(...node.children.map(c => c.x)) - 100;
      node.y = node.children.reduce((sum, c) => sum + c.y, 0) / node.children.length;
    }
  }
  assignX(root);

  function drawNode(node) {
    if (node.children && node.children.length > 0) {
      node.children.forEach(child => {
        const lineH = document.createElementNS("http://www.w3.org/2000/svg", "line");
        lineH.setAttribute("x1", node.x);
        lineH.setAttribute("y1", child.y);
        lineH.setAttribute("x2", child.x);
        lineH.setAttribute("y2", child.y);
        lineH.setAttribute("class", "branch-line");
        svg.appendChild(lineH);

        const lineV = document.createElementNS("http://www.w3.org/2000/svg", "line");
        lineV.setAttribute("x1", node.x);
        lineV.setAttribute("y1", node.y);
        lineV.setAttribute("x2", node.x);
        lineV.setAttribute("y2", child.y);
        lineV.setAttribute("class", "branch-line");
        svg.appendChild(lineV);

        const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
        label.setAttribute("x", (node.x + child.x) / 2);
        label.setAttribute("y", child.y - 5);
        label.setAttribute("class", "branch-label");
        label.textContent = node.distance.toFixed(2);
        svg.appendChild(label);

        drawNode(child);
      });
    } else {
      const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
      label.setAttribute("x", node.x + 5);
      label.setAttribute("y", node.y + 4);
      label.setAttribute("class", "leaf-label");
      label.textContent = node.name;
      svg.appendChild(label);
    }
  }

  drawNode(root);
}

function buildTree() {
  const boxes = document.querySelectorAll(".species-entry");
  const parsed = [];

  boxes.forEach(entry => {
    const nameInput = entry.querySelector(".species-name");
    const seqInput = entry.querySelector(".species-seq");
    const parsedSeq = parseDNAInput(nameInput, seqInput);
    if (parsedSeq.sequence.length >= 5) {
      parsed.push(parsedSeq);
    }
  });

  if (parsed.length < 2) {
    alert("Please enter at least two valid sequences.");
    return;
  }

  const tree = buildUPGMATree(parsed);
  renderTree(tree);
}

function addSpeciesBox() {
  const container = document.getElementById("species-inputs");

  const div = document.createElement("div");
  div.className = "species-entry";

  const name = document.createElement("input");
  name.type = "text";
  name.className = "species-name";
  name.placeholder = "Species Name";

  const seq = document.createElement("textarea");
  seq.className = "species-seq";
  seq.placeholder = "ATGC...";

  div.appendChild(name);
  div.appendChild(seq);
  container.appendChild(div);
}

function exportSVG() {
  const svg = document.getElementById("tree-svg");
  const clone = svg.cloneNode(true);

  // Inline critical styles
  function inlineStyles(element) {
    const computed = window.getComputedStyle(element);
    const tagName = element.tagName.toLowerCase();

    // Always set fill, stroke, etc.
    if (tagName === "line" || tagName === "path") {
      element.setAttribute("stroke", computed.stroke || "#000");
      element.setAttribute("stroke-width", computed.strokeWidth || "2");
    } else if (tagName === "text") {
      element.setAttribute("fill", computed.fill || "#000");
      element.setAttribute("font-size", computed.fontSize || "12px");
      element.setAttribute("font-family", computed.fontFamily || "Arial, sans-serif");
    }

    for (const child of element.children) {
      inlineStyles(child);
    }
  }

  inlineStyles(clone);

  // Serialize cloned SVG
  const serializer = new XMLSerializer();
  const svgString = serializer.serializeToString(clone);

  // Create blob and image
  const svgBlob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(svgBlob);
  const img = new Image();

  img.onload = function () {
    const canvas = document.createElement("canvas");
    canvas.width = svg.width.baseVal.value;
    canvas.height = svg.height.baseVal.value;
    const ctx = canvas.getContext("2d");

    // Force white background
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0);

    // Export as PNG
    const pngLink = document.createElement("a");
    pngLink.download = "phylogenetic_tree.png";
    pngLink.href = canvas.toDataURL("image/png");
    pngLink.click();

    URL.revokeObjectURL(url);
  };

  img.src = url;
}
