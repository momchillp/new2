// Parse FASTA or raw DNA sequence input to uppercase DNA string
function parseFASTA(input) {
  const lines = input.trim().split(/\r?\n/);
  if (lines[0].startsWith(">")) lines.shift();
  return lines.join("").replace(/[^ACGTacgt]/g, "").toUpperCase();
}

// Calculate GC content (%) of a DNA sequence
function calculateGCContent(seq) {
  const gcCount = (seq.match(/[GC]/gi) || []).length;
  return (gcCount / seq.length) * 100;
}

// Calculate melting temperature (Tm)
function calculateTm(seq) {
  const gcCount = (seq.match(/[GC]/gi) || []).length;
  const atCount = (seq.match(/[AT]/gi) || []).length;
  if (seq.length < 14) return 4 * gcCount + 2 * atCount;
  return 64.9 + 41 * (gcCount - 16.4) / seq.length;
}

// Get reverse complement of DNA sequence
function reverseComplement(seq) {
  return seq.split("").reverse().map(base => ({
    "A": "T", "T": "A", "C": "G", "G": "C"
  }[base] || base)).join("");
}

// Primer search with boundary checks and async processing
function findPrimersAsync(dna, regionStart, regionEnd, desiredTm, maxDiff, primerPairs, minProductLength, onProgress, onComplete) {
  const results = [];
  const maxPrimerLength = 30;
  const minPrimerLength = 18;

  const forwardSearchStart = Math.max(0, regionStart - 2000);
  const forwardSearchEnd = regionStart - minPrimerLength;
  const reverseSearchStart = regionEnd + minPrimerLength;
  const reverseSearchEnd = Math.min(dna.length - minPrimerLength, regionEnd + 2000);

  let fwdStart = forwardSearchStart;
  let fwdLen = minPrimerLength;
  let revStart = reverseSearchStart;
  let revLen = minPrimerLength;

  function next() {
    while (fwdStart <= forwardSearchEnd) {
      while (fwdLen <= maxPrimerLength) {
        if (fwdStart + fwdLen > dna.length) break;
        const fwdPrimer = dna.slice(fwdStart, fwdStart + fwdLen);
        const lastBaseFwd = fwdPrimer[fwdPrimer.length - 1];
        if (lastBaseFwd !== 'G' && lastBaseFwd !== 'C') {
        fwdLen++;
        continue;
        }
        const fwdTm = calculateTm(fwdPrimer);
        const fwdGC = calculateGCContent(fwdPrimer);
        if (fwdGC < 40 || fwdGC > 60 || Math.abs(fwdTm - desiredTm) > maxDiff) {
          fwdLen++;
          continue;
        }

        while (revStart <= reverseSearchEnd) {
          while (revLen <= maxPrimerLength) {
            if (revStart + revLen > dna.length) break;
            const revSeq = dna.slice(revStart, revStart + revLen);
            const revPrimer = reverseComplement(revSeq);
            const lastBaseRev = revPrimer[revPrimer.length - 1];
            if (lastBaseRev !== 'G' && lastBaseRev !== 'C') {
            revLen++;
            continue;
            }
            const revTm = calculateTm(revPrimer);
            const revGC = calculateGCContent(revPrimer);
            if (revGC < 40 || revGC > 60 || Math.abs(revTm - desiredTm) > maxDiff) {
              revLen++;
              continue;
            }

            const productLength = (revStart + revLen) - fwdStart;
            if ((fwdStart + fwdLen) > regionStart || revStart < regionEnd) {
            revLen++;
            continue;
            }

            if (fwdStart + fwdLen < regionStart && revStart > regionEnd) {
              results.push({
                fwdPrimer, revPrimer,
                fwdTm: fwdTm.toFixed(2), revTm: revTm.toFixed(2),
                fwdGC: fwdGC.toFixed(2), revGC: revGC.toFixed(2),
                fwdStart, revStart, productLength,
                fwdLen, revLen
              });
              if (results.length >= primerPairs) return finalize();
            }
            revLen++;
          }
          revLen = minPrimerLength;
          revStart++;
        }

        revStart = reverseSearchStart;
        revLen = minPrimerLength;
        fwdLen++;
      }
      fwdLen = minPrimerLength;
      fwdStart++;
      if (fwdStart % 100 === 0) {
        onProgress?.(results.length);
        return setTimeout(next, 0);
      }
    }
    finalize();
  }

  function finalize() {
    const unique = [];
    const seen = new Set();
    for (const r of results.sort((a, b) => a.productLength - b.productLength)) {
      const key = `${r.fwdStart}-${r.revStart}`;
      if (!seen.has(key)) {
        seen.add(key);
        unique.push(r);
      }
      if (unique.length >= primerPairs) break;
    }
    onComplete(unique);
  }

  next();
}

// Render primers to table
function renderResults(results) {
  const output = document.getElementById("primer-results");
  output.innerHTML = "";
  if (results.length === 0) return output.textContent = "No primers found.";

  const table = document.createElement("table");
  table.classList.add("primer-results-table");
  const header = document.createElement("tr");
  ["#", "Fwd Start", "Rev Start", "Fwd Primer", "Rev Primer", "Fwd Len", "Rev Len", "Fwd Tm", "Rev Tm", "Fwd GC", "Rev GC", "PCR Len"]
    .forEach(text => {
      const th = document.createElement("th");
      th.textContent = text;
      header.appendChild(th);
    });
  table.appendChild(header);

  results.forEach((r, i) => {
    const row = document.createElement("tr");
    [i + 1, r.fwdStart + 1, r.revStart + 1, r.fwdPrimer, r.revPrimer, r.fwdLen, r.revLen, r.fwdTm, r.revTm, r.fwdGC, r.revGC, r.productLength]
      .forEach(val => {
        const td = document.createElement("td");
        td.textContent = val;
        row.appendChild(td);
      });
    table.appendChild(row);
  });
  output.appendChild(table);
}

let zoom = 1;       // initial zoom factor
let pan = 0;        // initial pan offset

function drawPrimerMap(results, sequenceLength, regionStart = 0, regionEnd = 0) {
  const canvas = document.getElementById("sequence-map");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");

  const laneHeight = 30;
  const margin = 50;
  canvas.width = 1000;
  canvas.height = laneHeight * results.length + margin;

  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.font = "12px Arial";
  ctx.lineWidth = 2;

  // Helper to convert sequence coordinate to canvas coordinate
  const xScale = sequenceLength / (canvas.width / zoom);
  const toCanvasX = (pos) => ((pos - pan) / xScale);

  // Draw coordinate ruler
  ctx.fillStyle = "#000";
  ctx.beginPath();
  ctx.moveTo(0, 15);
  ctx.lineTo(canvas.width, 15);
  ctx.stroke();

  for (let bp = pan; bp <= pan + canvas.width * xScale; bp += 500) {
    const x = toCanvasX(bp);
    ctx.beginPath();
    ctx.moveTo(x, 10);
    ctx.lineTo(x, 20);
    ctx.stroke();
    ctx.fillText(bp.toString(), x + 2, 30);
  }

  // Highlight target region
  if (regionEnd > regionStart) {
    ctx.fillStyle = "rgba(255, 215, 0, 0.3)";
    const xStart = toCanvasX(regionStart);
    const xEnd = toCanvasX(regionEnd);
    ctx.fillRect(xStart, margin - 10, xEnd - xStart, canvas.height - margin);
  }

  // Draw each primer pair
  results.forEach((res, i) => {
    const y = i * laneHeight + margin;

    // Forward primer (green)
    const fx = toCanvasX(res.fwdStart);
    const fw = res.fwdLen / xScale;
    ctx.fillStyle = "#3aa95f";
    ctx.fillRect(fx, y, fw, 8);
    drawArrow(ctx, fx, y, fw, "forward");
    ctx.fillStyle = "#000";
    ctx.fillText(`Fwd(${res.fwdStart + 1})`, fx, y - 2);

    // Reverse primer (red)
    const rx = toCanvasX(res.revStart);
    const rw = res.revLen / xScale;
    ctx.fillStyle = "#ff4444";
    ctx.fillRect(rx, y + 10, rw, 8);
    drawArrow(ctx, rx, y + 10, rw, "reverse");
    ctx.fillStyle = "#000";
    ctx.fillText(`Rev(${res.revStart + 1})`, rx, y + 28);

    // PCR product (blue line)
    ctx.strokeStyle = "#0000ff";
    ctx.beginPath();
    ctx.moveTo(fx, y + 5);
    ctx.lineTo(toCanvasX(res.revStart + res.revLen), y + 5);
    ctx.stroke();
  });
}

// Draw directional arrow on primers
function drawArrow(ctx, x, y, width, direction) {
  ctx.fillStyle = direction === "forward" ? "#2e7d32" : "#b71c1c";
  const headSize = 5;
  ctx.beginPath();
  if (direction === "forward") {
    ctx.moveTo(x + width, y);
    ctx.lineTo(x + width, y + 8);
    ctx.lineTo(x + width + headSize, y + 4);
  } else {
    ctx.moveTo(x, y);
    ctx.lineTo(x, y + 8);
    ctx.lineTo(x - headSize, y + 4);
  }
  ctx.closePath();
  ctx.fill();
}

// Optional: Keyboard zoom and pan controls
document.addEventListener("keydown", (e) => {
  if (e.key === "+") zoom = Math.min(zoom * 1.2, 10);
  if (e.key === "-") zoom = Math.max(zoom / 1.2, 0.2);
  if (e.key === "ArrowRight") pan += 500;
  if (e.key === "ArrowLeft") pan -= 500;
  if (zoom < 1) pan = Math.min(pan, sequenceLength - (canvas.width * (1 / zoom)));
  if (pan < 0) pan = 0;

  // Redraw after zoom/pan (you must store the latest result set and call again)
  if (window.lastResults && window.lastSequenceLength) {
    drawPrimerMap(window.lastResults, window.lastSequenceLength, window.lastRegionStart, window.lastRegionEnd);
  }
});// Main event

document.getElementById("run-button").addEventListener("click", () => {
  const seqInput = document.getElementById("sequence-input").value;
  const sequence = parseFASTA(seqInput);
  if (!sequence) return alert("Enter a valid DNA sequence.");

  const desiredTm = parseFloat(document.getElementById("tm-input").value) || 60;
  const maxDiff = parseFloat(document.getElementById("tm-diff").value) || 3;
  const minProductLength = parseInt(document.getElementById("min-length").value) || 300;
  const primerPairs = parseInt(document.getElementById("num-pairs").value) || 5;

  let regionStart = parseInt(document.getElementById("region-start").value) - 1;
  let regionEnd = parseInt(document.getElementById("region-end").value) - 1;

  // Fallback to a centered 300 bp region if input is invalid or empty
  if (
  isNaN(regionStart) || isNaN(regionEnd) ||
  regionStart < 0 || regionEnd < 0 ||
  regionStart >= regionEnd || regionEnd >= sequence.length
) {
  const middle = Math.floor(sequence.length / 2);
  regionStart = middle - 150;
  regionEnd = middle + 150;
}

  const output = document.getElementById("primer-results");
  output.textContent = "Finding primers...";

  findPrimersAsync(sequence, regionStart, regionEnd, desiredTm, maxDiff, primerPairs, minProductLength,
    (progress) => output.textContent = `Searching... ${progress} found so far...`,
    (results) => {
      renderResults(results);
      drawPrimerMap(results, sequence.length);
    }
  );
});





