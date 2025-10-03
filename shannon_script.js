let speciesCounts = [0];

function adjustCount(index, delta) {
  const input = document.querySelector(`.species-row[data-index="${index}"] input`);
  const current = parseInt(input.value) || 0;
  const newVal = Math.max(0, current + delta);
  input.value = newVal;
  speciesCounts[index] = newVal;
}

function manualInput(index, value) {
  speciesCounts[index] = Math.max(0, parseInt(value) || 0);
}

function addSpecies() {
  const index = speciesCounts.length;
  speciesCounts.push(0);

  const row = document.createElement('div');
  row.className = 'species-row';
  row.setAttribute('data-index', index);
  row.innerHTML = `
    <span>Species ${index + 1}:</span>
    <div class="input-controls">
      <button onclick="adjustCount(${index}, -1)">−</button>
      <input type="number" min="0" value="0" oninput="manualInput(${index}, this.value)" />
      <button onclick="adjustCount(${index}, 1)">+</button>
    </div>
  `;
  document.getElementById('species-list').appendChild(row);
}

function calculateShannon() {
  const total = speciesCounts.reduce((sum, n) => sum + n, 0);

  if (total === 0) {
    document.getElementById('shannon-result').textContent = 'Please enter some counts.';
    document.getElementById('shannon-result1').textContent = '';
    document.getElementById('shannon-result2').textContent = '';
    return;
  }

   const H = -speciesCounts.reduce((sum, count) => {
    if (count === 0) return sum;
    const p = count / total;
    return sum + p * Math.log(p);
  }, 0);
  const result = H.toFixed(4);
  document.getElementById('shannon-result').textContent = `Shannon Diversity Index (H): ${result}`;

    const S = speciesCounts.filter(c => c > 0).length; // number of species with count > 0
  const result1 = (H / Math.log(S)).toFixed(4);
  document.getElementById('shannon-result1').textContent = `Pielou's Evenness Index (J): ${result1}`;

   document.getElementById('shannon-result2').textContent = `Total number of individuals: ${total}`;
}
