// ══════════════════════════════════════════════
// PSYCHROMETRIC FUNCTIONS (P = 101325 Pa)
// ══════════════════════════════════════════════
const P = 101325;

function pSat(T) {
  return 611.657 * Math.exp(17.2694 * T / (T + 238.29));
}

function wFromTHR(T, HR) {
  const ps = pSat(T);
  const pv = (HR / 100) * ps;
  return 1000 * 0.62198 * pv / (P - pv); // g/kg
}

function HRFromTW(T, w) {
  const pv = (w / 1000) * P / (0.62198 + w / 1000);
  const ps = pSat(T);
  return Math.min(100, (pv / ps) * 100);
}

function enthalpie(T, w) {
  return 1.006 * T + (w / 1000) * (2501 + 1.805 * T); // kJ/kg
}

function tRosee(T, HR) {
  const alpha = 17.2694, beta = 238.29;
  const gamma = Math.log(HR / 100) + alpha * T / (beta + T);
  return beta * gamma / (alpha - gamma);
}

function wFromTThum(T, Thum) {
  const ws = wFromTHR(Thum, 100);
  return ws - 6.6e-4 * (1 + 1.15e-3 * Thum) * (T - Thum) * (P / 100) / 1000;
}

function TFromHW(h, w) {
  return (h - 2501 * w / 1000) / (1.006 + 1.805 * w / 1000);
}

function wFromTH(T, h) {
  return 1000 * (h - 1.006 * T) / (2501 + 1.805 * T);
}

// ══════════════════════════════════════════════
// SOLVE STATE POINT
// ══════════════════════════════════════════════
function solvePoint() {
  const T    = parseFloat(document.getElementById('inp_T').value);
  const HR   = parseFloat(document.getElementById('inp_HR').value);
  const W    = parseFloat(document.getElementById('inp_W').value);
  const H    = parseFloat(document.getElementById('inp_H').value);
  const Thum = parseFloat(document.getElementById('inp_Thum').value);

  const hT = !isNaN(T), hHR = !isNaN(HR), hW = !isNaN(W), hH = !isNaN(H), hTh = !isNaN(Thum);

  let res = null;

  if (hT && hHR) {
    const w = wFromTHR(T, HR);
    if (w < 0) return null;
    res = { T, HR, w, h: enthalpie(T, w), Tr: tRosee(T, HR) };
  } else if (hT && hW) {
    const HR2 = HRFromTW(T, W);
    res = { T, HR: HR2, w: W, h: enthalpie(T, W), Tr: tRosee(T, HR2) };
  } else if (hT && hH) {
    const w = wFromTH(T, H);
    if (w < 0) return null;
    const HR2 = HRFromTW(T, w);
    res = { T, HR: HR2, w, h: H, Tr: tRosee(T, HR2) };
  } else if (hT && hTh) {
    const w = wFromTThum(T, Thum);
    if (w < 0) return null;
    const HR2 = HRFromTW(T, w);
    res = { T, HR: HR2, w, h: enthalpie(T, w), Tr: tRosee(T, HR2) };
  } else if (hH && hW) {
    const T2 = TFromHW(H, W);
    const HR2 = HRFromTW(T2, W);
    res = { T: T2, HR: HR2, w: W, h: H, Tr: tRosee(T2, HR2) };
  } else if (hW && hHR) {
    let Tg = 20;
    for (let i = 0; i < 100; i++) {
      const wt = wFromTHR(Tg, HR);
      if (Math.abs(wt - W) < 1e-6) break;
      const dwdT = (wFromTHR(Tg + 0.01, HR) - wt) / 0.01;
      Tg -= (wt - W) / dwdT;
    }
    res = { T: Tg, HR, w: W, h: enthalpie(Tg, W), Tr: tRosee(Tg, HR) };
  } else if (hH && hHR) {
    let Tg = 20;
    for (let i = 0; i < 100; i++) {
      const wt = wFromTHR(Tg, HR);
      const ht = enthalpie(Tg, wt);
      if (Math.abs(ht - H) < 1e-4) break;
      const wt2 = wFromTHR(Tg + 0.01, HR);
      const ht2 = enthalpie(Tg + 0.01, wt2);
      Tg -= (ht - H) / ((ht2 - ht) / 0.01);
    }
    const wt = wFromTHR(Tg, HR);
    res = { T: Tg, HR, w: wt, h: H, Tr: tRosee(Tg, HR) };
  } else if (hTh && hHR) {
    let Tg = Thum + 5;
    for (let i = 0; i < 100; i++) {
      const w2 = wFromTThum(Tg, Thum);
      const HR2 = HRFromTW(Tg, w2);
      if (Math.abs(HR2 - HR) < 1e-4) break;
      const w3 = wFromTThum(Tg + 0.01, Thum);
      const HR3 = HRFromTW(Tg + 0.01, w3);
      Tg -= (HR2 - HR) / ((HR3 - HR2) / 0.01);
    }
    const w2 = wFromTThum(Tg, Thum);
    res = { T: Tg, HR, w: w2, h: enthalpie(Tg, w2), Tr: tRosee(Tg, HR) };
  } else if (hTh && hW) {
    const ws = wFromTHR(Thum, 100);
    const A = 6.6e-4 * (P / 100) / 1000;
    const T2 = Thum + (ws - W) / A;
    const HR2 = HRFromTW(T2, W);
    res = { T: T2, HR: HR2, w: W, h: enthalpie(T2, W), Tr: tRosee(T2, HR2) };
  } else if (hTh && hH) {
    let Tg = 25;
    for (let i = 0; i < 100; i++) {
      const w2 = wFromTThum(Tg, Thum);
      const h2 = enthalpie(Tg, w2);
      if (Math.abs(h2 - H) < 1e-4) break;
      const w3 = wFromTThum(Tg + 0.01, Thum);
      const h3 = enthalpie(Tg + 0.01, w3);
      Tg -= (h2 - H) / ((h3 - h2) / 0.01);
    }
    const w2 = wFromTThum(Tg, Thum);
    const HR2 = HRFromTW(Tg, w2);
    res = { T: Tg, HR: HR2, w: w2, h: H, Tr: tRosee(Tg, HR2) };
  }

  if (!res) return null;
  if (res.HR < 0 || res.HR > 100.5 || res.w < 0 || res.w > 60 || res.T < -20 || res.T > 80) return null;
  return res;
}

// ══════════════════════════════════════════════
// AUTO-COMPUTE ON INPUT
// ══════════════════════════════════════════════
document.querySelectorAll('input[type=number]').forEach(inp => {
  inp.addEventListener('input', () => {
    const res = solvePoint();
    const cv  = document.getElementById('computed_vals');
    const err = document.getElementById('error_msg');
    if (res) {
      cv.style.display = 'flex';
      err.style.display = 'none';
      document.getElementById('cv_T').textContent  = res.T.toFixed(1)  + ' °C';
      document.getElementById('cv_HR').textContent = res.HR.toFixed(1) + ' %';
      document.getElementById('cv_W').textContent  = res.w.toFixed(2)  + ' g/kg';
      document.getElementById('cv_H').textContent  = res.h.toFixed(1)  + ' kJ/kg';
      document.getElementById('cv_Tr').textContent = res.Tr.toFixed(1) + ' °C';
    } else {
      cv.style.display = 'none';
      const count = ['inp_T','inp_HR','inp_W','inp_H','inp_Thum']
        .filter(id => document.getElementById(id).value !== '').length;
      err.style.display = count >= 2 ? 'block' : 'none';
    }
  });
});

// ══════════════════════════════════════════════
// POINTS STORAGE & STATE
// ══════════════════════════════════════════════
const COLORS = ['#58a6ff','#3fb950','#f78166','#d2a8ff','#ffa657','#79c0ff','#56d364','#ff7b72'];
let selectedColor = COLORS[0];
let points = [];
let pointCounter = 0;
let connectPoints = false;
let highlightedId = null;

// Build color dots
const colorRow = document.getElementById('color_row');
COLORS.forEach((c, i) => {
  const d = document.createElement('div');
  d.className = 'color-dot' + (i === 0 ? ' active' : '');
  d.style.background = c;
  d.onclick = () => {
    selectedColor = c;
    colorRow.querySelectorAll('.color-dot').forEach(x => x.classList.remove('active'));
    d.classList.add('active');
  };
  colorRow.appendChild(d);
});

function addPoint() {
  const res = solvePoint();
  if (!res) {
    document.getElementById('error_msg').style.display = 'block';
    return;
  }
  const lbl = document.getElementById('inp_label').value.trim() || ('P' + (++pointCounter));
  points.push({ ...res, color: selectedColor, label: lbl, id: Date.now() });
  renderPointsList();
  drawDiagram();
  clearInputs();
}

function clearInputs() {
  ['inp_T','inp_HR','inp_W','inp_H','inp_Thum'].forEach(id => {
    document.getElementById(id).value = '';
  });
  document.getElementById('inp_label').value = '';
  document.getElementById('computed_vals').style.display = 'none';
  document.getElementById('error_msg').style.display = 'none';
}

function removePoint(id) {
  points = points.filter(p => p.id !== id);
  renderPointsList();
  drawDiagram();
}

function renderPointsList() {
  const list = document.getElementById('points_list');
  if (points.length === 0) {
    list.innerHTML = '<div class="no-points">Aucun point placé</div>';
    return;
  }
  list.innerHTML = points.map(p => `
    <div class="point-item" onclick="highlightPoint(${p.id})">
      <div class="dot" style="background:${p.color}"></div>
      <span class="pt-label">${p.label}</span>
      <span class="pt-coords">${p.T.toFixed(1)}°C / ${p.HR.toFixed(0)}%</span>
      <button class="del-btn" onclick="event.stopPropagation();removePoint(${p.id})">✕</button>
    </div>
  `).join('');
}

function highlightPoint(id) {
  highlightedId = id;
  drawDiagram();
  setTimeout(() => { highlightedId = null; drawDiagram(); }, 1200);
}

function toggleConnect() {
  connectPoints = !connectPoints;
  const wrap = document.getElementById('toggle_wrap');
  wrap.classList.toggle('on', connectPoints);
  document.getElementById('toggle_label').textContent = connectPoints ? 'Points reliés' : 'Relier les points';
  drawDiagram();
}

// ══════════════════════════════════════════════
// DIAGRAM DRAWING
// ══════════════════════════════════════════════
const canvas = document.getElementById('psychro_canvas');
const ctx = canvas.getContext('2d');

const T_MIN = -10, T_MAX = 50;
const W_MIN = 0,   W_MAX = 30;
const MARGIN = { top: 30, right: 60, bottom: 50, left: 60 };

function setupCanvas() {
  const area = document.querySelector('.canvas-area');
  const aw = area.clientWidth  - 16;
  const ah = area.clientHeight - 16;
  canvas.width  = Math.max(280, aw);
  canvas.height = Math.max(200, ah);
  drawDiagram();
}

function diagramWidth()  { return canvas.width  - MARGIN.left - MARGIN.right; }
function diagramHeight() { return canvas.height - MARGIN.top  - MARGIN.bottom; }

function toCanvasX(T) { return MARGIN.left + (T - T_MIN) / (T_MAX - T_MIN) * diagramWidth(); }
function toCanvasY(w) { return MARGIN.top  + (1 - (w - W_MIN) / (W_MAX - W_MIN)) * diagramHeight(); }

function drawDiagram() {
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  ctx.fillStyle = '#0d1117';
  ctx.fillRect(0, 0, W, H);

  const dw = diagramWidth(), dh = diagramHeight();

  // ── CLIP ──
  ctx.save();
  ctx.beginPath();
  ctx.rect(MARGIN.left, MARGIN.top, dw, dh);
  ctx.clip();

  // ── GRID ──
  ctx.lineWidth = 0.5;
  for (let T = T_MIN; T <= T_MAX; T += 5) {
    ctx.strokeStyle = (T % 10 === 0) ? 'rgba(88,166,255,0.18)' : 'rgba(88,166,255,0.07)';
    ctx.beginPath();
    ctx.moveTo(toCanvasX(T), MARGIN.top);
    ctx.lineTo(toCanvasX(T), MARGIN.top + dh);
    ctx.stroke();
  }
  for (let w = 0; w <= W_MAX; w += 2) {
    ctx.strokeStyle = (w % 10 === 0) ? 'rgba(88,166,255,0.18)' : 'rgba(88,166,255,0.07)';
    ctx.beginPath();
    ctx.moveTo(MARGIN.left, toCanvasY(w));
    ctx.lineTo(MARGIN.left + dw, toCanvasY(w));
    ctx.stroke();
  }

  // ── HR CURVES ──
  [10, 20, 30, 40, 50, 60, 70, 80, 90, 100].forEach(hr => {
    ctx.beginPath();
    ctx.strokeStyle = hr === 100 ? 'rgba(88,166,255,0.7)' : 'rgba(88,166,255,0.25)';
    ctx.lineWidth = hr === 100 ? 1.5 : 0.8;
    let first = true;
    for (let T = T_MIN; T <= T_MAX; T += 0.5) {
      const w = wFromTHR(T, hr);
      if (w < W_MIN || w > W_MAX) { first = true; continue; }
      const cx = toCanvasX(T), cy = toCanvasY(w);
      first ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
      first = false;
    }
    ctx.stroke();
  });

  // HR labels
  ctx.font = '10px IBM Plex Mono';
  ctx.fillStyle = 'rgba(88,166,255,0.6)';
  [10, 20, 30, 40, 50, 60, 70, 80, 90].forEach(hr => {
    let Tlbl = T_MAX - 2;
    for (let T = T_MIN; T <= T_MAX; T += 0.5) {
      if (wFromTHR(T, hr) > W_MAX * 0.92) { Tlbl = T - 1; break; }
    }
    const w = wFromTHR(Tlbl, hr);
    if (w >= W_MIN && w <= W_MAX) ctx.fillText(hr + '%', toCanvasX(Tlbl) + 3, toCanvasY(w) - 3);
  });
  ctx.fillStyle = 'rgba(88,166,255,0.85)';
  ctx.fillText('100%', toCanvasX(18) + 3, toCanvasY(wFromTHR(18, 100)) - 3);

  // ── ENTHALPY LINES ──
  ctx.lineWidth = 0.6;
  for (let h = -10; h <= 120; h += 10) {
    ctx.beginPath();
    ctx.strokeStyle = 'rgba(63,185,80,0.2)';
    let first = true;
    for (let T = T_MIN; T <= T_MAX; T += 1) {
      const w = wFromTH(T, h);
      if (w < W_MIN || w > W_MAX) { first = true; continue; }
      const cx = toCanvasX(T), cy = toCanvasY(w);
      first ? ctx.moveTo(cx, cy) : ctx.lineTo(cx, cy);
      first = false;
    }
    ctx.stroke();
  }

  ctx.restore(); // end clip

  // ── AXES ──
  ctx.strokeStyle = 'rgba(88,166,255,0.4)';
  ctx.lineWidth = 1;
  ctx.strokeRect(MARGIN.left, MARGIN.top, dw, dh);

  ctx.font = '11px IBM Plex Mono';
  ctx.fillStyle = 'rgba(230,237,243,0.7)';
  ctx.textAlign = 'center';
  for (let T = T_MIN; T <= T_MAX; T += 5) {
    const cx = toCanvasX(T);
    ctx.fillText(T + '°', cx, MARGIN.top + dh + 18);
    ctx.strokeStyle = 'rgba(88,166,255,0.3)';
    ctx.lineWidth = 0.5;
    ctx.beginPath(); ctx.moveTo(cx, MARGIN.top + dh); ctx.lineTo(cx, MARGIN.top + dh + 5); ctx.stroke();
  }

  ctx.textAlign = 'right';
  for (let w = 0; w <= W_MAX; w += 2) {
    const cy = toCanvasY(w);
    ctx.fillText(w, MARGIN.left - 8, cy + 4);
    ctx.lineWidth = 0.5;
    ctx.beginPath(); ctx.moveTo(MARGIN.left, cy); ctx.lineTo(MARGIN.left - 4, cy); ctx.stroke();
  }

  ctx.fillStyle = 'rgba(88,166,255,0.8)';
  ctx.font = '12px IBM Plex Mono';
  ctx.textAlign = 'center';
  ctx.fillText('Température sèche T [°C]', MARGIN.left + dw / 2, canvas.height - 8);

  ctx.save();
  ctx.translate(14, MARGIN.top + dh / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('Humidité spécifique w [g/kg]', 0, 0);
  ctx.restore();

  ctx.fillStyle = 'rgba(63,185,80,0.5)';
  ctx.font = '10px IBM Plex Mono';
  ctx.textAlign = 'left';
  ctx.fillText('— h [kJ/kg]', MARGIN.left + 4, MARGIN.top + 14);

  // ── CONNECTING LINES ──
  if (connectPoints && points.length >= 2) {
    ctx.save();
    ctx.beginPath();
    ctx.rect(MARGIN.left, MARGIN.top, dw, dh);
    ctx.clip();

    for (let i = 0; i < points.length - 1; i++) {
      const a = points[i], b = points[i + 1];
      const ax = toCanvasX(a.T), ay = toCanvasY(a.w);
      const bx = toCanvasX(b.T), by = toCanvasY(b.w);

      const grad = ctx.createLinearGradient(ax, ay, bx, by);
      grad.addColorStop(0, a.color + 'cc');
      grad.addColorStop(1, b.color + 'cc');

      ctx.beginPath();
      ctx.moveTo(ax, ay);
      ctx.lineTo(bx, by);
      ctx.strokeStyle = grad;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([6, 4]);
      ctx.stroke();
      ctx.setLineDash([]);

      const mx2 = (ax + bx) / 2, my2 = (ay + by) / 2;
      const angle = Math.atan2(by - ay, bx - ax);
      const aSize = 6;
      ctx.beginPath();
      ctx.moveTo(mx2 + Math.cos(angle) * aSize, my2 + Math.sin(angle) * aSize);
      ctx.lineTo(mx2 + Math.cos(angle + 2.4) * aSize, my2 + Math.sin(angle + 2.4) * aSize);
      ctx.lineTo(mx2 + Math.cos(angle - 2.4) * aSize, my2 + Math.sin(angle - 2.4) * aSize);
      ctx.closePath();
      ctx.fillStyle = grad;
      ctx.fill();
    }
    ctx.restore();
  }

  // ── PLACED POINTS ──
  ctx.save();
  ctx.beginPath();
  ctx.rect(MARGIN.left - 20, MARGIN.top - 20, dw + 40, dh + 40);
  ctx.clip();

  points.forEach(p => {
    const cx = toCanvasX(p.T);
    const cy = toCanvasY(p.w);
    const isHighlighted = (p.id === highlightedId);
    const r = isHighlighted ? 9 : 6;

    if (isHighlighted) {
      ctx.beginPath();
      ctx.arc(cx, cy, 16, 0, Math.PI * 2);
      ctx.fillStyle = p.color + '20';
      ctx.fill();
    }

    ctx.strokeStyle = p.color + '60';
    ctx.lineWidth = 0.5;
    ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(MARGIN.left, cy); ctx.lineTo(cx, cy); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(cx, MARGIN.top + dh); ctx.lineTo(cx, cy); ctx.stroke();
    ctx.setLineDash([]);

    ctx.beginPath();
    ctx.arc(cx, cy, r, 0, Math.PI * 2);
    ctx.fillStyle = p.color;
    ctx.fill();
    ctx.strokeStyle = '#0d1117';
    ctx.lineWidth = 1.5;
    ctx.stroke();

    ctx.fillStyle = p.color;
    ctx.font = 'bold 11px IBM Plex Mono';
    ctx.textAlign = 'left';
    ctx.fillText(p.label, cx + r + 4, cy - 4);
  });

  ctx.restore();
}

// ══════════════════════════════════════════════
// TOOLTIP ON HOVER
// ══════════════════════════════════════════════
const tooltip = document.getElementById('tooltip');

canvas.addEventListener('mousemove', e => {
  const rect = canvas.getBoundingClientRect();
  const mx = e.clientX - rect.left;
  const my = e.clientY - rect.top;

  let found = null;
  points.forEach(p => {
    const cx = toCanvasX(p.T), cy = toCanvasY(p.w);
    if (Math.hypot(mx - cx, my - cy) < 10) found = p;
  });

  if (found) {
    tooltip.style.display = 'block';
    tooltip.style.left = (e.clientX + 12) + 'px';
    tooltip.style.top  = (e.clientY - 10) + 'px';
    tooltip.innerHTML =
      `<b style="color:${found.color}">${found.label}</b><br>` +
      `T = ${found.T.toFixed(1)} °C<br>` +
      `φ = ${found.HR.toFixed(1)} %<br>` +
      `w = ${found.w.toFixed(2)} g/kg<br>` +
      `h = ${found.h.toFixed(1)} kJ/kg<br>` +
      `T<sub>rosée</sub> = ${found.Tr.toFixed(1)} °C`;
  } else {
    tooltip.style.display = 'none';
    const T = T_MIN + (mx - MARGIN.left) / diagramWidth() * (T_MAX - T_MIN);
    const w = W_MAX - (my - MARGIN.top) / diagramHeight() * (W_MAX - W_MIN);
    if (T >= T_MIN && T <= T_MAX && w >= W_MIN && w <= W_MAX) {
      const HR = HRFromTW(T, w);
      if (HR >= 0 && HR <= 100) {
        const h = enthalpie(T, w);
        tooltip.style.display = 'block';
        tooltip.style.left = (e.clientX + 12) + 'px';
        tooltip.style.top  = (e.clientY - 10) + 'px';
        tooltip.innerHTML =
          `T = ${T.toFixed(1)} °C | φ = ${HR.toFixed(0)} %<br>` +
          `w = ${w.toFixed(2)} g/kg | h = ${h.toFixed(1)} kJ/kg`;
      }
    }
  }
});

canvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });

// ══════════════════════════════════════════════
// DRAWER (mobile)
// ══════════════════════════════════════════════
function toggleDrawer() {
  const panel   = document.querySelector('.panel');
  const overlay = document.getElementById('overlay');
  const btn     = document.getElementById('drawer_btn');
  const isOpen  = panel.classList.toggle('open');
  overlay.classList.toggle('visible', isOpen);
  btn.textContent = isOpen ? '✕' : '⊕';
}

function closeDrawer() {
  const panel   = document.querySelector('.panel');
  const overlay = document.getElementById('overlay');
  const btn     = document.getElementById('drawer_btn');
  panel.classList.remove('open');
  overlay.classList.remove('visible');
  btn.textContent = '⊕';
}

// ══════════════════════════════════════════════
// INIT
// ══════════════════════════════════════════════
window.addEventListener('resize', setupCanvas);
setupCanvas();
