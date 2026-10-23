/* ===================================================================
   Lattice Pricer — moteur de calcul
   Binomial (Cox-Ross-Rubinstein), trinomial (Boyle / Kamrad-Ritchken)
   et Black-Scholes-Merton, options européennes et américaines,
   avec rendement de dividende continu q.
   =================================================================== */
const Engine = (() => {
  // Fonction de répartition normale, précision ~1e-15 (algorithme de West, 2005)
  function ncdf(x) {
    const a = Math.abs(x);
    let c;
    if (a > 37) c = 0;
    else {
      const e = Math.exp(-a * a / 2);
      if (a < 7.07106781186547) {
        let b = 3.52624965998911e-02 * a + 0.700383064443688;
        b = b * a + 6.37396220353165; b = b * a + 33.912866078383;
        b = b * a + 112.079291497871; b = b * a + 221.213596169931;
        b = b * a + 220.206867912376;
        c = e * b;
        b = 8.83883476483184e-02 * a + 1.75566716318264;
        b = b * a + 16.064177579207; b = b * a + 86.7807322029461;
        b = b * a + 296.564248779674; b = b * a + 637.333633378831;
        b = b * a + 793.826512519948; b = b * a + 440.413735824752;
        c = c / b;
      } else {
        let b = a + 0.65; b = a + 4 / b; b = a + 3 / b; b = a + 2 / b; b = a + 1 / b;
        c = e / b / 2.506628274631;
      }
    }
    return x > 0 ? 1 - c : c;
  }
  const npdf = x => Math.exp(-x * x / 2) / 2.5066282746310002;

  /* ---------- Black-Scholes-Merton (européenne) ---------- */
  function bs(p) {
    const { S, K, T, r, q, sigma, call } = p;
    const sT = sigma * Math.sqrt(T);
    const d1 = (Math.log(S / K) + (r - q + sigma * sigma / 2) * T) / sT, d2 = d1 - sT;
    const dq = Math.exp(-q * T), dr = Math.exp(-r * T);
    const price = call ? S * dq * ncdf(d1) - K * dr * ncdf(d2)
                       : K * dr * ncdf(-d2) - S * dq * ncdf(-d1);
    const delta = call ? dq * ncdf(d1) : dq * (ncdf(d1) - 1);
    const gamma = dq * npdf(d1) / (S * sT);
    const vega = S * dq * npdf(d1) * Math.sqrt(T);
    const theta = call
      ? -S * dq * npdf(d1) * sigma / (2 * Math.sqrt(T)) - r * K * dr * ncdf(d2) + q * S * dq * ncdf(d1)
      : -S * dq * npdf(d1) * sigma / (2 * Math.sqrt(T)) + r * K * dr * ncdf(-d2) - q * S * dq * ncdf(-d1);
    const rho = call ? K * T * dr * ncdf(d2) : -K * T * dr * ncdf(-d2);
    return { price, delta, gamma, theta, vega, rho, d1, d2 };
  }

  /* ---------- Arbre binomial CRR ---------- */
  // opts.boundary = true : renvoie la frontière d'exercice S*(t)
  function binomial(p, N, opts = {}) {
    const { S, K, T, r, q, sigma, call, american } = p;
    const dt = T / N, u = Math.exp(sigma * Math.sqrt(dt)), d = 1 / u;
    const disc = Math.exp(-r * dt), g = Math.exp((r - q) * dt);
    const pu = (g - d) / (u - d), pd = 1 - pu;
    if (!(pu > 0 && pu < 1)) return { price: NaN, error: 'arbitrage' };
    const V = new Float64Array(N + 1);
    const pw = new Float64Array(2 * N + 1);
    pw[N] = S; for (let m = 1; m <= N; m++) { pw[N + m] = pw[N + m - 1] * u; pw[N - m] = pw[N - m + 1] * d; }
    const sAt = (k, j) => pw[N + 2 * j - k];
    for (let j = 0; j <= N; j++) {
      const s = sAt(N, j); V[j] = call ? Math.max(s - K, 0) : Math.max(K - s, 0);
    }
    const bnd = opts.boundary ? new Float64Array(N + 1).fill(NaN) : null;
    let V1 = null, V2 = null;
    for (let k = N - 1; k >= 0; k--) {
      for (let j = 0; j <= k; j++) {
        let v = disc * (pu * V[j + 1] + pd * V[j]);
        if (american) {
          const s = sAt(k, j), ex = call ? s - K : K - s;
          if (ex > v) {
            v = ex;
            if (bnd) bnd[k] = call ? (isNaN(bnd[k]) ? s : Math.min(bnd[k], s)) : (isNaN(bnd[k]) ? s : Math.max(bnd[k], s));
          }
        }
        V[j] = v;
      }
      if (k === 2) V2 = [V[0], V[1], V[2]];
      if (k === 1) V1 = [V[0], V[1]];
    }
    const out = { price: V[0], u, d, pu, pm: 0, pd, disc };
    if (V1) out.delta = (V1[1] - V1[0]) / (S * u - S * d);
    if (V2) {
      const Suu = S * u * u, Sdd = S * d * d;
      out.gamma = ((V2[2] - V2[1]) / (Suu - S) - (V2[1] - V2[0]) / (S - Sdd)) / ((Suu - Sdd) / 2);
      out.theta = (V2[1] - V[0]) / (2 * dt);
    }
    if (bnd) out.boundary = bnd;
    return out;
  }

  /* ---------- Arbre trinomial (calibrage par lambda) ---------- */
  function trinomial(p, N, lambda = Math.sqrt(1.5), opts = {}) {
    const { S, K, T, r, q, sigma, call, american } = p;
    const dt = T / N, u = Math.exp(lambda * sigma * Math.sqrt(dt)), d = 1 / u;
    const disc = Math.exp(-r * dt), g = Math.exp((r - q) * dt);
    const pm = 1 - 1 / (lambda * lambda);
    const pu = (g - pm - (1 - pm) * d) / (u - d), pd = 1 - pm - pu;
    if (!(pu > 0 && pd > 0 && pm >= 0)) return { price: NaN, error: 'proba' };
    const M = 2 * N + 1, V = new Float64Array(M);
    const pw = new Float64Array(M);
    pw[N] = S; for (let m = 1; m <= N; m++) { pw[N + m] = pw[N + m - 1] * u; pw[N - m] = pw[N - m + 1] * d; }
    const sAt = (k, i) => pw[N + i - k]; // i = 0..2k
    for (let i = 0; i < M; i++) {
      const s = sAt(N, i); V[i] = call ? Math.max(s - K, 0) : Math.max(K - s, 0);
    }
    const bnd = opts.boundary ? new Float64Array(N + 1).fill(NaN) : null;
    let V1 = null;
    for (let k = N - 1; k >= 0; k--) {
      for (let i = 0; i <= 2 * k; i++) {
        let v = disc * (pu * V[i + 2] + pm * V[i + 1] + pd * V[i]);
        if (american) {
          const s = sAt(k, i), ex = call ? s - K : K - s;
          if (ex > v) {
            v = ex;
            if (bnd) bnd[k] = call ? (isNaN(bnd[k]) ? s : Math.min(bnd[k], s)) : (isNaN(bnd[k]) ? s : Math.max(bnd[k], s));
          }
        }
        V[i] = v;
      }
      if (k === 1) V1 = [V[0], V[1], V[2]];
    }
    const out = { price: V[0], u, d, pu, pm, pd, disc };
    if (V1) {
      const Su = S * u, Sd = S * d;
      out.delta = (V1[2] - V1[0]) / (Su - Sd);
      out.gamma = ((V1[2] - V1[1]) / (Su - S) - (V1[1] - V1[0]) / (S - Sd)) / ((Su - Sd) / 2);
      out.theta = (V1[1] - V[0]) / dt;
    }
    if (bnd) out.boundary = bnd;
    return out;
  }

  /* ---------- Arbre complet (affichage, petit nombre de pas) ---------- */
  function tree(model, p, n, lambda) {
    const { S, K, T, r, q, sigma, call, american } = p;
    const dt = T / n, disc = Math.exp(-r * dt), g = Math.exp((r - q) * dt);
    let u, pu, pm, pd;
    if (model === 'bino') { u = Math.exp(sigma * Math.sqrt(dt)); pm = 0; pu = (g - 1 / u) / (u - 1 / u); pd = 1 - pu; }
    else { u = Math.exp(lambda * sigma * Math.sqrt(dt)); pm = 1 - 1 / (lambda * lambda); pu = (g - pm - (1 - pm) / u) / (u - 1 / u); pd = 1 - pm - pu; }
    const step = model === 'bino' ? 2 : 1; // écart d'indice entre deux nœuds voisins
    const levels = [];
    for (let k = 0; k <= n; k++) {
      const nodes = [];
      const jmin = -k, jmax = k;
      for (let j = jmin; j <= jmax; j += step) nodes.push({ j, S: S * Math.pow(u, j), V: 0, ex: false });
      levels.push(nodes);
    }
    const pay = s => call ? Math.max(s - K, 0) : Math.max(K - s, 0);
    levels[n].forEach(nd => { nd.V = pay(nd.S); });
    for (let k = n - 1; k >= 0; k--) {
      const nxt = new Map(levels[k + 1].map(nd => [nd.j, nd]));
      levels[k].forEach(nd => {
        const up = nxt.get(nd.j + (model === 'bino' ? 1 : 1));
        const dn = nxt.get(nd.j - 1);
        const mid = nxt.get(nd.j);
        let cont = model === 'bino'
          ? disc * (pu * up.V + pd * dn.V)
          : disc * (pu * up.V + pm * mid.V + pd * dn.V);
        nd.cont = cont;
        const ex = pay(nd.S);
        if (american && ex > cont + 1e-12) { nd.V = ex; nd.ex = true; } else nd.V = cont;
        nd.delta = (up.V - dn.V) / (up.S - dn.S);
      });
    }
    return { levels, u, pu, pm, pd };
  }

  return { ncdf, npdf, bs, binomial, trinomial, tree };
})();
if (typeof module !== 'undefined') module.exports = Engine;
