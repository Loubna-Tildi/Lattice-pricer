// node src/test.js : vérifie le moteur contre des valeurs de référence
const E = require('./engine.js');
const p = { S:100, K:100, T:1, r:0.05, q:0, sigma:0.2, call:true, american:false };
const ok = (name, v, ref, tol) => console.log((Math.abs(v - ref) < tol ? 'OK  ' : 'FAIL') + ` ${name}: ${v.toFixed(6)} (réf ${ref})`);
ok('BS call', E.bs(p).price, 10.4506, 1e-4);
ok('BS put', E.bs({ ...p, call:false }).price, 5.5735, 1e-4);
ok('Binomial N=1000 call', E.binomial(p, 1000).price, 10.4506, 3e-3);
ok('Trinomial N=1000 call', E.trinomial(p, 1000).price, 10.4506, 1e-3);
ok('Put américain trinomial N=2500', E.trinomial({ ...p, call:false, american:true }, 2500).price, 6.0903, 5e-4);
ok('λ = 1 redonne le binomial', E.trinomial(p, 100, 1).price, E.binomial(p, 100).price, 1e-9);
