// Follow a per-proposal permalink the way a curator arriving from a cell-ontology issue
// would, and assert it lands on the right card.
//
// This exists because the permalink shipped broken once already. Every card carried a
// "link to this proposal" anchor to its own #fragment and none of them worked: the cards
// are rendered from proposals.json after load, so a browser resolves an inbound fragment
// before the element exists and scrolls nowhere. render_smoke.js did not catch it because
// it asserted the anchor was PRESENT, which is the appearance of the feature rather than
// the behaviour. This asserts the behaviour.
//
// Two things it does that render_smoke.js deliberately does not. It exposes `location` as
// a global, since that is what focusFromHash reads and the other harness leaves undefined.
// And its getElementById returns null for ids the page has not emitted -- a document that
// invents an element for every lookup would make "the card was found" vacuously true.
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const DOCS = path.join(__dirname, "..", "docs");
const html = fs.readFileSync(path.join(DOCS, "index.html"), "utf8");
const proposals = JSON.parse(fs.readFileSync(path.join(DOCS, "proposals.json"), "utf8"));

const m = html.match(/<script>([\s\S]*?)<\/script>/);
if (!m) { console.error("FAIL: no <script> block in docs/index.html"); process.exit(1); }

// the fragment to follow: derived the way pid() does, from the first proposal, so the test
// keeps working when the queue changes and never hard-codes a proposal that may be withdrawn
const slug = p => (p.kind + "-" + p.organ + "-" + p.label).toLowerCase()
  .replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "");
const target = slug(proposals.proposals[proposals.proposals.length - 1]);

const log = [];
const known = {};
function mkEl(id) {
  return { id, innerHTML: "", textContent: "", style: {}, dataset: {}, open: false,
           children: [], _details: null,
           setAttribute(k, v) { this[k] = v; }, appendChild(c) { this.children.push(c); },
           querySelector(sel) {
             if (sel === "details" && /<details/.test(this.innerHTML)) {
               this._details = this._details || { open: false }; return this._details; }
             return null; },
           querySelectorAll() { return Object.values(known).filter(x => x._isButton); },
           scrollIntoView() { log.push("scroll:" + this.id); },
           classList: { add: c => log.push("class:" + c), remove() {}, toggle() {} } };
}
for (const id of ["stats", "bar", "list", "cycle"]) known[id] = mkEl(id);

const doc = {
  getElementById: id => known[id] || null,        // null for anything not emitted
  querySelector: () => null, querySelectorAll: () => [],
  createElement: () => { const b = mkEl("_btn"); b._isButton = true; return b; },
  addEventListener() {},
  documentElement: { setAttribute() {}, getAttribute: () => null, dataset: {} },
  body: mkEl("_body")
};
const fetchStub = f => Promise.resolve({
  ok: fs.existsSync(path.join(DOCS, f)),
  json: () => Promise.resolve(JSON.parse(fs.readFileSync(path.join(DOCS, f), "utf8")))
});
const location = { hash: "#" + target };
const ctx = { document: doc, fetch: fetchStub, console, URLSearchParams, setTimeout,
              location,
              window: { location, addEventListener: (e) => log.push("listen:" + e),
                        matchMedia: () => ({ matches: false, addEventListener() {} }) },
              localStorage: { getItem: () => null, setItem() {} } };
ctx.window.document = doc;
vm.createContext(ctx);
try {
  vm.runInContext(m[1], ctx, { filename: "docs/index.html <script>" });
} catch (e) {
  console.error("FAIL: the page's script threw at load: " + e.message);
  process.exit(1);
}

// register every card the render emitted, so getElementById behaves like a real document
function reindex() {
  for (const k of Object.keys(known)) if (k.includes("-")) delete known[k];
  const out = known.list.innerHTML || "";
  let n = 0, mm; const re = /<li class="p" id="([^"]+)"/g;
  while ((mm = re.exec(out))) { known[mm[1]] = mkEl(mm[1]);
                                known[mm[1]].innerHTML = "<details>"; n++; }
  return n;
}

setTimeout(() => {
  const problems = [];
  const cards = reindex();
  if (!cards) problems.push("no cards rendered, so the link cannot be tested");

  // 1. arriving cold on the fragment
  log.length = 0;
  ctx.focusFromHash();
  if (!known[target]) problems.push(`#${target} matches no card`);
  if (!log.includes("scroll:" + target)) problems.push(`#${target} did not scroll into view`);
  const det = known[target] && known[target]._details;
  if (!det || !det.open) problems.push(`#${target} did not open its evidence`);
  if (!log.some(l => l.startsWith("class:"))) problems.push(`#${target} was not highlighted`);

  // 2. a filter is hiding the target: the link must still resolve
  if (typeof ctx.setFilter === "function") {
    ctx.setFilter("marker-condition");
    const hidden = reindex();
    log.length = 0;
    ctx.focusFromHash();
    reindex();
    if (!known[target])
      problems.push(`#${target} unreachable while a filter is active (${hidden} cards shown)`);
    ctx.setFilter("all");
    reindex();
  } else {
    problems.push("setFilter is not exposed, so the filtered case cannot be tested");
  }

  // 3. a fragment matching nothing is a no-op, never a throw
  ctx.location.hash = "#no-such-proposal";
  try { ctx.focusFromHash(); }
  catch (e) { problems.push("an unknown fragment threw: " + e.message); }

  // 4. every anchor on the page points at an id the page emits
  const emitted = new Set(Object.keys(known).filter(k => k.includes("-")));
  const hrefs = [...(known.list.innerHTML.matchAll(/href="#([^"]+)">link to this proposal/g))]
    .map(x => x[1]);
  const dangling = hrefs.filter(h => !emitted.has(h));
  if (hrefs.length !== cards) problems.push(`${hrefs.length} permalinks for ${cards} cards`);
  if (dangling.length) problems.push(`permalink points at no card: ${dangling[0]}`);

  if (problems.length) { problems.forEach(p => console.error("FAIL: " + p)); process.exit(1); }
  console.log(`deep links OK: ${cards} cards, #${target} resolves cold and under a filter, ` +
              `${hrefs.length} anchors all point at emitted ids`);
}, 40);
