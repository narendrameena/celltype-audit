// Render the review queue with the real proposals.json and assert it produced cards.
//
// This exists because index.html shipped broken: a helper was declared inside another
// function, so calling it from the renderer threw a ReferenceError, and the page's catch
// handler reported that as "Could not load proposals.json" -- blaming the data for a bug
// in the code. No JavaScript runtime is available where the page is authored, so nothing
// caught it before deploy. GitHub's runners have node; this runs there, before the deploy
// step, and fails the build rather than publishing a page that renders nothing.
//
// It is not a browser. It provides the smallest DOM and fetch the page actually touches,
// which is enough to execute the render path end to end against real data.
const fs = require("fs");
const path = require("path");

const DOCS = path.join(__dirname, "..", "docs");
const html = fs.readFileSync(path.join(DOCS, "index.html"), "utf8");

const m = html.match(/<script>([\s\S]*?)<\/script>/);
if (!m) { console.error("FAIL: no <script> block in docs/index.html"); process.exit(1); }

const els = {};
function el(id) {
  if (!els[id]) els[id] = { innerHTML: "", textContent: "", style: {}, dataset: {},
                            classList: { add(){}, remove(){}, toggle(){} },
                            addEventListener(){}, querySelectorAll: () => [],
                            setAttribute(){}, appendChild(){} };
  return els[id];
}
const doc = {
  getElementById: el,
  querySelector: () => el("_q"),
  querySelectorAll: () => [],
  createElement: () => el("_c"),
  addEventListener(){},
  documentElement: { setAttribute(){}, getAttribute: () => null, dataset: {} },
  body: el("_body")
};

const json = f => JSON.parse(fs.readFileSync(path.join(DOCS, f), "utf8"));
const fetchStub = f => Promise.resolve({
  ok: fs.existsSync(path.join(DOCS, f)),
  json: () => Promise.resolve(json(f)),
  text: () => Promise.resolve(fs.readFileSync(path.join(DOCS, f), "utf8"))
});

const ctx = { document: doc, fetch: fetchStub, window: { matchMedia: () => ({ matches: false,
              addEventListener(){} }), location: { hash: "", search: "" }, addEventListener(){} },
              console, URLSearchParams, setTimeout, localStorage: { getItem: () => null,
              setItem(){} } };
ctx.window.document = doc;

const vm = require("vm");
vm.createContext(ctx);
try {
  vm.runInContext(m[1], ctx, { filename: "docs/index.html <script>" });
} catch (e) {
  console.error("FAIL: the page's script threw at load: " + e.message);
  process.exit(1);
}

setTimeout(() => {
  const out = el("list").innerHTML;
  const proposals = json("proposals.json").proposals;
  // matches the card open tag whether or not it carries attributes; the previous form
  // was an exact string and reported zero cards the moment an id was added to it
  const cards = (out.match(/<li class="p"[ >]/g) || []).length;
  const problems = [];

  if (/Could not load/.test(out)) problems.push("page reports it could not load proposals.json");
  if (/failed while rendering/.test(out)) problems.push("page reports a render failure: " + out.slice(0, 300));
  if (cards === 0) problems.push("no proposal cards were rendered");
  if (cards && cards !== proposals.length)
    problems.push(`rendered ${cards} cards for ${proposals.length} proposals`);
  // every card should offer the background link the proposals page promises
  const bg = (out.match(/class="bg"/g) || []).length;
  if (bg !== cards) problems.push(`${bg} background links for ${cards} cards`);
  // every card must carry its own derived reasoning and a stable anchor to link to
  const why = (out.match(/<details class="why">/g) || []).length;
  if (why !== cards) problems.push(`${why} reasoning blocks for ${cards} cards`);
  const ids = (out.match(/<li class="p" id="[^"]+"/g) || []).length;
  if (ids !== cards) problems.push(`${ids} cards have a stable id, expected ${cards}`);
  const perma = (out.match(/link to this proposal/g) || []).length;
  if (perma !== cards) problems.push(`${perma} permalinks for ${cards} cards`);
  if (/undefined|\[object Object\]|NaN/.test(out))
    problems.push("rendered output contains undefined/NaN - a field name is wrong");
  // proposals.json carries a lexical overlap score and no subsumption relation, so the
  // page must not describe a nearest match as broader or narrower than the label. It did,
  // and named a narrower term as broader.
  // Target the CLAIM, not the vocabulary. "Neither subsumes the other" is a verified
  // negative computed from the graph and must stay; "broader than" asserted about a
  // lexical near-match is the thing that was wrong.
  const bad = out.match(/\b(broader|narrower) than\b|\bwhich subsumes\b/gi);
  if (bad) problems.push("prose claims a subclass direction proposals.json does not carry: "
                         + bad.join(", "));

  if (problems.length) {
    problems.forEach(p => console.error("FAIL: " + p));
    process.exit(1);
  }
  console.log(`ok: ${cards} cards, ${bg} background links, ${why} reasoning blocks, ` +
              `${ids} anchored ids`);
}, 50);
