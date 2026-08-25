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
  // a background link only where it is specific to the kind; a card's own graph does the
  // rest. Asserting one per card is what let 17 identical #anchors links look intentional.
  const bg = (out.match(/class="bg"/g) || []).length;
  const kindsWithBg = proposals.filter(p => ["marker-condition", "missing-axiom"].includes(p.kind)).length;
  if (bg !== kindsWithBg) problems.push(`${bg} background links, expected ${kindsWithBg}`);
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
  // V10 and V12 are the checks that make a proposal falsifiable rather than merely
  // stated. They were hardcoded "not-run" on every proposal for weeks; if they regress to
  // that, the page is claiming a verification stack it did not run.
  const ranV = proposals.filter(p => ["pass", "fail"].includes((p.checks || {}).V10) ||
                                 ["pass", "fail"].includes((p.checks || {}).V12)).length;
  if (ranV === 0) problems.push("no proposal has V10 or V12 actually run");
  // every proposal must carry a verdict for every check in the stack -- an absent key is
  // a check that silently stopped applying, which reads on the page as if it never existed
  const STACK = ["V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V12"];
  const gaps = proposals.filter(p => STACK.some(v => !(p.checks || {})[v]));
  if (gaps.length) problems.push(`${gaps.length} proposals missing a verdict for some check`);
  // a proposal that failed exclusivity must say so where a curator will see it
  // every proposal must carry its own ontology graph, read from the release rather than
  // a shared explainer link: the generic one told a curator what an anchor set is and
  // nothing about the proposal in front of them
  // every verdict must carry a reason, and every reason must reach the page. A bare n/a
  // is indistinguishable from a skipped check, which is the ambiguity this removes.
  const STACK2 = ["V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V12"];
  const unexplained = proposals.reduce((n, p) =>
    n + STACK2.filter(v => p.checks[v] && !(p.check_detail || {})[v]).length, 0);
  if (unexplained) problems.push(`${unexplained} verdicts carry no reason`);
  const tables = (out.match(/<table class="vtab">/g) || []).length;
  if (tables !== proposals.length)
    problems.push(`${tables} check tables for ${proposals.length} proposals`);

  const graphed = proposals.filter(p => p.graph).length;
  if (graphed !== proposals.length)
    problems.push(`${graphed}/${proposals.length} proposals carry an ontology graph`);
  const svgs = (out.match(/<figure class="cg">/g) || []).length;
  if (svgs !== proposals.length) problems.push(`${svgs} graphs rendered for ${proposals.length} proposals`);

  const weak = proposals.filter(p => p.readiness === "weakened").length;
  const marks = (out.match(/Not ready to submit/g) || []).length;
  if (marks !== weak) problems.push(`${marks} "not ready" notices for ${weak} weakened proposals`);
  const detailed = proposals.filter(p => p.check_detail && (p.check_detail.V10 || p.check_detail.V12)).length;
  if (detailed !== proposals.length)
    problems.push(`${detailed}/${proposals.length} proposals carry falsification detail`);

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
