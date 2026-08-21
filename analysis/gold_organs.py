#!/usr/bin/env python3
"""Which organs have a hand-curated gold. One source of truth, read from disk.

This list was duplicated in six modules and had drifted in three of them: calibrate.py
scored four organs when seven existed and reported an AUC of 0.563 instead of 0.707;
gold_limits.py and prune.py still name four. Adding an eighth gold meant editing six
files, and forgetting one produced a silently narrower result rather than an error.

    from gold_organs import curated
    for organ in curated(): ...
"""
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _title(stem):
    return "_".join(w.capitalize() if i == 0 else w
                    for i, w in enumerate(stem.split("_")))


def curated(where=None):
    """Organ names, in the capitalisation the results files use, sorted."""
    d = where or HERE
    return sorted(_title(os.path.basename(p)[:-len("_gold.json")])
                  for p in glob.glob(os.path.join(d, "*_gold.json")))


if __name__ == "__main__":
    got = curated()
    print("%d curated organs: %s" % (len(got), ", ".join(got)))


# The organs curated most recently, after the mapper and the queue were fixed. They are
# never tuned on, so they act as a prospective held-out set. Declared here rather than in
# each figure script: fig2_performance.py carried its own copy naming Lung, Kidney and
# Heart -- the PREVIOUS held-out set -- and went on marking them as held out, and quoting
# their accuracy, for a month after skin, spleen and muscle superseded them.
HELD_OUT = ("Skin", "Spleen", "Muscle")
