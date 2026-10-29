"""Live block IDs of the reader: chapters 0-7 keep the IDs b-NNNNN of the original import; blocks of chapter N >= 8
of the full conversion become cN-NNNNN (appended later as separate chapters)."""
from __future__ import annotations
import re

FIRST_PREFIXED = 8


def chapter_of(blocks):
    """Map block ID -> chapter number, read off the level-1 headings 'N. ...' in document order."""
    chap, out = None, {}
    for b in blocks:
        if b['type'] == 'heading' and b.get('level') == 1 and not b.get('parentId'):
            m = re.match(r'(\d+)\.\s', b['source'])
            if m:
                chap = int(m.group(1))
        out[b['id']] = chap
    return out


def make_remap(blocks):
    ch = chapter_of(blocks)

    def remap(bid):
        c = ch.get(bid)
        if c is not None and c >= FIRST_PREFIXED and re.fullmatch(r'b-\d{5}', bid):
            return f'c{c}-{bid[2:]}'
        return bid
    return remap


def fix_links(s, remap):
    return re.sub(r'\(#(b-\d{5})\)', lambda m: f'(#{remap(m.group(1))})', s)
