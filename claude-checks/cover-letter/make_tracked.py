"""Write a copy of the cover letter with word-level tracked changes.

Usage (2026-10-08): python3 make_tracked.py \
    ~/Box/_sweetandsour/_drafts/cover-Nell-nectar-microbes.docx \
    ~/Box/_sweetandsour/_drafts/cover-Nell-nectar-microbes-PNAS-tracked.docx

Each edited paragraph is rebuilt from its original runs: unchanged text keeps
its original formatting, deletions become <w:del>, insertions <w:ins>.
In new text, *...* marks italics.
"""
import difflib
import html
import re
import shutil
import sys
import zipfile

SRC = sys.argv[1]
DST = sys.argv[2]
AUTHOR = "Claude"
DATE = "2026-10-08T00:00:00Z"

# paragraph index -> new text (None = unchanged)
NEW = {
    1: "8 October 2026",
    3: "Editorial Board",
    4: "*Proceedings of the National Academy of Sciences*",
    6: "Dear Editors,",
    8: ("We would like the manuscript entitled “Regional species coexistence "
        "despite local priority effects: the overlooked role of "
        "dispersal–community feedback” to be considered for publication as a "
        "Research Report in *PNAS*. The manuscript has not been published and "
        "is not under consideration for publication elsewhere."),
    10: ("Many organisms rely on animals to move between habitats, and those "
         "animals often choose where to go based on which species are already "
         "there. This feedback between dispersal and community composition has "
         "received little attention in theories of species coexistence. In many "
         "spatial coexistence mechanisms, local interactions generate spatial "
         "heterogeneity that allows coexistence at larger scales. Local priority "
         "effects, where the order in which species arrive determines the "
         "outcome of community assembly, also create spatial heterogeneity. "
         "However, existing theory predicts that local priority effects "
         "destabilize regional coexistence, so it is unclear whether regional "
         "coexistence can result from them."),
    12: ("Here, we use a mathematical model to show that dispersal–community "
         "feedback can turn local priority effects into regional coexistence. "
         "We model yeast and bacteria that compete in floral nectar and are "
         "carried between flowers by pollinators, which leave bacteria-dominated "
         "flowers more readily because bacteria sour the nectar. On each plant, "
         "this feedback creates alternative stable states: yeast-dominated "
         "plants keep their pollinators and bacteria-dominated plants lose them, "
         "so whichever microbe is common on a plant reinforces its own "
         "dominance. In a closed metacommunity of plants sharing a regional pool "
         "of pollinators, the same feedback lets each microbe invade when rare, "
         "and the two coexist across a nearly four-fold range of regional "
         "pollinator abundance. Without the feedback, coexistence is limited to "
         "a narrow range. Positive feedback within plants thus becomes negative "
         "feedback across the landscape, giving each microbe a set of plants "
         "where it is competitively dominant."),
    14: ("Our work shows how feedback between community states and dispersal "
         "can change the outcome of competition, and identifies an overlooked "
         "mechanism of species coexistence. Because many organisms depend on "
         "animals for dispersal, including plants with animal-dispersed seeds, "
         "vector-borne pathogens, and organisms living in ephemeral habitats, we "
         "expect the result to interest a broad readership in community "
         "ecology, disease ecology, and plant–animal interactions. We developed "
         "a new stochastic model for this system, and all code, including an "
         "open-source R package that reproduces every result, is archived on "
         "Zenodo."),
}

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
rev_id = [9000]


def next_id():
    rev_id[0] += 1
    return rev_id[0]


def tokens(text):
    return re.findall(r"\s+|\w+(?:[’']\w+)*|[^\w\s]", text)


def run_xml(text, rpr, deleted=False):
    tag = "w:delText" if deleted else "w:t"
    return ('<w:r>' + (f'<w:rPr>{rpr}</w:rPr>' if rpr else '') +
            f'<{tag} xml:space="preserve">{html.escape(text, quote=False)}</{tag}></w:r>')


def italic(rpr):
    if '<w:i/>' in rpr:
        return rpr
    return rpr + '<w:i/><w:iCs/>'


def plain(rpr):
    return rpr.replace('<w:i/>', '').replace('<w:iCs/>', '')


def runs_of(chars):
    """Group (char, rPr) into runs of equal formatting."""
    out = []
    for ch, rpr in chars:
        if out and out[-1][1] == rpr:
            out[-1][0] += ch
        else:
            out.append([ch, rpr])
    return out


def parse_new(text, base_rpr):
    """New text with *italic* markup -> list of (char, rPr)."""
    chars = []
    for i, seg in enumerate(re.split(r'\*', text)):
        r = italic(base_rpr) if i % 2 else plain(base_rpr)
        chars.extend((c, r) for c in seg)
    return chars


def rebuild(p_xml, new_text):
    ppr = re.search(r'<w:pPr>.*?</w:pPr>', p_xml, flags=re.S)
    ppr = ppr.group(0) if ppr else ''
    start = re.match(r'<w:p[^>]*>', p_xml).group(0)
    old_chars = []
    for r in re.findall(r'<w:r[ >].*?</w:r>', p_xml, flags=re.S):
        rpr = re.search(r'<w:rPr>(.*?)</w:rPr>', r, flags=re.S)
        rpr = rpr.group(1) if rpr else ''
        for t in re.findall(r'<w:t[^>]*>(.*?)</w:t>', r, flags=re.S):
            old_chars.extend((c, rpr) for c in html.unescape(t))
    base = plain(old_chars[0][1]) if old_chars else ''
    new_chars = parse_new(new_text, base)
    old_text = ''.join(c for c, _ in old_chars)
    new_plain = ''.join(c for c, _ in new_chars)
    a, b = tokens(old_text), tokens(new_plain)
    # map token index -> char offset
    def offsets(toks):
        o, acc = [], 0
        for t in toks:
            o.append(acc)
            acc += len(t)
        o.append(acc)
        return o
    oa, ob = offsets(a), offsets(b)
    body = []
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == 'equal':
            # keep original formatting unless italics changed
            for k in range(i2 - i1):
                ca = old_chars[oa[i1 + k]:oa[i1 + k + 1]]
                cb = new_chars[ob[j1 + k]:ob[j1 + k + 1]]
                if [r for _, r in ca] == [r for _, r in cb] or '*' not in new_text:
                    body.extend(run_xml(t, r) for t, r in runs_of(ca))
                else:
                    body.append(f'<w:del w:id="{next_id()}" w:author="{AUTHOR}" w:date="{DATE}">' +
                                ''.join(run_xml(t, r, True) for t, r in runs_of(ca)) + '</w:del>')
                    body.append(f'<w:ins w:id="{next_id()}" w:author="{AUTHOR}" w:date="{DATE}">' +
                                ''.join(run_xml(t, r) for t, r in runs_of(cb)) + '</w:ins>')
            continue
        if i2 > i1:
            ca = old_chars[oa[i1]:oa[i2]]
            body.append(f'<w:del w:id="{next_id()}" w:author="{AUTHOR}" w:date="{DATE}">' +
                        ''.join(run_xml(t, r, True) for t, r in runs_of(ca)) + '</w:del>')
        if j2 > j1:
            cb = new_chars[ob[j1]:ob[j2]]
            body.append(f'<w:ins w:id="{next_id()}" w:author="{AUTHOR}" w:date="{DATE}">' +
                        ''.join(run_xml(t, r) for t, r in runs_of(cb)) + '</w:ins>')
    return start + ppr + ''.join(body) + '</w:p>'


with zipfile.ZipFile(SRC) as z:
    files = {n: z.read(n) for n in z.namelist()}
doc = files['word/document.xml'].decode('utf-8')
paras = list(re.finditer(r'<w:p[ >].*?</w:p>', doc, flags=re.S))
out, last = [], 0
for i, m in enumerate(paras):
    out.append(doc[last:m.start()])
    out.append(rebuild(m.group(0), NEW[i]) if i in NEW else m.group(0))
    last = m.end()
out.append(doc[last:])
files['word/document.xml'] = ''.join(out).encode('utf-8')
with zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED) as z:
    for n, data in files.items():
        z.writestr(n, data)
print("wrote", DST)
