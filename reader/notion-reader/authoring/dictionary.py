"""Read-only local English–Chinese dictionary; selection never leaves the reader."""
import re
import sqlite3
from pathlib import Path

DEFAULT_PATH = Path.home() / 'Library/Application Support/AcademicAudio/Dictionaries/ecdict.sqlite3'
GLOSSARY = {
    'two-point function': '两点函数',
    'three-point function': '三点函数',
    'n-point function': 'n 点函数',
    'finite-size scaling': '有限尺寸标度',
    'boundary condition': '边界条件',
    'correlation function': '关联函数',
    'conformal field theory': '共形场论',
}

def lookup(term, configured_path=None):
    if not isinstance(term, str) or len(term) > 160:
        raise ValueError('Select an English word or short phrase.')
    query = ' '.join(term.strip().split()).replace('’', "'").replace('‐', '-').replace('‑', '-')
    query = query.strip('.,;:!?“”"()[]{}')
    if not query or len(query) > 80 or len(query.split()) > 8 or not re.fullmatch(r"[A-Za-z][A-Za-z '\-]*", query):
        raise ValueError('Select an English word or a phrase of up to eight words, without formulas.')
    key = query.casefold()
    for phrase, translation in GLOSSARY.items():
        if key in (phrase, phrase.replace('-', ' ')):
            return {'query': query, 'entries': [{'word': phrase, 'translation': translation,
                    'source': 'Reader physics glossary', 'match': 'phrase'}]}
    path = Path(configured_path).expanduser() if configured_path else DEFAULT_PATH
    if not path.is_file():
        return {'query': query, 'entries': [], 'unavailable': True,
                'message': 'The local dictionary is not installed.'}
    with sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True, timeout=3) as db:
        db.row_factory = sqlite3.Row
        rows = db.execute('SELECT word,phonetic,translation,definition FROM entries WHERE word=?', (key,)).fetchall()
        match = 'exact'
        if not rows:
            rows = db.execute('SELECT e.word,e.phonetic,e.translation,e.definition FROM forms f JOIN entries e ON e.word=f.word WHERE f.form=? LIMIT 4', (key,)).fetchall()
            match = 'inflection'
    return {'query': query, 'entries': [dict(row, source='ECDICT · MIT', match=match) for row in rows],
            'message': '' if rows else 'No local entry for this selection.'}
