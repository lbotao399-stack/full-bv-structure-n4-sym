"""Explicit clipboard handoff. Exact content matches only; never a nearest passage."""
from html.parser import HTMLParser
from html import unescape
import hashlib, re, unicodedata

class ClipboardHTML(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.root={'tag':'root','attrs':{},'children':[]};self.stack=[self.root];self.math_count=0
    def handle_starttag(self,tag,attrs):
        node={'tag':tag,'attrs':dict(attrs),'children':[]};self.stack[-1]['children'].append(node)
        if tag not in {'br','hr','img','meta','link','input','wbr'}:self.stack.append(node)
    def handle_startendtag(self,tag,attrs):self.handle_starttag(tag,attrs);self.handle_endtag(tag)
    def handle_endtag(self,tag):
        for i in range(len(self.stack)-1,0,-1):
            if self.stack[i]['tag']==tag:self.stack=self.stack[:i];break
    def handle_data(self,data):self.stack[-1]['children'].append(data)
    def annotation(self,node):
        if isinstance(node,str):return None
        if node['tag']=='annotation' and node['attrs'].get('encoding') in {'application/x-tex','text/latex'}:return ''.join(self.literal(c) for c in node['children'])
        for c in node['children']:
            result=self.annotation(c)
            if result is not None:return result
    def literal(self,node):return node if isinstance(node,str) else ''.join(self.literal(c) for c in node['children'])
    def render(self,node=None):
        node=self.root if node is None else node
        if isinstance(node,str):return node
        tag=node['tag'];attrs=node['attrs']
        if tag in {'script','style','noscript','button','meta','link'}:return ''
        if tag=='math' or 'katex' in attrs.get('class','').split():
            tex=self.annotation(node)
            if tex is not None:self.math_count+=1;return '$'+tex+'$'
        if attrs.get('data-equation'):
            self.math_count+=1;return '$'+attrs['data-equation']+'$'
        result=''.join(self.render(c) for c in node['children'])
        if tag in {'p','div','blockquote','h1','h2','h3','h4','li','tr'}:return '\n'+result+'\n'
        if tag=='br':return '\n'
        if tag in {'td','th'}:return result+'; '
        return result

MATH=re.compile(r'\$`([\s\S]*?)`\$|\$\$([\s\S]*?)\$\$|\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)|\$([^$]+)\$')

def plain_markdown(text):
    protected=[]
    def hold(m):
        protected.append(next(x for x in m.groups() if x is not None));return f'MATHCLIPTOKEN{len(protected)-1}END'
    text=MATH.sub(hold,text)
    parser=ClipboardHTML();parser.feed(text);text=parser.render()
    text=re.sub(r'MATHCLIPTOKEN(\d+)END',lambda m:protected[int(m[1])],text)
    text=re.sub(r'\[([^\]]+)\]\([^)]+\)',r'\1',text)
    return text

COMMANDS={'cup':'∪','cap':'∩','in':'∈','notin':'∉','setminus':'∖','times':'×','cdot':'·','pm':'±','mp':'∓','le':'≤','leq':'≤','ge':'≥','geq':'≥','ne':'≠','neq':'≠','mapsto':'↦','longmapsto':'↦','to':'→','rightarrow':'→','longrightarrow':'→','cong':'≅','equiv':'≡','sigma':'σ','mu':'μ','nu':'ν','lambda':'λ','chi':'χ','pi':'π','theta':'θ','alpha':'α','beta':'β','gamma':'γ','delta':'δ','epsilon':'ε','varepsilon':'ε','phi':'ϕ','varphi':'φ','ell':'ℓ','infty':'∞','ldots':'…','cdots':'⋯','quad':'','qquad':'','left':'','right':'','middle':'','mathbb':'','mathrm':'','mathbf':'','mathit':'','text':'','operatorname':''}
def key(text):
    text=plain_markdown(text)
    # Preserve variance and exponent grouping while normalizing typographic scripts.
    lower=dict(zip('₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓ', '0123456789+-=()aehijklmnoprstuvx'))
    upper=dict(zip('⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿⁱ', '0123456789+-=()ni'))
    for alphabet,prefix in [(lower,'_'),(upper,'^')]:
        text=re.sub('['+''.join(alphabet)+']+',lambda m:prefix+'{'+''.join(alphabet[c] for c in m[0])+'}',text)
    for _ in range(4):text=re.sub(r'\\(?:mathrm|mathbf|mathit|mathbb|text|operatorname)\{([^{}]*)\}',r'\1',text)
    text=re.sub(r'([_^])\s*\{(\\[A-Za-z]+|[A-Za-z0-9])\}',r'\1\2',text)
    text=re.sub(r'\\([A-Za-z]+)',lambda m:COMMANDS.get(m[1],m[0]),text)
    text=text.replace('−','-').replace('–','-').replace('’',"'")
    text=unicodedata.normalize('NFKC',text)
    # Ignore typographic delimiters, retaining operation symbols and their order.
    text=re.sub(r'\\([ ,;!])',r'\1',text)
    return re.sub(r'[\s$`*#]','',unescape(text)).strip()

def selection(payload,document):
    text=payload.get('text','');html=payload.get('html','')
    if not isinstance(text,str) or not isinstance(html,str):raise ValueError('Copy text or a formula from the Notion page first.')
    if len(text)>100000 or len(html)>2000000:raise ValueError('Select a smaller passage to read.')
    parser=ClipboardHTML();parser.feed(html);rich=parser.render().strip()
    source=rich if parser.math_count else text.strip() or rich
    source=re.sub(r'\n[ \t]*\n(?:[ \t]*\n)+','\n\n',source).strip()
    if not source:raise ValueError('Copy a passage in Codex Notion, then click Read copied selection.')
    target=key(source)
    blocks=document['blocks']
    keys=[key(b.get('caption') or b.get('alt') or '') if b['type']=='image' else key(b['source']) for b in blocks]
    # Match complete blocks or a complete consecutive range. Partial selections are read literally.
    for start in range(len(keys)):
        accum=''
        for end in range(start,len(keys)):
            accum+=keys[end]
            if accum==target and accum:
                ids=[b['id'] for b in blocks[start:end+1]]
                return {'kind':'source','startBlockId':ids[0],'blockIds':ids,'selectionText':source,'match':'exact normalized complete block range'}
            if len(accum)>=len(target) or not target.startswith(accum):break
    kind='paragraph'
    match=MATH.fullmatch(source)
    if match:
        source=next(x for x in match.groups() if x is not None);kind='equation'
    elif re.search(r'\\(?:frac|sum|begin|mathrm|sigma|chi|lambda|mathbb|left|sqrt)\b',source) and '$' not in source and not re.search(r'\b(?:the|and|for|is|are|with|let|of|in)\b',re.sub(r'\\[A-Za-z]+','',source),re.I):kind='equation'
    return {'kind':'selection','block':{'id':'selection-'+hashlib.sha256(source.encode()).hexdigest()[:16],'type':kind,'source':source},'selectionText':source,'match':'literal copied selection; no source match'}
