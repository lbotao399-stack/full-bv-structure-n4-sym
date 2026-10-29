"""Resolve narration from source prose; formula-only blocks inherit nearby prose."""
import re


def source_text(block):
    if block.get('type') == 'image':
        return str(block.get('caption') or block.get('alt') or '')
    return str(block.get('source') or '')


def prose_language(block):
    if block.get('type') == 'equation':
        return None
    text = re.sub(r'\$\$[\s\S]*?\$\$|\$`[\s\S]*?`\$|\$[^$]*?\$|\\\[[\s\S]*?\\\]|\\\([\s\S]*?\\\)', ' ', source_text(block))
    text = re.sub(r'<[^>]*>', ' ', text)
    if re.search(r'[\u3400-\u9fff]', text):
        return 'zh-CN'
    if re.search(r'[A-Za-z]{2,}', text):
        return 'en-US'
    return None


def source_language(block, document):
    direct = prose_language(block)
    if direct:
        return direct
    blocks = document.get('blocks', [])
    at = next((i for i, b in enumerate(blocks) if b.get('id') == block.get('id')), None)
    if at is not None:
        section = block.get('section')
        # A contained formula follows its own callout's language.
        if block.get('parentId'):
            parent = next((b for b in blocks if b.get('id') == block['parentId']), None)
            if parent and prose_language(parent):
                return prose_language(parent)
        for distance in range(1, len(blocks)):
            for i in (at-distance, at+distance):
                if 0 <= i < len(blocks):
                    candidate = blocks[i]
                    if section and candidate.get('section') not in (None, section):
                        continue
                    language = prose_language(candidate)
                    if language:
                        return language
    return 'zh-CN' if re.search(r'[\u3400-\u9fff]', document.get('title', '')) else 'en-US'
