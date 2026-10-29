"""Validated, content-addressed PNG storage for a local manuscript reader."""
from pathlib import Path
import base64
import binascii
import hashlib
import os
import struct
import tempfile
import zlib

MAX_IMAGE_BYTES = 16 * 1024 * 1024
MAX_JSON_BYTES = 24 * 1024 * 1024
MAX_PIXELS = 32_000_000
PNG_SIGNATURE = b'\x89PNG\r\n\x1a\n'


def capabilities():
    return {'supported': True, 'acceptedMimeTypes': ['image/png'],
            'maxImageBytes': MAX_IMAGE_BYTES, 'maxPixels': MAX_PIXELS}


def validate_png(data):
    """Accept noninterlaced 8-bit RGB/RGBA PNGs emitted by browser canvas.

    Check all chunk CRCs and the entire bounded pixel stream. No imaging library
    or image executable is loaded. Auxiliary chunks are not interpreted.
    """
    if not 1 <= len(data) <= MAX_IMAGE_BYTES or not data.startswith(PNG_SIGNATURE):
        raise ValueError('Paste a PNG image of at most 16 MiB.')
    pos = 8
    width = height = channels = None
    compressed = []
    saw_idat = ended_idat = saw_end = False
    while pos < len(data):
        if pos + 12 > len(data):
            raise ValueError('Incomplete PNG chunk.')
        length = struct.unpack_from('>I', data, pos)[0]
        kind = data[pos + 4:pos + 8]
        end = pos + 12 + length
        if end > len(data) or not all(65 <= x <= 90 or 97 <= x <= 122 for x in kind):
            raise ValueError('Invalid PNG chunk.')
        payload = data[pos + 8:pos + 8 + length]
        crc = struct.unpack_from('>I', data, pos + 8 + length)[0]
        if binascii.crc32(kind + payload) & 0xffffffff != crc:
            raise ValueError('The PNG checksum is invalid.')
        if width is None and kind != b'IHDR':
            raise ValueError('The PNG header is missing.')
        if kind == b'IHDR':
            if width is not None or length != 13:
                raise ValueError('Invalid PNG header.')
            width, height, depth, color, compression, filtering, interlace = struct.unpack('>IIBBBBB', payload)
            if not width or not height or width * height > MAX_PIXELS:
                raise ValueError('Images must contain between 1 and 32 million pixels.')
            if depth != 8 or color not in (2, 6) or compression or filtering or interlace:
                raise ValueError('Convert the image to an 8-bit RGB or RGBA PNG before pasting.')
            channels = 3 if color == 2 else 4
        elif kind == b'IDAT':
            if ended_idat:
                raise ValueError('Invalid PNG image-data order.')
            saw_idat = True
            compressed.append(payload)
        elif kind == b'IEND':
            if length or not saw_idat or end != len(data):
                raise ValueError('Invalid PNG ending.')
            saw_end = True
        elif kind == b'PLTE':
            if saw_idat or length < 3 or length > 768 or length % 3:
                raise ValueError('Invalid PNG palette.')
        elif not kind[0] & 32:
            raise ValueError('Unsupported PNG critical chunk.')
        if saw_idat and kind != b'IDAT':
            ended_idat = True
        pos = end
    if not saw_end:
        raise ValueError('Incomplete PNG image.')
    # Browser canvas emits standard row filters. Stream the decode so a
    # compressed bomb cannot allocate more than one small output chunk at once.
    stride = 1 + width * channels
    expected = stride * height
    offset = 0
    decoder = zlib.decompressobj()
    try:
        for item in compressed:
            pending = item
            while pending:
                pixels = decoder.decompress(pending, min(65536, expected - offset + 1))
                pending = decoder.unconsumed_tail
                if offset + len(pixels) > expected:
                    raise ValueError('PNG decoded data exceeds its declared dimensions.')
                first = (-offset) % stride
                if any(pixels[i] > 4 for i in range(first, len(pixels), stride)):
                    raise ValueError('Invalid PNG row filter.')
                offset += len(pixels)
                if decoder.unused_data:
                    raise ValueError('Extra PNG compressed data.')
                if not pending:
                    break
        if not decoder.eof or offset != expected:
            raise ValueError('PNG pixel data does not match its dimensions.')
    except zlib.error as exc:
        raise ValueError('Invalid PNG compression.') from exc
    return width, height


def store_image(reader_root, payload):
    if not isinstance(payload, dict) or not isinstance(payload.get('dataUrl'), str):
        raise ValueError('Expected an image data URL.')
    value = payload['dataUrl']
    prefix = 'data:image/png;base64,'
    if not value.startswith(prefix):
        raise ValueError('Only normalized PNG images may be uploaded.')
    encoded = value[len(prefix):]
    if len(encoded) > 4 * ((MAX_IMAGE_BYTES + 2) // 3):
        raise ValueError('Paste a PNG image of at most 16 MiB.')
    try:
        data = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as exc:
        raise ValueError('Invalid base64 image data.') from exc
    width, height = validate_png(data)
    digest = hashlib.sha256(data).hexdigest()
    root = Path(reader_root).resolve()
    directory = root / 'figures' / 'pasted'
    directory.mkdir(parents=True, exist_ok=True)
    if not directory.resolve().is_relative_to(root / 'figures'):
        raise ValueError('The image directory must stay inside this reader.')
    target = directory / (digest + '.png')
    if target.is_symlink():
        raise ValueError('An image target must be a regular local file.')
    if not target.exists():
        fd, temp_name = tempfile.mkstemp(prefix='.image-', suffix='.tmp', dir=directory)
        try:
            with os.fdopen(fd, 'wb') as stream:
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temp_name, target)
        finally:
            if os.path.exists(temp_name):
                os.unlink(temp_name)
    elif not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
        raise ValueError('An existing image file has unexpected content.')
    return {'source': 'figures/pasted/' + target.name, 'width': width,
            'height': height, 'bytes': len(data), 'sha256': digest}
