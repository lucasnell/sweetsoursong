"""Helpers to mark old LaTeX text as deleted with the `changes` package.

Text paragraphs become \\deleted{...}; display equations and other
environments are hidden with \\iffalse ... \\fi behind a visible marker.
"""
import re

def is_env_start(line):
    return re.match(r'\s*\\begin\{(equation|figure|table)', line) is not None

def wrap_deleted(text, marker="equation removed"):
    """Mark a block of old LaTeX as deleted, paragraph by paragraph."""
    lines = text.split('\n')
    out, para, i = [], [], 0
    def flush():
        body = [l for l in para]
        while body and body[-1].strip() == '':
            body.pop()
        if any(l.strip() and not l.strip().startswith('%') for l in body):
            # keep leading \noindent outside the deletion
            prefix = ''
            if body[0].lstrip().startswith('\\noindent'):
                prefix = '\\noindent\n'
                body[0] = body[0].replace('\\noindent', '', 1)
            out.append(prefix + '\\deleted{' + '\n'.join(body) + '\n}')
        elif body:
            out.append('\n'.join(body))
        para.clear()
    while i < len(lines):
        l = lines[i]
        if is_env_start(l):
            flush()
            env = re.match(r'\s*\\begin\{(\w+)', l).group(1)
            j = i
            while not re.match(r'\s*\\end\{' + env + r'\}', lines[j]):
                j += 1
            out.append('\\deleted{[' + marker + ']}\n\\iffalse\n' +
                       '\n'.join(lines[i:j + 1]) + '\n\\fi')
            i = j + 1
            continue
        m = re.match(r'(\s*)\\(sub)*section\*\{(.*)\}\s*$', l)
        if m:
            flush()
            out.append(l.replace('{' + m.group(3) + '}', '{\\deleted{' + m.group(3) + '}}', 1))
            i += 1
            continue
        if l.strip() == '':
            flush()
            out.append('')
        else:
            para.append(l)
        i += 1
    flush()
    return re.sub(r'\n{4,}', '\n\n\n', '\n'.join(out))
