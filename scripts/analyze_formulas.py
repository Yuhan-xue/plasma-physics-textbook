import re
import os

base_dir = r"C:\Users\53150\Documents\kimi\workspace\plasma_physics_textbook"

files_to_process = [
    "plasma_physics_textbook_v2.tex",
    "ch1_supplement.tex", "ch2_supplement.tex", "ch3_supplement.tex",
    "ch4_supplement.tex", "ch5_supplement.tex", "ch6_supplement.tex",
    "ch7_supplement.tex", "ch8_supplement.tex", "ch9_supplement.tex",
    "ch10_supplement.tex", "ch11_supplement.tex", "ch12_supplement.tex"
]

# tcolorbox environments that should cause skipping
tcolorbox_envs = ['theorem', 'definition', 'corollary', 'example', 'insight', 
                     'warning', 'history', 'review', 'tip', 'miniquiz', 
                     'chapterreview', 'formula', 'checklist']

def find_matching_end(content, start_pos, env_name):
    """Find the matching \end{env_name} for a \begin{env_name} starting at start_pos."""
def find_matching_end(content, start_pos, env_name):
    pattern = re.compile(rf'\\(begin|end)\\{{{env_name}\\}}}')
    depth = 1
    pos = start_pos
    while depth > 0 and pos < len(content):
        m = pattern.search(content, pos)
        if not m:
            return None
        if m.group(1) == 'begin':
            depth += 1
        else:
            depth -= 1
        if depth == 0:
            return m.end()
        pos = m.end()
    return None

def get_env_ranges(content, env_name):
    """Return list of (start, end) ranges for a given environment."""
    ranges = []
    pattern = re.compile(rf'\\begin\{{{env_name}\}}')
    for m in pattern.finditer(content):
        start = m.start()
        end = find_matching_end(content, m.end(), env_name)
        if end:
            ranges.append((start, end))
    return ranges

def is_inside_any_env(pos, ranges):
    """Check if position is inside any of the given ranges."""
    for start, end in ranges:
        if start < pos < end:
            return True
    return False

def get_line_num(content, pos):
    return content[:pos].count('\n') + 1

def extract_context(content, pos, window=150):
    start = max(0, pos - window)
    end = min(len(content), pos + window)
    return content[start:end].replace('\n', ' ')

def find_all_unnumbered_math(content):
    """Find all unnumbered display math environments."""
    results = []
    
    # Pattern 1: \[ ... \]
    for m in re.finditer(r'\\\[', content):
        start_pos = m.start()
        end_match = re.search(r'\\\]', content[start_pos+2:])
        if end_match:
            end_pos = start_pos + 2 + end_match.end()
            inner = content[start_pos+2:end_pos-2]
            results.append({
                'type': 'displaymath',
                'start': start_pos,
                'end': end_pos,
                'inner': inner,
                'line': get_line_num(content, start_pos)
            })
    
    # Pattern 2: $$ ... $$
    for m in re.finditer(r'(?<!\$)\$\$(?!\$)', content):
        start_pos = m.start()
        end_match = re.search(r'(?<!\$)\$\$(?!\$)', content[start_pos+2:])
        if end_match:
            end_pos = start_pos + 2 + end_match.end()
            inner = content[start_pos+2:end_pos-2]
            results.append({
                'type': 'dollars',
                'start': start_pos,
                'end': end_pos,
                'inner': inner,
                'line': get_line_num(content, start_pos)
            })
    
    # Pattern 3: equation*
    for m in re.finditer(r'\\begin\{equation\*\}', content):
        start_pos = m.start()
        end_match = re.search(r'\\end\{equation\*\}', content[start_pos:])
        if end_match:
            end_pos = start_pos + end_match.end()
            inner = content[start_pos + len(m.group()):end_pos - len(end_match.group())]
            results.append({
                'type': 'equation*',
                'start': start_pos,
                'end': end_pos,
                'inner': inner,
                'line': get_line_num(content, start_pos)
            })
    
    # Pattern 4: align*
    for m in re.finditer(r'\\begin\{align\*\}', content):
        start_pos = m.start()
        end_match = re.search(r'\\end\{align\*\}', content[start_pos:])
        if end_match:
            end_pos = start_pos + end_match.end()
            inner = content[start_pos + len(m.group()):end_pos - len(end_match.group())]
            results.append({
                'type': 'align*',
                'start': start_pos,
                'end': end_pos,
                'inner': inner,
                'line': get_line_num(content, start_pos)
            })
    
    # Pattern 5: gather*
    for m in re.finditer(r'\\begin\{gather\*\}', content):
        start_pos = m.start()
        end_match = re.search(r'\\end\{gather\*\}', content[start_pos:])
        if end_match:
            end_pos = start_pos + end_match.end()
            inner = content[start_pos + len(m.group()):end_pos - len(end_match.group())]
            results.append({
                'type': 'gather*',
                'start': start_pos,
                'end': end_pos,
                'inner': inner,
                'line': get_line_num(content, start_pos)
            })
    
    # Pattern 6: multline*
    for m in re.finditer(r'\\begin\{multline\*\}', content):
        start_pos = m.start()
        end_match = re.search(r'\\end\{multline\*\}', content[start_pos:])
        if end_match:
            end_pos = start_pos + end_match.end()
            inner = content[start_pos + len(m.group()):end_pos - len(end_match.group())]
            results.append({
                'type': 'multline*',
                'start': start_pos,
                'end': end_pos,
                'inner': inner,
                'line': get_line_num(content, start_pos)
            })
    
    results.sort(key=lambda x: x['start'])
    return results

def analyze_file(fpath):
    with open(fpath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    fname = os.path.basename(fpath)
    
    # Get all tcolorbox ranges
    all_tcb_ranges = []
    for env_name in tcolorbox_envs:
        ranges = get_env_ranges(content, env_name)
        all_tcb_ranges.extend(ranges)
    
    # Also get standard theorem-like environments that might be present
    for env_name in ['equation', 'align', 'gather', 'multline', 'split']:
        ranges = get_env_ranges(content, env_name)
        all_tcb_ranges.extend(ranges)
    
    # Find all unnumbered math
    unnumbered = find_all_unnumbered_math(content)
    
    inside_tcb = []
    outside_tcb = []
    
    for item in unnumbered:
        # Check if inside any tcolorbox
        in_tcb = is_inside_any_env(item['start'], all_tcb_ranges)
        if in_tcb:
            inside_tcb.append(item)
        else:
            outside_tcb.append(item)
    
    return {
        'fname': fname,
        'content': content,
        'inside_tcb': inside_tcb,
        'outside_tcb': outside_tcb,
        'tcb_ranges': all_tcb_ranges
    }

all_results = {}
for fname in files_to_process:
    fpath = os.path.join(base_dir, fname)
    all_results[fname] = analyze_file(fpath)
    
    r = all_results[fname]
    print(f"\n{'='*60}")
    print(f"{fname}")
    print(f"  Total unnumbered math: {len(r['inside_tcb']) + len(r['outside_tcb'])}")
    print(f"  Inside tcolorbox: {len(r['inside_tcb'])}")
    print(f"  Outside tcolorbox: {len(r['outside_tcb'])}")
    
    if r['outside_tcb']:
        print("  OUTSIDE formulas:")
        for item in r['outside_tcb']:
            ctx = extract_context(r['content'], item['start'])
            print(f"    Line {item['line']} ({item['type']}): {ctx[:100]}...")
    
    if r['inside_tcb'] and len(r['inside_tcb']) < 30:
        print("  INSIDE tcolorbox formulas:")
        for item in r['inside_tcb']:
            ctx = extract_context(r['content'], item['start'])
            print(f"    Line {item['line']} ({item['type']}): {ctx[:100]}...")
    elif r['inside_tcb']:
        print(f"  INSIDE tcolorbox: {len(r['inside_tcb'])} formulas (too many to list)")

