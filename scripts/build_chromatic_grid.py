from pathlib import Path
import hashlib
import json
import subprocess

import sys,os
ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'build/grid'
HERE.mkdir(parents=True,exist_ok=True)
BOOST = Path(sys.argv[1] if len(sys.argv)>1 else os.environ['BOOST_ROOT'])
header_source = ROOT/'cpp/include/chord/gl_weight.hpp'
cpp_source = ROOT/'cpp/src/gl_weight.cpp'
header = header_source.read_text()
assert header.count('kResidueLanes = 3;') == 1
header = header.replace('kResidueLanes = 3;', 'kResidueLanes = 110;')
local_include = HERE/'generated/include/chord'
local_include.mkdir(parents=True, exist_ok=True)
(local_include/'gl_weight.hpp').write_text(header)
cpp = cpp_source.read_text()
needle = ": primes_{1'000'000'007ULL, 1'000'000'009ULL, 998'244'353ULL} {"
assert cpp.count(needle) == 1
cpp = cpp.replace(needle, ': primes_{} {\n    primes_.fill(1000000007ULL);')
needle_mul = 'const auto product = static_cast<unsigned __int128>(lhs.value[lane]) * rhs.value[lane];'
assert cpp.count(needle_mul) == 1
cpp = cpp.replace(needle_mul, 'const auto product = lhs.value[lane] * rhs.value[lane];')
needle_mod = 'result.value[lane] = static_cast<std::uint64_t>(product % primes_[lane]);'
assert cpp.count(needle_mod) == 1
cpp = cpp.replace(needle_mod, 'result.value[lane] = product % 1000000007ULL;')
(HERE/'generated/gl_weight.cpp').write_text(cpp)
manifest = {'purpose': 'Local 110-lane copy; recurrence unchanged; fixed prime 1000000007 permits uint64 multiplication (<2^64) and constant modular reduction; original files untouched.',
            'sources': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (header_source, cpp_source)}}
(HERE/'grid_sources.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
subprocess.run(['clang++', '-std=c++20', '-O3', '-DNDEBUG',
                '-I'+str(HERE/'generated/include'), '-I'+str(ROOT/'cpp/include'),
                '-I'+str(BOOST), str(ROOT/'cpp/apps/chromatic_grid.cpp'),
                str(HERE/'generated/gl_weight.cpp'), str(ROOT/'cpp/src/chord_diagram.cpp'),
                '-o', str(HERE/'grid_probe')], check=True)
