"""
Run tests for st_exact_deRham (Paper A1).

Usage:
    python3 run_tests.py                    # all tests
    python3 run_tests.py setup              # only setup
    python3 run_tests.py recurrence voronoi # only these two

Dependencies: numpy, scipy
"""
import sys
import os
import time

# Setup: add src/ and tests/ to path so imports work
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'tests'))

# Suppress OpenBLAS threading for reproducibility
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')

# All registered modules: name -> (module_path, has_runtime)
ALL_MODULES = {
    # src internal tests
    'structures':      ('core_math.spec.structures', False),
    'incidence':       ('core_math.operators.incidence', False),
    'kelvin':          ('core_math.builders.kelvin', False),
    'c15':             ('core_math.builders.c15_periodic', True),
    'wp_periodic':     ('core_math.builders.weaire_phelan_periodic', True),
    'solids':          ('core_math.builders.solids', False),
    'solids_periodic': ('core_math.builders.solids_periodic', False),
    'polyhedra':       ('core_math.builders.polyhedra', False),
    'hodge':           ('physics.hodge', False),
    'bloch':           ('physics.bloch', False),
    'gauge_bloch':     ('physics.gauge_bloch', False),
    # paper tests (Groups 1-6)
    'setup':           ('test_01_setup', False),
    'standard_fails':  ('test_02_standard_fails', False),
    'recurrence':      ('test_03_recurrence', False),
    'spectral_gains':  ('test_04_spectral_gains', False),
    'inexactness':     ('test_05_inexactness', False),
    'voronoi':         ('test_06_voronoi', False),
    'curvature':       ('test_07_curvature', False),
}

# Filter: positional args that aren't flags
filters = [a for a in sys.argv[1:] if not a.startswith('--')]


def run_module_tests(module_path, has_runtime=False):
    """Import module and run its _test_init() or main()."""
    mod = __import__(module_path, fromlist=['_test_init', 'main'])
    if hasattr(mod, '_test_init'):
        mod._test_init()
    elif hasattr(mod, 'main'):
        mod.main()
    if has_runtime and hasattr(mod, '_test_runtime'):
        mod._test_runtime()


if __name__ == '__main__':
    # Determine which modules to run
    if filters:
        to_run = {}
        for f in filters:
            matched = [k for k in ALL_MODULES if f in k]
            if not matched:
                print(f"Unknown module: {f}. Available: {list(ALL_MODULES.keys())}")
                sys.exit(1)
            for m in matched:
                to_run[m] = ALL_MODULES[m]
    else:
        to_run = ALL_MODULES

    print("=" * 60)
    print(f"st_exact_deRham — tests ({', '.join(to_run.keys())})")
    print("=" * 60)
    t0 = time.time()

    for name, (module_path, has_rt) in to_run.items():
        run_module_tests(module_path, has_runtime=has_rt)

    elapsed = time.time() - t0
    print("\n" + "=" * 60)
    print(f"DONE ({elapsed:.1f}s)")
    print("=" * 60)
