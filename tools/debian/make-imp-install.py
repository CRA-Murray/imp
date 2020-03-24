#!/usr/bin/python3

import sys
import os

all_modules = sorted(d for d in os.listdir('.') if os.path.isdir(d))
exclude_modules = frozenset(('mpi', 'spb', 'scratch', 'kernel'))

non_mpi_modules = [m for m in all_modules if m not in exclude_modules]
data_modules = [m for m in all_modules
                if os.path.exists(os.path.join(m, 'data'))]

non_mpi_binaries = []
for m in non_mpi_modules + ['rmf/dependency/RMF']:
    bindir = os.path.join(m, 'bin')
    if os.path.exists(bindir):
        for b in os.listdir(bindir):
            fname, ext = os.path.splitext(b)
            if not fname.startswith('__') and ext in ('', '.cpp', '.py'):
                non_mpi_binaries.append(fname)

def subst_line(line, template, modules):
    for m in modules:
        sys.stdout.write(line.replace(template, m))

for line in sys.stdin:
    if "@NON_MPI_MODULES@" in line:
        subst_line(line, '@NON_MPI_MODULES@', non_mpi_modules)
    elif "@DATA_MODULES@" in line:
        subst_line(line, '@DATA_MODULES@', data_modules)
    elif "@NON_MPI_BINARIES@" in line:
        subst_line(line, '@NON_MPI_BINARIES@', non_mpi_binaries)
    else:
        sys.stdout.write(line)