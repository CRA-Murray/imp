#!/usr/bin/env python

"""
   Make an alias for the kernel module headers.
"""

import tools
import os

def main():
    tools.link(os.path.join('include', 'IMP'),
               os.path.join('include', 'IMP', 'kernel'))

if __name__ == '__main__':
    main()
