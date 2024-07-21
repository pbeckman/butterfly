#!/usr/bin/env python

import numpy as np
import sys

from pathlib import Path

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'usage: {sys.argv[0]} <path.bin>')
        sys.exit(0)

    np.fromfile(Path(sys.argv[1]), dtype=np.int32).astype(
