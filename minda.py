#!/usr/bin/env python3

#(c) 2023 by Authors
#This file is a part of Minda.
#Released under the BSD license (see LICENSE file)

"""
This script sets up environment paths
and invokes Minda without installation.
"""

import os
import sys

def main():
    #Setting executable paths
    minda_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, minda_root)

    #Minda entry point
    from minda.main import main
    sys.exit(main())


if __name__ == "__main__":
    main()