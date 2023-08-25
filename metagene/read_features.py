#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-03-08 20:42


import os
import pickle

# Get cache directory
if "XDG_CACHE_HOME" in os.environ:
    cache_dir = os.path.join(os.environ["XDG_CACHE_HOME"], "metagene")
else:
    home_dir = os.path.expanduser("~")
    cache_dir = os.path.join(home_dir, ".metagene")

# Check if cache directory exists and create it if it doesn't
if not os.path.exists(cache_dir):
    os.makedirs(cache_dir)

# Define cache file name
cache_file = os.path.join(cache_dir, "ref.pickle")


def parse_gtf(input_file):
    return input_file


def parse():
    # Check if cache file exists and load it if it does
    if os.path.exists(cache_file):
        with open(cache_file, "rb") as f:
            ref = pickle.load(f)
    else:
        # Process input file
        ref = parse_gtf("ref.gtf")

        # Save processed data to cache file
        with open(cache_file, "wb") as f:
            pickle.dump(ref, f)


if __name__ == "__main__":
    parse()
