#!/bin/bash

which zenodo_get >/dev/null || { echo "zenodo_get: command not found"; exit 1; }

zenodo_get -d 10.5281/zenodo.7276971

wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-022-01128-z/MediaObjects/41593_2022_1128_MOESM3_ESM.xlsx