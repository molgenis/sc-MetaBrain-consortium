#!/bin/bash

which zenodo_get >/dev/null || { echo "zenodo_get: command not found"; exit 1; }

zenodo_get -d 10.5281/zenodo.7276971