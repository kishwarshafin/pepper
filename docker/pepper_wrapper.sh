#!/usr/bin/env bash

set -e

# Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

# RUN HELEN SCRIPT
python3 /opt/pepper/$@