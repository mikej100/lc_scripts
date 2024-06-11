#!/bin/bash
echo
echo "Script file information from git"
echo "================================"

echo "Here are repository; branch; most recent commit reference, date, author and comment"

#git -C "$SCRIPT_DIR" remote -v | grep fetch
git -C $(dirname ${BASH_SOURCE}) rev-parse --abbrev-ref HEAD
git -C $(dirname ${BASH_SOURCE}) show -s --format=%h%x09%ci%x09%an%x09%sk