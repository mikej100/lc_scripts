#!/bin/bash
remote=$(git -C $(dirname ${BASH_SOURCE}) remote -v | grep -m 1 -o 'http[^ ]*')
branch=$(git -C $(dirname ${BASH_SOURCE}) rev-parse --abbrev-ref HEAD)
hash=$(git -C $(dirname ${BASH_SOURCE}) show -s --format=%H)
info=$(git -C $(dirname ${BASH_SOURCE}) show -s --format=" Date: %ci. Author: %an. %nTitle: %s")

echo
echo "Script version control: latest commit information"
echo "github repository: ${remote}/tree/${hash}"
echo "Branch: $branch $info"
echo