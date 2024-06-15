#!/bin/bash
remote=$(git -C $(dirname ${BASH_SOURCE}) remote -v | grep -m 1 -o '.*' )
link_stem=$(echo ${remote} | sed -r 's/origin (http.*)\.git \(fetch\)/\1/' )

branch=$(git -C $(dirname ${BASH_SOURCE}) rev-parse --abbrev-ref HEAD)
hash=$(git -C $(dirname ${BASH_SOURCE}) show -s --format=%H)
info=$(git -C $(dirname ${BASH_SOURCE}) show -s --format=" Date: %ci. Author: %an. %nTitle: %s")

echo
echo "Script version control: latest commit information"
echo "github link: ${link_stem}/tree/${hash}"
echo "Branch: $branch $info"
echo