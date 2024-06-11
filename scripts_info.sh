#!/bin/bash
echo
echo "Script version control: latest commit information"
remote=$(git -C $(dirname ${BASH_SOURCE}) remote -v | grep -m 1 -o 'http.*\s')
branch=$(git -C $(dirname ${BASH_SOURCE}) rev-parse --abbrev-ref HEAD)
echo "github repository: ${remote} Branch: ${branch}"
git -C $(dirname ${BASH_SOURCE}) show -s --format="Short hash: %h. Date: %ci. Author: %an. %nTitle: %s"
echo