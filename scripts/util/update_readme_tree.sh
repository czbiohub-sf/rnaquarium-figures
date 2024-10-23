#!/bin/bash

PAT_START='\[comment\]: \# {REPO_TREE}\n\`\`\`\n'

IGNORE='__init__.py
pyproject.toml
environment-dev.yml
scripts
tests
scratch
intermediate'

# construct directory ignore parameter list from .gitignore entries
GITIGNORE=$(grep -vE '^(#.*|\s*)$' .gitignore | sed 's/\(.*\)/-I \1/g' | paste -s -d" ")
OURIGNORE=$(echo "${IGNORE}" | sed 's/\(.*\)/-I \1/g' | paste -s -d" ")

# call tree and escape forward slashes for sed
REPOTREE=$(tree -n --noreport ${GITIGNORE} ${OURIGNORE} | sed 's/\//\\\//g')
if [[ "$1" != "-q" ]]; then echo "${REPOTREE}"; fi

# replace between code block after {REPO_TREE} with the repository tree contents
perl -pni -e "BEGIN{undef $/;} s/${PAT_START}[^\`]*\`\`\`/${PAT_START}${REPOTREE}\n\`\`\`/g" README.md
