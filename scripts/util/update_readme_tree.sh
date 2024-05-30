#!/bin/bash

PAT_START='\[comment\]: \# {REPO_TREE}\n\`\`\`\n'

# construct directory ignore parameter list from .gitignore entries
TREEARGS=$(grep -vE '^(#.*|\s*)$' .gitignore | sed 's/\(.*\)/-I \1/g' | paste -s -d" ")
# call tree and escape forward slashes for sed
REPOTREE=$(tree -n --noreport ${TREEARGS} | sed 's/\//\\\//g')
if [[ "$1" != "-q" ]]; then echo "${REPOTREE}"; fi

# replace between code block after {REPO_TREE} with the repository tree contents
perl -pni -e "BEGIN{undef $/;} s/${PAT_START}[^\`]*\`\`\`/${PAT_START}${REPOTREE}\n\`\`\`/g" README.md