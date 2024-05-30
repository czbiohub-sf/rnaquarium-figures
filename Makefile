SHELL:=/bin/bash
DATAPROMPT:=data source directory to link: 

# bold teal
T_HL:=[36;1m
# reset
T_0:=[m

.PHONY: data lint lint_py lint_sh

data:
#	@:don't echo, readline completion  bold teal prompt,                 initial path  store to
	@read         -e                   -p '$(T_HL)$(DATAPROMPT)$(T_0)' -i `pwd`        DATASRC; \
	ln -s "$$DATASRC" data


lint_py:
	-ruff check

lint_sh:
	-shellcheck --extended-analysis=true -o all -e SC2312 -S style scripts/*.sh

lint: lint_py lint_sh
