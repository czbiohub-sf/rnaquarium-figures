SHELL := /bin/bash
DATAPROMPT :=data source directory to link: 

# bold teal
T_HL :=[36;1m
# reset
T_0 :=[0m

.PHONY: data

data:
#	@:don't echo, readline completion  bold teal prompt,                 initial path  store to
	@read         -e                   -p '$(T_HL)$(DATAPROMPT)$(T_0)' -i `pwd`      DATASRC
	ln -s $DATASRC data
