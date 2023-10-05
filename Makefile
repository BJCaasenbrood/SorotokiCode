PYTHON = python2
MKDOCS = mkdocs
MATLAB_FILES = interpolateSE3.m

SRC = $(filter $(MATLAB_FILES), $(wildcard *.m))
TAR = $(SRC:.m=.md)

mfiledir = docs/mfiles

$(info SRC is $(SRC))
$(info TAR is $(TAR))
$(info mfiledir is $(mfiledir))

.PHONY: all clean

all: $(TAR)

$(MATLAB_FILES:.m=.md): $(MATLAB_FILES) docs/matdoc.py docs/matdocparser.py
	$(info $@)
	$(info $(@D))
	mkdir -p $(mfiledir)
	$(PYTHON) ./docs/matdoc.py "$<" > "$(mfiledir)/$@"

doc-serve: mkdocs.yml
	$(MKDOCS) serve

clean:
	rm -f $(TAR)
	rm -rf $(mfiledir)