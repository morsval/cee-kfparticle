.DEFAULT_GOAL := all
.PHONY: all lite analysis legacy check clean help
all lite analysis legacy check clean help:
	+$(MAKE) -C CEErco $@
