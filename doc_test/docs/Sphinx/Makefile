SRC_DIR = src
SPHINX_MK_FILE=sphinx.mk
RM=rm
RMFLAGS=-rf

.PHONY: all
all : html

.PHONY: html
html :
	@make -f $(SPHINX_MK_FILE) html

.PHONY: clean
clean :
	$(RM) $(RMFLAGS) source/_doxygen/ source/api/
	@make  -f $(SPHINX_MK_FILE) clean

.PHONY: distclean
distclean : clean
	$(RM) $(RMFLAGS) build
