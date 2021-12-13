DOXYGEN =doxygen
DOXY_FILE=doxyfile.cfg
DOXY_CFG_FILE=doxyfile-prj.cfg

all : $(DOXY_FILE) $(DOXY_CFG_FILE)
	$(DOXYGEN) $(DOXY_CFG_FILE)

$(DOXY_FILE) :
	$(DOXYGEN) -g $@
