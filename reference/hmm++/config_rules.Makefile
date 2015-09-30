########################
#Rules
########################

# Ensure the build sub-directories exist:
$(ALL_BLD_DIRS):
	mkdir -p $@

define ENSURE_BUILD_DIRS_PREREQS
| $(ALL_BLD_DIRS)
endef


#Replace implicit rules:
$(BLD_DIR)/%.o: $(SRC_DIR)/%.cpp $(ENSURE_BUILD_DIRS_PREREQS)
	$(CXX) -c $(ALL_CXXFLAGS) $< -o $@

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c $(ENSURE_BUILD_DIRS_PREREQS)
	$(GCC) -c $(ALL_CFLAGS)   $< -o $@


#explicit rules for making libraries and files:
$(LIB_FULL_PATH): $(LIB_OBJECTS_FULL_PATH) $(ENSURE_BUILD_DIRS_PREREQS)
	ar rv $@ $?
	ranlib $@


#depend rules:
define MAKE_DEPEND_PREREQS
./Makefile $(DEPEND_HEADERS) $(ENSURE_BUILD_DIRS_PREREQS)
endef

define MAKE_DEPEND_RULE
	@echo making depend $<; \
	$(CXX) -MM $(ALL_CXXFLAGS) -MT $(BLD_DIR)/$*.o $< > $@.$$$$; \
	sed 's,\(.*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$
endef

$(BLD_DIR)/%.dep: $(SRC_DIR)/%.cpp $(MAKE_DEPEND_PREREQS)
	$(MAKE_DEPEND_RULE)

$(BLD_DIR)/%.dep: $(SRC_DIR)/%.c   $(MAKE_DEPEND_PREREQS)
	$(MAKE_DEPEND_RULE)


#maintenance rules:
.PHONY: clean
clean:
	$(RM) -r $(REMOVES)


.PHONY: zip
zip:
	$(RM) $(MAKE_NAME).tgz && \
	tar -cvzf $(MAKE_NAME).tgz \
	--exclude=*.o --exclude=*.dep --exclude=*~ --exclude=*.tgz --exclude=*.a --exclude=*.snapshot \
	Makefile $(SRC_DIR) $(ZIP_INCLUDES) && \
	$(RM) -r tmp_zip && mkdir tmp_zip && cd tmp_zip && mkdir $(MAKE_NAME) && cd $(MAKE_NAME) && tar -xvzf ../../$(MAKE_NAME).tgz && \
	cd ../ && tar -cvzf ../$(MAKE_NAME).tgz $(MAKE_NAME) && cd ../ && $(RM) -r tmp_zip


FORCE:



# Include generated dependency files:
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),zip)
include $(DEPENDENCIES_FILES)
endif
endif
