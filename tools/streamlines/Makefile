# Makefile

CC = gfortran
GFLAGS = -O3 -ffree-form -Wno-unused-parameter -ffree-line-length-0 -fdefault-double-8 -fdefault-real-8 -c

SOURCES_MODULES = global_variables.f
SOURCES_SUBPROGRAMS =  subprograms.f
SOURCES_MAIN = main.f
SOURCES = $(SOURCES_MODULES) $(SOURCES_SUBPROGRAMS) $(SOURCES_MAIN)
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = streamline

PROGRESS = \033[1;91mProgress:\033[0m

PERCENT = 0.01

execname: $(SOURCES)
	@# @rm -rf Make
	@# @mkdir Make
	@# @$(CC) $(GFLAGS) $(SOURCES_MODULES)
	@# @$(CC) $(GFLAGS) $(SOURCES_CREATE)
	@# @$(CC) $(OBJECTS_CREATE) -o $(EXECUTABLE_CREATE) 
	@# @./$(EXECUTABLE_CREATE)
	@# @rm -rf *.ο $(EXECUTABLE_CREATE)
	@# @$(CC) $(GFLAGS) $(SOURCES_MODULES)
	@# @$(CC) $(GFLAGS) $(SOURCES_SUBPROGRAMS)
	@# @$(CC) $(GFLAGS) $(SOURCES_MAIN)
	@# @$(CC) $(OBJECTS) -o $(EXECUTABLE)
	@# @mv *.o Make
	@# @mv *.mod Make
	

	@echo "\033[1mStart compiling"
	@echo -ne '$(PROGRESS) \033[1m  ..........                                                                                              \033[91m[10%]\033[0m\r'
	@sleep $(PERCENT)
		@rm -rf Make
	@echo -ne '$(PROGRESS) \033[1m  ...............                                                                                         \033[91m[10%]\033[0m\r'
	@sleep $(PERCENT)
		@mkdir Make
	@echo -ne '$(PROGRESS) \033[1m  ....................                                                                                    \033[91m[20%]\033[0m\r'
	@sleep $(PERCENT)
		@$(CC) $(GFLAGS) $(SOURCES_MODULES)
	@echo -ne '$(PROGRESS) \033[1m  .........................                                                                               \033[91m[20%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ..............................                                                                          \033[91m[30%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ...................................                                                                     \033[91m[30%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ........................................                                                                \033[91m[40%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  .............................................                                                           \033[91m[40%]\033[0m\r'
	@sleep $(PERCENT)
		@$(CC) $(GFLAGS) $(SOURCES_MODULES)
	@echo -ne '$(PROGRESS) \033[1m  ..................................................                                                      \033[91m[50%]\033[0m\r'
	@sleep $(PERCENT)
		@$(CC) $(GFLAGS) $(SOURCES_SUBPROGRAMS)
	@echo -ne '$(PROGRESS) \033[1m  .......................................................                                                 \033[91m[50%]\033[0m\r'
	@sleep $(PERCENT)
		@$(CC) $(GFLAGS) $(SOURCES_MAIN)
	@echo -ne '$(PROGRESS) \033[1m  ............................................................                                            \033[91m[60%]\033[0m\r'
	@sleep $(PERCENT)
	
	@echo -ne '$(PROGRESS) \033[1m  .................................................................                                       \033[91m[60%]\033[0m\r'
	@sleep $(PERCENT)
	
	@echo -ne '$(PROGRESS) \033[1m  ......................................................................                                  \033[91m[70%]\033[0m\r'
	@sleep $(PERCENT)
	
	@echo -ne '$(PROGRESS) \033[1m  ...........................................................................                             \033[91m[70%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ................................................................................                        \033[91m[80%]\033[0m\r'
	@sleep $(PERCENT)
		@$(CC) $(OBJECTS) -o $(EXECUTABLE)
	@echo -ne '$(PROGRESS) \033[1m  .....................................................................................                   \033[91m[80%]\033[0m\r'
	@sleep $(PERCENT)
	@echo -ne '$(PROGRESS) \033[1m  ..........................................................................................              \033[91m[90%]\033[0m\r'
	@sleep $(PERCENT)
		@mv *.o Make
	@echo -ne '$(PROGRESS) \033[1m  ...............................................................................................         \033[91m[90%]\033[0m\r'
	@sleep $(PERCENT)
		@mv *.mod Make
	@echo -ne '$(PROGRESS) \033[1m  ....................................................................................................    \033[91m[100%]\033[0m\r'
	@echo -ne '\n'
	@echo "\033[1mExecutable created"
	@echo "Compilation was successful!"

clean:
	rm -rf *.o *.mod $(EXECUTABLE)
	rm -rf Make
	rm -rf ../results
	rm -rf ../mesh
	@echo "\033[1mFiles removed!"
