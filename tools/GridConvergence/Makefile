# Makefile

CC = gcc
GFLAGS = -lm

SOURCES = GridConvergence.c
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = GridConvergence

PROGRESS = \033[1;91mProgress:\033[0m

PERCENT = 0.01

execname: $(SOURCES)
	@# @rm -rf Make
	@# @mkdir Make
	@# @$(CC) $(GFLAGS) $(SOURCES_MODULES)
	@# @$(CC) $(GFLAGS) $(SOURCES_SUBPROGRAMS)
	@# @$(CC) $(GFLAGS) $(SOURCES_MAIN)
	@# @$(CC) $(OBJECTS) -o $(EXECUTABLE)
	@# @mv *.o Make
	@# @mv *.mod Make
	

	@echo "\033[1mStart compiling"
	@echo -ne '$(PROGRESS) \033[1m  ..........                                                                                              \033[91m[10%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ...............                                                                                         \033[91m[10%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ....................                                                                                    \033[91m[20%]\033[0m\r'
	@sleep $(PERCENT)
		
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

	@echo -ne '$(PROGRESS) \033[1m  ..................................................                                                      \033[91m[50%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  .......................................................                                                 \033[91m[50%]\033[0m\r'
	@sleep $(PERCENT)
		@$(CC) $(SOURCES) $(GFLAGS) -o $(EXECUTABLE) 
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

	@echo -ne '$(PROGRESS) \033[1m  .....................................................................................                   \033[91m[80%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ..........................................................................................              \033[91m[90%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ...............................................................................................         \033[91m[90%]\033[0m\r'
	@sleep $(PERCENT)

	@echo -ne '$(PROGRESS) \033[1m  ....................................................................................................    \033[91m[100%]\033[0m\r'
	@echo -ne '\n'
	@echo "\033[1mExecutable created"
	@echo "Compilation was successful!"

clean:
	rm -rf *.o *.mod $(EXECUTABLE)
	@echo "\033[1mFiles removed!"