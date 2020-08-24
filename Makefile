# Debug option
DEBUG := 1

# Compiler selection
CC := gcc
# CC := icc
CXX := g++
#CXX := icpc

# Paths to directories
INCLUDE_DIR := ./include
#INCLUDES := $(addprefix -I,$(INCLUDES))
OBJ_DIR := ./obj
OBJ_SUB_DIR := $(patsubst $(INCLUDE_DIR)/%, ./obj/%, \
$(shell find $(INCLUDE_DIR) -type d))
SRC_DIR := ./src
EXE_DIR = .

# Files to compile
# SRC0 := $(wildcard $(SRC_DIR)/*.cc)
DRIVER := driver
DRIVER_WANG := driver_wang
DRIVER_BB := driver_bb
DRIVER_BB2 := driver_bb2
DRIVER_BB3 := driver_bb3
DRIVER_BB5 := driver_bb5
DRIVER_BB6 := driver_bb6
DRIVER_BB7 := driver_bb7
DRIVER_BB8 := driver_bb8
DRIVER_BB9 := driver_bb9
DRIVER_BB10 := driver_bb10
EXE := $(DRIVER) $(DRIVER_WANG) $(DRIVER_BB) $(DRIVER_BB2) $(DRIVER_BB3) $(DRIVER_BB5) $(DRIVER_BB6) $(DRIVER_BB7) $(DRIVER_BB8) $(DRIVER_BB9) $(DRIVER_BB10)
SRC_C := $(shell find . -name "*.c" -print)
SRC_CXX := $(shell find . -name "*.cpp" -print)
SRC := $(SRC_C) $(SRC_CXX)
DEPENDS := $(SRC_C:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.d) \
$(SRC_CXX:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.d)
OBJ_C := $(SRC_C:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
OBJ_CXX := $(SRC_CXX:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
OBJ := $(OBJ_C) $(OBJ_CXX)
OBJ_DRIVER := $(filter $(OBJ_DIR)/$(DRIVER)/%,$(OBJ))
OBJ_DRIVER_WANG := $(filter $(OBJ_DIR)/$(DRIVER_WANG)/%,$(OBJ))
OBJ_DRIVER_BB := $(filter $(OBJ_DIR)/$(DRIVER_BB)/%,$(OBJ))
OBJ_DRIVER_BB2 := $(filter $(OBJ_DIR)/$(DRIVER_BB2)/%,$(OBJ))
OBJ_DRIVER_BB3 := $(filter $(OBJ_DIR)/$(DRIVER_BB3)/%,$(OBJ))
OBJ_DRIVER_BB5 := $(filter $(OBJ_DIR)/$(DRIVER_BB5)/%,$(OBJ))
OBJ_DRIVER_BB6 := $(filter $(OBJ_DIR)/$(DRIVER_BB6)/%,$(OBJ))
OBJ_DRIVER_BB7 := $(filter $(OBJ_DIR)/$(DRIVER_BB7)/%,$(OBJ))
OBJ_DRIVER_BB8 := $(filter $(OBJ_DIR)/$(DRIVER_BB8)/%,$(OBJ))
OBJ_DRIVER_BB9 := $(filter $(OBJ_DIR)/$(DRIVER_BB9)/%,$(OBJ))
OBJ_DRIVER_BB10 := $(filter $(OBJ_DIR)/$(DRIVER_BB10)/%,$(OBJ))
#$(DRIVER): $(OBJ) := $(filter-out $(OBJ_DIR)/$(DRIVER_WANG)/%,$(OBJ))
#$(DRIVER_WANG): $(OBJ) := $(filter-out $(OBJ_DIR)/$(DRIVER)/%,$(OBJ))
CLEAN_LOG_TARGETS := $(OBJ_DIR)/*.log $(EXE:%=%.log)
CLEAN_TARGETS := $(OBJ) $(OBJ_DIR)
DISTCLEAN_TARGETS := $(CLEAN_TARGETS) $(EXE)

# Compiler options
CPPFLAGS :=
$(OBJ_DIR)/%.o: CPPFLAGS += -MMD -MP
CPPFLAGS += $(patsubst %,-I%,$(shell find $(INCLUDE_DIR) -type d))
ifeq ($(DEBUG),0)
	ifeq ($(CC),gcc)
		CFLAGS := -Wall -Wextra -std=c11 -O3 -fopenmp -mavx2
		CXXFLAGS := -Wall -Wextra -std=c++11 -O3 -fopenmp -mavx2
	else ifeq ($(CC),icc)
		CFLAGS := -Wall -Wextra -std=c11 -O3 -qopenmp -xCORE-AVX2
		CXXFLAGS := -Wall -Wextra -std=c++11 -O3 -qopenmp -xCORE-AVX2
	endif
else
	ifeq ($(CC),gcc)
		# CFLAGS := -Wall -Wextra -std=c11 -g -O3 -fopenmp -mavx2
		# CXXFLAGS := -Wall -Wextra -std=c++11 -g -O3 -fopenmp -mavx2
		CFLAGS := -Wall -Wextra -std=c11 -g -O2 -ftree-vectorize -fopenmp -mavx2
		CXXFLAGS := -Wall -Wextra -std=c++11 -g -O2 -ftree-vectorize -fopenmp -mavx2
	else ifeq ($(CC),icc)
		CFLAGS := -Wall -Wextra -std=c11 -g -O3 -qopenmp -xCORE-AVX2
		CXXFLAGS := -Wall -Wextra -std=c++11 -g -O3 -qopenmp -xCORE-AVX2
	endif
endif
CPPFLAGS += -D'CFLAGS="$(CC) $(CFLAGS)"' -D'CXXFLAGS="$(CXX) $(CXXFLAGS)"'

# Linker
LDFLAGS :=
ifeq ($(CC),gcc)
	LDLIBS := -fopenmp
else ifeq ($(CC),icc)
	LDLIBS := -qopenmp
endif
LDLIBS += -lm

# Intel MKL
# USE_MKL :=
ifdef USE_MKL
ifeq ($(CC),gcc)
	CPPFLAGS += -DHAVE_MKL=\"$(USE_MKL)\" -m64 -I${MKLROOT}/include
	LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed
	LDLIBS += -lmkl_intel_lp64

	# Parallel MKL using TBB
	ifeq ($(USE_MKL),tbb)
		LDLIBS += -lmkl_tbb_thread -lmkl_core -ltbb -lstdc++
	endif # Parallel MKL using TBB

	# Parallel MKL using Gomp
	ifeq ($(USE_MKL),gomp)
		LDLIBS += -lmkl_gnu_thread -lmkl_core -lgomp
	endif # Parallel MKL using Gomp

	# Parallel MKL using iomp5
	ifeq ($(USE_MKL),iomp)
		LDLIBS += -lmkl_intel_thread -lmkl_core -liomp5
	endif # Parallel MKL using iomp5

	# Sequential MKL
	ifeq ($(USE_MKL),serial)
		LDLIBS += -lmkl_sequential -lmkl_core
	endif # Sequential MKL

	LDLIBS += -lpthread -lm -ldl
else ifeq ($(CC),icc)

	CPPFLAGS += -DHAVE_MKL
	# Parallel MKL using TBB
	ifeq ($(USE_MKL),tbb)
		CPPFLAGS += -mkl=parallel
		LDLIBS += -mkl=parallel -ltbb -lstdc++ -lpthread -lm -ldl
	endif # Parallel MKL using TBB

	# Parallel MKL using iomp5
	ifeq ($(USE_MKL),iomp)
		CPPFLAGS += -mkl=parallel
		LDLIBS += -mkl=parallel -liomp5 -lpthread -lm -ldl
	endif # Parallel MKL using iomp5

	# Sequential MKL
	ifeq ($(USE_MKL),sequential)
		CPPFLAGS += -mkl=sequential -DHAVE_MKL_SEQUENTIAL
		LDLIBS += -mkl=sequential -lpthread -lm -ldl
	endif # Sequential MKL

endif
endif # Intel MKL

# TBB
# USE_TBB :=
ifdef USE_TBB
CPPFLAGS += -DHAVE_TBB
LDLIBS += -ltbb -lstdc++ -lm -ldl
endif # TBB

# Colors
COMPILE_COLOR := \033[0;34m
LINK_COLOR := \033[0;35m
INFO_COLOR := \033[0;36m
OK_COLOR := \033[1;32m
ERROR_COLOR := \033[1;31m
WARN_COLOR := \033[1;33m
WHITE := \033[m
B_WHITE := \033[0;1m

# Messages
ARROW := "$(B_WHITE)├$(WHITE)"
ARROW_END := "$(B_WHITE)└$(WHITE)"
OK_STRING := "$(OK_COLOR)[OK]$(WHITE)\n"
ERROR_STRING := "$(ERROR_COLOR)[ERROR]$(WHITE)\n"
WARN_STRING := "$(WARN_COLOR)[WARNING]$(WHITE)\n"
CLEAN_STRING := $(ERROR_COLOR)Cleaning$(WHITE)
COMPILE_STRING = "$(COMPILE_COLOR)Compiling$(WHITE) $(@F)"
LINK_STRING = "$(LINK_COLOR)Linking$(WHITE) $(@F)"

# Functions
RM := rm -rf
RM_VERBOSE := $(RM) -v
VALGRIND = valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all

define execute
@printf "%-84b" $(2);\
if [ $(DEBUG) -eq 1 ]; then\
	printf "\n%b " $(ARROW);\
	echo $(1);\
	$(1) 2> $@.log;\
	RESULT=$$?;\
	if [ $$RESULT -ne 0 ]; then \
		printf "%b %b" $(ARROW_END) $(ERROR_STRING); \
	elif [ -s $@.log ]; then \
		printf "%b %b" $(ARROW_END) $(WARN_STRING); \
	else \
		printf "%b %b" $(ARROW_END) $(OK_STRING); \
	fi;\
	cat -n $@.log;\
else \
	$(1) 2> $@.log;\
	RESULT=$$?;\
	if [ $$RESULT -ne 0 ]; then \
		printf "  %b" $(ERROR_STRING); \
	elif [ -s $@.log ]; then \
		printf "%b" $(WARN_STRING); \
	else \
		printf "     %b" $(OK_STRING); \
	fi; \
	cat -n $@.log;\
	rm -f $@.log;\
fi;\
exit $$RESULT
endef # execute

define clean_files
@printf "%-84b" $(2);\
if [ $(DEBUG) -eq 1 ]; then\
	printf "\n%b %s\n" $(ARROW) "$(RM_VERBOSE) $(1)";\
	$(RM_VERBOSE) $(1) | column -t 2> $@.log;\
	RESULT=$$?;\
	if [ $$RESULT -ne 0 ]; then \
		printf "%b %b" $(ARROW_END) $(ERROR_STRING); \
	elif [ -s $@.log ]; then \
		printf "%b %b" $(ARROW_END) $(WARN_STRING); \
	else \
		printf "%b %b" $(ARROW_END) $(OK_STRING); \
	fi;\
	cat -n $@.log;\
	rm -f $@.log;\
else \
	$(RM) $(1) 2> $@.log;\
	RESULT=$$?;\
	if [ $$RESULT -ne 0 ]; then \
		printf "  %b" $(ERROR_STRING); \
	elif [ -s $@.log ]; then \
		printf "%b" $(WARN_STRING); \
	else \
		printf "     %b" $(OK_STRING); \
	fi; \
	cat -n $@.log;\
	rm -f $@.log;\
fi;\
exit $$RESULT
endef # clean_files

# Targets
.PHONY: all clean_log clean distclean

all: $(EXE)
	@printf "%b\n" "$(OK_COLOR)Successful compilation$(WHITE)"

benchmarks: all
	@./$(DRIVER) ; ./$(DRIVER_WANG)

docs: all
	@mkdir -p ./docs;
	$(call execute,doxywizard Doxyfile,"Making documentation files")

memcheck_driver: $(DRIVER)
	@$(VALGRIND) ./$^

memcheck_driver_wang: $(DRIVER_WANG)
	@$(VALGRIND) ./$^

memcheck_driver_bb: $(DRIVER_BB)
	@$(VALGRIND) ./$^

memcheck_driver_bb2: $(DRIVER_BB2)
	@$(VALGRIND) ./$^

memcheck_driver_bb3: $(DRIVER_BB3)
	@$(VALGRIND) ./$^

memcheck_driver_bb5: $(DRIVER_BB5)
	@$(VALGRIND) ./$^

memcheck_driver_bb6: $(DRIVER_BB6)
	@$(VALGRIND) ./$^

memcheck_driver_bb7: $(DRIVER_BB7)
	@$(VALGRIND) ./$^

memcheck_driver_bb8: $(DRIVER_BB8)
	@$(VALGRIND) ./$^

memcheck_driver_bb9: $(DRIVER_BB9)
	@$(VALGRIND) ./$^

memcheck_driver_bb10: $(DRIVER_BB10)
	@$(VALGRIND) ./$^

$(DRIVER): $(OBJ_DRIVER)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_WANG): $(OBJ_DRIVER_WANG)
	$(call execute,$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB): $(OBJ_DRIVER_BB)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB2): $(OBJ_DRIVER_BB2)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB3): $(OBJ_DRIVER_BB3)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB5): $(OBJ_DRIVER_BB5)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB6): $(OBJ_DRIVER_BB6)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB7): $(OBJ_DRIVER_BB7)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB8): $(OBJ_DRIVER_BB8)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB9): $(OBJ_DRIVER_BB9)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

$(DRIVER_BB10): $(OBJ_DRIVER_BB10)
	$(call execute,$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@,$(LINK_STRING))

# $(OBJ): | $(OBJ_DIR) $(OBJ_SUB_DIR)

$(OBJ_DIR) $(OBJ_SUB_DIR):
	@mkdir -p $@

-include $(DEPENDS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c  | $(OBJ_DIR) $(OBJ_SUB_DIR)
	$(call execute,$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<,$(COMPILE_STRING))

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR) $(OBJ_SUB_DIR)
	$(call execute,$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ -c $<,$(COMPILE_STRING))

clean_log:
	$(call clean_files,$(CLEAN_LOG_TARGETS),"$(CLEAN_STRING) log files")

clean: clean_log
	$(call clean_files,$(CLEAN_TARGETS),"$(CLEAN_STRING) object and dependency \
	files")

distclean: clean
	$(call clean_files,$(DISTCLEAN_TARGETS),"$(CLEAN_STRING) executable files")
