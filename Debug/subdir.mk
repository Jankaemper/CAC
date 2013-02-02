################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../MyWindow_files.cpp \
../MyWindow_operations.cpp \
../MyWindow_overlay.cpp \
../arrangement_2.cpp \
../cartograms.cpp \
../demo_tab.cpp \
../forms.cpp \
../qt_layer.cpp 

OBJS += \
./MyWindow_files.o \
./MyWindow_operations.o \
./MyWindow_overlay.o \
./arrangement_2.o \
./cartograms.o \
./demo_tab.o \
./forms.o \
./qt_layer.o 

CPP_DEPS += \
./MyWindow_files.d \
./MyWindow_operations.d \
./MyWindow_overlay.d \
./arrangement_2.d \
./cartograms.d \
./demo_tab.d \
./forms.d \
./qt_layer.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


