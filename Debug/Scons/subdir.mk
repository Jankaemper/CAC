################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../Scons/tester2.o 

CPP_SRCS += \
../Scons/tester.cpp \
../Scons/tester2.cpp 

OBJS += \
./Scons/tester.o \
./Scons/tester2.o 

CPP_DEPS += \
./Scons/tester.d \
./Scons/tester2.d 


# Each subdirectory must supply rules for building sources it contributes
Scons/%.o: ../Scons/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


