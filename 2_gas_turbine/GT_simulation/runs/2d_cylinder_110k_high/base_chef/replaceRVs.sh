#!/bin/bash -l

### Boundary conditions are functions of random variables. Use conditions from constantine et al for now. 
pi=3.14159265

sigma_1=0.7
sigma_2=0.2
sigma_3=0.09

# read inputs from file
# run number is first input to identify correct line

Run_num="$1"

echo $Run_num
RAND_PARAM_FILE="setting_param_file.dat"

# Extract variables from file

line_string=$(cat setting_param_file.dat | sed -n "$((Run_num+1))""p")

IFS=', ' read -r -a line_array <<< "$line_string"

Y1="${line_array[0]}"
Y2="${line_array[1]}"
Y3="${line_array[2]}"
Y4="${line_array[3]}"

cos_theta="(-x1/sqrt(x1^2+x2^2))"

curve_max=100.0
curve_mid=`echo $Y4\*0.05 | bc -l`
curve_sd=0.01

#replace_inflow=0.00905 # Re 50
#replace_inflow=0.181 #Re 1000

vel_1=`echo $sigma_1\*$Y1 | bc -l`
vel_2=`echo 2\*$pi\*10 | bc -l`
vel_3=`echo $sigma_2\*$Y2 | bc -l`
vel_4=`echo 10\*$pi\*10 | bc -l`

replace_inflow_velocity="156.8*(1+($vel_1)*cos($vel_2*x2)+($vel_3)*cos($vel_4*x2))"

temp_inflow_1=`echo 1/\(2\*$curve_sd^2\) | bc -l`

replace_inflow_temperature="300+$curve_max*exp(-((x2-($curve_mid))^2)*$temp_inflow_1)"

replace_initial_temperature="300"

temp_cy=`echo 0.1+$sigma_3\*$Y3/2 | bc -l`
replace_cylinder="-50*exp(-$temp_cy*$cos_theta)"

char_inflow_velocity="${#replace_inflow_velocity} "
char_inflow_temperature="${#replace_inflow_temperature} "
char_initial_temperature="${#replace_initial_temperature} "
char_cylinder="${#replace_cylinder} "

#account for change in characters of character count too. ie 2 to 1. 
char_char_inflow_velocity="${#char_inflow_velocity}"
char_char_inflow_temperature="${#char_inflow_temperature}"
char_char_initial_temperature="${#char_initial_temperature}"
char_char_cylinder="${#char_cylinder}"

# this line appears to work for 1-99 change in characters... hopefully is fullproof. May need revision
char_change=$((char_inflow_velocity+char_inflow_temperature+char_initial_temperature+char_cylinder-23-2-16-2-26-2-27-2+char_char_inflow_velocity+char_char_inflow_temperature+char_char_initial_temperature+char_char_cylinder))

# need to adjust line after header to reflect change in num characters. 

echo "Updating inflow BC and cylinder BC"
sed "s,23 replace_inflow_velocity,$char_inflow_velocity$replace_inflow_velocity,g" cylinder.smd > geom_link.smd
sed -i "s,16 replace_cylinder,$char_cylinder$replace_cylinder,g" geom_link.smd
sed -i "s,26 replace_inflow_temperature,$char_inflow_temperature$replace_inflow_temperature,g" geom_link.smd
sed -i "s,27 replace_initial_temperature,$char_initial_temperature$replace_initial_temperature,g" geom_link.smd

header_line=$(awk 'END{print}' geom_link.smd)
array_1=($header_line)
num_1=${array_1[0]}
num_2=${array_1[1]}
num_1=$((num_1+char_change))
num_2=$((num_2+char_change))
header_line_2="$num_1 $num_2                                                                                                                    "

sed -i "s,$header_line,$header_line_2,g" geom_link.smd
