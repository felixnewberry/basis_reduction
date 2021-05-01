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

# makes sense to calculate things once here and simiplify expression in geom_bc file


vel_1=`echo $sigma_1\*$Y1 | bc -l`
vel_2=`echo 2\*$pi\*10 | bc -l`
vel_3=`echo $sigma_2\*$Y2 | bc -l`
vel_4=`echo 10\*$pi\*10 | bc -l`

#replace_inflow_velocity="156.8*(1+$sigma_1*$Y1*cos(2*$pi*x2*10)+$sigma_2*$Y2*cos(10*$pi*x2*10))"
replace_inflow_velocity="156.8*(1+($vel_1)*cos($vel_2*x2)+($vel_3)*cos($vel_4*x2))"

#replace_inflow_velocity="exp(5)" # works

#replace_inflow_velocity="200*exp(-((x2)^2))" # works

#replace_inflow_velocity="200*exp(-((x2)^2)*5000)" #works!
#replace_inflow_velocity="200*exp(-(x2^2)/(2*0.01^2))" #blows up


# 181 or 161 (depending on viscosity 1.81e-5 or 1.61e-4) - Re=1000000

# This breaks the velocity inlet too
#replace_inflow_velocity="156.8+$curve_max*exp(-(x2-$curve_mid)^2/(2*$curve_sd^2))"

temp_inflow_1=`echo 1/\(2\*$curve_sd^2\) | bc -l`

#replace_inflow_temperature="300+$curve_max*exp(-((x2)^2)*5000)" #works!
replace_inflow_temperature="300+$curve_max*exp(-((x2-($curve_mid))^2)*$temp_inflow_1)"

#replace_inflow_temperature="300"
#replace_inflow_temperature="305"
#replace_inflow_temperature="300+5*(cos(x2*$pi*5))^2"
#replace_initial_temperature="300+5*(cos(x2*$pi*5))^2"

replace_initial_temperature="300"
#replace_inflow_temperature="300+$curve_max*exp(-(x2-$curve_mid)^2/(2*$curve_sd^2))"
#replace_initial_temperature="300+$curve_max*exp(-(x2-$curve_mid)^2/(2*$curve_sd^2))"

temp_cy=`echo 0.1+$sigma_3\*$Y3/2 | bc -l`
#replace_cylinder="-50*exp(-(0.1+$sigma_3*$Y3)*$cos_theta/2)"
replace_cylinder="-50*exp(-$temp_cy*$cos_theta)"


#replace_cylinder="10*$cos_theta"
# checked that cos_theta does indeed work

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
sed "s,23 replace_inflow_velocity,$char_inflow_velocity$replace_inflow_velocity,g" 2d_cylinder.smd > geom_link.smd
sed -i "s,16 replace_cylinder,$char_cylinder$replace_cylinder,g" geom_link.smd
sed -i "s,26 replace_inflow_temperature,$char_inflow_temperature$replace_inflow_temperature,g" geom_link.smd
sed -i "s,27 replace_initial_temperature,$char_initial_temperature$replace_initial_temperature,g" geom_link.smd

#</header>
#header_line=tail -1 geom_link.smd
#header_line=awk 'END{print}' geom_link.smd

header_line=$(awk 'END{print}' geom_link.smd)
array_1=($header_line)
num_1=${array_1[0]}
num_2=${array_1[1]}
num_1=$((num_1+char_change))
num_2=$((num_2+char_change))
header_line_2="$num_1 $num_2                                                                                                                    "

sed -i "s,$header_line,$header_line_2,g" geom_link.smd
